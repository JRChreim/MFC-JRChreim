!> energies (6-eqn to 4-eqn) equilibrium through an infinitely fast (algebraic)
!> procedure.

#:include 'macros.fpp'

module m_phase_change

#ifndef MFC_POST_PROCESS

    use m_derived_types        !< Definitions of the derived types

    use m_helper               !< calling module to convert

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use ieee_arithmetic

    implicit none

    private; public :: s_initialize_phasechange_module, &
 s_infinite_relaxation_k

    !> @name Parameters for the first order transition phase change
    !> @{
    integer, parameter :: max_iter = 1e8_wp        !< max # of iterations
    real(wp), parameter :: pCr = 4.94e7_wp   !< Critical water pressure
    real(wp), parameter :: TCr = 385.05_wp + 273.15_wp  !< Critical water temperature
    real(wp), parameter :: mixM = 1.0e-8_wp !< threshold for 'mixture cell'. If Y < mixM, phase change does not happen
    integer, parameter :: lp = 1    !< index for the liquid phase of the reacting fluid
    integer, parameter :: vp = 2    !< index for the vapor phase of the reacting fluid
    !> @}

    !> @name Gibbs free energy phase change parameters
    !> @{
    real(wp) :: A, B, C, D
    !> @}

    !$acc declare create(max_iter,pCr,TCr,mixM,lp,vp,A,B,C,D)

contains

    !>  The purpose of this subroutine is to initialize the phase change module
        !!      by setting the parameters needed for phase change and
        !!      selecting the phase change module that will be used
        !!      (pT- or pTg-equilibrium)
    subroutine s_initialize_phasechange_module()

        ! variables used in the calculation of the saturation curves for fluids 1 and 2
        A = (gs_min(lp)*cvs(lp) - gs_min(vp)*cvs(vp) &
             + qvps(vp) - qvps(lp))/((gs_min(vp) - 1.0_wp)*cvs(vp))

        B = (qvs(lp) - qvs(vp))/((gs_min(vp) - 1.0_wp)*cvs(vp))

        C = (gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
            /((gs_min(vp) - 1.0_wp)*cvs(vp))

        D = ((gs_min(lp) - 1.0_wp)*cvs(lp)) &
            /((gs_min(vp) - 1.0_wp)*cvs(vp))

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to calculate the solution to the selected relaxation system
        if ( .not. any((/1, 4, 5, 6/) == relax_model)) then
            call s_mpi_abort('relaxation solver was not set!')
        end if

    end subroutine s_initialize_phasechange_module !-------------------------------

    !>  This subroutine is created to activate either the pT- (N fluids) or the
        !!      pTg-equilibrium (2 fluids for g-equilibrium)
        !!      model, also considering mass depletion, depending on the incoming
        !!      state conditions.
        !!  @param q_cons_vf Cell-average conservative variables
    subroutine s_infinite_relaxation_k(q_cons_vf) ! ----------------
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(wp) :: pS, pSOV, pSSL !< equilibrium pressure for mixture, overheated vapor, and subcooled liquid
        real(wp) :: TS, TSatOV, TSatSL, TSOV, TSSL !< equilibrium temperature for mixture, overheated vapor, and subcooled liquid. Saturation Temperatures at overheated vapor and subcooled liquid
        real(wp) :: rhoe, rhoeT, dynE, rhos !< total internal energies (different calculations), kinetic energy, and total entropy
        real(wp) :: rho, rM, m1, m2 !< total density, total reacting mass, individual reacting masses
        logical :: TR

        !$acc declare create(pS, pSOV, pSSL, TS, TSatOV, TSatSL, TSOV, TSSL, rhoe, rhoeT, dynE, rhos, rho, rM, m1, m2, TR)

        real(wp), dimension(num_fluids) :: p_infOV, p_infpT, p_infSL, alphak, alpharhoe0k, m0k, rhok, Tk

        !< Generic loop iterators
        integer :: i, j, k, l

        !$acc declare create(p_infOV, p_infpT, p_infSL, m0k, rhok, Tk)

        ! assigning value to the global parameter
        max_iter_pc_ts = 0

        ! starting equilibrium solver
        !$acc parallel loop collapse(3) gang vector default(present) private(p_infOV, p_infpT, p_infSL, alphak, m0k, rhok, Tk, pS, pSOV, pSSL, TS, TSatOV, TSatSL, TSOV, TSSL, rhoe, rhoeT, dynE, rhos, rho, rM, m1, m2, TR)
        do j = 0, m
            do k = 0, n
                do l = 0, p
                    ! trigger for phase change. This will be used for checking many conditions through the code
                    TR = .true.

                    ! computing mixture density, volume fraction, and internal energy, so as saving original variables
                    ! in case phase change is cancelled.
                    rhoeT = 0.0_wp
                    !$acc loop seq
                    do i = 1, num_fluids
                        ! initial volume fraction
                        alphak(i) = q_cons_vf(i + advxb - 1)%sf(j, k, l)
                        
                        ! initial partial density
                        m0k(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)

                        ! calculating the total internal energy such that the energy-fraction for each of the
                        ! fluids can be proportionally distributed when the sum of the internal energies differs from the
                        if ( ieee_is_nan( m0k(i) ) ) then
                            print *, 'partial densities', m0k(i)
                        end if

                        if (model_eqns == 3) then
                            ! initial volume fraction
                            alpharhoe0k(i) = q_cons_vf(i + intxb - 1)%sf(j, k, l)
                            
                            ! total mixture energy by the summation of the internal energy equations
                            rhoeT = rhoeT + alpharhoe0k(i)
                        end if
                    end do
                    
                    ! Mixture density
                    rho = sum(m0k)

                    ! kinetic energy as an auxiliary variable to the calculation of the total internal energy
                    dynE = 0.0_wp
                    !$acc loop seq
                    do i = momxb, momxe
                        dynE = dynE + 5.0e-1_wp*q_cons_vf(i)%sf(j, k, l)**2 / max(rho, sgm_eps)
                    end do

                    ! calculating the internal mixture energy that MUST be preserved throughout pT- and pTg-relaxation procedures
                    ! This calulation is performed as the total energy minus the kinetic one as energy it is preserved at discontinuities
                    rhoe = q_cons_vf(E_idx)%sf(j, k, l) - dynE

                    ! calculating the total reacting mass for the phase change process. By hypothesis, this should not change
                    ! throughout the phase-change process. In ANY HYPOTHESIS, UNLESS correcting them artifically
                    rM = m0k(lp) + m0k(vp)

                    ! correcting possible negative mass fraction values
                    ! at this time, TR is updated if phase change needs to be stopped
                    call s_correct_partial_densities(2, alphak, alpharhoe0k, m0k, rM, rho, TR, i, j, k, l)

                    ! if phase change is still necessary
                    if (TR) then
                        ! (old) p-equilibrium
                        if (relax_model == 1) then
                            call s_old_infinite_p_relaxation_k(j, k, l, alphak, alpharhoe0k, m0k, pS, rhoe, Tk)                            
                        ! p-equilibrium
                        elseif (relax_model == 4) then
                            call s_infinite_p_relaxation_k(j, k, l, alphak, alpharhoe0k, m0k, pS, rho, rhoe, rhoeT, Tk)
                        ! pT-equilibrium
                        elseif (relax_model == 5) then
                            ! for this case, MFL cannot be either 0 or 1, so I chose it to be 2
                            call s_infinite_pt_relaxation_k(j, k, l, m0k, 2, pS, p_infpT, rhoe, rM, TS)
                            Tk = spread(TS, 1, num_fluids)
                        ! pT-pTg equilibrium
                        elseif (relax_model == 6) then

                            ! pT-equilibrium as rhe initial condition
                            call s_infinite_pt_relaxation_k(j, k, l, m0k, 2, pS, p_infpT, rhoe, rM, TS)
                            Tk = spread(TS, 1, num_fluids)

                            ! updating the densities and volume fractions used for thresholds
                            rhok = (pS + ps_inf)/((gs_min - 1)*cvs*Tk)
            
                            ! new volume fractions, after partial densities and p- or pT-equilibrium
                            alphak = m0k / rhok

                            ! 1 - model activation, 1st order transition (p,T) <= (pCr, TCr)
                            ! if ( ( pS < pCr ) .and. ( pS > 0 ) .and. &
                            if ( ( pS < pCr ) .and. &
                            ! 2.1 Homogeneous pTg-equilibrium criterium
                            ( ( ( pS < 0 ) .and. ( pS + minval(p_infpT) > 0.0_wp ) ) &
                            .or. &
                            ! 2.2. Heterogeneous pTg-equilibrium.
                            ( (alphak(lp) > palpha_eps) .and. (alphak(vp) > palpha_eps) ) ) &
                            ) then
                            ! then
                                ! updating m1 and m2 AFTER correcting the partial densities. These values are 
                                ! stored to be retrieved in case the final state is a mixture of fluids
                                m1 = m0k(lp)
                                m2 = m0k(vp)

                                ! checking if fluid is either subcoooled liquid or overheated vapor (NOT metastability)

                                ! overheated vapor
                                ! depleting the mass of liquid and tranferring the total mass to vapor
                                m0k(lp) = mixM*rM
                                m0k(vp) = (1.0_wp - mixM)*rM

                                ! calling pT-equilibrium for overheated vapor, which is MFL = 0
                                call s_infinite_pt_relaxation_k(j, k, l, m0k, 0, pSOV, p_infOV, rhoe, rM, TSOV)

                                ! calculating Saturation temperature
                                call s_TSat(pSOV, TSatOV, TSOV)

                                ! subcooled liquid 
                                ! tranferring the total mass to liquid and depleting the mass of vapor
                                m0k(lp) = (1.0_wp - mixM)*rM
                                m0k(vp) = mixM*rM

                                ! calling pT-equilibrium for subcooled liquid, which is MFL = 1                       
                                call s_infinite_pt_relaxation_k(j, k, l, m0k, 1, pSSL, p_infSL, rhoe, rM, TSSL)

                                ! calculating Saturation temperature
                                call s_TSat(pSSL, TSatSL, TSSL)

                                ! checking the conditions for overheated vapor
                                if (TSOV > TSatOV) then

                                    ! Assigning pressure and temperature
                                    pS = pSOV ; TS = TSOV

                                    ! correcting the liquid and vapor partial densities
                                    m0k(lp) = mixM*rM
                                    m0k(vp) = (1.0_wp - mixM)*rM

                                ! checking the conditions for subcooled liquid
                                elseif (TSSL < TSatSL) then

                                    ! Assigning pressure and temperature
                                    pS = pSSL ; TS = TSSL

                                    ! correcting the liquid and vapor partial densities
                                    m0k(lp) = (1.0_wp - mixM)*rM
                                    m0k(vp) = mixM*rM

                                ! if not, mixture of fluids. Starting phase change (pTg)
                                else
                                    ! returning partial pressures to what they were after the partial density correction 
                                    m0k(lp) = m1
                                    m0k(vp) = m2

                                    ! pTg-relaxation
                                    call s_infinite_ptg_relaxation_k(j, k, l, alphak, alpharhoe0k, m0k, pS, p_infpT, rho, rhoe, rM, TR, TS)
                                    ! if no pTg happens, the solver will return to the hyperbolic state variables
                                    if ( TR .eqv. .false. ) then
                                        !$acc loop seq
                                        do i = 1, num_fluids
                                            ! returning partial densities to what they were previous to any relaxation scheme.
                                            m0k(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)
                                        end do 
                                        ! cycles the innermost loop to the next iteration
                                        cycle 
                                    end if

                                end if
                                Tk = spread(TS, 1, num_fluids)
                            else
                                !$acc loop seq
                                do i = 1, num_fluids
                                    ! returning partial densities to what they were previous to any relaxation scheme.
                                    m0k(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l) 
                                end do
                                ! cycles the innermost loop to the next iteration
                                cycle 
                            end if
                        end if
                    else
                        !$acc loop seq
                        do i = 1, num_fluids
                            ! returning partial densities to what they were previous to any relaxation scheme.
                            m0k(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)
                        end do
                    end if
                    ! updating conservative variables after the any relaxation procedures
                    call update_conservative_vars( j, k, l, m0k, pS, q_cons_vf, Tk )
                end do
            end do
        end do
    end subroutine s_infinite_relaxation_k ! ----------------

    !>  This auxiliary subroutine is created to activate the pT-equilibrium for N fluids
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param pS equilibrium pressure at the interface
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rhoe mixture energy
    subroutine s_infinite_p_relaxation_k(j, k, l, alphaik, alpharhoeik, mik, pS, rho, rhoe, rhoeT, Tk)
        !$acc routine seq

        ! initializing variables
        real(wp), intent(IN) :: rho, rhoe, rhoeT
        real(wp), intent(OUT) :: pS
        real(wp), dimension(num_fluids), intent(OUT) :: Tk
        real(wp), dimension(num_fluids), intent(in)  :: alphaik, alpharhoeik, mik
        integer, intent(IN) :: j, k, l

        real(wp) :: fp, fpp, pO !< variables for the Newton Solver
        real(wp) :: Econst, gamma, pi_inf, mQ, TS !< auxiliary variables
        real(wp), dimension(num_fluids) :: alpha0k, alpharhoe0k, m0k
        real(wp), dimension(num_fluids) :: alphak, alpharhoek, rhok
        character(20) :: nss, pSs, Econsts

        integer :: i, na, ns, nsL !< generic loop iterators

        ! think how to solve pK, alphak, alpharhoek, given the energy mismatechs at the discontinuities 
        alpha0k = alphaik

        ! Re-distributing the initial internal energy values such that rhoe = rhoeT, eventually.
        ! alpharhoe0k does the paper of alpharhoeik, but it is its value adjusted to match rhoe = rhoeT
        alpharhoe0k = alpharhoeik ! * rhoe / rhoeT

        ! distributing the partial density
        m0k = mik

        ! Numerical correction of the volume fractions
        IF (mpp_lim) THEN
            DO i = 1, num_fluids
                IF ((m0k(i) < 0.0_wp) .OR. (alpha0k(i) < 0.0_wp)) THEN
                    alpha0k(i)       = 0.0_wp
                    alpharhoe0k(i)   = 0.0_wp
                    m0k(i)           = 0.0_wp
                END IF
                IF (alpha0k(i) > 1.0_wp) alpha0k(i) = 1.0_wp
                ! IF (m0k(i)/sum(m0k) > 1.0_wp) m0k(i) = sum(m0k)
                ! IF (alpharhoe0k(i)/sum(alpharhoe0k) > 1.0_wp) alpharhoe0k(i) = sum(alpharhoe0k)
            END DO
            alpha0k = alpha0k / sum(alpha0k)
            ! alpharhoe0k = alpharhoe0k / sum(alpharhoe0k)
        END IF

        ! initial conditions for starting the solver. For pressure, as long as the initial guess
        ! is in (-min(gs_min*ps_inf), +infty), a solution should be able to be found.
        pS = 1.0e4_wp

        ! internal energies - first estimate
        alpharhoek = alpharhoe0k

        ! counter for the outer loop
        nsL = 0

        do while  (( ( abs(   sum( alpharhoek ) - sum( alpharhoe0k ) ) > ptgalpha_eps ) &
               .and. ( abs( ( sum( alpharhoek ) - sum( alpharhoe0k ) ) / sum( alpharhoe0k ) ) > ptgalpha_eps ) ) &
               .or.  ( nSL == 0 ) )
            ! increasing counter
            nsL = nsL + 1

            ! Variable to check the energy constraint before initializing the p-relaxation procedure. This ensures
            ! global convergence will be estabilished
            Econst = sum( (gs_min - 1.0_wp)*(alpharhoek - m0k*qvs) / (gs_min*ps_inf - minval(ps_inf)) )

#ifndef MFC_OpenACC
            ! energy constraint for the p-equilibrium
            if ((minval(ps_inf) > 0) .and. (Econst <= 1.0_wp) .or. (nsL > max_iter)) then

                call s_real_to_str(Econst, Econsts)
                call s_mpi_abort('Solver for the p-relaxation solver failed (m_phase_change, s_infinite_p_relaxation_k) &
&                   . Please, check energy constraint. Econst ~'//Econsts//'. Aborting!')

            end if           
#endif

            ! if the energy constraint is satisfied, we start the newton solver
            ! counter
            ns = 0

            ! Newton solver for p-equilibrium. ns <= 1, is to ensure the internal energy correction happens at least once.
            ! in the loosely coupled algorithm. This is different than the pT-equilibrium case, in which no energy correction is needed.
            ! A solution is found when f(p) = 1
            fp = 0.0_wp
            do while (((abs(fp - 1.0_wp) > ptgalpha_eps)) .or. (ns <= 1))
                ! increasing counter
                ns = ns + 1

                ! updating old pressure
                pO = pS

                ! updating functions used in the Newton's solver. f(p)
                fp = sum( (gs_min - 1.0_wp)*(alpharhoek - m0k*qvs) / (pO + gs_min*ps_inf) )

                ! updating functions used in the Newton's solver. f'(p)
                fpp = sum( - 1.0_wp * (gs_min - 1.0_wp)*(alpharhoek - m0k*qvs) / ((pO + gs_min*ps_inf)**2) )

                ! updating the relaxed pressure
                pS = pO + ((1.0_wp - fp)/fpp)/(1.0_wp - (1.0_wp - fp + abs(1.0_wp - fp))/(2.0_wp*fpp*(pO + minval(gs_min*ps_inf))))

                ! checking if pS is within expected bounds
                if ( ((pS <= -1.0_wp*minval(gs_min*ps_inf)) .or. (ieee_is_nan(pS)) .or. (ns > max_iter)) ) then

                    if (proc_rank == 0) then

                        ! call s_tattletale((/ 0.0_wp,  0.0_wp/), reshape((/ 0.0_wp,  0.0_wp,  0.0_wp,  0.0_wp/), (/2, 2/)) &
                        !                 , j, (/ 0.0_wp,  0.0_wp,  0.0_wp,  0.0_wp/), k, l, rhoeT, ps_inf, pS, (/pS - pO, pS + pO/) &
                        !                 , rhoe, q_cons_vf, 0.0_wp)
                    end if

                    call s_real_to_str(pS, pSs)
                    call s_int_to_str(ns, nss)
                    call s_mpi_abort('Solver for the p-relaxation failed (m_phase_change, s_infinite_p_relaxation_k). &
                    &   pS ~'//pSs//'. ns = '//nss//'. Aborting!')

                end if
            end do

            ! updating fluid variables, together with the relaxed pressure, in a loosely coupled procedure
            ! volume fractions
            alphak = (gs_min - 1.0_wp)*(alpharhoek - m0k*qvs) / (pS + gs_min*ps_inf)

            ! internal energies               
            alpharhoek = alpharhoe0k - ( pS + pS ) * (alphak - alpha0k) / 2
        end do

        ! Mixture-total-energy correction ==================================

        ! The mixture-total-energy correction of the mixture pressure P is not necessary here
        ! because the primitive variables are directly recovered later on by the conservative
        ! variables (see s_convert_conservative_to_primitive_variables called in s_compute_rhs).
        ! However, the internal-energy equations should be reset with the corresponding mixture
        ! pressure from the correction. This step is carried out below.
        gamma  = sum( alphak * 1.0_wp * 1.0_wp / ( gs_min - 1.0_wp ) )
        pi_inf = sum( alphak * gs_min * ps_inf / ( gs_min - 1.0_wp ) )
        mQ     = sum( m0k * qvs )    

        ! restarted pressure
        pS = (rhoe - pi_inf - mQ) / gamma

        ! densities
        rhok = m0k / alphak

        ! (NOT common) temperature
        Tk = (pS + ps_inf)/((gs_min - 1)*cvs*rhok)

        ! updating maximum number of iterations
        max_iter_pc_ts = maxval((/max_iter_pc_ts, ns/))
    end subroutine s_infinite_p_relaxation_k ! -----------------------

    subroutine s_old_infinite_p_relaxation_k(j, k, l, alpha0k, alpharhoe0k, m0k, pS, rhoe, Tk)
        ! Description: The purpose of this procedure is to infinitely relax
        !              the pressures from the internal-energy equations to a
        !              unique pressure, from which the corresponding volume
        !              fraction of each phase are recomputed. For conservation
        !              purpose, this pressure is finally corrected using the
        !              mixture-total-energy equation.

        ! make sure to clean all p_relaxation codes, removing unnecessary variables. It would be goo0d not to have q_cons_vf as an input variable anywhere

        ! Relaxed pressure, initial partial pressures, function f(p) and its partial
        ! derivative df(p), isentropic partial density, sum of volume fractions,
        ! mixture density, dynamic pressure, surface energy, specific heat ratio
        ! function, liquid stiffness function (two variations of the last two
        ! initializing variables
        real(wp), intent(IN) :: rhoe
        real(wp), intent(OUT) :: pS
        real(wp), dimension(num_fluids), intent(OUT) :: Tk
        real(wp), dimension(num_fluids), intent(in) :: alpha0k, alpharhoe0k, m0k

        integer, intent(IN) :: j, k, l

        real(wp) :: fp, fpp, mQ, drhodp, den, num, pO !< variables for the Newton Solver
        real(wp) :: gamma, pi_inf !< auxiliary variables
        real(wp), dimension(num_fluids) :: alphak, alpharhoek, mk, pk, rhok 
        character(20) :: nss, pSs

        integer :: i, ns !< generic loop iterators

        ! initializing the partial energies
        alphak       = alpha0k
        alpharhoek   = alpharhoe0k
        mk           = m0k

        ! Numerical correction of the volume fractions
        IF (mpp_lim) THEN
            DO i = 1, num_fluids
                IF ((mk(i) < 0.0_wp) .OR. (alphak(i) < 0.0_wp)) THEN
                    alphak(i)       = 0.0_wp
                    alpharhoek(i)   = 0.0_wp
                    mk(i)           = 0.0_wp
                END IF
                IF (alphak(i) > 1.0_wp) alphak(i) = 1.0_wp
            END DO
            alphak = alphak / sum(alphak)
        END IF

        ! Initial state
        DO i = 1, num_fluids
            IF (alphak(i) > sgm_eps) THEN
                pk(i) = ( ( alpharhoek(i) - mk(i) * qvs(i) )/ alphak(i) - fluid_pp(i)%pi_inf) / fluid_pp(i)%gamma               
                ! Physical pressure?
                IF (pk(i) <= -(1.0_wp - ptgalpha_eps)*ps_inf(i) + ptgalpha_eps) pk(i) = -(1.0_wp - ptgalpha_eps)*ps_inf(i) + ptgalpha_eps
            ELSE
                pk(i) = 0.0_wp
            END IF
        END DO

        pS = sum( alphak * pk )

        ! Iterative process for relaxed pressure determination
        ns    = 0
        fp    = 2.0_wp * ptgalpha_eps
        fpp   = 1d9
        rhok  = 0.0_wp

        ! Note that a solution is found when f(p) = 1
        ! keep an eye on this
        DO WHILE (abs(fp) > ptgalpha_eps .or. (ns == 0))
            
            ! setting old pressure
            pO = pS
            
            ! updating pressure
            pS = pO - fp / fpp

            ! Newton-Raphson method
            fp  = -1.0_wp
            fpp = 0.0_wp
            DO i = 1, num_fluids
                ! Physical pressure?
                IF (pS <= -(1.0_wp - ptgalpha_eps)*ps_inf(i) + ptgalpha_eps) pS = -(1.0_wp - ptgalpha_eps)*ps_inf(i) + ptgalpha_eps
                ! calulating fp and fpp
                IF (alphak(i) > sgm_eps) THEN
                    num       = gs_min(i) * (pS + ps_inf(i))
                    den       = num + pk(i) - pS
                    rhok(i)   = mk(i) / MAX(alphak(i),sgm_eps) * num / den
                    drhodp    = mk(i) / MAX(alphak(i),sgm_eps) * gs_min(i) * ( pk(i) + ps_inf(i) ) / ( den**2 )
                    fp        = fp  + mk(i) / rhok(i)
                    fpp       = fpp - mk(i) * drhodp / (rhok(i)**2)
                END IF
            END DO

            ! Convergence?
            ns = ns + 1
#ifndef MFC_OpenACC
                ! energy constraint for the p-equilibrium
                ! checking if pressure is within expected bounds
                if ((pS <= -1.0_wp*minval(gs_min*ps_inf)) .or. (ieee_is_nan(pS)) .or. (ns > max_iter)) then

                    if (proc_rank == 0) then

                        ! call s_tattletale((/ 0.0_wp,  0.0_wp/), reshape((/ 0.0_wp,  0.0_wp,  0.0_wp,  0.0_wp/), (/2, 2/)) &
                        !                 , j, (/ 0.0_wp,  0.0_wp,  0.0_wp,  0.0_wp/), k, l, 0.0_wp, ps_inf, pS, (/pS - pO, pS + pO/) &
                        !                 , rhoe, q_cons_vf, 0.0_wp)
                    end if

                    call s_real_to_str(pS, pSs)
                    call s_int_to_str(ns, nss)
                    call s_mpi_abort('Solver for the p-relaxation failed (m_phase_change, s_old_infinite_p_relaxation_k). &
                    &   pS ~'//pSs//'. ns = '//nss//'. Aborting!')

                end if
#endif
        end do

        ! Cell update of the volume fraction - this is not needed to be outputed, as it will be later fixed by a more general subroutine
        DO i = 1, num_fluids
            IF (alphak(i) .GT. sgm_eps) alphak(i) = mk(i) / rhok(i)
        END DO                    

        ! Mixture-total-energy correction ==================================

        ! The mixture-total-energy correction of the mixture pressure P is not necessary here
        ! because the primitive variables are directly recovered later on by the conservative
        ! variables (see s_convert_conservative_to_primitive_variables called in s_compute_rhs).
        ! However, the internal-energy equations should be reset with the corresponding mixture
        ! pressure from the correction. This step is carried out below.
        gamma  = sum( alphak * 1.0_wp * 1.0_wp / ( gs_min - 1.0_wp ) )
        pi_inf = sum( alphak * gs_min * ps_inf / ( gs_min - 1.0_wp ) )
        mQ     = sum( mk * qvs )    

        pS = (rhoe - pi_inf - mQ) / gamma

        ! (NOT common) temperature
        Tk = (pS + ps_inf)/((gs_min - 1)*cvs*rhok)

        ! updating maximum number of iterations
        max_iter_pc_ts = maxval((/max_iter_pc_ts, ns/))

    end subroutine s_old_infinite_p_relaxation_k ! -----------------------

    !>  This auxiliary subroutine is created to activate the pT-equilibrium for N fluids
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param MFL flag that tells whether the fluid is pure gas (0), pure liquid (1), or a mixture (2)
        !!  @param pS equilibrium pressure at the interface
        !!  @param p_infpT stiffness for the participating fluids under pT-equilibrium
        !!  @param rM sum of the reacting masses
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rhoe mixture energy
        !!  @param TS equilibrium temperature at the interface
    subroutine s_infinite_pt_relaxation_k(j, k, l, m0k, MFL, pS, p_infpT, rhoe, rM, TS)

#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_infinite_pt_relaxation_k
#else
        !$acc routine seq
#endif

        ! initializing variables
        real(wp), intent(out) :: pS, TS
        real(wp), dimension(num_fluids), intent(out) :: p_infpT
        real(wp), intent(in) :: rhoe, rM
        real(wp), intent(in), dimension(num_fluids) :: m0k
        integer, intent(in) :: j, k, l, MFL
        integer, dimension(num_fluids) :: ig !< flags to toggle the inclusion of fluids for the pT-equilibrium
        real(wp) :: gp, gpp, hp, pO, mCP, mQ !< variables for the Newton Solver
        character(20) :: nss, pSs, Econsts

        integer :: i, ns !< generic loop iterators

        ! auxiliary variables for the pT-equilibrium solver
        mCP = 0.0_wp; mQ = 0.0_wp; p_infpT = ps_inf;

        ! these are slowing the computations significantly. Think about a workaround
        ig = 0

        ! Performing tests before initializing the pT-equilibrium
        !$acc loop seq
        do i = 1, num_fluids

            ! check if all alpha(i)*rho(i) are negative. If so, abort
#ifndef MFC_OpenACC
            ! check which indices I will ignore (no need to abort the solver in this case). Adjust this sgm_eps value for mixture cells
            if( ( m0k(i) .ge.  0.0_wp ) .and. ( m0k(i) - rM * mixM .le. sgm_eps ) ) then

                ig(i) = i

                ! PRINT *, ig

                ! this value is rather arbitrary, as I am interested in MINVAL( ps_inf ) for the solver.
                ! This way, I am ensuring this value will not be selected.
                p_infpT(i) = 2 * MAXVAL( ps_inf )

            end if
#endif

            if ( ( bubbles_euler .eqv. .false. ) .or. ( bubbles_euler .and. (i /= num_fluids) ) ) then
                ! sum of the total alpha*rho*cp of the system                
                mCP = mCP + m0k(i)*cvs(i)*gs_min(i)

                ! sum of the total alpha*rho*q of the system
                mQ = mQ + m0k(i)*qvs(i)
            end if
        end do

        ! Checking energy constraint. In the case we are calculating the possibility of having subcooled liquid or
        ! overheated vapor, the energy constraint might not be satisfied, as are hypothetically transferring all the 
        ! mass from one phase to the other. When this is the case, we simply ignore this possibility, set pS = TS = 0,
        ! and discard the hypothesis. The solver can thus move forward.
        if ((rhoe - mQ - minval(p_infpT)) < 0.0_wp) then

            if ((MFL == 0) .or. (MFL == 1)) then

                ! Assigning zero values for mass depletion cases
                ! pressure
                pS = 0.0_wp

                ! temperature
                TS = 0.0_wp

                return
#ifndef MFC_OpenACC
            else
                if (proc_rank == 0) then

                    ! call s_tattletale((/ 0.0_wp,  0.0_wp/), reshape((/ 0.0_wp,  0.0_wp,  0.0_wp,  0.0_wp/), (/2, 2/)) &
                    !                   , j, (/ 0.0_wp,  0.0_wp,  0.0_wp,  0.0_wp/), k, l, mQ, p_infpT, pS, (/abs(pS - pO), abs(pS - pO)/) &
                    !                   , rhoe, q_cons_vf, TS)

                end if

                call s_real_to_str(rhoe - mQ - minval(p_infpT), Econsts)
                call s_mpi_abort('Solver for the pT-relaxation solver failed (m_phase_change, s_infinite_pt_relaxation_k) &
                &    . Energy constraint~'//Econsts//'. Aborting!')

#endif
            end if
        end if

        ! Maybe improve this condition afterwards. As long as the initial guess is in between -min(ps_inf)
        ! and infinity, a solution should be able to be found.
        pS = 1.0e4_wp

        ! starting counter for the Newton solver
        ns = 0
  
        ! arbitrary value for g(p), used for the newton solver
        gp = 2.0_wp

        ! Newton solver for pT-equilibrium. 1d1 is arbitrary, and ns == 0, to the loop is entered at least once.
        ! Note that a solution is found when gp(p) = 1
        ! keep an eye on this
        ! do while (((abs(gp - 1.0_wp) > ptgalpha_eps) .and. (abs((gp - 1.0_wp)/gp) > ptgalpha_eps/1.0e6_wp)) .or. (ns == 0))
        do while (((abs(gp - 1.0_wp) > ptgalpha_eps)) .or. (ns == 0))

            ! increasing counter
            ns = ns + 1

            ! updating old pressure
            pO = pS

            ! updating functions used in the Newton's solver
            gpp = 0.0_wp; gp = 0.0_wp; hp = 0.0_wp
            !$acc loop seq
            do i = 1, num_fluids

                ! given pS always change, I need ig( i ) and gp to be in here, as it dynamically updates.
                ! Note that I do not need to use p_infpT here, but I will do it for consistency
                if (i /= ig(i)) then

                    gp = gp + (gs_min(i) - 1.0_wp)*m0k(i)*cvs(i) &
                        *(rhoe + pO - mQ)/(mCP*(pO + p_infpT(i)))

                    gpp = gpp + (gs_min(i) - 1.0_wp)*m0k(i)*cvs(i) &
                        *(p_infpT(i) - rhoe + mQ)/(mCP*(pO + p_infpT(i))**2)

                end if

            end do

            hp = 1.0_wp/(rhoe + pO - mQ) + 1.0_wp/(pO + minval(p_infpT))

            ! updating common pressure for the newton solver
            pS = pO + ((1.0_wp - gp)/gpp)/(1.0_wp - (1.0_wp - gp + abs(1.0_wp - gp)) &
                                          /(2.0_wp*gpp)*hp)

            ! check if solution is out of bounds (which I believe it won`t happen given the solver is gloabally convergent.
#ifndef MFC_OpenACC
            if ((pS <= -1.0_wp*minval(p_infpT)) .or. (ieee_is_nan(pS)) .or. (ns > max_iter)) then
                if (proc_rank == 0) then
                    ! call s_tattletale((/0.0_wp, 0.0_wp/), reshape((/0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp/), (/2, 2/)) &
                    !                   , j, (/0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp/), k, l, mQ, p_infpT, pS, (/pS - pO, pS + pO/) &
                    !                   , rhoe, q_cons_vf, TS)
                end if

                call s_real_to_str(pS, pSs)
                call s_int_to_str(nS, nss)
                call s_mpi_abort('Solver for the pT-relaxation failed (m_phase_change, s_infinite_pt_relaxation_k). &
                &   pS ~'//pSs//'. ns = '//nss//'. Aborting!')

            end if
#endif
        end do

        ! common temperature
        TS = (rhoe + pS - mQ)/mCP

        ! updating maximum number of iterations
        max_iter_pc_ts = maxval((/max_iter_pc_ts, ns/))

    end subroutine s_infinite_pt_relaxation_k ! -----------------------

    !>  This auxiliary subroutine is created to activate the pTg-equilibrium for N fluids under pT
        !!      and 2 fluids under pTg-equilibrium. There is a final common p and T during relaxation
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param pS equilibrium pressure at the interface
        !!  @param p_infpT stiffness for the participating fluids under pT-equilibrium
        !!  @param rhoe mixture energy
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param TS equilibrium temperature at the interface
    subroutine s_infinite_ptg_relaxation_k(j, k, l, alphak, alpharhoe0k, m0k, pS, p_infpT, rho, rhoe, rM, TR, TS)

#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_infinite_ptg_relaxation_k
#else
        !$acc routine seq
#endif

        real(wp), dimension(num_fluids), intent(INOUT) :: alphak, alpharhoe0k, m0k
        real(wp), dimension(num_fluids), intent(IN) :: p_infpT
        real(wp), intent(INOUT) :: pS, TS, rM
        real(wp), intent(IN) :: rho, rhoe
        integer, intent(IN) :: j, k, l
        logical, intent(INOUT) :: TR
        real(wp), dimension(num_fluids) :: p_infpTg
        real(wp), dimension(2, 2) :: Jac, InvJac, TJac
        real(wp), dimension(2) :: R2D, DeltamP
        real(wp), dimension(3) :: Oc
        real(wp) :: Om, OmI ! underrelaxation factor
        real(wp) :: mCP, mCPD, mCVGP, mCVGP2, mQ, mQD, TSat ! auxiliary variables for the pTg-solver
        character(20) :: nss, pSs, Econsts

        !< Generic loop iterators
        integer :: i, ns

        ! checking if homogeneous cavitation is expected. If yes, transfering an amount of mass to the depleted (liquid,
        ! for the moment) phase, and then let the algorithm run. 
        ! checking if homogeneous cavitation is possible

        ! is the fluid at a metastable state with enough 'energy' for phase change to happen?
        if ((pS < 0.0_wp) .and. (rM > (rhoe - gs_min(lp)*ps_inf(lp)/(gs_min(lp) - 1.0e-1_wp))/qvs(lp))) then

            ! transfer a bit of mass to the deficient phase, enforce phase change
            call s_correct_partial_densities(1, alphak, alpharhoe0k, m0k, rM, rho, TR, i, j, k, l)

        ! the metastable state is not enough to sustain phase change
        elseif (pS < 0.0_wp) then
            
            ! cancel any phase-change updates.
            TR = .false.

            ! ends the execution of this function and returns control to the calling function
            return

        ! if not homogeneous, then heterogeneous. Thus, setting up an arbitrary initial condition in case the one from
        ! the p(T)-equilibrium solver could lead to numerical issues
        elseif ((pS < 1.0d-1) .and. (pS >= 0.0_wp)) then
            ! improve this initial condition
            pS = 1.0e4_wp

        end if

        ! Relaxation factor. This value is rather arbitrary, with a certain level of self adjustment.
        OmI = 1.0e-1_wp
        ! Critical relaxation factors, for variable sub-relaxation
        Oc(1) = OmI; Oc(2) = OmI; Oc(3) = OmI

        R2D(1) = 0.0_wp; R2D(2) = 0.0_wp
        DeltamP(1) = 0.0_wp; DeltamP(2) = 0.0_wp
        ! starting counter for the Newton solver
        ns = 0

        ! Newton solver for pTg-equilibrium. 1d6 is arbitrary, and ns == 0, to the loop is entered at least once.
        do while (((sqrt(R2D(1)**2 + R2D(2)**2) > ptgalpha_eps) &
                   .and. ((sqrt(R2D(1)**2 + R2D(2)**2)/rhoe) > (ptgalpha_eps/1d6))) &
                  .or. (ns == 0))

            ! Updating counter for the iterative procedure
            ns = ns + 1

            ! Auxiliary variables to help in the calculation of the residue
            mCP = 0.0_wp; mCPD = 0.0_wp; mCVGP = 0.0_wp; mCVGP2 = 0.0_wp; mQ = 0.0_wp; mQD = 0.0_wp
            ! Those must be updated through the iterations, as they either depend on
            ! the partial masses for all fluids, or on the equilibrium pressure
            !$acc loop seq
            do i = 1, num_fluids
                if ( ( bubbles_euler .eqv. .false. ) .or. ( bubbles_euler .and. (i /= num_fluids) ) ) then

                    ! sum of the total alpha*rho*cp of the system
                    mCP = mCP + m0k(i)*cvs(i)*gs_min(i)

                    ! sum of the total alpha*rho*q of the system
                    mQ = mQ + m0k(i)*qvs(i)

                    ! These auxiliary variables now need to be updated, as the partial densities now
                    ! vary at every iteration.
                    if ((i /= lp) .and. (i /= vp)) then

                        mCVGP  = mCVGP + m0k(i)*cvs(i)*(gs_min(i) - 1)/(pS + ps_inf(i))

                        mCVGP2 = mCVGP2 + m0k(i)*cvs(i)*(gs_min(i) - 1)/((pS + ps_inf(i))**2)

                        mQD = mQD + m0k(i)*qvs(i)

                        ! sum of the total alpha*rho*cp of the system
                        mCPD = mCPD + m0k(i)*cvs(i)*gs_min(i)

                    end if
                end if
            end do

            ! Checking pressure and energy criteria for the (pT) solver to find a solution
#ifndef MFC_OpenACC
            if ((pS <= -1.0_wp*minval(ps_inf)) .or. ((rhoe - mQ - minval(ps_inf)) < 0.0_wp)) then
                if (proc_rank == 0) then
                    ! call s_tattletale(DeltamP, InvJac, j, Jac, k, l, mQ, ps_inf, pS &
                    !                   , R2D, rhoe, q_cons_vf, TS)
                end if

                call s_real_to_str(rhoe - mQ - minval(ps_inf), Econsts)
                call s_real_to_str(pS, pSs)
                call s_mpi_abort('Solver for the pTg-relaxation failed (m_phase_change, s_infinite_ptg_relaxation_k). &
                &   pS ~'//pSs//'. Econst = '//Econsts//'. Aborting!')

            end if
#endif
            ! calculating the (2D) Jacobian Matrix used in the solution of the pTg-quilibrium model
            call s_compute_jacobian_matrix(InvJac, j, Jac, k, l, m0k, mCPD, mCVGP, mCVGP2, pS, rM, TJac)

            ! calculating correction array for Newton's method
            DeltamP = matmul(InvJac, R2D)

            ! checking if the correction in the mass/pressure will lead to negative values for those quantities
            ! If so, adjust the underrelaxation parameter Om

#ifndef MFC_OpenACC
            ! creating criteria for variable underrelaxation factor
            if (m0k(lp) - Om*DeltamP(1) <= 0.0_wp) then
                Oc(1) = m0k(lp)/(2*DeltamP(1))
            else
                Oc(1) = OmI
            end if
            if (m0k(vp) + Om*DeltamP(1) <= 0.0_wp) then
                Oc(2) = -m0k(vp)/(2*DeltamP(1))
            else
                Oc(2) = OmI
            end if
            if (pS + minval(ps_inf) - Om*DeltamP(2) <= 0.0_wp) then
                Oc(3) = (pS + minval(ps_inf))/(2*DeltamP(2))
            else
                Oc(3) = OmI
            end if
            ! choosing amonst the minimum relaxation maximum to ensure solver will not produce unphysical values
            Om = minval(Oc)
#else
            Om = OmI
#endif
            ! updating two reacting 'masses'. Recall that inert 'masses' do not change during the phase change
            ! liquid
            m0k(lp) = m0k(lp) - Om*DeltamP(1)

            ! gas
            m0k(vp) = m0k(vp) + Om*DeltamP(1)

            ! updating pressure
            pS = pS - Om*DeltamP(2)

            ! calculating residuals, which are (i) the difference between the Gibbs Free energy of the gas and the liquid
            ! and (ii) the energy before and after the phase-change process.
            call s_compute_pTg_residue(j, k, l, m0k, mCPD, mCVGP, mQD, pS, rhoe, rM, R2D)

            ! checking if the residue returned any NaN values
#ifndef MFC_OpenACC
            if ((ieee_is_nan(R2D(1))) .or. (ieee_is_nan(R2D(2))) .or. (ns > max_iter)) then

                print *, 'homogeneous cavitation criteria', pS + minval(p_infpT)

                print *, 'Oc', Oc

                print *, 'Om', Om

                print *, 'pO', pS + Om*DeltamP(2)

                print *, 'pOOmI', pS + OmI*DeltamP(2)

                print *, 'pS', pS

                print *, 'TS', (rhoe + pS - mQ)/mCP

                print *, 'R2D', R2D

                print *, 'DeltamP', DeltamP

                print *, 'jkl', j,k,l

                call s_TSat(pS, TSat, (rhoe + pS - mQ)/mCP)

                print *, 'TSat', TSat 

                call s_int_to_str(ns, nss)
                call s_mpi_abort('Residual for the pTg-relaxation possibly returned NaN values. ns = ' &
                                 //nss//' (m_phase_change, s_infinite_ptg_relaxation_k). Aborting!')

            end if
#endif
        end do

        ! common temperature
        TS = (rhoe + pS - mQ)/mCP

        ! updating maximum number of iterations
        max_iter_pc_ts = maxval((/max_iter_pc_ts, ns/))
    end subroutine s_infinite_ptg_relaxation_k ! -----------------------

    !>  This auxiliary subroutine corrects the partial densities of the REACTING fluids in case one of them is negative
        !!      but their sum is positive. Inert phases are not corrected at this moment
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rM sum of the reacting masses
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
    subroutine s_correct_partial_densities(CT, alpha0k, alpharhoe0k, m0k, rM, rho, TR, i, j, k, l)

#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_correct_partial_densities
#else
        !$acc routine seq
#endif

        !> @name variables for the correction of the reacting partial densities
        !> @{
        real(wp), dimension(num_fluids), intent(INOUT) :: alpha0k, alpharhoe0k, m0k
        real(wp), intent(INOUT) :: rM
        real(wp), intent(IN) :: rho
        logical, intent(INOUT) :: TR
        integer, intent(IN) :: CT, j, k, l
        integer :: i
        !> @}

        ! CT = 0: No Mass correction; ! CT = 1: Reacting Mass correction;
        ! CT = 2: crude correction; else: Total Mass correction
        if (CT == 0) then
            !$acc loop seq
            do i = 1, num_fluids
                ! analysis in terms of mass fraction because this is what is needed for the relaxation process. Volume
                ! fraction does not matter
                if (m0k(i)/rho < mixM) then
                    ! do not continue relaxation
                    TR = .false.
                end if
            end do
        elseif (CT == 1) then
            if (rM < 0.0_wp) then
                ! reacting masses are very negative so as to affect the physics of the problem, so phase change will not be activated
                if ((m0k(lp)/rM < mixM) .or. &
                    (m0k(vp)/rM < mixM)) then
                    
                    ! do not continue relaxation
                    TR = .false.
                ! reacting masses are not as negative so I can disregard them,
                ! expecting no significant changes in the physics of the simulation
                else
                    m0k(lp) = mixM*rM
                    m0k(vp) = mixM*rM

                    ! continue relaxation
                    TR = .true.
                end if
            ! correcting the partial densities of the reacting fluids. In case liquid is negative
            elseif (m0k(lp)/rM < mixM) then

                m0k(lp) = mixM*rM
                m0k(vp) = (1.0_wp - mixM)*rM
                
                ! continue relaxation
                TR = .true.

            ! correcting the partial densities of the reacting fluids. In case vapor is negative
            elseif (m0k(vp)/rM < mixM) then

                m0k(lp) = (1.0_wp - mixM)*rM
                m0k(vp) = mixM*rM

                ! continue relaxation
                TR = .true.

            end if
        elseif (CT == 2) then
            do i = 1, num_fluids
                if (alpha0k(i) < 0 .or. m0k(i)  < 0 ) then
                    alpha0k(i) = 0.0_wp
                    m0k(i) = 0.0_wp
                    if (model_eqns .eq. 3) then
                        alpharhoe0k(i) = 0.0_wp
                    end if 
                end if
            end do
            ! continue relaxation
            TR = .true.
        else
            !$acc loop seq
            do i = 1, num_fluids
                if ((m0k(i)/rho) < mixM) then
                    m0k(i) = mixM*rho
                end if
            end do
            ! continue relaxation
            TR = .true.
        end if

        rM = m0k(lp) + m0k(vp)

    end subroutine s_correct_partial_densities

    !>  This auxiliary subroutine calculates the 2 x 2 Jacobian and, its inverse and transpose
        !!      to be used in the pTg-equilibirium procedure
        !!  @param InvJac Inverse of the Jacobian Matrix
        !!  @param j generic loop iterator for x direction
        !!  @param Jac Jacobian Matrix
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param mCPD  sum of the total alpha*rho*cp
        !!  @param mCVGP auxiliary variable for the calculation of the matrices: alpha*rho*cv*(g-1)/press
        !!  @param mCVGP2 auxiliary variable for the calculation of the matrices: alpha*rho*cv*(g-1)/press^2
        !!  @param pS equilibrium pressure at the interface
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param TJac Transpose of the Jacobian Matrix
    subroutine s_compute_jacobian_matrix(InvJac, j, Jac, k, l, m0k, mCPD, mCVGP, mCVGP2, pS, rM, TJac)

#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_compute_jacobian_matrix
#else
        !$acc routine seq
#endif

        real(wp), dimension(num_fluids), intent(IN) :: m0k
        real(wp), intent(IN) :: pS, mCPD, mCVGP, mCVGP2, rM
        integer, intent(IN) :: j, k, l
        real(wp), dimension(2, 2), intent(OUT) :: Jac, InvJac, TJac
        real(wp) :: TS, dFdT, dTdm, dTdp ! mass of the reacting fluid, total reacting mass, and auxiliary variables

        TS = 1/(rM*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) &
                + m0k(lp)*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                + mCVGP)

        dFdT = &
            -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*log(TS) &
            - (qvps(lp) - qvps(vp)) &
            + cvs(lp)*(gs_min(lp) - 1)*log(pS + ps_inf(lp)) &
            - cvs(vp)*(gs_min(vp) - 1)*log(pS + ps_inf(vp))

        dTdm = -(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                 - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)))*TS**2

        dTdp = (rM*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2 &
                + m0k(lp)*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp))**2 &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2) &
                + mCVGP2)*TS**2

        ! F = (F1,F2) is the function whose roots we are looking for
        ! x = (m1, p) are the independent variables. m1 = mass of the first participant fluid, p = pressure
        ! F1 = 0 is the Gibbs free energy quality
        ! F2 = 0 is the enforcement of the thermodynamic (total - kinectic) energy
        ! dF1dm
        Jac(1, 1) = dFdT*dTdm

        ! dF1dp
        Jac(1, 2) = dFdT*dTdp + TS &
                    *(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)))

        ! dF2dm
        Jac(2, 1) = (qvs(vp) - qvs(lp) &
                     + (cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                     /(m0k(lp)*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                       + rM*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP) &
                     - (m0k(lp)*(cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                        - rM*cvs(vp)*gs_min(vp) - mCPD) &
                     *(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                       - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                     /((m0k(lp)*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                            - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                        + rM*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP)**2))/1
        ! dF2dp
        Jac(2, 2) = (1 + (m0k(lp)*(cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                          - rM*cvs(vp)*gs_min(vp) - mCPD) &
                     *(m0k(lp)*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp))**2 &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2) &
                       + rM*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2 + mCVGP2) &
                     /(m0k(lp)*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                       + rM*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP)**2)/1

        ! intermediate elements of J^{-1}
        InvJac(1, 1) = Jac(2, 2)
        InvJac(1, 2) = -1.0_wp*Jac(1, 2)
        InvJac(2, 1) = -1.0_wp*Jac(2, 1)
        InvJac(2, 2) = Jac(1, 1)

        ! elements of J^{T}
        TJac(1, 1) = Jac(1, 1)
        TJac(1, 2) = Jac(2, 1)
        TJac(2, 1) = Jac(1, 2)
        TJac(2, 2) = Jac(2, 2)

        ! dividing by det(J)
        InvJac = InvJac/(Jac(1, 1)*Jac(2, 2) - Jac(1, 2)*Jac(2, 1))

    end subroutine s_compute_jacobian_matrix

    !>  This auxiliary subroutine computes the residue of the pTg-equilibrium procedure
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param mCPD  sum of the total alpha*rho*cp
        !!  @param mCVGP auxiliary variable for the calculation of the matrices: alpha*rho*cv*(g-1)/press
        !!  @param mQD sum of the total alpha*rho*qv
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param pS equilibrium pressure at the interface
        !!  @param rhoe mixture energy
        !!  @param R2D (2D) residue array
    subroutine s_compute_pTg_residue(j, k, l, m0k, mCPD, mCVGP, mQD, pS, rhoe, rM, R2D)

#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_compute_pTg_residue
#else
        !$acc routine seq
#endif

        real(wp), dimension(num_fluids), intent(IN) :: m0k
        real(wp), intent(IN) :: pS, rhoe, mCPD, mCVGP, mQD, rM
        integer, intent(IN) :: j, k, l
        real(wp), dimension(2), intent(OUT) :: R2D
        real(wp) :: TS !< mass of the reacting liquid, total reacting mass, equilibrium temperature

        ! relaxed temperature
        TS = 1/(rM*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) &
                + m0k(lp)*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                + mCVGP)

        ! Gibbs Free Energy Equality condition (DG)
        R2D(1) = TS*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                     *(1 - log(TS)) - (qvps(lp) - qvps(vp)) &
                     + cvs(lp)*(gs_min(lp) - 1)*log(pS + ps_inf(lp)) &
                     - cvs(vp)*(gs_min(vp) - 1)*log(pS + ps_inf(vp))) &
                 + qvs(lp) - qvs(vp)

        ! Constant Energy Process condition (DE)
        R2D(2) = (rhoe + pS &
                  + m0k(lp)*(qvs(vp) - qvs(lp)) - rM*qvs(vp) - mQD &
                  + (m0k(lp)*(gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
                     - rM*gs_min(vp)*cvs(vp) - mCPD) * TS ) / 1

    end subroutine s_compute_pTg_residue

    ! SUBROUTINE CREATED TO TELL ME WHERE THE ERROR IN THE PT- AND PTG-EQUILIBRIUM SOLVERS IS
    subroutine s_tattletale(DeltamP, InvJac, j, Jac, k, l, mQ, p_infA, pS, R2D, rhoe, q_cons_vf, TS) ! ----------------

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(wp), dimension(2, 2), intent(IN) :: Jac, InvJac
        real(wp), dimension(num_fluids), intent(IN) :: p_infA
        real(wp), dimension(2), intent(IN) :: R2D, DeltamP
        real(wp), intent(IN) :: pS, TS
        real(wp), intent(IN) :: rhoe, mQ
        integer, intent(IN) :: j, k, l
        real(wp) :: rho
        !< Generic loop iterator
        integer :: i

        print *, 'j, k, l', j, k, l

        print *, 'rhoe', rhoe

        print *, 'mQ', mQ

        print *, 'Energy constrain', (rhoe - mQ - minval(p_infA))

        print *, 'R2D', R2D

        print *, 'l2(R2D)', sqrt(R2D(1)**2 + R2D(2)**2)

        print *, 'DeltamP', DeltamP

        print *, 'pS', pS

        print *, '-min(ps_inf)', -minval(p_infA)

        print *, 'TS', TS

        do i = 1, num_fluids

            rho = rho + q_cons_vf(i + contxb - 1)%sf(j, k, l)

        end do

        print *, 'rho', rho

        do i = 1, num_fluids

            print *, 'i', i

            print *, 'alpha_i', q_cons_vf(i + advxb - 1)%sf(j, k, l)

            print *, 'alpha_rho_i', q_cons_vf(i + contxb - 1)%sf(j, k, l)

            print *, 'mq_i', q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

            if (model_eqns == 3) then
                print *, 'internal energies', q_cons_vf(i + intxb - 1)%sf(j, k, l)
            end if

            print *, 'Y_i', q_cons_vf(i + contxb - 1)%sf(j, k, l)/rho

        end do

        print *, 'J', Jac, 'J-1', InvJac

    end subroutine s_tattletale

    ! Newton Solver for the finding the Saturation temperature TSat for a given saturation pressure
    subroutine s_TSat(pSat, TSat, TSIn)

#ifdef _CRAYFTN
        !DIR$ INLINEALWAYS s_TSat
#else
        !$acc routine seq
#endif

        real(wp), intent(OUT) :: TSat
        real(wp), intent(IN) :: pSat, TSIn
        real(wp) :: dFdT, FT, Om !< auxiliary variables
        character(20) :: nss, pSatS, TSatS

        ! Generic loop iterators
        integer :: ns

        ! in case of fluid under tension (p - p_inf > 0, T > 0), or, when subcooled liquid/overheated vapor cannot be
        ! phisically sustained (p = 0, T = 0)
        if ((pSat == 0.0_wp) .and. (TSIn == 0.0_wp)) then

            ! assigning Saturation temperature
            TSat = 0.0_wp

        else

            ! calculating initial estimate for temperature in the TSat procedure. I will also use this variable to
            ! iterate over the Newton's solver
            TSat = TSIn

            ! underrelaxation factor
            Om = 1d-1
            ! starting counter for the Newton solver
            ns = 0

            ! Newton solver for finding the saturation temperature as function of pressure. ns == 0, so the loop is
            ! entered at least once.
            do while ((abs(FT) > ptgalpha_eps) .or. (ns == 0))

                ! Updating counter for the iterative procedure
                ns = ns + 1

                ! calculating residual
                ! FT = A + B / TSat + C * DLOG( TSat ) + D * DLOG( ( pSat + ps_inf( lp ) ) ) - DLOG( pSat + ps_inf( vp ) )
                
                FT = TSat*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                           *(1 - log(TSat)) - (qvps(lp) - qvps(vp)) &
                           + cvs(lp)*(gs_min(lp) - 1)*log(pSat + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)*log(pSat + ps_inf(vp))) &
                     + qvs(lp) - qvs(vp)

                ! calculating the jacobian
                ! dFdT = - B / ( TSat ** 2) + C / TSat

                dFdT = &
                    -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*log(TSat) &
                    - (qvps(lp) - qvps(vp)) &
                    + cvs(lp)*(gs_min(lp) - 1)*log(pSat + ps_inf(lp)) &
                    - cvs(vp)*(gs_min(vp) - 1)*log(pSat + ps_inf(vp))

                ! updating saturation temperature
                TSat = TSat - Om*FT/dFdT

#ifndef MFC_OpenACC
                ! Checking if TSat returns a NaN
                if ((ieee_is_nan(TSat)) .or. (ns > max_iter)) then

                    ! PRINT *, 'FT', FT
                    ! PRINT *, 'TSat, pSat', TSat, pSat

                    call s_int_to_str(ns, nss)
                    call s_real_to_str(TSat, TSatS)
                    call s_real_to_str(pSat, pSatS)
                    call s_mpi_abort('TSat = '//TSatS//', pSat = '// pSatS //' (by assumption of first order transition). &
&                     ns = '//nss//'. m_phase_change, s_TSat. Aborting!')

                end if
#endif
            end do

        end if

        ! PRINT *, 'ps_inf', ps_inf( lp ), ps_inf( vp )

        ! PRINT *, 'ABCD', A, B, C, D
        ! PRINT *, 'FT', FT
        ! PRINT *, 'TSat, pSat', TSat, pSat
    end subroutine s_TSat

    subroutine update_conservative_vars( j, k, l, m0k, pS, q_cons_vf, Tk )
        
        !$acc routine seq
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(wp), intent(IN) :: pS
        real(wp), dimension(num_fluids), intent(IN) :: m0k, Tk
        integer, intent(IN) :: j, k, l
        real(wp), dimension(num_fluids) :: sk, hk, gk, ek, rhok
        integer :: i

        ! Thermodynamic state calculations
        ! entropy
        sk = cvs*DLOG((Tk**gs_min)/((pS + ps_inf)**(gs_min - 1.0_wp))) + qvps

        ! enthalpy
        hk = gs_min*cvs*Tk + qvs

        ! Gibbs-free energy
        gk = hk - Tk*sk

        ! densities
        rhok = (pS + ps_inf)/((gs_min - 1)*cvs*Tk)

        ! internal energy
        ek = (pS + gs_min*ps_inf)/(pS + ps_inf)*cvs*Tk + qvs

        ! calculating volume fractions, internal energies, and total entropy
        ! rhos =  0.0_wp
        !$acc loop seq
        do i = 1, num_fluids

            if ( ( bubbles_euler .eqv. .false. ) .or. ( bubbles_euler .and. (i /= num_fluids) ) ) then

                ! mass factions. Note that, at most, only liquid and vapor masses should change
                q_cons_vf(i + contxb - 1)%sf(j, k, l) = m0k(i)

                ! volume fractions
                q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/rhok(i)

                ! alpha*rho*e
                if (model_eqns == 3) then
                    q_cons_vf(i + intxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)*ek(i)
                end if

                ! Total entropy
                ! rhos = rhos + q_cons_vf(i + contxb - 1)%sf(j, k, l)*sk(i)
            end if

        end do
        
    end subroutine update_conservative_vars

    subroutine s_real_to_str(rl, res)
        character(len=*) :: res
        real(wp), intent(IN) :: rl
        write (res, '(F10.4)') rl
        res = trim(res)
    end subroutine s_real_to_str

#endif
end module m_phase_change