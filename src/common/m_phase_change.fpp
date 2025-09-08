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

    use m_helper_basic         !< Functions to compare floating point numbers

    implicit none

    private; public :: s_initialize_phasechange_module, &
 s_infinite_relaxation_k

    !> @name Parameters for the first order transition phase change
    !> @{
    integer, parameter :: max_iter = 1e8_wp             !< max # of iterations
    real(wp), parameter :: pCr = 4.94e7_wp              !< Critical water pressure
    real(wp), parameter :: TCr = 385.05_wp + 273.15_wp  !< Critical water temperature
    real(wp), parameter :: mixM = sgm_eps               !< threshold for 'mixture cell'. If Y < mixM, phase change does not happen
    integer, parameter :: lp = 1                        !< index for the liquid phase of the reacting fluid
    integer, parameter :: vp = 2                        !< index for the vapor phase of the reacting fluid
    !> @}

    !> @name Gibbs free energy phase change parameters
    !> @{
    real(wp) :: A, B, C, D
    !> @}

    $:GPU_DECLARE(create='[max_iter,pCr,TCr,mixM,lp,vp,A,B,C,D]')

contains

    !>  The purpose of this subroutine is to initialize the phase change module
        !!      by setting the parameters needed for phase change and
        !!      selecting the phase change module that will be used
        !!      (pT- or pTg-equilibrium)
    impure subroutine s_initialize_phasechange_module
        ! variables used in the calculation of the saturation curves for fluids 1 and 2
        A = (gs_min(lp)*cvs(lp) - gs_min(vp)*cvs(vp) + qvps(vp) - qvps(lp)) &
            /((gs_min(vp) - 1.0_wp)*cvs(vp))

        B = (qvs(lp) - qvs(vp))/((gs_min(vp) - 1.0_wp)*cvs(vp))

        C = (gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
            /((gs_min(vp) - 1.0_wp)*cvs(vp))

        D = ((gs_min(lp) - 1.0_wp)*cvs(lp)) &
            /((gs_min(vp) - 1.0_wp)*cvs(vp))

    end subroutine s_initialize_phasechange_module

    !>  This subroutine is created to activate either the pT- (N fluids) or the
        !!      pTg-equilibrium (2 fluids for g-equilibrium)
        !!      model, also considering mass depletion, depending on the incoming
        !!      state conditions.
        !!  @param q_cons_vf Cell-average conservative variables
    impure subroutine s_infinite_relaxation_k(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp) :: pS, pSOV, pSSL !< equilibrium pressure for mixture, overheated vapor, and subcooled liquid
        real(wp) :: TS, TSatOV, TSatSL, TSOV, TSSL !< equilibrium temperature for mixture, overheated vapor, and subcooled liquid. Saturation Temperatures at overheated vapor and subcooled liquid
        real(wp) :: rhoe, rhoe6E, dynE, rhos !< total internal energies (different calculations), kinetic energy, and total entropy
        real(wp) :: rho, rM, m1, m2 !< total density, total reacting mass, individual reacting masses
        logical :: TR

        $:GPU_DECLARE(create='[pS,pSOV,pSSL,TS,TSatOV,TSatSL,TSOV,TSSL]')
        $:GPU_DECLARE(create='[rhoe,rhoe6E,dynE,rhos,rho,rM,m1,m2,TR]')

        real(wp), dimension(num_fluids) :: p_infOV, p_infpT, p_infSL, alphak, alpharhoe0k, m0k, rhok, Tk
        $:GPU_DECLARE(create='[p_infOV,p_infpT,p_infSL,alphak,alpharhoe0k,m0k,rhok,Tk]')

        !< Generic loop iterators
        integer :: i, j, k, l

        ! assigning value to the global parameter
        max_iter_pc_ts = 0

        ! starting equilibrium solver
        $:GPU_PARALLEL_LOOP(collapse=3, private='[pS,pSOV,pSSL,TS,TSatOV,TSatSL,TSOV,TSSL, &
            & rhoe,rhoeT,dynE,rhos,rho,rM,m1,m2,TR, &
            & p_infOV,p_infpT,p_infSL,alphak,alpharhoe0k,m0k,rhok,Tk]')
        do j = 0, m
            do k = 0, n
                do l = 0, p
                    ! trigger for phase change. This will be used for checking many conditions through the code
                    TR = .true.

                    ! mixture internal energy for the 6-equation model. Note that, in sharp interfacesm rhoe \neq rhoe6E
                    rhoe6E = 0.0_wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        ! initial volume fraction
                        alphak(i) = q_cons_vf(i + advxb - 1)%sf(j, k, l)

                        ! initial partial density
                        m0k(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)

                        if ( ieee_is_nan( m0k(i) ) ) then
                            ! print *, 'partial density i = ', i, 'is NaN'
                        end if

                        ! calculating the total internal energy such that the energy-fraction for each of the
                        ! fluids can be proportionally distributed when the sum of the internal energies differs from
                        ! either approach
                        if (model_eqns == 3) then
                            ! initial volume fraction
                            alpharhoe0k(i) = q_cons_vf(i + intxb - 1)%sf(j, k, l)

                            ! total mixture energy by the summation of the internal energy equations
                            rhoe6E = rhoe6E + alpharhoe0k(i)
                        end if
                    end do

                    ! Mixture density
                    rho = sum(m0k)

                    ! calculating the total reacting mass for the phase change process. By hypothesis, this should not change
                    ! throughout the phase-change process. In ANY HYPOTHESIS, UNLESS correcting them artifically
                    rM = m0k(lp) + m0k(vp)

                    ! correcting possible negative mass fraction values
                    ! at this time, TR is updated if phase change needs to be stopped
                    call s_correct_partial_densities(2, alphak, alpharhoe0k, m0k, rM, rho, TR, i, j, k, l)

                    ! kinetic energy as an auxiliary variable to the calculation of the total internal energy
                    dynE = 0.0_wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = momxb, momxe
                        dynE = dynE + 5.0e-1_wp*q_cons_vf(i)%sf(j, k, l)**2 / max(rho, sgm_eps)
                    end do

                    ! calculating the internal mixture energy that MUST be preserved throughout pT- and pTg-relaxation procedures
                    ! This calulation is performed as the total energy minus the kinetic one as energy it is preserved at discontinuities
                    rhoe = q_cons_vf(E_idx)%sf(j, k, l) - dynE

                    ! if phase change is still necessary
                    if (TR) then
                        select case (relax_model)
                        case (1) ! (old) p-equilibrium
                            call s_old_infinite_p_relaxation_k(j, k, l, alphak, alpharhoe0k, m0k, pS, rhoe, Tk)                            
                        case (4) ! p-equilibrium
                            call s_infinite_p_relaxation_k(j, k, l, alphak, alpharhoe0k, m0k, pS, rho, rhoe, rhoe6E, Tk)
                        case (5) ! pT-equilibrium
                            ! for this case, MFL cannot be either 0 or 1, so I chose it to be 2
                            call s_infinite_pt_relaxation_k(j, k, l, m0k, 2, pS, p_infpT, rhoe, rM, TS)
                            Tk = spread(TS, 1, num_fluids)
                        case (6) ! pT-pTg equilibrium
                            ! pT-equilibrium as rhe initial condition
                            call s_infinite_pt_relaxation_k(j, k, l, m0k, 2, pS, p_infpT, rhoe, rM, TS)
                            Tk = spread(TS, 1, num_fluids)

                            ! updating the densities and volume fractions used for thresholds
                            rhok = (pS + ps_inf)/((gs_min - 1)*cvs*Tk)
            
                            ! new volume fractions, after partial densities and p- or pT-equilibrium
                            alphak = m0k / rhok

                            ! 1 - model activation, 1st order transition (p,T) <= (pCr, TCr)
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
                                m1 = m0k(lp) ; m2 = m0k(vp)

                                ! checking if fluid is either subcoooled liquid or overheated vapor (NOT metastability)

                                ! overheated vapor
                                ! depleting the mass of liquid and tranferring the total mass to vapor
                                m0k(lp) = mixM*rM ; m0k(vp) = (1.0_wp - mixM)*rM

                                ! calling pT-equilibrium for overheated vapor, which is MFL = 0
                                call s_infinite_pt_relaxation_k(j, k, l, m0k, 0, pSOV, p_infOV, rhoe, rM, TSOV)

                                ! calculating Saturation temperature
                                call s_TSat(pSOV, TSatOV, TSOV)

                                ! subcooled liquid 
                                ! tranferring the total mass to liquid and depleting the mass of vapor
                                m0k(lp) = (1.0_wp - mixM)*rM ; m0k(vp) = mixM*rM

                                ! calling pT-equilibrium for subcooled liquid, which is MFL = 1                       
                                call s_infinite_pt_relaxation_k(j, k, l, m0k, 1, pSSL, p_infSL, rhoe, rM, TSSL)

                                ! calculating Saturation temperature
                                call s_TSat(pSSL, TSatSL, TSSL)

                                ! checking the conditions for overheated vapor
                                if (TSOV > TSatOV) then

                                    ! Assigning pressure and temperature
                                    pS = pSOV ; TS = TSOV

                                    ! correcting the liquid and vapor partial densities
                                    m0k(lp) = mixM*rM ; m0k(vp) = (1.0_wp - mixM)*rM

                                ! checking the conditions for subcooled liquid
                                elseif (TSSL < TSatSL) then

                                    ! Assigning pressure and temperature
                                    pS = pSSL ; TS = TSSL

                                    ! correcting the liquid and vapor partial densities
                                    m0k(lp) = (1.0_wp - mixM)*rM ; m0k(vp) = mixM*rM

                                ! if not, mixture of fluids. Starting phase change (pTg)
                                else
                                    ! returning partial pressures to what they were after the partial density correction 
                                    m0k(lp) = m1 ; m0k(vp) = m2

                                    ! pTg-relaxation
                                    call s_infinite_ptg_relaxation_k(j, k, l, alphak, alpharhoe0k, m0k, pS, p_infpT, rho, rhoe, rM, TR, TS)
                                    ! if no pTg happens, the solver will return to the hyperbolic state variables
                                    if ( TR .eqv. .false. ) then
                                        $:GPU_LOOP(parallelism='[seq]')
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
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_fluids
                                    ! returning partial densities to what they were previous to any relaxation scheme.
                                    m0k(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)
                                end do
                                ! cycles the innermost loop to the next iteration
                                cycle 
                            end if
                        end select
                    else
                        $:GPU_LOOP(parallelism='[seq]')
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
    impure subroutine s_infinite_p_relaxation_k(j, k, l, alphaik, alpharhoeik, mik, pS, rho, rhoe, rhoe6E, Tk)
        !$acc routine seq

        ! initializing variables
        real(wp), intent(in) :: rho, rhoe, rhoe6E
        real(wp), intent(out) :: pS
        real(wp), dimension(num_fluids), intent(out) :: Tk
        real(wp), dimension(num_fluids), intent(in)  :: alphaik, alpharhoeik, mik
        integer, intent(in) :: j, k, l

        real(wp) :: fp, fpp, pO !< variables for the Newton Solver
        real(wp) :: Econst, gamma, pi_inf, mQ, TS !< auxiliary variables
        real(wp), dimension(num_fluids) :: alpha0k, alpharhoe0k, m0k
        real(wp), dimension(num_fluids) :: alphak, alpharhoek, rhok
        character(20) :: nss, pSs, Econsts
        integer, dimension(num_fluids) :: iFix, iVar !< auxiliary index for choosing appropiate values for conditional sums

        integer :: i, na, ns, nsL !< generic loop iterators

        iFix = (/ (i, i=1,num_fluids) /) 

        ! think how to solve pK, alphak, alpharhoek, given the energy mismatechs at the discontinuities 
        alpha0k = alphaik

        ! Re-distributing the initial internal energy values such that rhoe = rhoe6E, eventually.
        ! alpharhoe0k does the paper of alpharhoeik, but it is its value adjusted to match rhoe = rhoe6E
        alpharhoe0k = alpharhoeik ! * rhoe / rhoe6E

        ! distributing the partial density
        m0k = mik

        ! Numerical correction of the volume fractions, partial densities, and internal energies
        IF (mpp_lim) THEN          
          ! auxiliry variable to avoid do loops
          iVar = iFix
          iVar( pack( iFix, ( alphaik >= 0 ) .and. ( mik >= 0 ) ) ) = 0
          ! if alpha < 0 or m < 0, set everything to 0
          m0k( pack( iVar, iVar /= 0 ) ) = 0.0_wp
          alpharhoe0k( pack( iVar, iVar /= 0 ) ) = 0.0_wp
          alpha0k( pack( iVar, iVar /= 0 ) ) = 0.0_wp
          ! if alpha > 1, set it to 1
          alpha0k( pack( iFix, alphaik > 1.0_wp ) ) = 1.0_wp

          ! renormalize alpha
          alpha0k = alpha0k / sum(alpha0k)
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

                        ! call s_whistleblower((/ 0.0_wp,  0.0_wp/), reshape((/ 0.0_wp,  0.0_wp,  0.0_wp,  0.0_wp/), (/2, 2/)) &
                        !                 , j, (/ 0.0_wp,  0.0_wp,  0.0_wp,  0.0_wp/), k, l, rhoe6E, ps_inf, pS, (/pS - pO, pS + pO/) &
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

    impure subroutine s_old_infinite_p_relaxation_k(j, k, l, alpha0k, alpharhoe0k, m0k, pS, rhoe, Tk)
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
        real(wp), intent(in) :: rhoe
        real(wp), intent(out) :: pS
        real(wp), dimension(num_fluids), intent(out) :: Tk
        real(wp), dimension(num_fluids), intent(in) :: alpha0k, alpharhoe0k, m0k

        integer, intent(in) :: j, k, l

        real(wp) :: fp, fpp, mQ, pO !< variables for the Newton Solver
        real(wp) :: gamma, pi_inf !< auxiliary variables
        real(wp), dimension(num_fluids) :: alphak, alpharhoek, mk, pk, rhok, num, den, drhodp
        integer, dimension(num_fluids) :: iVar, iFix !< auxiliary index for choosing appropiate values for conditional sums
        character(20) :: nss, pSs
        !> @}

        integer :: i, ns !< generic loop iterators
        
        iFix = (/ (i, i=1,num_fluids) /)

        ! initializing the partial energies
        alphak       = alpha0k
        alpharhoek   = alpharhoe0k
        mk           = m0k

        ! Numerical correction of the volume fractions, partial densities, and internal energies
        IF (mpp_lim) THEN  
          ! auxiliry variable to avoid do loops
          iVar = iFix
          iVar( pack( iFix, ( alpha0k >= 0 ) .and. ( m0k >= 0 ) ) ) = 0
          
          ! if alpha < 0 or m < 0, set everything to 0
          mk( pack( iVar, iVar /= 0 ) ) = 0.0_wp
          alpharhoek( pack( iVar, iVar /= 0 ) ) = 0.0_wp
          alphak( pack( iVar, iVar /= 0 ) ) = 0.0_wp
          ! if alpha > 1, set it to 1
          alphak( pack( iFix, alpha0k > 1.0_wp ) ) = 1.0_wp

          ! renormalize alpha
          alphak = alphak / sum(alphak)           
        END IF

        ! Initial state
        pk = 0.0_wp
        
        ! auxiliry variable to avoid do loops
        iVar = iFix
        iVar( pack( iFix, alphak <= sgm_eps ) ) = 0

        pk( pack( iVar, iVar /= 0 ) ) = ( ( alpharhoek( pack( iVar, iVar /= 0 ) ) &
        - mk( pack( iVar, iVar /= 0 ) ) * qvs( pack( iVar, iVar /= 0 ) ) ) /      &
        alphak( pack( iVar, iVar /= 0 ) ) - fluid_pp( pack( iVar, iVar /= 0 ) )%pi_inf) &
        / fluid_pp( pack( iVar, iVar /= 0 ) )%gamma               

        ! auxiliry variable to avoid do loops
        iVar = iFix
        iVar( pack( iFix, .not. ( pk < -(1.0_wp - ptgalpha_eps)*ps_inf + ptgalpha_eps ) ) ) = 0

        pk( pack( iVar, iVar /= 0 ) ) = -(1.0_wp - ptgalpha_eps)*ps_inf( pack( iVar, iVar /= 0 ) ) &
        + ptgalpha_eps
        
        ! initial guess for the relaxed pressure
        pS = sum( alphak * pk )

        ! Iterative process for relaxed pressure determination
        ns    = 0
        fp    = 2.0_wp * ptgalpha_eps
        fpp   = 2.0_wp * ptgalpha_eps

        ! Note that a solution is found when f(p) = 1
        ! keep an eye on this
        DO WHILE (abs(fp) > ptgalpha_eps .or. (ns == 0))
            
            ! setting old pressure
            pO = pS

            ! Newton-Raphson method
            ! Physical pressure?
            if ( pO < minval(-(1.0_wp - ptgalpha_eps)*ps_inf + ptgalpha_eps) ) then
              pO = minval(-(1.0_wp - ptgalpha_eps)*ps_inf + ptgalpha_eps)
            end if

            ! calulating fp and fpp
            iVar = iFix
            iVar( pack( iFix, alphak <= sgm_eps ) ) = 0
            
            num       = gs_min * (pO + ps_inf)
            den       = num + pk - pO
            
            rhok = 0.0_wp  
            rhok( pack( iVar, iVar /= 0 ) ) = mk( pack( iVar, iVar /= 0 ) ) &
            / alphak( pack( iVar, iVar /= 0 ) ) * num( pack( iVar, iVar /= 0 ) ) &
            / den( pack( iVar, iVar /= 0 ) )

            drhodp = 0.0_wp
            drhodp( pack( iVar, iVar /= 0 ) ) = mk( pack( iVar, iVar /= 0 ) ) &
            / alphak( pack( iVar, iVar /= 0 ) ) * gs_min( pack( iVar, iVar /= 0 ) ) &
            * ( pk( pack( iVar, iVar /= 0 ) ) + ps_inf( pack( iVar, iVar /= 0 ) ) ) &
            / ( den( pack( iVar, iVar /= 0 ) ) ** 2 )

            fp        = sum( mk( pack( iVar, iVar /= 0 ) ) / rhok( pack( iVar, iVar /= 0 ) ) ) - 1.0_wp
            fpp       = sum( - mk( pack( iVar, iVar /= 0 ) ) * drhodp( pack( iVar, iVar /= 0 ) ) &
            / ( rhok( pack( iVar, iVar /= 0 ) ) **2 ) )
            
            ! updating pressure
            pS = pO - fp / fpp

            ! Convergence?
            ns = ns + 1
#ifndef MFC_OpenACC
                ! energy constraint for the p-equilibrium
                ! checking if pressure is within expected bounds
                if ((pS <= -1.0_wp*minval(gs_min*ps_inf)) .or. (ieee_is_nan(pS)) .or. (ns > max_iter)) then
                    if (proc_rank == 0) then
                        ! call s_whistleblower((/ 0.0_wp,  0.0_wp/), reshape((/ 0.0_wp,  0.0_wp,  0.0_wp,  0.0_wp/), (/2, 2/)) &
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
        alphak( pack( iFix, alphak > sgm_eps ) ) = &
                mk( pack( iFix, alphak > sgm_eps ) ) / rhok( pack( iFix, alphak > sgm_eps ) )

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
    impure subroutine s_infinite_pt_relaxation_k(j, k, l, m0k, MFL, pS, p_infpT, rhoe, rM, TS)

        $:GPU_ROUTINE(function_name='s_infinite_pt_relaxation_k', &
            & parallelism='[seq]', cray_inline=True)

        ! initializing variables
        real(wp), intent(out) :: pS, TS
        real(wp), dimension(num_fluids), intent(out) :: p_infpT
        real(wp), intent(in) :: rhoe, rM
        real(wp), intent(in), dimension(num_fluids) :: m0k
        integer, intent(in) :: j, k, l, MFL
        integer, dimension(num_fluids) :: iVar, iFix !< auxiliary index for choosing appropiate values for conditional sums
        real(wp) :: gp, gpp, hp, pO, mCP, mQ !< variables for the Newton Solver
        character(20) :: nss, pSs, Econsts

        integer :: i, ns !< generic loop iterators

        ! auxiliary variables for the pT-equilibrium solver
        p_infpT = ps_inf

        ! auxiliary variables to skip index looping
        iFix = (/ (i, i=1,num_fluids) /)
        iVar = iFix

        ! dismissing fluids that do not participate into pT-relaxation due to the small amount of mass fraction
        iVar( pack( iFix, .not. ( ( m0k - rM * mixM <= sgm_eps ) .and. ( m0k >= 0 ) ) ) ) = 0
        
        ! this value is rather arbitrary, as I am interested in MINVAL( ps_inf ) for the solver.
        ! This way, I am ensuring this value will not be selected.
        p_infpT( pack(iVar, iVar /= 0) ) = 2 * MAXVAL( ps_inf )
        
        ! if ( ( bubbles_euler .eqv. .false. ) .or. ( bubbles_euler .and. (i /= num_fluids) ) ) then
          ! sum of the total alpha*rho*cp of the system
          mCP = sum( m0k * cvs * gs_min )

          ! sum of the total alpha*rho*q of the system
          mQ = sum( m0k * qvs )
        ! end if

        ! Checking energy constraint. In the case we are calculating the possibility of having subcooled liquid or
        ! overheated vapor, the energy constraint might not be satisfied, as are hypothetically transferring all the 
        ! mass from one phase to the other. When this is the case, we simply ignore this possibility, set pS = TS = 0,
        ! and discard the hypothesis. The solver can thus move forward.
        if ((rhoe - mQ - minval(p_infpT)) < 0.0_wp) then

            if ( any((/ 0, 1 /) == MFL ) ) then

                ! Assigning zero values for mass depletion cases
                ! pressure and temperature
                pS = 0.0_wp ; TS = 0.0_wp

                return
#ifndef MFC_OpenACC
            else
                if (proc_rank == 0) then

                    ! call s_whistleblower((/ 0.0_wp,  0.0_wp/), reshape((/ 0.0_wp,  0.0_wp,  0.0_wp,  0.0_wp/), (/2, 2/)) &
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
        do while ( ( ( abs( gp - 1.0_wp ) > ptgalpha_eps ) ) .or. ( ns == 0 ) )

            ! increasing counter
            ns = ns + 1

            ! updating old pressure
            pO = pS

            ! updating functions used in the Newton's solver
            gp = sum( ( gs_min - 1.0_wp ) * m0k * cvs * ( rhoe + pO - mQ ) &
              / ( mCP * ( pO + p_infpT ) ) ) &
               - sum( ( gs_min( pack(iVar, iVar /= 0) ) - 1.0_wp ) * m0k( pack(iVar, iVar /= 0) ) &
               * cvs( pack(iVar, iVar /= 0) ) * ( rhoe + pO - mQ ) &
              / ( mCP * ( pO + p_infpT( pack(iVar, iVar /= 0) ) ) ) )

            gpp = sum( ( gs_min - 1.0_wp ) * m0k * cvs * ( p_infpT - rhoe + mQ ) &
              / ( mCP * ( pO + p_infpT ) **2 ) ) &
              - sum( ( gs_min( pack(iVar, iVar /= 0) ) - 1.0_wp ) * m0k( pack(iVar, iVar /= 0) ) &
              * cvs( pack(iVar, iVar /= 0) ) * ( p_infpT( pack(iVar, iVar /= 0) ) - rhoe + mQ )  &
              / ( mCP * ( pO + p_infpT( pack(iVar, iVar /= 0) ) ) **2 ) )

            hp = 1.0_wp/(rhoe + pO - mQ) + 1.0_wp/(pO + minval(p_infpT))

            ! updating common pressure for the newton solver
            pS = pO + ((1.0_wp - gp)/gpp)/(1.0_wp - (1.0_wp - gp + abs(1.0_wp - gp)) &
                                          /(2.0_wp*gpp)*hp)

            ! check if solution is out of bounds (which I believe it won`t happen given the solver is gloabally convergent.
#ifndef MFC_OpenACC
            if ((pS <= -1.0_wp*minval(p_infpT)) .or. (ieee_is_nan(pS)) .or. (ns > max_iter)) then
                if (proc_rank == 0) then
                    ! call s_whistleblower((/0.0_wp, 0.0_wp/), reshape((/0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp/), (/2, 2/)) &
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
    impure subroutine s_infinite_ptg_relaxation_k(j, k, l, alphak, alpharhoe0k, m0k, pS, p_infpT, rho, rhoe, rM, TR, TS)

        $:GPU_ROUTINE(function_name='s_infinite_ptg_relaxation_k', &
            & parallelism='[seq]', cray_inline=True)

        real(wp), dimension(num_fluids), intent(inout) :: alphak, alpharhoe0k, m0k
        real(wp), dimension(num_fluids), intent(in) :: p_infpT
        real(wp), intent(inout) :: pS, TS, rM
        real(wp), intent(in) :: rho, rhoe
        integer, intent(in) :: j, k, l
        logical, intent(inout) :: TR
        real(wp), dimension(num_fluids) :: p_infpTg
        real(wp), dimension(2, 2) :: Jac, InvJac, TJac
        real(wp), dimension(2) :: R2D, DeltamP
        real(wp), dimension(3) :: Oc
        real(wp) :: Om, OmI ! underrelaxation factor
        real(wp) :: mCP, mCPD, mCVGP, mCVGP2, mQ, mQD, TSat ! auxiliary variables for the pTg-solver
        character(20) :: nss, pSs, Econsts, R2D1s, R2D2s 

        !< Generic loop iterators
        integer :: i, ns

        ! assigning the relexant pi_infs based on the previous pT-equilibrium
        p_infpTg = p_infpT

        ! checking if homogeneous cavitation is expected. If yes, transfering a small amount of mass to the depleted
        ! phase, and then let the algorithm run.
        
        ! is the fluid at a metastable state with enough 'energy' for phase change to happen?
        if ((pS < 0.0_wp) .and. (rM > (rhoe - gs_min(lp)*ps_inf(lp)/(gs_min(lp) - 1.0e-1_wp))/qvs(lp))) then

            ! transfer a bit of mass to the deficient phase, enforce phase change
            call s_correct_partial_densities(1, alphak, alpharhoe0k, m0k, rM, rho, TR, i, j, k, l)

            ! The following changes do not have a strong physicial rationale. So keep an eye on them,
            ! as they might be a weakness of the solver

            ! reestabilish the fluid parameters that are now important for the phase change
            p_infpTg(lp) = ps_inf(lp) ; p_infpTg(vp) = ps_inf(vp)

            ! give an arbitrary value 'positive' value for the pressure as,
            ! since now both vapor and liquid exist pS > -min(pi_inf) for the 
            ! solver to converge (at least) within the pR-relaxation context
            pS = 1.0e4_wp

        ! the metastable state is not enough to sustain phase change
        elseif (pS < 0.0_wp) then
            
            ! cancel any phase-change updates.
            TR = .false.

            ! ends the execution of this function and returns control to the calling function
            return

        ! if not homogeneous, then heterogeneous. Thus, setting up an arbitrary initial condition in case the one from
        ! the p(T)-equilibrium solver could lead to numerical issues
        elseif ((pS < 1.e-1_wp) .and. (pS >= 0.0_wp)) then
            ! improve this initial condition
            pS = 1.0e4_wp
        end if

        ! Relaxation factor. This value is rather arbitrary, with a certain level of self adjustment.
        OmI = 1.0e-1_wp
        ! Critical relaxation factors, for variable sub-relaxation
        Oc = OmI;

        R2D = 0.0_wp ; DeltamP = 0.0_wp;
        ! starting counter for the Newton solver
        ns = 0

        ! Newton solver for pTg-equilibrium. 1d6 is arbitrary, and ns == 0, to the loop is entered at least once.
        do while (((sqrt(R2D(1)**2 + R2D(2)**2) > ptgalpha_eps) &
                    .and. ((sqrt(R2D(1)**2 + R2D(2)**2)/rhoe) > (ptgalpha_eps/1.e6_wp))) &
                  .or. (ns == 0))

            ! Updating counter for the iterative procedure
            ns = ns + 1

            ! Auxiliary variables to help in the calculation of the residue
            ! Those must be updated through the iterations, as they either depend on
            ! the partial masses for all fluids, or on the equilibrium pressure

            ! sum of the total alpha*rho*cp of the system
            mCP = sum( m0k * cvs * gs_min )

            ! mCP - the contribution from the reacting phases
            mCPD = mCP - m0k(lp) * cvs(lp) * gs_min(lp) &
                       - m0k(vp) * cvs(vp) * gs_min(vp)

            ! sum of the total alpha*rho*q of the system
            mQ = sum( m0k * qvs )

            ! mQ - the contribution from the reacting phases
            mQD = mQ - m0k(lp) * qvs(lp) &
                     - m0k(vp) * qvs(vp)

            ! mCVG - the contribution from the reacting phases
            mCVGP = sum( m0k * cvs * ( gs_min - 1 ) / ( pS + ps_inf ) ) &
                  - m0k(lp) * cvs(lp) * ( gs_min(lp) - 1 ) / ( pS + ps_inf(lp) ) &
                  - m0k(vp) * cvs(vp) * ( gs_min(vp) - 1 ) / ( pS + ps_inf(vp) )
              
            ! mCVG2 - the contribution from the reacting phases
            mCVGP2 = sum( m0k * cvs * ( gs_min - 1 ) / ( ( pS + ps_inf ) ** 2 ) ) &
                  - m0k(lp) * cvs(lp) * ( gs_min(lp) - 1 ) / ( ( pS + ps_inf(lp) ) ** 2 ) &
                  - m0k(vp) * cvs(vp) * ( gs_min(vp) - 1 ) / ( ( pS + ps_inf(vp) ) ** 2 )

            ! calculating the (2D) Jacobian Matrix used in the solution of the pTg-quilibrium model
            call s_compute_jacobian_matrix(InvJac, j, Jac, k, l, m0k, mCPD, mCVGP, mCVGP2, pS, rM, TJac)

            ! print *, InvJac, Jac, TJac

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
            if (pS + minval(p_infpTg) - Om*DeltamP(2) <= 0.0_wp) then
                Oc(3) = (pS + minval(p_infpTg))/(2*DeltamP(2))
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

                print *, Om
                call s_int_to_str(ns, nss)
                call s_real_to_str(R2D(1), R2D1s)
                call s_real_to_str(R2D(2), R2D2s)
                call s_mpi_abort('Residual for the pTg-relaxation possibly returned NaN values. ns = ' &
                                  //nss//', R2D(1) = '//R2D1s//', R2D(2) = '//R2D2s// &
                                  ' (m_phase_change, s_infinite_ptg_relaxation_k). Aborting!')

            end if
            ! if not,
            ! Checking pressure and energy criteria for the (pT) solver to find a solution
            if ((pS <= -1.0_wp*minval(p_infpTg)) .or. ((rhoe - mQ - minval(p_infpTg)) < 0.0_wp)) then
                if (proc_rank == 0) then
                    ! call s_whistleblower(DeltamP, InvJac, j, Jac, k, l, mQ, p_infpTg, pS &
                    !                   , R2D, rhoe, q_cons_vf, TS)
                end if

                call s_real_to_str(rhoe - mQ - minval(p_infpTg), Econsts)
                call s_real_to_str(pS, pSs)
                call s_mpi_abort('Solver for the pTg-relaxation failed (m_phase_change, s_infinite_ptg_relaxation_k). &
                &   pS ~'//pSs//'. Econst = '//Econsts//'. Aborting!')

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
    impure subroutine s_correct_partial_densities(CT, alpha0k, alpharhoe0k, m0k, rM, rho, TR, i, j, k, l)
        $:GPU_ROUTINE(function_name='s_correct_partial_densities', &
            & parallelism='[seq]', cray_inline=True)

        !> @name variables for the correction of the reacting partial densities
        !> @{
        real(wp), dimension(num_fluids), intent(inout) :: alpha0k, alpharhoe0k, m0k
        real(wp), intent(inout) :: rM
        real(wp), intent(in) :: rho
        logical, intent(inout) :: TR
        integer, intent(in) :: CT, j, k, l
        integer, dimension(num_fluids) :: iFix, iVar !< auxiliary index for choosing appropiate values for conditional sums
        integer :: i
        !> @}

        iFix = (/ (i, i=1,num_fluids) /)
        iVar = iFix

        ! CT = 0: No Mass correction; ! CT = 1: Reacting Mass correction;
        ! CT = 2: crude correction; else: Total Mass correction
        if (CT == 0) then
          if ( sum( pack(m0k, m0k < rho * mixM ) ) /= 0 ) then
            TR = .false.
          end if
        elseif (CT == 1) then
            if (rM < 0.0_wp) then
                ! reacting masses are very negative so as to affect the physics of the problem, so phase change will not be activated
                if ( any( (/ m0k(lp), m0k(vp) /) < rM * mixM ) ) then
                    ! do not continue relaxation
                    TR = .false.
                    ! reacting masses are not as negative so I can disregard them,
                    ! expecting no significant changes in the physics of the simulation
                else
                    m0k(lp) = mixM*rM ; m0k(vp) = mixM*rM

                    ! continue relaxation
                    TR = .true.
                end if
            ! correcting the partial densities of the reacting fluids. In case liquid is negative
            elseif (m0k(lp) < rM * mixM) then

                m0k(lp) = mixM*rM ; m0k(vp) = (1.0_wp - mixM)*rM
                
                ! continue relaxation
                TR = .true.

            ! correcting the partial densities of the reacting fluids. In case vapor is negative
            elseif (m0k(vp) < rM * mixM) then

                m0k(lp) = (1.0_wp - mixM)*rM ; m0k(vp) = mixM*rM

                ! continue relaxation
                TR = .true.

            end if
        elseif (CT == 2) then
            ! auxiliry variable to avoid do loops - use Morgan's Law
            iVar( pack( iFix, ( alpha0k >= 0 ) .and. ( m0k >= 0 ) ) ) = 0

            ! if either the volume fraction or the partial density is negative, make them positive
            alpha0k(pack(iVar, iVar /= 0) ) = 0.0_wp
            m0k(pack(iVar, iVar /= 0) ) = 0.0_wp
            if (model_eqns .eq. 3) then
                alpharhoe0k(pack(iVar, iVar /= 0) ) = 0.0_wp
            end if

            ! continue relaxation
            TR = .true.
        else
            ! if there are any insignificant values, make them significant 
            m0k( pack( iFix, m0k < rho * mixM ) ) = rho * mixM

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
    impure subroutine s_compute_jacobian_matrix(InvJac, j, Jac, k, l, m0k, mCPD, mCVGP, mCVGP2, pS, rM, TJac)
        $:GPU_ROUTINE(function_name='s_compute_jacobian_matrix', &
            & parallelism='[seq]', cray_inline=True)

        real(wp), dimension(num_fluids), intent(in) :: m0k
        real(wp), intent(in) :: pS, mCPD, mCVGP, mCVGP2, rM
        integer, intent(in) :: j, k, l
        real(wp), dimension(2, 2), intent(out) :: Jac, InvJac, TJac
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
    impure subroutine s_compute_pTg_residue(j, k, l, m0k, mCPD, mCVGP, mQD, pS, rhoe, rM, R2D)
        $:GPU_ROUTINE(function_name='s_compute_pTg_residue', &
            & parallelism='[seq]', cray_inline=True)

        real(wp), dimension(num_fluids), intent(in) :: m0k
        real(wp), intent(in) :: pS, rhoe, mCPD, mCVGP, mQD, rM
        integer, intent(in) :: j, k, l
        real(wp), dimension(2), intent(out) :: R2D
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
    impure subroutine s_whistleblower(DeltamP, InvJac, j, Jac, k, l, mQ, p_infA, pS, R2D, rhoe, q_cons_vf, TS) ! ----------------

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        real(wp), dimension(2, 2), intent(in) :: Jac, InvJac
        real(wp), dimension(num_fluids), intent(in) :: p_infA
        real(wp), dimension(2), intent(in) :: R2D, DeltamP
        real(wp), intent(in) :: pS, TS
        real(wp), intent(in) :: rhoe, mQ
        integer, intent(in) :: j, k, l
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

    end subroutine s_whistleblower

        !>  This auxiliary subroutine finds the Saturation temperature for a given
        !!      saturation pressure through a newton solver
        !!  @param pSat Saturation Pressure
        !!  @param TSat Saturation Temperature
        !!  @param TSIn equilibrium Temperature
    impure subroutine s_TSat(pSat, TSat, TSIn)
        $:GPU_ROUTINE(function_name='s_TSat',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), intent(out) :: TSat
        real(wp), intent(in) :: pSat, TSIn
        real(wp) :: dFdT, FT, Om !< auxiliary variables
        character(20) :: nss, pSatS, TSatS

        ! Generic loop iterators
        integer :: ns

        ! in case of fluid under tension (p - p_inf > 0, T > 0), or, when subcooled liquid/overheated vapor cannot be
        ! phisically sustained (p = 0, T = 0)
        if ((pSat <= 0.0_wp) .and. (TSIn >= 0.0_wp)) then

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

                    call s_int_to_str(ns, nss)
                    call s_real_to_str(TSat, TSatS)
                    call s_real_to_str(pSat, pSatS)
                    call s_mpi_abort('TSat = '//TSatS//', pSat = '// pSatS //' (by assumption of first order transition). &
&                     ns = '//nss//'. m_phase_change, s_TSat. Aborting!')

                end if
#endif
            end do

        end if

    end subroutine s_TSat

    impure subroutine update_conservative_vars( j, k, l, m0k, pS, q_cons_vf, Tk )
        
        !$acc routine seq
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(in) :: pS
        real(wp), dimension(num_fluids), intent(in) :: m0k, Tk
        integer, intent(in) :: j, k, l
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
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids
            ! if ( ( bubbles_euler .eqv. .false. ) .or. ( bubbles_euler .and. (i /= num_fluids) ) ) then
          
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
            ! end if
        end do
        
    end subroutine update_conservative_vars

    impure subroutine s_real_to_str(rl, res)
        real(wp), intent(in) :: rl
        character(len=*), intent(out) :: res
        write (res, '(F10.4)') rl
        res = trim(res)
    end subroutine s_real_to_str

#endif
end module m_phase_change