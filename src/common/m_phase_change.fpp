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

    private; public :: pVap_sf, &
                       s_initialize_phasechange_module, &
                       s_finalize_phasechange_module, &
                       s_compute_bubbles_euler_vapor_pressure, &
                       s_infinite_relaxation_k

    !> @name Parameters for the first order transition phase change
    !> @{
    integer, parameter  :: max_iter = 50                       !< max # of iterations
    real(wp), parameter :: pCr      = 4.94e7_wp                !< Critical water pressure
    real(wp), parameter :: TCr      = 385.05_wp + 273.15_wp    !< Critical water temperature
    real(wp), parameter :: mixM     = 0*sgm_eps                !< threshold for 'mixture cell'. If Y < mixM, phase change does not happen
    integer, parameter  :: lp       = 1                        !< index for the liquid phase of the reacting fluid
    integer, parameter  :: vp       = 2                        !< index for the vapor phase of the reacting fluid
    !> @}

    real(wp), allocatable, dimension(:, :, :) :: pVap_sf !< Relaxed cellwise vapor pressure used by Eulerian bubbles
    $:GPU_DECLARE(create='[pVap_sf]')

contains

    impure subroutine s_initialize_phasechange_module()

        if (.not. bubbles_euler) return

        @:ALLOCATE(pVap_sf(0:m, 0:n, 0:p))
        pVap_sf = pv
        $:GPU_UPDATE(device='[pVap_sf]')

    end subroutine s_initialize_phasechange_module

    impure subroutine s_finalize_phasechange_module()

        if (allocated(pVap_sf)) then
            @:DEALLOCATE(pVap_sf)
        end if

    end subroutine s_finalize_phasechange_module

    !>  Compute the common initial vapor pressure used by Eulerian bubbles from
        !!      the reference bubble temperature. At initialization, all Eulerian
        !!      bubbles share the same temperature, so they also share the same
        !!      saturation pressure. This assumes fluid 1 is the host liquid and
        !!      fluid 2 is its vapor, matching the saturation-property solver.
    impure subroutine s_compute_bubbles_euler_vapor_pressure()

        real(wp) :: pGuess, pVap, TSatRef
        
        ! A dflt_real input signals that no vapor pressure should be used. Additionally
        ! if the Eulerian bubble model is not activated, there is no need to compute the vapor pressure for the bubbles.
        if ( ( bub_pp%pv == dflt_real .or. bub_pp%pv == 0.0_wp ) .or. (.not. bubbles_euler) .or. (num_fluids < 2) ) return

        ! if this is not the case, then we update the vapor pressure pv to its saturation value at the reference temperature. This is the initial vapor pressure for all the bubbles, and it will be updated during the simulation with the pT- or pTg-equilibrium solver.
        call s_Saturation_Properties(pVap, TSatRef, bub_pp%p0ref, 2)
        pv = pVap
        bub_pp%pv = pVap

        if (allocated(pVap_sf)) then
            pVap_sf = pVap
            $:GPU_UPDATE(device='[pVap_sf]')
        end if

    end subroutine s_compute_bubbles_euler_vapor_pressure

    logical elemental function f_is_negligible_phase_mass(mass, reacting_mass) result(is_negligible)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: mass, reacting_mass

        ! check if this should be a strict inequality or not. In some cases, for whatever reason, the nonstrict inequality did not work.
        is_negligible = (0.0_wp <= mass) .and. (mass <= reacting_mass * mixM + sgm_eps)
    end function f_is_negligible_phase_mass

    !>  This subroutine is created to activate either the pT- (N fluids) or the
        !!      pTg-equilibrium (2 fluids for g-equilibrium)
        !!      model, also considering mass depletion, depending on the incoming
        !!      state conditions.
        !!  @param q_cons_vf Cell-average conservative variables
    subroutine s_infinite_relaxation_k(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp) :: pS, pSOV, pSSL !< equilibrium pressure for mixture, overheated vapor, and subcooled liquid
        real(wp) :: TS, TSatOV, TSatSL, TSOV, TSSL !< equilibrium temperature for mixture, overheated vapor, and subcooled liquid. Saturation Temperatures at overheated vapor and subcooled liquid
        real(wp) :: pVapSG !< Vapor pressure from the intermediate pT state, used only for the Blake/subgrid trigger
        real(wp) :: rhoe, dynE !< total internal energies (different calculations), kinetic energy, and total entropy
        real(wp) :: rho, rM !< total density, total reacting mass
        real(wp) :: alpha_b !< volume fraction of the subgrid component, used in case icsg is activated
        logical :: TR, TIC, TSG

        $:GPU_DECLARE(create='[pS,pSOV,pSSL,TS,TSatOV,TSatSL,TSOV,TSSL,pVapSG]')
        $:GPU_DECLARE(create='[rhoe,dynE,rho,rM,TR]')

        real(wp), dimension(nb) :: mass_b, R_b !< subgrid variables, used in case icsg is activated

        real(wp), dimension(num_fluids) :: p_infOV, p_infpT, p_infSL, alphak, me0k, m0k, rhok, Tk
        $:GPU_DECLARE(create='[p_infOV,p_infpT,p_infSL,alphak,me0k,m0k,rhok,Tk]')

        real(wp), dimension(2) :: mOr !< Individual reacting masses
        $:GPU_DECLARE(create='[mOr]')

        !< Generic loop iterators
        integer :: cb, i, j, k, l

        ! assigning value to the global parameter
        max_iter_pc_ts = 0

        ! starting equilibrium solver
        $:GPU_PARALLEL_LOOP(collapse=3, private='[pS,pSOV,pSSL,TS,TSatOV,TSatSL,TSOV,TSSL,pVapSG,rhoe,rhoeT,dynE,rho,rM,TR,p_infOV,p_infpT,p_infSL,alphak,me0k,m0k,mOr,rhok,Tk]')
        do j = 0, m
            do k = 0, n
                do l = 0, p
                    ! trigger for phase change. This will be used for checking many conditions through the code
                    TR = .true.

                    if (bubbles_euler) pVap_sf(j, k, l) = pv

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        ! initial volume fraction
                        alphak(i) = q_cons_vf(i + advxb - 1)%sf(j, k, l)

                        ! initial partial density
                        m0k(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)

                        ! calculating the total internal energy such that the energy-fraction for each of the
                        ! fluids can be proportionally distributed when the sum of the internal energies differs from
                        ! either approach
                        if (model_eqns == 3) then
                            ! initial volume fraction
                            me0k(i) = q_cons_vf(i + intxb - 1)%sf(j, k, l)
                        end if
                    end do

                    ! check is the total number of NaN partial densities is == num_fluids. In this case, we assume
                    ! the solver will fail irrespective of the relaxation. There is no magic
                    if ( count( ieee_is_nan( m0k ) ) == num_fluids ) then
                      TR = .false.
                    end if

                    call s_correct_partial_densities(2, alphak, me0k, m0k, rM, rho, TR, i, j, k, l)

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
                            call s_old_infinite_p_relaxation(j, k, l, alphak, me0k, m0k, pS, rhoe, Tk)
                        case (4) ! p-equilibrium
                            call s_infinite_p_relaxation(j, k, l, alphak, me0k, m0k, pS, rhoe, rM, Tk)
                        case (5) ! pT-equilibrium
                            ! for this case, MFL cannot be either 0 or 1, so I chose it to be 2
                            call s_infinite_pt_relaxation(j, k, l, m0k, 2, pS, p_infpT, rhoe, rM, TS)
                            Tk = spread(TS, 1, num_fluids)
                        case (6) ! pT-pTg equilibrium
                            ! pT-equilibrium as rhe initial condition
                            call s_infinite_pt_relaxation(j, k, l, m0k, 2, pS, p_infpT, rhoe, rM, TS)
                            Tk = spread(TS, 1, num_fluids)

                            ! updating the densities and volume fractions used for thresholds
                            rhok = (pS + ps_inf)/((gs_min - 1)*cvs*Tk)

                            ! new volume fractions, after partial densities and p- or pT-equilibrium
                            alphak = m0k / rhok

                            !! phase change triggers
                            ! For Interface Capturing (enough alphak)
                            TIC = (alphak(lp) > palpha_eps) .and. (alphak(vp) > palpha_eps)

                            ! For Subgrid (enough alpha_b)
                            ! in case interface capturing and subgrid are activated. Subgrid trigger
                            alpha_b = q_cons_vf(alf_idx)%sf(j, k, l) ; TSG = .false.
                            if (bubbles_euler .and. (alphak(lp) > alpha_b) .and. ( alphak(lp) > 0.5 )) then
                              ! Vapor pressure from the intermediate pT state is
                              ! only needed here to evaluate the Blake/subgrid
                              ! trigger before deciding whether pTg relaxation is
                              ! activated.
                              call s_Saturation_Properties(pVapSG, TS, pS, 2)

                              do cb = 1, nb

                                ! this is true for the monodisperse case, for the moment. I need to expand this to 'R0ref(cb)'
                                ! mass_b(cb) = rho0ref * 4.0_wp * pi * R0ref ** 3.0_wp / 3.0_wp

                                R_b(cb) = q_cons_vf(bub_idx%rs(cb))%sf(j, k, l) / q_cons_vf(n_idx)%sf(j, k, l)

                                ! print *, 'Volume fraction: ', alphak(lp)
                                ! print *, 'Volume fraction bubble: ', alpha_b

                                call s_SG_trigger( alpha_b, mass_b(cb), pS, R_b(cb), pVapSG, TSG )

                              end do
                            end if

                            ! 1 - model activation, 1st order transition (p,T) <= (pCr, TCr)
                            if ( ( pS < pCr ) .and. &
                            ! 2.1 Homogeneous pTg-equilibrium criterium
                            ( ( ( pS < 0 ) .and. ( pS + minval(p_infpT) > 0.0_wp ) ) &
                            .or. &
                            ! 2.2. Heterogeneous pTg-equilibrium (either IC or SG activated).
                            TIC .or. TSG &
                            ) ) then
                                ! updating m1 and m2 AFTER correcting the partial densities. These values are
                                ! stored to be retrieved in case the final state is a mixture of fluids
                                mOr = (/ m0k(lp), m0k(vp) /)

                                ! checking if fluid is either subcoooled liquid or overheated vapor (NOT metastability)

                                ! overheated vapor
                                ! depleting the mass of liquid and tranferring the total mass to vapor
                                m0k(lp) = mixM*rM ; m0k(vp) = (1.0_wp - mixM)*rM

                                ! calling pT-equilibrium for overheated vapor, which is MFL = 0
                                call s_infinite_pt_relaxation(j, k, l, m0k, 0, pSOV, p_infOV, rhoe, rM, TSOV)

                                ! computing the saturation temperature for the overheated vapor state
                                call s_Saturation_Properties(pSOV, TSatOV, TSOV, 1)

                                ! subcooled liquid
                                ! tranferring the total mass to liquid and depleting the mass of vapor
                                m0k(lp) = (1.0_wp - mixM)*rM ; m0k(vp) = mixM*rM

                                ! calling pT-equilibrium for subcooled liquid, which is MFL = 1
                                call s_infinite_pt_relaxation(j, k, l, m0k, 1, pSSL, p_infSL, rhoe, rM, TSSL)

                                ! computing the saturation temperature for the subcooled liquid state
                                call s_Saturation_Properties(pSSL, TSatSL, TSSL, 1)

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
! #ifndef MFC_OpenACC
!                                     print *, 'case 6 fallback to pTg: j,k,l = ', j, k, l
!                                     print *, 'pS = ', pS, ' pCr = ', pCr
!                                     print *, 'TIC = ', TIC, ' TSG = ', TSG
!                                     print *, 'TSOV = ', TSOV
!                                     print *, 'TSatOV = ', TSatOV
!                                     print *, 'TSOV - TSatOV = ', TSOV - TSatOV
!                                     print *, 'TSSL = ', TSSL
!                                     print *, 'TSatSL = ', TSatSL
!                                     print *, 'TSSL - TSatSL = ', TSSL - TSatSL
!                                     print *, 'm0k(lp) = ', m0k(lp), ' m0k(vp) = ', m0k(vp)
! #endif
                                    ! returning partial pressures to what they were after the partial density correction
                                    m0k(lp) = mOr(1) ; m0k(vp) = mOr(2)

                                    ! pTg-relaxation
                                    call s_infinite_ptg_relaxation(j, k, l, alphak, me0k, m0k, pS, p_infpT, rho, rhoe, rM, TR, TS, TSG)
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

                                    ! if ( m0k(vp) > 0.0_wp ) then
                                    !     print *, 'alphak(lp):', alphak(lp)
                                    !     print *, 'alphak(vp):', alphak(vp)
                                    !     print *, 'alpha_b:', alpha_b
                                    ! end if

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
                    ! Update conservative variables only when the relaxation
                    ! path remains active for this cell.
                    if (TR) then 
                        call update_conservative_vars( j, k, l, m0k, pS, q_cons_vf, Tk )
                        ! Store the final relaxed-cell vapor pressure
                        ! for the next Euler bubble update.
                        if (bubbles_euler .and. pv > 0.0_wp) then
                            call s_Saturation_Properties(pVap_sf(j, k, l), TS, pS, 2)
                        end if
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
    end subroutine s_infinite_relaxation_k ! ----------------

    !>  This auxiliary subroutine is created to activate the pT-equilibrium for N fluids
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param pS equilibrium pressure at the interface
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rhoe mixture energy
    impure subroutine s_infinite_p_relaxation(j, k, l, alpha0k, me0k, m0k, pS, rhoe, rM, Tk)
        $:GPU_ROUTINE(function_name='s_infinite_p_relaxation', &
            & parallelism='[seq]', cray_inline=True)

        ! initializing variables
        real(wp), intent(in) :: rhoe, rM
        real(wp), intent(out) :: pS
        real(wp), dimension(num_fluids), intent(out) :: Tk
        real(wp), dimension(num_fluids), intent(in)  :: alpha0k, me0k, m0k
        integer, intent(in) :: j, k, l

        real(wp) :: fp, fpp !< variables for the Newton Solver
        real(wp) :: Econst, Om, TS !< auxiliary variables
        real(wp), dimension(num_fluids) :: alphak, mek, meik
        logical, dimension(num_fluids) :: is_negligible_mass
        character(20) :: nss, pSs, Econsts
        integer, dimension(num_fluids) :: iFix, iAuxSP, iAuxZP !< auxiliary index for choosing appropiate values for conditional sums
        real(wp), dimension(3) :: Oc
        integer, dimension(:), allocatable :: iSP, iZP

        integer :: mF !< multiplying factor for the tolerance of the solver
        integer :: i, na, ns, nsL !< generic loop iterators

        ! indices for all the fluids/phases
        iFix = (/ (i, i=1,num_fluids) /)

        is_negligible_mass = f_is_negligible_phase_mass(m0k, rM)

        ! indices for zero-mass phases (negligible amount of partial density). Fluids with negative partial densities should
        ! not be present at this point, since they have already been corrected at the first call of s_correct_partial_densities
        iAuxZP = merge(iFix, 0, is_negligible_mass)
        iZP = pack(iAuxZP, iAuxZP /= 0)

        ! indices for phases that have a significant partial density
        iAuxSP = merge(iFix, 0, .not. is_negligible_mass)
        iSP = pack(iAuxSP, iAuxSP /= 0)

        ! Re-distributing the initial internal energies such that rhoe = rhoe6E. This step might be unecessary in the future
        meik = me0k * rhoe / sum(me0k)

        ! initial conditions for starting the solver. For pressure, as long as the initial guess
        ! is in (-min(gs_min*ps_inf), +infty), a solution should be found.
        pS = (rhoe - sum( m0k(iSP) * qvs(iSP) ) - sum( alpha0k(iSP) * pi_infs(iSP) ) ) / sum( alpha0k(iSP) * gammas(iSP) )

        if ( pS < -1.0_wp * minval( gs_min(iSP) * ps_inf(iSP) ) ) then
          pS = -1.0_wp * minval( gs_min(iSP) * ps_inf(iSP) ) + ptgalpha_eps * 101325
        end if

        ! internal energies - first estimate
        mek = meik * ( 1 + ptgalpha_eps )

        ! volume fractions - first estimate
        alphak = alpha0k * ( 1 + ptgalpha_eps ) / sum( alpha0k * ( 1 + ptgalpha_eps ) )

        ! counter for the outer loop
        nsL = 0

        ! Relaxation factor. Although this is not needed for Newton Solver for finding p, it seems to be needed to update
        ! the internal energies after finding pS.
        Om = under_relax

        do while ( ( ( abs( sum( mek(iSP) ) - rhoe ) > ptgalpha_eps ) .and. ( abs( ( sum( mek(iSP) ) - rhoe ) / rhoe ) > ptgalpha_eps ) ) .or.  ( nSL == 0 ) )
            ! increasing counter
            nsL = nsL + 1

            ! Variable to check the energy constraint before initializing the p-relaxation procedure. This ensures
            ! global convergence will be estabilished
            Econst = sum( (gs_min(iSP) - 1.0_wp) * ( mek(iSP) - m0k(iSP) * qvs(iSP) ) / ( gs_min(iSP) * ps_inf(iSP) - minval( ps_inf(iSP) ) ) )

#ifndef MFC_OpenACC
            ! energy constraint for the p-equilibrium
            if ((minval( ps_inf(iSP) ) > 0) .and. (Econst <= 1.0_wp) .or. (nsL > max_iter)) then

              call s_whistleblower((/ 0.0_wp,  0.0_wp/), (/ (/1/fpp, 0.0_wp/), (/0.0_wp, 0.0_wp/) /), j &
                                , (/ (/fpp, 0.0_wp/), (/0.0_wp, 0.0_wp/) /), k, l, m0k, nsL, ps_inf &
                                , pS, (/ sum( mek ) - rhoe, 0.0_wp/), rhoe, alphak * (pS + ps_inf) / ( (gs_min - 1.0_wp) * m0k * cvs ))

              call s_real_to_str(Econst, Econsts)
              call s_mpi_abort('Solver for the p-relaxation solver failed (m_phase_change, s_infinite_p_relaxation) &
&                   . Please, check energy constraint. Econst ~'//Econsts//'. Aborting!')
            end if
#endif

            ! if the energy constraint is satisfied, we start the newton solver counter, and the multiplier
            ns = 0 ; mF = 1

            ! Newton solver for p-equilibrium. ns <= 1, is to ensure the internal energy correction happens at least once.
            ! in the loosely coupled algorithm. This is different than the pT-equilibrium case, in which no energy correction is needed.
            ! A solution is found when f(p) = 1
            fp = 0.0_wp
            do while ( ( ( abs(fp - 1.0_wp) > mF * ptgalpha_eps ) ) .or. ( ns <= 1 ) )
                ! increasing counter
                ns = ns + 1

                ! updating functions used in the Newton's solver. f(p)
                fp = sum( alphak(iSP) )

                ! updating functions used in the Newton's solver. f'(p)
                fpp = sum( -1.0_wp * alphak(iSP) / ( pS + gs_min(iSP) * ps_inf(iSP) ) )

                ! updating the relaxed pressure
                pS = pS + ( ( 1.0_wp - fp ) / fpp ) / ( 1.0_wp - ( 1.0_wp - fp + abs( 1.0_wp - fp ) ) / ( 2.0_wp * fpp * ( pS + minval( gs_min(iSP) * ps_inf(iSP) ) ) ) )

                ! updating the underelaxation parameters.
                ! First restriction
                if ( any( pS + gs_min(iSP) * ps_inf(iSP) > 0 ) ) then
                  Oc(1) = minval( ( meik(iSP) - m0k(iSP) * qvs(iSP) ) / ( pS * ( alphak(iSP) - alpha0k(iSP) ) ) ) / 2
                end if

                ! second restriction
                if ( any( pS + gs_min(iSP) * ps_inf(iSP) < 0 ) ) then
                  Oc(2) = maxval( ( meik(iSP) - m0k(iSP) * qvs(iSP) ) / ( pS * ( alphak(iSP) - alpha0k(iSP) ) ) ) / 1
                end if

                ! updating internal energies. An underrelaxation factor is needed due to the closure for mek
                if ( ( Om >= minval( (meik(iSP) - m0k(iSP) * qvs(iSP) ) / ( pS * (alphak(iSP) - alpha0k(iSP)) ) ) ) &
                .and. ( minval( (meik(iSP) - m0k(iSP) * qvs(iSP) ) / ( pS * (alphak(iSP) - alpha0k(iSP)) ) ) > 0 ) ) then
                  Om = minval( (meik(iSP) - m0k(iSP) * qvs(iSP) ) / ( pS * (alphak(iSP) - alpha0k(iSP)) ) ) / 2
                else
                  Om = under_relax
                end if

                ! updating phase variables, together with the relaxed pressure, in a loosely coupled procedure
                ! internal energies
                mek(iSP) = meik(iSP) - Om * pS * ( alphak(iSP) - alpha0k(iSP) )

                ! volume fractions
                alphak(iSP) = ( gs_min(iSP) - 1.0_wp ) * ( mek(iSP) - m0k(iSP) * qvs(iSP) ) / ( pS + gs_min(iSP) * ps_inf(iSP) )

                ! checking if pS is within expected bounds
                if ( ((pS <= -1.0_wp*minval( gs_min(iSP) * ps_inf(iSP) ) ) .or. (ieee_is_nan(pS))) .and. ( ns <= max_iter ) ) then

                  ! In case the newton-Raphson procedure for pS makes it <= -1.0_wp*minval(gs_min*ps_inf) due to the
                  ! estimates for the fluid internal energies, restart the pressure so that the solver can continue.
                  ! keep an eye on this, as it has not been tested
                  pS = (rhoe - sum( m0k(iSP) * qvs(iSP) ) - sum( alpha0k(iSP) * pi_infs(iSP) ) ) / sum( alpha0k(iSP) * gammas(iSP) )
                  if ( pS < -1.0_wp * minval( gs_min(iSP) * ps_inf(iSP) ) ) then
                    pS = -1.0_wp * minval( gs_min(iSP) * ps_inf(iSP) ) + ptgalpha_eps * 101325
                  end if

                  print *, 'pS restarted due to unphysical values pressures during the Newton solver. ns = ', ns, 'Continuing...'

                else if (ns > max_iter) then

                    ! restarting solver by increasing the tolerance. This is completely ad hoc. Keep an eye on it
                    if ( (abs(fp - 1.0_wp) < 10 * mF * ptgalpha_eps) .and. (abs(fp - 1.0_wp) > mF * ptgalpha_eps) ) then

                      ns = 0 ; mF = mF + 1

                      print *, 'ptgalpha_eps increased from ', (mF - 1) * ptgalpha_eps, ' to ', mF * ptgalpha_eps, &
                      '. Newton solver restarted'

                    else

                      call s_whistleblower((/ 0.0_wp,  0.0_wp/), (/ (/1/fpp, 0.0_wp/), (/0.0_wp, 0.0_wp/) /), j &
                                        , (/ (/fpp, 0.0_wp/), (/0.0_wp, 0.0_wp/) /), k, l, m0k, ns, ps_inf &
                                        , pS, (/fp - 1.0_wp, 0.0_wp/), rhoe, alphak * (pS + ps_inf) / ( (gs_min - 1.0_wp) * m0k * cvs ) )

                      call s_real_to_str(pS, pSs)
                      call s_int_to_str(ns, nss)
                      call s_mpi_abort('Solver for the p-relaxation failed (m_phase_change, s_infinite_p_relaxation). &
                      &   pS ~'//pSs//'. ns = '//nss//'. Aborting!')
                    end if
                end if
            end do
        end do

        ! (NOT common) temperatures
        Tk(iSP) = alphak(iSP) * (pS + ps_inf(iSP)) / ( (gs_min(iSP) - 1.0_wp) * m0k(iSP) * cvs(iSP) )

        ! correcting nonphysical temperatures due to alphak/m0k = 0/0 division. Note that the state for these fluids
        ! need not be determined, as alpha = alpharho = alpharhoe = 0 (state components for the q array). They cannot only
        ! be either 0 or NaN for the sake of the algorithm.
        Tk(iZP) = (pS + ps_inf(iZP)) / ( (gs_min(iZP) - 1.0_wp) * cvs(iZP) )

        ! updating maximum number of iterations
        max_iter_pc_ts = maxval((/max_iter_pc_ts, ns/))

    end subroutine s_infinite_p_relaxation ! -----------------------

    ! Description: The purpose of this procedure is to infinitely relax
    !              the pressures from the internal-energy equations to a
    !              unique pressure, from which the corresponding volume
    !              fraction of each phase are recomputed. For conservation
    !              purpose, this pressure is finally corrected using the
    !              mixture-total-energy equation.
    ! Relaxed pressure, initial partial pressures, function f(p) and its partial
    ! derivative df(p), isentropic partial density, sum of volume fractions,
    ! mixture density, dynamic pressure, surface energy, specific heat ratio
    ! function, liquid stiffness function (two variations of the last two
    ! initializing variables
    impure subroutine s_old_infinite_p_relaxation(j, k, l, alpha0k, me0k, m0k, pS, rhoe, Tk)

        $:GPU_ROUTINE(function_name='s_old_infinite_p_relaxation', &
            & parallelism='[seq]', cray_inline=True)
        real(wp), intent(in) :: rhoe
        real(wp), intent(out) :: pS
        real(wp), dimension(num_fluids), intent(out) :: Tk
        real(wp), dimension(num_fluids), intent(in) :: alpha0k, me0k, m0k

        integer, intent(in) :: j, k, l

        real(wp) :: fp, fpp, mQ, pO !< variables for the Newton Solver
        real(wp), dimension(num_fluids) :: alphak, mek, mk, pk, rhok, num, den, drhodp
        integer, dimension(num_fluids) :: iVar, iFix !< auxiliary index for choosing appropiate values for conditional sums
        character(20) :: nss, pSs
        !> @}

        integer :: i, ns !< generic loop iterators

        ! indices for all the fluids/phases
        iFix = (/ (i, i=1,num_fluids) /)

        ! initializing the partial energies, volume fractions, and masses
        alphak = alpha0k ; mek = me0k ; mk = m0k

        ! Numerical correction of the volume fractions, partial densities, and internal energies
        if (mpp_lim) then
          ! auxiliry variable to avoid do loops. Correting only variables for fluids whose alpha < 0 OR m < 0
          iVar = iFix ; iVar( pack( iFix, ( alpha0k >= 0 ) .and. ( m0k >= 0 ) ) ) = 0

          ! if alpha < 0 or m < 0, set everything to 0
          mk( pack( iVar, iVar /= 0 ) ) = 0.0_wp
          mek( pack( iVar, iVar /= 0 ) ) = 0.0_wp
          alphak( pack( iVar, iVar /= 0 ) ) = 0.0_wp
          ! if alpha > 1, set it to 1
          alphak( pack( iFix, alpha0k > 1.0_wp ) ) = 1.0_wp

          ! renormalize alpha
          alphak = alphak / sum(alphak)
        end if

        ! Initial state
        pk = 0.0_wp

        ! auxiliry variable to avoid do loops. Disregarding variables when alpha <= sgm_eps
        iVar = iFix ; iVar( pack( iFix, alphak <= sgm_eps ) ) = 0

        pk( pack( iVar, iVar /= 0 ) ) = ( ( mek( pack( iVar, iVar /= 0 ) ) &
        - mk( pack( iVar, iVar /= 0 ) ) * qvs( pack( iVar, iVar /= 0 ) ) ) /      &
        alphak( pack( iVar, iVar /= 0 ) ) - pi_infs( pack( iVar, iVar /= 0 ) ) ) &
        / gammas( pack( iVar, iVar /= 0 ) )

        ! auxiliry variable to avoid do loops. Disregarding pressures that are not physical
        iVar = iFix ; iVar( pack( iFix, .not. ( pk < -(1.0_wp - ptgalpha_eps)*ps_inf + ptgalpha_eps ) ) ) = 0

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
        do while ( abs( fp ) > ptgalpha_eps .or. ( ns == 0 ) )

            ! setting old pressure
            pO = pS

            ! Newton-Raphson method
            ! Physical pressure?
            if ( pO < minval(-(1.0_wp - ptgalpha_eps)*ps_inf + ptgalpha_eps) ) then
              pO = minval(-(1.0_wp - ptgalpha_eps)*ps_inf + ptgalpha_eps)
            end if

            ! calulating fp and fpp. Disregarding fluids for alpha < sgm_eps
            iVar = iFix ; iVar( pack( iFix, alphak <= sgm_eps ) ) = 0

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

                  call s_whistleblower((/ -fp/fpp,  0.0_wp/), (/ (/1/fpp, 0.0_wp/), (/0.0_wp, 0.0_wp/) /), j &
                                    , (/ (/fpp, 0.0_wp/), (/0.0_wp, 0.0_wp/) /), k, l, m0k, ns, ps_inf &
                                    , pS, (/fp, 0.0_wp/), rhoe, Tk)

                  call s_real_to_str(pS, pSs)
                  call s_int_to_str(ns, nss)
                  call s_mpi_abort('Solver for the Old p-relaxation failed (m_phase_change, s_old_infinite_p_relaxation). &
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
        pS = (rhoe - sum( alphak * pi_infs ) - sum( mk * qvs ) ) / sum( alphak * gammas )

        ! (NOT common) temperature
        Tk = (pS + ps_inf)/((gs_min - 1)*cvs*rhok)

        ! updating maximum number of iterations
        max_iter_pc_ts = maxval((/max_iter_pc_ts, ns/))

    end subroutine s_old_infinite_p_relaxation ! -----------------------

    !>  This auxiliary subroutine is created to activate the pT-equilibrium for N fluids
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
        !!  @param MFL flag that tells whether the fluid is gas (0), liquid (1), or a mixture (2)
        !!  @param pS equilibrium pressure at the interface
        !!  @param p_infpT stiffness for the participating fluids under pT-equilibrium
        !!  @param rM sum of the reacting masses
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rhoe mixture energy
        !!  @param TS equilibrium temperature at the interface
    subroutine s_infinite_pt_relaxation(j, k, l, m0k, MFL, pS, p_infpT, rhoe, rM, TS)

        $:GPU_ROUTINE(function_name='s_infinite_pt_relaxation', &
            & parallelism='[seq]', cray_inline=True)

        ! initializing variables
        real(wp), intent(out) :: pS, TS
        real(wp), dimension(num_fluids), intent(out) :: p_infpT
        real(wp), intent(in) :: rhoe, rM
        real(wp), intent(in), dimension(num_fluids) :: m0k
        logical, dimension(num_fluids) :: is_negligible_mass
        integer, intent(in) :: j, k, l, MFL
        integer, dimension(num_fluids) :: iFix, iAuxSP, iAuxZP !< auxiliary index for choosing appropiate values for conditional sums
        integer, dimension(:), allocatable :: iSP, iZP
        real(wp) :: gp, gpp, hp, pO, mCP, mQ !< variables for the Newton Solver
        character(20) :: nss, pSs, Econsts

        integer :: i, ns !< generic loop iterators

        ! auxiliary variables for the pT-equilibrium solver
        p_infpT = ps_inf

        ! indices for all the fluids/phases
        iFix = (/ (i, i=1,num_fluids) /)

        is_negligible_mass = f_is_negligible_phase_mass(m0k, rM)

        ! indices for zero-mass phases (negligible amount of partial density). Fluids with negative partial densities should
        ! not be present at this point, since they have already been corrected at the first call of s_correct_partial_densities
        iAuxZP = merge(iFix, 0, is_negligible_mass)
        iZP = pack(iAuxZP, iAuxZP /= 0)

        ! indices for phases that have a significant partial density
        iAuxSP = merge(iFix, 0, .not. is_negligible_mass)
        iSP = pack(iAuxSP, iAuxSP /= 0)

        ! this value is rather arbitrary, as I am interested in MINVAL( ps_inf ) for the solver.
        ! This way, I am ensuring this value will not be selected.
        p_infpT(iZP) = 2 * maxval( ps_inf )

        ! sum of the total alpha*rho*cp of the system. Note that these variables already dismiss the subgrid fluid
        ! physical parameters
        mCP = sum( m0k(iSP) * cvs(iSP) * gs_min(iSP) )

        ! sum of the total alpha*rho*q of the system
        mQ = sum( m0k(iSP) * qvs(iSP) )

        ! Checking energy constraint. In the case we are calculating the possibility of having subcooled liquid or
        ! overheated vapor, the energy constraint might not be satisfied, as are hypothetically transferring all the
        ! mass from one phase to the other. When this is the case, we simply ignore this possibility, set pS = TS = 0,
        ! and discard the hypothesis. The solver can thus move forward.
        if ((rhoe - mQ - minval(ps_inf(iSP))) < 0.0_wp) then

            if ( any((/ 0, 1 /) == MFL ) ) then
! #ifndef MFC_OpenACC
!                 print *, 'pT trial rejected by energy constraint: j,k,l = ', j, k, l
!                 print *, 'MFL = ', MFL
!                 print *, 'm0k(lp) = ', m0k(lp), ' m0k(vp) = ', m0k(vp)
!                 print *, 'rhoe = ', rhoe
!                 print *, 'mQ = ', mQ
!                 print *, 'minval(ps_inf(iSP)) = ', minval(ps_inf(iSP))
!                 print *, 'energy_margin = ', rhoe - mQ - minval(ps_inf(iSP))
! #endif

                ! Assigning zero values for pressure and temperature in case of mass depletion cases
                pS = 0.0_wp ; TS = 0.0_wp

                return
#ifndef MFC_OpenACC
            else
                call s_whistleblower((/ 0.0_wp,  0.0_wp/), (/ (/0.0_wp, 0.0_wp/), (/0.0_wp, 0.0_wp/) /), j &
                                  , (/ (/0.0_wp, 0.0_wp/), (/0.0_wp, 0.0_wp/) /), k, l, m0k, ns, p_infpT &
                                  , 0.0_wp, (/0.0_wp, 0.0_wp/), rhoe, spread(0.0_wp, 1, num_fluids))

                call s_real_to_str(rhoe - mQ - minval(ps_inf(iSP)), Econsts)
                call s_mpi_abort('Solver for the pT-relaxation solver failed (m_phase_change, s_infinite_pt_relaxation) &
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
        gp = 0.0_wp

        ! Newton solver for pT-equilibrium. gdi = 0.0_wp is arbitrary, and ns == 0, to the loop is entered at least once.
        ! Note that a solution is found when gp(p) = 1
        do while ( ( ( abs( gp - 1.0_wp ) > ptgalpha_eps ) ) .or. ( ns == 0 ) )

            ! increasing counter
            ns = ns + 1

            ! updating old pressure
            pO = pS

            ! updating functions used in the Newton's solver
            gp = sum( ( gs_min(iSP) - 1.0_wp ) * m0k(iSP) * cvs(iSP) * ( rhoe + pO - mQ ) / ( mCP * ( pO + ps_inf(iSP) ) ) )

            gpp = sum( ( gs_min(iSP) - 1.0_wp ) * m0k(iSP) * cvs(iSP) * ( ps_inf(iSP) - rhoe + mQ ) / ( mCP * ( pO + ps_inf(iSP) ) **2 ) )

            hp = 1.0_wp/(rhoe + pO - mQ) + 1.0_wp/(pO + minval(ps_inf(iSP)))

            ! updating common pressure for the newton solver
            pS = pO + ((1.0_wp - gp) / gpp) / (1.0_wp - (1.0_wp - gp + abs(1.0_wp - gp)) * hp / (2.0_wp*gpp))

            ! common temperature
            TS = (rhoe + pS - mQ) / mCP

            ! check if solution is out of bounds (which I believe it won`t happen given the solver is gloabally convergent.
#ifndef MFC_OpenACC
            if ((pS <= -1.0_wp*minval(ps_inf(iSP))) .or. (ieee_is_nan(pS)) .or. (ns > max_iter)) then

              call s_whistleblower((/0.0_wp, 0.0_wp/), (/ (/1/gpp, 0.0_wp/), (/0.0_wp, 0.0_wp/) /), j &
                                , (/ (/gpp, 0.0_wp/), (/0.0_wp, 0.0_wp/) /), k, l, m0k, ns, p_infpT &
                                , pS, (/abs( gp - 1.0_wp ), 0.0_wp/), rhoe, spread(TS, 1, num_fluids))

              call s_real_to_str(pS, pSs); call s_int_to_str(nS, nss)
              call s_mpi_abort('Solver for the pT-relaxation failed (m_phase_change, s_infinite_pt_relaxation). &
              &   pS ~'//pSs//'. ns = '//nss//'. Aborting!')
            end if
#endif
        end do

        ! updating maximum number of iterations
        max_iter_pc_ts = maxval((/max_iter_pc_ts, ns/))

    end subroutine s_infinite_pt_relaxation ! -----------------------

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
    subroutine s_infinite_ptg_relaxation(j, k, l, alphak, me0k, m0k, pS, p_infpT, rho, rhoe, rM, TR, TS, TSG)

        $:GPU_ROUTINE(function_name='s_infinite_ptg_relaxation', &
            & parallelism='[seq]', cray_inline=True)

        real(wp), dimension(num_fluids), intent(inout) :: alphak, me0k, m0k
        real(wp), dimension(num_fluids), intent(in) :: p_infpT
        real(wp), intent(inout) :: pS, TS, rM
        real(wp), intent(in) :: rhoe
        integer, intent(in) :: j, k, l
        logical, intent(inout) :: TR, TSG ! triggering parameters
        real(wp), dimension(num_fluids) :: p_infpTg, hk, gk, sk
        real(wp), dimension(2, 2) :: Jac, InvJac, TJac
        real(wp), dimension(2) :: R2D, DeltamP
        real(wp), dimension(3) :: Oc
        real(wp), parameter :: Om_floor = 1.0e-12_wp ! minimum positive relaxation factor
        real(wp) :: Om ! underrelaxation factor
        real(wp) :: maxg, mCP, mCPD, mCVGP, mCVGP2, mQ, mQD, rho, TSat ! auxiliary variables for the pTg-solver
        character(20) :: nss, pSs, Econsts, R2D1s, R2D2s

        !< Generic loop iterators
        integer :: i, ns

        ! assigning the relexant pi_infs based on the previous pT-equilibrium
        p_infpTg = p_infpT

        ! checking if homogeneous cavitation is expected. If yes, transfering a small amount of mass to the depleted
        ! phase, and then let the algorithm run.

        ! is the fluid at a metastable state with enough 'energy' for phase change to happen? Or, is the subgrid bubble
        ! volume fraction large enough?
        if ( ( (pS < -1.47e10_wp) .and. (rM > (rhoe - gs_min(lp)*ps_inf(lp)/(gs_min(lp) - 1.0e-1_wp))/qvs(lp)) ) .or. &
        TSG ) then

            ! transfer a bit of mass to the deficient phase, enforce phase change
            call s_correct_partial_densities(1, alphak, me0k, m0k, rM, rho, TR, i, j, k, l)

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

        ! Relaxation factor. This value is initially user-defined, with a certain level of self adjustment.
        Oc = under_relax;

        R2D = 0.0_wp ; DeltamP = 0.0_wp;
        ! starting counter for the Newton solver
        ns = 0

        ! (initial) common temperature
        TS = (rhoe + pS - mQ)/mCP

        ! entropy
        sk = cvs*log((TS**gs_min)/((pS + ps_inf)**(gs_min - 1.0_wp))) + qvps

        ! enthalpy
        hk = gs_min*cvs*TS + qvs

        ! Gibbs-free energy.
        gk = hk - TS*sk

        ! maximum Gibbs Free Energy for the reacting phase, used as a relative criterion for the solver
        maxg = maxval([gk(lp),gk(vp)])

        ! Newton solver for pTg-equilibrium. 1d6 is arbitrary, and ns == 0, to the loop is entered at least once.
        do while ( ( ( norm2(R2D) > ptgalpha_eps ) .and. ( norm2( R2D*(/maxg,rhoe/) ) / norm2( (/maxg,rhoe/) ) > ptgalpha_eps ) ) &
          .or. ( ns == 0 ) )

            ! Updating counter for the iterative procedure
            ns = ns + 1

            ! Auxiliary variables to help in the calculation of the residue
            ! Those must be updated through the iterations, as they either depend on
            ! the partial masses for all fluids, or on the equilibrium pressure

            ! sum of the total alpha*rho*cp of the system
            mCP = sum( m0k * cvs * gs_min )

            ! mCP - the contribution from the reacting phases
            mCPD = mCP - m0k(lp) * cvs(lp) * gs_min(lp) - m0k(vp) * cvs(vp) * gs_min(vp)

            ! sum of the total alpha*rho*q of the system
            mQ = sum( m0k * qvs )

            ! mQ - the contribution from the reacting phases
            mQD = mQ - m0k(lp) * qvs(lp) - m0k(vp) * qvs(vp)

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

            ! calculating correction array for Newton's method
            DeltamP = matmul(InvJac, R2D)

            ! checking if the correction in the mass/pressure will lead to negative values for those quantities
            ! If so, adjust the underrelaxation parameter Om
#ifndef MFC_OpenACC
            ! reset the trial factor at the start of each Newton iteration
            Om = under_relax

            ! creating criteria for variable underrelaxation factor
            if (m0k(lp) - Om*DeltamP(1) <= 0.0_wp) then
                Oc(1) = m0k(lp)/(2*DeltamP(1))
            else
                Oc(1) = under_relax
            end if
            if (m0k(vp) + Om*DeltamP(1) <= 0.0_wp) then
                Oc(2) = -m0k(vp)/(2*DeltamP(1))
            else
                Oc(2) = under_relax
            end if
            if (pS + minval(p_infpTg) - Om*DeltamP(2) <= 0.0_wp) then
                Oc(3) = (pS + minval(p_infpTg))/(2*DeltamP(2))
            else
                Oc(3) = under_relax
            end if
            ! choosing amongst the minimum relaxation maximum to ensure solver will not produce unphysical values.
            ! If the limiter becomes nonpositive, fall back to a tiny positive step instead of reversing direction.
            if (minval(Oc) > 0.0_wp) then
                Om = min(under_relax, minval(Oc))
            else
                Om = Om_floor
            end if
#else
            Om = under_relax
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
            call s_compute_pTg_residual(j, k, l, m0k, mCPD, mCVGP, mQD, pS, rhoe, rM, R2D)

            ! updating common temperature
            TS = (rhoe + pS - mQ)/mCP

            ! entropy
            sk = cvs*log((TS**gs_min)/((pS + ps_inf)**(gs_min - 1.0_wp))) + qvps

            ! enthalpy
            hk = gs_min*cvs*TS + qvs

            ! Gibbs-free energy
            gk = hk - TS*sk

            ! maximum Gibbs Free Energy for the reacting phase, used as a relative criterion for the solver
            maxg = maxval([gk(lp),gk(vp)])

            ! print *, 'j,k,l = ', j, k, l
            ! print *, 'ns = ', ns
            ! print *, 'mCP = ', mCP
            ! print *, 'mQ = ', mQ
            ! print *, 'TS = ', TS
            ! print *, 'Om = ', Om
            ! print *, 'DeltamP = ', DeltamP
            ! print *, 'pS = ', pS
            ! print *, 'pS + minval(p_infpTg) = ', pS + minval(p_infpTg)
            ! print *, 'pS + ps_inf(lp) = ', pS + ps_inf(lp)
            ! print *, 'pS + ps_inf(vp) = ', pS + ps_inf(vp)
            ! print *, 'R2D = ', R2D
            ! print *, 'norm2(R2D) = ', norm2(R2D)
            ! print *, 'det(J) = ', Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1)
            ! print *, 'Jac = ', Jac
            ! print *, 'mCPD = ', mCPD
            ! print *, 'mQD = ', mQD
            ! print *, 'mCVGP = ', mCVGP
            ! print *, 'mCVGP2 = ', mCVGP2
            ! print *, 'm0k(lp) = ', m0k(lp)
            ! print *, 'm0k(vp) = ', m0k(vp)
            ! print *, 'maxg = ', maxg
            ! print *, 'relative residual = ', norm2(R2D*(/maxg,rhoe/))/norm2((/maxg,rhoe/))


          ! checking if the residue returned any NaN values
#ifndef MFC_OpenACC
          if (ieee_is_nan(norm2(R2D)) .or. (ns > max_iter)) then

            call s_whistleblower(DeltamP, InvJac, j, Jac, k, l, m0k, ns, p_infpTg &
                                , pS, R2D, rhoe, spread(TS, 1, num_fluids))

            call s_real_to_str(R2D(1), R2D1s) ; call s_real_to_str(R2D(2), R2D2s)
            call s_real_to_str(rhoe - mQ - minval(p_infpTg), Econsts)
            call s_int_to_str(ns, nss); call s_real_to_str(pS, pSs)
            call s_mpi_abort('Solver for the pTg-relaxation failed (m_phase_change, s_infinite_ptg_relaxation). &
            &. ns = ' //nss//', R2D(1) = '//R2D1s//', R2D(2) = '//R2D2s//'. pS ~'//pSs//' &
            . Econst = '//Econsts//'. Aborting!')
          end if
#endif
        end do

        ! updating maximum number of iterations
        max_iter_pc_ts = maxval((/max_iter_pc_ts, ns/))

    end subroutine s_infinite_ptg_relaxation ! -----------------------

    !>  This auxiliary subroutine corrects the partial densities of the REACTING fluids in case one of them is negative
        !!      but their sum is positive. Inert phases are not corrected at this moment
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param rM sum of the reacting masses
        !!  @param j generic loop iterator for x direction
        !!  @param k generic loop iterator for y direction
        !!  @param l generic loop iterator for z direction
    subroutine s_correct_partial_densities(CT, alpha0k, me0k, m0k, rM, rho, TR, i, j, k, l)
        $:GPU_ROUTINE(function_name='s_correct_partial_densities', &
            & parallelism='[seq]', cray_inline=True)

        !> @name variables for the correction of the reacting partial densities
        !> @{
        real(wp), dimension(num_fluids), intent(inout) :: alpha0k, me0k, m0k
        real(wp), intent(out) :: rM, rho
        logical, intent(inout) :: TR
        integer, intent(in) :: CT, j, k, l
        integer, dimension(num_fluids) :: iFix, iAuxZP !< auxiliary index for choosing appropiate values for conditional sums
        integer, dimension(:), allocatable :: iZP
        integer :: i
        !> @}

        ! indices for all the fluids/phases
        iFix = (/ (i, i=1,num_fluids) /)
        iAuxZP = iFix

        ! Mixture density
        rho = sum(m0k)

        ! calculating the total reacting mass for the phase change process. By hypothesis, this should not change
        ! throughout the phase-change process. In ANY HYPOTHESIS, UNLESS correcting them artifically
        rM = m0k(lp) + m0k(vp)

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
            ! zero phase indices. Auxiliary variable to avoid do loops - use Morgan's Law
            if ( any((/ 1, 4 /) == relax_model ) ) then
                ! this iAuxZP is only valid when we use either the old or new p-relaxations, as they are only
                ! used with the 6-equation model. Note that they test the phisical validity of the initial conditions
                ! iAuxZP( pack( iFix, ( alpha0k > 0 ) .and. ( m0k > 0 ) .and. ( me0k > m0k * qvs ) ) ) = 0
                iAuxZP( pack( iFix, ( m0k > 0 ) ) ) = 0
            else
                ! this is used for either pT- or pTg-relaxation, as regardless of the equation model, the phasic internal
                ! energies are not important
                iAuxZP( pack( iFix, ( alpha0k > 0 ) .and. ( m0k > 0 ) ) ) = 0
            end if
            iZP = pack(iAuxZP, iAuxZP /= 0)

            ! if either the volume fraction or the partial density is negative, make them positive
            alpha0k(iZP) = 0.0_wp

            ! the largest value of alpha0k must be one
            alpha0k( pack( iFix, alpha0k > 1.0_wp ) ) = 1.0_wp

            m0k(iZP) = 0.0_wp

            if (model_eqns == 3) then
              me0k( pack( iFix, me0k < 0.0_wp ) ) = 0.0_wp
              ! me0k(iZP) = 0.0_wp
              ! me0k = me0k / sum(alpha0k)
            end if

            ! renormalizing alpha. This should be the last correction, as the
            ! other variables would not benefit from the proportional adjustment
            alpha0k = alpha0k / sum(alpha0k)

            ! continue relaxation
            TR = .true.
        else
            ! if there are any insignificant values, make them significant
            m0k( pack( iFix, m0k < rho * mixM ) ) = rho * mixM

            ! continue relaxation
            TR = .true.
        end if

        ! Mixture density
        rho = sum(m0k)

        ! calculating the total reacting mass for the phase change process.
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
    subroutine s_compute_pTg_residual(j, k, l, m0k, mCPD, mCVGP, mQD, pS, rhoe, rM, R2D)
        $:GPU_ROUTINE(function_name='s_compute_pTg_residual', &
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
                  + (m0k(lp)*(cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                  - rM*cvs(vp)*gs_min(vp) - mCPD) * TS ) / 1

    end subroutine s_compute_pTg_residual

    ! SUBROUTINE CREATED TO TELL ME WHERE THE ERROR IN THE PT- AND PTG-EQUILIBRIUM SOLVERS IS
    impure subroutine s_whistleblower(DeltamP, InvJac, j, Jac, k, l, mk, ns, p_inf, pS, R2D, rhoe, Tk) ! ----------------

        real(wp), dimension(2, 2), intent(in) :: Jac, InvJac
        real(wp), dimension(num_fluids), intent(in) :: mk, p_inf, Tk
        real(wp), dimension(2), intent(in) :: R2D, DeltamP
        real(wp), intent(in) :: pS, rhoe
        integer, intent(in) :: j, k, l, ns
        real(wp), dimension(num_fluids) :: alphak, ek, hk, gk, sk, rhok
        integer, dimension(num_fluids) :: iFix
        real(wp) :: maxg, rho
        !< Generic loop iterator
        integer :: i

        ! auxiliary variable to correct
        iFix = (/ (i, i=1,num_fluids) /)

        ! auxiliary calculations
        ! Thermodynamic state
        ! entropy
        sk = cvs*log((Tk**gs_min)/((pS + ps_inf)**(gs_min - 1.0_wp))) + qvps

        ! enthalpy
        hk = gs_min*cvs*Tk + qvs

        ! Gibbs-free energy
        gk = hk - Tk*sk

        ! maximum Gibbs Free Energy for the reacting phase, used as a relative criterion for the solver
        maxg = maxval([gk(lp),gk(vp)])

        ! densities
        rhok = (pS + ps_inf)/((gs_min - 1)*cvs*Tk)

        ! internal energy
        ek = (pS + gs_min*ps_inf)/(pS + ps_inf)*cvs*Tk + qvs

        !! printing output !!
        print *, 'cell at which the error happened'

        print *, '[x,y,z]_index', j, k, l

        !! per fluid quantities !!
        print *, 'per fluid quantities'

        print *, 'mk', mk

        print *, 'rhok', rhok

        alphak = mk / rhok

        alphak( pack( iFix, mk < sgm_eps ) ) = 0._wp

        print *, 'alphak', alphak

        if (model_eqns == 3) then
            print *, 'mek', mk * ek
        end if

        print *, 'mqk', mk * qvs

        print *, 'Yk', mk / sum(mk)

        print *, 'Tk', Tk

        print *, 'ps_inf', ps_inf

        !! global quantities !!
        print *, 'global quantities'

        print *, 'mQ', sum(mk * qvs)

        print *, 'pS', pS

        print *, 'rho', sum(mk)

        print *, 'rhoe', rhoe

        if (model_eqns == 3) then
            print *, 'rhoe6E', sum( mk * ek )
            print *, 'abs( rhoe - rhoe6E )', abs( rhoe - sum( mk * ek ) )
        end if

        print *, 'expected TS', (rhoe + pS - sum(mk * qvs)) / sum(mk * cvs * gs_min)

        !! Solver quantities !!
        print *, 'solver quantities'

        print *, 'ns', ns

        print *, 'l2(R2D)', norm2(R2D)

        if ( any( (/ 1, 4 /) == relax_model ) ) then
          print *, 'l2(R2Dr)', norm2(R2D)/abs(rhoe)
          print *, 'Energy constrain', sum( (gs_min - 1.0_wp)*(mk*ek - mk*qvs) / (gs_min*p_inf - minval(p_inf)) )
        else
          print *, 'l2(R2Dr)', norm2(R2D*(/maxg,rhoe/))/norm2((/maxg,rhoe/))
          print *, 'Energy constrain', (rhoe - sum(mk * qvs) - minval(p_inf))
        end if

        print *, 'DeltamP', DeltamP

        print *, '-min(ps_inf)', -minval(p_inf)

        print *, 'J', Jac, 'J-1', InvJac

    end subroutine s_whistleblower

        !>  This auxiliary routine computes the requested saturation property
        !!      from the supplied saturation state and initial guess
        !!  @param pSat Saturation pressure. Input for iSatOut = 1 and output for iSatOut = 2
        !!  @param TSat Saturation temperature. Output for iSatOut = 1 and input for iSatOut = 2
        !!  @param SatIn Initial guess for the Newton solver: temperature when computing TSat, pressure when computing pSat
        !!  @param iSatOut Saturation-property selector chosen by the caller
    subroutine s_Saturation_Properties(pSat, TSat, SatIn, iSatOut)
        $:GPU_ROUTINE(function_name='s_Saturation_Properties',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), intent(inout) :: pSat, TSat
        real(wp), intent(in) :: SatIn
        integer, intent(in) :: iSatOut
        real(wp) :: dFdp, dFdT, FSatProp, Om, pMin
        character(20) :: iSatOutS, nss, SatInS, pSatS, TSatS

        ! Generic loop iterators
        integer :: ns

        ! Shared Newton-solver state. The selected branch below only changes the
        ! residual and jacobian formulas for the requested saturation property.
        Om = under_relax
        ns = 0
        FSatProp = 2.0_wp*ptgalpha_eps

        select case (iSatOut)
        case (1)
            ! Compute saturation temperature from a prescribed saturation pressure.
            ! in case of fluid under tension (p - p_inf > 0, T > 0), or, when subcooled liquid/overheated vapor cannot be
            ! phisically sustained (p = 0, T = 0)
            if ((pSat <= sgm_eps) .and. (SatIn >= 0.0_wp)) then

                ! assigning Saturation temperature
                TSat = 0.0_wp

            else

                ! calculating initial estimate for temperature in the TSat procedure. I will also use this variable to
                ! iterate over the Newton's solver
                TSat = SatIn

                ! Newton solver for finding the saturation temperature as function of pressure. ns == 0, so the loop is
                ! entered at least once.
                do while ( ( abs(FSatProp) > ptgalpha_eps ) .or. ( ns == 0 ) )

                    ! Updating counter for the iterative procedure
                    ns = ns + 1

                    ! residual for the saturation-temperature solve
                    FSatProp = TSat*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                                     *(1 - log(TSat)) - (qvps(lp) - qvps(vp)) &
                                     + cvs(lp)*(gs_min(lp) - 1)*log(pSat + ps_inf(lp)) &
                                     - cvs(vp)*(gs_min(vp) - 1)*log(pSat + ps_inf(vp))) &
                               + qvs(lp) - qvs(vp)

                    ! calculating the jacobian
                    dFdT = &
                        -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*log(TSat) &
                        - (qvps(lp) - qvps(vp)) &
                        + cvs(lp)*(gs_min(lp) - 1)*log(pSat + ps_inf(lp)) &
                        - cvs(vp)*(gs_min(vp) - 1)*log(pSat + ps_inf(vp))

                    ! updating saturation temperature
                    TSat = TSat - Om*FSatProp/dFdT

#ifndef MFC_OpenACC
                    ! Checking if TSat returns a NaN
                    if ((ieee_is_nan(TSat)) .or. (ns > max_iter)) then

                        call s_int_to_str(ns, nss)
                        call s_real_to_str(TSat, TSatS)
                        call s_real_to_str(pSat, pSatS)
                        call s_real_to_str(SatIn, SatInS)
                        call s_mpi_abort('pSat = '//pSatS//', TSat = '//TSatS//', SatIn = '//SatInS//'. &
                        & ns = '//nss//'. m_phase_change, s_Saturation_Properties. Aborting!')

                    end if
#endif
                end do

            end if

        case (2)
            ! Compute saturation pressure from a prescribed saturation temperature.
            ! minimum pressure that keeps the logarithms well-defined
            pMin = maxval((/ -(1.0_wp - ptgalpha_eps)*ps_inf(lp) + ptgalpha_eps, &
                            -(1.0_wp - ptgalpha_eps)*ps_inf(vp) + ptgalpha_eps /))

            ! if the prescribed saturation temperature is nonphysical or
            ! the phase change state cannot be sustained
            if (TSat <= sgm_eps) then

                ! assigning Saturation pressure
                pSat = 0.0_wp

            else

                ! calculating initial estimate for pressure in the pSat procedure. I will also use this variable to
                ! iterate over the Newton's solver
                pSat = max(SatIn, pMin)

                ! Newton solver for finding the saturation pressure as function of temperature. ns == 0, so the loop is
                ! entered at least once.
                do while ( ( abs(FSatProp) > ptgalpha_eps ) .or. ( ns == 0 ) )

                    ! Updating counter for the iterative procedure
                    ns = ns + 1

                    ! residual for the saturation-pressure solve, using the same
                    ! Gibbs-equality expression as in the saturation-temperature path
                    FSatProp = TSat*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                                     *(1 - log(TSat)) - (qvps(lp) - qvps(vp)) &
                                     + cvs(lp)*(gs_min(lp) - 1)*log(pSat + ps_inf(lp)) &
                                     - cvs(vp)*(gs_min(vp) - 1)*log(pSat + ps_inf(vp))) &
                               + qvs(lp) - qvs(vp)

                    ! jacobian of the Gibbs-equality residual with respect to pressure
                    dFdp = TSat*(cvs(lp)*(gs_min(lp) - 1)/(pSat + ps_inf(lp)) &
                                 - cvs(vp)*(gs_min(vp) - 1)/(pSat + ps_inf(vp)))

                    ! updating saturation pressure and keeping the logarithms well-defined
                    pSat = max(pSat - Om*FSatProp/dFdp, pMin)

#ifndef MFC_OpenACC
                    ! Checking if pSat returns a NaN
                    if ((ieee_is_nan(pSat)) .or. (ns > max_iter)) then

                        call s_int_to_str(ns, nss)
                        call s_real_to_str(TSat, TSatS)
                        call s_real_to_str(pSat, pSatS)
                        call s_real_to_str(SatIn, SatInS)
                        call s_mpi_abort('pSat = '//pSatS//', TSat = '//TSatS//', pIn = '//SatInS//'. &
                        & ns = '//nss//'. m_phase_change, s_Saturation_Properties. Aborting!')

                    end if
#endif
                end do

            end if

#ifndef MFC_OpenACC
        case default
            call s_int_to_str(iSatOut, iSatOutS)
            call s_mpi_abort('Unsupported saturation output choice = '//iSatOutS// &
                           & '. Use 1 for outputting temperature or 2 for outputting pressure. ' // &
                           & 'm_phase_change, s_Saturation_Properties. Aborting!')
#endif
        end select

    end subroutine s_Saturation_Properties

    subroutine update_conservative_vars(j, k, l, m0k, pS, q_cons_vf, Tk )

      $:GPU_ROUTINE(function_name='update_conservative_vars',parallelism='[seq]', &
            & cray_inline=True)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp), intent(in) :: pS
        real(wp), dimension(num_fluids), intent(in) :: m0k, Tk
        integer, intent(in) :: j, k, l
        real(wp), dimension(num_fluids) :: sk, hk, gk, ek, rhok
        integer :: i

        ! Thermodynamic state calculations
        ! entropy
        sk = cvs*log((Tk**gs_min)/((pS + ps_inf)**(gs_min - 1.0_wp))) + qvps

        ! enthalpy
        hk = gs_min*cvs*Tk + qvs

        ! Gibbs-free energy
        gk = hk - Tk*sk

        ! densities
        rhok = (pS + ps_inf)/((gs_min - 1)*cvs*Tk)

        ! internal energy
        ek = (pS + gs_min*ps_inf)/(pS + ps_inf)*cvs*Tk + qvs

        ! assigning volume fractions, internal energies, and total entropy
        $:GPU_LOOP(parallelism='[seq]')
        do i = 1, num_fluids

          ! mass factions. Note that, at most, only liquid and vapor masses should change
          q_cons_vf(i + contxb - 1)%sf(j, k, l) = m0k(i)

          ! volume fractions
          q_cons_vf(i + advxb - 1)%sf(j, k, l) = m0k(i)/rhok(i)

          ! alpha*rho*e
          if (model_eqns == 3) then
              q_cons_vf(i + intxb - 1)%sf(j, k, l) = m0k(i)*ek(i)
          end if

        end do
    end subroutine update_conservative_vars

        !>  This auxiliary subroutine enables phase change based on subgrid
        !!  criterium, if subgrid model is activated. This is based on Fuster's
        !!  work (Stability of bubbly liquids and its connection to the process
        !!  of cavitation inception)
    subroutine s_SG_trigger( alpha_b, massIn_b, pS, RIn_b, pVap, TSG )
        $:GPU_ROUTINE(function_name='s_SG_trigger',parallelism='[seq]', &
            & cray_inline=True)

        real(wp), intent(in)    :: alpha_b, massIn_b, pS, RIn_b, pVap
        logical, intent(inout)  :: TSG
        real(wp) :: RBlake

        !! first approximation: dilute limit - Blake's critical radius for
        !! either mono or polydisperse bubbles, since they are into the dilute
        !! limit
        ! RBlake = ( 3.0_wp * gam * R_g * rho0ref / ( 2.0_wp * ss * R0ref ** ( 3.0_wp * gam - 6.0_wp ) ) ) ** ( 1 / ( 5.0_wp - 3.0_wp * gam ) )
        ! below the boiling point, the subgrid model is not activated, as for any T, the equilibrium is stable
        if (pS > pVap) then
            TSG = .false.
        else
            RBlake = 2.0_wp * ss / ( pVap - pS ) * ( 1.0_wp - 1.0_wp / ( 3.0_wp * gam ) )
            TSG = RIn_b > RBlake
        end if

        ! if (TSG) then
        !   print *, 'RIn_b', RIn_b
        !   print *, 'RBlake', RBlake
        !   Print *, '( pVap - pS )', ( pVap - pS )
        !   print *, 'pVap', pVap
        !   print *, 'pS', pS
        ! end if
        ! TSG = alpha_b > 1.0e-4_wp

    end subroutine s_SG_trigger

    impure subroutine s_real_to_str(rl, res)
        real(wp), intent(in) :: rl
        character(len=*), intent(out) :: res
        write (res, '(ES20.8)') rl
        res = trim(res)
    end subroutine s_real_to_str

#endif
end module m_phase_change
