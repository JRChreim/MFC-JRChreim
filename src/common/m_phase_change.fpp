!>
!! @file m_phase_change.fpp
!! @brief Contains module m_phasechange

#:include 'macros.fpp'

!> @brief This module is used to relax the model equations (6-eqn model)
!> towards pressure and temperature (6-eqn to 4-eqn), and (if wanted) Gibbs free
!> energies (6-eqn to 4-eqn) equilibrium through an infinitely fast (algebraic)
!> procedure.
module m_phase_change

#ifndef MFC_POST_PROCESS

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_helper               !< calling module to convert

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use ieee_arithmetic

    ! ==========================================================================

    implicit none

    private; public :: s_initialize_phasechange_module, &
 s_relaxation_solver, &
 s_infinite_relaxation_k, &
 s_finalize_relaxation_solver_module

    !> @name Abstract interface for creating function pointers
    !> @{
    abstract interface

        !> @name Abstract subroutine for the infinite relaxation solver
        !> @{
        subroutine s_abstract_relaxation_solver(q_cons_vf) ! -------
            import :: scalar_field, sys_size
            type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        end subroutine
        !> @}

    end interface
    !> @}

    !> @name Parameters for the first order transition phase change
    !> @{
    integer, parameter :: max_iter = 1d6        !< max # of iterations
    real(kind(0d0)), parameter :: pCr = 4.94d7   !< Critical water pressure
    real(kind(0d0)), parameter :: TCr = 385.05 + 273.15  !< Critical water temperature
    real(kind(0d0)), parameter :: mixM = sgm_eps !< threshold for 'mixture cell'. If Y < mixM, phase change does not happen
    integer, parameter :: lp = 1    !< index for the liquid phase of the reacting fluid
    integer, parameter :: vp = 2    !< index for the vapor phase of the reacting fluid
    !> @}

    !> @name Gibbs free energy phase change parameters
    !> @{
    real(kind(0d0)) :: A, B, C, D
    !> @}

    !$acc declare create(max_iter,pCr,TCr,mixM,lp,vp,A,B,C,D)

    procedure(s_abstract_relaxation_solver), pointer :: s_relaxation_solver => null()

contains

    !>  The purpose of this subroutine is to initialize the phase change module
        !!      by setting the parameters needed for phase change and
        !!      selecting the phase change module that will be used
        !!      (pT- or pTg-equilibrium)
    subroutine s_initialize_phasechange_module()

        ! variables used in the calculation of the saturation curves for fluids 1 and 2
        A = (gs_min(lp)*cvs(lp) - gs_min(vp)*cvs(vp) &
             + qvps(vp) - qvps(lp))/((gs_min(vp) - 1.0d0)*cvs(vp))

        B = (qvs(lp) - qvs(vp))/((gs_min(vp) - 1.0d0)*cvs(vp))

        C = (gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
            /((gs_min(vp) - 1.0d0)*cvs(vp))

        D = ((gs_min(lp) - 1.0d0)*cvs(lp)) &
            /((gs_min(vp) - 1.0d0)*cvs(vp))

        ! Associating procedural pointer to the subroutine that will be
        ! utilized to calculate the solution to the selected relaxation system
        if (any((/4, 5, 6/) == relax_model)) then
            s_relaxation_solver => s_infinite_relaxation_k
        else
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
        real(kind(0.0d0)) :: pS, pSOV, pSSL !< equilibrium pressure for mixture, overheated vapor, and subcooled liquid
        real(kind(0.0d0)) :: TS, TSOV, TSSL, TSatOV, TSatSL !< equilibrium temperature for mixture, overheated vapor, and subcooled liquid. Saturation Temperatures at overheated vapor and subcooled liquid
        real(kind(0.0d0)) :: rhoe, rhoeT, dynE, rhos !< total internal energies (different calculations), kinetic energy, and total entropy
        real(kind(0.0d0)) :: rho, rM, m1, m2 !< total density, total reacting mass, individual reacting masses
        logical :: TR

        !$acc declare create(pS, pSOV, pSSL, TS, TSOV, TSatOV, TSatSL, TSSL, rhoe, rhoeT, dynE, rhos, rho, rM, m1, m2, TR)

        real(kind(0d0)), dimension(num_fluids) :: p_infOV, p_infpT, p_infSL, alpha0k, alphak, m0k, rhok, Tk

        !< Generic loop iterators
        integer :: i, j, k, l, w

        !$acc declare create(p_infOV, p_infpT, p_infSL, alpha0k, alphak, rhok, Tk)

        ! assigning value to the global parameter
        max_iter_pc_ts = 0

        ! starting equilibrium solver
        !$acc parallel loop collapse(3) gang vector default(present) private(p_infOV, p_infpT, p_infSL, alpha0k, alphak, rhok, Tk, pS, pSOV, pSSL, TS, TSOV, TSatOV, TSatSL, TSSL, rhoe, rhoeT, dynE, rhos, rho, rM, m1, m2, TR)
        do j = 0, m
            do k = 0, n
                do l = 0, p
                    ! trigger for phase change. This will be used for checking many conditions through the code
                    TR = .true.

                    ! computing mixture density, volume fraction, and internal energy, so as saving original variables
                    ! in case phase change is cancelled.
                    rho = 0.0d0; rhoeT = 0.0d0
                    !$acc loop seq
                    do i = 1, num_fluids
                        ! original volume fractions, before any relaxation
                        alpha0k(i) = q_cons_vf(i + advxb - 1)%sf(j, k, l)
                        
                        ! original mass fractions, before any relaxation
                        m0k(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)
                              
                        ! Mixture density
                        rho = rho + q_cons_vf(i + contxb - 1)%sf(j, k, l)

                        ! calculating the total internal energy such that the energy-fraction for each of the
                        ! fluids can be proportionally distributed when the sum of the internal energies differs from the
                        if (model_eqns .eq. 3) then
                            rhoeT = rhoeT + q_cons_vf(i + intxb - 1)%sf(j, k, l)
                        end if
                    end do
                    
                    ! kinetic energy as an auxiliary variable to the calculation of the total internal energy
                    dynE = 0.0d0
                    !$acc loop seq
                    do i = momxb, momxe
                        dynE = dynE + 5.0d-1*q_cons_vf(i)%sf(j, k, l)**2/rho
                    end do

                    ! calculating the internal mixture energy that MUST be preserved throughout pT- and pTg-relaxation procedures
                    ! This calulation is performed as the total energy minus the kinetic one as energy it is preserved at discontinuities
                    rhoe = q_cons_vf(E_idx)%sf(j, k, l) - dynE

                    ! calculating the total reacting mass for the phase change process. By hypothesis, this should not change
                    ! throughout the phase-change process.
                    rM = m0k(lp) + m0k(vp)

                    ! correcting possible negative mass fraction values
                    ! at this time, TR is updated if phase change needs to be stopped
                    call s_correct_partial_densities(2, q_cons_vf, rM, rho, TR, i, j, k, l)

                    ! if ( q_cons_vf(vp + advxb - 1)%sf(j, k, l) > 1E-8) then
                    !     PRINT *, 'phase change, m, vf, 176'
                    !     PRINT *, q_cons_vf(lp + contxb - 1)%sf(j, k, l), q_cons_vf(vp + contxb - 1)%sf(j, k, l), q_cons_vf(3 + contxb - 1)%sf(j, k, l)
                    !     PRINT *, q_cons_vf(lp + advxb - 1)%sf(j, k, l), q_cons_vf(vp + advxb - 1)%sf(j, k, l), q_cons_vf(3 + advxb - 1)%sf(j, k, l)
                    ! end if

                    if (TR) then
                        ! p-equilibrium, only
                        if (relax_model == 4) then
                            call s_infinite_p_relaxation_k(j, k, l, pS, q_cons_vf, rho, rhoe, rhoeT, Tk)
                        ! pT-equilibrium only
                        elseif (relax_model == 5) then
                            ! for this case, MFL cannot be either 0 or 1, so I chose it to be 2
                            call s_infinite_pt_relaxation_k(j, k, l, 2, pS, p_infpT, q_cons_vf, rhoe, rM, TS)
                            Tk = spread(TS, 1, num_fluids)
                        ! if pT-pTg equilibrium
                        elseif (relax_model == 6) then

                            call s_infinite_pt_relaxation_k(j, k, l, 2, pS, p_infpT, q_cons_vf, rhoe, rM, TS)
                            Tk = spread(TS, 1, num_fluids)

                            ! if ( q_cons_vf(vp + advxb - 1)%sf(j, k, l) > 1E-8) then
                            !     PRINT *, 'mpp_variables, m, vf, 176'
                            !     PRINT *, q_cons_vf(i + contxb - 1)%sf(j, k, l), q_cons_vf(vp + advxb - 1)%sf(j, k, l)
                            ! end if

                            ! updating the densities and volume fractions used for thresholds
                            rhok = (pS + ps_inf)/((gs_min - 1)*cvs*Tk)
            
                            !$acc loop seq
                            do i = 1, num_fluids
                                ! new volume fractions, after p- or pT-equilibrium
                                alphak(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/rhok(i)
                            end do
                            
                            ! 1 - model activation, 1st order transition (p,T) <= (pCr, TCr)
                            if ( ( pS < pCr ) .and. &
                            ! 2.1 Homogeneous pTg-equilibrium criterium
                            ( ( pS + minval(p_infpT) .gt. 0.0 ) &
                            .or. &
                            ! 2.2. Heterogeneous pTg-equilibrium.
                            ( (alphak(lp) > palpha_eps) .and. (alphak(vp) > palpha_eps) ) ) &
                            ) then

                                ! updating m1 and m2 AFTER correcting the partial densities. These values must be
                                ! stored in case the final state is a mixture of fluids
                                m1 = q_cons_vf(lp + contxb - 1)%sf(j, k, l)
                                m2 = q_cons_vf(vp + contxb - 1)%sf(j, k, l)

                                ! checking if fluid is either subcoooled liquid or overheated vapor (NOT metastability)

                                ! overheated vapor
                                ! depleting the mass of liquid and tranferring the total mass to vapor
                                q_cons_vf(lp + contxb - 1)%sf(j, k, l) = 0*mixM*rM
                                q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0d0 - 0*mixM)*rM

                                ! calling pT-equilibrium for overheated vapor, which is MFL = 0
                                call s_infinite_pt_relaxation_k(j, k, l, 0, pSOV, p_infOV, q_cons_vf, rhoe, rM, TSOV)

                                ! calculating Saturation temperature
                                call s_TSat(pSOV, TSatOV, TSOV)

                                ! subcooled liquid 
                                ! tranferring the total mass to liquid and depleting the mass of vapor
                                q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0d0 - 0*mixM)*rM
                                q_cons_vf(vp + contxb - 1)%sf(j, k, l) = 0*mixM*rM

                                ! calling pT-equilibrium for subcooled liquid, which is MFL = 1                       
                                call s_infinite_pt_relaxation_k(j, k, l, 1, pSSL, p_infSL, q_cons_vf, rhoe, rM, TSSL)

                                ! calculating Saturation temperature
                                call s_TSat(pSSL, TSatSL, TSSL)

                                ! checking the conditions for overheated vapor
                                if (TSOV > TSatOV) then

                                    ! Assigning pressure and temperature
                                    pS = pSOV
                                    TS = TSOV

                                    ! correcting the liquid and vapor partial densities
                                    q_cons_vf(lp + contxb - 1)%sf(j, k, l) = 0*mixM*rM
                                    q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0d0 - 0*mixM)*rM

                                    w = 0

                                ! checking the conditions for subcooled liquid
                                elseif (TSSL < TSatSL) then

                                    ! Assigning pressure and temperature
                                    pS = pSSL
                                    TS = TSSL

                                    ! correcting the liquid and vapor partial densities
                                    q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0d0 - 0*mixM)*rM
                                    q_cons_vf(vp + contxb - 1)%sf(j, k, l) = 0*mixM*rM

                                    w = 1

                                ! if not, allowing phase change (pTg)
                                else
                                    ! returning partial pressures to what they were after the partial density correction 
                                    q_cons_vf(lp + contxb - 1)%sf(j, k, l) = m1
                                    q_cons_vf(vp + contxb - 1)%sf(j, k, l) = m2

                                    ! pTg-equilibrium solver.
                                    call s_infinite_ptg_relaxation_k(j, k, l, pS, p_infpT, rho, rhoe, rM, q_cons_vf, TR, TS)

                                    ! if no pTg happens, the solver will return to the hyperbolic state variables
                                    if ( TR .eqv. .false. ) then
                                        !$acc loop seq
                                        do i = 1, num_fluids
                                            ! returning partial densities to what they were previous to any relaxation scheme
                                            q_cons_vf(i + contxb - 1)%sf(j, k, l) = m0k(i) 
                                        end do 
                                        ! ends the execution of this function and returns control to the calling function
                                        return
                                    end if

                                    w = 2

                                end if
                                Tk = spread(TS, 1, num_fluids)
                            else
                                !$acc loop seq
                                do i = 1, num_fluids
                                    ! returning partial densities to what they were for the hyperbolic solver
                                    q_cons_vf(i + contxb - 1)%sf(j, k, l) = m0k(i)
                                end do
                                ! ends the execution of this function and returns control to the calling function
                                return
                            end if
                        end if
                        
                        ! if ( q_cons_vf(vp + advxb - 1)%sf(j, k, l) > 1E-8) then
                        !     PRINT *, 'phase change, m, vf, 297'
                        !     PRINT *, w
                        !     PRINT *, q_cons_vf(lp + contxb - 1)%sf(j, k, l), q_cons_vf(vp + contxb - 1)%sf(j, k, l), q_cons_vf(3 + contxb - 1)%sf(j, k, l)
                        !     PRINT *, q_cons_vf(lp + advxb - 1)%sf(j, k, l), q_cons_vf(vp + advxb - 1)%sf(j, k, l), q_cons_vf(3 + advxb - 1)%sf(j, k, l)
                        ! end if
                        
                        ! updating conservative variables through either p- or pT-equilibrium
                        call update_conservative_vars( j, k, l, pS, q_cons_vf, Tk )

                        ! if ( q_cons_vf(vp + advxb - 1)%sf(j, k, l) > 1E-8) then
                        !     PRINT *, 'phase change, m, vf, 304'
                        !     PRINT *, w
                        !     PRINT *, q_cons_vf(lp + contxb - 1)%sf(j, k, l), q_cons_vf(vp + contxb - 1)%sf(j, k, l), q_cons_vf(3 + contxb - 1)%sf(j, k, l)
                        !     PRINT *, q_cons_vf(lp + advxb - 1)%sf(j, k, l), q_cons_vf(vp + advxb - 1)%sf(j, k, l), q_cons_vf(3 + advxb - 1)%sf(j, k, l)
                        ! end if

                    else
                        !$acc loop seq
                        do i = 1, num_fluids
                            ! returning partial densities to what they were for the hyperbolic solver
                            q_cons_vf(i + contxb - 1)%sf(j, k, l) = m0k(i) 
                        end do
                    end if
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
    subroutine s_infinite_p_relaxation_k(j, k, l, pS, q_cons_vf, rho, rhoe, rhoeT, Tk)
        !$acc routine seq

        ! initializing variables
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(IN) :: rho, rhoe, rhoeT
        real(kind(0.0d0)), intent(OUT) :: pS
        real(kind(0.0d0)), dimension(num_fluids), intent(OUT) :: Tk
        integer, intent(IN) :: j, k, l

        real(kind(0.0d0)) :: fp, fpp, pO !< variables for the Newton Solver
        real(kind(0.0d0)) :: Econst !< auxiliary variables
        real(kind(0.0d0)), dimension(num_fluids) :: alpha0k, alphak, alpharhoek, alpharhoe0k, rhok
        character(20) :: nss, pSs, Econsts

        integer :: i, ns !< generic loop iterators

        !$acc loop seq
        do i = 1, num_fluids

            ! Re-distributiong internal energy values such that rhoe = rhoeT, eventually.
            ! initial value for internal energy
            if (model_eqns .eq. 3) then
                alpharhoe0k(i) = q_cons_vf(i + intxb - 1)%sf(j, k, l)*rhoe/rhoeT
            else 
                alpharhoe0k(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)*rhoe/rho
            end if

            alphak(i) = q_cons_vf(i + advxb - 1)%sf(j, k, l)

        end do

        ! final value for internal energy, that will be updated with the iterative procedure
        alpharhoek = alpharhoe0k

        ! calculating initial estimate for pressure in the p-relaxation procedure. I will also use this variable to
        ! iterate over the Newton's solver
        pO = 1.0d0

        ! Maybe improve this condition afterwards. As long as the initial guess is in between -min(gs_min*ps_inf)
        ! and infinity, a solution should be able to be found.
        pS = 1.0d3

        ! starting counter for the Newton solver
        ns = 0

        ! Newton solver for p-equilibrium. 1d6 is arbitrary, and ns <= 1, is to ensure the internal energy correction
        ! happens at least once. This is different than the pT-equilibrium case, in which no energy correction is needed.
        do while (((DABS(pS - pO) > ptgalpha_eps) .and. (DABS((pS - pO)/pO) > ptgalpha_eps/1.0d6)) .or. (ns <= 1))

            ! increasing counter
            ns = ns + 1

            ! updating old pressure
            pO = pS

            ! updating functions used in the Newton's solver
            fp = 0.0d0; fpp = 0.0d0; 
            !$acc loop seq
            do i = 1, num_fluids

                fp = fp + (gs_min(i) - 1.0d0)*(alpharhoek(i) - q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)) &
                     /(pO + gs_min(i)*ps_inf(i))

                fpp = fpp - (gs_min(i) - 1.0d0)*(alpharhoek(i) - q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)) &
                      /((pO + gs_min(i)*ps_inf(i))**2)

            end do

            ! updating the relaxed pressure (solution)
            pS = pO + ((1.0d0 - fp)/fpp)/(1.0d0 - (1.0d0 - fp + DABS(1.0d0 - fp))/(2.0d0*fpp*(pO + minval(gs_min*ps_inf))))

            ! Variable to check energy constraint
            Econst = 0.0d0
            !$acc loop seq
            ! updating fluid variables, together with the relaxed pressure, in a loosely coupled procedure
            do i = 1, num_fluids

                ! volume fractions
                alpha0k(i) = alphak(i)

                alphak(i) = (gs_min(i) - 1.0d0)*(alpharhoek(i) - q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)) &
                            /(pS + gs_min(i)*ps_inf(i))

                ! internal energies
                ! alpharhoek(i) = alpharhoek(i) + 1 / (gs_min(i) - 1) &
                ! * (alphak(i)*(pS + gs_min(i)*ps_inf(i)) - alpha0k(i)*(pO + gs_min(i)*ps_inf(i)))

                alpharhoek(i) = alpharhoe0k(i) - pS*(alphak(i) - q_cons_vf(i + advxb - 1)%sf(j, k, l))

                ! densities
                rhok(i) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/alphak(i)

                ! Variable to check energy constraint
                Econst = Econst + (gs_min(i) - 1.0d0)*(alpharhoek(i) - q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)) &
                         /(gs_min(i)*ps_inf(i) - minval(ps_inf))

            end do

#ifndef MFC_OpenACC
            ! energy constraint for the p-equilibrium
            if ((minval(ps_inf) > 0) .and. (Econst <= 1.0d0)) then

                call s_real_to_str(Econst, Econsts)
                call s_mpi_abort('Solver for the p-relaxation solver failed (m_phase_change, s_infinite_p_relaxation_k) &
&                   . Please, check energy constraint. Econst ~'//Econsts//'. Aborting!')

                ! checking if pressure is within expected bounds
            elseif ((pS <= -1.0d0*minval(gs_min*ps_inf)) .or. (ieee_is_nan(pS)) .or. (ns > max_iter)) then

                if (proc_rank == 0) then

                    print *, 'ns', ns

                    print *, 'energies (N,O)', rhoe, rhoeT, 'DE', rhoe - rhoeT

                    print *, 'fp', 'fpp', fp, fpp

                    call s_tattletale((/0.0d0, 0.0d0/), reshape((/0.0d0, 0.0d0, 0.0d0, 0.0d0/), (/2, 2/)) &
                                      , j, (/0.0d0, 0.0d0, 0.0d0, 0.0d0/), k, l, rhoeT, ps_inf, pS, (/pS - pO, pS + pO/) &
                                      , rhoe, q_cons_vf, SUM(Tk)/num_fluids)
                end if

                call s_real_to_str(pS, pSs)
                call s_int_to_str(ns, nss)
                call s_mpi_abort('Solver for the p-relaxation failed (m_phase_change, s_infinite_p_relaxation_k). &
                &   pS ~'//pSs//'. ns = '//nss//'. Aborting!')

            end if
#endif
        end do

        ! (NOT common) temperature
        Tk = (pS + ps_inf)/((gs_min - 1)*cvs*rhok)

        ! updating maximum number of iterations
        max_iter_pc_ts = maxval((/max_iter_pc_ts, ns/))
    end subroutine s_infinite_p_relaxation_k ! -----------------------

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
    subroutine s_infinite_pt_relaxation_k(j, k, l, MFL, pS, p_infpT, q_cons_vf, rhoe, rM, TS)
        !$acc routine seq

        ! initializing variables
        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(OUT) :: pS, TS
        real(kind(0.0d0)), dimension(num_fluids), intent(OUT) :: p_infpT
        real(kind(0.0d0)), intent(IN) :: rhoe, rM
        integer, intent(IN) :: j, k, l, MFL
        integer, dimension(num_fluids) :: ig !< flags to toggle the inclusion of fluids for the pT-equilibrium
        real(kind(0.0d0)) :: gp, gpp, hp, pO, mCP, mQ !< variables for the Newton Solver
        character(20) :: nss, pSs, Econsts

        integer :: i, ns !< generic loop iterators

        ! auxiliary variables for the pT-equilibrium solver
        mCP = 0.0d0; mQ = 0.0d0; p_infpT = ps_inf;

        ! these are slowing the computations significantly. Think about a workaround
        ig(1:num_fluids) = 0

        ! Performing tests before initializing the pT-equilibrium
        !$acc loop seq
        do i = 1, num_fluids

            ! check if all alpha(i)*rho(i) are negative. If so, abort
! #ifndef MFC_OpenACC
            ! check which indices I will ignore (no need to abort the solver in this case). Adjust this sgm_eps value for mixture cells
            if( ( q_cons_vf( i + contxb - 1 )%sf( j, k, l ) .ge. 0.0D0 ) &
                    .and. ( q_cons_vf( i + contxb - 1 )%sf( j, k, l ) - rM * mixM .le. sgm_eps ) ) then

                ig(i) = i

                ! this value is rather arbitrary, as I am interested in MINVAL( ps_inf ) for the solver.
                ! This way, I am ensuring this value will not be selected.
                p_infpT(i) = 2 * MAXVAL( ps_inf )

            end if
! #endif

            ! sum of the total alpha*rho*cp of the system
            mCP = mCP + q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i)*gs_min(i)

            ! sum of the total alpha*rho*q of the system
            mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)
        end do

        ! Checking energy constraint. In the case we are calculating the possibility of having subcooled liquid or
        ! overheated vapor, the energy constraint might not be satisfied, as are hypothetically transferring all the 
        ! mass from one phase to the other. When this is the case, we simply ignore this possibility, set pS = TS = 0,
        ! and discard the hypothesis. The solver can thus move forward.
        if ((rhoe - mQ - minval(p_infpT)) < 0.0d0) then

            if ((MFL == 0) .or. (MFL == 1)) then

                ! Assigning zero values for mass depletion cases
                ! pressure
                pS = 0.0d0

                ! temperature
                TS = 0.0d0

                return
#ifndef MFC_OpenACC
            else
                if (proc_rank == 0) then

                    call s_tattletale((/0.0d0, 0.0d0/), reshape((/0.0d0, 0.0d0, 0.0d0, 0.0d0/), (/2, 2/)) &
                                      , j, (/0.0d0, 0.0d0, 0.0d0, 0.0d0/), k, l, mQ, p_infpT, pS, (/DABS(pS - pO), DABS(pS - pO)/) &
                                      , rhoe, q_cons_vf, TS)

                end if

                call s_real_to_str(rhoe - mQ - minval(p_infpT), Econsts)
                call s_mpi_abort('Solver for the pT-relaxation solver failed (m_phase_change, s_infinite_pt_relaxation_k) &
                &    . Energy constraint~'//Econsts//'. Aborting!')

#endif
            end if
        end if

        ! calculating initial estimate for pressure in the pT-relaxation procedure. I will also use this variable to
        ! iterate over the Newton's solver
        pO = 0.0d0

        ! Maybe improve this condition afterwards. As long as the initial guess is in between -min(ps_inf)
        ! and infinity, a solution should be able to be found.
        pS = 1.0d4

        ! starting counter for the Newton solver
        ns = 0

        ! Newton solver for pT-equilibrium. 1d6 is arbitrary, and ns == 0, to the loop is entered at least once.
        do while (((DABS(pS - pO) > ptgalpha_eps) .and. (DABS((pS - pO)/pO) > ptgalpha_eps/1.0d6)) .or. (ns == 0))

            ! increasing counter
            ns = ns + 1

            ! updating old pressure
            pO = pS

            ! updating functions used in the Newton's solver
            gpp = 0.0d0; gp = 0.0d0; hp = 0.0d0
            !$acc loop seq
            do i = 1, num_fluids

                ! given pS always change, I need ig( i ) and gp to be in here, as it dynamically updates.
                ! Note that I do not need to use p_infpT here, but I will do it for consistency
                if (i /= ig(i)) then

                    gp = gp + (gs_min(i) - 1.0d0)*q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                        *(rhoe + pO - mQ)/(mCP*(pO + p_infpT(i)))

                    gpp = gpp + (gs_min(i) - 1.0d0)*q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                        *(p_infpT(i) - rhoe + mQ)/(mCP*(pO + p_infpT(i))**2)

                end if

            end do

            hp = 1.0d0/(rhoe + pO - mQ) + 1.0d0/(pO + minval(p_infpT))

            ! updating common pressure for the newton solver
            pS = pO + ((1.0d0 - gp)/gpp)/(1.0d0 - (1.0d0 - gp + DABS(1.0d0 - gp)) &
                                          /(2.0d0*gpp)*hp)

            ! check if solution is out of bounds (which I believe it won`t happen given the solver is gloabally convergent.
#ifndef MFC_OpenACC
            if ((pS <= -1.0d0*minval(p_infpT)) .or. (ieee_is_nan(pS)) .or. (ns > max_iter)) then
                if (proc_rank == 0) then
                    call s_tattletale((/0.0d0, 0.0d0/), reshape((/0.0d0, 0.0d0, 0.0d0, 0.0d0/), (/2, 2/)) &
                                      , j, (/0.0d0, 0.0d0, 0.0d0, 0.0d0/), k, l, mQ, p_infpT, pS, (/pS - pO, pS + pO/) &
                                      , rhoe, q_cons_vf, TS)
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
    subroutine s_infinite_ptg_relaxation_k(j, k, l, pS, p_infpT, rho, rhoe, rM, q_cons_vf, TR, TS)

        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: p_infpT
        real(kind(0.0d0)), intent(INOUT) :: pS, TS, rM
        real(kind(0.0d0)), intent(IN) :: rho, rhoe
        integer, intent(IN) :: j, k, l
        logical, intent(INOUT) :: TR
        real(kind(0.0d0)), dimension(num_fluids) :: p_infpTg
        real(kind(0.0d0)), dimension(2, 2) :: Jac, InvJac, TJac
        real(kind(0.0d0)), dimension(2) :: R2D, DeltamP
        real(kind(0.0d0)), dimension(3) :: Oc
        real(kind(0.0d0)) :: Om, OmI ! underrelaxation factor
        real(kind(0.0d0)) :: mCP, mCPD, mCVGP, mCVGP2, mQ, mQD ! auxiliary variables for the pTg-solver
        character(20) :: nss, pSs, Econsts

        !< Generic loop iterators
        integer :: i, ns

        ! checking if homogeneous cavitation is expected. If yes, transfering an amount of mass to the depleted (liquid,
        ! for the moment) phase, and then let the algorithm run. 
        ! checking if homogeneous cavitation is possible
        
        ! PRINT *, 'pc happening'

        ! is the fluid at a metastable state with enough 'energy' for phase change to happen?
        if ((pS < -6.0d5) .and. (q_cons_vf(lp + contxb - 1)%sf(j, k, l) + q_cons_vf(vp + contxb - 1)%sf(j, k, l) &
                                    > (rhoe - gs_min(lp)*ps_inf(lp)/(gs_min(lp) - 1))/qvs(lp))) then

            ! transfer a bit of mass to the deficient phase, enforce phase0chane
            call s_correct_partial_densities(1, q_cons_vf, rM, rho, TR, i, j, k, l)

            PRINT *, 'is it pS < 7.5E-5'
            
        ! the metastable state is not enough to sustain phase change
        elseif (pS < 0.0d0) then
            
            PRINT *, 'is it pS < 0'
            PRINT *, 'pS', pS
            ! cancel any phase-change updates.
            TR = .false.

            ! ends the execution of this function and returns control to the calling function
            return

        ! if not homogeneous, then heterogeneous. Thus, setting up an arbitrary initial condition in case the one from
        ! the p(T)-equilibrium solver could lead to numerical issues
        elseif ((pS < 1.0d-1) .and. (pS >= 0.0d0)) then
            ! improve this initial condition
            pS = 1.0d4

            ! PRINT *, 'is it pS > 0'
        end if

        ! PRINT *, 'phase change is happening'
        ! PRINT *, q_cons_vf(lp + contxb - 1)%sf(j, k, l), q_cons_vf(vp + contxb - 1)%sf(j, k, l), q_cons_vf(3 + contxb - 1)%sf(j, k, l)
        ! if ( q_cons_vf(lp + advxb - 1)%sf(j, k, l) .le. 1.0d-8 .or. q_cons_vf(vp + advxb - 1)%sf(j, k, l) .le. 1.0d-8 ) then
            ! PRINT *, 'phase change is happening, mass, vf'
            ! PRINT *, q_cons_vf(lp + contxb - 1)%sf(j, k, l), q_cons_vf(vp + contxb - 1)%sf(j, k, l), q_cons_vf(3 + contxb - 1)%sf(j, k, l)
            ! PRINT *, q_cons_vf(lp + advxb - 1)%sf(j, k, l), q_cons_vf(vp + advxb - 1)%sf(j, k, l), q_cons_vf(3 + advxb - 1)%sf(j, k, l)
        ! end if

        ! Relaxation factor. This value is rather arbitrary, with a certain level of self adjustment.
        OmI = 1.0d-1
        ! Critical relaxation factors, for variable sub-relaxation
        Oc(1) = OmI; Oc(2) = OmI; Oc(3) = OmI

        R2D(1) = 0.0d0; R2D(2) = 0.0d0
        DeltamP(1) = 0.0d0; DeltamP(2) = 0.0d0
        ! starting counter for the Newton solver
        ns = 0

        ! Newton solver for pTg-equilibrium. 1d6 is arbitrary, and ns == 0, to the loop is entered at least once.
        do while (((DSQRT(R2D(1)**2 + R2D(2)**2) > ptgalpha_eps) &
                   .and. ((DSQRT(R2D(1)**2 + R2D(2)**2)/rhoe) > (ptgalpha_eps/1d6))) &
                  .or. (ns == 0))

            ! Updating counter for the iterative procedure
            ns = ns + 1

            ! Auxiliary variables to help in the calculation of the residue
            mCP = 0.0d0; mCPD = 0.0d0; mCVGP = 0.0d0; mCVGP2 = 0.0d0; mQ = 0.0d0; mQD = 0.0d0
            ! Those must be updated through the iterations, as they either depend on
            ! the partial masses for all fluids, or on the equilibrium pressure
            !$acc loop seq
            do i = 1, num_fluids
                ! sum of the total alpha*rho*cp of the system
                mCP = mCP + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                      *cvs(i)*gs_min(i)

                ! sum of the total alpha*rho*q of the system
                mQ = mQ + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

                ! These auxiliary variables now need to be updated, as the partial densities now
                ! vary at every iteration.
                if ((i /= lp) .and. (i /= vp)) then

                    mCVGP = mCVGP + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                            *cvs(i)*(gs_min(i) - 1)/(pS + ps_inf(i))

                    mCVGP2 = mCVGP2 + q_cons_vf(i + contxb - 1)%sf(j, k, l) &
                             *cvs(i)*(gs_min(i) - 1)/((pS + ps_inf(i))**2)

                    mQD = mQD + q_cons_vf(i + contxb - 1)%sf(j, k, l)*qvs(i)

                    ! sum of the total alpha*rho*cp of the system
                    mCPD = mCPD + q_cons_vf(i + contxb - 1)%sf(j, k, l)*cvs(i) &
                           *gs_min(i)

                end if
            end do

            ! Checking pressure and energy criteria for the (pT) solver to find a solution
#ifndef MFC_OpenACC
            if ((pS <= -1.0d0*minval(ps_inf)) .or. ((rhoe - mQ - minval(ps_inf)) < 0.0d0)) then
                if (proc_rank == 0) then
                    call s_tattletale(DeltamP, InvJac, j, Jac, k, l, mQ, ps_inf, pS &
                                      , R2D, rhoe, q_cons_vf, TS)
                end if

                call s_real_to_str(rhoe - mQ - minval(ps_inf), Econsts)
                call s_real_to_str(pS, pSs)
                call s_mpi_abort('Solver for the pTg-relaxation failed (m_phase_change, s_infinite_ptg_relaxation_k). &
                &   pS ~'//pSs//'. Econst = '//Econsts//'. Aborting!')

            end if
#endif
            ! calculating the (2D) Jacobian Matrix used in the solution of the pTg-quilibrium model
            call s_compute_jacobian_matrix(InvJac, j, Jac, k, l, mCPD, mCVGP, mCVGP2, pS, q_cons_vf, TJac)

            ! calculating correction array for Newton's method
            DeltamP = -1.0d0*matmul(InvJac, R2D)

            ! checking if the correction in the mass/pressure will lead to negative values for those quantities
            ! If so, adjust the underrelaxation parameter Om

#ifndef MFC_OpenACC
            ! creating criteria for variable underrelaxation factor
            if (q_cons_vf(lp + contxb - 1)%sf(j, k, l) + Om*DeltamP(1) < 0.0d0) then
                Oc(1) = -q_cons_vf(lp + contxb - 1)%sf(j, k, l)/(2*DeltamP(1))
            else
                Oc(1) = OmI
            end if
            if (q_cons_vf(vp + contxb - 1)%sf(j, k, l) - Om*DeltamP(1) < 0.0d0) then
                Oc(2) = q_cons_vf(vp + contxb - 1)%sf(j, k, l)/(2*DeltamP(1))
            else
                Oc(2) = OmI
            end if
            if (pS + 1.0d0*minval(ps_inf) + Om*DeltamP(2) < 0.0d0) then
                Oc(3) = (pS - 1.0d0*minval(ps_inf))/(2*DeltamP(2))
            else
                Oc(3) = OmI
            end if
            ! choosing amonst the minimum relaxation maximum to ensure solver will not produce unphysical values
            Om = minval(Oc)
#else
            Om = 1.0d-1
#endif
            ! updating two reacting 'masses'. Recall that inert 'masses' do not change during the phase change
            ! liquid
            q_cons_vf(lp + contxb - 1)%sf(j, k, l) = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + Om*DeltamP(1)

            ! gas
            q_cons_vf(vp + contxb - 1)%sf(j, k, l) = q_cons_vf(vp + contxb - 1)%sf(j, k, l) - Om*DeltamP(1)

            ! updating pressure
            pS = pS + Om*DeltamP(2)

            ! calculating residuals, which are (i) the difference between the Gibbs Free energy of the gas and the liquid
            ! and (ii) the energy before and after the phase-change process.
            call s_compute_pTg_residue(j, k, l, mCPD, mCVGP, mQD, q_cons_vf, pS, rhoe, R2D)

            ! checking if the residue returned any NaN values
#ifndef MFC_OpenACC
            if ((ieee_is_nan(R2D(1))) .or. (ieee_is_nan(R2D(2))) .or. (ns > max_iter)) then

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
    subroutine s_correct_partial_densities(CT, q_cons_vf, rM, rho, TR, i, j, k, l)
        !$acc routine seq

        !> @name variables for the correction of the reacting partial densities
        !> @{
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), intent(INOUT) :: rM
        real(kind(0.0d0)), intent(IN) :: rho
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
                if (q_cons_vf(i + contxb - 1)%sf(j, k, l)/rho < 0) then
                    ! do not continue relaxation
                    TR = .false.
                end if
            end do
        elseif (CT == 1) then
            if (rM < 0.0d0) then
                ! reacting masses are very negative so as to affect the physics of the problem, so phase change will not be activated
                if ((q_cons_vf(lp + contxb - 1)%sf(j, k, l)/rM < mixM) .or. &
                    (q_cons_vf(vp + contxb - 1)%sf(j, k, l)/rM < mixM)) then
                    
                    ! do not continue relaxation
                    TR = .false.
                ! reacting masses are not as negative so I can disregard them,
                ! expecting no significant changes in the physics of the simulation
                else
                    q_cons_vf(lp + contxb - 1)%sf(j, k, l) = mixM*rM
                    q_cons_vf(vp + contxb - 1)%sf(j, k, l) = mixM*rM

                    rM = q_cons_vf(lp + contxb - 1)%sf(j, k, l) + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

                    ! continue relaxation
                    TR = .true.
                end if
            ! correcting the partial densities of the reacting fluids. In case liquid is negative
            elseif (q_cons_vf(lp + contxb - 1)%sf(j, k, l)/rM < mixM) then

                q_cons_vf(lp + contxb - 1)%sf(j, k, l) = mixM*rM
                q_cons_vf(vp + contxb - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM
                
                ! continue relaxation
                TR = .true.

            ! correcting the partial densities of the reacting fluids. In case vapor is negative
            elseif (q_cons_vf(vp + contxb - 1)%sf(j, k, l)/rM < mixM) then

                q_cons_vf(lp + contxb - 1)%sf(j, k, l) = (1.0d0 - mixM)*rM
                q_cons_vf(vp + contxb - 1)%sf(j, k, l) = mixM*rM

                ! continue relaxation
                TR = .true.

            end if
        elseif (CT == 2) then
            do i = 1, num_fluids
                if (q_cons_vf(i + advxb - 1)%sf(j, k, l) < 0 .or. &
                    q_cons_vf(i + contxb - 1)%sf(j, k, l)  < 0 ) then
                    q_cons_vf(i + advxb - 1)%sf(j, k, l) = 0.0
                    q_cons_vf(i + contxb - 1)%sf(j, k, l) = 0.0
                    if (model_eqns .eq. 3) then
                        q_cons_vf(i + intxb - 1)%sf(j, k, l) = 0.0
                    end if 
                end if
            end do
            ! continue relaxation
            TR = .true.
        else
            !$acc loop seq
            do i = 1, num_fluids
                if ((q_cons_vf(i + contxb - 1)%sf(j, k, l)/rho) < mixM) then
                    q_cons_vf(i + contxb - 1)%sf(j, k, l) = mixM*rho
                end if
            end do
            ! continue relaxation
            TR = .true.
        end if
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
    subroutine s_compute_jacobian_matrix(InvJac, j, Jac, k, l, mCPD, mCVGP, mCVGP2, pS, q_cons_vf, TJac)
        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(IN) :: pS, mCPD, mCVGP, mCVGP2
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(2, 2), intent(OUT) :: Jac, InvJac, TJac
        real(kind(0.0d0)) :: ml, mT, TS, dFdT, dTdm, dTdp ! mass of the reacting fluid, total reacting mass, and auxiliary variables

        ! mass of the reacting liquid
        ml = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

        ! mass of the two participating fluids
        mT = q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
             + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

        TS = 1/(mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) &
                + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                + mCVGP)

        dFdT = &
            -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*DLOG(TS) &
            - (qvps(lp) - qvps(vp)) &
            + cvs(lp)*(gs_min(lp) - 1)*DLOG(pS + ps_inf(lp)) &
            - cvs(vp)*(gs_min(vp) - 1)*DLOG(pS + ps_inf(vp))

        dTdm = -(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                 - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)))*TS**2

        dTdp = (mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2 &
                + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp))**2 &
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
                     /(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                       + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP) &
                     - (ml*(cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                        - mT*cvs(vp)*gs_min(vp) - mCPD) &
                     *(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                       - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                     /((ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                            - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                        + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP)**2))/1
        ! dF2dp
        Jac(2, 2) = (1 + (ml*(cvs(vp)*gs_min(vp) - cvs(lp)*gs_min(lp)) &
                          - mT*cvs(vp)*gs_min(vp) - mCPD) &
                     *(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp))**2 &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2) &
                       + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))**2 + mCVGP2) &
                     /(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                       + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP)**2)/1

        ! intermediate elements of J^{-1}
        InvJac(1, 1) = Jac(2, 2)
        InvJac(1, 2) = -1.0d0*Jac(1, 2)
        InvJac(2, 1) = -1.0d0*Jac(2, 1)
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
    subroutine s_compute_pTg_residue(j, k, l, mCPD, mCVGP, mQD, q_cons_vf, pS, rhoe, R2D)
        !$acc routine seq

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), intent(IN) :: pS, rhoe, mCPD, mCVGP, mQD
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)), dimension(2), intent(OUT) :: R2D
        real(kind(0.0d0)) :: ml, mT, TS !< mass of the reacting liquid, total reacting mass, equilibrium temperature

        ! mass of the reacting liquid
        ml = q_cons_vf(lp + contxb - 1)%sf(j, k, l)

        ! mass of the two participating fluids
        mT = q_cons_vf(lp + contxb - 1)%sf(j, k, l) &
             + q_cons_vf(vp + contxb - 1)%sf(j, k, l)

        TS = 1/(mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) &
                + ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                      - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                + mCVGP)

        ! Gibbs Free Energy Equality condition (DG)
        R2D(1) = TS*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                     *(1 - DLOG(TS)) - (qvps(lp) - qvps(vp)) &
                     + cvs(lp)*(gs_min(lp) - 1)*DLOG(pS + ps_inf(lp)) &
                     - cvs(vp)*(gs_min(vp) - 1)*DLOG(pS + ps_inf(vp))) &
                 + qvs(lp) - qvs(vp)

        ! Constant Energy Process condition (DE)
        R2D(2) = (rhoe + pS &
                  + ml*(qvs(vp) - qvs(lp)) - mT*qvs(vp) - mQD &
                  + (ml*(gs_min(vp)*cvs(vp) - gs_min(lp)*cvs(lp)) &
                     - mT*gs_min(vp)*cvs(vp) - mCPD) &
                  /(ml*(cvs(lp)*(gs_min(lp) - 1)/(pS + ps_inf(lp)) &
                        - cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp))) &
                    + mT*cvs(vp)*(gs_min(vp) - 1)/(pS + ps_inf(vp)) + mCVGP))/1

    end subroutine s_compute_pTg_residue

    ! SUBROUTINE CREATED TO TELL ME WHERE THE ERROR IN THE PT- AND PTG-EQUILIBRIUM SOLVERS IS
    subroutine s_tattletale(DeltamP, InvJac, j, Jac, k, l, mQ, p_infA, pS, R2D, rhoe, q_cons_vf, TS) ! ----------------

        type(scalar_field), dimension(sys_size), intent(IN) :: q_cons_vf
        real(kind(0.0d0)), dimension(2, 2), intent(IN) :: Jac, InvJac
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: p_infA
        real(kind(0.0d0)), dimension(2), intent(IN) :: R2D, DeltamP
        real(kind(0.0d0)), intent(IN) :: pS, TS
        real(kind(0.0d0)), intent(IN) :: rhoe, mQ
        integer, intent(IN) :: j, k, l
        real(kind(0.0d0)) :: rho
        !< Generic loop iterator
        integer :: i

        print *, 'j, k, l', j, k, l

        print *, 'rhoe', rhoe

        print *, 'mQ', mQ

        print *, 'Energy constrain', (rhoe - mQ - minval(p_infA))

        print *, 'R2D', R2D

        print *, 'l2(R2D)', DSQRT(R2D(1)**2 + R2D(2)**2)

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

            if (model_eqns .eq. 3) then
                print *, 'internal energies', q_cons_vf(i + intxb - 1)%sf(j, k, l)
            end if

            print *, 'Y_i', q_cons_vf(i + contxb - 1)%sf(j, k, l)/rho

        end do

        print *, 'J', Jac, 'J-1', InvJac

    end subroutine s_tattletale

    ! Newton Solver for the finding the Saturation temperature TSat for a given saturation pressure
    subroutine s_TSat(pSat, TSat, TSIn)
        !$acc routine seq

        real(kind(0.0d0)), intent(OUT) :: TSat
        real(kind(0.0d0)), intent(IN) :: pSat, TSIn
        real(kind(0.0d0)) :: dFdT, FT, Om !< auxiliary variables
        character(20) :: nss, pSatS, TSatS

        ! Generic loop iterators
        integer :: ns

        ! in case of fluid under tension (p - p_inf > 0, T > 0), or, when subcooled liquid/overheated vapor cannot be
        ! phisically sustained (p = 0, T = 0)
        if ((pSat .le. 0.0d0) .and. (TSIn .ge. 0.0d0)) then

            ! assigning Saturation temperature
            TSat = 0.0d0

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
            do while ((DABS(FT) > ptgalpha_eps) .or. (ns == 0))

                ! Updating counter for the iterative procedure
                ns = ns + 1

                ! calculating residual
                ! FT = A + B / TSat + C * DLOG( TSat ) + D * DLOG( ( pSat + ps_inf( lp ) ) ) - DLOG( pSat + ps_inf( vp ) )

                FT = TSat*((cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp)) &
                           *(1 - DLOG(TSat)) - (qvps(lp) - qvps(vp)) &
                           + cvs(lp)*(gs_min(lp) - 1)*DLOG(pSat + ps_inf(lp)) &
                           - cvs(vp)*(gs_min(vp) - 1)*DLOG(pSat + ps_inf(vp))) &
                     + qvs(lp) - qvs(vp)

                ! calculating the jacobian
                ! dFdT = - B / ( TSat ** 2) + C / TSat

                dFdT = &
                    -(cvs(lp)*gs_min(lp) - cvs(vp)*gs_min(vp))*DLOG(TSat) &
                    - (qvps(lp) - qvps(vp)) &
                    + cvs(lp)*(gs_min(lp) - 1)*DLOG(pSat + ps_inf(lp)) &
                    - cvs(vp)*(gs_min(vp) - 1)*DLOG(pSat + ps_inf(vp))

                ! updating saturation temperature
                TSat = TSat - Om*FT/dFdT

#ifndef MFC_OpenACC
                ! Checking if TSat returns a NaN
                if ((ieee_is_nan(TSat)) .or. (ns > max_iter)) then

                    call s_int_to_str(ns, nss)
                    call s_real_to_str(TSat, TSatS)
                    call s_real_to_str(pSat, pSatS)
                    PRINT *, pSat
                    call s_mpi_abort('TSat = '//TSatS//', pSat = '// pSatS //' (by assumption of first order transition). &
&                     ns = '//nss//'. m_phase_change, s_TSat. Aborting!')

                end if
#endif
            end do

        end if

    end subroutine s_TSat

    subroutine update_conservative_vars( j, k, l, pS, q_cons_vf, Tk )
        
        !$acc routine seq
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        real(kind(0.0d0)), intent(IN) :: pS
        real(kind(0.0d0)), dimension(num_fluids), intent(IN) :: Tk
        integer, intent(IN) :: j, k, l
        real(kind(0d0)), dimension(num_fluids) :: sk, hk, gk, ek, rhok
        integer :: i

        ! Calculations AFTER then intended equilibrium
        ! entropy
        sk = cvs*DLOG((Tk**gs_min)/((pS + ps_inf)**(gs_min - 1.0d0))) + qvps

        ! enthalpy
        hk = gs_min*cvs*Tk + qvs

        ! Gibbs-free energy
        gk = hk - Tk*sk

        ! densities
        rhok = (pS + ps_inf)/((gs_min - 1)*cvs*Tk)

        ! internal energy
        ek = (pS + gs_min*ps_inf)/(pS + ps_inf)*cvs*Tk + qvs

        ! calculating volume fractions, internal energies, and total entropy
        ! rhos = 0.0d0
        !$acc loop seq
        do i = 1, num_fluids

            ! volume fractions
            q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)/rhok(i)

            ! if ( q_cons_vf(vp + advxb - 1)%sf(j, k, l) > 1.0d-8) then
            !     PRINT *, 'inside TR 1'
            !     PRINT *, q_cons_vf(vp + advxb - 1)%sf(j, k, l)
            !     PRINT *, 'inside TR 2'
            ! end if

            ! alpha*rho*e
            if (model_eqns .eq. 3) then
                q_cons_vf(i + intxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)*ek(i)
            end if

            ! Total entropy
            ! rhos = rhos + q_cons_vf(i + contxb - 1)%sf(j, k, l)*sk(i)

        end do
        
    end subroutine update_conservative_vars

    subroutine s_real_to_str(rl, res)
        character(len=*) :: res
        real(kind(0.0d0)), intent(IN) :: rl
        write (res, '(F10.4)') rl
        res = trim(res)
    end subroutine s_real_to_str

    !>  This subroutine finalizes the phase change module
    subroutine s_finalize_relaxation_solver_module()

        s_relaxation_solver => null()

    end subroutine s_finalize_relaxation_solver_module

#endif
end module m_phase_change