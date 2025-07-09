!> this module (de-)normalize the input variables such that the solver can be
!> solve equations normalized while the inputs and outputs are not.

#:include 'macros.fpp'

module m_normalize

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
    integer, parameter :: max_iter = 1e8_wp          !< max # of iterations
    real(wp), parameter :: p0 = 101325.0_wp             !< Reference Pressure in Pa
    real(wp), parameter :: rho0 = 1.0e3_wp           !< Reference density in kg / m3
    real(wp), parameter :: mixM = 1.0e-8_wp !< threshold for 'mixture cell'. If Y < mixM, phase change does not happen
    integer, parameter :: lp = 1    !< index for the liquid phase of the reacting fluid
    integer, parameter :: vp = 2    !< index for the vapor phase of the reacting fluid
    !> @}

    $:GPU_DECLARE(create='[pCr,TCr,mixM,lp,vp,A,B,C,D]')

contains

    !>  The purpose of this subroutine is to normalize and denormalize the code 
        !!      such that phase change bubble models can be merged easily.
        !!      This function will be used at the beginning of each time step
        !!      (normalize) and after m_phase_change.fpp (denormalize) module
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param ND string to toogle Normalize and Denormalize inputs
    subroutine s_norm_denorm(ND, q_cons_vf)
                
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        character(len=1), intent(in) :: ND

        select case(ND)
            case ("N")

            do j = 0, m
                do k = 0, n
                    do l = 0, p
                        do i = 1, num_fluids
                        ! mass factions
                        q_cons_vf(i + contxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l)

                        ! volume fractions
                        q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + advxb - 1)%sf(j, k, l)

                        ! mixture energies
                        if (model_eqns == 3) then
                            q_cons_vf(i + intxb - 1)%sf(j, k, l) = q_cons_vf(i + intxb - 1)%sf(j, k, l)
                        end if

                        ! total energy
                        q_cons_vf(E_idx)%sf(j, k, l) = q_cons_vf(E_idx)%sf(j, k, l)

                        end do

                        $:GPU_LOOP(parallelism='[seq]')
                        ! momentum
                        do i = momxb, momxe
                            q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do

            case ("D")

        end select

    end subroutine s_norm_denorm !-------------------------------

end module m_normalize