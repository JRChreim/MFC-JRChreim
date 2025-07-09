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

    private; public :: s_norm_denorm
    
    !> @name Parameters for the first order transition phase change
    !> @{
    real(wp), parameter :: p0 = 101325.0_wp             !< Reference Pressure in Pa
    real(wp), parameter :: rho0 = 1.0e3_wp              !< Reference density  in kg / m3
    real(wp), parameter :: u0 = 1.0e-8_wp               !< Reference velocity in m  / s
    !> @}

    $:GPU_DECLARE(create='[pCr,TCr,mixM,lp,vp,A,B,C,D]')

contains

    !>  The purpose of this subroutine is to normalize and denormalize the code 
        !!      such that phase change bubble models can be merged easily.
        !!      This function will be used at the beginning of each time step
        !!      (normalize) and after m_phase_change.fpp (denormalize) module
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param ND string to toogle Normalize and Denormalize inputs
    impure subroutine s_norm_denorm(ND, q_cons_vf, q_prim_vf, q_T_sf)
                
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf, q_prim_vf
        type(scalar_field), intent(inout) :: q_T_sf
        character(len=1), intent(in) :: ND

        ! choosing parameters based on the whether the variables are being
        ! (de)normalized
        select case(ND)
            case ("N")
                AlphaRhoND = rho0, AlphaRhoEND = , TotEND = , MomND = 
            case ("D")
                AlphaRhoND = 1.0e0_wp ./ rho0, AlphaRhoEND = 1.0e0_wp ./ , TotEND = 1.0e0_w / , MomND = 1.0e0_wp /
        end select

        $:GPU_PARALLEL_LOOP(collapse=3, private='[pS,pSOV,pSSL,TS,TSatOV,TSatSL,TSOV,TSSL, &
        & rhoe,rhoeT,dynE,rhos,rho,rM,m1,m2,TR,p_infOV,p_infpT,p_infSL,alphak,alpharhoe0k,m0k,rhok,Tk]')
        ! (de)normalizing the conservative variables
        do j = 0, m
            do k = 0, n
                do l = 0, p
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        ! mass factions
                        q_cons_vf(i + contxb - 1)%sf(j, k, l) = q_cons_vf(i + contxb - 1)%sf(j, k, l) / AlphaRhoND

                        ! volume fractions
                        q_cons_vf(i + advxb - 1)%sf(j, k, l) = q_cons_vf(i + advxb - 1)%sf(j, k, l)

                        ! mixture energies
                        if (model_eqns == 3) then
                            q_cons_vf(i + intxb - 1)%sf(j, k, l) = q_cons_vf(i + intxb - 1)%sf(j, k, l) / AlphaRhoEND
                        end if

                        ! total energy
                        q_cons_vf(E_idx)%sf(j, k, l) = q_cons_vf(E_idx)%sf(j, k, l) / TotEND

                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    ! momentum
                    do i = momxb, momxe
                        q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l) / MomND
                    end do

                    ! updating the respective primitive variables. As of 08/Jul/2025, I am not concerned with
                    ! q_T_sf
                    call s_convert_conservative_to_primitive_variables(q_cons_vf, q_T_sf, q_prim_vf, idwint)

                    ! parsing them to GPU
                    do i = 1, sys_size
                        $:GPU_UPDATE(host='[q_cons_vf(i)%sf(:,:,:)]')
                        $:GPU_UPDATE(host='[q_prim_vf(i)%sf(:,:,:)]')
                    end do
                end do
            end do
        end do

    end subroutine s_norm_denorm !-------------------------------

end module m_normalize