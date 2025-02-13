#:def Hardcoded1DVariables()
    ! Place any declaration of intermediate variables here

#:enddef

#:def Hardcoded1D()

    select case (patch_icpp(patch_id)%hcid)
    case (100) ! 1D Pressure Pulse - Solution to Wave Equation, 1D cartesian

        q_prim_vf(momxb)%sf(i, j, 0) =  patch_icpp(1)%vel(1) &
                                        * ( exp( - ( x_cc(i) - 1 ) ** 2 / 2 ) &
                                        -   exp( - ( x_cc(i) + 1 ) ** 2 / 2 ) )

    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
