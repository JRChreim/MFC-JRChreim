!>
!! @file m_grid.f90
!! @brief Contains module m_grid

!> @brief  This module takes care of creating the rectilinear grid on which
!!              the data for the initial condition will be laid out and on which
!!              the simulation will eventually be computed. The grid may either
!!              be uniform or non-uniform. Non-uniform grids can be generated
!!              using either the hyperbolic tangent mapping of Johnsen (2007)
!!              or a geometric-progression construction with a refined core.
!!              Alternatively to synthesizing a new grid, the user may select to
!!              read in a preexisting one. This is carried out through the module
!!              m_start_up.f90. In such a case, the responsibility of this module
!!              becomes only to allocate/deallocate the necessary grid variables
!!              for the cell-centers and cell-boundaries locations.
module m_grid

    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_mpi_proxy             ! Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

#ifdef MFC_MPI
    use mpi                     ! Message passing interface (MPI) module
#endif

    implicit none

    private; 
    public :: s_initialize_grid_module, &
              s_generate_grid, &
              s_generate_serial_grid, &
              s_generate_parallel_grid, &
              s_finalize_grid_module

    abstract interface

        impure subroutine s_generate_abstract_grid

        end subroutine s_generate_abstract_grid

    end interface

    procedure(s_generate_abstract_grid), pointer :: s_generate_grid => null()

    integer, parameter  :: geom_ratio_max_iter = 200
    real(wp), parameter :: geom_ratio_tol = 1.0e-12_wp
    real(wp), parameter :: geom_ratio_max = 1.0e5_wp

contains

    !> The following subroutine generates either a uniform or
        !!              non-uniform rectilinear grid in serial, defined by the parameters
        !!              inputted by the user. The grid information is stored in
        !!              the grid variables containing coordinates of the cell-
        !!              centers and cell-boundaries.
    impure subroutine s_generate_serial_grid

        ! Generic loop iterator
        integer :: i             !< generic loop operators

        ! Grid Generation in the x-direction
        dx = (x_domain%end - x_domain%beg)/real(m + 1, wp)

        do i = 0, m
            x_cb(i - 1) = x_domain%beg + dx*real(i, wp)
        end do

        x_cb(m) = x_domain%end

        call s_stretch_grid_by_type(x_cb, m, x_domain, stretch_x, stretch_type, a_x, x_a, x_b, loops_x, x_cc, dx)

        ! Grid Generation in the y-direction
        if (n == 0) return

        if (grid_geometry == 2 .and. f_approx_equal(y_domain%beg, 0.0_wp)) then
            !IF (grid_geometry == 2) THEN

            dy = (y_domain%end - y_domain%beg)/real(2*n + 1, wp)

            y_cb(-1) = y_domain%beg

            do i = 1, n
                y_cb(i - 1) = y_domain%beg + dy*real(2*i - 1, wp)
            end do

        else

            dy = (y_domain%end - y_domain%beg)/real(n + 1, wp)

            do i = 0, n
                y_cb(i - 1) = y_domain%beg + dy*real(i, wp)
            end do

        end if

        y_cb(n) = y_domain%end

        call s_stretch_grid_by_type(y_cb, n, y_domain, stretch_y, stretch_type, a_y, y_a, y_b, loops_y, y_cc, dy)

        ! Grid Generation in the z-direction
        if (p == 0) return

        dz = (z_domain%end - z_domain%beg)/real(p + 1, wp)

        do i = 0, p
            z_cb(i - 1) = z_domain%beg + dz*real(i, wp)
        end do

        z_cb(p) = z_domain%end

        call s_stretch_grid_by_type(z_cb, p, z_domain, stretch_z, stretch_type, a_z, z_a, z_b, loops_z, z_cc, dz)

    end subroutine s_generate_serial_grid

    !> The following subroutine generates either a uniform or
        !!              non-uniform rectilinear grid in parallel, defined by the parameters
        !!              inputted by the user. The grid information is stored in
        !!              the grid variables containing coordinates of the cell-
        !!              centers and cell-boundaries.
    impure subroutine s_generate_parallel_grid

#ifdef MFC_MPI

        ! Locations of cell boundaries
        real(wp), allocatable, dimension(:) :: x_cb_glb, y_cb_glb, z_cb_glb !<
            !! Locations of cell boundaries

        character(LEN=path_len + name_len) :: file_loc !<
            !! Generic string used to store the address of a file

        integer :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE) :: status

        integer :: i !< Generic loop integers

        allocate (x_cb_glb(-1:m_glb))
        allocate (y_cb_glb(-1:n_glb))
        allocate (z_cb_glb(-1:p_glb))

        ! Grid generation in the x-direction
        dx = (x_domain%end - x_domain%beg)/real(m_glb + 1, wp)
        do i = 0, m_glb
            x_cb_glb(i - 1) = x_domain%beg + dx*real(i, wp)
        end do
        x_cb_glb(m_glb) = x_domain%end
        call s_stretch_grid_by_type(x_cb_glb, m_glb, x_domain, stretch_x, stretch_type, a_x, x_a, x_b, loops_x)

        ! Grid generation in the y-direction
        if (n_glb > 0) then

            if (grid_geometry == 2 .and. f_approx_equal(y_domain%beg, 0.0_wp)) then
                dy = (y_domain%end - y_domain%beg)/real(2*n_glb + 1, wp)
                y_cb_glb(-1) = y_domain%beg
                do i = 1, n_glb
                    y_cb_glb(i - 1) = y_domain%beg + dy*real(2*i - 1, wp)
                end do
            else
                dy = (y_domain%end - y_domain%beg)/real(n_glb + 1, wp)
                do i = 0, n_glb
                    y_cb_glb(i - 1) = y_domain%beg + dy*real(i, wp)
                end do
            end if
            y_cb_glb(n_glb) = y_domain%end
            call s_stretch_grid_by_type(y_cb_glb, n_glb, y_domain, stretch_y, stretch_type, a_y, y_a, y_b, loops_y)

            ! Grid generation in the z-direction
            if (p_glb > 0) then
                dz = (z_domain%end - z_domain%beg)/real(p_glb + 1, wp)
                do i = 0, p_glb
                    z_cb_glb(i - 1) = z_domain%beg + dz*real(i, wp)
                end do
                z_cb_glb(p_glb) = z_domain%end
                call s_stretch_grid_by_type(z_cb_glb, p_glb, z_domain, stretch_z, stretch_type, a_z, z_a, z_b, loops_z)
            end if
        end if

        ! Write cell boundary locations to grid data files
        file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'x_cb.dat'
        data_size = m_glb + 2
        call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                           mpi_info_int, ifile, ierr)
        call MPI_FILE_WRITE(ifile, x_cb_glb, data_size, mpi_p, status, ierr)
        call MPI_FILE_CLOSE(ifile, ierr)

        if (n > 0) then
            file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'y_cb.dat'
            data_size = n_glb + 2
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                               mpi_info_int, ifile, ierr)
            call MPI_FILE_WRITE(ifile, y_cb_glb, data_size, mpi_p, status, ierr)
            call MPI_FILE_CLOSE(ifile, ierr)

            if (p > 0) then
                file_loc = trim(case_dir)//'/restart_data'//trim(mpiiofs)//'z_cb.dat'
                data_size = p_glb + 2
                call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                                   mpi_info_int, ifile, ierr)
                call MPI_FILE_WRITE(ifile, z_cb_glb, data_size, mpi_p, status, ierr)
                call MPI_FILE_CLOSE(ifile, ierr)
            end if
        end if

        deallocate (x_cb_glb, y_cb_glb, z_cb_glb)

#endif

    end subroutine s_generate_parallel_grid

    !> Dispatch the grid stretching routine according to the user-selected
        !! stretch type. 'stretch_type = 1' selects the hyperbolic tangent map
        !! and 'stretch_type = 2' selects the geometric-progression core.
        !! The geometric-progression count is passed in as a real value and
        !! converted explicitly here before reaching the integer-only helper.
        !! Cell-center and minimum-spacing outputs are optional so the MPI
        !! grid path can skip them.
    impure subroutine s_stretch_grid_by_type(cb, cell_end, domain, stretch, stretch_type, a, coord_a, coord_b, loops, cc, dS)

        integer, intent(in) :: cell_end
        real(wp), intent(inout) :: cb(-1:cell_end)
        type(bounds_info), intent(in) :: domain
        logical, intent(in) :: stretch
        integer, intent(in) :: stretch_type
        real(wp), intent(in) :: a, coord_a, coord_b
        integer, intent(in) :: loops
        real(wp), intent(out), optional :: cc(0:cell_end)
        real(wp), intent(out), optional :: dS

        integer :: num_refined

        if (stretch) then
            select case (stretch_type)
                case (1)
                    call s_stretch_grid_hyper_tan(cb, cell_end, domain, a, coord_a, coord_b, loops)
                case (2)
                    num_refined = nint(a)
                    call s_stretch_grid_geom_prog(cb, cell_end, domain, coord_a, coord_b, num_refined)
            end select
        end if

        if (have_cc .neqv. have_dS) then
            call s_mpi_abort('s_stretch_grid_by_type requires both cc and dS or neither.')
        end if

        if (present(cc)) then
            cc(0:cell_end) = (cb(0:cell_end) + cb(-1:cell_end - 1))/2._wp
            print *, 'Stretched grid: min/max [x,y,z] grid: ', minval(cc(:)), maxval(cc(:))
        end if
        if (present(dS)) then
            dS = minval(cb(0:cell_end) - cb(-1:cell_end - 1))
            if (num_procs > 1) call s_mpi_reduce_min(dS)
        end if            

    end subroutine s_stretch_grid_by_type

    !> Apply the current hyperbolic-tangent stretching formula to one
        !! coordinate direction.
        !! This helper centralizes the stretching logic so the serial and
        !! parallel grid generators can share the same implementation.
        !! Center recomputation and minimum-spacing reporting are optional so
        !! the same routine can be used both for local and global grids.
    impure subroutine s_stretch_grid_hyper_tan(cb, cell_end, domain, a, coord_a, coord_b, loops)

        integer, intent(in) :: cell_end
        real(wp), intent(inout) :: cb(-1:cell_end)
        type(bounds_info), intent(in) :: domain
        real(wp), intent(in) :: a, coord_a, coord_b
        integer, intent(in) :: loops

        real(wp) :: length, a_loc, coord_a_loc, coord_b_loc
        integer :: i, j


        length = abs(domain%end - domain%beg)
        a_loc = a
        coord_a_loc = coord_a/length
        coord_b_loc = coord_b/length

        cb = cb/length

        do j = 1, loops
            do i = -1, cell_end
                cb(i) = cb(i)/a_loc* &
                        (a_loc + log(cosh(a_loc*(cb(i) - coord_a_loc))) &
                          + log(cosh(a_loc*(cb(i) - coord_b_loc))) &
                          - 2._wp*log(cosh(a_loc*(coord_b_loc - coord_a_loc)/2._wp)))
            end do
        end do

        cb = cb*length

    end subroutine s_stretch_grid_hyper_tan

    !> Build a 1D mesh using a uniform refined core and geometric tails on the
        !! left and right sides. This mirrors the Python mesh_stretch_geom_prog
        !! workflow and is intentionally separate from the log-cosh stretching.
    impure subroutine s_stretch_grid_geom_prog(cb, cell_end, domain, coord_a, coord_b, num_refined)

        integer, intent(in) :: cell_end
        real(wp), intent(inout) :: cb(-1:cell_end)
        type(bounds_info), intent(in) :: domain
        real(wp), intent(in) :: coord_a, coord_b
        integer, intent(in) :: num_refined

        integer :: i, n_left, n_right, n_stretched, n_total
        real(wp) :: first_spacing, left_len, right_len, ratio_left, ratio_right
        real(wp) :: refined_beg_loc, refined_end_loc, spacing, uniform_length

        if (f_is_default(coord_a) .or. f_is_default(coord_b) .or. (num_refined <= 0)) then
            call s_mpi_abort('Geometric progression grid stretching requires '// &
                             'refined bounds and a positive number of refined cells.')
        end if

        if (coord_b <= coord_a) then
            call s_mpi_abort('Geometric progression grid stretching requires coord_b > coord_a.')
        end if

        n_total = cell_end + 1

        if (num_refined > n_total) then
            call s_mpi_abort('Geometric progression grid stretching requires num_refined <= total number of cells.')
        end if

        n_stretched = n_total - num_refined
        refined_beg_loc = coord_a
        refined_end_loc = coord_b

        left_len = refined_beg_loc - domain%beg
        right_len = domain%end - refined_end_loc

        if (n_stretched == 0) then
            if ((.not. f_approx_equal(refined_beg_loc, domain%beg)) .or. &
                (.not. f_approx_equal(refined_end_loc, domain%end))) then
                call s_mpi_abort('Geometric progression grid stretching with '// &
                                 'zero stretched cells requires the refined bounds '// &
                                 'to span the full domain.')
            end if

            n_left = 0
            n_right = 0
            refined_beg_loc = domain%beg
            refined_end_loc = domain%end

        elseif ((left_len <= 0._wp) .and. (right_len <= 0._wp)) then
            call s_mpi_abort('Geometric progression grid stretching leaves no room for stretched cells outside the refined region.')

        elseif (left_len <= 0._wp) then
            n_left = 0
            n_right = n_stretched
            refined_beg_loc = domain%beg

        elseif (right_len <= 0._wp) then
            n_left = n_stretched
            n_right = 0
            refined_end_loc = domain%end

        else
            n_left = nint(real(n_stretched, wp)*left_len/(left_len + right_len))
            n_left = max(0, min(n_stretched, n_left))
            n_right = n_stretched - n_left
        end if

        first_spacing = (refined_end_loc - refined_beg_loc)/real(num_refined, wp)

        cb(-1) = domain%beg

        if (n_left > 0) then
            ratio_left = f_solve_geometric_ratio(first_spacing, n_left, left_len)

            do i = 1, n_left
                spacing = first_spacing*ratio_left**real(n_left - i, wp)
                cb(i - 1) = cb(i - 2) + spacing
            end do
        end if

        do i = 1, num_refined
            cb(n_left + i - 1) = cb(n_left + i - 2) + first_spacing
        end do

        if (n_right > 0) then
            ratio_right = f_solve_geometric_ratio(first_spacing, n_right, right_len)

            do i = 1, n_right
                spacing = first_spacing*ratio_right**real(i - 1, wp)
                cb(n_left + num_refined + i - 1) = cb(n_left + num_refined + i - 2) + spacing
            end do
        end if

        cb(cell_end) = domain%end

    end subroutine s_stretch_grid_geom_prog

    !> Solve for the geometric ratio that matches a stretched length and a
        !! fixed first spacing.
    real(wp) function f_solve_geometric_ratio(first_spacing, count, stretched_length) result(ratio)

        real(wp), intent(in) :: first_spacing, stretched_length
        integer, intent(in) :: count

        real(wp) :: lower, upper, f_lower, f_upper, mid, f_mid, uniform_length
        integer :: iter

        if (count <= 0) then
            ratio = 1._wp
            return
        end if

        uniform_length = real(count, wp)*first_spacing

        if (f_approx_equal(stretched_length, uniform_length)) then
            ratio = 1._wp
            return
        end if

        if (stretched_length < uniform_length) then
            lower = geom_ratio_tol
            upper = 1._wp - geom_ratio_tol
        else
            lower = 1._wp + geom_ratio_tol
            upper = geom_ratio_max
        end if

        f_lower = f_geometric_ratio_residual(lower, first_spacing, count, stretched_length)
        f_upper = f_geometric_ratio_residual(upper, first_spacing, count, stretched_length)

        if (stretched_length > uniform_length) then
            do while ((f_lower*f_upper > 0._wp) .and. (upper < huge(1._wp)/10._wp))
                upper = min(upper*10._wp, huge(1._wp)/10._wp)
                f_upper = f_geometric_ratio_residual(upper, first_spacing, count, stretched_length)
            end do
        end if

        if (f_lower*f_upper > 0._wp) then
            call s_mpi_abort('Could not bracket the geometric progression ratio for the requested mesh stretching.')
        end if

        do iter = 1, geom_ratio_max_iter
            mid = 0.5_wp*(lower + upper)
            f_mid = f_geometric_ratio_residual(mid, first_spacing, count, stretched_length)

            if (abs(f_mid) <= 1.0e-12_wp*max(1._wp, abs(uniform_length))) exit
            if (abs(upper - lower) <= 1.0e-12_wp*max(1._wp, abs(mid))) exit

            if (f_lower*f_mid <= 0._wp) then
                upper = mid
                f_upper = f_mid
            else
                lower = mid
                f_lower = f_mid
            end if
        end do

        ratio = mid

    end function f_solve_geometric_ratio

    !> Residual used to solve the geometric ratio.
    real(wp) function f_geometric_ratio_residual(ratio, first_spacing, count, stretched_length) result(residual)

        real(wp), intent(in) :: ratio, first_spacing, stretched_length
        integer, intent(in) :: count

        real(wp) :: log_term, geometric_tail

        if (ratio <= 0._wp) then
            residual = huge(1._wp)
            return
        end if

        if (count == 0) then
            residual = stretched_length*(ratio - 1._wp)
            return
        end if

        log_term = real(count, wp)*log(ratio)

        if (log_term > 700._wp) then
            geometric_tail = huge(1._wp)
        else
            geometric_tail = exp(log_term) - 1._wp
        end if

        residual = stretched_length*(ratio - 1._wp) - first_spacing*geometric_tail

    end function f_geometric_ratio_residual

    !> Computation of parameters, allocation procedures, and/or
        !!              any other tasks needed to properly setup the module
    impure subroutine s_initialize_grid_module

        if (parallel_io .neqv. .true.) then
            s_generate_grid => s_generate_serial_grid
        else
            s_generate_grid => s_generate_parallel_grid
        end if

    end subroutine s_initialize_grid_module

    !> Deallocation procedures for the module
    impure subroutine s_finalize_grid_module

        s_generate_grid => null()

    end subroutine s_finalize_grid_module

end module m_grid
