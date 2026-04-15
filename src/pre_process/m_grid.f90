!>
!! @file m_grid.f90
!! @brief Contains module m_grid

!> @brief  This module takes care of creating the rectilinear grid on which
!!              the data for the initial condition will be laid out and on which
!!              the simulation will eventually be computed. The grid may either
!!              be uniform or non-uniform. Non-uniform grids are generated using
!!              the hyperbolic tangent function, see Johnsen (2007) for details.
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

        call s_stretch_grid(x_cb, m, x_domain, stretch_x, a_x, x_a, x_b, loops_x, x_cc, dx)

        if (stretch_x) then
            print *, 'Stretched grid: min/max x grid: ', minval(x_cc(:)), maxval(x_cc(:))
            if (num_procs > 1) call s_mpi_reduce_min(dx)
        end if

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

        call s_stretch_grid(y_cb, n, y_domain, stretch_y, a_y, y_a, y_b, loops_y, y_cc, dy)

        if (stretch_y .and. (num_procs > 1)) call s_mpi_reduce_min(dy)

        ! Grid Generation in the z-direction
        if (p == 0) return

        dz = (z_domain%end - z_domain%beg)/real(p + 1, wp)

        do i = 0, p
            z_cb(i - 1) = z_domain%beg + dz*real(i, wp)
        end do

        z_cb(p) = z_domain%end

        call s_stretch_grid(z_cb, p, z_domain, stretch_z, a_z, z_a, z_b, loops_z, z_cc, dz)

        if (stretch_z .and. (num_procs > 1)) call s_mpi_reduce_min(dz)

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
        call s_stretch_grid(x_cb_glb, m_glb, x_domain, stretch_x, a_x, x_a, x_b, loops_x)

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
            call s_stretch_grid(y_cb_glb, n_glb, y_domain, stretch_y, a_y, y_a, y_b, loops_y)

            ! Grid generation in the z-direction
            if (p_glb > 0) then
                dz = (z_domain%end - z_domain%beg)/real(p_glb + 1, wp)
                do i = 0, p_glb
                    z_cb_glb(i - 1) = z_domain%beg + dz*real(i, wp)
                end do
                z_cb_glb(p_glb) = z_domain%end
                call s_stretch_grid(z_cb_glb, p_glb, z_domain, stretch_z, a_z, z_a, z_b, loops_z)
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

    !> Apply the current 1D stretching formula to one coordinate direction.
        !! This helper centralizes the stretching logic so the serial and
        !! parallel grid generators can share the same implementation.
        !! Center recomputation and minimum-spacing reporting are optional so
        !! the same routine can be used both for local and global grids.
    impure subroutine s_stretch_grid(cb, cell_end, domain, stretch, a, coord_a, coord_b, loops, cc, min_spacing)

        integer, intent(in) :: cell_end
        real(wp), intent(inout) :: cb(-1:cell_end)
        type(bounds_info), intent(in) :: domain
        logical, intent(in) :: stretch
        real(wp), intent(in) :: a, coord_a, coord_b
        integer, intent(in) :: loops
        real(wp), intent(out), optional :: cc(0:cell_end)
        real(wp), intent(out), optional :: min_spacing

        real(wp) :: length, a_loc, coord_a_loc, coord_b_loc
        integer :: i, j

        if (stretch) then

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

        end if

        if (present(cc)) cc(0:cell_end) = (cb(0:cell_end) + cb(-1:cell_end - 1))/2._wp
        if (present(min_spacing)) min_spacing = minval(cb(0:cell_end) - cb(-1:cell_end - 1))

    end subroutine s_stretch_grid

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
