! ------------------------------------------------------------------------------
! ----- Module: mod_laplace_solver ---------------------------------------------
! ----- Author: Arnar Jonsson --------------------------------------------------
! ----- Date: 2025 -------------------------------------------------------------
! ------------------------------------------------------------------------------

module mod_laplace_solver
    
    ! Dependencies
    use mod_global
    use mod_verlet
    use mod_pair
    use mod_work_function
    use mkl_pardiso ! Required for pardiso solver
    use mkl_spblas  ! Required for sparse matrices
    ! use iso_c_binding, only: c_int, c_ptr, c_f_pointer

    implicit none

    ! Accessibility
    private
    public :: Init_Laplace_Solver, Calculate_Laplace_Field, Calculate_Laplace_Field_at, &
                Clean_Up_Laplace, Place_Electron, Write_Laplace_Data !, Sample_Field

    integer(kind=8), parameter :: pardiso_precision = 8

    ! Pardiso variables
    type(MKL_PARDISO_HANDLE), dimension(64) :: pt
    integer(kind=8), dimension(64) :: iparm

    integer(kind=8), parameter :: maxfct=1, mnum=1, mtype=11, nrhs=1, msglvl=1
    integer(kind=8), parameter :: analysis=11, factorization=22, solving=33, solving1=331, solving2=333

    integer(kind=8) :: ndiff
    integer(kind=8), allocatable, dimension(:) :: perm

    ! Matrix
    real(kind=pardiso_precision), allocatable, dimension(:) :: nnz_values
    integer(kind=8), allocatable, dimension(:) :: nnz_col_index, nnz_ia, nnz_center
    integer(kind=8) :: nnz

    ! System
    real(kind=pardiso_precision), allocatable, dimension(:) :: b
    real(kind=8), allocatable, dimension(:) :: iCharge, oCharge

    ! Solution
    real(kind=pardiso_precision), allocatable, dimension(:) :: x
    real(kind=8), allocatable, dimension(:) :: voltage
    real(kind=8), allocatable, dimension(:,:) :: laplace_field

    ! Discretization
    integer(kind=8), allocatable, dimension(:) :: gridPoints, gridPointsActive
    real(kind=8) :: hx, hy, hz
    real(kind=8) :: emitter_radius, min_d
    real(kind=8), dimension(3) :: h
    real(kind=8), dimension(3) :: div_h, div_h2, div_h3
    real(kind=8), dimension(2) :: lim_x, lim_y, lim_z
    integer(kind=8) :: Nx, Ny, Nz, nrGrid, nrGridActive, NxNy

contains

    subroutine Init_Laplace_Solver() ! DONE
        ! Initialize environment based on input file

        print *, 'Laplace: initializing grid'
        call init_grid()

        print *, 'Laplace: initializing matrix'
        call init_matrix()

        print *, 'Laplace: initializing pardiso'
        call init_pardiso()

        ! call write_grid()

    end subroutine Init_Laplace_Solver

    subroutine Calculate_Laplace_Field() ! TODO
        ! Calculate the electric field
        ! print *, 'Laplace: finding electrons'
        call find_electrons()
        ! print *, 'Laplace: updating matrix'
        call update_matrix()
        ! print *, 'Laplace: solving matrix'
        call solve_system()
        ! print *, 'Laplace: allocate voltage'
        call allocate_voltage()
        ! print *, 'Laplace: calculating field'
        call calculate_field()
        ! print *, 'Laplace: done'

        ! call write_average_field()
        
    end subroutine Calculate_Laplace_Field

! ------------------------------------------------------------------------------
! ----- Solution ---------------------------------------------------------------
! ------------------------------------------------------------------------------

    subroutine solve_system() ! TODO
        ! Solve the system
        integer(kind=8) :: i
        real(kind=8) :: start, end

        if (allocated(laplace_field)) then
            deallocate(laplace_field)
        end if

        call pardiso_phase(33,0)

    end subroutine solve_system

! -------------------------------------------------------------------------------------------------------------
! ----- Initialization ----------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------------------

    subroutine init_grid() ! DONE
        ! Initialize parameters and store points that are outside of the emitter or one layer within it
        integer(kind=8) :: i, a
        integer(kind=8), dimension(3) :: boundaries
        real(kind=8) :: Lx, Ly, Lz

        Nx = laplace_intervals(1); Ny = laplace_intervals(2); Nz = laplace_intervals(3)

        Lx = laplace_dim(1); Ly = laplace_dim(2); Lz = laplace_dim(3)

        NxNy = Nx*Ny

        nrGrid = Nx*Ny*Nz

        hx = Lx / (Nx-1)
        hy = Ly / (Ny-1)
        hz = Lz / (Nz-1)

        h = (/hx, hy, hz/)

        min_d = length_scale**2

        div_h = (/1.0d0 / hx, 1.0d0 / hy, 1.0d0 / hz/)
        div_h2 = (/1.0d0 / (hx**2), 1.0d0 / (hy**2), 1.0d0 / (hz**2)/)
        div_h3 = (/1.0d0 / (hx**3), 1.0d0 / (hy**3), 1.0d0 / (hz**3)/)

        lim_x = (/laplace_pos(1), laplace_pos(1)+Lx/)
        lim_y = (/laplace_pos(2), laplace_pos(2)+Ly/)
        lim_z = (/laplace_pos(3), laplace_pos(3)+Lz/)

        emitter_radius = emitters_dim(1,1)

        allocate(gridPoints(nrGrid), gridPointsActive(nrGrid))

        nrGridActive = 0
        nnz = 0

        do i=1,nrGrid
            if (is_emitter(i) == 1) then
                if (is_first_layer(i) == 1) then
                    nrGridActive = nrGridActive + 1
                    gridPointsActive(nrGridActive) = i
                    gridPoints(i) = nrGridActive

                    nnz = nnz + 1
                else
                    gridPoints(i) = 0
                end if
            ! else if (is_anode(i) == 1) then
            !     nrGridActive = nrGridActive + 1
            !     gridPointsActive(nrGridActive) = i
            !     gridPoints(i) = nrGridActive

            !     nnz = nnz + 1
            else
                nrGridActive = nrGridActive + 1
                gridPointsActive(nrGridActive) = i
                gridPoints(i) = nrGridActive

                boundaries = is_boundary(i)

                nnz = nnz+1

                do a=1,3
                    select case (boundaries(a))
                        case (0)
                            nnz = nnz+2
                        case (1)
                            nnz = nnz+3
                        case (2)
                            nnz = nnz+3
                    end select
                end do
            end if
        end do

        allocate(iCharge(nrGridActive), oCharge(nrGridActive), x(nrGridActive))
        iCharge = 0
        oCharge = 0

    end subroutine init_grid

    subroutine init_matrix() ! DONE
        ! Initialize matrix without boundary conditions
        integer(kind=8) :: g, i, a, center, next_center
        integer(kind=8), dimension(3) :: boundaries
        integer(kind=8), dimension(3,3) :: indices
        integer(kind=8) :: next, next_next, next_next_next, prev, prev_prev, prev_prev_prev

        allocate(nnz_values(nnz), nnz_col_index(nnz), nnz_ia(nrGridActive+1))
        allocate(b(nrGridActive), nnz_center(nrGridActive))

        nnz_values = 0.0d0
        b = 0.0d0

        center = 1

        do g=1,nrGridActive

            i = gridPointsActive(g)

            if (is_emitter(i) == 1) then

                ! Place the center
                nnz_center(g) = center
                next_center = center + 1

                ! Place the column indices
                nnz_ia(g) = center
                nnz_col_index(center) = g

                ! Place the values
                call insert_emitter_boundary(g)

            ! else if (is_anode(i) == 1) then
            !     ! Place the center
            !     nnz_center(g) = center
            !     next_center = center + 1

            !     ! Place the column indices
            !     nnz_ia(g) = center
            !     nnz_col_index(center) = g

            !     ! Place the values
            !     call insert_anode_boundary(g)

            else

                boundaries = is_boundary(i)

                ! Place the center
                do a=1,3
                    select case (boundaries(a))
                        case (0)
                            center = center + 1 ! make space for one to the left
                        case (1)
                            center = center ! do nothing
                        case (2)
                            center = center + 3 ! make space for three to the left
                    end select
                end do
                nnz_center(g) = center
                next_center = center + 1
                nnz_ia(g) = center

                ! Place the column indices
                nnz_col_index(center) = g

                call get_finite_indices(g,indices)
            
                do a=1,3
                    select case (boundaries(a))
                        case (0)
                            prev = indices(a,1)
                            next = indices(a,3)

                            nnz_col_index(prev) = move(g,-1,a)
                            nnz_col_index(next) = move(g,1,a)

                            next_center = next + 1
                            nnz_ia(g) = prev
                        case (1)
                            next = indices(a,1)
                            next_next = indices(a,2)
                            next_next_next = indices(a,3)

                            nnz_col_index(next) = move(g,1,a)
                            nnz_col_index(next_next) = move(g,2,a)
                            nnz_col_index(next_next_next) = move(g,3,a)

                            next_center = next_next_next + 1
                        case (2)
                            prev = indices(a,3)
                            prev_prev = indices(a,2)
                            prev_prev_prev = indices(a,1)

                            nnz_col_index(prev) = move(g,-1,a)
                            nnz_col_index(prev_prev) = move(g,-2,a)
                            nnz_col_index(prev_prev_prev) = move(g,-3,a)      
                            
                            nnz_ia(g) = prev_prev_prev
                    end select
                end do

                ! Place the values
                call insert_finite_difference(g)
            end if

            center = next_center

        end do

        nnz_ia(nrGridActive+1) = nnz+1
    end subroutine init_matrix

    subroutine init_pardiso() ! DONE
        real(kind=8) :: start, end

        ! Initialize pardiso solver and do analysis, symbolic factorization, and numerical factorization
        if (pardiso_precision == 4) then
            iparm(28) = 1
        end if
        call pardisoinit(pt,mtype,iparm)
        ! Solver parameters
        iparm(1) = 1 ! use user-defined parameters
        iparm(2) = 3 ! 1 = minimum degree algorithm, 2 = metis, 3 = parallel metis
        ! iparm(3) = 0 ! reserved
        iparm(4) = 0 ! >0 = preconditioned CGS
        iparm(5) = 0 ! 0 = ignored; 1 = user supplied fill-in permutation; 2 = returns permutation
        iparm(6) = 0 ! 0 = write solution on x; 1 = write solution on b
        ! iparm(7) = 0 ! output number of iterative refinement steps
        ! iparm(8) =  1e3 ! maximum number of iterative refinement steps
        ! iparm(9) = 12 ! tolerance = 10^(-iparm(9))	
        iparm(10) = 13 ! pivoting perturbation: 13 = 10^-13; 8 = 10^-8
        iparm(11) = 1 ! 0 = no scaling; 1 = scaling
        iparm(12) = 0 ! 0 = solve linear system; 1 = solve conjugate transposed system, 2 = solve transposed system
        iparm(13) = 1 ! 0 = no symmetric weighted matching; 1 = symmetric weighted matching
        ! iparm(14) = 0 ! output number of perturbed pivots
        ! iparm(15) = 0 ! output peak memory on symbolic factorization
        ! iparm(16) = 0 ! output permanent memory on symbolic factorization
        ! iparm(17) = 0 ! output size of factors on numerical factorization and solving
        ! iparm(18) = 0 ! output number of non-zero elements in the factors
        ! iparm(19) = 0 ! output number of floating point operations in factorization
        ! iparm(20) = 0 ! output cgs diagnostics
        iparm(21) = 0 ! 0 = 1x1 diagonal pivoting; 1 = 1x1 and 2x2 pivoting
        ! iparm(22) = 0 ! output
        ! iparm(23) = 0 ! output
        iparm(24) = 0 ! 0 = classic factorization algorithm; 1 = two-level factorization; 10 = improved two-level factorization
        iparm(25) = 0 ! 0 = parallel forward/backward solve; 1 = sequential forward/backward solve
        ! iparm(26) = 0 ! reserved
        iparm(27) = 1 ! 0 = check matrix off; 1 = check matrix on
        if (pardiso_precision == 4) then
            iparm(28) = 1 ! 0 = double precision; 1 = single precision
        end if
        ! iparm(29) = 0 ! reserved
        ! iparm(30) = 0 ! output
        iparm(31) = 2 ! 0 = no partial solve
        ! iparm(32) = 0 ! reserved
        ! iparm(33) = 0 ! reserved
        ! iparm(34) = 0 ! 0 = CNR mode and oneapi determines optimal number of threads; 1 = CNR mode and user determines number of threads; 2 = CNR mode and user determines number of threads; 3 = CNR mode and user determines number of threads
        iparm(35) = 0 ! 0 = one-based indexing; 1 = zero-based indexing
        iparm(36) = 0 ! 0 = not compute shur; 1 = compute shur
        iparm(37) = 0 ! 0 = CSR format; 1 = BSR format
        ! iparm(38) = 0 ! reserved
        iparm(39) = 0 ! 0 = full factorization; 1 = low-rank update
        ! iparm(40) = 0 ! reserved
        ! iparm(41) = 0 ! reserved
        ! iparm(42) = 0 ! reserved
        iparm(43) = 0 ! 0 = do not compute diagonal inverse; 1 = compute diagonal inverse
        ! iparm(44) = 0 ! reserved
        ! iparm(45) = 0 ! reserved
        ! iparm(46) = 0 ! reserved
        ! iparm(47) = 0 ! reserved
        ! iparm(48) = 0 ! reserved
        ! iparm(49) = 0 ! reserved
        ! iparm(50) = 0 ! reserved
        ! iparm(51) = 0 ! reserved
        ! iparm(52) = 0 ! reserved
        ! iparm(53) = 0 ! reserved
        ! iparm(54) = 0 ! reserved
        ! iparm(55) = 0 ! reserved
        iparm(56) = 0 ! 0 = internal function used to work with pivot; 1 = use mkl_pardiso_pivot callback routine
        ! iparm(57) = 0 ! reserved
        ! iparm(58) = 0 ! reserved
        ! iparm(59) = 0 ! reserved
        iparm(60) = 0 ! 0 = in-core-memory; 1 = mix; 2 = out-of-core memory
        ! iparm(61) = 0 ! reserved
        ! iparm(62) = 0 ! reserved
        ! iparm(63) = 0 ! size of the minimum OOC memory for numerical factorization
        ! iparm(64) = 0 ! reserved

        allocate(perm(nrGridActive))

        ! call cpu_time(start)
        ! print *, 'Laplace: symbolic factorization'
        call pardiso_phase(11,1)
        iparm(27) = 0
        print *, 'Laplace: required memory:', iparm(17)*1e-6, 'Gb'
        ! print *, 'Laplace: numerical factorization'
        call pardiso_phase(22,1)
        ! print *, 'Laplace: solving'
        perm = 0
        call pardiso_phase(33,1)
        ! call cpu_time(end)
        
    end subroutine init_pardiso

    subroutine pardiso_phase(phase,msg)
        integer(kind=8), intent(in) :: phase, msg
        integer(kind=8) :: error

        call pardiso_64(pt,maxfct,mnum,mtype,phase,nrGridActive,nnz_values,nnz_ia,nnz_col_index,perm,nrhs,iparm,msg,b,x,error)


        if (error == 0) then
            ! print *, '  Laplace: analysis and symbolic factorization successful'
        else
            print *, '  Laplace: phase', phase, 'failed with error code ', error
        end if

    end subroutine pardiso_phase

! -------------------------------------------------------------------------------------------------------------
! ----- Voltage and field calculations ------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------------------

    subroutine allocate_voltage() ! DONE & PARALLEL
        integer(kind=8) :: i, g

        allocate(voltage(nrGrid))

        voltage = 0.0d0

        ! $OMP PARALLEL DO DEFAULT(NONE) &
        ! $OMP& SHARED(nrGrid, voltage, x, gridPoints) &
        ! $OMP& PRIVATE(i, g)

        do i=1,nrGrid
            if (gridPoints(i) /= 0) then
                g = gridPoints(i)
                voltage(i) = real(x(g),kind=8)

            end if
        end do

        ! $OMP END PARALLEL DO

    end subroutine allocate_voltage

    subroutine calculate_field() ! DONE & PARALLEL
        integer(kind=8) :: i, a
        integer(kind=8), dimension(3) :: boundaries
        real(kind=8), dimension(3) :: pos
        real(kind=8) :: l, l_const = q_0 / (4.0d0*pi*epsilon_0)

        allocate(laplace_field(nrGrid,3))

        laplace_field = 0.0d0

        ! $OMP PARALLEL DO DEFAULT(NONE) &
        ! $OMP& SHARED(nrGrid, voltage, laplace_field) &
        ! $OMP& PRIVATE(i, boundaries, boundary)

        do i=1,nrGrid
            ! print *, 'Laplace: i:', i, 'nrGrid:', nrGrid
            boundaries = is_boundary(i)
            do a=1,3
                select case (boundaries(a))
                    case (0) ! Inner point
                        ! print *, '  Laplace: central difference'
                        laplace_field(i,a) = - central_difference(voltage,i,a)
                        if (a==3) then
                            laplace_field(i,a) = laplace_field(i,a) + applied_field(i)
                        end if
                    case (1) ! Startpoint
                        ! print *, '  Laplace: forward difference'
                        laplace_field(i,a) = - forward_difference(voltage,i,a)
                        if (a==3) then
                            laplace_field(i,a) = laplace_field(i,a) + applied_field(i)
                        end if
                    case (2) ! Endpoint
                        ! print *, '  Laplace: backward difference'
                        laplace_field(i,a) = - backward_difference(voltage,i,a)
                        if (a==3) then
                            laplace_field(i,a) = laplace_field(i,a) + applied_field(i)
                        end if
                end select
            end do
        end do

        ! $OMP END PARALLEL DO

        ! Deallocate the voltage vector
        deallocate(voltage)

    end subroutine calculate_field

    function average_field() ! DONE & PARALLEL
        integer(kind=8) :: i, j, k, num
        real(kind=8) :: sum, average_field

        sum = 0.0d0
        num = 0

        ! $OMP PARALLEL DO DEFAULT(NONE) &
        ! $OMP& SHARED(Nx, Ny, laplace_field) &
        ! $OMP& PRIVATE(i, j, k) &
        ! $OMP& REDUCTION(+:sum, num)

        do i=0,Nx-1
            do j=0,Ny-1
                k = i + j*Nx+1
                sum = sum + laplace_field(k,3)
                num = num + 1
            end do
        end do

        ! $OMP END PARALLEL DO

        if (num > 0) then
            average_field = sum / num
        else
            average_field = 0.0d0
        end if
    end function average_field

    subroutine write_average_field() ! DONE
        integer(kind=8) :: IFAIL
        real(kind=8) :: avg_field

        avg_field = average_field()

        write (ud_laplace_average_field, "(E16.8)", iostat=IFAIL) avg_field
    end subroutine write_average_field

    subroutine write_grid() ! DONE
        integer(kind=8) :: IFAIL, i

        do i=1,NxNy
            write (ud_grid, "(i8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, i8)", iostat=IFAIL) is_emitter(i), cart_coord(i)
        end do
    end subroutine write_grid

    subroutine get_finite_indices(g,indices)
        integer(kind=8), intent(in) :: g
        integer(kind=8), dimension(3,3), intent(out) :: indices
        integer(kind=8) :: i, a, center, start, end
        integer(kind=8), dimension(3) :: boundaries

        i = gridPointsActive(g)
        boundaries = is_boundary(i)
        center = nnz_center(g)

        start = center
        end = center

        do a=1,3
            select case (boundaries(a))
                case (0)
                    start = start - 1
                    end = end + 1
                    indices(a,1) = start
                    indices(a,3) = end
                case (1)
                    end = end + 1
                    indices(a,1) = end
                    end = end + 1
                    indices(a,2) = end
                    end = end + 1
                    indices(a,3) = end
                case (2)
                    start = start - 1
                    indices(a,3) = start
                    start = start - 1
                    indices(a,2) = start
                    start = start - 1
                    indices(a,1) = start
            end select
        end do

    end subroutine

! -----------------------------------------------------------------------------
! ----- Matrix operations -----------------------------------------------------
! -----------------------------------------------------------------------------

    subroutine update_matrix()
        integer(kind=8) :: g, nrchanges
        nrchanges = 0
        do g=1,nrGridActive
            if (oCharge(g) /= iCharge(g)) then
                print *, 'oCharge:', oCharge(g), 'iCharge:', iCharge(g)
                b(g) = iCharge(g)
                oCharge(g) = iCharge(g)
                if (b(g) == 0.0d0) then
                    perm(g) = 0
                else
                    perm(g) = 1
                end if
                nrchanges = nrchanges + 1
            end if
        end do
        print *, 'Laplace: number of changes:', nrchanges

        ! print *, 'Laplace: number of changes:', ndiff

    end subroutine update_matrix

    subroutine insert_electron_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = nnz_center(g)

        nnz_values(center) = div_h2(3)
        b(g) = iCharge(g)*div_h2(3)
    end subroutine insert_electron_boundary

    ! subroutine remove_electron_boundary(g)
    !     integer(kind=8), intent(in) :: g
    !     integer(kind=8) :: center

    !     center = nnz_center(g)

    !     nnz_values(center) = 0.0d0
    !     b(g) = 0.0d0
    ! end subroutine remove_electron_boundary

    subroutine insert_emitter_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = nnz_center(g)

        nnz_values(center) = div_h2(3)
        b(g) = 0.0d0
    end subroutine insert_emitter_boundary

    subroutine insert_anode_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = nnz_center(g)

        nnz_values(center) = div_h2(3)
        b(g) = V_s*div_h2(3)
    end subroutine insert_anode_boundary

    subroutine insert_finite_difference(g) ! DONE
        ! Insert finite difference equation into matrix
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center
        integer(kind=8), dimension(3) :: boundaries
        integer(kind=8), dimension(3,3) :: indices
        integer(kind=8) :: i, a, next, next_next, next_next_next, prev, prev_prev, prev_prev_prev

        center = nnz_center(g)

        i = gridPointsActive(g)
        boundaries = is_boundary(i)
        call get_finite_indices(g,indices)

        ! Insert finite difference equation for each axis
        do a=1,3          
            select case (boundaries(a))
                case (0) ! Inner point (values for central difference)
                    prev = indices(a,1)
                    next = indices(a,3)

                    nnz_values(prev) = nnz_values(prev) + div_h2(a)
                    nnz_values(center) = nnz_values(center) - 2.0d0*div_h2(a)
                    nnz_values(next) = nnz_values(next) + div_h2(a)

                case (1) ! Startpoint (values for forward difference)
                    next = indices(a,1)
                    next_next = indices(a,2)
                    next_next_next = indices(a,3)

                    nnz_values(center) = nnz_values(center) + 2.0d0*div_h2(a)
                    nnz_values(next) = nnz_values(next) - 5.0d0*div_h2(a)
                    nnz_values(next_next) = nnz_values(next_next) + 4.0d0*div_h2(a)
                    nnz_values(next_next_next) = nnz_values(next_next_next) - div_h2(a)

                case (2) ! Endpoint (values for backward difference)
                    prev = indices(a,3)
                    prev_prev = indices(a,2)
                    prev_prev_prev = indices(a,1)

                    nnz_values(prev_prev_prev) = nnz_values(prev_prev_prev) - div_h2(a)
                    nnz_values(prev_prev) = nnz_values(prev_prev) + 4.0d0*div_h2(a)
                    nnz_values(prev) = nnz_values(prev) - 5.0d0*div_h2(a)
                    nnz_values(center) = nnz_values(center) + 2.0d0*div_h2(a)
            end select
        end do

        b(g) = 0.0d0
    end subroutine insert_finite_difference

!     subroutine remove_finite_difference(g) ! DONE
!         ! Remove finite difference equation from matrix
!         integer(kind=8), intent(in) :: g
!         integer(kind=8), dimension(3) :: boundaries
!         integer(kind=8), dimension(3,3) :: indices
!         integer(kind=8) :: i, a, center
!         integer(kind=8) :: prev_prev_prev, prev_prev, prev, next, next_next, next_next_next

!         center = nnz_center(g)
!         i = gridPointsActive(g)
!         boundaries = is_boundary(i)	
!         call get_finite_indices(g,indices)

!         nnz_values(center) = 0.0d0
!         do a=1,3
!             select case (boundaries(a))
!                 case (0)
!                     prev = indices(a,1)
!                     next = indices(a,3)

!                     nnz_values(prev) = 0.0d0
!                     nnz_values(next) = 0.0d0
!                 case (1)
!                     next = indices(a,1)
!                     next_next = indices(a,2)
!                     next_next_next = indices(a,3)

!                     nnz_values(next) = 0.0d0
!                     nnz_values(next_next) = 0.0d0
!                     nnz_values(next_next_next) = 0.0d0
!                 case (2)
!                     prev = indices(a,3)
!                     prev_prev = indices(a,2)
!                     prev_prev_prev = indices(a,1)

!                     nnz_values(prev) = 0.0d0
!                     nnz_values(prev_prev) = 0.0d0
!                     nnz_values(prev_prev_prev) = 0.0d0
!             end select
!         end do

!         b(g) = 0.0d0
!     end subroutine remove_finite_difference

! -----------------------------------------------------------------------------
! ----- Finite difference -----------------------------------------------------
! -----------------------------------------------------------------------------

    function central_difference(f,x,a)
        real(kind=8), dimension(nrGrid), intent(in) :: f
        integer(kind=8), intent(in) :: x, a
        integer(kind=8) :: x_m1h, x_m2h, x_p1h, x_p2h
        real(kind=8) :: central_difference

        x_m1h = move(x,-1,a)
        x_m2h = move(x,-2,a)
        x_p1h = move(x,1,a)
        x_p2h = move(x,2,a)

        ! central_difference = (f(x_m2h)-8.0d0*f(x_m1h)+8.0d0*f(x_p1h)-f(x_p2h))/(12.0d0*h(a))
        central_difference = (f(x_p1h)-f(x_m1h))/(2.0d0*h(a))
    end function central_difference

    function forward_difference(f,x,a)
        real(kind=8), dimension(nrGrid), intent(in) :: f
        integer(kind=8), intent(in) :: x, a
        integer(kind=8) :: x_p1h, x_p2h, x_p3h, x_p4h
        real(kind=8) :: forward_difference, p2_forward_difference, p3_forward_difference, p5_forward_difference
        real(kind=8) :: temp_forward_difference, temp_temp_forward_difference

        x_p1h = move(x,1,a)
        x_p2h = move(x,2,a)
        x_p3h = move(x,3,a)
        x_p4h = move(x,4,a)

        p2_forward_difference = (f(x)-f(x_p1h))/(h(a))
        p3_forward_difference = (-3.0d0*f(x)+4.0d0*f(x_p1h)-1.0d0*f(x_p2h))/(2.0d0*h(a))
        p5_forward_difference = (-25.0d0*f(x)+48.0d0*f(x_p1h)-36.0d0*f(x_p2h)+16.0d0*f(x_p3h)-3.0d0*f(x_p4h))/(12.0d0*h(a))

        forward_difference = min(p2_forward_difference,p3_forward_difference, p5_forward_difference)
        if (forward_difference == p2_forward_difference) then
            temp_forward_difference = min(p3_forward_difference, p5_forward_difference)
            temp_temp_forward_difference = max(p3_forward_difference, p5_forward_difference)
        else if (forward_difference == p3_forward_difference) then
            temp_forward_difference = min(p2_forward_difference, p5_forward_difference)
            temp_temp_forward_difference = max(p2_forward_difference, p5_forward_difference)
        else
            temp_forward_difference = min(p2_forward_difference, p3_forward_difference)
            temp_temp_forward_difference = max(p2_forward_difference, p3_forward_difference)
        end if

        if (temp_forward_difference*forward_difference > 0.0d0) then
            forward_difference = temp_forward_difference
        end if
        if (temp_temp_forward_difference*forward_difference > 0.0d0) then
            forward_difference = temp_temp_forward_difference
        end if

        ! forward_difference = p3_forward_difference

    end function forward_difference

    function backward_difference(f,x,a)
        real(kind=8), dimension(nrGrid), intent(in) :: f
        integer(kind=8), intent(in) :: x, a
        integer(kind=8) :: x_m1h, x_m2h, x_m3h, x_m4h
        real(kind=8) :: backward_difference

        x_m1h = move(x,-1,a)
        x_m2h = move(x,-2,a)
        x_m3h = move(x,-3,a)
        x_m4h = move(x,-4,a)

        backward_difference = (25.0d0*f(x)-48.0d0*f(x_m1h)+36.0d0*f(x_m2h)-16.0d0*f(x_m3h)+3.0d0*f(x_m4h))/(12.0d0*h(a))
    end function backward_difference

! -----------------------------------------------------------------------------
! ----- Boundary conditions ---------------------------------------------------
! -----------------------------------------------------------------------------

    function applied_field(i) ! DONE
        integer(kind=8), intent(in) :: i
        real(kind=8), dimension(3) :: grid_coord
        real(kind=8) :: applied_field, r

        grid_coord = cart_coord(i)
        r = emitters_pos(3,1) + d - grid_coord(3) + min_d

        applied_field = -V_d / r
    end function applied_field

    function charge_voltage(charge_pos, i) ! DONE
        ! Calculate the boundary value for an electron
        real(kind=8), dimension(3), intent(in) :: charge_pos
        integer(kind=8), intent(in) :: i
        real(kind=8), dimension(3) :: grid_coord
        real(kind=8) :: charge_voltage, r
        
        grid_coord = cart_coord(i)
        r = norm2(grid_coord - charge_pos) + min_d

        charge_voltage = -q_0 / (4.0d0*pi*epsilon_0*r)
    end function charge_voltage

    subroutine find_electrons() ! TODO
        ! Aggregate boundary values from electrons
        real(kind=8), dimension(3) :: elec_pos, new_pos, node_pos
        real(kind=8), dimension(2,4) :: weights
        integer(kind=8), dimension(2,4) :: nodes
        real(kind=8) :: weight, weight_sum=0.0d0, new_weight_sum
        integer(kind=8) :: g, i, k, u, l, m
        real(kind=8) :: dx, dy, dz
        real(kind=8) :: wx, wy, wz

        do k=1,nrPart
            if (particles_mask(k) .eqv. .true.) then
                elec_pos = particles_cur_pos(:,k)
                if (is_inside(elec_pos) == 1) then
                    i = disc_coord(elec_pos(1), elec_pos(2), elec_pos(3))

                    ! -----------------------------------------------
                    u = i
                    node_pos = cart_coord(u)

                    dx = elec_pos(1) - node_pos(1)
                    dy = elec_pos(2) - node_pos(2)
                    dz = elec_pos(3) - node_pos(3)

                    wx = dx / hx
                    wy = dy / hy
                    wz = dz / hz

                    weight = wx*wy*wz
                    if ((is_emitter(u) == 1) .or. (is_anode(u) == 1)) then
                        weights(1,1) = 0.0d0
                    else
                        weights(1,1) = weight
                        weight_sum = weight_sum + weight
                    end if
                    nodes(1,1) = u

                    ! --
                    u = move(i,1,1)
                    weight = (1.0d0 - wx)*wy*wz
                    if ((is_emitter(u) == 1) .or. (is_anode(u) == 1)) then
                        weights(1,2) = 0.0d0
                    else
                        weights(1,2) = weight
                        weight_sum = weight_sum + weight
                    end if
                    nodes(1,2) = u

                    ! --
                    u = move(i,1,2)
                    weight = wx*(1.0d0 - wy)*wz
                    if ((is_emitter(u) == 1) .or. (is_anode(u) == 1)) then
                        weights(1,3) = 0.0d0
                    else
                        weights(1,3) = weight
                        weight_sum = weight_sum + weight
                    end if
                    nodes(1,3) = u

                    ! --
                    u = move(move(i,1,2),1,1)
                    weight = (1.0d0 - wx)*(1.0d0 - wy)*wz
                    if ((is_emitter(u) == 1) .or. (is_anode(u) == 1)) then
                        weights(1,4) = 0.0d0
                    else
                        weights(1,4) = weight
                        weight_sum = weight_sum + weight
                    end if
                    nodes(1,4) = u
                    g = gridPoints(u)
                    
                    ! -------------------------------------------------
                    u = move(i,1,3)
                    weight = wx*wy*(1.0d0 - wz)
                    if ((is_emitter(u) == 1) .or. (is_anode(u) == 1)) then
                        weights(2,1) = 0.0d0
                    else
                        weights(2,1) = weight
                        weight_sum = weight_sum + weight
                    end if
                    nodes(2,1) = u

                    ! --
                    u = move(move(i,1,1),1,3)
                    weight = (1.0d0 - wx)*wy*(1.0d0 - wz)
                    if ((is_emitter(u) == 1) .or. (is_anode(u) == 1)) then
                        weights(2,2) = 0.0d0
                    else
                        weights(2,2) = weight
                        weight_sum = weight_sum + weight
                    end if
                    nodes(2,2) = u

                    ! --
                    u = move(move(i,1,2),1,3)
                    weight = wx*(1.0d0 - wy)*(1.0d0 - wz)
                    if ((is_emitter(u) == 1) .or. (is_anode(u) == 1)) then
                        weights(2,3) = 0.0d0
                    else
                        weights(2,3) = weight
                        weight_sum = weight_sum + weight
                    end if
                    nodes(2,3) = u

                    ! --
                    u = move(move(move(i,1,2),1,1),1,3)
                    weight = (1.0d0 - wx)*(1.0d0 - wy)*(1.0d0 - wz)
                    if ((is_emitter(u) == 1) .or. (is_anode(u) == 1)) then
                        weights(2,4) = 0.0d0
                    else
                        weights(2,4) = weight
                        weight_sum = weight_sum + weight
                    end if
                    nodes(2,4) = u

                    ! -------------------------------------------------

                    do l=1,2
                        do m=1,4
                            weights(l,m) = weights(l,m) / weight_sum
                            if (weights(l,m) < 0.0d0) then
                                print*, 'Laplace: negative weight:', weights(l,m)
                            end if
                            g = gridPoints(nodes(l,m))
                            iCharge(g) = iCharge(g) - weights(l,m)*particles_charge(k)/(hx*hy*hz*epsilon_0)
                        end do
                    end do

                end if
            end if
        end do
    end subroutine find_electrons

! -----------------------------------------------------------------------------
! ----- Indexing functions ----------------------------------------------------
! -----------------------------------------------------------------------------

    function disc_coord(x,y,z) ! DONE
        ! Change cartesian coordinates to discrete coordinates
        real(kind=8), intent(in) :: x, y, z
        integer(kind=8) :: x_coord,y_coord,z_coord,disc_coord

        x_coord = (x-laplace_pos(1)) / hx
        y_coord = (y-laplace_pos(2)) / hy
        z_coord = (z-laplace_pos(3)) / hz

        ! print *, 'Laplace: disc_coord: Lx=',laplace_dim(1), ' Ly=',laplace_dim(2), ' Lz=',laplace_dim(3)
        ! print *, 'Laplace: disc_coord: hx=',hx, ' hy=',hy, ' hz=',hz
        ! print *, 'Laplace: disc_coord: Nx=',Nx, ' y=',Ny, ' z=',Nz	

        disc_coord = x_coord + y_coord*Nx + z_coord*NxNy + 1
    end function disc_coord

    function cart_coord(i) ! DONE
        ! Change discrete coordinates to cartesian coordinates
        integer(kind=8), intent(in) :: i
        integer(kind=8) :: k, x_coord, y_coord, z_coord
        real(kind=8), dimension(3) :: cart_coord
        k = i-1

        x_coord = mod(k,Nx)
        y_coord = mod(k,NxNy) / Nx
        z_coord = k / NxNy

        cart_coord(1) = x_coord*hx + laplace_pos(1)
        cart_coord(2) = y_coord*hy + laplace_pos(2)
        cart_coord(3) = z_coord*hz + laplace_pos(3)
    end function cart_coord

    function move(i,nr_steps,axis)
        integer(kind=8), intent(in) :: i, nr_steps, axis
        integer(kind=8) :: move

        select case (axis)
            case (1)
                move = i + nr_steps
            case (2)
                move = i + nr_steps*Nx
            case (3)
                move = i + nr_steps*NxNy
        end select
    end function move

    function is_boundary(i) ! DONE
        ! Returns an array with the boundary type for each axis. 0 = innerpoint, 1 = startpoint, 2 = endpoint
        integer(kind=8), intent(in) :: i
        integer(kind=8), dimension(3) :: is_boundary

        is_boundary = (/0,0,0/)
        if (mod(i-1,Nx) == 0) then
            is_boundary(1) = 1
        else if (mod(i-1,Nx) == Nx-1) then
            is_boundary(1) = 2
        end if

        if (mod(i-1,NxNy)/Nx == 0) then
            is_boundary(2) = 1
        else if (mod(i-1,NxNy)/Nx == Ny-1) then
            is_boundary(2) = 2
        end if

        if ((i-1)/NxNy == 0) then
            is_boundary(3) = 1
        else if ((i-1)/NxNy == Nz-1) then
            is_boundary(3) = 2
        end if

    end function is_boundary

    function is_anode(i)
        integer(kind=8), intent(in) :: i
        integer(kind=8) :: is_anode
        real(kind=8), dimension(3) :: point_coord

        is_anode = 0

        point_coord = cart_coord(i)

        if (point_coord(3) >= box_dim(3)) then
            is_anode = 1
        end if

    end function is_anode

    function is_emitter(i)
        ! Returns 1 if point is within the emitter and 0 otherwise
        integer(kind=8), intent(in) :: i
        integer(kind=8) :: is_emitter
        real(kind=8), dimension(3) :: point_coord
        real(kind=8) :: point_radius

        is_emitter = 0

        point_coord = cart_coord(i)

        if (point_coord(3) <= emitters_pos(3,1)) then
            select case (emitters_type(1))
                case (1) ! Circle
                    point_radius = sqrt(point_coord(1)**2 + point_coord(2)**2)
                    if (point_radius <= emitter_radius) then
                        is_emitter = 1
                    else
                        is_emitter = 0
                    end if

                case (2) ! Rectangle
                    if (emitters_pos(1,1) <= point_coord(1) .and. point_coord(1) <= emitters_pos(1,1)+emitters_dim(1,1) &
                        .and. emitters_pos(2,1) <= point_coord(2) .and. point_coord(2) <= emitters_pos(2,1)+emitters_dim(2,1)) then
                        is_emitter = 1
                    else
                        is_emitter = 0
                    end if
                case (4) ! Ring
                    point_radius = sqrt(point_coord(1)**2 + point_coord(2)**2)
                    if ((point_radius <= emitters_ring(1,1)) .and. (point_radius >= emitters_ring(2,1))) then
                        is_emitter = 1
                    else
                        is_emitter = 0
                    end if
            end select
        end if

    end function

    function is_first_layer(i)
        ! Returns 1 if point is within the first layer of the emitter and 0 otherwise
        integer(kind=8), intent(in) :: i
        integer(kind=8) :: axis, is_first_layer, nr_surrounding


        nr_surrounding = 0
        is_first_layer = 0

        if (is_emitter(i) == 1) then
            do axis=1,3
                if (is_emitter(move(i,1,axis)) == 1) then
                    nr_surrounding = nr_surrounding + 1
                end if
                if (is_emitter(move(i,-1,axis)) == 1) then
                    nr_surrounding = nr_surrounding + 1
                end if
            end do

            if (nr_surrounding < 6) then
                is_first_layer = 1
            end if
        end if

    end function is_first_layer

    function is_inside(pos)
        real(kind=8), dimension(3), intent(in) :: pos
        integer(kind=8) :: is_inside

        if (lim_x(1) <= pos(1) .and. pos(1) <= lim_x(2) &
            .and. lim_y(1) <= pos(2) .and. pos(2) <= lim_y(2) &
            .and. lim_z(1) <= pos(3) .and. pos(3) <= lim_z(2)) then
            is_inside = 1
        else
            is_inside = 0
        end if
    end function

! -----------------------------------------------------------------------------
! ----- Clean Up --------------------------------------------------------------
! -----------------------------------------------------------------------------
    
    subroutine Clean_Up_Laplace() ! DONE
        integer(kind=8) :: error, phase

        ! Clean up memory
        print '(a)', 'Laplace: Clean up'

        ! Clear pardiso memory
        phase = -1
        if (.not. allocated(perm)) then
            allocate(perm(nrGridActive))
        end if

        print *, '  Laplace: clearing pardiso memory'
        call pardiso_64(pt,maxfct,mnum,mtype,phase,nrGridActive,nnz_values,nnz_ia,nnz_col_index,perm,nrhs,iparm,msglvl,b,x,error)
        
        ! call pardiso_32(pt,maxfct,mnum,mtype,phase,nrGridActive,nnz_values,nnz_ia,nnz_col_index,perm,nrhs,iparm,msglvl,b,x,error)
        if (error == 0) then
            print *, '  Laplace: clearing pardiso memory successful'
        else
            print *, '  Laplace: clearing pardiso memory failed with error code ', error
        end if
        ! Deallocate CSR matrix
        deallocate(nnz_values)
        deallocate(nnz_col_index)
        deallocate(nnz_ia)
        deallocate(nnz_center)
        ! Deallocate other arrays
        deallocate(b)
        deallocate(oCharge)
        deallocate(gridPoints)
        deallocate(gridPointsActive)
        deallocate(x)
        deallocate(perm)

        ! These should have been deallocated already in the last iteration
        if (allocated(iCharge)) then 
            deallocate(iCharge)
        end if
        if (allocated(laplace_field)) then
            deallocate(laplace_field) 
        end if
        if (allocated(voltage)) then
            deallocate(voltage)
        end if
    end subroutine Clean_Up_Laplace

    subroutine check_matrix()
        integer(kind=8) :: i,count,index,end

        do i=1,nrGridActive
            print*, 'Laplace: row=',i
            count = 0
            index = nnz_ia(i)
            end = nnz_ia(i+1)
            do while(index < end)
                print*,'    Laplace: column=',nnz_col_index(index)
                index = index + 1
                count = count + 1
            end do
        end do

    end subroutine check_matrix

    subroutine Place_Electron(step)
        integer(kind=8), intent(in) :: step
        real(kind=8), dimension(3) :: par_pos, par_vel

        if (step == 1) then
            par_pos(1) = 0.0d0*length_scale
            par_pos(2) = 0.0d0*length_scale
            par_pos(3) = 1.05d0*length_scale 
            par_vel = 0.0d0
            call Add_Particle(par_pos,par_vel,species_elec,step,1,-1)
            ! par_pos(1) = 5.0d0*length_scale
            ! par_pos(2) = 5.0d0*length_scale
            ! par_pos(3) = 1.05d0*length_scale 
            ! par_vel = 0.0d0
            ! call Add_Particle(par_pos,par_vel,species_elec,step,1,-1)
        end if
    end subroutine Place_Electron

    subroutine Write_Laplace_Data()
        integer :: ud_lp_field, ud_ic_field, IFAIL, lvl, i, j, k
        character(len=128) :: filename_lp_field, filename_ic_field
        real(kind=8), dimension(3) :: cur_pos, ic_field, lp_field
        
        do lvl=0,2
            ! write(filename_lp_voltage, '(a15,i0,a4)') 'out/lp_voltage_',lvl,'.bin'
            write(filename_lp_field, '(a13,i0,a4)') 'out/lp_field_',lvl,'.bin'
            write(filename_ic_field, '(a13,i0,a4)') 'out/ic_field_',lvl,'.bin'

            ! ! Open the voltage file
            ! open(newunit=ud_lp_voltage, iostat=IFAIL, file=filename_lp_voltage, status='REPLACE', action='WRITE', access='STREAM')
            ! if (IFAIL /= 0) then
            !     print *, 'RUMDEED: Failed to open the laplace voltage file.'
            !     return 
            ! end if

            ! Open the field file
            open(newunit=ud_lp_field, iostat=IFAIL, file=filename_lp_field, status='REPLACE', action='WRITE', access='STREAM')
            if (IFAIL /= 0) then
                print *, 'RUMDEED: Failed to open the laplace field file.'
                return 
            end if

            ! Open the field file
            open(newunit=ud_ic_field, iostat=IFAIL, file=filename_ic_field, status='REPLACE', action='WRITE', access='STREAM')
            if (IFAIL /= 0) then
                print *, 'RUMDEED: Failed to open the image charge field file.'
                return 
            end if

            do i=0,99
                do j=0,99
                    cur_pos(1) = i*laplace_dim(1)/99 + laplace_pos(1)
                    cur_pos(2) = j*laplace_dim(2)/99 + laplace_pos(2)
                    cur_pos(3) = 0.0d0

                    lp_field = Calculate_Laplace_Field_at(cur_pos)
                    ic_field = Calc_Field_at(cur_pos)
                    ! print *, 'Laplace: writing: x=',test_pos(1), ' y=',test_pos(2), ' z=',test_pos(3)
                    ! write(unit=ud_lp_voltage,iostat=IFAIL) i, j, voltage(k), is_emitter(k)
                    write(unit=ud_lp_field,iostat=IFAIL) i, j, lp_field(3), is_emitter(k)
                    write(unit=ud_ic_field,iostat=IFAIL) i, j, ic_field(3), is_emitter(k)
                end do
            end do

            ! close(unit=ud_lp_voltage, iostat=IFAIL, status='keep')
            close(unit=ud_lp_field, iostat=IFAIL, status='keep')
            close(unit=ud_ic_field, iostat=IFAIL, status='keep')
        end do
        
    end subroutine Write_Laplace_Data

    ! subroutine Sample_Field(step)
    !     integer, intent(in) :: step
    !     integer             :: IFAIL, x, y, ud_field_pos, Nx, Ny
    !     character(len=128)  :: filename
    !     real(kind=precision)    :: hx, hy, cur_rad
    !     real(kind=precision), dimension(3) :: cur_pos, cur_field


    !     if (sample_field_file .eqv. .true.) then
    !       if (mod(step,sample_field_rate) == 0) then
    !         ! Create the file
    !         write(filename, '(a10, i0, a4)') 'out/field-', step, '.bin'
    !         ! Open the file
    !         open(newunit=ud_field_pos, iostat=IFAIL, file=filename, status='REPLACE', action='WRITE', access='STREAM')
    !         if (IFAIL /= 0) then
    !         print *, 'RUMDEED: Failed to open the field position file.'
    !         return 
    !         end if

    !         Nx = 1000
    !         Ny = 1000
    !         hx = 2.0d0*emitters_dim(1,1) / Nx
    !         hy = 2.0d0*emitters_dim(2,1) / Ny

    !         do x=1,Nx
    !           do y=1,Ny
    
    !             cur_pos(1) = emitters_pos(1,1) + (x-1)*hx
    !             cur_pos(2) = emitters_pos(2,1) + (y-1)*hy
    !             cur_pos(3) = 0.0d0

    !             cur_rad = sqrt(cur_pos(1)**2 + cur_pos(2)**2)

    !             if (emitters_ring(2,1) <= cur_rad .and. cur_rad <= emitters_ring(1,1)) then
    
    !                 if (laplace .eqv. .true.) then
    !                     cur_field = Calculate_Laplace_Field_at(cur_pos)
    !                 else
    !                     cur_field = Calc_Field_at(cur_pos)
    !                 end if

    !                 ! print *, cur_field(3)

    !                 ! Write data
    !                 write(unit=ud_field_pos, iostat=IFAIL) cur_pos(:), cur_field(3)
    !             end if
    !           end do
    !         end do 

    !         ! Close the file
    !         close(ud_field_pos, iostat=IFAIL)
    !       end if
    !     end if
    ! end subroutine Sample_Field

    function Calculate_Laplace_Field_at(pos)
        real(kind=8), dimension(3), intent(in) :: pos
        real(kind=8), dimension(3) :: Calculate_Laplace_Field_at
        real(kind=8), dimension(3) :: pos_bot1, pos_bot2, pos_bot3, pos_bot4, pos_bot5, pos_bot6
        real(kind=8), dimension(3) :: pos_top1, pos_top2, pos_top3, pos_top4, pos_top5, pos_top6
        real(kind=8), dimension(3) :: field_bot1, field_bot2, field_bot3, field_bot4, field_top1, field_top2, field_top3, field_top4
        real(kind=8), dimension(3) :: field_bot5, field_bot6, field_bot, field_top5, field_top6, field_top
        integer(kind=8) :: i,k

        if (is_inside(pos) == 0) then
            print *, 'Laplace: calculating field outside the domain'
            return
        end if

        i = disc_coord(pos(1), pos(2), pos(3))

        ! Bottom side
        pos_bot1 = cart_coord(i)
        field_bot1 = laplace_field(i,:)

        ! if (pos_bot1(1) == -5.0d-9) then
        !     print *, 'Laplace: k=',i
        !     print *, 'Laplace: pos_bot1=',pos_bot1
        ! end if

        k = move(i,1,1)
        pos_bot2 = cart_coord(k)
        field_bot2 = laplace_field(k,:)

        k = move(i,1,2)
        pos_bot3 = cart_coord(k)
        field_bot3 = laplace_field(k,:)

        ! if (pos_bot1(1) == -5.0d-9) then
        !     print *, 'Laplace: k=',k
        !     print *, 'Laplace: pos_bot3=',pos_bot3
        ! end if

        k = move(move(i,1,2),1,1)
        pos_bot4 = cart_coord(k)
        field_bot4 = laplace_field(k,:)

        ! Top side
        k = move(i,1,3)
        pos_top1 = cart_coord(k)
        field_top1 = laplace_field(k,:)

        k = move(move(i,1,3),1,1)
        pos_top2 = cart_coord(k)
        field_top2 = laplace_field(k,:)

        k = move(move(i,1,3),1,2)
        pos_top3 = cart_coord(k)
        field_top3 = laplace_field(k,:)

        k = move(move(move(i,1,3),1,2),1,1)
        pos_top4 = cart_coord(k)
        field_top4 = laplace_field(k,:)

        ! Bottom interpolation
        field_bot5 = linear_interpolate_3d(pos_bot1(1),pos_bot2(1),pos(1),field_bot1,field_bot2)
        field_bot6 = linear_interpolate_3d(pos_bot3(1),pos_bot4(1),pos(1),field_bot3,field_bot4)
        field_bot = linear_interpolate_3d(pos_bot1(2),pos_bot3(2),pos(2),field_bot5,field_bot6)

        ! Top interpolation
        field_top5 = linear_interpolate_3d(pos_top1(1),pos_top2(1),pos(1),field_top1,field_top2)
        field_top6 = linear_interpolate_3d(pos_top3(1),pos_top4(1),pos(1),field_top3,field_top4)
        field_top = linear_interpolate_3d(pos_top1(2),pos_top3(2),pos(2),field_top5,field_top6)

        ! Final interpolation
        Calculate_Laplace_Field_at = linear_interpolate_3d(pos_bot1(3),pos_top1(3),pos(3),field_bot,field_top)

        ! if ((pos(1) /= 5.0d-9) .and. (pos(2) /= 5.0d-9)) then
        ! if (isnan(Calculate_Laplace_Field_at(3)) .eqv. .true.) then

        !     print *, ''
        !     print *, 'Laplace: pos_bot1=',pos_bot1
        !     print *, 'Laplace: pos_bot2=',pos_bot2
        !     print *, 'Laplace: pos_bot3=',pos_bot3
        !     print *, 'Laplace: pos_bot4=',pos_bot4
        !     print *, ''
        !     print *, 'Laplace: pos=',pos
        !     print *, 'Laplace: field_bot1=',field_bot1(3)
        !     print *, 'Laplace: field_bot2=',field_bot2(3)
        !     print *, 'Laplace: field_bot3=',field_bot3(3)
        !     print *, 'Laplace: field_bot4=',field_bot4(3)
        !     print *, 'Laplace: field_bot5=',field_bot5(3)
        !     print *, 'Laplace: field_bot6=',field_bot6(3)
        !     print *, 'Laplace: field_bot=',field_bot(3)

        !     print *, ''
        !     print *, 'Laplace: pos_top1=',pos_top1
        !     print *, 'Laplace: pos_top2=',pos_top2
        !     print *, 'Laplace: pos_top3=',pos_top3
        !     print *, 'Laplace: pos_top4=',pos_top4

        !     print *, ''

        !     print *, 'Laplace: field_top1=',field_top1(3)
        !     print *, 'Laplace: field_top2=',field_top2(3)
        !     print *, 'Laplace: field_top3=',field_top3(3)
        !     print *, 'Laplace: field_top4=',field_top4(3)
        !     print *, 'Laplace: field_top5=',field_top5(3)
        !     print *, 'Laplace: field_top6=',field_top6(3)
        !     print *, 'Laplace: field_top=',field_top(3)

        !     print *, 'Laplace: Calculate_Laplace_Field_at: NaN'
        ! end if
        ! end if

        ! print*, 'Laplace: field_top1=',field_top1(3)
        ! print*, 'Laplace: field_top2=',field_top2(3)
        ! print*, 'Laplace: field_top3=',field_top3(3)
        ! print*, 'Laplace: field_top4=',field_top4(3)

        ! print*, 'Laplace: field_bot1=',field_bot1(3)
        ! print*, 'Laplace: field_bot2=',field_bot2(3)
        ! print*, 'Laplace: field_bot3=',field_bot3(3)
        ! print*, 'Laplace: field_bot4=',field_bot4(3)

        ! print*, 'Laplace: Calculate_Laplace_Field_at=',Calculate_Laplace_Field_at(3)

    end function Calculate_Laplace_Field_at

    function linear_interpolate_3d(x1,x2,xi,f1,f2)
        real(kind=8), intent(in) :: x1,x2,xi
        real(kind=8), dimension(3), intent(in) :: f1, f2
        real(kind=8), dimension(3) :: linear_interpolate_3d
        integer(kind=8) :: i

        do i=1,3
            linear_interpolate_3d(i) = f1(i) + (f2(i)-f1(i))/(x2-x1+min_d)*(xi-x1)
        end do

    end function linear_interpolate_3d

end module mod_laplace_solver