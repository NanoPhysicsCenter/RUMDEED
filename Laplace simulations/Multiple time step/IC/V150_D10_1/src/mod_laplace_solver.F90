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
    use mkl_pardiso ! Required for pardiso solver
    use mkl_spblas  ! Required for sparse matrices
    ! use iso_c_binding, only: c_int, c_ptr, c_f_pointer

    implicit none

    ! Accessibility
    private
    public :: Init_Laplace_Solver, Calculate_Laplace_Field, Calculate_Laplace_Field_at, Clean_Up_Laplace, Place_Electron, Write_Laplace_Data

    ! Pardiso variables
    type(MKL_PARDISO_HANDLE), dimension(64) :: pt
    integer(kind=8), dimension(64) :: iparm

    integer(kind=8), parameter :: maxfct=1, mnum=1, mtype=11, nrhs=1, msglvl=1
    integer(kind=8), parameter :: analysis=11, factorization=22, solving=33, solving1=331, solving2=333

    integer(kind=8) :: ndiff
    integer(kind=8), allocatable, dimension(:) :: perm, permdiff

    ! Matrix
    double precision, allocatable, dimension(:) :: nnz_values
    integer(kind=8), allocatable, dimension(:) :: nnz_col_index, nnz_ia, nnz_center
    integer(kind=8) :: nnz

    ! System
    double precision, allocatable, dimension(:) :: b
    double precision, allocatable, dimension(:) :: iCharge, oCharge

    ! Solution
    double precision, allocatable, dimension(:) :: voltage, x
    double precision, allocatable, dimension(:,:) :: laplace_field

    ! Discretization
    integer(kind=8), allocatable, dimension(:) :: gridPoints, gridPointsActive
    double precision :: hx, hy, hz, emitter_radius
    double precision, dimension(3) :: div_h, div_h2, div_h3
    double precision, dimension(2) :: lim_x, lim_y, lim_z
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

    subroutine Calculate_Laplace_Field(step) ! TODO
        integer(kind=8), intent(in) :: step
        ! Calculate the electric field
        print *, 'Laplace: finding electrons'
        call find_electrons()
        print *, 'Laplace: updating matrix'
        call update_matrix()
        print *, 'Laplace: solving matrix'
        call solve_system()
        print *, 'Laplace: allocate voltage'
        call allocate_voltage()
        print *, 'Laplace: calculating field'
        call calculate_field()
        print *, 'Laplace: done'

        ! call write_average_field()
        
    end subroutine Calculate_Laplace_Field

! ------------------------------------------------------------------------------
! ----- Solution ---------------------------------------------------------------
! ------------------------------------------------------------------------------

    subroutine solve_system() ! TODO
        ! Solve the system
        integer(kind=8) :: i
        double precision :: start, end

        if (allocated(laplace_field)) then
            deallocate(laplace_field)
        end if

        if (ndiff /= 0) then
            ! iparm(4) = 31+
            ! iparm(9) = 12 ! tolerance = 10^(-iparm(9))	
            ! iparm(11) = 0
            ! iparm(13) = 0
            ! iparm(24) = 10
            call pardiso_phase(11) 
            ! iparm(39) = 1
            call pardiso_phase(22)
        end if
        call pardiso_phase(33)

        ! print *, 'Laplace: iterative solving time:', end-start

        ! print*, 'Laplace: solve_system', iparm(20)

    end subroutine solve_system

! -------------------------------------------------------------------------------------------------------------
! ----- Initialization ----------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------------------

    subroutine init_grid() ! DONE
        ! Initialize parameters and store points that are outside of the emitter or one layer within it
        integer(kind=8) :: i, a
        integer(kind=8), dimension(3) :: boundaries
        double precision :: Lx, Ly, Lz

        Nx = laplace_intervals(1); Ny = laplace_intervals(2); Nz = laplace_intervals(3)

        Lx = laplace_dim(1); Ly = laplace_dim(2); Lz = laplace_dim(3)

        NxNy = Nx*Ny

        nrGrid = Nx*Ny*Nz

        hx = Lx / (Nx-1)
        hy = Ly / (Ny-1)
        hz = Lz / (Nz-1)

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

        allocate(oCharge(nrGridActive), x(nrGridActive))

    end subroutine init_grid

    subroutine init_matrix() ! DONE
        ! Initialize matrix without boundary conditions
        integer(kind=8) :: g, i, a, center, next_center
        integer(kind=8), dimension(3) :: boundaries
        integer(kind=8), dimension(3,3) :: indices
        integer(kind=8) :: next, next_next, next_next_next, prev, prev_prev, prev_prev_prev

        allocate(nnz_values(nnz), nnz_col_index(nnz), nnz_ia(nrGridActive+1))
        allocate(b(nrGridActive), nnz_center(nrGridActive))

        nnz_values = nnz_values*0.0d0
        b = b*0.0d0

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
        double precision :: start, end

        ! Initialize pardiso solver and do analysis, symbolic factorization, and numerical factorization
        call pardisoinit(pt,mtype,iparm)
        ! Solver parameters
        iparm(1) = 1 ! use user-defined parameters
        iparm(2) = 3 ! 1 = minimum degree algorithm, 2 = metis, 3 = parallel metis
        ! iparm(3) = 0 ! reserved
        iparm(4) = 0 ! >0 = preconditioned CGS
        iparm(5) = 0 ! 0 = ignored; 1 = user supplied fill-in permutation; 2 = returns permutation
        iparm(6) = 0 ! 0 = write solution on x; 1 = write solution on b
        ! iparm(7) = 0 ! output number of iterative refinement steps
        iparm(8) = 100 ! maximum number of iterative refinement steps
        iparm(9) = 6 ! tolerance = 10^(-iparm(9))	
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
        iparm(28) = 0 ! 0 = double precision; 1 = single precision
        ! iparm(29) = 0 ! reserved
        ! iparm(30) = 0 ! output
        iparm(31) = 0 ! 0 = no partial solve
        ! iparm(32) = 0 ! reserved
        ! iparm(33) = 0 ! reserved
        ! iparm(34) = 0 ! 0 = CNR mode and oneapi determines optimal number of threads; 1 = CNR mode and user determines number of threads; 2 = CNR mode and user determines number of threads; 3 = CNR mode and user determines number of threads
        iparm(35) = 0 ! 0 = one-based indexing; 1 = zero-based indexing
        iparm(36) = 0 ! 0 = not compute shur; 1 = compute shur
        iparm(37) = 0 ! 0 = CSR format; 1 = BSR format
        ! iparm(38) = 0 ! reserved
        iparm(39) = 1 ! 0 = full factorization; 1 = low-rank update
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
        call pardiso_phase(11)
        iparm(27) = 0
        call pardiso_phase(22)
        call pardiso_phase(33)
        ! call cpu_time(end)

        ! print *, 'Laplace: pardiso initialization time:', end-start

        ! deallocate(perm)

        
    end subroutine init_pardiso

    subroutine pardiso_phase(phase)
        integer(kind=8), intent(in) :: phase
        double precision, dimension(nrGridActive) :: ddum
        integer(kind=8) :: error

        ! print *, '  Laplace: analysis and symbolic factorization'
        call pardiso_64(pt,maxfct,mnum,mtype,phase,nrGridActive,nnz_values,nnz_ia,nnz_col_index,perm,nrhs,iparm,msglvl,b,x,error)
        if (error == 0) then
            ! print *, '  Laplace: analysis and symbolic factorization successful'
        else
            print *, '  Laplace: phase', phase, 'failes with error code ', error
        end if
        ! print *, '  Laplace:', iparm(17)*1e-6, 'Gb needed for numerical factorization'

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
                voltage(i) = x(g)
            end if
        end do

        ! $OMP END PARALLEL DO

    end subroutine allocate_voltage

    subroutine calculate_field() ! DONE & PARALLEL
        integer(kind=8) :: i, a
        integer(kind=8), dimension(3) :: boundaries

        allocate(laplace_field(nrGrid,3))

        laplace_field = 0.0d0

        ! $OMP PARALLEL DO DEFAULT(NONE) &
        ! $OMP& SHARED(nrGrid, voltage, laplace_field) &
        ! $OMP& PRIVATE(i, boundaries, boundary)

        do i=1,nrGrid
            boundaries = is_boundary(i)
            do a=1,3
                select case (boundaries(a))
                    case (0) ! Inner point
                        laplace_field(i,a) = - central_difference(voltage,i,a)
                        if (a==3) then
                            laplace_field(i,a) = laplace_field(i,a) + applied_field(i)
                        end if
                    case (1) ! Startpoint
                        laplace_field(i,a) = - forward_difference(voltage,i,a)
                        if (a==3) then
                            laplace_field(i,a) = laplace_field(i,a) + applied_field(i)
                        end if
                    case (2) ! Endpoint
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
        double precision :: sum, average_field

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
        double precision :: avg_field

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
        integer(kind=8) :: g

        ! allocate(permdiff(nrGridActive*9))
        ndiff = 0

        do g=1,nrGridActive
            if (oCharge(g) /= iCharge(g)) then
                if (iCharge(g) == 0.0d0) then
                    call remove_electron_boundary(g)
                    call insert_finite_difference(g)
                else if (oCharge(g) == 0.0d0) then
                    call remove_finite_difference(g)
                    call insert_electron_boundary(g)
                else
                    call remove_electron_boundary(g)
                    call insert_electron_boundary(g)
                    ! if (allocated(permdiff)) then
                    !     ndiff = ndiff + 1
                    !     permdiff((ndiff-1)*2+1) = g
                    !     permdiff((ndiff-1)*2+2) = g
                    ! else
                    !     ndiff = ndiff + 1
                    ! end if
                end if

                oCharge(g) = iCharge(g)
            end if
        end do

        ! print *, 'Laplace: number of changes:', ndiff

        deallocate(iCharge)

        ! allocate(perm(1+2*ndiff))
        ! call reduce_perm()
        ! deallocate(permdiff)

        ! print *, 'Laplace: perm:', perm(2:ndiff*2+1)

    end subroutine update_matrix

    subroutine insert_electron_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = nnz_center(g)

        nnz_values(center) = 1.0d0
        b(g) = iCharge(g)
    end subroutine insert_electron_boundary

    subroutine remove_electron_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = nnz_center(g)

        nnz_values(center) = 0.0d0
        b(g) = 0.0d0
    end subroutine remove_electron_boundary

    subroutine insert_emitter_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = nnz_center(g)

        nnz_values(center) = 1.0d0
        b(g) = 0.0d0
    end subroutine insert_emitter_boundary

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

                    if (allocated(permdiff)) then
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = prev
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = next
                    else
                        ndiff = ndiff + 2
                    end if
                case (1) ! Startpoint (values for forward difference)
                    next = indices(a,1)
                    next_next = indices(a,2)
                    next_next_next = indices(a,3)

                    nnz_values(center) = nnz_values(center) + 2.0d0*div_h2(a)
                    nnz_values(next) = nnz_values(next) - 5.0d0*div_h2(a)
                    nnz_values(next_next) = nnz_values(next_next) + 4.0d0*div_h2(a)
                    nnz_values(next_next_next) = nnz_values(next_next_next) - 1.0d0*div_h2(a)

                    if (allocated(permdiff)) then
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = next
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = next_next
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = next_next_next
                    else
                        ndiff = ndiff + 3
                    end if
                case (2) ! Endpoint (values for backward difference)
                    prev = indices(a,3)
                    prev_prev = indices(a,2)
                    prev_prev_prev = indices(a,1)

                    nnz_values(prev_prev_prev) = nnz_values(prev_prev_prev) - 1.0d0*div_h2(a)
                    nnz_values(prev_prev) = nnz_values(prev_prev) + 4.0d0*div_h2(a)
                    nnz_values(prev) = nnz_values(prev) - 5.0d0*div_h2(a)
                    nnz_values(center) = nnz_values(center) + 2.0d0*div_h2(a)

                    if (allocated(permdiff)) then
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = prev
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = prev_prev
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = prev_prev_prev
                    else
                        ndiff = ndiff + 3
                    end if
            end select
        end do

        if (allocated(permdiff)) then
            ndiff = ndiff + 1
            permdiff((ndiff-1)*2+1) = g
            permdiff((ndiff-1)*2+2) = center
        end if

        b(g) = 0.0d0
    end subroutine insert_finite_difference

    subroutine remove_finite_difference(g) ! DONE
        ! Remove finite difference equation from matrix
        integer(kind=8), intent(in) :: g
        integer(kind=8), dimension(3) :: boundaries
        integer(kind=8), dimension(3,3) :: indices
        integer(kind=8) :: i, a, center
        integer(kind=8) :: prev_prev_prev, prev_prev, prev, next, next_next, next_next_next

        center = nnz_center(g)
        i = gridPointsActive(g)
        boundaries = is_boundary(i)	
        call get_finite_indices(g,indices)

        nnz_values(center) = 0.0d0
        do a=1,3
            select case (boundaries(a))
                case (0)
                    prev = indices(a,1)
                    next = indices(a,3)

                    nnz_values(prev) = 0.0d0
                    nnz_values(next) = 0.0d0

                    if (allocated(permdiff)) then
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = prev
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = next
                    else
                        ndiff = ndiff + 2
                    end if
                case (1)
                    next = indices(a,1)
                    next_next = indices(a,2)
                    next_next_next = indices(a,3)

                    nnz_values(next) = 0.0d0
                    nnz_values(next_next) = 0.0d0
                    nnz_values(next_next_next) = 0.0d0

                    if (allocated(permdiff)) then
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = next
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = next_next
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = next_next_next
                    else
                        ndiff = ndiff + 3
                    end if
                case (2)
                    prev = indices(a,3)
                    prev_prev = indices(a,2)
                    prev_prev_prev = indices(a,1)

                    nnz_values(prev) = 0.0d0
                    nnz_values(prev_prev) = 0.0d0
                    nnz_values(prev_prev_prev) = 0.0d0

                    if (allocated(permdiff)) then
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = prev
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = prev_prev
                        ndiff = ndiff + 1
                        permdiff((ndiff-1)*2+1) = g
                        permdiff((ndiff-1)*2+2) = prev_prev_prev
                    else
                        ndiff = ndiff + 3
                    end if
            end select
        end do

        if (allocated(permdiff)) then
            ndiff = ndiff + 1
            permdiff((ndiff-1)*2+1) = g
            permdiff((ndiff-1)*2+2) = center
        end if

        b(g) = 0.0d0
    end subroutine remove_finite_difference

    subroutine reduce_perm()
        integer(kind=8) :: i

        perm(1) = ndiff
        do i=1,ndiff*2
            perm(i+1) = permdiff(i) - 1
        end do

    end subroutine reduce_perm

! -----------------------------------------------------------------------------
! ----- Finite difference -----------------------------------------------------
! -----------------------------------------------------------------------------

    function central_difference(f,x,a)
        double precision, dimension(nrGrid), intent(in) :: f
        integer(kind=8), intent(in) :: x, a
        integer(kind=8) :: x_m1h, x_p1h
        double precision :: central_difference

        x_m1h = move(x,-1,a)
        x_p1h = move(x,1,a)

        central_difference = (f(x_p1h) - f(x_m1h))*div_h(a)/2.0d0
    end function central_difference

    function forward_difference(f,x,a)
        double precision, dimension(nrGrid), intent(in) :: f
        integer(kind=8), intent(in) :: x, a
        integer(kind=8) :: x_p1h, x_p2h, x_3h, x_4h
        double precision :: forward_difference

        x_p1h = move(x,1,3)
        x_p2h = move(x,2,3)

        forward_difference = (-3.0d0*f(x)+4.0d0*f(x_p1h)-f(x_p2h))*div_h(a)/2.0d0
        ! forward_difference = (-25.0d0*f(x)+48.0d0*f(x_1h)-36*f(x_2h)+16.0d0*f(x_3h)-3.0d0*f(x_4h))/(12.0d0*hz)
    end function forward_difference

    function backward_difference(f,x,a)
        double precision, dimension(nrGrid), intent(in) :: f
        integer(kind=8), intent(in) :: x, a
        integer(kind=8) :: x_m1h, x_m2h
        double precision :: backward_difference

        x_m1h = move(x,-1,a)
        x_m2h = move(x,-2,a)

        backward_difference = (f(x_m2h)-4.0d0*f(x_m1h)+3.0d0*f(x))*div_h(a)/2.0d0
    end function backward_difference

! -----------------------------------------------------------------------------
! ----- Boundary conditions ---------------------------------------------------
! -----------------------------------------------------------------------------

    function applied_field(i) ! DONE
        integer(kind=8), intent(in) :: i
        double precision, dimension(3) :: grid_coord
        double precision :: applied_field

        grid_coord = cart_coord(i)

        applied_field = -V_d / abs(emitters_pos(3,1)+box_dim(3)-grid_coord(3))
    end function applied_field

    function charge_voltage(charge_pos, i) ! DONE
        ! Calculate the boundary value for an electron
        double precision, dimension(3), intent(in) :: charge_pos
        integer(kind=8), intent(in) :: i
        double precision, dimension(3) :: grid_coord
        double precision :: charge_voltage
        
        grid_coord = cart_coord(i)
        charge_voltage = -q_0 / ((4.0d0*pi*epsilon_0) * norm2(grid_coord - charge_pos))
    end function charge_voltage

    subroutine find_electrons() ! DONE
        ! Aggregate boundary values from electrons
        double precision, dimension(3) :: elec_pos, new_pos
        integer(kind=8) :: g, i, k, u
        double precision :: boundary_voltage

        allocate(iCharge(nrGridActive))
        
        iCharge = 0.0d0

        ! $OMP PARALLEL DO DEFAULT(NONE) &
        ! $OMP& SHARED(particles_mask, particles_cur_pos, iCharge, lim_x, lim_y, lim_z, nrPart, gridPoints) &
        ! $OMP& PRIVATE(elec_pos, g, i, k)
        do k=1,nrPart
            if (particles_mask(k) .eqv. .true.) then
                elec_pos = particles_cur_pos(:,k)
                ! print *, 'Laplace: checking a particle at x=',elec_pos(1), ' y=',elec_pos(2), ' z=',elec_pos(3)
                if (is_inside(elec_pos) == 1) then

                    ! print *, 'Laplace: particle position: x=',elec_pos(1), ' y=',elec_pos(2), ' z=',elec_pos(3)
                    i = disc_coord(elec_pos(1), elec_pos(2), elec_pos(3))
                    ! print *, 'Laplace: particle index: i=',i
                    ! new_pos = cart_coord2(i)
                    ! print *, 'Laplace: new particle position: x=',new_pos(1), ' y=',new_pos(2), ' z=',new_pos(3)

                    ! Bottom side boundary values
                    u = i
                    g = gridPoints(u)
                    if (g > nrGridActive) then
                        print *, 'Laplace: find_electrons: g > nrGridActive'
                        print *, 'Laplace: find_electrons: pos =',elec_pos
                    end if
                    if (is_emitter(u) == 0) then
                        if (g /= 0) then 
                            ! new_pos = cart_coord(u)
                            ! print *, 'Laplace: adding boundary: x=',new_pos(1), ' y=',new_pos(2), ' z=',new_pos(3)
                            ! $OMP CRITICAL
                            boundary_voltage = charge_voltage(elec_pos,u) 
                            ! print *, 'Laplace: adding boundary voltage:', boundary_voltage
                            iCharge(g) = iCharge(g) + boundary_voltage
                            ! $OMP END CRITICAL
                        end if
                    end if

                    u = move(i,1,1)
                    g = gridPoints(u)
                    if (g > nrGridActive) then
                        print *, 'Laplace: find_electrons: g > nrGridActive'
                        print *, 'Laplace: find_electrons: pos =',elec_pos
                    end if
                    if (is_emitter(u) == 0) then
                        if (g /= 0) then 
                            ! new_pos = cart_coord(u)
                            ! print *, 'Laplace: adding boundary: x=',new_pos(1), ' y=',new_pos(2), ' z=',new_pos(3)
                            ! $OMP CRITICAL
                            boundary_voltage = charge_voltage(elec_pos,u) 
                            ! print *, 'Laplace: adding boundary voltage:', boundary_voltage
                            iCharge(g) = iCharge(g) + boundary_voltage
                            ! $OMP END CRITICAL
                        end if
                    end if

                    u = move(i,1,2)
                    g = gridPoints(u)
                    if (g > nrGridActive) then
                        print *, 'Laplace: find_electrons: g > nrGridActive'
                        print *, 'Laplace: find_electrons: pos =',elec_pos
                    end if
                    if (is_emitter(u) == 0) then
                        if (g /= 0) then 
                            ! new_pos = cart_coord(u)
                            ! print *, 'Laplace: adding boundary: x=',new_pos(1), ' y=',new_pos(2), ' z=',new_pos(3)
                            ! $OMP CRITICAL
                            boundary_voltage = charge_voltage(elec_pos,u) 
                            ! print *, 'Laplace: adding boundary voltage:', boundary_voltage
                            iCharge(g) = iCharge(g) + boundary_voltage
                            ! $OMP END CRITICAL
                        end if
                    end if

                    u = move(move(i,1,2),1,1)
                    g = gridPoints(u)
                    if (g > nrGridActive) then
                        print *, 'Laplace: find_electrons: g > nrGridActive'
                        print *, 'Laplace: find_electrons: pos =',elec_pos
                    end if
                    if (is_emitter(u) == 0) then
                        if (g /= 0) then 
                            ! new_pos = cart_coord(u)
                            ! print *, 'Laplace: adding boundary: x=',new_pos(1), ' y=',new_pos(2), ' z=',new_pos(3)
                            ! $OMP CRITICAL
                            boundary_voltage = charge_voltage(elec_pos,u) 
                            ! print *, 'Laplace: adding boundary voltage:', boundary_voltage
                            iCharge(g) = iCharge(g) + boundary_voltage
                            ! $OMP END CRITICAL
                        end if
                    end if
                    
                    ! Top side boundary values
                    u = move(i,1,3)
                    g = gridPoints(u)
                    if (g > nrGridActive) then
                        print *, 'Laplace: find_electrons: g > nrGridActive'
                        print *, 'Laplace: find_electrons: pos =',elec_pos
                    end if
                    if (is_emitter(u) == 0) then
                        if (g /= 0) then 
                            ! new_pos = cart_coord(u)
                            ! print *, 'Laplace: adding boundary: x=',new_pos(1), ' y=',new_pos(2), ' z=',new_pos(3)
                            ! $OMP CRITICAL
                            boundary_voltage = charge_voltage(elec_pos,u) 
                            ! print *, 'Laplace: adding boundary voltage:', boundary_voltage
                            iCharge(g) = iCharge(g) + boundary_voltage
                            ! $OMP END CRITICAL
                        end if
                    end if

                    u = move(move(i,1,1),1,3)
                    g = gridPoints(u)
                    if (g > nrGridActive) then
                        print *, 'Laplace: find_electrons: g > nrGridActive'
                        print *, 'Laplace: find_electrons: pos =',elec_pos
                    end if
                    if (is_emitter(u)==0) then
                        if (g /= 0 ) then 
                            ! new_pos = cart_coord(u)
                            ! print *, 'Laplace: adding boundary: x=',new_pos(1), ' y=',new_pos(2), ' z=',new_pos(3)
                            ! $OMP CRITICAL
                            boundary_voltage = charge_voltage(elec_pos,u) 
                            ! print *, 'Laplace: adding boundary voltage:', boundary_voltage
                            iCharge(g) = iCharge(g) + boundary_voltage
                            ! $OMP END CRITICAL
                        end if
                    end if

                    u = move(move(i,1,2),1,3)
                    g = gridPoints(u)
                    if (g > nrGridActive) then
                        print *, 'Laplace: find_electrons: g > nrGridActive'
                        print *, 'Laplace: find_electrons: pos =',elec_pos
                    end if
                    if (is_emitter(u) == 0) then
                        if (g /= 0) then 
                            ! new_pos = cart_coord(u)
                            ! print *, 'Laplace: adding boundary: x=',new_pos(1), ' y=',new_pos(2), ' z=',new_pos(3)
                            ! $OMP CRITICAL
                            boundary_voltage = charge_voltage(elec_pos,u) 
                            ! print *, 'Laplace: adding boundary voltage:', boundary_voltage
                            iCharge(g) = iCharge(g) + boundary_voltage
                            ! $OMP END CRITICAL
                        end if
                    end if

                    u = move(move(move(i,1,2),1,1),1,3)
                    g = gridPoints(u)
                    if (g > nrGridActive) then
                        print *, 'Laplace: find_electrons: g > nrGridActive'
                        print *, 'Laplace: find_electrons: pos =',elec_pos
                    end if
                    if (is_emitter(u) == 0) then
                        if (g /= 0) then 
                            ! new_pos = cart_coord(u)
                            ! print *, 'Laplace: adding boundary: x=',new_pos(1), ' y=',new_pos(2), ' z=',new_pos(3)
                            ! $OMP CRITICAL
                            boundary_voltage = charge_voltage(elec_pos,u) 
                            ! print *, 'Laplace: adding boundary voltage:', boundary_voltage
                            iCharge(g) = iCharge(g) + boundary_voltage
                            ! $OMP END CRITICAL
                        end if
                    end if

                end if
            end if
        end do
        ! $OMP END PARALLEL DO
    end subroutine find_electrons

! -----------------------------------------------------------------------------
! ----- Indexing functions ----------------------------------------------------
! -----------------------------------------------------------------------------

    function disc_coord(x,y,z) ! DONE
        ! Change cartesian coordinates to discrete coordinates
        double precision, intent(in) :: x, y, z
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
        double precision, dimension(3) :: cart_coord
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

    function is_top(i)
        integer(kind=8), intent(in) :: i
        integer(kind=8) :: is_top
        
        if (i >= (Nz-1)*NxNy) then
            is_top = 1
        else
            is_top = 0
        end if

    end function is_top

    function is_emitter(i)
        ! Returns 1 if point is within the emitter and 0 otherwise
        integer(kind=8), intent(in) :: i
        integer(kind=8) :: is_emitter
        double precision, dimension(3) :: point_coord
        double precision :: point_radius

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
        double precision, dimension(3), intent(in) :: pos
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
        double precision, dimension(3) :: par_pos, par_vel

        if (step == 1) then
            par_pos(1) = 0.0d0*length_scale
            par_pos(2) = 0.0d0*length_scale
            par_pos(3) = 1.0d0*length_scale 
            par_vel = 0.0d0
            call Add_Particle(par_pos,par_vel,species_elec,step,1,-1)
            ! par_pos(1) = 3.0d0*length_scale
            ! par_pos(2) = 0.0d0*length_scale
            ! par_pos(3) = 1.0d0*length_scale 
            ! par_vel = 0.0d0
            ! call Add_Particle(par_pos,par_vel,species_elec,step,1,-1)
        end if
    end subroutine Place_Electron

    subroutine Write_Laplace_Data()
        integer :: ud_lp_field, ud_ic_field, IFAIL, lvl, i, j, k
        character(len=128) :: filename_lp_field, filename_ic_field
        double precision, dimension(3) :: ic_pos, ic_field
        
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

            do i=0,Nx-1
                do j=0,Ny-1
                    k = i + j*Nx+1 + lvl*NxNy
                    ic_pos = cart_coord(k)
                    ic_field = Calc_Field_at(ic_pos)
                    ! print *, 'Laplace: writing: x=',test_pos(1), ' y=',test_pos(2), ' z=',test_pos(3)
                    ! write(unit=ud_lp_voltage,iostat=IFAIL) i, j, voltage(k), is_emitter(k)
                    write(unit=ud_lp_field,iostat=IFAIL) i, j, laplace_field(k,3), is_emitter(k)
                    write(unit=ud_ic_field,iostat=IFAIL) i, j, ic_field(3), is_emitter(k)
                end do
            end do

            ! close(unit=ud_lp_voltage, iostat=IFAIL, status='keep')
            close(unit=ud_lp_field, iostat=IFAIL, status='keep')
            close(unit=ud_ic_field, iostat=IFAIL, status='keep')
        end do
        
    end subroutine Write_Laplace_Data

    function Calculate_Laplace_Field_at(pos)
        double precision, dimension(3), intent(in) :: pos
        double precision, dimension(3) :: Calculate_Laplace_Field_at
        double precision, dimension(3) :: pos_bot1, pos_bot2, pos_bot3, pos_bot4, pos_bot5, pos_bot6
        double precision, dimension(3) :: pos_top1, pos_top2, pos_top3, pos_top4, pos_top5, pos_top6
        double precision, dimension(3) :: field_bot1, field_bot2, field_bot3, field_bot4, field_top1, field_top2, field_top3, field_top4
        double precision, dimension(3) :: field_bot5, field_bot6, field_bot, field_top5, field_top6, field_top
        integer(kind=8) :: i,k

        if (is_inside(pos) == 0) then
            Calculate_Laplace_Field_at = Calc_Field_at(pos)
            return
        end if

        i = disc_coord(pos(1), pos(2), pos(3))

        pos_bot1 = cart_coord(k)
        field_bot1 = laplace_field(k,:)
        k = move(i,1,1)
        pos_bot2 = cart_coord(k)
        field_bot2 = laplace_field(k,:)
        k = move(i,1,2)
        pos_bot3 = cart_coord(k)
        field_bot3 = laplace_field(k,:)
        k = move(move(i,1,2),1,1)
        pos_bot4 = cart_coord(k)
        field_bot4 = laplace_field(k,:)

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

    end function Calculate_Laplace_Field_at

    function linear_interpolate_3d(x1,x2,xi,f1,f2)
        double precision, intent(in) :: x1,x2,xi
        double precision, dimension(3), intent(in) :: f1, f2
        double precision, dimension(3) :: linear_interpolate_3d
        integer(kind=8) :: i

        do i=1,3
            linear_interpolate_3d(i) = f1(i) + (f2(i)-f1(i))/(x2-x1)*(xi-x1)
        end do

    end function linear_interpolate_3d

end module mod_laplace_solver