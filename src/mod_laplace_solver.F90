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
    use iso_c_binding, only: c_int, c_ptr, c_f_pointer

    implicit none

    ! Accessibility
    private
    public :: Init_Laplace_Solver, Calculate_Laplace_Field, Calculate_Laplace_Field_at, Clean_Up_Laplace, Place_Electron, Write_Laplace_Data

    ! Pardiso variables
    type(MKL_PARDISO_HANDLE), dimension(64) :: pt
    integer(kind=8), dimension(64) :: iparm
    integer(kind=8), allocatable, dimension(:) :: perm
    integer(kind=8) :: maxfct=1, mnum=1, mtype=11, phase, n, nrhs=1
    integer(kind=8) :: msglvl=0, error, analysis=11, factorization=22, solving=33, solving1=331, solving2=333

    double precision, allocatable, dimension(:) :: voltage, laplace_field, x

    ! Matrix arrays
    double precision, allocatable, dimension(:) :: values, iCharge, oCharge, b
    integer(kind=8), allocatable, dimension(:) :: col_index

    double precision, allocatable, dimension(:) :: nnz_values
    integer(kind=8), allocatable, dimension(:) :: nnz_col_index, nnz_ia

    ! System parameters
    double precision :: hx, hy, hz, emitter_radius
    double precision, dimension(3) :: div_h2, div_h3
    double precision, dimension(2) :: lim_x, lim_y, lim_z
    integer(kind=8) :: Nx, Ny, Nz, nrGrid, nrGridActive, NxNy, n7, n19, nnz

contains

! -----------------------------------------------------------------------------
! ----- Solver ----------------------------------------------------------------
! -----------------------------------------------------------------------------

    subroutine Calculate_Laplace_Field(step) ! TODO
        integer(kind=8), intent(in) :: step
        ! Calculate the electric field
        print *, 'Laplace: finding electrons'
        call find_electrons()
        print *, 'Laplace: updating matrix'
        call update_matrix()
        print *, 'Laplace: solving matrix'
        call solve_matrix()
        print *, 'Laplace: allocate voltage'
        call allocate_voltage()
        print *, 'Laplace: calculating field'
        call calculate_field()
        print *, 'Laplace: done'

        ! call write_average_field()
        
    end subroutine Calculate_Laplace_Field

    subroutine solve_matrix() ! TODO
        ! Solve the system
        integer :: i
        ! double precision, dimension(nrGridActive) :: ddum

        print *, '  Laplace: solve_matrix: start'
        ! Solver parameters
        ! iparm(2) = 3 ! parallel nested dissection algorithm
        ! iparm(24) = 1 ! 10 = parallel factorization
        ! iparm(25) = 0 ! 2 = parallel solve; must be coupled with iparm(60)=0
        ! iparm(31) = 1 ! 0 = use all rhs; 2 = use sparsity of rhs
        ! iparm(35) = 0 ! one-based indexing (required)
        ! iparm(37) = 0 ! CSR format (required)
        ! iparm(60) = 1 ! 0 = in-core-memory (fast); 1 = mix as necessary; 2 = out-of-core memory (slow)

        ! iparm(1) = 1 ! Use customized iparm
        ! iparm(2) = 3 ! Use parallel METIS algorithm
        ! iparm(10) = 13 ! Perturb pivots
        ! iparm(11) = 1 ! Enable scaling
        ! iparm(24) = 1 ! Two-level parallel factorization
        ! iparm(31) = 0 ! perm(i) = 1 if b(i) is nonzero
        iparm(35) = 0 ! one-based indexing
        iparm(37) = 0 ! CSR format
        iparm(60) = 1 ! 0 = in-core-memory (fast); 1 = mix as necessary; 2 = out-of-core memory (slow)

        ! iparm(31) = 1 ! perm(i) = 1 if b(i) is nonzero
        ! iparm(35) = 0 ! one-based indexing
        ! iparm(37) = 0 ! CSR format
        ! iparm(60) = 1 ! 0 = in-core-memory; 1 = mix as necessary; 2 = out-of-core memory


        print *, '  Laplace: allocating reduced CSR memory'
        deallocate(laplace_field)
        allocate(nnz_values(nnz), nnz_col_index(nnz), nnz_ia(nrGridActive+1), x(nrGridActive))
        print *, '  Laplace: allocating reduced CSR memory successful'
        print *, '  Laplace: reducing matrix'
        call reduce_matrix()
        print *, '  Laplace: reducing matrix successful'

        ! print *, '  Laplace: solve_matrix: nnz=',nnz
        ! print *, '  Laplace: solve_matrix: size(nnz_values)=',size(nnz_values)
        ! print *, '  Laplace: solve_matrix: size(nnz_col_index)=',size(nnz_col_index)
        ! print *, '  Laplace: solve_matrix: size(nnz_ia)=',size(nnz_ia)


        ! call check_matrix()


        ! Analysis ----------------------------
        print *, '  Laplace: solve_matrix: analysis and symbolic factorization'
        phase = analysis  
        iparm(27) = 1 ! Check matrix on for analysis
        call pardiso_64(pt,maxfct,mnum,mtype,phase,n,nnz_values,nnz_ia,nnz_col_index,perm,nrhs,iparm,msglvl,b,x,error)
        iparm(27) = 0 ! Check matrix off for factorization and solving
        if (error == 0) then
            print *, '  Laplace: solve_matrix: analysis and symbolic factorization successful'
        else
            print *, '  Laplace: solve_matrix: analysis and symbolic factorization with error code ', error
        end if
        print *, '  Laplace: solve_matrix:', iparm(17)*1e-6, 'Gb needed for numerical factorization'

        ! Factorization ----------------------------
        print *, '  Laplace: solve_matrix: numerical factorization'
        phase = factorization     
        call pardiso_64(pt,maxfct,mnum,mtype,phase,n,nnz_values,nnz_ia,nnz_col_index,perm,nrhs,iparm,msglvl,b,x,error)
        if (error == 0) then
            print *, '  Laplace: solve_matrix: numerical factorization successful'
        else
            print *, '  Laplace: solve_matrix: numerical factorization failed with error code ', error
        end if

        ! Solving ----------------------------
        if (solving /= 0) then
            print *, '  Laplace: solve_matrix: solving'
            phase = solving    
            call pardiso_64(pt,maxfct,mnum,mtype,phase,nrGridActive,nnz_values,nnz_ia,nnz_col_index,perm,nrhs,iparm,msglvl,b,x,error)
            if (error == 0) then
                print *, '  Laplace: solve_matrix: solving successful'
            else
                print *, '  Laplace: solve_matrix: solving failed with error code ', error
            end if
        else
            print *, '  Laplace: solve_matrix: forward substitution'
            phase = solving1    
            call pardiso_64(pt,maxfct,mnum,mtype,phase,nrGridActive,nnz_values,nnz_ia,nnz_col_index,perm,nrhs,iparm,msglvl,b,x,error)
            if (error == 0) then
                print *, '  Laplace: solve_matrix: forward substitution successful'
            else
                print *, '  Laplace: solve_matrix: forward substitution failed with error code ', error
            end if

            print *, '  Laplace: solve_matrix: backward substitution'
            phase = solving2
            call pardiso_64(pt,maxfct,mnum,mtype,phase,nrGridActive,nnz_values,nnz_ia,nnz_col_index,perm,nrhs,iparm,msglvl,b,x,error)
            if (error == 0) then
                print *, '  Laplace: solve_matrix: backward substitution successful'
            else
                print *, '  Laplace: solve_matrix: backward substitution failed with error code ', error
            end if
        end if

        if (norm2(x)==0.0d0) then
            print *, 'Laplace: solve_matrix: x is zero'
        end if
        if (norm2(b)==0.0d0) then
            print *, 'Laplace: solve_matrix: b is zero'
        end if

        ! Clear memory ----------------------------
        print *, '  Laplace: solve_matrix: clearing memory'
        phase = -1
        call pardiso_64(pt,maxfct,mnum,mtype,phase,n,nnz_values,nnz_ia,nnz_col_index,perm,nrhs,iparm,msglvl,b,x,error)
        if (error == 0) then
            print *, '  Laplace: solve_matrix: clearing memory successful'
        else
            print *, '  Laplace: solve_matrix: clearing memory failed with error code ', error
        end if

        ! Deallocate the reduced CSR matrix
        deallocate(nnz_values)
        deallocate(nnz_ia)
        deallocate(nnz_col_index) 
    end subroutine solve_matrix

! -------------------------------------------------------------------------------------------------------------
! ----- Voltage and field calculations ------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------------------

    subroutine allocate_voltage() ! DONE & PARALLEL
        integer(kind=8) :: i, g

        voltage = 0.0d0

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP& SHARED(nrGrid, voltage, x, gridPoints) &
        !$OMP& PRIVATE(i, g)

        do i=1,nrGrid
            if (gridPoints(i) /= 0) then
                g = gridPoints(i)
                voltage(i) = x(g)
            end if
        end do

        !$OMP END PARALLEL DO

        ! Deallocate the solution vector
        deallocate(x)
    end subroutine allocate_voltage

    subroutine calculate_field() ! DONE & PARALLEL
        integer(kind=8) :: i, boundary
        integer(kind=8), dimension(3) :: boundaries

        allocate(laplace_field(nrGrid))

        laplace_field = 0.0d0

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP& SHARED(nrGrid, voltage, laplace_field) &
        !$OMP& PRIVATE(i, boundaries, boundary)

        do i=1,nrGrid
            boundaries = is_boundary(i)
            boundary = boundaries(3)

            select case (boundary)
                case (0) ! Inner point
                    laplace_field(i) = - central_difference(voltage,i) + applied_field(i)
                case (1) ! Startpoint
                    laplace_field(i) = - forward_difference(voltage,i) + applied_field(i)
                case (2) ! Endpoint
                    laplace_field (i) = - backward_difference(voltage,i) + applied_field(i)
            end select
        end do

        !$OMP END PARALLEL DO

        ! Deallocate the voltage vector
        deallocate(voltage)

    end subroutine calculate_field

    function average_field() ! DONE & PARALLEL
        integer(kind=8) :: i, j, k, num
        double precision :: sum, average_field

        sum = 0.0d0
        num = 0

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP& SHARED(Nx, Ny, laplace_field) &
        !$OMP& PRIVATE(i, j, k) &
        !$OMP& REDUCTION(+:sum, num)

        do i=0,Nx-1
            do j=0,Ny-1
                k = i + j*Nx+1
                sum = sum + laplace_field(k)
                num = num + 1
            end do
        end do

        !$OMP END PARALLEL DO

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

! -----------------------------------------------------------------------------
! ----- Initialization --------------------------------------------------------
! -----------------------------------------------------------------------------

    subroutine Init_Laplace_Solver(nthreads) ! DONE
        ! Initialize environment based on input file
        double precision :: Lx, Ly, Lz
        integer(kind=8), intent(in) :: nthreads

        ! print *, 'Laplace: initializing solver'

        Nx = laplace_intervals(1); Ny = laplace_intervals(2); Nz = laplace_intervals(3)
        Lx = laplace_dim(1); Ly = laplace_dim(2); Lz = laplace_dim(3)
        ! print *, '  Nx=',Nx,'   Ny=',Ny,'   Nz=',Nz

        NxNy = Nx*Ny
        nrGrid = Nx * Ny * Nz
        nrGridActive = 0


        hx = Lx / (Nx-1)
        hy = Ly / (Ny-1)
        hz = Lz / (Nz-1)
        div_h2 = (/1.0d0 / (hx**2), 1.0d0 / (hy**2), 1.0d0 / (hz**2)/)
        div_h3 = (/1.0d0 / (hx**3), 1.0d0 / (hy**3), 1.0d0 / (hz**3)/)

        lim_x = (/laplace_pos(1), laplace_pos(1)+Lx/)
        lim_y = (/laplace_pos(2), laplace_pos(2)+Ly/)
        lim_z = (/laplace_pos(3), laplace_pos(3)+Lz/)


        emitter_radius = emitters_dim(1,1)

        allocate(gridPoints(nrGrid), gridPointsActive(nrGrid))
        print *, 'Laplace: initializing gridpoints'
        call init_grid()

        n = nrGridActive
        n7 = 7*nrGridActive
        n19 = 19*nrGridActive

        print *, 'Laplace: allocating arrays'
        allocate(values(n19), col_index(n19), b(nrGridActive), oCharge(nrGridActive))
        allocate(voltage(nrGrid), laplace_field(nrGrid))
        allocate(perm(nrGridActive))

        print *, 'Laplace: initializing matrix'
        call init_matrix()

        print *, 'Laplace: initializing Pardiso'
        call init_pardiso(nthreads)

        call write_grid()
        
    end subroutine Init_Laplace_Solver

    subroutine init_pardiso(nthreads) ! DONE
        integer(kind=8), intent(in) :: nthreads
        ! call mkl_set_num_threads(nthreads)
        call pardisoinit(pt,mtype,iparm)
    end subroutine init_pardiso

    subroutine init_grid() ! DONE
        ! Store points that are outside of the emitter or one layer within it
        integer(kind=8) :: i


        do i=1,nrGrid
            ! print *, 'Laplace: gridpoint nr', i
            ! print *, '      coord=', cart_coord(i)
            if (is_emitter(i) == 1) then
                ! print *, '  it is emitter'
                if (is_first_layer(i) == 1) then
                    ! print *, '      it is first_layer'
                    nrGridActive = nrGridActive + 1
                    gridPointsActive(nrGridActive) = i
                    gridPoints(i) = nrGridActive
                else
                    ! print *, 'Laplace: checking not first layer'
                    gridPoints(i) = 0
                end if
            else
                ! print *, '  it is not emitter'
                ! print *, 'Laplace: checking not emitter'
                ! print *, 'Laplace: updating nrGridActive'
                nrGridActive = nrGridActive + 1
                ! print *, 'Laplace: updating gridPointsActive at i=',nrGridActive
                gridPointsActive(nrGridActive) = i
                ! print *, 'Laplace: updating gridPoints at i=', i
                gridPoints(i) = nrGridActive
            end if
        end do

    end subroutine init_grid

    subroutine init_matrix() ! DONE
        ! Initialize matrix without boundary conditions
        integer(kind=8) :: g, i, a, center
        integer(kind=8), dimension(3) :: boundaries
        integer(kind=8) :: next, next_next, next_next_next, prev, prev_prev, prev_prev_prev

        values = values*0.0d0
        b = b*0.0d0
        ocharge = 0.0d0
        icharge = 0.0d0

        center = 10

        do g=1,nrGridActive

            col_index(center) = g
            do a=1,3
                next = center + (a-1)*3 + 1
                next_next = next + 1
                next_next_next = next + 2

                prev = center - (a-1)*3 - 1
                prev_prev = prev - 1
                prev_prev_prev = prev - 2

                col_index(prev_prev_prev) = move(g,-3,a)
                col_index(prev_prev) = move(g,-2,a)
                col_index(prev) = move(g,-1,a)
                col_index(next) = move(g,1,a)
                col_index(next_next) = move(g,2,a)
                col_index(next_next_next) = move(g,3,a)
            end do

            ! rows_start(g) = g
            ! rows_end(g) = g

            i = gridPointsActive(g)
            
            if (is_emitter(i) == 1) then
                call insert_emitter_boundary(g)
            else
                call insert_finite_difference(g)
            end if

            center = center+19

        end do
    end subroutine init_matrix

! -----------------------------------------------------------------------------
! ----- Matrix operations -----------------------------------------------------
! -----------------------------------------------------------------------------

    subroutine update_matrix() ! DONE & PARALLEL
        ! Update boundary conditions in matrix
        integer(kind=8) :: g

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP& SHARED(nrGridActive, iCharge, oCharge) &
        !$OMP& PRIVATE(g)

        do g=1,nrGridActive
            if (oCharge(g) /= 0.0d0) then

                call remove_electron_boundary(g)

                if (iCharge(g) == 0.0d0) then
                    call insert_finite_difference(g)
                end if

                oCharge(g) = 0.0d0
            else
                if (iCharge(g) /= 0.0d0) then
                    call remove_finite_difference(g)
                end if
            end if

            if (iCharge(g) /= 0.0d0) then
                call insert_electron_boundary(g)
                oCharge(g) = iCharge(g)
            end if
        end do

        !$OMP END PARALLEL DO

        ! iCharge has been written in oCharge and is no longer needed
        deallocate(iCharge)
    end subroutine update_matrix

    subroutine reduce_matrix()
        integer(kind=8) :: g, j, index
        integer(kind=8) :: cur_ind, cur_col
        logical :: start

        ! nnz_rows_start = [rows_start, nnz]
        ! nnz_rows_end = [rows_end, nnz]

        index = 0
        cur_ind = 0

        ! $OMP PARALLEL DO DEFAULT(NONE) &
        ! $OMP& SHARED(nrGridActive, values, col_index, nnz_values, nnz_col_index, nnz_ia, nnz_rows_start, nnz_rows_end) &
        ! $OMP& PRIVATE(g, j, cur_col, start) &
        ! $OMP& REDUCTION(+:index, cur_ind)

        do g=1,nrGridActive
            start = .true.

            ! if (g>100 .and. g<=200) then
            !     print *, 'Laplace: reduce_matrix: g=',g
            ! end if
            do j=1,19
                cur_ind = cur_ind + 1
                cur_col = col_index(cur_ind)
                ! if (g>100 .and. g<=200) then
                !     print *, '  	values(cur_ind)=',values(cur_ind), ' cur_col=',cur_col, ' index=',index+1
                ! end if
                if ((values(cur_ind) < 0.0d0) .or. (values(cur_ind) > 0.0d0)) then
                    index = index+1
                    nnz_values(index) = values(cur_ind)
                    nnz_col_index(index) = cur_col

                    if (start .eqv. .true.) then
                        nnz_ia(g) = index
                        start = .false.
                    end if
                    
                    ! if (cur_col < nnz_rows_start(g)) then
                    !     nnz_rows_start(g) = cur_col
                    ! end if
                    ! if (cur_col > nnz_rows_end(g)) then
                    !     nnz_rows_end(g) = cur_col
                    ! end if
                end if
            end do
            ! if (g>100 .and. g<=200) then
            !     print *, '   nnz_ia(g)=', nnz_ia(g)
            ! end if
        end do

        ! $OMP END PARALLEL DO

        nnz_ia(nrGridActive+1) = nnz


    end subroutine reduce_matrix

    subroutine insert_electron_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = (g-1)*19 + 10

        values(center) = 1.0d0
        b(g) = iCharge(g)
        perm(g) = 1

        nnz = nnz + 1
    end subroutine insert_electron_boundary

    subroutine remove_electron_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = (g-1)*19 + 10

        if (values(center) /= 0.0d0) then
            values(center) = 0.0d0
            nnz = nnz - 1
        end if

        b(g) = 0.0d0
        perm(g) = 0
    end subroutine remove_electron_boundary

    subroutine insert_emitter_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = (g-1)*19 + 10

        values(center) = 1.0d0
        b(g) = 0.0d0
        perm(g) = 0

        nnz = nnz + 1
    end subroutine insert_emitter_boundary

    subroutine insert_finite_difference(g) ! DONE
        ! Insert finite difference equation into matrix
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center
        integer(kind=8), dimension(3) :: boundaries
        integer(kind=8) :: i, a, next, next_next, next_next_next, prev, prev_prev, prev_prev_prev

        center = (g-1)*19 + 10

        i = gridPointsActive(g)

        boundaries = is_boundary(i)

        ! Insert finite difference equation for each axis
        do a=1,3

            next = center + (a-1)*3 + 1
            next_next = next + 1
            next_next_next = next + 2
            prev = center - (a-1)*3 - 1
            prev_prev = prev - 1
            prev_prev_prev = prev - 2

            select case (boundaries(a))
                case (0) ! Inner point (values for central difference)

                    values(prev) = values(prev) + div_h2(a)

                    values(center) = values(center) - 2.0d0*div_h2(a)

                    values(next) = values(next) + div_h2(a)

                    nnz = nnz + 2 ! prev and next

                case (1) ! Startpoint (values for forward difference)

                    values(center) = values(center) + 2.0d0*div_h2(a)

                    values(next) = values(next) - 5.0d0*div_h2(a)

                    values(next_next) = values(next_next) + 4.0d0*div_h2(a)

                    values(next_next_next) = values(next_next_next) - 1.0d0*div_h2(a)

                    nnz = nnz + 3 ! three next

                case (2) ! Endpoint (values for backward difference)

                    values(prev_prev_prev) = values(prev_prev_prev) - 1.0d0*div_h2(a)

                    values(prev_prev) = values(prev_prev) + 4.0d0*div_h2(a)

                    values(prev) = values(prev) - 5.0d0*div_h2(a)

                    values(center) = values(center) + 2.0d0*div_h2(a)

                    nnz = nnz + 3 ! three prev

            end select
        end do

        b(g) = 0.0d0
        perm(g) = 0

        if (values(center) /= 0.0d0) then
            nnz = nnz + 1 ! center
        end if
    end subroutine insert_finite_difference

    subroutine remove_finite_difference(g) ! DONE
        ! Remove finite difference equation from matrix
        integer(kind=8), intent(in) :: g
        integer(kind=8), dimension(3) :: boundaries
        integer(kind=8) :: i, center

        center = (g-1)*19 + 10

        if (values(center) /= 0.0d0) then
            values(center) = 0.0d0
            nnz = nnz - 1
        end if

        do i=1,8
            if (values(center+i) /= 0.0d0) then
                values(center+i) = 0.0d0
                nnz = nnz - 1
            end if
            if (values(center-i) /= 0.0d0) then
                values(center-i) = 0.0d0
                nnz = nnz - 1
            end if
        end do

    end subroutine remove_finite_difference

! -----------------------------------------------------------------------------
! ----- Finite difference -----------------------------------------------------
! -----------------------------------------------------------------------------

    function central_difference(array,index)
        double precision, dimension(nrGrid), intent(in) :: array
        integer(kind=8), intent(in) :: index
        integer(kind=8) :: index_prev, index_next
        double precision :: central_difference

        index_prev = move(index,-1,3)
        index_next = move(index,1,3)

        central_difference = (array(index_next) - array(index_prev))/(2.0d0*hz)
    end function central_difference

    function forward_difference(f,x)
        double precision, dimension(nrGrid), intent(in) :: f
        integer(kind=8), intent(in) :: x
        integer(kind=8) :: x_1h, x_2h, x_3h, x_4h
        double precision :: forward_difference

        x_1h = move(x,1,3)
        x_2h = move(x,2,3)
        x_3h = move(x,3,3)
        x_4h = move(x,4,3)

        forward_difference = (-3.0d0*f(x)+4.0d0*f(x_1h)-f(x_2h))/(2.0d0*hz)
        ! forward_difference = (-25.0d0*f(x)+48.0d0*f(x_1h)-36*f(x_2h)+16.0d0*f(x_3h)-3.0d0*f(x_4h))/(12.0d0*hz)
    end function forward_difference

    function backward_difference(array,index)
        double precision, dimension(nrGrid), intent(in) :: array
        integer(kind=8), intent(in) :: index
        integer(kind=8) :: index_prev, index_prev_prev
        double precision :: backward_difference

        index_prev = move(index,-1,3)
        index_prev_prev = move(index,-2,3)

        backward_difference = (array(index_prev_prev)-4.0d0*array(index_prev)+3.0d0*array(index))/(2.0d0*hz)
    end function backward_difference

! -----------------------------------------------------------------------------
! ----- Boundary conditions ---------------------------------------------------
! -----------------------------------------------------------------------------

    function applied_field(i)
        integer(kind=8), intent(in) :: i
        double precision, dimension(3) :: grid_coord
        double precision :: applied_field

        grid_coord = cart_coord(i)

        applied_field = -V_d / abs(emitters_pos(3,1)+box_dim(3)-grid_coord(3))
    end function applied_field

    function elec_voltage(elec_pos, i) ! DONE
        ! Calculate the boundary value for an electron
        double precision, dimension(3), intent(in) :: elec_pos
        integer(kind=8), intent(in) :: i
        double precision, dimension(3) :: grid_coord
        double precision :: elec_voltage
        
        grid_coord = cart_coord(i)
        elec_voltage = -q_0 / ((4.0d0*pi*epsilon_0) * norm2(grid_coord - elec_pos))
    end function elec_voltage

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
                if ((lim_x(1)<=elec_pos(1)) .and. (elec_pos(1)<=lim_x(2)) &
                    .and. (lim_y(1)<=elec_pos(2)) .and. (elec_pos(2)<=lim_y(2)) &
                    .and. (lim_z(1)<=elec_pos(3)) .and. (elec_pos(3)<=lim_z(2))) then

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
                            boundary_voltage = elec_voltage(elec_pos,u) 
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
                            boundary_voltage = elec_voltage(elec_pos,u) 
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
                            boundary_voltage = elec_voltage(elec_pos,u) 
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
                            boundary_voltage = elec_voltage(elec_pos,u) 
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
                            boundary_voltage = elec_voltage(elec_pos,u) 
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
                            boundary_voltage = elec_voltage(elec_pos,u) 
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
                            boundary_voltage = elec_voltage(elec_pos,u) 
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
                            boundary_voltage = elec_voltage(elec_pos,u) 
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
        ! Clean up memory
        print '(a)', 'Laplace: Clean up'
        ! Deallocate CSR matrix
        deallocate(values)
        deallocate(col_index)
        deallocate(b)
        ! Deallocate other arrays
        deallocate(oCharge)
        deallocate(gridPoints)
        deallocate(gridPointsActive)
        deallocate(perm)

        ! These should have been deallocated already in the last iteration
        if (allocated(iCharge)) then 
            deallocate(iCharge)
        end if
        if (allocated(laplace_field)) then
            deallocate(laplace_field) 
        end if
        if (allocated(nnz_values)) then
            deallocate(nnz_values)
        end if
        if (allocated(nnz_col_index)) then
            deallocate(nnz_col_index)
        end if
        if (allocated(nnz_ia)) then
            deallocate(nnz_ia)
        end if
        if (allocated(voltage)) then
            deallocate(voltage)
        end if
        if (allocated(x)) then
            deallocate(x)
        end if
    end subroutine Clean_Up_Laplace

    ! subroutine check_matrix()
    !     integer(kind=8) :: i,count,low,high,col,index


    !     index = 1
    !     do i=1,nrGridActive
    !         count = 0
    !         low = nnz_rows_start(i)
    !         high = nnz_rows_end(i)
    !         col = nnz_col_index(index)

    !         do while ((low <= col) .and. (col < high))
    !             index = index + 1
    !             col = nnz_col_index(index)
    !             count = count + 1
    !         end do
    !         count  = count + 1
    !         index = index + 1

    !         print *, 'Laplace: row=',i,' count=',count

    !     end do


    ! end subroutine check_matrix

    subroutine Place_Electron(step)
        integer(kind=8), intent(in) :: step
        double precision, dimension(3) :: par_pos, par_vel

        if (step == 1) then
            par_pos(1) = 0.0d0*length_scale
            par_pos(2) = 0.0d0*length_scale
            par_pos(3) = 0.5d0*length_scale 
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
                    write(unit=ud_lp_field,iostat=IFAIL) i, j, laplace_field(k), is_emitter(k)
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
        double precision, dimension(3) :: pos_bot1, pos_bot2, pos_bot3, pos_bot4, pos_top1, pos_top2, pos_top3, pos_top4
        double precision, dimension(3) :: field_bot1, field_bot2, field_bot3, field_bot4, field_top1, field_top2, field_top3, field_top4
        double precision, dimension(3) :: field_bot5, field_bot6, field_bot, field_top5, field_top6, field_top
        integer(kind=8) :: i,k

        if (is_inside(pos) == 0) then
            Calculate_Laplace_Field_at = (/0.0d0,0.0d0,0.0d0/)
            return
        end if

        ! i = disc_coord(pos(1), pos(2), pos(3))

        ! pos_bot1 = cart_coord(k)
        ! field_bot1 = laplace_field(k)
        ! k = move(i,1,1)
        ! pos_bot2 = cart_coord(k)
        ! field_bot2 = laplace_field(k)
        ! k = move(i,1,2)
        ! pos_bot3 = cart_coord(k)
        ! field_bot3 = laplace_field(k)
        ! i = move(move(i,1,2),1,1)
        ! pos_bot4 = cart_coord(k)
        ! field_bot4 = laplace_field(k)

        ! k = move(i,1,3)
        ! pos_top1 = cart_coord(k)
        ! field_top1 = laplace_field(k)
        ! k = move(move(i,1,3),1,1)
        ! pos_top2 = cart_coord(k)
        ! field_top2 = laplace_field(k)
        ! k = move(move(i,1,3),1,2)
        ! pos_top3 = cart_coord(k)
        ! field_top3 = laplace_field(k)
        ! k = move(move(move(i,1,3),1,2),1,1)
        ! pos_top4 = cart_coord(k)
        ! field_top4 = laplace_field(k)

        ! ! Bottom interpolation

        ! field_bot5 = linear_interpolate_3d(pos_bot1(1),pos_bot2(1),pos(1),field_bot1,field_bot2)
        ! field_bot6 = linear_interpolate_3d(pos_bot3(1),pos_bot4(1),pos(1),field_bot3,field_bot4)
        ! field_bot = linear_interpolate_3d(pos_bot5(2),pos_bot6(2),pos(2),field_bot5,field_bot6)

        ! ! Top interpolation
        ! field_top5 = linear_interpolate_3d(pos_top1(1),pos_top2(1),pos(1),field_top1,field_top2)
        ! field_top6 = linear_interpolate_3d(pos_top3(1),pos_top4(1),pos(1),field_top3,field_top4)
        ! field_top = linear_interpolate_3d(pos_top5(2),pos_top6(2),pos(2),field_top5,field_top6)

        ! ! Final interpolation
        ! Calculate_Laplace_Field_at = linear_interpolate_3d(pos_bot(3),pos_top(3),pos(3),field_bot,field_top)

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