! ------------------------------------------------------------------------------
! ----- Module: mod_laplace_solver ---------------------------------------------
! ----- Author: Arnar Jonsson --------------------------------------------------
! ----- Date: 2025 -------------------------------------------------------------
! ------------------------------------------------------------------------------

module mod_laplace_solver
    
    ! Dependencies
    use mod_global
    use mod_pair
    use mkl_pardiso ! Required for pardiso solver
    use mkl_spblas  ! Required for sparse matrices
    use iso_c_binding, only: c_int, c_ptr, c_f_pointer

    implicit none

    ! Accessibility
    private
    public :: Init_Laplace_Solver, Calculate_Laplace_Field, Clean_Up_Laplace, Place_Electron, Write_Laplace_Data

    ! Local variables
    type(MKL_PARDISO_HANDLE), dimension(64) :: pt
    ! integer*8, dimension(64) :: pt
    integer, dimension(64) :: iparm

    double precision, allocatable, dimension(:) :: voltage, laplace_field

    ! Matrix arrays
    double precision, allocatable, dimension(:) :: values, iCharge, oCharge, b
    integer(kind=8), allocatable, dimension(:) :: rows_start, rows_end, col_index

    double precision, allocatable, dimension(:) :: nnz_values
    integer(kind=8), allocatable, dimension(:) :: nnz_rows_start, nnz_rows_end, nnz_col_index, nnz_ia

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
        double precision, dimension(nrGridActive) :: x
        integer(kind=8), intent(in) :: step
        ! Calculate the electric field
        print *, 'Laplace: finding electrons'
        call find_electrons()
        print *, 'Laplace: updating matrix'
        call update_matrix()
        print *, 'Laplace: solving matrix'
        call solve_matrix(x)
        print *, 'Laplace: allocate voltage'
        call allocate_voltage(x)
        print *, 'Laplace: calculating field'
        call calculate_field()
        print *, 'Laplace: done'

        call write_average_field()
        
    end subroutine Calculate_Laplace_Field

    subroutine solve_matrix(x) ! TODO
        ! Solve the system
        integer :: i
        double precision, dimension(nrGridActive), intent(inout) :: x
        double precision, dimension(nrGridActive) :: ddum
        integer(kind=8), dimension(nrGridActive) :: perm
        integer(kind=8) :: maxfct,mnum,mtype,phase, n, nrhs
        integer(kind=8) :: msglvl, error, reordering, factorization, solving

        ! Matrix parameters
        n = nrGridActive
        nrhs = 1
        maxfct = 1
        mnum = 1
        mtype = 11
        ! Solver parameters
        msglvl = 0
        reordering = 11
        factorization = 12
        solving = 33

        iparm(27) = 1


        allocate(nnz_values(nnz), nnz_col_index(nnz), nnz_ia(nrGridActive+1), nnz_rows_start(nrGridActive+1), nnz_rows_end(nrGridActive+1))
        call reduce_matrix()

        ! call check_matrix()


        print *, '  Laplace: solve_matrix: reordering and symbolic factorization'
        phase = reordering  
        call pardiso_d(pt,maxfct,mnum,mtype,phase,n,nnz_values,nnz_ia,nnz_col_index,perm,nrhs,iparm,msglvl,ddum,ddum,error)
        if (error == 0) then
            print *, '  Laplace: solve_matrix: symbolic factorization successful'
        else
            print *, '  Laplace: solve_matrix: symbolic factorization failed with error code ', error
        end if

        iparm(27) = 0

        print *, '  Laplace: solve_matrix: factorization'
        phase = factorization     
        call pardiso_d(pt,maxfct,mnum,mtype,phase,n,nnz_values,nnz_ia,nnz_col_index,perm,nrhs,iparm,msglvl,ddum,ddum,error)
        if (error == 0) then
            print *, '  Laplace: solve_matrix: factorization successful'
        else
            print *, '  Laplace: solve_matrix: factorization failed with error code ', error
        end if

        print *, '  Laplace: solve_matrix: solving'
        phase = solving    
        call pardiso_d(pt,maxfct,mnum,mtype,phase,nrGridActive,nnz_values,nnz_ia,nnz_col_index,perm,nrhs,iparm,msglvl,b,x,error)
        

        select case (error)
            case (0)
                print *, '  Laplace: solve_matrix: solving successful'
            case (-1)
                print *, 'Laplace: solve_matrix: input inconsistent'
            case (-2)
                print *, 'Laplace: solve_matrix: not enough memory'
            case (-3)
                print *, 'Laplace: solve_matrix: reordering problem'
            case (-4)
                print *, 'Laplace: solve_matrix: zero pivot'
            case (-5)
                print *, 'Laplace: solve_matrix: unclassified error'
            case (-6)
                print *, 'Laplace: solve_matrix: reordering failed'
            case (-7)
                print *, 'Laplace: solve_matrix: diagonal matrix is singular'
            case (-8)
                print *, 'Laplace: solve_matrix: 32-bit overflow'
            case (-9)
                print *, 'Laplace: solve_matrix: not enough memory for OOC'
            case (-10)
                print *, 'Laplace: solve_matrix: error opening OOC files'
            case (-11)
                print *, 'Laplace: solve_matrix: read/write error with OOC files'
            case (-12)
                print *, 'Laplace: solve_matrix: pardiso_64 called from 32-bit library'
            case (-13)
                print *, 'Laplace: solve_matrix: interrupted by the mkl_progress function'
            case (-15)
                print *, 'Laplace: solve_matrix: internal error'
        end select

        if (norm2(x)==0.0d0) then
            print *, 'Laplace: solve_matrix: x is zero'
        end if
        if (norm2(b)==0.0d0) then
            print *, 'Laplace: solve_matrix: b is zero'
        end if

        deallocate(nnz_values, nnz_ia, nnz_col_index, nnz_rows_start, nnz_rows_end) 
    end subroutine solve_matrix

    ! subroutine solve_matrix_v2(x)
    !     integer :: i, index
    !     double precision, dimension(nrGridActive), intent(out) :: x
    !     integer, dimension(nnz) :: nnz_rowA, nnz_colA
    !     double precision, dimension(nnz) :: nnz_valA
    !     integer, dimension(nrGridActive) :: perm
    !     integer :: maxfct=1,mnum=1,mtype=11,phase
    !     integer :: msglvl, error
    !     x = 0.0d0

    !     i = 1
    !     index = 1
    !     ! $OMP PARALLEL DO DEFAULT(NONE) &
    !     ! $OMP& SHARED(nnz, nnz_rowA, nnz_colA, nnz_valA, rowA, colA, valA, linkA) &
    !     ! $OMP& PRIVATE(i,index)
    !     do i=1,nnz
    !         if (index > 7*nrGridActive) then
    !             print *, 'Laplace: index out of bounds'
    !             print *, 'Laplace: index=',index
    !             exit
    !         end if
    !         nnz_rowA(i) = rowA(index)
    !         nnz_colA(i) = colA(index)
    !         nnz_valA(i) = valA(index)
    !         index = linkA(index)
    !     end do
    !     phase = 23

    !     call pardiso_d(pt,maxfct,mnum,mtype,phase,nrGridActive,nnz_valA, nnz_rowA, nnz_colA, perm,nrGridActive,iparm,msglvl,b,x,error)
        
        
    ! end subroutine solve_matrix_v2

    subroutine allocate_voltage(x)
        double precision, dimension(nrGridActive), intent(in) :: x
        integer(kind=8) :: i, g

        voltage = 0.0d0

        do i=1,nrGrid
            if (gridPoints(i) /= 0) then
                g = gridPoints(i)
                voltage(i) = x(g)
            end if
        end do
    end subroutine allocate_voltage

    subroutine calculate_field()
        integer(kind=8) :: i, g, boundary
        integer(kind=8), dimension(3) :: boundaries

        laplace_field = 0.0d0

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

    end subroutine calculate_field

    function average_field()
        integer(kind=8) :: i, j, k, num
        double precision :: sum, average_field

        sum = 0.0d0
        num = 0

        do i=0,Nx-1
            do j=0,Ny-1
                k = i + j*Nx+1
                sum = sum + laplace_field(k)
                num = num + 1
            end do
        end do

        if (num > 0) then
            average_field = sum / num
        else
            average_field = 0.0d0
        end if
    end function average_field

    subroutine write_average_field()
        integer(kind=8) :: IFAIL
        double precision :: avg_field

        avg_field = average_field()

        write (ud_laplace_average_field, "(E16.8)", iostat=IFAIL) avg_field
    end subroutine write_average_field

    subroutine write_grid()
        integer(kind=8) :: IFAIL, i

        do i=1,NxNy
            write (ud_grid, "(i8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, i8)", iostat=IFAIL) is_emitter(i), cart_coord(i)
        end do
    end subroutine write_grid

! -----------------------------------------------------------------------------
! ----- Initialization --------------------------------------------------------
! -----------------------------------------------------------------------------

    subroutine Init_Laplace_Solver() ! DONE
        ! Initialize environment based on input file
        double precision :: Lx, Ly, Lz

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

        print *, 'Laplace: limits: x=',lim_x(1),lim_x(2), ' y=',lim_y(1),lim_y(2), ' z=',lim_z(1),lim_z(2)

        emitter_radius = emitters_dim(1,1)

        allocate(gridPoints(nrGrid), gridPointsActive(nrGrid))
        print *, 'Laplace: initializing gridpoints'
        call init_grid()

        n7 = 7*nrGridActive
        n19 = 19*nrGridActive

        print *, 'Laplace: allocating arrays'
        allocate(values(n19), col_index(n19), rows_start(nrGridActive), rows_end(nrGridActive), b(nrGridActive), iCharge(nrGridActive), oCharge(nrGridActive))

        allocate(voltage(nrGrid), laplace_field(nrGrid))

        print *, 'Laplace: initializing matrix'
        call init_matrix()

        print *, 'Laplace: initializing Pardiso'
        call init_pardiso()

        call write_grid()
        
    end subroutine Init_Laplace_Solver

    subroutine init_pardiso()
        integer :: mtype=11
        call pardisoinit(pt,mtype,iparm)
    end subroutine init_pardiso

    subroutine init_grid()
        ! Store points that are outside of the emitter or one layer within it
        integer(kind=8) :: i
        do i=1,nrGrid
            ! print *, 'Laplace: gridpoint'
            ! print *, '      coord=', cart_coord(i)
            if (is_emitter(i) == 1) then
                ! print *, '      it is emitter'
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

            rows_start(g) = g
            rows_end(g) = g

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

    subroutine update_matrix() ! DONE
        ! Update boundary conditions in matrix
        integer(kind=8) :: g

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
    end subroutine update_matrix

    subroutine reduce_matrix()
        integer(kind=8) :: g, j, index
        integer(kind=8) :: cur_ind, cur_col
        logical :: start

        nnz_rows_start = [rows_start, nnz]
        nnz_rows_end = [rows_end, nnz]

        index = 0
        cur_ind = 0
     
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
                    
                    if (cur_col < nnz_rows_start(g)) then
                        nnz_rows_start(g) = cur_col
                    end if
                    if (cur_col > nnz_rows_end(g)) then
                        nnz_rows_end(g) = cur_col
                    end if
                end if
            end do
            ! if (g>100 .and. g<=200) then
            !     print *, '   nnz_ia(g)=', nnz_ia(g)
            ! end if
        end do

        nnz_ia(nrGridActive+1) = nnz

    end subroutine reduce_matrix

    subroutine insert_electron_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = (g-1)*19 + 10

        values(center) = 1.0d0
        b(g) = iCharge(g)

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
    end subroutine remove_electron_boundary

    subroutine insert_emitter_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = (g-1)*19 + 10

        values(center) = 1.0d0
        b(g) = 0.0d0

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

                    values(center) = values(center) + 2.0d0*div_h3(a)

                    values(next) = values(next) - 5.0d0*div_h3(a)

                    values(next_next) = values(next_next) + 4.0d0*div_h3(a)

                    values(next_next_next) = values(next_next_next) - div_h3(a)

                    nnz = nnz + 3 ! three next

                case (2) ! Endpoint (values for backward difference)

                    values(prev_prev_prev) = values(prev_prev_prev) - div_h3(a)

                    values(prev_prev) = values(prev_prev) + 4.0d0*div_h3(a)

                    values(prev) = values(prev) - 5.0d0*div_h3(a)

                    values(center) = values(center) + 2.0d0*div_h3(a)

                    nnz = nnz + 3 ! three prev

            end select
        end do

        b(g) = 0.0d0

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
        integer(kind=8) :: x_1h, x_2h
        double precision :: forward_difference

        x_1h = move(x,1,3)
        x_2h = move(x,2,3)

        forward_difference = (-3.0d0*f(x)+4.0d0*f(x_1h)-f(x_2h))/(2.0d0*hz)
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

! -----------------------------------------------------------------------------
! ----- Clean Up --------------------------------------------------------------
! -----------------------------------------------------------------------------
    
    subroutine Clean_Up_Laplace() ! DONE
        ! Clean up memory
        deallocate(values, col_index, rows_start, rows_end, b, iCharge, oCharge, gridPoints, gridPointsActive)
        deallocate(voltage, laplace_field)
    end subroutine Clean_Up_Laplace

    subroutine check_matrix()
        integer(kind=8) :: i,count,low,high,col,index


        index = 1
        do i=1,nrGridActive
            count = 0
            low = nnz_rows_start(i)
            high = nnz_rows_end(i)
            col = nnz_col_index(index)

            do while ((low <= col) .and. (col < high))
                index = index + 1
                col = nnz_col_index(index)
                count = count + 1
            end do
            count  = count + 1
            index = index + 1

            print *, 'Laplace: row=',i,' count=',count

        end do


    end subroutine check_matrix

    subroutine Place_Electron(step)
        integer(kind=8), intent(in) :: step
        double precision, dimension(3) :: par_pos, par_vel

        if (step == 1) then
            par_pos(1) = -2.0d0*length_scale
            par_pos(2) = -2.0d0*length_scale
            par_pos(3) = 4.0d0*length_scale 
            par_vel = 0.0d0
            call Add_Particle(par_pos,par_vel,species_elec,step,1,-1)
            ! par_pos(1) = 4.0d0*length_scale
            ! par_pos(2) = 4.0d0*length_scale
            ! par_pos(3) = 1.5d0*length_scale 
            ! par_vel = 0.0d0
            ! call Add_Particle(par_pos,par_vel,species_elec,step,1,-1)
        end if
    end subroutine Place_Electron

    subroutine Write_Laplace_Data()
        integer :: ud_laplace_voltage, ud_laplace_field, IFAIL, lvl, i, j, k
        character(len=128) :: filename_field, filename_voltage
        integer :: x, y
        double precision, dimension(3) :: test_pos
        
        do lvl=0,2
            write(filename_voltage, '(a20,i0,a4)') 'out/laplace_voltage_',lvl,'.bin'
            write(filename_field, '(a18,i0,a4)') 'out/laplace_field_',lvl,'.bin'

            ! Open the voltage file
            open(newunit=ud_laplace_voltage, iostat=IFAIL, file=filename_voltage, status='REPLACE', action='WRITE', access='STREAM')
            if (IFAIL /= 0) then
                print *, 'RUMDEED: Failed to open the laplace voltage file.'
                return 
            end if

            ! Open the field file
            open(newunit=ud_laplace_field, iostat=IFAIL, file=filename_field, status='REPLACE', action='WRITE', access='STREAM')
            if (IFAIL /= 0) then
                print *, 'RUMDEED: Failed to open the laplace field file.'
                return 
            end if

            do i=0,Nx-1
                do j=0,Ny-1
                    k = i + j*Nx+1 + lvl*NxNy
                    test_pos = cart_coord(k)
                    ! print *, 'Laplace: writing: x=',test_pos(1), ' y=',test_pos(2), ' z=',test_pos(3)
                    write(unit=ud_laplace_voltage,iostat=IFAIL) i, j, voltage(k), is_emitter(k)
                    write(unit=ud_laplace_field,iostat=IFAIL) i, j, laplace_field(k), is_emitter(k)
                end do
            end do

            close(unit=ud_laplace_voltage, iostat=IFAIL, status='keep')
            close(unit=ud_laplace_field, iostat=IFAIL, status='keep')
        end do
        
    end subroutine Write_Laplace_Data

end module mod_laplace_solver