! ------------------------------------------------------------------------------
! ----- Module: mod_laplace_solver ---------------------------------------------
! ----- Author: Arnar Jonsson --------------------------------------------------
! ----- Date: 2025 -------------------------------------------------------------
! ------------------------------------------------------------------------------

module mod_laplace_solver
    
    ! Dependencies
    use mod_global
    use mkl_pardiso ! Required for pardiso solver
    use mkl_spblas  ! Required for sparse matrices
    use iso_c_binding, only: c_int, c_ptr, c_f_pointer

    implicit none

    ! Accessibility
    private
    public :: Init_Laplace_Solver, Calculate_Laplace_Field, Clean_Up_Laplace

    ! Local variables
    type(MKL_PARDISO_HANDLE), dimension(64) :: pt
    ! integer(8), dimension(64) :: pt
    integer, dimension(64) :: iparm

    ! Matrix arrays
    double precision, allocatable, dimension(:) :: values, iCharge, oCharge, b
    integer(kind=8), allocatable, dimension(:) :: rows_start, rows_end, col_index

    double precision, allocatable, dimension(:) :: nnz_values
    integer(kind=8), allocatable, dimension(:) :: nnz_rows_start, nnz_rows_end, nnz_col_index, nnz_ia

    ! System parameters
    double precision :: hx, hy, hz, emitter_radius
    double precision, dimension(3) :: div_h2, div_2_h2
    double precision, dimension(2) :: lim_x, lim_y, lim_z
    integer :: Nx, Ny, Nz, nrGrid, nrGridActive, NxNy, n7, n13, nnz

contains

! -----------------------------------------------------------------------------
! ----- Solver ----------------------------------------------------------------
! -----------------------------------------------------------------------------

    subroutine Calculate_Laplace_Field() ! TODO
        double precision, dimension(nrGridActive) :: x
        double precision, dimension(nrGrid) :: voltage, field
        ! Calculate the electric field
        print *, 'Laplace: finding electrons'
        call find_electrons()
        print *, 'Laplace: updating matrix'
        call update_matrix()
        print *, 'Laplace: solving matrix'
        call solve_matrix(x)
        print *, 'Laplace: allocate voltage'
        call allocate_voltage(x,voltage)
        print *, 'Laplace: calculating field'
        call calculate_field(voltage,field)
        print *, 'Laplace: done'
        
    end subroutine Calculate_Laplace_Field

    subroutine solve_matrix(x) ! TODO
        ! Solve the system
        integer :: i
        type(sparse_matrix_t) :: cooA
        type(sparse_matrix_t) :: csrA
        double precision, dimension(nrGridActive), intent(inout) :: x
        integer(kind=8), dimension(nrGridActive) :: perm
        integer(kind=8) :: maxfct=1,mnum=1,mtype=11,phase=11
        integer(kind=8) :: msglvl=0, error


        x = 0.0d0
        ! print *, 'Laplace: creating COO'
        ! call create_coo(cooA)
        ! print *, 'Laplace: converting COO to CSR'
        ! call coo_to_csr(cooA,csrA)
        ! print *, 'Extracting CSR'
        ! call export_csr(csrA, nnz_valA, nnz_rowA, nnz_colA)
        ! print *, 'Solving matrix'
        ! iparm(60) = 2

        print *, 'Laplace: reduce matrix'
        allocate(nnz_values(nnz), nnz_col_index(nnz), nnz_ia(nrGridActive+1), nnz_rows_start(nrGridActive+1), nnz_rows_end(nrGridActive+1))
        call reduce_matrix()

        ! call check_matrix()
        ! iparm(60) = 2

        print *, 'Laplace: calling pardiso'
        call pardiso_64(pt,maxfct,mnum,mtype,phase,nrGridActive,nnz_values,nnz_ia,nnz_col_index,perm,nrGridActive,iparm,msglvl,b,x,error)
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

    subroutine allocate_voltage(x,voltage)
        double precision, dimension(nrGridActive), intent(in) :: x
        double precision, dimension(nrGrid), intent(out) :: voltage
        integer(kind=8) :: i, g

        voltage = 0.0d0

        do i=1,nrGrid
            if (gridPoints(i) /= 0) then
                g = gridPoints(i)
                voltage(i) = x(g)
            end if
        end do
    end subroutine allocate_voltage

    subroutine calculate_field(voltage, field)
        double precision, dimension(nrGrid), intent(out) :: field
        double precision, dimension(nrGrid), intent(in) :: voltage
        integer(kind=8) :: i, g, boundary
        integer(kind=8), dimension(3) :: boundaries

        field = 0.0d0

        do i=1,nrGrid
            boundaries = is_boundary(i)
            boundary = boundaries(3)
            select case (boundary)
                case (0) ! Inner point
                    field(i) = central_difference(voltage,i)
                case (1) ! Startpoint
                    field(i) = forward_difference(voltage,i)
                case (2) ! Endpoint
                    field (i) = backward_difference(voltage,i)
            end select
        end do

    end subroutine calculate_field

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
        div_2_h2 = (/2.0d0 / (hx**2), 2.0d0 / (hy**2), 2.0d0 / (hz**2)/)

        lim_x = (/laplace_pos(1), laplace_pos(1)+Lx/)
        lim_y = (/laplace_pos(2), laplace_pos(2)+Ly/)
        lim_z = (/laplace_pos(3), laplace_pos(3)+Lz/)

        emitter_radius = emitters_dim(1,1)

        allocate(gridPoints(nrGrid), gridPointsActive(nrGrid))
        print *, 'Laplace: initializing gridpoints'
        call init_grid()
        n7 = 7*nrGridActive
        n13 = 13*nrGridActive

        print *, 'Laplace: allocating arrays'
        allocate(values(n13), col_index(n13), rows_start(nrGridActive), rows_end(nrGridActive), b(nrGridActive), iCharge(nrGridActive), oCharge(nrGridActive))

        print *, 'Laplace: initializing matrix'
        call init_matrix()

        print *, 'Laplace: initializing Pardiso'
        call init_pardiso()
        
    end subroutine Init_Laplace_Solver

    subroutine init_pardiso()
        integer :: mtype=11
        call pardisoinit(pt,mtype,iparm)
    end subroutine init_pardiso

    subroutine init_grid()
        ! Store points that are outside of the emitter or one layer within it
        integer(kind=8) :: i
        do i=1,nrGrid
            if (is_emitter(i) == 1) then
                ! print *, 'Laplace: checking emitter'
                if (is_first_layer(i) == 1) then
                    ! print *, 'Laplace: checking first layer'
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
        integer(kind=8) :: next, next_next, prev, prev_prev

        values = values*0.0d0
        b = b*0.0d0

        center = 7

        do g=1,nrGridActive

            col_index(center) = g
            do a=1,3
                next = center + (a-1)*2 + 1
                next_next = next + 1
                prev = center - (a-1)*2 - 1
                prev_prev = prev - 1

                col_index(prev_prev) = move(g,-2,a)
                col_index(prev) = move(g,-1,a)
                col_index(next) = move(g,1,a)
                col_index(next_next) = move(g,2,a)
            end do

            rows_start(g) = g
            rows_end(g) = g

            i = gridPointsActive(g)
            
            if (is_emitter(i) == 1) then
                call insert_emitter_boundary(g)
                nnz = nnz + 1
            else
                call insert_finite_difference(g)
                nnz = nnz + 7
            end if

            center = center+13

        end do
    end subroutine init_matrix

! -----------------------------------------------------------------------------
! ----- Matrix operations -----------------------------------------------------
! -----------------------------------------------------------------------------

    subroutine update_matrix() ! DONE
        ! Update boundary conditions in matrix
        integer(kind=8) :: g

        ! $OMP PARALLEL DO DEFAULT(NONE) &
        ! $OMP& SHARED(nrGridActive, oCharge, iCharge, b, valA) &
        ! $OMP& PRIVATE(g) &
        ! $OMP& REDUCTION(+:nnz)
        do g=1,nrGridActive
            if (oCharge(g) /= 0.0d0) then
                if (iCharge(g) == 0.0d0) then
                    ! $OMP CRITICAL
                    call insert_finite_difference(g)
                    nnz = nnz + 6
                    ! $OMP END CRITICAL
                    b(g) = 0.0d0
                end if
            else
                if (iCharge(g) /= 0.0d0) then
                    ! $OMP CRITICAL
                    call remove_finite_difference(g)
                    nnz = nnz - 6
                    ! $OMP END CRITICAL
                end if
            end if

            if (iCharge(g) /= 0.0d0) then

                call insert_electron_boundary(g)
            end if
        end do
        ! $OMP END PARALLEL DO

        oCharge = iCharge
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

            do j=1,13
                cur_ind = cur_ind + 1
                cur_col = col_index(cur_ind)
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
        end do

        nnz_ia(nrGridActive+1) = nnz

    end subroutine reduce_matrix

    subroutine insert_electron_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = (g-1)*13 + 7

        values(center) = 1.0d0
        b(g) = iCharge(g)
    end subroutine

    subroutine insert_emitter_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = (g-1)*13 + 7

        values(center) = 1.0d0
        b(g) = 0.0d0

    end subroutine

    subroutine insert_finite_difference(g) ! DONE
        ! Insert finite difference equation into matrix
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center
        integer(kind=8), dimension(3) :: boundaries
        integer(kind=8) :: i, a, next, next_next, prev, prev_prev

        center = (g-1)*13 + 7

        i = gridPointsActive(g)

        boundaries = is_boundary(i)

        ! Insert finite difference equation for each axis
        do a=1,3

            next = center + (a-1)*2 + 1
            next_next = next + 1
            prev = center - (a-1)*2 - 1
            prev_prev = prev - 1


            select case (boundaries(a))
                case (0) ! Inner point (values for central difference)
                    values(center) = values(center) - div_2_h2(a)

                    values(next) = values(next) + div_2_h2(a)

                    values(prev) = values(prev) - div_2_h2(a)

                case (1) ! Startpoint (values for forward difference)

                    values(center) = values(center) + div_h2(a)

                    values(next) = values(next) - div_h2(a)

                    values(next_next) = values(next_next) + div_h2(a)

                case (2) ! Endpoint (values for backward difference)
                    values(center) = values(center) + div_h2(a)

                    values(prev) = values(prev) + div_h2(a)

                    values(prev_prev) = values(prev_prev) - div_h2(a)

            end select
        end do

        b(g) = 0.0d0
    end subroutine insert_finite_difference

    subroutine remove_finite_difference(g) ! DONE
        ! Remove finite difference equation from matrix
        integer(kind=8), intent(in) :: g
        integer(kind=8), dimension(3) :: boundaries
        integer(kind=8) :: i, center

        center = (g-1)*13 + 7

        values(center) = 0.0d0
        do i=1,7
            values(center+i) = 0.0d0
            values(center-i) = 0.0d0
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

    central_difference = (array(index_next) - array(index_prev))/(2*hz)
end function central_difference

function forward_difference(array,index)
    double precision, dimension(nrGrid), intent(in) :: array
    integer(kind=8), intent(in) :: index
    integer(kind=8) :: index_next, index_next_next
    double precision :: forward_difference

    index_next = move(index,1,3)
    index_next_next = move(index,2,3)

    forward_difference = (-3*array(index)+4*array(index_next)-array(index_next_next))/(2*hz)
end function forward_difference

function backward_difference(array,index)
    double precision, dimension(nrGrid), intent(in) :: array
    integer(kind=8), intent(in) :: index
    integer(kind=8) :: index_prev, index_prev_prev
    double precision :: backward_difference

    index_prev = move(index,-1,3)
    index_prev_prev = move(index,-2,3)

    backward_difference = (array(index_prev_prev)-4*array(index_prev)+3*array(index))/(2*hz)
end function backward_difference



! -----------------------------------------------------------------------------
! ----- Boundary conditions ---------------------------------------------------
! -----------------------------------------------------------------------------

    function elec_boundary(elec_pos, i) ! DONE
        ! Calculate the boundary value for an electron
        double precision, dimension(3), intent(in) :: elec_pos
        integer(kind=8), intent(in) :: i
        double precision, dimension(3) :: grid_coord
        double precision :: elec_boundary
        
        grid_coord = cart_coord(i)
        elec_boundary = q_0/((4*pi*epsilon_0) * 1/norm2(grid_coord - elec_pos))
    end function elec_boundary

    subroutine find_electrons() ! DONE
        ! Aggregate boundary values from electrons
        double precision, dimension(3) :: elec_pos
        integer(kind=8) :: g, i, k
        
        iCharge = iCharge*0.0d0

        ! $OMP PARALLEL DO DEFAULT(NONE) &
        ! $OMP& SHARED(particles_mask, particles_cur_pos, iCharge, lim_x, lim_y, lim_z, nrPart, gridPoints) &
        ! $OMP& PRIVATE(elec_pos, g, i, k)
        do k=1,nrPart
            if (particles_mask(k) .eqv. .true.) then
                elec_pos = particles_cur_pos(k,:)
                if ((lim_x(1)>=elec_pos(1)) .and. (elec_pos(1)<=lim_x(2)) &
                    .and. (lim_y(1)>=elec_pos(2)) .and. (elec_pos(2)<=lim_y(2)) &
                    .and. (lim_z(1)>=elec_pos(2)) .and. (elec_pos(3)<=lim_z(2))) then

                    i = disc_coord(elec_pos(1), elec_pos(2), elec_pos(3))

                    ! Bottom side boundary values
                    g = gridPoints(i)
                    if (g /= 0) then 
                        ! $OMP CRITICAL
                        iCharge(g) = iCharge(g) + elec_boundary(elec_pos,i) 
                        ! $OMP END CRITICAL
                    end if
                    g = gridPoints(move(i,1,1))
                    if (g /= 0) then 
                        ! $OMP CRITICAL
                        iCharge(g) = iCharge(g) + elec_boundary(elec_pos,move(i,1,1)) 
                        ! $OMP END CRITICAL
                    end if
                    g = gridPoints(move(i,1,2))
                    if (g /= 0) then 
                        ! $OMP CRITICAL
                        iCharge(g) = iCharge(g) + elec_boundary(elec_pos,move(i,1,2)) 
                        ! $OMP END CRITICAL
                    end if
                    g = gridPoints(move(move(i,1,2),1,1))
                    if (g /= 0) then 
                        ! $OMP CRITICAL
                        iCharge(g) = iCharge(g) + elec_boundary(elec_pos,move(move(i,1,2),1,1)) 
                        ! $OMP END CRITICAL
                    end if
                    
                    ! Top side boundary values
                    g = gridPoints(move(i,1,3))
                    if (g /= 0) then 
                        ! $OMP CRITICAL
                        iCharge(g) = iCharge(g) + elec_boundary(elec_pos,move(i,1,3)) 
                        ! $OMP END CRITICAL
                    end if
                    g = gridPoints(move(move(i,1,1),1,3))
                    if (g /= 0 ) then 
                        ! $OMP CRITICAL
                        iCharge(g) = iCharge(g) + elec_boundary(elec_pos,move(move(i,1,1),1,3)) 
                        ! $OMP END CRITICAL
                    end if
                    g = gridPoints(move(move(i,1,2),1,3))
                    if (g /= 0) then 
                        ! $OMP CRITICAL
                        iCharge(g) = iCharge(g) + elec_boundary(elec_pos,move(move(i,1,2),1,3)) 
                        ! $OMP END CRITICAL
                    end if
                    g = gridPoints(move(move(move(i,1,2),1,1),1,3))
                    if (g /= 0) then 
                        ! $OMP CRITICAL
                        iCharge(g) = iCharge(g) + elec_boundary(elec_pos,move(move(move(i,1,2),1,1),1,3)) 
                        ! $OMP END CRITICAL
                    end if

                end if
            end if
        end do
        ! $OMP END PARALLEL DO
    end subroutine find_electrons

! -----------------------------------------------------------------------------
! ----- Matrix conversion -----------------------------------------------------
! -----------------------------------------------------------------------------

    subroutine create_csr(csrA)
        type(sparse_matrix_t), intent(inout) :: csrA
        integer :: status

        status = mkl_sparse_d_create_csr(&
            csrA, &
            SPARSE_INDEX_BASE_ONE, &
            nrGridActive, &
            nrGridActive, &
            nnz_rows_start, &
            nnz_rows_end, &
            nnz_col_index, &
            nnz_values)

        select case (status)
            case (SPARSE_STATUS_SUCCESS)
                print *, 'Laplace: create_csr: successful'
            case (SPARSE_STATUS_NOT_INITIALIZED)
                print *, 'Laplace: create_csr: not initialized'
            case (SPARSE_STATUS_ALLOC_FAILED)
                print *, 'Laplace: create_csr: allocation failed'
            case (SPARSE_STATUS_INVALID_VALUE)
                print *, 'Laplace: create_csr: invalid value'
            case (SPARSE_STATUS_EXECUTION_FAILED)
                print *, 'Laplace: create_csr: execution failed'
            case (SPARSE_STATUS_INTERNAL_ERROR)
                print *, 'Laplace: create_csr: internal error'
            case (SPARSE_STATUS_NOT_SUPPORTED)
                print *, 'Laplace: create_csr: not supported'
            case default
                print *, 'Laplace: create_csr: unknown error'
        end select
    end subroutine create_csr

    ! subroutine create_coo(cooA) ! TODO
    !     ! Create COO format matrix from matrix arrays
    !     type(sparse_matrix_t), intent(inout) :: cooA
    !     integer(4) :: status
    !     integer :: i, index
    !     integer, dimension(nnz) :: nnz_rowA, nnz_colA
    !     double precision, dimension(nnz) :: nnz_valA

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
    !     ! $OMP END PARALLEL DO  
    !     ! print *, 'Laplace: creating COO'

    !     status = mkl_sparse_d_create_coo(cooA, SPARSE_INDEX_BASE_ONE, nrGridActive, nrGridActive, nnz, nnz_rowA, nnz_colA, nnz_valA)

    !     select case (status)
    !         case (SPARSE_STATUS_SUCCESS)
    !             print *, 'Laplace: create_coo: successful'
    !         case (SPARSE_STATUS_NOT_INITIALIZED)
    !             print *, 'Laplace: create_coo: not initialized'
    !         case (SPARSE_STATUS_ALLOC_FAILED)
    !             print *, 'Laplace: create_coo: allocation failed'
    !         case (SPARSE_STATUS_INVALID_VALUE)
    !             print *, 'Laplace: create_coo: invalid value'
    !         case (SPARSE_STATUS_EXECUTION_FAILED)
    !             print *, 'Laplace: create_coo: execution failed'
    !         case (SPARSE_STATUS_INTERNAL_ERROR)
    !             print *, 'Laplace: create_coo: internal error'
    !         case (SPARSE_STATUS_NOT_SUPPORTED)
    !             print *, 'Laplace: create_coo: not supported'
    !         case default
    !             print *, 'Laplace: create_coo: unknown error'
    !     end select

    ! end subroutine create_coo

    ! subroutine coo_to_csr(cooA, csrA) ! TODO
    !     ! Convert COO to CSR format
    !     type(sparse_matrix_t), intent(in) :: cooA
    !     type(sparse_matrix_t), intent(out) :: csrA
    !     integer :: status

    !     status = mkl_sparse_convert_csr(cooA, SPARSE_OPERATION_NON_TRANSPOSE, csrA)

    !     select case (status)
    !         case (SPARSE_STATUS_SUCCESS)
    !             print *, 'Laplace: coo_to_csr: successful'
    !         case (SPARSE_STATUS_NOT_INITIALIZED)
    !             print *, 'Laplace: coo_to_csr: not initialized'
    !         case (SPARSE_STATUS_ALLOC_FAILED)
    !             print *, 'Laplace: coo_to_csr: allocation failed'
    !         case (SPARSE_STATUS_INVALID_VALUE)
    !             print *, 'Laplace: coo_to_csr: invalid value'
    !         case (SPARSE_STATUS_EXECUTION_FAILED)
    !             print *, 'Laplace: coo_to_csr: execution failed'
    !         case (SPARSE_STATUS_INTERNAL_ERROR)
    !             print *, 'Laplace: coo_to_csr: internal error'
    !         case (SPARSE_STATUS_NOT_SUPPORTED)
    !             print *, 'Laplace: coo_to_csr: not supported'
    !         case default
    !             print *, 'Laplace: coo_to_csr: unknown error'
    !     end select

    ! end subroutine coo_to_csr

    subroutine export_csr(csrA)
        type(sparse_matrix_t), intent(in) :: csrA
        integer :: status, i
        integer(c_int) :: indexing

        integer :: nrows, ncols

        type(c_ptr) :: values_ptr, rows_start_ptr, rows_end_ptr, col_index_ptr

        print *, 'Laplace: calling mkl_sparse_d_export_csr'
        status = mkl_sparse_d_export_csr(&
            &   csrA,           &
            &   indexing,       &
            &   nrows,          &
            &   ncols,          &
            &   rows_start_ptr, &
            &   rows_end_ptr,   &
            &   col_index_ptr,  &
            &   values_ptr)

        select case (status)
            case (SPARSE_STATUS_SUCCESS)
                print *, 'Laplace: export_csr: successful'
            case (SPARSE_STATUS_NOT_INITIALIZED)
                print *, 'Laplace: export_csr: not initialized'
            case (SPARSE_STATUS_ALLOC_FAILED)
                print *, 'Laplace: export_csr: allocation failed'
            case (SPARSE_STATUS_INVALID_VALUE)
                print *, 'Laplace: export_csr: invalid value'
            case (SPARSE_STATUS_EXECUTION_FAILED)
                print *, 'Laplace: export_csr: execution failed'
            case (SPARSE_STATUS_INTERNAL_ERROR)
                print *, 'Laplace: export_csr: internal error'
            case (SPARSE_STATUS_NOT_SUPPORTED)
                print *, 'Laplace: export_csr: not supported'
            case default
                print *, 'Laplace: export_csr: unknown error'
        end select

        call c_f_pointer(rows_start_ptr, nnz_rows_start, (/nrows+1/))
        call c_f_pointer(rows_end_ptr, nnz_rows_end, (/nrows+1/))

        nnz = rows_end(nrows-1) - 1

        call c_f_pointer(col_index_ptr, nnz_col_index, (/nnz/))
        call c_f_pointer(values_ptr, nnz_values, (/nnz/))

    end subroutine export_csr

! -----------------------------------------------------------------------------
! ----- Indexing functions ----------------------------------------------------
! -----------------------------------------------------------------------------

    function disc_coord(x,y,z) ! DONE
        ! Change cartesian coordinates to discrete coordinates
        double precision, intent(in) :: x, y, z
        integer(kind=8) :: disc_coord

        disc_coord = x / hx + (y / hy)*Nx + (z / hz)*NxNy + 1
    end function disc_coord

    function cart_coord(i) ! DONE
        ! Change discrete coordinates to cartesian coordinates
        integer(kind=8), intent(in) :: i
        integer(kind=8) :: k
        double precision, dimension(3) :: cart_coord
        k = i-1
        cart_coord(1) = mod(k,Nx) * hx + lim_x(1)
        cart_coord(2) = mod(k,NxNy)/Nx * hy + lim_y(1)
        cart_coord(3) = (k / (NxNy)) * hz + lim_z(1)
    end function cart_coord

    function move(i,nr_steps,axis) ! DONE
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
        if (mod(mod(i-1,NxNy),Ny) == 0) then
            is_boundary(1) = 1
        else if (mod(mod(i-1,NxNy),Ny) == Nx-1) then
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

    function is_emitter(i)
        ! Returns 1 if point is within the emitter and 0 otherwise
        integer(kind=8), intent(in) :: i
        integer(kind=8) :: is_emitter
        double precision, dimension(3) :: point_coord
        double precision :: point_radius

        point_coord = cart_coord(i)
        point_radius = sqrt(point_coord(1)**2 + point_coord(2)**2)
        if ((point_radius <= emitter_radius) .and. (point_coord(3) <= emitters_pos(3,1))) then
            is_emitter = 1
        else
            is_emitter = 0
        end if

    end function

    function is_first_layer(i)
        ! Returns 1 if point is within the first layer of the emitter and 0 otherwise
        integer(kind=8), intent(in) :: i
        integer(kind=8) :: axis, is_first_layer, nr_surrounding


        nr_surrounding = 0

        if (is_emitter(i) == 1) then
            do axis=1,3
                if (is_emitter(move(i,1,axis)) == 1) then
                    nr_surrounding = nr_surrounding + 1
                end if
                if (is_emitter(move(i,-1,axis)) == 1) then
                    nr_surrounding = nr_surrounding + 1
                end if
            end do
        end if

        if (nr_surrounding == 6) then
            is_first_layer = 0
        else
            is_first_layer = 1
        end if

    end function is_first_layer

! -----------------------------------------------------------------------------
! ----- Clean Up --------------------------------------------------------------
! -----------------------------------------------------------------------------
    
    subroutine Clean_Up_Laplace() ! DONE
        ! Clean up memory
        deallocate(values, col_index, rows_start, rows_end, b, iCharge, oCharge, gridPoints, gridPointsActive)
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

end module mod_laplace_solver