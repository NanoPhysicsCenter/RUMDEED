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

    implicit none

    ! Accessibility
    private
    public :: Init_Laplace_Solver, Calculate_Laplace_Field, Clean_Up_Laplace

    ! Local variables
    integer, dimension(64) :: pt
    integer, dimension(64) :: iparm

    double precision :: hx, hy, hz, emitter_radius
    double precision, dimension(3) :: div_h2, div_2_h2
    double precision, dimension(2) :: lim_x, lim_y, lim_z
    integer :: Nx, Ny, Nz, nrGrid, nrGridActive, nnz, NxNy, n7

contains

! -----------------------------------------------------------------------------
! ----- Solver ----------------------------------------------------------------
! -----------------------------------------------------------------------------

    subroutine Calculate_Laplace_Field() ! TODO
        double precision, dimension(nrGridActive) :: x
        ! Calculate the electric field
        print *, 'Laplace: finding electrons'
        call find_electrons()
        print *, 'Laplace: updating matrix'
        call update_matrix()
        print *, 'Laplace: solving matrix'
        call solve_matrix(x)
        
    end subroutine Calculate_Laplace_Field

    subroutine solve_matrix(x) ! TODO
        ! Solve the system
        type(sparse_matrix_t) :: cooA
        ! type(sparse_matrix_t) :: csrA
        double precision, dimension(nrGridActive), intent(out) :: x
        integer, dimension(nrGridActive) :: perm
        integer :: maxfct=1,mnum=1,mtype=11,phase=11
        integer :: msglvl, error
        x = 0.0d0
        print *, 'Laplace: creating COO'
        call create_coo(cooA)
        print *, 'Laplace: converting COO to CSR'
        ! call coo_to_csr(cooA,csrA)

        ! ! print *, 'Laplace: solving matrix'
        ! call pardiso_d(pt,maxfct,mnum,mtype,phase,nrGridActive,csrA, perm,nrGridActive,iparm,msglvl,b,x,error)
    end subroutine solve_matrix

! -----------------------------------------------------------------------------
! ----- Initialization --------------------------------------------------------
! -----------------------------------------------------------------------------

    subroutine Init_Laplace_Solver() ! DONE
        ! Initialize environment based on input file
        double precision :: Lx, Ly, Lz

        Nx = laplace_intervals(1); Ny = laplace_intervals(2); Nz = laplace_intervals(3)
        Lx = laplace_dim(1); Ly = laplace_dim(2); Lz = laplace_dim(3)

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

        print *, 'Laplace: allocating arrays'
        allocate(linkA(n7), valA(n7),  rowA(n7), colA(n7), b(nrGridActive), iCharge(nrGridActive), oCharge(nrGridActive))

        nnz = 0
        print *, 'Laplace: initializing matrix'
        call init_matrix()

        print *, 'Laplace: initializing Pardiso'
        call init_pardiso()

    end subroutine Init_Laplace_Solver

    subroutine init_pardiso()
        integer :: mtype=11
        ! call pardisoinit(pt,mtype,iparm)
    end subroutine init_pardiso

    subroutine init_grid()
        ! Store points that are outside of the emitter or one layer within it
        integer :: i
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
        integer :: g, i, axis, index
        integer, dimension(3) :: boundaries

        rowA = rowA*0
        colA = colA*0
        valA = valA*0.0d0
        b = b*0.0d0
        index = 1

        do g=1,nrGridActive
            i = gridPointsActive(g)
            
            ! print *, "Laplace: inputting gridpoint",i
            if (is_emitter(i) == 1) then
                rowA(index) = g
                colA(index) = g

                valA(index) = 1.0d0
                nnz = nnz + 1 ! Add 1 because emitter needs only 1 value

                linkA(index) = index+6
                index = index+6
            else
                ! Define 7 cells with this row index
                rowA(index:index+6) = (/g,g,g,g,g,g,g/)
                colA(index) = g

                linkA(index) = index+1
                index = index+1

                boundaries = is_boundary(i)

                ! Define the columns for each cell in accordance with point type 
                do axis=1,3
                    select case (boundaries(axis))
                        case (0) ! Inner point (central difference)
                            colA(index) = gridPoints(move(i,-1,axis))
                            linkA(index) = index+1
                            index = index + 1
                            colA(index) = gridPoints(move(i,1,axis))
                            linkA(index) = index+1
                            index = index+1
                        case (1) ! Startpoint (forward difference)
                            colA(index) = gridPoints(move(i,1,axis))
                            linkA(index) = index+1
                            index = index + 1
                            colA(index) = gridPoints(move(i,2,axis))
                            linkA(index) = index+1
                            index = index + 1
                        case (2) ! Endpoint (backward difference)
                            colA(index) = gridPoints(move(i,-2,axis))
                            linkA(index) = index+1
                            index = index + 1
                            colA(index) = gridPoints(move(i,-1,axis))  
                            linkA(index) = index+1  
                            index = index + 1                    
                    end select
                end do
                ! Define values for each cell in accordance with point type
                call insert_finite_difference(g)
                nnz = nnz + 7 ! Add 7 because this is the number of values added by finite difference
            end if
        end do
    end subroutine init_matrix

! -----------------------------------------------------------------------------
! ----- Matrix operations -----------------------------------------------------
! -----------------------------------------------------------------------------

    subroutine update_matrix() ! DONE
        ! Update boundary conditions in matrix
        integer :: g
        oCharge = iCharge

        ! $OMP PARALLEL DO DEFAULT(NONE) &
        ! $OMP& SHARED(nrGridActive, oCharge, iCharge, b, valA) &
        ! $OMP& PRIVATE(g) &
        ! $OMP& REDUCTION(+:nnz)
        do g=1,nrGridActive
            if (oCharge(g) /= 0.0d0) then
                if (iCharge(g) == 0.0d0) then
                    ! $OMP CRITICAL
                    call insert_finite_difference(g)
                    ! $OMP END CRITICAL
                    nnz = nnz + 6 ! Add only 6 because one was being accounted for by oCharge
                    b(g) = 0.0d0
                end if
            else
                if (iCharge(g) /= 0.0d0) then
                    ! $OMP CRITICAL
                    call remove_finite_difference(g)
                    ! $OMP END CRITICAL
                    nnz = nnz - 6 ! Remove only 6 because one will be used for by iCharge
                end if
            end if

            if (iCharge(g) /= 0.0d0) then
                valA(1+(g-1)*7) = 1.0d0
                b(g) = iCharge(g)
            end if
        end do
        ! $OMP END PARALLEL DO
    end subroutine update_matrix

    subroutine insert_finite_difference(g) ! DONE
        ! Insert finite difference equation into matrix
        integer, intent(in) :: g
        integer, dimension(3) :: boundaries
        integer :: i, a, index

        i = gridPointsActive(g)

        index = 1 + (g-1)*7
        boundaries = is_boundary(i)

        ! Insert finite difference equation for each axis
        do a=1,3
            select case (boundaries(a))
                case (0) ! Inner point (values for central difference)
                    valA(index) = valA(index) - div_2_h2(a)
                    valA(index+a) = valA(index+1) + div_h2(a)
                    valA(index+a+1) = valA(index+2) - div_2_h2(a)
                case (1) ! Startpoint (values for forward difference)
                    valA(index) = valA(index) + div_h2(a)
                    valA(index+a) = valA(index+1) - div_2_h2(a)
                    valA(index+a+1) = valA(index+2) + div_h2(a)
                case (2) ! Endpoint (values for backward difference)
                    valA(index) = valA(index) + div_h2(a)
                    valA(index+a) = valA(index+1) + div_h2(a)
                    valA(index+a+1) = valA(index+2) - div_2_h2(a)
            end select
        end do
    end subroutine insert_finite_difference

    subroutine remove_finite_difference(g) ! DONE
        ! Remove finite difference equation from matrix
        integer, intent(in) :: g
        integer, dimension(3) :: boundaries
        integer :: i, a, index

        i = gridPointsActive(g)

        index = 1 + (g-1)*7
        boundaries = is_boundary(i)

        ! remove finite difference equation for each axis
        do a=1,3
            select case (boundaries(a))
            case (0) ! Inner point (values for central difference)
                valA(index) = 0.0d0
                valA(index+a) = 0.0d0
                valA(index+a+1) = 0.0d0
            case (1) ! Startpoint (values for forward difference)
                valA(index) = 0.0d0
                valA(index+a) = 0.0d0
                valA(index+a+1) = 0.0d0
            case (2) ! Endpoint (values for backward difference)
                valA(index) = 0.0d0
                valA(index+a) = 0.0d0
                valA(index+a+1) = 0.0d0
            end select
        end do
    end subroutine remove_finite_difference

! -----------------------------------------------------------------------------
! ----- Boundary conditions ---------------------------------------------------
! -----------------------------------------------------------------------------

    function elec_boundary(elec_pos, i) ! DONE
        ! Calculate the boundary value for an electron
        double precision, dimension(3), intent(in) :: elec_pos
        integer, intent(in) :: i
        double precision, dimension(3) :: grid_coord
        double precision :: elec_boundary
        
        grid_coord = cart_coord(i)
        elec_boundary = q_0/((4*pi*epsilon_0) * 1/norm2(grid_coord - elec_pos))
    end function elec_boundary

    subroutine find_electrons() ! DONE
        ! Aggregate boundary values from electrons
        double precision, dimension(3) :: elec_pos
        integer :: g, i, k
        
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

    subroutine create_coo(cooA) ! TODO
        ! Create COO format matrix from matrix arrays
        type(sparse_matrix_t), intent(inout) :: cooA
        integer(4) :: stat
        integer :: i, index
        integer, dimension(nnz) :: nnz_rowA, nnz_colA
        double precision, dimension(nnz) :: nnz_valA

        print *, 'Laplace: aggregating COO'
        i = 1
        index = 1
        ! $OMP PARALLEL DO DEFAULT(NONE) &
        ! $OMP& SHARED(nnz, nnz_rowA, nnz_colA, nnz_valA, rowA, colA, valA, linkA) &
        ! $OMP& PRIVATE(i,index)
        do i=1,nnz
            if (index > 7*nrGridActive) then
                print *, 'Laplace: index out of bounds'
                print *, 'Laplace: index=',index
                exit
            end if
            nnz_rowA(i) = rowA(index)
            nnz_colA(i) = colA(index)
            nnz_valA(i) = valA(index)
            index = linkA(index)
        end do
        ! $OMP END PARALLEL DO  
        ! print *, 'Laplace: creating COO'
        print *, 'Laplace: nrGridActive=',nrGridActive
        print *, 'Laplace: nnz=',nnz
        ! print *, 'Laplace: size of nnz_rowA=',size(nnz_rowA)
        ! print *, 'Laplace: size of nnz_colA=',size(nnz_colA)
        ! print *, 'Laplace: size of nnz_valA=',size(nnz_valA)
        ! print *, 'Laplace: max of nnz_rowA=',maxval(nnz_rowA)
        ! print *, 'Laplace: max of nnz_colA=',maxval(nnz_colA)
        ! print *, 'Laplace: min of nnz_rowA=',minval(nnz_rowA)
        ! print *, 'Laplace: min of nnz_colA=',minval(nnz_colA)
        ! stat = mkl_sparse_d_create_coo(cooA, SPARSE_INDEX_BASE_ONE, nrGridActive, nrGridActive, nnz, nnz_rowA, nnz_colA, nnz_valA)

    end subroutine create_coo

    subroutine coo_to_csr(cooA, csrA) ! TODO
        ! Convert COO to CSR format
        type(sparse_matrix_t), intent(in) :: cooA
        type(sparse_matrix_t), intent(out) :: csrA
        integer :: status
        status = mkl_sparse_convert_csr(cooA, 1, csrA)
    end subroutine coo_to_csr

! -----------------------------------------------------------------------------
! ----- Indexing functions ----------------------------------------------------
! -----------------------------------------------------------------------------

    function disc_coord(x,y,z) ! DONE
        ! Change cartesian coordinates to discrete coordinates
        double precision, intent(in) :: x, y, z
        integer :: disc_coord

        disc_coord = x / hx + (y / hy)*Nx + (z / hz)*NxNy + 1
    end function disc_coord

    function cart_coord(i) ! DONE
        ! Change discrete coordinates to cartesian coordinates
        integer, intent(in) :: i
        integer :: k
        double precision, dimension(3) :: cart_coord
        k = i-1
        cart_coord(1) = mod(k,Nx) * hx + lim_x(1)
        cart_coord(2) = mod(k,NxNy)/Nx * hy + lim_y(1)
        cart_coord(3) = (k / (NxNy)) * hz + lim_z(1)
    end function cart_coord

    function move(i,nr_steps,axis) ! DONE
        integer, intent(in) :: i, nr_steps, axis
        integer :: move

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
        integer, intent(in) :: i
        integer, dimension(3) :: is_boundary

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
        integer, intent(in) :: i
        integer :: is_emitter
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
        integer, intent(in) :: i
        integer :: axis, is_first_layer, nr_surrounding


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
        deallocate(linkA, valA, rowA, colA, b, iCharge, oCharge, gridPoints, gridPointsActive)
    end subroutine Clean_Up_Laplace


end module mod_laplace_solver