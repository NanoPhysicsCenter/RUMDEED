! ------------------------------------------------------------------------------
! ----- Module: mod_polarso ---------------------------------------------
! ----- Author: Arnar Jonsson --------------------------------------------------
! ----- Date: 2025 -------------------------------------------------------------
! ------------------------------------------------------------------------------

! The POLARSO Laplace solver is a compile-time option because it needs the
! Intel MKL library (PARDISO). Build with 'make POLARSO=yes' to enable it.
! Without it this file provides a stub module with the same interface (see
! the #else branch at the bottom), so the callers compile unchanged; the
! stubs can never be reached because use_polarso can only be set to .true.
! by Read_Polarso_Variables, which in the stub refuses to enable it.
#if defined(_POLARSO_)

module mod_polarso
    
    ! Dependencies
    use mod_global
    use mod_verlet
    use mod_pair
    use mod_work_function
    use mkl_pardiso ! Required for pardiso solver
    use mkl_spblas  ! Required for some sparse matrix operations
    use mkl_rci     ! Required for iterative solver

    implicit none

    ! Accessibility
    private
    public :: PL_Init_Solver, PL_Calculate_Field, PL_Update_Field, PL_Calculate_Field_At, &
                PL_Clean_Up, Place_Electron, Write_Polarso_Data, Read_Polarso_Variables

    integer(kind=8), parameter :: pardiso_precision = 8

    ! fgmres variables
    integer(kind=8), parameter :: fgmres_tol = 12
    integer(kind=8), parameter :: max_iter = 200
    integer(kind=8), dimension(128) :: ipar
    real(kind=8), dimension(128) :: dpar
    real(kind=8), allocatable, dimension(:) :: tmp
    integer(kind=8) :: rci_request

    ! Pardiso variables
    integer(kind=8), parameter :: switch_tol = 8
    integer(kind=8), parameter :: pardiso_tol = 12
    type(MKL_PARDISO_HANDLE), dimension(64) :: pt
    integer(kind=8), dimension(64) :: iparm

    integer(kind=8), parameter :: maxfct=1, mnum=1, mtype=11, nrhs=1, msglvl=1
    integer(kind=8), parameter :: analysis=11, factorization=22, solving=33, solving1=331, solving2=333

    integer(kind=8) :: ndiff
    integer(kind=8), allocatable, dimension(:) :: perm

    ! Matrix
    real(kind=pardiso_precision), allocatable, dimension(:) :: nnz_values
    integer(kind=8), allocatable, dimension(:) :: nnz_col_index, nnz_ia, nnz_center
    real(kind=pardiso_precision), allocatable, dimension(:) :: b
    integer(kind=8) :: nnz

    ! System
    integer(kind=8) :: newNrCharge, oldNrCharge
    integer(kind=8), allocatable, dimension(:,:) :: newCharge_index, oldCharge_index
    real(kind=8), allocatable, dimension(:,:) :: newCharge_density

    ! Solution
    real(kind=pardiso_precision), allocatable, dimension(:) :: x
    real(kind=8), allocatable, dimension(:) :: voltage
    real(kind=8), allocatable, dimension(:,:) :: polarso_field

    ! Discretization
    integer(kind=8), allocatable, dimension(:) :: gridPoints, gridPointsActive
    real(kind=8) :: hx, hy, hz
    real(kind=8) :: emitter_radius, min_d
    real(kind=8), dimension(3) :: h_grid ! Grid spacing (renamed from h to avoid clash with Planck's constant in mod_global)
    real(kind=8), dimension(3) :: div_h, div_h2, div_h3
    real(kind=8), dimension(3,2) :: elec_lim, grid_lim, emit_lim
    integer(kind=8) :: Nx, Ny, Nz, nrGrid, nrGridActive, NxNy
    logical :: anode_in_grid, cathode_in_grid

    ! Input file
    namelist /polarso/ use_polarso, polarso_dim, polarso_pos, &
                   polarso_step, polarso_padding, write_field_files

contains

    subroutine PL_Init_Solver() ! DONE
        ! Initialize environment based on input file

        print *, 'Polarso: initializing grid'
        call init_grid()

        print *, 'Polarso: initializing matrix'
        call init_matrix()

        print *, 'Polarso: initializing pardiso'
        call init_pardiso()

        ! call write_grid()

    end subroutine PL_Init_Solver

    subroutine Read_Polarso_Variables()
        integer :: ud_polarso, IFAIL, emit

        !Open the file 'polarso' for reading and check for errors.
        !The file is optional: without it the defaults apply (use_polarso = .false.).
        open(newunit=ud_polarso, iostat=IFAIL, file='polarso', status='OLD', action='read')
        if (IFAIL /= 0) then
        print '(a)', 'RUMDEED: No polarso file found, POLARSO disabled'
        return
        end if

        read(ud_polarso, NML=polarso)

        ! Close the 'polarso' file
        close(unit=ud_polarso, iostat=IFAIL, status='keep')

        polarso_dim = polarso_dim * length_scale
        polarso_pos = polarso_pos * length_scale
        polarso_padding = polarso_padding * length_scale
        polarso_step = polarso_step * length_scale
    end subroutine Read_Polarso_Variables

    subroutine PL_Calculate_Field()
        ! Note: write_field_files only writes diagnostic output; it must not
        ! change the simulation. Call Place_Electron manually (e.g. from a
        ! test driver) if a single test electron is wanted.
        print *, 'Polarso: finding electrons'
        call full_update_charge()
        print *, 'Polarso: updating matrix'
        call full_update_rhs()
        print *, 'Polarso: solving system'
        call solve_system()
        print *, 'Polarso: allocate voltage'
        call allocate_voltage()
        print *, 'Polarso: calculating field'
        call calculate_field()
        print *, 'Polarso: done'
        if (write_field_files .eqv. .true.) then
            call Write_Polarso_Data()
        end if

        ! call write_average_field()
        
    end subroutine PL_Calculate_Field

    subroutine PL_Update_Field(index)
        integer(kind=8), intent(in) :: index ! Index of the new particle in particle arrays
        integer(kind=8) :: charge_slot

        ! partial_update_charge claims a slot in the newCharge_* arrays; the
        ! RHS update must use that slot, not the particle index (the two only
        ! coincide when every particle so far was deposited).
        call partial_update_charge(index, charge_slot)
        if (charge_slot == 0) then
            return ! Particle outside the solver window: field unchanged
        end if
        call partial_update_rhs(charge_slot)
        call solve_system()
        call allocate_voltage()
        call calculate_field()

    end subroutine PL_Update_Field

! ------------------------------------------------------------------------------
! ----- Solution ---------------------------------------------------------------
! ------------------------------------------------------------------------------

    subroutine solve_system()
        ! Allocatable rather than automatic: nrGridActive can be large and
        ! automatic arrays live on the stack.
        real(kind=8), allocatable, dimension(:) :: sol_b
        real(kind=8) :: error, norm_b

        ! Attempt direct solution with pardiso
        ! print *, 'Polarso: attempting direct solver'
        call direct_solve_system()

        ! Compute relative error (absolute when b = 0, i.e. no space charge)
        sol_b = sparse_dot(nrGridActive, nnz_values, nnz_ia, nnz_col_index, x)
        norm_b = norm2(b)
        if (norm_b > 0.0d0) then
            error = norm2(sol_b-b)/norm_b
        else
            error = norm2(sol_b)
        end if

        ! Possibly use iterative solver
        if (error > 10.0d0**(-switch_tol)) then
            print *, 'Polarso: switching to iterative solver'
            call iterative_solve_system()
        end if
        
    end subroutine solve_system

    subroutine direct_solve_system()
        ! Solve the system
        call pardiso_phase(33,0)

    end subroutine direct_solve_system

    subroutine iterative_solve_system() ! TODO
        real(kind=8), allocatable, dimension(:) :: sol_b
        real(kind=8) :: error, norm_b
        integer(kind=8) :: par_error
        integer(kind=8) :: itercount
        integer(kind=8) :: n_restart
        ! Solve the system iteratively
        ! dfgmres_init sets ipar(15), the restart length, to min(150, n), so
        ! tmp can be sized to match before the init call (ipar(15) itself is
        ! undefined until dfgmres_init has run).
        n_restart = min(150_8, nrGridActive)
        allocate(tmp((2*n_restart+1)*nrGridActive + n_restart*(n_restart+9)/2 + 1))
        call init_dfgmres()
        ipar(5) = max_iter
        ipar(8) = 0
        ipar(9) = 1
        ipar(11) = 1

        ! print *, 'Polarso: checking dfgmres'
        call dfgmres_check(nrGridActive, x, b, rci_request, ipar, dpar, tmp)
        if (rci_request /= 0) then
            print *, '  Polarso: dfgmres_check failed with error code ', rci_request
        end if

        iterate: do while (.true.)
            ! print *, 'Polarso: DFGMRES iteration', ipar(4), 'of', ipar(5)
            call dfgmres(nrGridActive, x, b, rci_request, ipar, dpar, tmp)

            ! print *, 'size of tmp:', size(tmp)
            ! print *, 'ipar(23)+nrGridActive-1', ipar(23)+nrGridActive-1
            ! print *, 'ipar(22)+nrGridActive-1', ipar(22)+nrGridActive-1

            ! print *, 'Polarso: receiving feedback'
            select case(rci_request)
                case (0)
                    exit iterate
                case (1) ! Residual computation
                    ! print *, 'Polarso: dot product'
                    tmp(ipar(23):ipar(23)+nrGridActive-1) = sparse_dot(nrGridActive, nnz_values, nnz_ia, nnz_col_index, tmp(ipar(22):ipar(22)+nrGridActive-1))
                    cycle
                case (2) ! Stop condition
                    if (ipar(4) < ipar(5)) then
                        sol_b = sparse_dot(nrGridActive, nnz_values, nnz_ia, nnz_col_index, x)
                        norm_b = norm2(b)
                        if (norm_b > 0.0d0) then
                            error = norm2(sol_b-b)/norm_b
                        else
                            error = norm2(sol_b)
                        end if
                        if (error <= 10.0d0**(-fgmres_tol)) then
                            exit iterate
                        else
                            cycle
                        end if
                    else
                        exit iterate
                    end if
                case (3) ! Preconditioner
                    ! print *, 'Polarso: dfgmres: solving'
                    call pardiso_64(pt,maxfct,mnum,mtype,33,nrGridActive,nnz_values,nnz_ia,nnz_col_index, &
                    & perm,nrhs,iparm,0,tmp(ipar(22):ipar(22)+nrGridActive-1),tmp(ipar(23):ipar(23)+nrGridActive-1),par_error)
                case (4)
                    if (dpar(7) <= 10.0d0**(-16)) then
                        ! print *, 'Polarso: Reached minimum tolerance'
                        exit iterate
                    else
                        cycle
                    end if
                case default ! Error
                    print *, '  Polarso: dfgmres failed with error code ', rci_request
                    exit iterate
            end select
        end do iterate

        ! print *, 'Polarso: getting solution'
        call dfgmres_get(nrGridActive, x, b, rci_request, ipar, dpar, tmp, itercount)

        deallocate(tmp)

    end subroutine iterative_solve_system

! -------------------------------------------------------------------------------------------------------------
! ----- Initialization ----------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------------------


    subroutine init_grid() ! DONE
        ! Initialize parameters and store points that are outside of the emitter or one layer within it
        integer(kind=8) :: i, axis
        integer(kind=8), dimension(3) :: boundaries
        real(kind=8) :: Lx, Ly, Lz

        allocate(newCharge_index(8,max_particles), newCharge_density(8,max_particles), oldCharge_index(8,max_particles))


        Lx = polarso_dim(1)+2.0d0*polarso_padding; Ly = polarso_dim(2)+2.0d0*polarso_padding; Lz = polarso_dim(3)

        Nx = ceiling(Lx / polarso_step)+1; Ny = ceiling(Ly / polarso_step)+1; Nz = ceiling(Lz / polarso_step)+1
        Lx = (Nx-1)*polarso_step; Ly = (Ny-1)*polarso_step; Lz = (Nz-1)*polarso_step

        ! hx = Lx / (Nx-1)
        ! hy = Ly / (Ny-1)
        ! hz = Lz / (Nz-1)
        hx = Lx/(Nx-1)
        hy = Ly/(Ny-1)
        hz = Lz/(Nz-1)
        h_grid = (/hx, hy, hz/)

        ! Physical limits of the grid
        elec_lim(1,:) = (/polarso_pos(1), polarso_pos(1)+polarso_dim(1)/)
        elec_lim(2,:) = (/polarso_pos(2), polarso_pos(2)+polarso_dim(2)/)
        elec_lim(3,:) = (/polarso_pos(3), polarso_pos(3)+polarso_dim(3)/)
        ! print *, 'Polarso: physical limits:', elec_lim(1,:), elec_lim(2,:), elec_lim(3,:)

        ! print *, 'Polarso: physical limits:', elec_lim(1,:), elec_lim(2,:), elec_lim(3,:)

        ! Grid limits
        grid_lim(1,:) = (/polarso_pos(1)-polarso_padding, polarso_pos(1)-polarso_padding+Lx/)
        grid_lim(2,:) = (/polarso_pos(2)-polarso_padding, polarso_pos(2)-polarso_padding+Ly/)
        grid_lim(3,:) = (/polarso_pos(3), polarso_pos(3)+Lz/)

        ! print *, 'Polarso: grid limits:', grid_lim(1,:), grid_lim(2,:), grid_lim(3,:)

        ! Emitter limits
        emit_lim(1,:) = (/emitters_pos(1,1), emitters_pos(1,1)+emitters_dim(1,1)/)
        emit_lim(2,:) = (/emitters_pos(2,1), emitters_pos(2,1)+emitters_dim(2,1)/)
        emit_lim(3,:) = (/emitters_pos(3,1), emitters_pos(3,1)+emitters_dim(3,1)/)
        ! print *, 'Polarso: emitter limits:', emit_lim(1,:), emit_lim(2,:), emit_lim(3,:)

        ! The physical anode is the plane z = polarso_pos(3) + box_dim(3). The
        ! top grid plane carries anode (Dirichlet) rows only when it coincides
        ! with that plane (within half a cell, to be safe against roundoff in
        ! hz and against polarso_dim(3) being rounded up to a whole number of
        ! steps). Otherwise the solver window ends below the anode and the top
        ! plane gets open-boundary rows like the lateral boundaries.
        anode_in_grid = (abs(grid_lim(3,2) - (polarso_pos(3) + box_dim(3))) < 0.5d0*hz)

        ! Likewise the bottom grid plane is the grounded cathode plane (the
        ! conductor the emitter sits on, z = 0 in the planar model) only when
        ! the window starts there.
        cathode_in_grid = (abs(grid_lim(3,1)) < 0.5d0*hz)

        ! print *, 'Polarso: emitter limits:', emit_lim(1,:), emit_lim(2,:), emit_lim(3,:)

        ! emit_lim = grid_lim

        NxNy = Nx*Ny
        nrGrid = Nx*Ny*Nz

        div_h = (/1.0d0, 1.0d0, 1.0d0/)
        ! Stencil weights of the finite-difference Laplacian with the rows
        ! multiplied through by hz**2, so cubic grids keep O(1) matrix entries
        ! next to the O(1) Dirichlet rows. The RHS deposition in
        ! partial_update_charge carries the matching hz**2 factor.
        div_h2 = (/(hz*hz)/(hx*hx), (hz*hz)/(hy*hy), 1.0d0/)

        emitter_radius = emitters_dim(1,1)

        allocate(gridPoints(nrGrid), gridPointsActive(nrGrid))

        nrGridActive = 0
        nnz = 0
        min_d = length_scale**4

        ! $OMP PARALLEL DO DEFAULT (NONE) &
        ! $OMP& SHARED(nrGrid, gridPoints, gridPointsActive) &
        ! $OMP& PRIVATE(i,a,boundaries) &
        ! $OMP& REDUCTION(+:nrGridActive, nnz)
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
            else if ((is_anode(i) == 1) .or. (is_cathode(i) == 1)) then
                nrGridActive = nrGridActive + 1
                gridPointsActive(nrGridActive) = i
                gridPoints(i) = nrGridActive

                nnz = nnz + 1
            else
                nrGridActive = nrGridActive + 1
                gridPointsActive(nrGridActive) = i
                gridPoints(i) = nrGridActive

                boundaries = is_boundary(i)

                nnz = nnz+1

                do axis=1,3
                    select case (boundaries(axis))
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
        ! $OMP END PARALLEL DO

        ! print *, 'nrGrid', nrGrid
        ! print *, 'nrGridActive', nrGridActive
        allocate(polarso_field(nrGrid,3))

    end subroutine init_grid

    subroutine init_matrix() ! DONE
        ! Initialize matrix without boundary conditions
        integer(kind=8) :: g, i, axis, center, next_center
        integer(kind=8), dimension(3) :: boundaries
        integer(kind=8), dimension(3,3) :: indices
        integer(kind=8) :: next1, next2, next3, prev1, prev2, prev3

        allocate(x(nrGridActive), b(nrGridActive))
        allocate(nnz_values(nnz), nnz_col_index(nnz), nnz_ia(nrGridActive+1))
        allocate(nnz_center(nrGridActive))

        nnz_values = 0.0d0
        b = 0.0d0
        newNrCharge = 0
        oldNrCharge = 0

        center = 1

        ! $OMP PARALLEL DO DEFAULT (NONE) &
        ! $OMP& SHARED(nrGridActive, gridPointsActive) &
        ! $OMP& SHARED(nnz_center, nnz_ia, nnz_col_index) &
        ! $OMP& PRIVATE(g, i, a, boundaries, indices) &
        ! $OMP& PRIVATE(next_center, next1, next2, next3, prev1, prev2, prev3) &
        ! $OMP& REDUCTION(+:nnz_values, b, center)
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
            else if ((is_anode(i) == 1) .or. (is_cathode(i) == 1)) then
                ! Place the center
                nnz_center(g) = center
                next_center = center + 1
                ! Place the column indices
                nnz_ia(g) = center
                nnz_col_index(center) = g
                ! Place the values
                call insert_anode_boundary(g)
            else
                boundaries = is_boundary(i)
                ! Place the center
                do axis=1,3
                    select case (boundaries(axis))
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
                do axis=1,3
                    select case (boundaries(axis))
                        case (0)
                            prev1 = indices(axis,1)
                            next1 = indices(axis,3)
            
                            nnz_col_index(prev1) = gridPoints(move(i,-1,axis))
                            nnz_col_index(next1) = gridPoints(move(i,1,axis))
                            next_center = next1 + 1

                            nnz_ia(g) = prev1
                        case (1)
                            next1 = indices(axis,1)
                            next2 = indices(axis,2)
                            next3 = indices(axis,3)
                            
                            nnz_col_index(next1) = gridPoints(move(i,1,axis))
                            nnz_col_index(next2) = gridPoints(move(i,2,axis))
                            nnz_col_index(next3) = gridPoints(move(i,3,axis))

                            next_center = next3 + 1
                        case (2)
                            prev1 = indices(axis,3)
                            prev2 = indices(axis,2)
                            prev3 = indices(axis,1)
                            
                            nnz_col_index(prev1) = gridPoints(move(i,-1,axis))
                            nnz_col_index(prev2) = gridPoints(move(i,-2,axis))
                            nnz_col_index(prev3) = gridPoints(move(i,-3,axis))  
 
                            nnz_ia(g) = prev3
                    end select
                end do
                ! Place the values
                call insert_laplace_equation(g)
            end if
            center = next_center
        end do
        ! $OMP END PARALLEL DO
        nnz_ia(nrGridActive+1) = nnz+1

        deallocate(nnz_center)

        ! print *, 'Polarso: nnz:', nnz
        ! print *, 'Polarso: size(nnz_values):', size(nnz_values)
        ! print *, 'Polarso: size(nnz_col_index):', size(nnz_col_index)
        ! print *, 'Polarso: size(nnz_ia):', size(nnz_ia)
        ! print *, 'Polarso: size(b):', size(b)
        
    end subroutine init_matrix

    subroutine init_pardiso() ! DONE
        allocate(perm(nrGridActive))
        perm = 0 ! Not used (iparm(5) = 0) but must not be undefined

        ! Initialize pardiso solver and do analysis, symbolic factorization, and numerical factorization
        if (pardiso_precision == 4) then
            iparm(28) = 1
        end if
        call pardisoinit(pt,mtype,iparm)

        ! Solver parameters
        iparm(1) = 1 ! use user-defined parameters
        iparm(2) = 2 ! 1 = minimum degree algorithm, 2 = metis, 3 = parallel metis
                     ! (1 segfaults inside MKL 2026.1's reordering; 2 is the MKL
                     ! default and the right choice for 3D grids anyway)
        ! iparm(3) = 0 ! reserved
        iparm(4) = 0 ! >0 = preconditioned CGS
        iparm(5) = 0 ! 0 = ignored; 1 = user supplied fill-in permutation; 2 = returns permutation
        iparm(6) = 0 ! 0 = write solution on x; 1 = write solution on b
        ! iparm(7) = 0 ! output number of iterative refinement steps
        iparm(8) =  20 ! maximum number of iterative refinement steps
        iparm(9) = pardiso_tol ! tolerance = 10**(-iparm(9))	
        iparm(10) = 18 ! pivoting perturbation: 13 = 10**-13; 8 = 10**-8
        iparm(11) = 1 ! 0 = no scaling; 1 = scaling
        iparm(12) = 0 ! 0 = solve linear system; 1 = solve conjugate transposed system, 2 = solve transposed system
        iparm(13) = 1 ! 0 = no symmetric weighted matching; 1 = symmetric weighted matching
        ! iparm(14) = 0 ! output number of perturbed pivots
        ! iparm(15) = 0 ! output peak memory on symbolic factorization
        ! iparm(16) = 0 ! output permanent memory on symbolic factorization
        ! iparm(17) = 0 ! output size of factors on numerical factorization and solving
        ! iparm(18) = 0 ! output number of non-zero elements in the factors
        ! iparm(19) = 0 ! output number of floating point operations in factorization
        iparm(20) = 0 ! output cgs diagnostics
        iparm(21) = 1 ! 0 = 1x1 diagonal pivoting; 1 = 1x1 and 2x2 pivoting
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
        iparm(31) = 0 ! 0 = no partial solve
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

        ! call cpu_time(start)
        ! print *, 'Polarso: symbolic factorization'
        call pardiso_phase(11,1)
        iparm(27) = 0
        print *, 'Polarso: required memory:', iparm(17)*1e-6, 'Gb'
        ! print *, 'Polarso: numerical factorization'
        call pardiso_phase(22,1)
        ! print *, 'Polarso: solving'
        call pardiso_phase(33,1)

    end subroutine init_pardiso

    subroutine init_dfgmres()
        call dfgmres_init(nrGridActive, x, b, rci_request, ipar, dpar, tmp)

        if (rci_request /= 0) then
            print *, '  Polarso: dfgmres_init failed with error code ', rci_request
        end if
    end subroutine init_dfgmres

    subroutine pardiso_phase(phase,msg)
        integer(kind=8), intent(in) :: phase, msg
        integer(kind=8) :: error

        call pardiso_64(pt,maxfct,mnum,mtype,phase,nrGridActive,nnz_values,nnz_ia,nnz_col_index,perm,nrhs,iparm,msg,b,x,error)

        if (error == 0) then
            ! print *, '  Polarso: analysis and symbolic factorization successful'
        else
            print *, '  Polarso: phase', phase, 'failed with error code ', error
        end if

    end subroutine pardiso_phase

! -------------------------------------------------------------------------------------------------------------
! ----- Voltage and field calculations ------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------------------

    subroutine allocate_voltage() ! DONE & PARALLEL
        integer(kind=8) :: i, g

        allocate(voltage(nrGrid))

        voltage = 0.0d0

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP& SHARED(nrGrid, voltage, x, gridPoints) &
        !$OMP& PRIVATE(i, g)
        ! print *, nrGrid
        do i=1,nrGrid
            ! print *, 'Polarso: i:', i
            if (gridPoints(i) /= 0) then
                ! print *, 'x(g)', x(g)
                g = gridPoints(i)
                ! print *, 'Polarso: g:', g
                voltage(i) = x(g)
            end if
        end do
        !$OMP END PARALLEL DO

    end subroutine allocate_voltage

    subroutine calculate_field() ! DONE & PARALLEL
        integer(kind=8) :: i, axis
        integer(kind=8), dimension(3) :: boundaries
        real(kind=8) :: E_applied

        ! Uniform vacuum field of the planar gap, superposed on the
        ! space-charge field from the solve (grounded electrodes, see
        ! insert_anode_boundary). Same expression as Set_Voltage in
        ! mod_verlet; V_d can change every time step.
        E_applied = -V_d/d

        polarso_field = 0.0d0

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP& SHARED(nrGrid, voltage, polarso_field, E_applied) &
        !$OMP& PRIVATE(i, axis, boundaries)
        do i=1,nrGrid
            boundaries = is_boundary(i)
            do axis=1,3
                select case (boundaries(axis))
                    case (0) ! Inner point
                        polarso_field(i,axis) = - central_difference(voltage,i,axis)
                    case (1) ! Startpoint
                        polarso_field(i,axis) = - forward_difference(voltage,i,axis)
                    case (2) ! Endpoint
                        polarso_field(i,axis) = - backward_difference(voltage,i,axis)
                end select
            end do
            polarso_field(i,3) = polarso_field(i,3) + E_applied
        end do
        !$OMP END PARALLEL DO

        ! Deallocate the voltage vector
        deallocate(voltage)

    end subroutine calculate_field

    subroutine get_finite_indices(g,indices)
        integer(kind=8), intent(in) :: g
        integer(kind=8), dimension(3,3), intent(out) :: indices
        integer(kind=8) :: i, axis, center, start, end
        integer(kind=8), dimension(3) :: boundaries

        i = gridPointsActive(g)
        boundaries = is_boundary(i)
        center = nnz_center(g)

        start = center
        end = center

        do axis=1,3
            select case (boundaries(axis))
                case (0)
                    start = start - 1
                    end = end + 1
                    indices(axis,1) = start
                    indices(axis,3) = end
                case (1)
                    end = end + 1
                    indices(axis,1) = end
                    end = end + 1
                    indices(axis,2) = end
                    end = end + 1
                    indices(axis,3) = end
                case (2)
                    start = start - 1
                    indices(axis,3) = start
                    start = start - 1
                    indices(axis,2) = start
                    start = start - 1
                    indices(axis,1) = start
            end select
        end do

    end subroutine

! -----------------------------------------------------------------------------
! ----- Matrix operations -----------------------------------------------------
! -----------------------------------------------------------------------------

    subroutine full_update_rhs()
        integer(kind=8) :: i
        ! Rebuild the RHS from scratch each time step. All boundary rows are
        ! homogeneous (emitter and anode grounded, see insert_anode_boundary),
        ! so only the deposited space charge contributes.
        b = 0.0d0
        oldNrCharge = newNrCharge
        oldCharge_index = 0
        do i=1,newNrCharge
            call partial_update_rhs(i)
        end do
    end subroutine full_update_rhs

    subroutine partial_update_rhs(i)
        integer(kind=8), intent(in) :: i
        integer(kind=8) :: g, j
        do j=1,8
            g = newCharge_index(j,i)
            b(g) = b(g) + newCharge_density(j,i)
            oldCharge_index(j,i) = g
        end do
    end subroutine partial_update_rhs

    subroutine partial_remove_rhs(i)
        integer(kind=8), intent(in) :: i
        integer(kind=8) :: g, j
        do j=1,8
            g = oldCharge_index(j,i)
            b(g) = 0.0d0
        end do
    end subroutine partial_remove_rhs

    subroutine insert_emitter_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = nnz_center(g)

        nnz_values(center) = 1.0d0
        b(g) = 0.0d0
    end subroutine insert_emitter_boundary

    subroutine insert_anode_boundary(g)
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center

        center = nnz_center(g)

        nnz_values(center) = 1.0d0
        ! Grounded electrode planes (anode, and the cathode plane outside the
        ! emitter): the solver computes only the space-charge potential, with
        ! every conductor at 0 V. The applied vacuum field -V_d/d is
        ! superposed analytically in calculate_field, which also keeps the
        ! field correct when the solver window is shorter than the gap. For
        ! the flat planar emitters this solver supports, the superposition is
        ! exact.
        b(g) = 0.0d0
    end subroutine insert_anode_boundary

    subroutine insert_laplace_equation(g)
        ! Interior points get the 7-point finite-difference Laplacian. Points
        ! on an open boundary (the lateral faces, and the top plane when the
        ! solver window ends below the anode) get a homogeneous Neumann
        ! closure instead: the sum over the boundary axes of the second-order
        ! one-sided first derivatives, set to zero. The previous scheme
        ! imposed the PDE itself with one-sided second differences on the
        ! boundary, which made the matrix numerically singular
        ! (cond ~ 1e12, with the near-null mode concentrated on the
        ! boundary edges); the Neumann closure gives cond ~ 1e3.
        integer(kind=8), intent(in) :: g
        integer(kind=8) :: center
        integer(kind=8), dimension(3) :: boundaries
        integer(kind=8), dimension(3,3) :: indices
        integer(kind=8) :: i, axis, next, next_next, prev, prev_prev

        center = nnz_center(g)

        i = gridPointsActive(g)
        boundaries = is_boundary(i)
        call get_finite_indices(g,indices)

        if (any(boundaries /= 0)) then
            ! Neumann closure row. The sign of each axis contribution is
            ! chosen so that every axis adds +3 to the diagonal. Slots of
            ! non-boundary axes keep their explicit zeros.
            do axis=1,3
                select case (boundaries(axis))
                    case (1) ! Startpoint: 3*f(x) - 4*f(x+h) + f(x+2h) = 0
                        next = indices(axis,1)
                        next_next = indices(axis,2)

                        nnz_values(center) = nnz_values(center) + 3.0d0
                        nnz_values(next) = nnz_values(next) - 4.0d0
                        nnz_values(next_next) = nnz_values(next_next) + 1.0d0

                    case (2) ! Endpoint: 3*f(x) - 4*f(x-h) + f(x-2h) = 0
                        prev = indices(axis,3)
                        prev_prev = indices(axis,2)

                        nnz_values(center) = nnz_values(center) + 3.0d0
                        nnz_values(prev) = nnz_values(prev) - 4.0d0
                        nnz_values(prev_prev) = nnz_values(prev_prev) + 1.0d0
                end select
            end do
        else
            ! Interior point: central second difference along every axis
            do axis=1,3
                prev = indices(axis,1)
                next = indices(axis,3)

                nnz_values(prev) = nnz_values(prev) + div_h2(axis)
                nnz_values(center) = nnz_values(center) - 2.0d0*div_h2(axis)
                nnz_values(next) = nnz_values(next) + div_h2(axis)
            end do
        end if

        b(g) = 0.0d0
    end subroutine insert_laplace_equation

! -----------------------------------------------------------------------------
! ----- Finite difference -----------------------------------------------------
! -----------------------------------------------------------------------------

    function first_finite_difference_3d(f,i)
        real(kind=8), dimension(nrGrid), intent(in) :: f
        integer(kind=8), intent(in) :: i
        integer(kind=8), dimension(3) ::  boundaries
        integer(kind=8) :: axis
        real(kind=8), dimension(3) :: first_finite_difference_3d

        boundaries = is_boundary(i)
        
        do axis=1,3
            select case (boundaries(axis))
                case (0) ! Inner point
                    ! print *, '  Polarso: central difference'
                    first_finite_difference_3d(axis) = central_difference(f,i,axis)
                case (1) ! Startpoint
                    ! print *, '  Polarso: forward difference'
                    first_finite_difference_3d(axis) = forward_difference(f,i,axis)
                case (2) ! Endpoint
                    ! print *, '  Polarso: backward difference'
                    first_finite_difference_3d(axis) = backward_difference(f,i,axis)
            end select
        end do

    end function first_finite_difference_3d

    function central_difference(f,x,a)
        real(kind=8), dimension(nrGrid), intent(in) :: f
        integer(kind=8), intent(in) :: x, a
        integer(kind=8) :: x_m1h, x_p1h
        real(kind=8) :: central_difference

        x_m1h = move(x,-1,a)
        x_p1h = move(x,1,a)

        central_difference = (f(x_p1h)-f(x_m1h))/(2.0d0*h_grid(a))
    end function central_difference

    function forward_difference(f,x,a)
        real(kind=8), dimension(nrGrid), intent(in) :: f
        integer(kind=8), intent(in) :: x, a
        integer(kind=8) :: x_p1h, x_p2h
        real(kind=8) :: forward_difference

        x_p1h = move(x,1,a)
        x_p2h = move(x,2,a)

        forward_difference = (-3.0d0*f(x)+4.0d0*f(x_p1h)-1.0d0*f(x_p2h))/(2.0d0*h_grid(a))

    end function forward_difference

    function backward_difference(f,x,a)
        ! Second-order backward difference, the mirror of forward_difference
        real(kind=8), dimension(nrGrid), intent(in) :: f
        integer(kind=8), intent(in) :: x, a
        integer(kind=8) :: x_m1h, x_m2h
        real(kind=8) :: backward_difference

        x_m1h = move(x,-1,a)
        x_m2h = move(x,-2,a)

        backward_difference = (3.0d0*f(x)-4.0d0*f(x_m1h)+f(x_m2h))/(2.0d0*h_grid(a))
    end function backward_difference

! -----------------------------------------------------------------------------
! ----- Boundary conditions ---------------------------------------------------
! -----------------------------------------------------------------------------

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

    subroutine full_update_charge()
        integer(kind=8) :: k, i, charge_slot
        if (allocated(newCharge_density)) then
            deallocate(newCharge_density)
        end if
        if (allocated(newCharge_index)) then
            deallocate(newCharge_index)
        end if
        allocate(newCharge_density(8,max_particles))
        allocate(newCharge_index(8,max_particles))
        newNrCharge = 0

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP& SHARED(nrElec, particles_elec_pointer, particles_mask) &
        !$OMP& PRIVATE(k, i, charge_slot)
        do k=1,nrElec
            i = particles_elec_pointer(k)
            if (particles_mask(i) .eqv. .true.) then
                call partial_update_charge(i, charge_slot)
            end if
        end do
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP& SHARED(nrIon, particles_ion_pointer, particles_mask) &
        !$OMP& PRIVATE(k, i, charge_slot)
        do k=1,nrIon
            i = particles_ion_pointer(k)
            if (particles_mask(i) .eqv. .true.) then
                call partial_update_charge(i, charge_slot)
            end if
        end do
        !$OMP END PARALLEL DO
    end subroutine full_update_charge

    subroutine partial_update_charge(k, charge_slot)
        ! Deposit the charge of particle k onto the grid (cloud-in-cell).
        ! Returns in charge_slot the column of newCharge_index /
        ! newCharge_density this charge occupies, or 0 if the particle was
        ! not deposited (outside the solver window, or no free grid nodes
        ! around it).
        integer, intent(in) :: k
        integer(kind=8), intent(out) :: charge_slot
        integer(kind=8), dimension(2,4) :: cube_id
        real(kind=8), dimension(2,4) :: cube_weight
        real(kind=8), dimension(3) :: elec_pos, node_pos
        real(kind=8) :: weight, weight_sum
        integer(kind=8) :: g, i, l, m
        real(kind=8) :: dx, dy, dz
        real(kind=8) :: wx, wy, wz

        charge_slot = 0

        elec_pos = particles_cur_pos(:,k)
        if (is_inside(elec_pos) == 0) then
            return
        end if

        i = disc_coord(elec_pos(1), elec_pos(2), elec_pos(3))
        cube_id = cube(i)

        weight_sum = 0.0d0
        do l=1,2
            do m=1,4
                g = cube_id(l,m)

                if ((is_emitter(g) == 0) .and. (is_anode(g) == 0) .and. (is_cathode(g) == 0)) then
                    node_pos = cart_coord(g)

                    dx = abs(elec_pos(1) - node_pos(1))
                    dy = abs(elec_pos(2) - node_pos(2))
                    dz = abs(elec_pos(3) - node_pos(3))
                    wx = 1.0d0 - dx / hx
                    wy = 1.0d0 - dy / hy
                    wz = 1.0d0 - dz / hz
                    weight = wx*wy*wz

                    if (weight < 0.0d0) then
                        weight = 0.0d0
                    end if

                    cube_weight(l,m) = weight
                    weight_sum = weight_sum + weight
                else
                    cube_weight(l,m) = 0.0d0
                end if
            end do
        end do

        ! All eight surrounding nodes are electrode nodes: the charge cannot
        ! be deposited (avoids division by zero below)
        if (weight_sum <= 0.0d0) then
            return
        end if

        ! full_update_charge calls this from inside an OpenMP parallel loop, so
        ! the counter must be incremented and read as one atomic operation.
        ! Otherwise two threads take the same column of newCharge_index /
        ! newCharge_density and one of the charges is lost.
        !$OMP ATOMIC CAPTURE
        newNrCharge = newNrCharge + 1
        charge_slot = newNrCharge
        !$OMP END ATOMIC

        do l=1,2
            do m=1,4
                cube_weight(l,m) = cube_weight(l,m) / weight_sum
                g = gridPoints(cube_id(l,m))
                newCharge_index((l-1)*4+m,charge_slot) = g
                ! RHS of the Poisson rows: -rho/epsilon_0 with the rows
                ! multiplied by hz**2 (see div_h2 in init_grid), and
                ! rho = weight*charge/(hx*hy*hz). Reduces to the cubic-cell
                ! form -weight*charge/(hz*epsilon_0) when hx = hy = hz.
                newCharge_density((l-1)*4+m,charge_slot) = - cube_weight(l,m)*particles_charge(k)*hz/(hx*hy*epsilon_0)
            end do
        end do

    end subroutine partial_update_charge
   
! -----------------------------------------------------------------------------
! ----- Indexing functions ----------------------------------------------------
! -----------------------------------------------------------------------------

    function disc_coord(x,y,z) ! DONE
        ! Change cartesian coordinates to discrete coordinates. The cell
        ! index is clamped to [0, N-2] per axis so that the deposition /
        ! interpolation cube (this corner plus one step in each axis, see
        ! cube()) always stays inside the grid, also for points exactly on
        ! the upper grid boundary.
        real(kind=8), intent(in) :: x, y, z
        integer(kind=8) :: x_coord,y_coord,z_coord,disc_coord

        x_coord = min(max(int((x-grid_lim(1,1)) / hx, kind=8), 0_8), Nx-2)
        y_coord = min(max(int((y-grid_lim(2,1)) / hy, kind=8), 0_8), Ny-2)
        z_coord = min(max(int((z-grid_lim(3,1)) / hz, kind=8), 0_8), Nz-2)

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

        cart_coord(1) = x_coord*hx + grid_lim(1,1)
        cart_coord(2) = y_coord*hy + grid_lim(2,1)
        cart_coord(3) = z_coord*hz + grid_lim(3,1)
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

    function cube(i)
        integer(kind=8), intent(in) :: i
        integer(kind=8), dimension(2,4) :: cube

        cube(1,1) = i
        cube(1,2) = move(i,1,1)
        cube(1,3) = move(i,1,2)
        cube(1,4) = move(move(i,1,2),1,1)

        cube(2,1) = move(i,1,3)
        cube(2,2) = move(move(i,1,1),1,3)
        cube(2,3) = move(move(i,1,2),1,3)
        cube(2,4) = move(move(move(i,1,2),1,1),1,3)

    end function cube

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
        ! Returns 1 if the point lies on the top grid plane and that plane is
        ! the physical anode (see anode_in_grid in init_grid). Index-based so
        ! it cannot miss the plane through floating-point roundoff.
        integer(kind=8), intent(in) :: i
        integer(kind=8) :: is_anode

        is_anode = 0
        if (anode_in_grid .and. ((i-1)/NxNy == Nz-1)) then
            is_anode = 1
        end if

    end function is_anode

    function is_cathode(i)
        ! Returns 1 if the point lies on the bottom grid plane and that plane
        ! is the grounded cathode plane the emitter sits on (see
        ! cathode_in_grid in init_grid). The emitter itself is handled
        ! separately by is_emitter / is_first_layer.
        integer(kind=8), intent(in) :: i
        integer(kind=8) :: is_cathode

        is_cathode = 0
        if (cathode_in_grid .and. ((i-1)/NxNy == 0)) then
            is_cathode = 1
        end if

    end function is_cathode

    function is_emitter(i)
        ! Returns 1 if point is within the emitter and 0 otherwise
        integer(kind=8), intent(in) :: i
        integer(kind=8) :: is_emitter
        real(kind=8), dimension(3) :: point_coord
        real(kind=8) :: point_radius

        is_emitter = 0

        point_coord = cart_coord(i)

        if (point_coord(3) <= emit_lim(3,2)) then
            select case (emitters_type(1))
                case (1) ! Circle
                    point_radius = sqrt(point_coord(1)**2 + point_coord(2)**2)
                    if (point_radius <= emitter_radius) then
                        is_emitter = 1
                    else
                        is_emitter = 0
                    end if

                case (2) ! Rectangle
                    if ((emit_lim(1,1) <= point_coord(1)) .and. (point_coord(1) <= emit_lim(1,2)) &
                        .and. (emit_lim(2,1) <= point_coord(2)) .and. (point_coord(2) <= emit_lim(2,2))) then
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

        if ((elec_lim(1,1) <= pos(1)) .and. (pos(1) <= elec_lim(1,2)) &
            .and. (elec_lim(2,1) <= pos(2)) .and. (pos(2) <= elec_lim(2,2)) &
            .and. (elec_lim(3,1) < pos(3)) .and. (pos(3) < elec_lim(3,2))) then
            is_inside = 1
        else
            is_inside = 0
        end if
    end function

! -----------------------------------------------------------------------------
! ----- Clean Up --------------------------------------------------------------
! -----------------------------------------------------------------------------
    
    subroutine PL_Clean_Up() ! DONE
        integer(kind=8) :: error, phase

        ! Nothing to clean up when the solver was never initialized
        if (use_polarso .eqv. .false.) then
            return
        end if

        ! Clean up memory
        print '(a)', 'Polarso: Clean up'

        ! Clear pardiso memory
        print *, '  Polarso: clearing pardiso memory'
        call pardiso_phase(-1,0)
        
        ! Deallocate CSR matrix
        deallocate(nnz_values)
        deallocate(nnz_col_index)
        deallocate(nnz_ia)
        ! Deallocate other arrays
        deallocate(b)
        deallocate(newCharge_index)
        deallocate(newCharge_density)
        deallocate(gridPoints)
        deallocate(gridPointsActive)
        deallocate(x)
        deallocate(perm)

        ! These should have been deallocated already in the last iteration
        if (allocated(polarso_field)) then
            deallocate(polarso_field) 
        end if
        if (allocated(voltage)) then
            deallocate(voltage)
        end if
    end subroutine PL_Clean_Up

    subroutine check_matrix()
        integer(kind=8) :: i,count,index,end

        do i=1,nrGridActive
            print*, 'Polarso: row=',i
            count = 0
            index = nnz_ia(i)
            end = nnz_ia(i+1)
            do while(index < end)
                print*,'    Polarso: column=',nnz_col_index(index)
                index = index + 1
                count = count + 1
            end do
        end do

    end subroutine check_matrix

    subroutine Place_Electron()
        real(kind=8), dimension(3) :: par_pos, par_vel

        par_pos = (/0.0d0, 0.0d0, 2.0d0/)*length_scale
        par_vel = 0.0d0
        call Add_Particle(par_pos,par_vel,species_elec,1,1,-1)

    end subroutine Place_Electron

    subroutine Write_Polarso_Data()
        integer :: ud_pl_field, ud_ic_field, IFAIL, lvl, i, j, k
        character(len=128) :: filename_pl_field, filename_ic_field
        real(kind=8), dimension(3) :: cur_pos, ic_field, pl_field
        real(kind=8), dimension(3,2) :: lim
        real(kind=8) :: Lx, Ly, hx_temp, hy_temp
        
        ! write(filename_lp_voltage, '(a15,i0,a4)') 'out/lp_voltage_',lvl,'.bin'
        write(filename_pl_field, '(a16)') 'out/pl_field.bin'
        write(filename_ic_field, '(a16)') 'out/ic_field.bin'

        ! ! Open the voltage file
        ! open(newunit=ud_lp_voltage, iostat=IFAIL, file=filename_lp_voltage, status='REPLACE', action='WRITE', access='STREAM')
        ! if (IFAIL /= 0) then
        !     print *, 'RUMDEED: Failed to open the Polarso voltage file.'
        !     return 
        ! end if

        ! Open the field file
        open(newunit=ud_pl_field, iostat=IFAIL, file=filename_pl_field, status='REPLACE', action='WRITE', access='STREAM')
        if (IFAIL /= 0) then
            print *, 'RUMDEED: Failed to open the polarso field file.'
            return 
        end if

        ! Open the field file
        open(newunit=ud_ic_field, iostat=IFAIL, file=filename_ic_field, status='REPLACE', action='WRITE', access='STREAM')
        if (IFAIL /= 0) then
            print *, 'RUMDEED: Failed to open the image charge field file.'
            return 
        end if

        lim = emit_lim

        Lx = lim(1,2) - lim(1,1)
        Ly = lim(2,2) - lim(2,1)
        hx_temp = Lx / 199
        hy_temp = Ly / 199

        do i=0,199
            do j=0,199
                cur_pos(1) = lim(1,1) + i*hx_temp
                cur_pos(2) = lim(2,1) + j*hy_temp
                cur_pos(3) = 0.0d0

                pl_field = PL_Calculate_Field_At(cur_pos)
                ic_field = Calc_Field_at(cur_pos)
                ! print *, 'Polarso: writing: x=',test_pos(1), ' y=',test_pos(2), ' z=',test_pos(3)
                ! write(unit=ud_lp_voltage,iostat=IFAIL) i, j, voltage(k), is_emitter(k)
                write(unit=ud_pl_field,iostat=IFAIL) i, j, pl_field(3)
                write(unit=ud_ic_field,iostat=IFAIL) i, j, ic_field(3)
            end do
        end do

        ! close(unit=ud_lp_voltage, iostat=IFAIL, status='keep')
        close(unit=ud_pl_field, iostat=IFAIL, status='keep')
        close(unit=ud_ic_field, iostat=IFAIL, status='keep')
        
    end subroutine Write_Polarso_Data

    function PL_Calculate_Field_At(point_pos)
        real(kind=8), dimension(3), intent(in) :: point_pos
        real(kind=8), dimension(3) :: PL_Calculate_Field_At
        integer(kind=8) :: i, k, l, m

        integer(kind=8), dimension(2,4) :: grid_id
        real(kind=8), dimension(2,4,3) :: grid_pos, grid_field
        real(kind=8), dimension(2,3) :: temp_field_front, temp_field_back, temp_field

        real(kind=8), dimension(3) :: point_field

        ! Points outside the solver grid (the window can be smaller than the
        ! simulation box) see only the applied vacuum field; interpolating
        ! there would index polarso_field out of bounds.
        if ((point_pos(1) < grid_lim(1,1)) .or. (point_pos(1) > grid_lim(1,2)) .or. &
            (point_pos(2) < grid_lim(2,1)) .or. (point_pos(2) > grid_lim(2,2)) .or. &
            (point_pos(3) < grid_lim(3,1)) .or. (point_pos(3) > grid_lim(3,2))) then
            PL_Calculate_Field_At = (/0.0d0, 0.0d0, -V_d/d/)
            return
        end if

        i = disc_coord(point_pos(1), point_pos(2), point_pos(3))

        grid_id = cube(i)

        do l=1,2
            do m=1,4
                k = grid_id(l,m)
                grid_pos(l,m,:) = cart_coord(k)
                grid_field(l,m,:) = polarso_field(k,:)
            end do
        end do

        ! Interpolate derivative
        do l=1,2
            temp_field_front(l,:) = trilinear_interpolate(grid_pos(l,1,1),grid_pos(l,2,1),point_pos(1),grid_field(l,1,:),grid_field(l,2,:))
            temp_field_back(l,:) = trilinear_interpolate(grid_pos(l,3,1),grid_pos(l,4,1),point_pos(1),grid_field(l,3,:),grid_field(l,4,:))
            temp_field(l,:) = trilinear_interpolate(grid_pos(l,1,2),grid_pos(l,3,2),point_pos(2),temp_field_front(l,:),temp_field_back(l,:))
        end do
        PL_Calculate_Field_At = trilinear_interpolate(grid_pos(1,1,3),grid_pos(2,1,3),point_pos(3),temp_field(1,:),temp_field(2,:))

    end function PL_Calculate_Field_At

    

    function trilinear_interpolate(x1,x2,xi,f1,f2)
        real(kind=8), intent(in) :: x1,x2,xi
        real(kind=8), dimension(3), intent(in) :: f1, f2
        real(kind=8), dimension(3) :: trilinear_interpolate
        integer(kind=8) :: i

        do i=1,3
            trilinear_interpolate(i) = f1(i) + (f2(i)-f1(i))/(x2-x1+min_d)*(xi-x1)
        end do

    end function trilinear_interpolate

    function lagrange_interpolation(x,x1,x2,x3,x4,y1,y2,y3,y4)
        real(kind=8), intent(in) :: x,x1,x2,x3,x4,y1,y2,y3,y4
        real(kind=8) :: lagrange_interpolation, L1, L2, L3, L4
        integer(kind=8) :: i

        lagrange_interpolation = 0.0d0

        L1 = ((x-x2)*(x-x3)*(x-x4))/((x1-x2)*(x1-x3)*(x1-x4))
        L2 = ((x-x1)*(x-x3)*(x-x4))/((x2-x1)*(x2-x3)*(x2-x4))
        L3 = ((x-x1)*(x-x2)*(x-x4))/((x3-x1)*(x3-x2)*(x3-x4))
        L4 = ((x-x1)*(x-x2)*(x-x3))/((x4-x1)*(x4-x2)*(x4-x3))

        lagrange_interpolation = L1*y1 + L2*y2 + L3*y3 + L4*y4
    end function lagrange_interpolation

    function sparse_dot(n,a,ia,ja,v)
        integer(kind=8), intent(in) :: n
        real(kind=8), dimension(:), intent(in) :: a,v
        integer(kind=8), dimension(:), intent(in) :: ia,ja
        real(kind=8), allocatable, dimension(:) :: sparse_dot
        integer(kind=8) :: i, j, col_start, col_end

        allocate(sparse_dot(n))

        ! Each row is written by exactly one thread, so a plain shared write
        ! is correct (a reduction over the whole result vector is not needed)
        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(n,a,ia,ja,v,sparse_dot) &
        !$OMP PRIVATE(i,j,col_start,col_end)
        do i=1,n
            col_start = ia(i)
            col_end = ia(i+1)-1

            sparse_dot(i) = 0.0d0
            do j=col_start,col_end
                sparse_dot(i) = sparse_dot(i) + a(j)*v(ja(j))
            end do
        end do
        !$OMP END PARALLEL DO

    end function sparse_dot
end module mod_polarso

#else

! ------------------------------------------------------------------------------
! Stub version of mod_polarso, used when RUMDEED is built without the POLARSO
! Laplace solver (the default). It has the same public interface as the real
! module but no MKL dependency. Read_Polarso_Variables stops the program if a
! 'polarso' input file is present, so a run that expects the solver cannot
! silently proceed without it; the other stubs are unreachable because
! use_polarso stays .false. (they stop the program if called anyway).
! ------------------------------------------------------------------------------
module mod_polarso

    use mod_global, only: use_polarso

    implicit none

    private
    public :: PL_Init_Solver, PL_Calculate_Field, PL_Update_Field, PL_Calculate_Field_At, &
                PL_Clean_Up, Place_Electron, Write_Polarso_Data, Read_Polarso_Variables

contains

    subroutine Read_Polarso_Variables()
        logical :: polarso_file_exists

        inquire(file='polarso', exist=polarso_file_exists)
        if (polarso_file_exists .eqv. .true.) then
            print '(a)', 'RUMDEED: ERROR a polarso input file is present, but this binary was'
            print '(a)', 'RUMDEED: built without the POLARSO solver. Rebuild with make POLARSO=yes.'
            stop
        end if
    end subroutine Read_Polarso_Variables

    subroutine Polarso_Not_Available()
        print '(a)', 'RUMDEED: ERROR POLARSO was called but this binary was built without it.'
        print '(a)', 'RUMDEED: Rebuild with make POLARSO=yes.'
        stop
    end subroutine Polarso_Not_Available

    subroutine PL_Init_Solver()
        call Polarso_Not_Available()
    end subroutine PL_Init_Solver

    subroutine PL_Calculate_Field()
        call Polarso_Not_Available()
    end subroutine PL_Calculate_Field

    subroutine PL_Update_Field(index)
        integer, intent(in) :: index

        call Polarso_Not_Available()
    end subroutine PL_Update_Field

    function PL_Calculate_Field_At(point_pos)
        real(kind=8), dimension(3), intent(in) :: point_pos
        real(kind=8), dimension(3) :: PL_Calculate_Field_At

        PL_Calculate_Field_At = 0.0d0
        call Polarso_Not_Available()
    end function PL_Calculate_Field_At

    subroutine PL_Clean_Up()
        ! Nothing to clean up, the solver was never initialized
    end subroutine PL_Clean_Up

    subroutine Place_Electron()
        call Polarso_Not_Available()
    end subroutine Place_Electron

    subroutine Write_Polarso_Data()
        call Polarso_Not_Available()
    end subroutine Write_Polarso_Data

end module mod_polarso

#endif