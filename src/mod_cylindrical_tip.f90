!--------------------------------------------!
! Module for emission from a cylindrical tip !
! Kristinn Torfason                          !
! 14.09.23                                   !
!--------------------------------------------!

Module mod_cylindrical_tip
    use mod_global
    use mod_verlet
    use mod_pair
    use mod_work_function
    use kdtree2_precision_module
    use kdtree2_module
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan

    implicit none

    PRIVATE
    PUBLIC :: Init_Cylindrical_Tip, Clean_Up_Cylindrical_Tip

    ! ----------------------------------------------------------------------------
    ! KD-tree variables
    type(kdtree2), pointer :: kd_tree => null() ! KD-tree for the mesh
    double precision, dimension(:, :), allocatable :: kd_mesh ! Mesh for the KD-tree
    double precision, dimension(:, :), allocatable :: kd_data ! Data for the KD-tree

    ! ----------------------------------------------------------------------------
    ! Variables
    integer, dimension(:), allocatable :: nrEmitted_emitters

    ! ----------------------------------------------------------------------------
    ! Parameters for the cylindrical tip
    double precision, parameter :: length_scale_cyl = 1.0d-6 ! Length scale for the cylindrical tip
    double precision, parameter :: radius_cyl = 12.0d0*length_scale_cyl ! Radius of the cylinder
    double precision, parameter :: height_cyl = 100.0d0*length_scale_cyl ! Height of the cylinder
    double precision, parameter :: radius_cor = 1.0d0*length_scale_cyl ! Radius of the corner

    integer, parameter          :: sec_top = 1 ! Top section
    integer, parameter          :: sec_side = 2 ! Side section
    integer, parameter          :: sec_corner = 3 ! Corner section

    ! ----------------------------------------------------------------------------
    ! Area of each section of the cylindrical tip
    double precision, parameter :: area_top = pi * (radius_cyl - radius_cor)**2
    double precision, parameter :: area_side = 2.0d0 * pi * radius_cyl * (height_cyl - radius_cor)
    double precision, parameter :: area_corner = pi**2 * radius_cor * (2.0d0 * radius_cor/pi + radius_cyl - radius_cor)
    double precision, parameter :: area_total = area_top + area_side + area_corner

    ! ----------------------------------------------------------------------------
    ! Probability of emitting from each section
    double precision            :: prob_top = area_top / area_total
    double precision            :: prob_side = area_side / area_total
    double precision            :: prob_corner = area_corner / area_total

    ! ----------------------------------------------------------------------------
    ! Other variables
    double precision :: residual = 0.0d0 ! Residual from the previous time step
    ! Constant used in MC integration (function Elec_Supply_V2)
    double precision :: time_step_div_q0
    ! Average field
    double precision, dimension(1:3)   :: F_avg = 0.0d0
    integer                            :: ud_cyl_debug ! Unit for debugging the cylindrical tip

    ! ----------------------------------------------------------------------------
    ! Constants for field emission
    ! First Fowler-Nordheim constant in units [ A eV V^{-2} ]
    double precision, parameter :: a_FN = q_02/(16.0d0*pi**2*h_bar)

    ! Second Fowler-Nodheim constant in units [ eV^{-3/2} V m^{-1} ]
    double precision, parameter :: b_FN = -4.0d0/(3.0d0*h_bar) * sqrt(2.0d0*m_0*q_0)
  
    ! Constant used for calculation of the l in v_y and t_y.
    ! The units are [ eV^{2} V^{-1} m ]
    ! See Forbes, R. G., & Deane, J. H. (2007, November).
    ! "Reformulation of the standard theory of Fowler–Nordheim tunnelling and cold field electron emission."
    ! In Proceedings of the Royal Society of London A: Mathematical,
    ! Physical and Engineering Sciences (Vol. 463, No. 2087, pp. 2907-2927). The Royal Society.
    double precision, parameter :: l_const = q_0 / (4.0d0*pi*epsilon_0)

    ! ----------------------------------------------------------------------------
    ! Variables for the Metropolis-Hastings algorithm
    integer, parameter                 :: N_MH_step = 25 ! Number of steps to do in the MH algorithm

contains

!-------------------------------------------!
! Initialize the cylinder emission
! TODO: Change this for the new system
subroutine Init_Cylindrical_Tip()
    ! Local variables
    double precision :: int_res
    double precision, dimension(1:3) :: E_test, pos_test
    integer :: IFAIL

    ! Allocate the number of emitters
    allocate(nrEmitted_emitters(1:nrEmit))

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_Planar ! Check_Boundary_Cylinder

    ! Function for the electric field in the system
    ptr_field_E => field_E_cylinder

    ! The function that does the emission
    ptr_Do_Emission => Do_Field_Emission_Cylinder

    ! The function to do image charge effects
    ptr_Image_Charge_effect => Force_Image_charges_v2 ! Force_Image_charges for the cylinder

    call Read_work_function()

    ! Parameters used in the module
    time_step_div_q0 = time_step / q_0

    ! Create the KD-tree
    call Create_KD_Tree()

    ! ! Open the file for debugging
    ! open(newunit=ud_cyl_debug, file='cyl_debug.dt', status='replace', action='write')

    ! ! Debugging the kd-tree, its data and interpolation
    ! pos_test(1) = -0.011920996665774017d0*1.0d-3
    ! pos_test(2) = -2.6827637010898934d-5*1.0d-3
    ! pos_test(3) =  0.0993894990837496d0*1.0d-3
    ! E_test = field_E_cylinder(pos_test)
    ! print *, 'pos_test = ', pos_test
    ! print *, 'E_test = ', E_test
    ! print *, 'E_test = ', sqrt(E_test(1)**2 + E_test(2)**2 + E_test(3)**2)
    ! ! Results should be close to 285521.7478622014 V/m
    ! print *, '? 285521.7478622014 V/m ?'
    ! print *, ''

    ! print *, 'Calling Calc_E_edge_cyl'
    ! call Calc_E_edge_cyl()
    ! print *, 'Calling Calc_E_side_cyl'
    ! call Calc_E_top_cyl()
    ! print *, 'Calling Calc_E_top_cyl'
    ! call Calc_E_corner_cyl()

    ! ! Debugging the integration
    ! print *, 'Calling Do_Cuba_Suave'
    ! call Do_Cuba_Suave(int_res)
    ! !print *, 'Calling Do_Cuba_Divonne'
    ! !call Do_Cuba_Divonne(int_res)

    ! !print *, 'RUMDEED: Integration result = ', int_res
    ! !print *, 'Top area = ', area_top
    ! !print *, 'Side area = ', area_side
    ! !print *, 'Corner area = ', area_corner

    ! close(unit=ud_cyl_debug, iostat=IFAIL, status='keep')
    ! stop
end subroutine Init_Cylindrical_Tip

!-------------------------------------------!
! Create the K-dimensional tree using the kdtree2 module
subroutine Create_KD_Tree()
    character (len=*), parameter :: filename_meshdata="Cyl_mesh_data.txt"
    integer                      :: IFAIL
    integer                      :: ud_meshdata
    integer                      :: n_points, k

    ! Read the mesh and data from file
    ! Open the file for reading
    open(newunit=ud_meshdata, iostat=IFAIL, file=filename_meshdata, status='OLD', action='READ')
    if (IFAIL /= 0) then
      print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file ', filename_meshdata
      stop
    end if

    ! Count the number of lines in the file
    n_points = 0
    do
      read(ud_meshdata, *, iostat=IFAIL)
      if (IFAIL /= 0) exit
      n_points = n_points + 1
    end do

    ! Allocate the array to hold the file
    allocate(kd_mesh(1:3, 1:n_points))
    allocate(kd_data(1:3, 1:n_points))

    ! Rewind to the start of the file
    rewind(ud_meshdata)

    ! Read the file into the array and scale the coordinates from mm to m
    do k = 1, n_points
        read(ud_meshdata, *) kd_mesh(1, k), kd_mesh(2, k), kd_mesh(3, k), kd_data(1, k), kd_data(2, k), kd_data(3, k)
        ! Convert from mm to m
        kd_mesh(:, k) = kd_mesh(:, k) * 1.0d-3
    end do

    ! Close the file
    close(unit=ud_meshdata, iostat=IFAIL, status='keep')

    ! Create the kd-tree
    kd_tree => kdtree2_create(kd_mesh, sort=.true., rearrange=.true.)
end subroutine Create_KD_Tree

function field_E_cylinder(pos) result(field_E)
    ! Input variables
    double precision, dimension(1:3), intent(in) :: pos ! Position to calculate the electric field at

    ! Output variables
    double precision, dimension(1:3)             :: field_E ! Electric field at the point

    ! Local variables
    integer, parameter                        :: nn = 4 ! number of nearest neighbors to find
    integer                                   :: k ! loop variable
    type(kdtree2_result)                      :: kd_results(1:nn) ! results from the kd-tree
    double precision, dimension(1:nn)         :: w ! weights for the interpolation
    double precision, parameter               :: eps = length_scale**3 ! Small number to avoid division by zero
    integer, parameter                        :: p = 1 ! Power for the inverse distance weighting

    !! For debugging start with planar field
    !field_E = field_E_planar(pos)

    ! Find the nearest neighbors
    !print *, 'pos = ', pos
    ! Note: The kd-tree is not thread safe, so we have to use a critical section
    !$omp critical (kdtree2)
    call kdtree2_n_nearest(tp=kd_tree, qv=pos, nn=nn, results=kd_results)
    !$omp end critical (kdtree2)

    ! Interpolate the electric field at the point using the nearest neighbors and the inverse distance to each as weights
    ! Calculate the weights as 1/distance^p
    !print *, 'kd_results(1:nn)%idx = ', kd_results(1:nn)%idx
    !print *, 'kd_results(1:nn)%dis = ', kd_results(1:nn)%dis
    w(1:nn) = 1.0d0/(kd_results(1:nn)%dis**p + eps)

    ! Normalize the weights
    w = w / sum(w)
    
    ! Calculate the interpolated electric field
    field_E = 0.0d0
    do k = 1, nn
        field_E = field_E + w(k) * kd_data(:, kd_results(k)%idx)
    end do
    field_E = field_E * 0.25d2 ! Debug to make field larger
end function field_E_cylinder

subroutine Do_Field_Emission_Cylinder(step)
    integer, intent(in) :: step
    integer             :: i, IFAIL
    integer, parameter  :: emit = 1

    ! Integration
    double precision                 :: N_sup
    integer                          :: N_round

    ! Emission variables
    double precision                 :: D_f, Df_avg, F_norm, rnd
    double precision, dimension(1:3) :: F, n_vec
    integer                          :: s, sec
    integer                          :: nrElecEmit
    double precision, dimension(1:3) :: par_pos, par_vel

    nrEmitted_emitters = 0

    call Do_Surface_Integration(N_sup)
    !N_sup = 1.0d0 ! Start with one electron per time step
    N_round = nint(N_sup + residual)
    residual = N_sup - N_round

    nrElecEmit = 0
    nrEmitted_emitters(emit) = 0

    ! Loop over the electrons to be emitted.
    do s = 1, N_round
        ! Get the position of the particle
        ! ndim_in, F_out, F_norm_out, pos_xyz_out, sec_out
        call Metropolis_Hastings_cyl_tip(N_MH_step, F, F_norm, par_pos, sec, n_vec)
  
        ! Check if the field is favorable for emission or not
        if (F_norm >= 0.0d0) then
          D_f = -huge(1.0d0)
          print *, 'Warning: F > 0.0d0 (Do_Field_Emission_Cylinder)'
          print *, 'F_norm = ', F_norm
          print *, 'F = ', F
          print *, 'par_pos = ', par_pos
          print *, 'sec = ', sec
          stop
        else
          D_f = Escape_Prob_log(F_norm, par_pos)
        end if
        df_avg = df_avg + exp(D_f)
  
        CALL RANDOM_NUMBER(rnd)
        if (log(rnd) <= D_f) then
          ! Calculate the position of the particle (1 nm above the surface)
          par_pos = par_pos + n_vec*length_scale
          ! Calculate the velocity of the particle
          par_vel = 0.0d0

          ! Add the particle to the system
          call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1, sec)

          nrElecEmit = nrElecEmit + 1
          nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
        end if
      end do

    write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, nrElecEmit, nrElec, &
                                                       & (nrEmitted_emitters(i), i = 1, nrEmit)
end subroutine Do_Field_Emission_Cylinder

subroutine Metropolis_Hastings_cyl_tip(ndim_in, F_out, F_norm_out, pos_xyz_out, sec_out, n_vec)
    ! Input variables
    integer, intent(in)                           :: ndim_in
    ! Output variables
    double precision, intent(out)                 :: F_norm_out
    double precision, intent(out), dimension(1:3) :: F_out, n_vec
    double precision, intent(out), dimension(1:3) :: pos_xyz_out
    integer, intent(out)                          :: sec_out

    ! Local variables
    integer                          :: IFAIL, k
    double precision, dimension(1:3) :: pos_cur, pos_new, pos_xyz_cur, pos_xyz_new
    double precision, dimension(1:3) :: F_cur, F_new
    double precision                 :: F_norm_cur, F_norm_new
    integer                          :: sec_cur, sec_new
    double precision                 :: alpha, rnd
    !double precision                 :: ratio
    integer                          :: i, accepted = 0, rejected = 0

    k = 0
    ! Get a random initial position on the surface
    do
      k = k + 1
      call Get_Random_Surface_Pos(pos_cur, sec_cur)

      F_cur = Calc_Field_at(pos_cur)
      pos_xyz_cur = Convert_to_xyz(pos_cur, sec_cur)
      F_norm_cur = Field_normal(pos_xyz_cur, pos_cur, F_cur, sec_cur)
      if (F_norm_cur < 0.0d0) then ! Check that the field is negative
        exit ! Exit the loop
      end if

      if (k > 100000) then
        print *, 'RUMDEED: ERROR Cannot find a favorable starting position (Metropolis_Hastings_cyl_tip)'
        stop
      end if
    end do

    ! Do the Metropolis-Hastings algorithm
    do i = 1, ndim_in ! Loop over all the steps

        ! Propose a new position
        call Jump_MH(pos_cur, sec_cur, pos_new, sec_new)
        pos_xyz_new = Convert_to_xyz(pos_new, sec_new)

        ! Find the jump probability at the current position
        alpha = Get_Jump_Probability(F_norm_cur, F_new, F_norm_new, pos_new, pos_xyz_new, sec_new)

        ! Accept or reject the jump
        call random_number(rnd)
        if (rnd < alpha) then
            ! Accept the jump
            accepted = accepted + 1

            ! Update the position
            pos_cur = pos_new
            sec_cur = sec_new
            pos_xyz_cur = pos_xyz_new

            ! Update the electric field
            F_cur = F_new
            F_norm_cur = F_norm_new
            
            write (unit=ud_mh, iostat=IFAIL) cur_time, pos_xyz_cur/length_scale, sec_cur, pos_cur
        else
            ! Reject the jump
            rejected = rejected + 1
        end if
    end do

    ! Calculate the acceptance ratio
    !ratio = dble(accepted) / dble(accepted + rejected) * 100.d0
    !print *, 'RUMDEED: Metropolis-Hastings acceptance ratio = ', ratio, '%'

    ! Set the output variables
    F_out = F_cur
    F_norm_out = F_norm_cur
    pos_xyz_out = pos_xyz_cur
    sec_out = sec_cur
    n_vec = surface_normal(pos_xyz_cur, pos_cur, sec_cur)
end subroutine Metropolis_Hastings_cyl_tip

subroutine Jump_MH(pos_cur, sec_cur, pos_new, sec_new)
    ! Input variables
    double precision, dimension(1:3), intent(in) :: pos_cur
    integer, intent(in)                          :: sec_cur

    ! Output variables
    double precision, dimension(1:3), intent(out) :: pos_new
    integer, intent(out)                          :: sec_new

    sec_new = sec_cur ! Set the section to be the same as the current section

    ! Check which section the point is in
    if (sec_cur == sec_top) then
        call Jump_MH_Top(pos_cur, pos_new)
        call Check_Boundary_Top(pos_new, sec_new)
    else if (sec_cur == sec_side) then
        call Jump_MH_Side(pos_cur, pos_new)
        call Check_Boundary_Side(pos_new, sec_new)
    else if (sec_cur == sec_corner) then
        call Jump_MH_Corner(pos_cur, pos_new)
        call Check_Boundary_Corner(pos_new, sec_new)
    else
        print *, 'RUMDEED: Error Unknown section (Jump_MH)'
        print *, 'sec_cur = ', sec_cur
        stop
    end if
end subroutine Jump_MH

! Jump from the current position to a new position on the top surface
! Coordinates are cartesian (x, y, z)
subroutine Jump_MH_Top(pos_cur, pos_new)
    ! Input variables
    double precision, dimension(1:3), intent(in) :: pos_cur
    
    ! Output variables
    double precision, dimension(1:3), intent(out) :: pos_new

    ! Local variables
    double precision, dimension(1:2) :: std_xy

    ! Jump in x and y
    std_xy(:) = (radius_cyl - radius_cor) * 0.05d0 ! Standard deviation for the jump in x and y
    pos_new(1:2) = box_muller(pos_cur(1:2), std_xy) ! Jump in x and y
    pos_new(3) = height_cyl ! z is fixed to be at the top of the cylinder
end subroutine Jump_MH_Top

! Check the boundary conditions for the top surface
! Boundary is a circle with radius (radius_cyl - radius_cor)
! x^2 + y^2 should be less than (radius_cyl - radius_cor)^2
! If x^2 + y^2 is more than (radius_cyl - radius_cor)^2, then the point is moved to the corner
! z should be at the top of the cylinder
subroutine Check_Boundary_Top(pos_new, sec_new)
    ! Input/output variables
    double precision, dimension(1:3), intent(inout) :: pos_new
    integer, intent(inout)                          :: sec_new
    
    ! Local variables
    double precision                 :: r2
    double precision                 :: phi_new, theta_new
    !double precision                 :: dx, dy
    double precision, dimension(1:3) :: p_c, p_vec
    !double precision, dimension(1:3) :: R_unit

    ! Check if the point is inside the boundary
    r2 = pos_new(1)**2 + pos_new(2)**2
    ! Reflect the point back to the top surface
    ! if (r2 > (radius_cyl - radius_cor)**2) then
    !     R_unit = pos_new / sqrt(r2)

    !     dx = (radius_cyl - radius_cor) * R_unit(1) - pos_new(1)
    !     dy = (radius_cyl - radius_cor) * R_unit(2) - pos_new(2)

    !     pos_new(1) = pos_new(1) + 2.0d0 * dx
    !     pos_new(2) = pos_new(2) + 2.0d0 * dy
    ! end if
    if (r2 > (radius_cyl - radius_cor)**2) then
        ! Move the point to the corner
        sec_new = sec_corner

        ! Convert to toroidal coordinates and move the point to the surface of the torus
        ! Calculate the phi coordinate of the point
        phi_new = atan2(pos_new(2), pos_new(1))

        ! Calculate the theta coordinate of the point
        ! Find the center point inside the torus closest to p
        p_c(1) = (radius_cyl - radius_cor) * cos(phi_new)
        p_c(2) = (radius_cyl - radius_cor) * sin(phi_new)
        p_c(3) = height_cyl - radius_cor

        ! Vector from the center inside of the torus to the point
        p_vec = pos_new - p_c
        
        ! Find the angle this vector makes with the xy-plane
        ! This is the same theta angle as in the toroidal coordinates
        theta_new = atan2(p_vec(3), sqrt(p_vec(1)**2 + p_vec(2)**2))

        pos_new(1) = radius_cor
        pos_new(2) = phi_new
        pos_new(3) = theta_new
    else ! The point is inside the boundary
        ! Keep the point on the top surface
        sec_new = sec_top
    end if
end subroutine Check_Boundary_Top

! Jump from the current position to a new position on the side surface
! Coordinates are cylindrical (rho, phi, z)
subroutine Jump_MH_Side(pos_cur, pos_new)
    ! Input variables
    double precision, dimension(1:3), intent(in) :: pos_cur
    
    ! Output variables
    double precision, dimension(1:3), intent(out) :: pos_new

    ! Local variables
    double precision, dimension(1:2) :: std_phi_z

    ! Jump in phi and z
    std_phi_z(1) = 2.0d0*pi * 0.005d0 ! Standard deviation for the jump in phi
    std_phi_z(2) = (height_cyl - radius_cor) * 0.05d0 ! Standard deviation for the jump in z

    pos_new(2:3) = box_muller(pos_cur(2:3), std_phi_z) ! Jump in phi and z
    pos_new(1)   = radius_cyl ! rho is fixed to be radius_cyl
end subroutine Jump_MH_Side

! Check the boundary conditions for the side surface
! Boundary is a cylinder with radius radius_cyl and height height_cyl
! phi should be between 0 and 2*pi
! z should be between 0 and (height_cyl - radius_cor)
! If z is more than (height_cyl - radius_cor), then the point is moved to the corner.
! If z is less than 0, then the point if reflected back to the side surface.
! rho should equal radius_cyl
subroutine Check_Boundary_Side(pos_new, sec_new)
    ! Input/output variables
    double precision, dimension(1:3), intent(inout) :: pos_new
    integer, intent(inout)                          :: sec_new

    ! Local variables
    double precision                 :: phi_new, theta_new
    double precision, dimension(1:3) :: p_c, p_vec

    if (pos_new(3) > (height_cyl - radius_cor)) then
        ! Move the point to the corner
        sec_new = sec_corner

        ! Convert to toroidal coordinates and move the point to the surface of the torus
        ! Calculate the phi coordinate of the point
        phi_new = atan2(pos_new(2), pos_new(1))

        ! Calculate the theta coordinate of the point
        ! Find the center point inside the torus closest to p
        p_c(1) = (radius_cyl - radius_cor) * cos(phi_new)
        p_c(2) = (radius_cyl - radius_cor) * sin(phi_new)
        p_c(3) = height_cyl - radius_cor

        ! Vector from the center inside of the torus to the point
        p_vec = pos_new - p_c
        
        ! Find the angle this vector makes with the xy-plane
        ! This is the same theta angle as in the toroidal coordinates
        theta_new = atan2(p_vec(3), sqrt(p_vec(1)**2 + p_vec(2)**2))

        pos_new(1) = radius_cor
        pos_new(2) = phi_new
        pos_new(3) = theta_new

        ! ! Reflect the point back to the side surface
        ! pos_new(3) = 2.0d0 * (height_cyl - radius_cor) - pos_new(3)
    else if (pos_new(3) < 0.0d0) then
        ! Reflect the point back to the side surface
        pos_new(3) = -pos_new(3)

        ! Check that phi is between 0 and 2*pi
        pos_new(2) = modulo(pos_new(2), 2.0d0*pi)

        ! Keep the point on the side surface
        sec_new = sec_side
    else
        ! Keep the point on the side surface
        sec_new = sec_side

        ! Check that phi is between 0 and 2*pi
        pos_new(2) = modulo(pos_new(2), 2.0d0*pi)
    end if
end subroutine Check_Boundary_Side

! Jump from the current position to a new position on the corner surface
! Coordinates are toroidal (rho, phi, theta)
subroutine Jump_MH_Corner(pos_cur, pos_new)
    ! Input variables
    double precision, dimension(1:3), intent(in) :: pos_cur
    
    ! Output variables
    double precision, dimension(1:3), intent(out) :: pos_new

    ! Local variables
    double precision, dimension(1:2) :: std_phi_theta

    ! Jump in phi and theta
    std_phi_theta(1) = 2.0d0*pi * 0.005d0 ! Standard deviation for the jump in phi
    std_phi_theta(2) = pi/2.0d0 * 0.005d0 ! Standard deviation for the jump in theta

    pos_new(2:3) = box_muller(pos_cur(2:3), std_phi_theta) ! Jump in phi and theta
    pos_new(1)   = radius_cor ! rho is fixed to be radius_cor
end subroutine Jump_MH_Corner

! Check the boundary conditions for the corner surface
! Boundary is a torus with major radius (radius_cyl - radius_cor) and minor radius radius_cor
! phi should be between 0 and 2*pi
! theta should be between 0 and pi/2
! If theta is more than pi/2, then the point is moved to the top surface.
! If theta is less then 0, then the point is moved to the side surface.
! rho should equal radius_cor
subroutine Check_Boundary_Corner(pos_new, sec_new)
    ! Input/output variables
    double precision, dimension(1:3), intent(inout) :: pos_new
    integer, intent(inout)                          :: sec_new

    ! Local variables
    double precision                 :: z_new

    if (pos_new(3) > pi/2.0d0) then
        ! Convert to cartesian coordinates and move the point to the surface of the cylinder
        pos_new = Convert_to_xyz(pos_new, sec_new)
        pos_new(3) = height_cyl

        ! Move the point to the top surface
        sec_new = sec_top
    else if (pos_new(3) < 0.0d0) then
        ! Move the point to the side surface
        sec_new = sec_side

        ! Convert to cylindrical coordinates and move the point to the surface of the cylinder
        z_new = (height_cyl - radius_cor) + radius_cor * sin(pos_new(3))

        pos_new(1) = radius_cor
        ! Same phi, check that it is between 0 and 2*pi
        pos_new(2) = modulo(pos_new(2), 2.0d0*pi)
        pos_new(3) = z_new
    else
        ! Keep the point on the corner surface
        sec_new = sec_corner

        ! rho is fixed to be radius_cor
        pos_new(1) = radius_cor

        ! Check that phi is between 0 and 2*pi
        pos_new(2) = modulo(pos_new(2), 2.0d0*pi)

        ! Check that theta is between 0 and pi/2
        pos_new(3) = modulo(pos_new(3), pi/2.0d0)
    end if
end subroutine Check_Boundary_Corner

!-------------------------------------------!
! Calculate the jump probability for the Metropolis-Hastings algorithm
function Get_Jump_Probability(F_norm_cur, F_new, F_norm_new, pos_new, pos_xyz_new, sec_new) result(alpha)
    ! Input variables
    double precision, intent(in)                  :: F_norm_cur ! Component of the electric field normal to the surface at the current position
    double precision, dimension(1:3), intent(in)  :: pos_xyz_new, pos_new ! Position of the particle in cartesian coordinates and section coordinates
    integer, intent(in)                           :: sec_new ! Section the particle is in

    ! Output variables
    double precision, intent(out)                 :: F_norm_new ! Component of the electric field normal to the surface at the new position
    double precision, dimension(1:3), intent(out) :: F_new ! Electric field at the new position
    double precision                              :: alpha ! Jump probability

    ! Calculate the electric field at the new position
    F_new = Calc_Field_at(pos_xyz_new)
    F_norm_new = Field_normal(pos_xyz_new, pos_new, F_new, sec_new)

    ! Calculate the jump probability
    if (F_norm_new >= 0.0d0) then ! If the field is positive, then the jump probability is zero
      alpha = 0.0d0
    else
      alpha = F_norm_new / F_norm_cur
    end if
end function Get_Jump_Probability

! -------------------------------------------!
! Get a random position on the top surface
subroutine Get_Random_Surface_Pos(pos_out, sec_out)
    ! Input variables
    double precision, intent(out), dimension(1:3) :: pos_out
    integer, intent(out)                          :: sec_out

    ! Local variables
    double precision :: rnd

    ! Get a random number
    call random_number(rnd)
    !rnd = 0.0d0 ! Always top section for debugging
    !rnd = prob_top + 0.01d0 ! Always side section for debugging
    !rnd = prob_top + prob_side + 0.01d0 ! Always corner section for debugging

    ! Check which section to emit from
    if (rnd < prob_top) then
        ! Top surface
        sec_out = sec_top
        pos_out = Get_Top_Surface_Pos()
    else if (rnd < prob_top + prob_side) then
        ! Side surface
        sec_out = sec_side
        pos_out = Get_Side_Surface_Pos()
    else
        ! Corner surface
        sec_out = sec_corner
        pos_out = Get_Corner_Surface_Pos()
    end if
end subroutine Get_Random_Surface_Pos

! -------------------------------------------!
! Get a random position on the top surface
! Coordinates are cartesian (x, y, z)
function Get_Top_Surface_Pos() result(pos_out)
    double precision, dimension(1:3) :: pos_out

    ! Generate a random position on the top surface
    ! x and y are random numbers between -(radius_cyl - radius_cor) and (radius_cyl - radius_cor)
    call random_number(pos_out(1:2))
    pos_out(1:2) = (2.0d0 * pos_out(1:2) - 1.0d0) * (radius_cyl - radius_cor)

    ! They have to satisfy x^2 + y^2 < (radius_cyl - radius_cor)^2
    do while (pos_out(1)**2 + pos_out(2)**2 > (radius_cyl - radius_cor)**2)
        call random_number(pos_out(1:2))
        pos_out(1:2) = (2.0d0 * pos_out(1:2) - 1.0d0) * (radius_cyl - radius_cor)
    end do

    ! z is fixed to be at the top of the cylinder
    pos_out(3) = height_cyl
end function Get_Top_Surface_Pos

! -------------------------------------------!
! Get a random position on the side surface
! Coordinates are cylindrical (rho, phi, z)
function Get_Side_Surface_Pos() result(pos_out)
    double precision, dimension(1:3) :: pos_out

    ! Generate a random position on the side surface
    call random_number(pos_out(2:3))

    ! phi is a random number between 0 and 2*pi
    pos_out(2) = pos_out(2) * 2.0d0 * pi

    ! z is a random number between 0 and (height_cyl - radius_cor)
    pos_out(3) = pos_out(3) * (height_cyl - radius_cor)

    ! rho is fixed to be radius_cyl
    pos_out(1) = radius_cyl
end function Get_Side_Surface_Pos

! -------------------------------------------!
! Get a random position on the corner surface
! Coordinates are toroidal (rho, phi, theta)
function Get_Corner_Surface_Pos() result(pos_out)
    double precision, dimension(1:3) :: pos_out

    ! Generate a random position on the corner surface
    call random_number(pos_out(2:3))

    ! phi is a random number between 0 and 2*pi
    pos_out(2) = pos_out(2) * 2.0d0 * pi

    ! theta is a random number between 0 and pi/2
    pos_out(3) = pos_out(3) * pi / 2.0d0

    ! rho is fixed to be radius_cor
    pos_out(1) = radius_cor
end function Get_Corner_Surface_Pos

! -------------------------------------------!
! Convert to cartesian coordinates
pure function Convert_to_xyz(pos_in, sec_in) result(pos_out)
    ! Input variables
    double precision, dimension(1:3), intent(in) :: pos_in
    integer, intent(in)                          :: sec_in

    ! Output variables
    double precision, dimension(1:3) :: pos_out

    ! Check which section the point is in
    if (sec_in == sec_top) then ! Top section
        ! Coordinates are cartesian (x, y, z)
        pos_out = pos_in

    else if (sec_in == sec_side) then ! Side section
        ! Coordinates are cylindrical (rho, phi, z)
        pos_out(1) = radius_cyl * cos(pos_in(2))
        pos_out(2) = radius_cyl * sin(pos_in(2))
        pos_out(3) = pos_in(3)

    else if (sec_in == sec_corner) then ! Corner section
        ! Coordinates are toroidal (rho, phi, theta)
        pos_out(1) = (radius_cyl - radius_cor + radius_cor*cos(pos_in(3))) * cos(pos_in(2))
        pos_out(2) = (radius_cyl - radius_cor + radius_cor*cos(pos_in(3))) * sin(pos_in(2))
        pos_out(3) = (height_cyl - radius_cor) + radius_cor * sin(pos_in(3))

    else ! Error
        !print *, 'Error: Unknown section (Convert_to_xyz)'
        !stop
        pos_out = ieee_value( pos_out, ieee_signaling_nan )
    end if
end function Convert_to_xyz

!-------------------------------------------!
! Clean up the cylinder emission
subroutine Clean_Up_Cylindrical_Tip()
    deallocate(nrEmitted_emitters)
end subroutine Clean_Up_Cylindrical_Tip

!-------------------------------------------!
! Surface integration

  ! ----------------------------------------------------------------------------
  ! This function is called to do the surface integration.
  ! Here we select the method to do it.
  !
subroutine Do_Surface_Integration(N_sup)
    double precision, intent(out) :: N_sup ! Number of electrons

    SELECT CASE (cuba_method)
    case(cuba_method_suave)
      call Do_Cuba_Suave(N_sup)
    case(cuba_method_divonne)
      call Do_Cuba_Divonne(N_sup)
    case default
      print '(a)', 'RUMDEED: ERROR UNKNOWN INTEGRATION METHOD'
      print *, cuba_method
      stop
    end select
  end subroutine Do_Surface_Integration

!----------------------------------------------------------------------------------------
! The functions v_y and t_y are because of the image charge effect in the FN equation.
! The approximation for v_y and t_y are taken from
! Forbes, R. G., & Deane, J. H. (2007, November).
! "Reformulation of the standard theory of Fowler–Nordheim tunnelling and cold field electron emission."
! In Proceedings of the Royal Society of London A: Mathematical,
! Physical and Engineering Sciences (Vol. 463, No. 2087, pp. 2907-2927). The Royal Society.
!
  double precision function v_y(F, pos)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos
    integer                                      :: emit = 1
    double precision                             :: l

    if (image_charge .eqv. .true.) then
      l = l_const * (-1.0d0*F) / w_theta_xy(pos, emit)**2 ! l = y^2, y = 3.79E-4 * sqrt(F_old) / w_theta
      if (l > 1.0d0) then
        l = 1.0d0
      end if
      v_y = 1.0d0 - l + 1.0d0/6.0d0 * l * log(l)
    else
      v_y = 1.0d0
    end if
  end function v_y

  double precision function t_y(F, pos)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos
    integer                                      :: emit = 1
    double precision                             :: l

    if (image_charge .eqv. .true.) then
      l = l_const * (-1.0d0*F) / w_theta_xy(pos, emit)**2 ! l = y^2, y = 3.79E-4 * sqrt(F_old) / w_theta
      if (l > 1.0d0) then
        print *, 'Error: l > 1.0'
        print *, 'l = ', l, ', F = ', F, ', t_y = ', t_y
        print *, 'x = ', pos(1)/length_scale, 'y = ', pos(2)/length_scale, ' z = ,', pos(3)/length_scale
        l = 1.0d0
        !call Write_Current_Position()
        !stop
      end if

      t_y = 1.0d0 + l*(1.0d0/9.0d0 - 1.0d0/18.0d0*log(l))
    else
      t_y = 1.0d0
    end if
  end function t_y

  !-----------------------------------------------------------------------------
  ! This function returns the escape probability of the Electrons.
  ! Escape_prob = exp(-b_FN*w_theta^(3/2)*v_y/F) .
  double precision function Escape_Prob(F, pos)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos
    integer                                      :: emit = 1

    Escape_Prob = exp(b_FN * (sqrt(w_theta_xy(pos, emit)))**3 * v_y(F, pos) / (-1.0d0*F))

  end function Escape_Prob

  ! Returns the log of the escape probability.
  double precision function Escape_Prob_log(F, pos)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos
    integer                                      :: emit = 1

    Escape_Prob_log = b_FN * (sqrt(w_theta_xy(pos, emit)))**3 * v_y(F, pos) / (-1.0d0*F)

  end function Escape_Prob_log

  !-----------------------------------------------------------------------------
  ! A simple function that calculates
  ! A_FN/(t**2(l)*w_theta(x,y)) F**2(x,y)
  ! pos: Position to calculate the function
  ! F: The z-component of the field at par_pos, it should be F < 0.0d0.
  double precision function Elec_Supply_V2(F, pos)
    double precision, dimension(1:3), intent(in) :: pos
    double precision,                 intent(in) :: F
    integer                                      :: emit = 1

    Elec_Supply_V2 = time_step_div_q0 * a_FN/(t_y(F, pos)**2*w_theta_xy(pos, emit)) * F**2
  end function Elec_Supply_V2

  !-----------------------------------------------------------------------------
  ! Gives the component of the field normal to the surface
  pure function Field_normal(par_xyz_pos, par_pos_org, field_E, sec_in)
    double precision, dimension(1:3), intent(in) :: par_xyz_pos, par_pos_org, field_E
    integer                         , intent(in) :: sec_in
    double precision                             :: Field_normal
    double precision, dimension(1:3)             :: unit_vec

    unit_vec = surface_normal(par_xyz_pos, par_pos_org, sec_in)
    Field_normal = dot_product(unit_vec, field_E)
  end function Field_normal

  !-----------------------------------------------------------------------------
  ! Find a unit vector normal to the surface
  pure function surface_normal(par_pos, par_pos_org, sec_in)
    double precision, dimension(1:3)             :: surface_normal
    double precision, dimension(1:3), intent(in) :: par_pos, par_pos_org
    integer                         , intent(in) :: sec_in

    if (sec_in == sec_top) then
        ! Top surface
        ! The normal vector is (0, 0, 1)
        surface_normal(1) = 0.0d0
        surface_normal(2) = 0.0d0
        surface_normal(3) = 1.0d0
    else if (sec_in == sec_side) then
        ! Side surface
        ! The normal vector is (cos(phi), sin(phi), 0)
        surface_normal(1) = cos(par_pos(2))
        surface_normal(2) = sin(par_pos(2))
        surface_normal(3) = 0.0d0
    else if (sec_in == sec_corner) then
        ! Corner surface
        ! The normal vector is (cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta))
        surface_normal(1) = cos(par_pos(2)) * cos(par_pos(3))
        surface_normal(2) = sin(par_pos(2)) * cos(par_pos(3))
        surface_normal(3) = sin(par_pos(3))
    else
        ! Error
        surface_normal(:) = ieee_value( surface_normal(:), ieee_signaling_nan )
    end if
  end function surface_normal

  subroutine Do_Cuba_Divonne(N_sup)
    implicit none
    double precision, intent(out) :: N_sup
    double precision              :: N_top, N_side, N_corner
    integer                       :: i
    integer                       :: IFAIL
  
    
    ! Cuba integration variables (common)
    integer, parameter :: ndim = 2 ! Number of dimensions
    integer, parameter :: ncomp = 1 ! Number of components in the integrand
    integer            :: userdata = 0 ! User data passed to the integrand
    integer, parameter :: nvec = 1 ! Number of points given to the integrand function
    !double precision   :: epsrel = 1.0d-4 ! Requested relative error
    !double precision   :: epsabs = 0.25d0 ! Requested absolute error
    integer            :: flags = 0+4 ! Flags
    integer            :: seed = 0 ! Seed for the rng. Zero will use Sobol.
    !integer            :: mineval = 75000 ! Minimum number of integrand evaluations
    !integer            :: maxeval = 10000000 ! Maximum number of integrand evaluations
  
    ! Divonne specific
    integer :: key1 = 47 ! 〈in〉, determines sampling in the partitioning phase:
                    ! key1 = 7,9,11,13 selects the cubature rule of degree key1.
                    ! Note that the degree-11 rule is available only in 3 dimensions, the degree-13 rule only in 2 dimensions.
                    ! For other values of key1, a quasi-random sample of n_1=|key1| points is used,
                    ! where the sign of key1 determines the type of sample,
                    ! – key1 > 0, use a Korobov quasi-random sample,
                    ! – key1 < 0, use a “standard” sample (a Sobol quasi-random sample if seed= 0, otherwise a pseudo-random sample).
    integer :: key2 = 1 !<in>, determines sampling in the final integration phase:
                    ! key2 = 7, 9, 11, 13 selects the cubature rule of degree key2. Note that the degree-11
                    ! rule is available only in 3 dimensions, the degree-13 rule only in 2 dimensions.
                    ! For other values of key2, a quasi-random sample is used, where the sign of key2
                    ! determines the type of sample,
                    ! - key2 > 0, use a Korobov quasi-random sample,
                    ! - key2 < 0, use a "standard" sample (see description of key1 above),
                    ! and n_2 = |key2| determines the number of points,
                    ! - n_2 > 40, sample n2 points,
                    ! - n_2 < 40, sample n2 nneed points, where nneed is the number of points needed to
                    ! reach the prescribed accuracy, as estimated by Divonne from the results of the
                    ! partitioning phase.
    integer :: key3 = -1 ! <in>, sets the strategy for the refinement phase:
                    ! key3 = 0, do not treat the subregion any further.
                    ! key3 = 1, split the subregion up once more.
                    ! Otherwise, the subregion is sampled a third time with key3 specifying the sampling
                    ! parameters exactly as key2 above.
    integer :: maxpass = 2! <in>, controls the thoroughness of the partitioning phase: The
  ! partitioning phase terminates when the estimated total number of integrand evaluations
  ! (partitioning plus final integration) does not decrease for maxpass successive iterations.
  ! A decrease in points generally indicates that Divonne discovered new structures of
  ! the integrand and was able to find a more effective partitioning. maxpass can be
  ! understood as the number of `safety' iterations that are performed before the partition
  ! is accepted as final and counting consequently restarts at zero whenever new structures are found.
    double precision :: border = 0.0d0 ! <in>, the width of the border of the integration region.
  ! Points falling into this border region will not be sampled directly, but will be extrapolated
  ! from two samples from the interior. Use a non-zero border if the integrand
  ! subroutine cannot produce values directly on the integration boundary.
    double precision :: maxchisq = 10.0d0 !<in>, the maximum \chi^2 value a single subregion is allowed
  ! to have in the final integration phase. Regions which fail this \chi^2 test and whose
  ! sample averages differ by more than mindeviation move on to the refinement phase.
    double precision :: mindeviation = 0.25d0 !<in>, a bound, given as the fraction of the requested
  ! error of the entire integral, which determines whether it is worthwhile further
  ! examining a region that failed the \chi^2 test. Only if the two sampling averages
  ! obtained for the region differ by more than this bound is the region further treated.
    integer :: ngiven = 0 ! <in>, the number of points in the xgiven array.
    integer :: ldxgiven = ndim ! <in>, the leading dimension of xgiven, i.e. the offset between one point and the next in memory.
    !double precision, allocatable :: xgiven(:,:) ! xgiven(ldxgiven,ngiven) <in>, a list of points where the integrand
  ! might have peaks. Divonne will consider these points when partitioning the
  ! integration region. The idea here is to help the integrator find the extrema of the integrand
  ! in the presence of very narrow peaks. Even if only the approximate location
  ! of such peaks is known, this can considerably speed up convergence.
    integer :: nextra = 0 ! <in>, the maximum number of extra points the peak-finder subroutine
  ! will return. If nextra is zero, peakfinder is not called and an arbitrary object
  ! may be passed in its place, e.g. just 0.
  
  
    ! Output
    character          :: statefile = "" ! File to save the state in. Empty string means don't do it.
    integer            :: spin = -1 ! Spinning cores
    integer            :: nregions ! <out> The actual number of subregions needed
    integer            :: neval ! <out> The actual number of integrand evaluations needed
    integer            :: fail ! <out> Error flag (0 = Success, -1 = Dimension out of range, >0 = Accuracy goal was not met)
    double precision, dimension(1:ncomp) :: integral ! <out> The integral of the integrand over the unit hybercube
    double precision, dimension(1:ncomp) :: error ! <out> The presumed absolute error
    double precision, dimension(1:ncomp) :: prob ! <out> The chi-square probability
  
    ! Initialize the average field to zero
    F_avg = 0.0d0

    !-----------------------------------------------------------------------------
    ! Integrate over the top surface
    ! Pass the number of the emitter being integrated over to the integrand as userdata
    userdata = sec_top
   
    call divonne(ndim, ncomp, integrand_cuba_cyl, userdata, nvec,&
                cuba_epsrel, cuba_epsabs, flags, seed, cuba_mineval, cuba_maxeval,&
                key1, key2, key3, maxpass,&
                border, maxchisq, mindeviation,&
                ngiven, ldxgiven, 0, nextra, 0,&
                statefile, spin,&
                nregions, neval, fail, integral, error, prob)
  
    if (fail /= 0) then
      if (abs(error(1) - cuba_epsabs) > 1.0d-2) then
        print '(a)', 'RUMDEED: WARNING Cuba did not return 0'
        print *, 'Fail = ', fail
        print *, 'nregions = ', nregions
        print *, 'neval = ', neval, ' max is ', cuba_maxeval
        print *, 'error(1) = ', error(1)
        print *, 'integral(1)*epsrel = ', integral(1)*cuba_epsrel
        print *, 'epsabs = ', cuba_epsabs
        print *, 'prob(1) = ', prob(1)
        print *, 'integral(1) = ', integral(1)
        call Flush_Data()
      end if
    end if
  
    ! Store the results
    N_top = integral(1)

    ! Write the output variables of the integration to a file
    write(ud_integrand, '(i8, tr2, i8, tr2, i4, tr2, ES12.4, tr2, ES12.4, tr2, ES12.4)', iostat=IFAIL) &
                            & nregions, neval, fail, integral(1), error(1), prob(1)

    !-----------------------------------------------------------------------------
    ! Integrate over the side surface
    userdata = sec_side
   
    call divonne(ndim, ncomp, integrand_cuba_cyl, userdata, nvec,&
                 cuba_epsrel, cuba_epsabs, flags, seed, cuba_mineval, cuba_maxeval,&
                 key1, key2, key3, maxpass,&
                 border, maxchisq, mindeviation,&
                 ngiven, ldxgiven, 0, nextra, 0,&
                 statefile, spin,&
                 nregions, neval, fail, integral, error, prob)
   
    if (fail /= 0) then
      if (abs(error(1) - cuba_epsabs) > 1.0d-2) then
        print '(a)', 'RUMDEED: WARNING Cuba did not return 0'
        print *, 'Fail = ', fail
        print *, 'nregions = ', nregions
        print *, 'neval = ', neval, ' max is ', cuba_maxeval
        print *, 'error(1) = ', error(1)
        print *, 'integral(1)*epsrel = ', integral(1)*cuba_epsrel
        print *, 'epsabs = ', cuba_epsabs
        print *, 'prob(1) = ', prob(1)
        print *, 'integral(1) = ', integral(1)
        call Flush_Data()
      end if
    end if
   
    ! Store the results
    N_side = integral(1)

    ! Write the output variables of the integration to a file
    write(ud_integrand, '(i8, tr2, i8, tr2, i4, tr2, ES12.4, tr2, ES12.4, tr2, ES12.4)', iostat=IFAIL) &
                            & nregions, neval, fail, integral(1), error(1), prob(1)

    !-----------------------------------------------------------------------------
    ! Integrate over the corner surface
    userdata = sec_corner
   
    call divonne(ndim, ncomp, integrand_cuba_cyl, userdata, nvec,&
                 cuba_epsrel, cuba_epsabs, flags, seed, cuba_mineval, cuba_maxeval,&
                 key1, key2, key3, maxpass,&
                 border, maxchisq, mindeviation,&
                 ngiven, ldxgiven, 0, nextra, 0,&
                 statefile, spin,&
                 nregions, neval, fail, integral, error, prob)
   
    if (fail /= 0) then
      if (abs(error(1) - cuba_epsabs) > 1.0d-2) then
        print '(a)', 'RUMDEED: WARNING Cuba did not return 0'
        print *, 'Fail = ', fail
        print *, 'nregions = ', nregions
        print *, 'neval = ', neval, ' max is ', cuba_maxeval
        print *, 'error(1) = ', error(1)
        print *, 'integral(1)*epsrel = ', integral(1)*cuba_epsrel
        print *, 'epsabs = ', cuba_epsabs
        print *, 'prob(1) = ', prob(1)
        print *, 'integral(1) = ', integral(1)
        call Flush_Data()
      end if
    end if
   
    ! Store the results
    N_side = integral(1)

    ! Finish calculating the average field
    F_avg = F_avg / neval

    ! Write the output variables of the integration to a file
    write(ud_integrand, '(i8, tr2, i8, tr2, i4, tr2, ES12.4, tr2, ES12.4, tr2, ES12.4)', iostat=IFAIL) &
                          & nregions, neval, fail, integral(1), error(1), prob(1)

    ! Calculate the total number of electrons emitted
    !print *, 'N_top = ', N_top
    !print *, 'N_side = ', N_side
    !print *, 'N_corner = ', N_corner
    N_sup = N_top + N_side + N_corner
    !print *, 'N_sup = ', N_sup
    
    ! Set probabilities of emitting from different surfaces
    prob_top = N_top / N_sup
    prob_side = N_side / N_sup
    prob_corner = N_corner / N_sup
  end subroutine Do_Cuba_Divonne

  subroutine Do_Cuba_Suave(N_sup)
    implicit none
   ! Input / output variables
    double precision, intent(out) :: N_sup
    double precision              :: N_top, N_side, N_corner
    integer                       :: IFAIL

    ! Cuba integration variables
    integer, parameter :: ndim = 2 ! Number of dimensions
    integer, parameter :: ncomp = 1 ! Number of components in the integrand
    integer            :: userdata = 0 ! User data passed to the integrand
    integer, parameter :: nvec = 1 ! Number of points given to the integrand function
    double precision   :: epsrel = 1.0d-1 ! Requested relative error
    double precision   :: epsabs = 1.0d-3 ! Requested absolute error
    integer            :: flags = 0+4 ! Flags
    integer            :: seed = 0 ! Seed for the rng. Zero will use Sobol.
    integer            :: mineval = 100000 ! Minimum number of integrand evaluations
    integer            :: maxeval = 5000000 ! Maximum number of integrand evaluations
    integer            :: nnew = 5000 ! Number of integrand evaluations in each subdivision
    integer            :: nmin = 1000 ! Minimum number of samples a former pass must contribute to a subregion to be considered in the region's compound integral value.
    double precision   :: flatness = 10.0d0 ! Determine how prominently out-liers, i.e. samples with a large fluctuation, 
                                           ! figure in the total fluctuation, which in turn determines how a region is split up.
                                           ! As suggested by its name, flatness should be chosen large for 'flat" integrand and small for 'volatile' integrands
                                           ! with high peaks.
    character          :: statefile = "" ! File to save the state in. Empty string means don't do it.
    integer            :: spin = -1 ! Spinning cores
    integer            :: nregions ! <out> The actual number of subregions needed
    integer            :: neval ! <out> The actual number of integrand evaluations needed
    integer            :: fail ! <out> Error flag (0 = Success, -1 = Dimension out of range, >0 = Accuracy goal was not met)
    double precision, dimension(1:ncomp) :: integral ! <out> The integral of the integrand over the unit hybercube
    double precision, dimension(1:ncomp) :: error ! <out> The presumed absolute error
    double precision, dimension(1:ncomp) :: prob ! <out> The chi-square probability


    ! Initialize the average field to zero
    F_avg = 0.0d0

    ! ------------------------------------------------------------------------------
    ! Integrate over the top surface
    userdata = sec_top

    call suave(ndim, ncomp, integrand_cuba_cyl, userdata, nvec, &
     & epsrel, epsabs, flags, seed, &
     & mineval, maxeval, nnew, nmin, flatness, &
     & statefile, spin, &
     & nregions, neval, fail, integral, error, prob)

    if (fail /= 0) then
      print '(a)', 'RUMDEED: WARNING Cuba did not return 0 (TOP)'
      print *, fail
      print *, error
      print *, integral(1)*cuba_epsrel
      print *, cuba_epsabs
      print *, prob
      print *, integral(1)
      call Flush_Data()
    end if


    !! Round the results to the nearest integer
    !N_sup = nint( integral(1) )
    N_top = integral(1)

    ! Write the output variables of the integration to a file
    write(ud_integrand, '(i8, tr2, i8, tr2, i4, tr2, ES12.4, tr2, ES12.4, tr2, ES12.4)', iostat=IFAIL) &
                          & nregions, neval, fail, integral(1), error(1), prob(1)

    ! ------------------------------------------------------------------------------
    ! Integrate over the side surface
    userdata = sec_side

    call suave(ndim, ncomp, integrand_cuba_cyl, userdata, nvec, &
      & epsrel, epsabs, flags, seed, &
      & mineval, maxeval, nnew, nmin, flatness, &
      & statefile, spin, &
      & nregions, neval, fail, integral, error, prob)

    if (fail /= 0) then
      print '(a)', 'RUMDEED: WARNING Cuba did not return 0 (SIDE)'
      print *, fail
      print *, error
      print *, integral(1)*cuba_epsrel
      print *, cuba_epsabs
      print *, prob
      print *, integral(1)
      call Flush_Data()
    end if


    !! Round the results to the nearest integer
    !N_sup = nint( integral(1) )
    N_side = integral(1)

    ! Write the output variables of the integration to a file
    write(ud_integrand, '(i8, tr2, i8, tr2, i4, tr2, ES12.4, tr2, ES12.4, tr2, ES12.4)', iostat=IFAIL) &
                          & nregions, neval, fail, integral(1), error(1), prob(1)

    ! ------------------------------------------------------------------------------
    ! Integrate over the side surface
    userdata = sec_corner

    call suave(ndim, ncomp, integrand_cuba_cyl, userdata, nvec, &
      & epsrel, epsabs, flags, seed, &
      & mineval, maxeval, nnew, nmin, flatness, &
      & statefile, spin, &
      & nregions, neval, fail, integral, error, prob)

    if (fail /= 0) then
      print '(a)', 'RUMDEED: WARNING Cuba did not return 0 (CORNER)'
      print *, fail
      print *, error
      print *, integral(1)*cuba_epsrel
      print *, cuba_epsabs
      print *, prob
      print *, integral(1)
      call Flush_Data()
    end if


    !! Round the results to the nearest integer
    !N_sup = nint( integral(1) )
    N_corner = integral(1)

    ! Write the output variables of the integration to a file
    write(ud_integrand, '(i8, tr2, i8, tr2, i4, tr2, ES12.4, tr2, ES12.4, tr2, ES12.4)', iostat=IFAIL) &
                          & nregions, neval, fail, integral(1), error(1), prob(1)

    ! Finish calculating the average field
    F_avg = F_avg / neval

    ! Calculate the total number of electrons emitted
    !print *, 'N_top = ', N_top
    !print *, 'N_side = ', N_side
    !print *, 'N_corner = ', N_corner
    N_sup = N_top + N_side + N_corner
    !print *, 'N_sup = ', N_sup

    ! Set probabilities of emitting from different surfaces
    prob_top = N_top / N_sup
    prob_side = N_side / N_sup
    prob_corner = N_corner / N_sup
  end subroutine Do_Cuba_Suave

  ! ----------------------------------------------------------------------------
  ! The integration function for the Cuba library
  !
  integer function integrand_cuba_cyl(ndim, xx, ncomp, ff, userdata)
    ! Input / output variables
    integer, intent(in) :: ndim ! Number of dimensions (Should be 2)
    integer, intent(in) :: ncomp ! Number of vector-components in the integrand (Always 1 here)
    integer, intent(in) :: userdata ! Additional data passed to the integral function (In our case the number of the emitter)
    double precision, intent(in), dimension(1:ndim)   :: xx ! Integration points, between 0 and 1
    double precision, intent(out), dimension(1:ncomp) :: ff ! Results of the integrand function

    ! Variables used for calculations
    double precision, dimension(1:3) :: par_pos, par_pos_org, field
    double precision                 :: field_norm

    ! Variables used for the transformation from the hypercube to the surface
    double precision, dimension(1:2) :: upper, lower, range
    double precision                 :: jacobian ! Jacobian of the transformation from the hypercube to the surface

    !print *, 'Integrand called'
    ! Surface position
    ! Cuba does the integration over the hybercube.
    ! It gives us coordinates between 0 and 1.
    ! Check which section the emitter is in
    if (userdata == sec_top) then
      !print *, 'Top surface'
        ! Use cylindrical coordinates for the top surface here to avoid rejecting points
        ! r [0, radius_cyl - radius_cor]
        ! phi [0, 2*pi]
        upper(1:2) = (/radius_cyl - radius_cor, 2.0d0*pi/)
        lower(1:2) = (/0.0d0, 0.0d0/)

        !upper(1:2) = 1.0d0*(radius_cyl - radius_cor)/length_scale
        !lower(1:2) = -1.0d0*(radius_cyl - radius_cor)/length_scale

        ! Range
        range(1:2) = upper(1:2) - lower(1:2)

        ! r and phi coordinates of the point on the surface
        par_pos_org(1:2) = lower(1:2) + xx(1:2)*range(1:2)
        par_pos_org(3) = height_cyl ! z

        ! Convert to cartesian coordinates
        par_pos(1) = xx(1)*upper(1)*cos(xx(2)*upper(2))
        par_pos(2) = xx(1)*(radius_cyl-radius_cor)*sin(xx(2)*2.0d0*pi)
        par_pos(3) = height_cyl

        ! Jacobian of the transformation from the hypercube
        jacobian = par_pos_org(1)*range(1)*range(2) ! dx dy = r dr dphi
        
        !ff(1) = 1.0d0 ! Integrand results debug
        
        ! ! Check if the point is inside the circle on the top surface
        ! if ((scaledx(1)**2 + scaledx(2)**2) <= ((radius_cyl-radius_cor)/length_scale)**2) then
        !     ff(1) = 1.0d0
        ! else
        !     ff(1) = 0.0d0 ! Point is outside the circle reject it
        ! end if
        integrand_cuba_cyl = 0 ! Return value to Cuba, 0 = success

        !write(ud_cyl_debug, *) par_pos(1:3), par_pos_org(1:3), userdata
    else if (userdata == sec_side) then
      !print *, 'Side surface'
        ! Use cylindrical coordinates for the side surface
        ! phi [0, 2*pi]
        ! z [0, height_cyl - radius_cor]
        upper(1:2) = (/2.0d0*pi, height_cyl - radius_cor/)
        lower(1:2) = (/0.0d0, 0.0d0/)

        ! Range
        range(1:2) = upper(1:2) - lower(1:2)

        ! rho, phi and z coordinates of the point on the surface
        par_pos_org(1) = radius_cyl
        par_pos_org(2:3) = lower(1:2) + xx(1:2)*range(1:2)

        ! Convert to cartesian coordinates
        par_pos(1) = radius_cyl*cos(par_pos_org(2))
        par_pos(2) = radius_cyl*sin(par_pos_org(2))
        par_pos(3) = par_pos_org(3)

        ! Jacobian of the transformation from the hypercube
        jacobian = radius_cyl*range(1)*range(2) ! r dphi dz

        !ff(1) = 1.0d0 ! Integrand results debug

        integrand_cuba_cyl = 0 ! Return value to Cuba, 0 = success

        !write(ud_cyl_debug, *) par_pos(1:3), par_pos_org(1:3), userdata
    else if (userdata == sec_corner) then
      !print *, 'Corner surface'
        ! Use toroidal coordinates for the corner surface
        ! phi [0, 2*pi]
        ! theta [0, pi/2]
        upper(1:2) = (/2.0d0*pi, pi/2.0d0/)
        lower(1:2) = (/0.0d0, 0.0d0/)

        ! Range
        range(1:2) = upper(1:2) - lower(1:2)

        ! rho, phi and theta coordinates of the point on the surface
        par_pos_org(1) = radius_cor
        par_pos_org(2:3) = lower(1:2) + xx(1:2)*range(1:2)

        ! Convert to cartesian coordinates
        par_pos(1) = (radius_cyl - radius_cor + radius_cor*cos(par_pos_org(3)))*cos(par_pos_org(2))
        par_pos(2) = (radius_cyl - radius_cor + radius_cor*cos(par_pos_org(3)))*sin(par_pos_org(2))
        par_pos(3) = (height_cyl - radius_cor) + radius_cor*sin(par_pos_org(3))

        ! Jacobian of the transformation from the hypercube
        jacobian = radius_cor*(radius_cyl - radius_cor + radius_cor*cos(par_pos_org(3)))*range(1)*range(2) ! rho*(R_0 + rho*cos(theta)) dphi dtheta

        !ff(1) = 1.0d0 ! Integrand results debug

        integrand_cuba_cyl = 0 ! Return value to Cuba, 0 = success

        !write(ud_cyl_debug, *) par_pos(1:3), par_pos_org(1:3), userdata
    else
        integrand_cuba_cyl = 100 ! Error
        !print *, 'RUMDEED: Error Unknown section (integrand_cuba)'
        !stop
    end if

    ! Calculate the electric field on the surface
    field = Calc_Field_at(par_pos)
    field_norm = Field_normal(par_pos, par_pos_org, field, userdata)

    ! Add to the average field
    F_avg = F_avg + field

    ! Check if the field is favorable for emission
    if (field_norm < 0.0d0) then
      ! The field is favorable for emission
      ! Calculate the electron supply at this point
      ff(1) = Elec_Supply_V2(field_norm, par_pos)
      !print *, 'Electron supply = ', ff(1)
    else
    !   ! The field is NOT favorable for emission
    !   ! This point does not contribute
      !print *, 'Field not favorable for emission'
      !print *, 'field_norm = ', field_norm
      !print *, 'par_pos = ', par_pos/length_scale_cyl

      ff(1) = 0.0d0
      !pause
    end if

    ! We multiply with the jacobian because Cuba does the 
    ! integration over the hybercube, i.e. from 0 to 1.
    ff(1) = jacobian*ff(1)
    !if (userdata == sec_corner) then
    !  print *, 'ff(1) = ', ff(1)
    !end if
  end function integrand_cuba_cyl

!-------------------------------------------!
! Tests

!-----------------------------------------------------------------------
! Calculate the electric field along the edge of the cylinder
! Coordinates are cylindrical (rho, phi, z)
subroutine Calc_E_edge_cyl()
    implicit none
    ! Local
    integer :: k
    integer, parameter :: N_p = 100
    real(kdkind), dimension(1:3) :: p, E_vec
    double precision :: phi, rho
    integer :: ud_cyl_side_left, ud_cyl_side_right, IFAIL
    character (len=*), parameter :: filename_cyl_side_left="cyl_side_left.dt"
    character (len=*), parameter :: filename_cyl_side_right="cyl_side_right.dt"

    ! Data arrays
    double precision, dimension(:, :), allocatable :: data_cyl_side_left, data_cyl_side_right
    double precision, dimension(:, :), allocatable :: p_cyl_side_left, p_cyl_side_right
    double precision, dimension(:), allocatable    :: len_cyl_side_left, len_cyl_side_right

    print *, 'Calc_E_edge_cyl started'

    allocate(data_cyl_side_right(1:3, 1:N_p))
    allocate(data_cyl_side_left(1:3, 1:N_p))
    allocate(len_cyl_side_right(1:N_p))
    allocate(len_cyl_side_left(1:N_p))
    allocate(p_cyl_side_right(1:3, 1:N_p))
    allocate(p_cyl_side_left(1:3, 1:N_p))

    rho = radius_cyl

    ! Calculate the electric field along the edges of the cylinder
    print *, 'Calculating the electric field along the edges of the cylinder'
    do k = 1, N_p
        ! Right edge
        phi = 0.0d0
        p(1) = rho*cos(phi)
        p(2) = rho*sin(phi)
        p(3) = (k-1)*(height_cyl-radius_cor)/(N_p-1)

        p_cyl_side_right(:, k) = p(:)
        len_cyl_side_right(k) = p(3)

        E_vec = field_E_cylinder(p)
        ! Store data in array
        data_cyl_side_right(:, k) = E_vec(:)

        ! Left edge
        phi = pi
        p(1) = rho*cos(phi)
        p(2) = rho*sin(phi)
        p(3) = (k-1)*(height_cyl-radius_cor)/(N_p-1)

        p_cyl_side_left(:, k) = p(:)
        len_cyl_side_left(k) = (height_cyl - radius_cor) - p(3)

        E_vec = field_E_cylinder(p)
        ! Store data in array
        data_cyl_side_left(:, k) = E_vec(:)
    end do

    ! Open data files
    print *, 'Opening data files'
    open(newunit=ud_cyl_side_left, iostat=IFAIL, file=filename_cyl_side_left, status='REPLACE', action='WRITE')
    if (IFAIL /= 0) then
      print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file ', filename_cyl_side_left
      stop
    end if
    open(newunit=ud_cyl_side_right, iostat=IFAIL, file=filename_cyl_side_right, status='REPLACE', action='WRITE')
    if (IFAIL /= 0) then
      print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file ', filename_cyl_side_right
      stop
    end if

    ! Write data to files
    print *, 'Writing data to files'
    do k = 1, N_p
      write(ud_cyl_side_left, *) data_cyl_side_left(:, k), len_cyl_side_left(k), p_cyl_side_left(:, k)
      write(ud_cyl_side_right, *) data_cyl_side_right(:, k), len_cyl_side_right(k), p_cyl_side_right(:, k)
    end do

    ! Close data files
    close(unit=ud_cyl_side_left, iostat=IFAIL, status='keep')
    close(unit=ud_cyl_side_right, iostat=IFAIL, status='keep')

    deallocate(data_cyl_side_left)
    deallocate(data_cyl_side_right)
    deallocate(p_cyl_side_left)
    deallocate(p_cyl_side_right)
    deallocate(len_cyl_side_left)
    deallocate(len_cyl_side_right)

    print *, 'Calc_E_edge_cyl finished'
end subroutine Calc_E_edge_cyl

    !-----------------------------------------------------------------------
    ! Calculate the electric field along the corner of the cylinder
    ! Coordinates are torodial (rho, phi, theta)
subroutine Calc_E_corner_cyl()
    implicit none
    ! Local
    integer :: k
    integer, parameter :: N_p = 100
    real(kdkind), dimension(1:3) :: p, E_vec
    double precision :: phi, theta, rho
    integer :: ud_cyl_corner_left, ud_cyl_corner_right, IFAIL
    
    character (len=*), parameter :: filename_cyl_corner_left="cyl_corner_left.dt"
    character (len=*), parameter :: filename_cyl_corner_right="cyl_corner_right.dt"

    ! Data arrays
    double precision, dimension(:, :), allocatable :: data_cyl_corner_left, data_cyl_corner_right
    double precision, dimension(:, :), allocatable :: p_cyl_corner_left, p_cyl_corner_right
    double precision, dimension(:), allocatable    :: len_cyl_corner_left, len_cyl_corner_right

    allocate(data_cyl_corner_left(1:3, 1:N_p))
    allocate(data_cyl_corner_right(1:3, 1:N_p))
    allocate(len_cyl_corner_left(1:N_p))
    allocate(len_cyl_corner_right(1:N_p))
    allocate(p_cyl_corner_left(1:3, 1:N_p))
    allocate(p_cyl_corner_right(1:3, 1:N_p))

    rho = radius_cor

    ! Calculate the electric field along the corner of the cylinder
    do k = 1, N_p
        ! Right corner
        phi = 0.0d0
        theta = (k-1)*pi/(2.0d0*(N_p-1))
        p(1) = (radius_cyl - radius_cor + rho*cos(theta)) * cos(phi)
        p(2) = (radius_cyl - radius_cor + rho*cos(theta)) * sin(phi)
        p(3) = (height_cyl - radius_cor) + rho*sin(theta)

        p_cyl_corner_right(:, k) = p(:)
        len_cyl_corner_right(k) = radius_cor*theta

        E_vec = field_E_cylinder(p)
        ! Store data in array
        data_cyl_corner_right(:, k) = E_vec(:)

        ! Left corner
        phi = pi
        theta = (k-1)*pi/(2.0d0*(N_p-1))
        p(1) = (radius_cyl - radius_cor + rho*cos(theta)) * cos(phi)
        p(2) = (radius_cyl - radius_cor + rho*cos(theta)) * sin(phi)
        p(3) = (height_cyl - radius_cor) + rho*sin(theta)

        p_cyl_corner_left(:, k) = p(:)
        len_cyl_corner_left(k) = radius_cor*(pi/2.0d0 - theta)

        E_vec = field_E_cylinder(p)
        ! Store data in array
        data_cyl_corner_left(:, k) = E_vec(:)
    end do

        ! Open data files
    open(newunit=ud_cyl_corner_left, iostat=IFAIL, file=filename_cyl_corner_left, status='REPLACE', action='WRITE')
    if (IFAIL /= 0) then
      print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file ', filename_cyl_corner_left
      stop
    end if
    open(newunit=ud_cyl_corner_right, iostat=IFAIL, file=filename_cyl_corner_right, status='REPLACE', action='WRITE')
    if (IFAIL /= 0) then
      print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file ', filename_cyl_corner_right
      stop
    end if

    ! Write data to files
    do k = 1, N_p
        write(ud_cyl_corner_left, *) data_cyl_corner_left(:, k), len_cyl_corner_left(k), p_cyl_corner_left(:, k)
        write(ud_cyl_corner_right, *) data_cyl_corner_right(:, k), len_cyl_corner_right(k), p_cyl_corner_right(:, k)
    end do

    ! Close data files
    close(unit=ud_cyl_corner_left, iostat=IFAIL, status='keep')
    close(unit=ud_cyl_corner_right, iostat=IFAIL, status='keep')

    deallocate(data_cyl_corner_left)
    deallocate(data_cyl_corner_right)
    deallocate(p_cyl_corner_left)
    deallocate(p_cyl_corner_right)
    deallocate(len_cyl_corner_left)
    deallocate(len_cyl_corner_right)
end subroutine Calc_E_corner_cyl

!-----------------------------------------------------------------------
! Calculate the electric field along the top of the cylinder
! Coordinates are cylindrical (rho, phi, z)
subroutine Calc_E_top_cyl()
    implicit none
    ! Local
    integer :: k
    integer, parameter :: N_p = 100
    real(kdkind), dimension(1:3) :: p, E_vec
    double precision :: rho_max
    integer :: ud_cyl_top_left, ud_cyl_top_right, IFAIL
    character (len=*), parameter :: filename_cyl_top_left="cyl_top_left.dt"
    character (len=*), parameter :: filename_cyl_top_right="cyl_top_right.dt"

    ! Data arrays
    double precision, dimension(:, :), allocatable :: data_cyl_top_left, data_cyl_top_right
    double precision, dimension(:, :), allocatable :: p_cyl_top_left, p_cyl_top_right
    double precision, dimension(:), allocatable    :: len_cyl_top_left, len_cyl_top_right

    allocate(data_cyl_top_left(1:3, 1:N_p))
    allocate(data_cyl_top_right(1:3, 1:N_p))
    allocate(len_cyl_top_left(1:N_p))
    allocate(len_cyl_top_right(1:N_p))
    allocate(p_cyl_top_left(1:3, 1:N_p))
    allocate(p_cyl_top_right(1:3, 1:N_p))

    rho_max = radius_cyl - radius_cor

    ! Calculate the electric field along the top of the cylinder
    do k = 1, N_p
        ! Right side
        !phi = 0.0d0
        !rho = (k-1)*rho_max/(N_p-1)
        !p(1) = rho*cos(phi)
        !p(2) = rho*sin(phi)
        !p(3) = height_cyl
        p(1) = (k-1)*(radius_cyl-radius_cor)/(N_p-1)
        p(2) = 0.0d0
        p(3) = height_cyl

        p_cyl_top_right(:, k) = p(:)
        len_cyl_top_right(k) = (radius_cyl - radius_cor) - p(1)

        E_vec = field_E_cylinder(p)
        ! Store data in array
        data_cyl_top_left(:, k) = E_vec(:)

        ! Left side
        !phi = pi
        !rho = (k-1)*rho_max/(N_p-1)
        !p(1) = rho*cos(phi)
        !p(2) = rho*sin(phi)
        !p(3) = height_cyl
        p(1) = (k-1)*(radius_cyl-radius_cor)/(N_p-1)*(-1.0d0)
        p(2) = 0.0d0
        p(3) = height_cyl

        p_cyl_top_left(:, k) = p(:)
        len_cyl_top_left(k) = abs(p(1))

        E_vec = field_E_cylinder(p)
        ! Store data in array
        data_cyl_top_right(:, k) = E_vec(:)
    end do

    ! Open data files
    open(newunit=ud_cyl_top_left, iostat=IFAIL, file=filename_cyl_top_left, status='REPLACE', action='WRITE')
    if (IFAIL /= 0) then
      print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file ', filename_cyl_top_left
      stop
    end if
    open(newunit=ud_cyl_top_right, iostat=IFAIL, file=filename_cyl_top_right, status='REPLACE', action='WRITE')
    if (IFAIL /= 0) then
      print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file ', filename_cyl_top_right
      stop
    end if

    ! Write data to files
    do k = 1, N_p
        write(ud_cyl_top_left, *) data_cyl_top_left(:, k), len_cyl_top_left(k), p_cyl_top_left(:, k)
        write(ud_cyl_top_right, *) data_cyl_top_right(:, k), len_cyl_top_right(k), p_cyl_top_right(:, k)
    end do

    ! Close data files
    close(unit=ud_cyl_top_left, iostat=IFAIL, status='keep')
    close(unit=ud_cyl_top_right, iostat=IFAIL, status='keep')

    deallocate(data_cyl_top_left)
    deallocate(data_cyl_top_right)
    deallocate(p_cyl_top_left)
    deallocate(p_cyl_top_right)
    deallocate(len_cyl_top_left)
    deallocate(len_cyl_top_right)
end subroutine Calc_E_top_cyl

end module mod_cylindrical_tip