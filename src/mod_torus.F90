!--------------------------------------------!
! Module for emission from a torus           !
! Kristinn Torfason                          !
! 11.03.25                                   !
!--------------------------------------------!

Module mod_torus
    use mod_global
#if defined(_OPENMP)
    use omp_lib
#endif
    use mod_verlet
    use mod_pair
    !use mod_work_function
    use kdtree2_precision_module
    use kdtree2_module
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan

    implicit none
    PRIVATE
    PUBLIC :: Init_Torus, Clean_Up_Torus

    ! ----------------------------------------------------------------------------
    ! KD-tree variables
    type(kdtree2), pointer :: kd_tree => null() ! KD-tree for the mesh
    double precision, dimension(:, :), allocatable :: kd_mesh ! Mesh for the KD-tree
    double precision, dimension(:, :), allocatable :: kd_data ! Data for the KD-tree

    double precision :: time_step_div_q0

    ! ----------------------------------------------------------------------------
    ! Variables
    integer, dimension(:), allocatable :: nrEmitted_emitters

    double precision, dimension(1:3)   :: F_avg = 0.0d0

    ! Parameters for the torus
    double precision :: rho = 0.012d0 * 1.0d-3 ! mm -> m
    double precision :: R_y = 0.05d0 * 1.0d-3 + 0.012d0*1.0d-3 ! mm -> m
    double precision :: R_z = 0.1d0 * 1.0d-3 + 0.012d0*1.0d-3 ! mm -> m

    ! ----------------------------------------------------------------------------
    ! Variables for the Metropolis-Hastings algorithm
    integer, parameter                 :: N_MH_step = 55 ! Number of steps to do in the MH algorithm


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

contains
!-------------------------------------------!
! Initialize the cylinder emission
! TODO: Change this for the new system
subroutine Init_Torus()
    ! Local variables
    !double precision :: int_res
    !double precision, dimension(1:3) :: E_test, pos_test, vel_test
    !integer :: IFAIL

    ! Allocate the number of emitters
    allocate(nrEmitted_emitters(1:nrEmit))

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_Torus

    ! Function for the electric field in the system
    ptr_field_E => field_E_Torus

    ! The function that does the emission
    ptr_Do_Emission => Do_Field_Emission_Torus_simple

    ! The function to do image charge effects
    ptr_Image_Charge_effect => Image_Charge_Torus ! Force_Image_charges for the cylinder
    !ptr_Image_Charge_effect => Force_Image_charges_v2

    !call Read_work_function()

    ! Parameters used in the module
    time_step_div_q0 = time_step / q_0

    ! Create the KD-tree
    call Create_KD_Tree()

    ! Output field along a curve on top of the looped CNT
    call Calc_E_Along_Top()
end subroutine Init_Torus

!-------------------------------------------!
! Clean up the cylinder emission
subroutine Clean_Up_Torus()
    deallocate(nrEmitted_emitters)

    !close(unit=ud_cyl_debug)
end subroutine Clean_Up_Torus

!-------------------------------------------!
! Create the K-dimensional tree using the kdtree2 module
subroutine Create_KD_Tree()
    character (len=*), parameter :: filename_meshdata = "Torus_mesh_data.txt"
    character (len=*), parameter :: filename_imagedata = "image_charge_data.txt"
    integer                      :: IFAIL
    integer                      :: ud_meshdata, ud_imagedata
    integer                      :: n_points, k
    double precision             :: max_z

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

    ! Rewind to the start of the file
    rewind(ud_meshdata)

    ! Allocate the arrays to hold the file
    allocate(kd_mesh(1:3, 1:n_points))
    allocate(kd_data(1:3, 1:n_points))

    ! Read the file into the array and scale the coordinates from mm to m
    do k = 1, n_points
      read(ud_meshdata, *) kd_mesh(1, k), kd_mesh(2, k), kd_mesh(3, k), kd_data(1, k), kd_data(2, k), kd_data(3, k)
    
      ! Convert from mm to m
      kd_mesh(:, k) = kd_mesh(:, k) * 1.0d-3
    end do

    ! Find the maximum z-coordinate
    max_z = maxval(kd_mesh(3, :))

    ! Compare with the height of the simulation box
    if (abs(max_z - box_dim(3))/length_scale > 1.0d-6) then
      print *, 'RUMDEED: ERROR The maximum z-coordinate in the mesh data is not equal to the height of the simulation box'
      print *, 'max_z = ', max_z, ' m'
      print *, 'box_dim(3) = ', box_dim(3), ' m'
      stop
    end if

    ! Close the file
    close(unit=ud_meshdata, iostat=IFAIL, status='keep')

    ! Create the kd-tree
    kd_tree => kdtree2_create(kd_mesh, sort=.false., rearrange=.true.)

    ! Check if the kd-tree was created successfully
    if (.not. associated(kd_tree)) then
      print *, 'RUMDEED: ERROR The kd-tree was not created successfully'
      stop
    end if

    ! To do: Image charge data
end subroutine Create_KD_Tree

!-------------------------------------------!
! Output field along a curve on top of the looped cnt
subroutine Calc_E_Along_Top()
  implicit none
  integer :: k
  integer, parameter :: N_p = 10000

  real(kdkind), dimension(1:3) :: p, E_vec, E_vec_image
  double precision :: phi, rho, theta
  integer :: ud_cyl_around, IFAIL
  !character (len=*), parameter :: filename_cyl_circle="cyl_E_circle.dt"
  character (len=100) :: filename_cyl_circle

  ! Data arrays
  double precision, dimension(:, :, :), allocatable :: data_cyl_around
  double precision, dimension(:, :), allocatable :: p_cyl_around
  double precision, dimension(:), allocatable    :: len_cyl_around

  !print *, 'Calc_E_circle_cyl started'

  allocate(data_cyl_around(1:3, 1:3, 1:N_p))
  allocate(len_cyl_around(1:N_p))
  allocate(p_cyl_around(1:3, 1:N_p))

  ! Calculate the electric field around the cylinder
  !print *, 'Calculating the electric field around the cylinder'
  do k = 1, N_p
    ! Right corner
    phi = (k-1)*1.0d0*pi/(N_p-1)
    theta = 0.0d0
    p(1) = 0.0d0
    p(2) = (R_y + rho*cos(theta))*cos(phi)
    p(3) = (R_z + rho*cos(theta))*sin(phi)

    p_cyl_around(:, k) = p(:)
    len_cyl_around(k) = phi

    E_vec = field_E_Torus(p) ! Vacuum field 
    E_vec_image = Calc_Field_at(p) ! Total field
    ! Store data in array
    data_cyl_around(1, :, k) = E_vec_image ! Total field
    data_cyl_around(2, :, k) = E_vec(:) ! Vacuum field
    data_cyl_around(3, :, k) = E_vec_image(:) - E_vec ! Electric field due to the image charges and electrons
  end do

  ! Open data file
  !print *, 'Opening data file'

  filename_cyl_circle = "out/loop_E_top.dt"

  open(newunit=ud_cyl_around, iostat=IFAIL, file=filename_cyl_circle, status='REPLACE', action='WRITE')
  if (IFAIL /= 0) then
    print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file ', filename_cyl_circle
    stop
  end if

  ! Write data to file
  !print *, 'Writing data to file'
  do k = 1, N_p
    write(ud_cyl_around, *) data_cyl_around(1, :, k), data_cyl_around(2, :, k), data_cyl_around(3, :, k), &
                            len_cyl_around(k), p_cyl_around(:, k)
  end do

  ! Close data file
  close(unit=ud_cyl_around, iostat=IFAIL, status='keep')

  deallocate(data_cyl_around)
  deallocate(p_cyl_around)
  deallocate(len_cyl_around)

  !print *, 'Calc_E_circle_cyl finished'
end subroutine Calc_E_Along_Top

!-------------------------------------------!
! Check the boundary conditions for the torus
! Coordinates are in the form of (x, y, z)
! x = (R + rho*cos(theta))*cos(phi)
! y = rho*sin(theta)
! z = (R + rho*cos(theta))*sin(phi)
! where R is the radius of the torus, rho is the distance from the center of the torus
! and theta and phi are the angles in the toroidal coordinates.
subroutine Check_Boundary_Torus(i)
    integer, intent(in) :: i
    integer             :: sec

    ! Local variables
    double precision, dimension(1:3) :: pos
    double precision                 :: phi, x_c, y_c, z_c

    pos = particles_cur_pos(:, i)

    ! Check if we are below or above the planes at z = 0 and z = d
    if (pos(3) <= 0.0d0) then
        call Mark_Particles_Remove(i, remove_bot)
    else if (pos(3) > box_dim(3)) then
        call Mark_Particles_Remove(i, remove_top)
    else
        ! Check if we are inside the torus

        ! Find the closest point on the torus center curve
        phi = atan2(pos(2), pos(1)) ! See notes (Check this! after flip)
        x_c = 0.0d0
        y_c = R_y*cos(phi)
        z_c = R_z*sin(phi)

        ! Calculate the distance from the particle to the center of the torus
        d = sqrt((pos(1) - x_c)**2 + (pos(2) - y_c)**2 + (pos(3) - z_c)**2)
        if (d <= rho) then
            ! The particle is outside the torus
            call Mark_Particles_Remove(i, remove_bot)
        end if
    end if


    ! To do
    ! Check the boundary conditions for the particle
    !if (Check_Boundary_Cylinder_pos(particles_cur_pos(:, i), sec) .eqv. .true.) then
      ! Remove the particle from the system
    !  call Mark_Particles_Remove(i, sec)
    !end if
  
end subroutine Check_Boundary_Torus

!-------------------------------------------!
! Calculate the electric field at a point in the system
function field_E_Torus(pos) result(field_E)
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
    field_E = field_E * 3.5d0 ! Debug to make field larger
end function field_E_Torus

subroutine Do_Field_Emission_Torus_simple(step)
    integer, intent(in) :: step
    integer             :: i, IFAIL, per
    integer, parameter  :: emit = 1
  
    ! Integration
    double precision                 :: N_sup
    integer                          :: N_round
  
    ! Emission variables
    double precision                 :: D_f, Df_avg, F_norm, rnd
    double precision, dimension(1:3) :: F, n_vec
    integer                          :: s, sec, sec_b
    integer                          :: nrElecEmit
    double precision, dimension(1:3) :: par_pos, par_vel, old_pos
    !logical                          :: to_pause = .false.
  
    nrEmitted_emitters = 0
  
    ! Check if the step is 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90% or 100% of the total steps
    !if (mod(step, steps/10) == 0) then
      !per = nint(dble(step)/dble(steps)*100.0d0) ! Calculate the percentage of the simulation
      !print *, 'Step = ', step
      !print *, 'Percentage = ', per, '%'
  
      ! Calculate the field at the surface
      !print *, 'Calling Calc_E_edge_cyl'
      !call Calc_E_edge_cyl(per)
      !print *, 'Calling Calc_E_side_cyl'
      !call Calc_E_top_cyl(per)
      !print *, 'Calling Calc_E_top_cyl'
      !call Calc_E_corner_cyl(per)
  
      !call Calc_E_cyl(per)
      !call Calc_E_circle_cyl(per)
    !end if
  
    call Do_Surface_Integration_simple(N_sup)
    N_round = Rand_Poisson(N_sup)
    !print *, 'N_round = ', N_round
    !pause

    nrElecEmit = 0
    nrEmitted_emitters(emit) = 0
  
    Df_avg = 0.0d0
    i = 0
    ! Loop over the electrons to be emitted.
    do s = 1, N_round
        ! Get the position of the particle
        ! ndim_in, F_out, F_norm_out, pos_xyz_out, sec_out
        call Metropolis_Hastings_Torus(N_MH_step, F, F_norm, par_pos, n_vec)
  
        ! Check if the field is favorable for emission or not
        if (F_norm >= 0.0d0) then
          D_f = -huge(1.0d0)
          print *, 'Warning: F > 0.0d0 (Do_Field_Emission_Torus)'
          print *, 'F_norm = ', F_norm
          print *, 'F = ', F
          print *, 'par_pos = ', par_pos
          !print *, 'sec = ', sec
          stop
        else
  
          par_pos = par_pos + n_vec*length_scale
          call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1, sec)
  
          nrElecEmit = nrElecEmit + 1
          nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
        end if
      end do
  
      ! Df_avg = Df_avg / dble(i)
      ! !print *, 'df_avg = ', Df_avg
      ! if (Df_avg > 1.0d-4) then
      !   print *, 'RUMDEED: Df_avg > 1.0d-4 (Do_Field_Emission_Cylinder)'
      !   print *, step
      !   write_position_file = .true.
      !   call Write_Position(0)
      !   cought_stop_signal = .true.
      ! end if
  
    write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, nrElecEmit, nrElec, &
                                                       & (nrEmitted_emitters(i), i = 1, nrEmit)
end subroutine Do_Field_Emission_Torus_simple

!-------------------------------------------!
! Image charge effect
function Image_Charge_Torus(pos_1, pos_2)
    double precision, dimension(1:3)             :: Image_Charge_Torus ! Image charge effect
    double precision, intent(in), dimension(1:3) :: pos_1 ! Position of the particle we are calculating the force/acceleration on
    double precision, intent(in), dimension(1:3) :: pos_2 ! Position of the particle that is acting on the particle at pos_1

    Image_Charge_Torus = 0.0d0
end function Image_Charge_Torus

!-------------------------------------------!
! Generate a random position on the surface of the torus
subroutine Get_Random_Surface_Pos(pos, pos_torus)
    double precision, dimension(1:3), intent(out) :: pos
    double precision, dimension(1:3), intent(out)  :: pos_torus ! Toroidal coordinates
    double precision :: phi, theta

    double precision, dimension(1:2) :: rnd

    ! Get random numbers
    call random_number(rnd)

    ! Generate random angles
    phi = 1.0d0 * pi * rnd(1)
    theta = 2.0d0 * pi * rnd(2)

    pos_torus(1) = rho
    pos_torus(2) = phi
    pos_torus(3) = theta

    ! Calculate the xyz position on the surface of the torus
    pos = Convert_to_xyz(pos_torus)
end subroutine Get_Random_Surface_Pos

function Field_Normal(pos, pos_torus, F)
    double precision                             :: Field_Normal
    double precision, dimension(1:3), intent(in) :: pos ! Position to calculate the field normal at
    double precision, dimension(1:3), intent(in) :: pos_torus ! Position of in toroidal coordinates
    double precision, dimension(1:3), intent(in) :: F ! Electric field at the point
    double precision, dimension(1:3)             :: unit_vector ! Unit vector from the surface at pos

    double precision :: phi, theta

    ! Calculate the angles from the position
    !phi = pos_torus(2)
    !theta = pos_torus(3)

    ! Get a unit vector from the surface at pos
    unit_vector = Get_Unit_Vector(pos_torus)

    ! Calculate the normal vector to the surface
    Field_Normal = dot_product(unit_vector, F)
end function Field_Normal

!-------------------------------------------!
! Calculate the unit vector from the surface at pos
! Check this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Get_Unit_Vector(pos_torus)
    double precision, dimension(1:3)             :: Get_Unit_Vector
    double precision, dimension(1:3), intent(in) :: pos_torus ! Position of in toroidal coordinates

    double precision :: phi, theta

    ! Calculate the angles from the position
    ! rho is fixed, should be rho = pos_torus(1)
    !phi = pos_torus(2)
    !theta = pos_torus(3)

    ! Get a unit vector from the surface at pos
    ! Unit vector is [ sin(theta), cos(theta)*cos(phi), cos(theta)*sin(phi) ]
    Get_Unit_Vector(1) = sin(pos_torus(3))
    Get_Unit_Vector(2) = cos(pos_torus(3))*cos(pos_torus(2))
    Get_Unit_Vector(3) = cos(pos_torus(3))*sin(pos_torus(2))
end function Get_Unit_Vector

!-------------------------------------------!
! Convert from toroidal coordinates to Cartesian coordinates
function Convert_to_xyz(pos_torus)
    double precision, dimension(1:3)             :: Convert_to_xyz
    double precision, dimension(1:3), intent(in) :: pos_torus ! Position of in toroidal coordinates

    ! Calculate the position on the surface of the torus
    Convert_to_xyz(1) = rho * sin(pos_torus(3))
    Convert_to_xyz(2) = (R_y + rho * cos(pos_torus(3))) * cos(pos_torus(2))
    Convert_to_xyz(3) = (R_z + rho * cos(pos_torus(3))) * sin(pos_torus(2))
end function Convert_to_xyz

!-------------------------------------------!
! Convert from Cartesian coordinates to toroidal coordinates
! Check this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Convert_to_torus(pos)
    double precision, dimension(1:3)             :: Convert_to_torus
    double precision, dimension(1:3), intent(in) :: pos ! Position to calculate the toroidal coordinates at

    ! Calculate the toroidal coordinates
    Convert_to_torus(1) = rho
    Convert_to_torus(2) = atan2(pos(1), pos(2))
    Convert_to_torus(3) = atan2(pos(3), pos(2))
end function Convert_to_torus

!-------------------------------------------!
! Metropolis-Hastings algorithm for the torus
subroutine Metropolis_Hastings_Torus(ndim_in, F_out, F_norm_out, pos_xyz_out, n_vec)
    integer, intent(in)                           :: ndim_in
    double precision, dimension(1:3), intent(out) :: F_out
    double precision, intent(out)                 :: F_norm_out
    double precision, dimension(1:3), intent(out) :: pos_xyz_out
    double precision, dimension(1:3), intent(out) :: n_vec

    ! Local variables
    double precision, dimension(1:3) :: F_cur, F_new
    double precision, dimension(1:3) :: pos_cur, pos_new
    double precision, dimension(1:3) :: pos_cur_torus, pos_new_torus
    !double precision, dimension(1:3) :: n_vec_cur, n_vec_new
    double precision :: F_norm, F_norm_new, F_norm_cur
    double precision :: alpha, rnd

    integer :: IFAIL, k, i
    integer :: accepted, rejected

    ! Initialize the output variables
    F_out = 0.0d0
    F_norm_out = 0.0d0
    pos_xyz_out = 0.0d0

    ! Get an initial position for the particle
    k = 0
    do
        call Get_Random_Surface_Pos(pos_cur, pos_cur_torus)

        F_cur = Calc_Field_at(pos_cur)
        F_norm_cur = Field_Normal(pos_cur, pos_cur_torus, F_cur)

        if (F_norm_cur < 0.0d0) then
            ! The field is favorable for emission
            exit
        end if

        if (k > 100000) then
            print *, 'RUMDEED: ERROR: Unable to find a position for the particle (Metropolis_Hastings_Torus)'
            stop
        else
            k = k + 1
        end if
    end do
    
    ! Do the Metropolis-Hastings algorithm
    do i = 1, ndim_in
        call Jump_MH(pos_cur, pos_cur_torus, pos_new, pos_new_torus)

        alpha = Get_Jump_Probability(pos_cur, pos_cur_torus, F_norm_cur, F_norm_new)

        ! Accept or reject the jump
        call random_number(rnd)
        if (rnd < alpha) then
            ! Accept the jump
            accepted = accepted + 1

            ! Update the position
            pos_cur = pos_new
            pos_cur_torus = pos_new_torus

            ! Update the electric field
            F_cur = F_new
            F_norm_cur = F_norm_new
            
            !write (unit=ud_mh, iostat=IFAIL) cur_time, pos_xyz_cur/length_scale, sec_cur, pos_cur
        else
            ! Reject the jump
            rejected = rejected + 1
        end if
    end do

    ! Set the output variables
    pos_xyz_out = pos_cur
    n_vec = Get_Unit_Vector(pos_cur_torus)
    F_out = F_cur
    F_norm_out = F_norm_cur
end subroutine Metropolis_Hastings_Torus

function Get_Jump_Probability(pos, pos_torus, F_norm_cur, F_norm_new) result(alpha)
    double precision, dimension(1:3), intent(in) :: pos ! Position to calculate the jump probability at
    double precision, dimension(1:3), intent(in) :: pos_torus ! Position of in toroidal coordinates
    double precision, intent(in)                 :: F_norm_cur ! Electric field at the current position
    double precision, intent(out)                 :: F_norm_new ! Electric field at the new position
    double precision :: alpha

    ! Local variables
    double precision, dimension(1:3) :: F_new

    ! Calculate the electric field at the new position
    F_new = Calc_Field_at(pos)
    F_norm_new = Field_Normal(pos, pos_torus, F_new)

    ! Calculate the jump probability
    ! Calculate the jump probability
    if (F_norm_new >= 0.0d0) then ! If the field is positive, then the jump probability is zero
      alpha = 0.0d0
    else
      alpha = F_norm_new / F_norm_cur
    end if

end function Get_Jump_Probability

!-------------------------------------------!
! Jump in the Metropolis-Hastings algorithm
subroutine Jump_MH(pos_cur, pos_torus, pos_new, pos_new_torus)
    double precision, dimension(1:3), intent(in)  :: pos_cur
    double precision, dimension(1:3), intent(in)  :: pos_torus
    double precision, dimension(1:3), intent(out) :: pos_new
    double precision, dimension(1:3), intent(out) :: pos_new_torus

    ! Local
    double precision, dimension(1:2) :: std_xy

    ! Generate a new position on the surface from the current position
    std_xy(1) = 1.0d0*pi*0.05d0
    std_xy(2) = 2.0d0*pi*0.05d0
    pos_new_torus(2:3) = box_muller(pos_torus(2:3), std_xy) ! Jump in x and y

    ! Check the angles are within the limits and rebound if necessary
    ! phi is between 0 and pi
    if (pos_new_torus(2) < 0.0d0) then
      pos_new_torus(2) = -pos_new_torus(2)
    else if (pos_new_torus(2) >= pi) then
      pos_new_torus(2) = 2.0d0*pi - pos_new_torus(2) ! Reflect back, d = P - L, P' = L - d = 2L - P
    end if
    ! theta is between 0 and 2*pi
    if (pos_new_torus(3) < 0.0d0) then
      pos_new_torus(3) = -pos_new_torus(3)
    else if (pos_new_torus(3) >= 2.0d0*pi) then
      pos_new_torus(3) = 4.0d0*pi - pos_new_torus(3)
    end if

    pos_new_torus(1) = rho ! rho is fixed

    ! Calculate the x and y coordinates
    pos_new = Convert_to_xyz(pos_new_torus)

    ! Calculate the normal vector to the surface at the new position
    !n_vec_new = Get_Unit_Vector(pos_new_torus)
end subroutine Jump_MH


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
  !integer                                      :: emit = 1
  double precision                             :: l

  if (image_charge .eqv. .true.) then
    l = l_const * (-1.0d0*F) / w_theta_pos_tip(pos)**2 ! l = y^2, y = 3.79E-4 * sqrt(F_old) / w_theta
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
  !integer                                      :: emit = 1
  double precision                             :: l

  if (image_charge .eqv. .true.) then
    l = l_const * (-1.0d0*F) / w_theta_pos_tip(pos)**2 ! l = y^2, y = 3.79E-4 * sqrt(F_old) / w_theta
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
  !integer                                      :: emit = 1

  !print *, 'F = ', F
  !print *, 'pos = ', pos
  !print *, 'v_y = ', v_y(F, pos)
  !print *, 'w_theta = ', w_theta_pos_tip(pos)
  !print *, b_FN * (sqrt(w_theta_pos_tip(pos)))**3 * v_y(F, pos) / abs(F)
  !print *, 'Escape_prob = ', exp(b_FN * (sqrt(w_theta_pos_tip(pos)))**3 * v_y(F, pos) / abs(F))
  !print *, ''
  Escape_Prob = exp(b_FN * (sqrt(w_theta_pos_tip(pos)))**3 * v_y(F, pos) / abs(F))
  if (Escape_Prob > 1.0d0) then
    print *, 'Escape_prob is larger than 1.0'
    print *, 'Escape_prob = ', Escape_Prob
    print *, ''
  else if (Escape_Prob < 0.0d0) then
    print *, 'Escape_prob is smaller than 0.0'
    print *, 'Escape_prob = ', Escape_Prob
    print *, ''
  end if
end function Escape_Prob


!-----------------------------------------------------------------------------
! A simple function that calculates
! A_FN/(t**2(l)*w_theta(x,y)) F**2(x,y)
! pos: Position to calculate the function
! F: The z-component of the field at par_pos, it should be F < 0.0d0.
double precision function Elec_Supply_V2(F, pos)
  double precision, dimension(1:3), intent(in) :: pos
  double precision,                 intent(in) :: F
  !integer                                      :: emit = 1

  Elec_Supply_V2 = time_step_div_q0 * a_FN/(t_y(F, pos)**2*w_theta_pos_tip(pos)) * F**2
end function Elec_Supply_V2

double precision function w_theta_pos_tip(pos)
  double precision, dimension(1:3), intent(in) :: pos

  w_theta_pos_tip = 4.65d0
end function w_theta_pos_tip


!-------------------------------------------!
! Do the surface integration for the torus
subroutine Do_Surface_Integration_simple(N_sup)
    double precision, intent(out) :: N_sup ! Number of electrons

    SELECT CASE (cuba_method)
    case(cuba_method_suave)
      !call Do_Cuba_Suave(N_sup)
      print '(a)', 'RUMDEED: ERROR SIMPLE METHOD NOT IMPLEMENTED FOR SUAVE'
      stop
    case(cuba_method_divonne)
      call Do_Cuba_Divonne_simple(N_sup)
    case default
      print '(a)', 'RUMDEED: ERROR UNKNOWN INTEGRATION METHOD'
      print *, cuba_method
      stop
    end select
end subroutine Do_Surface_Integration_simple

  subroutine Do_Cuba_Divonne_simple(N_sup)
    implicit none
    double precision, intent(out) :: N_sup
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
    ! Set seed to non-zero value based on the current time
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
    integer :: key2 = -1 !<in>, determines sampling in the final integration phase:
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

    ! Set seed
    !seed = my_seed(1)

    !-----------------------------------------------------------------------------
    ! Integrate over the surface
   
    call divonne(ndim, ncomp, integrand_cuba_torus_simple, userdata, nvec,&
                cuba_epsrel, cuba_epsabs, flags, seed, cuba_mineval, cuba_maxeval,&
                key1, key2, key3, maxpass,&
                border, maxchisq, mindeviation,&
                ngiven, ldxgiven, 0, nextra, 0,&
                statefile, spin,&
                nregions, neval, fail, integral, error, prob)
  
    if (fail /= 0) then
      if (abs(error(1) - cuba_epsabs) > 1.0d-2) then
        print '(a)', 'RUMDEED: WARNING Cuba did not return 0'
        print *, 'Top integration'
        print *, 'time_step = ', time_step
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
    N_sup = integral(1)
    !print *, 'N_sup = ', N_sup

    ! Write the output variables of the integration to a file
    write(ud_integrand, '(i8, tr2, i8, tr2, i4, tr2, ES12.4, tr2, ES12.4, tr2, ES12.4)', iostat=IFAIL) &
                            & nregions, neval, fail, integral(1), error(1), prob(1)

    ! Finish calculating the average field
    F_avg = F_avg / neval

  end subroutine Do_Cuba_Divonne_simple

  ! ----------------------------------------------------------------------------
  ! The integration function for the Cuba library
  !
  integer function integrand_cuba_torus_simple(ndim, xx, ncomp, ff, userdata)
    ! Input / output variables
    integer, intent(in) :: ndim ! Number of dimensions (Should be 2)
    integer, intent(in) :: ncomp ! Number of vector-components in the integrand (Always 1 here)
    integer, intent(in) :: userdata ! Additional data passed to the integral function (In our case the number of the emitter)
    double precision, intent(in), dimension(1:ndim)   :: xx ! Integration points, between 0 and 1
    double precision, intent(out), dimension(1:ncomp) :: ff ! Results of the integrand function

    ! Variables used for calculations
    double precision, dimension(1:3) :: par_pos, par_pos_org, field
    double precision                 :: field_norm
    double precision                 :: h_theta, h_phi

    ! Variables used for the transformation from the hypercube to the surface
    double precision, dimension(1:2) :: upper, lower, range
    double precision                 :: jacobian ! Jacobian of the transformation from the hypercube to the surface

    ! Other variables
    integer :: IFAIL

    ! Use toroidal coordinates for the corner surface (r, phi, theta)
    ! phi [0, pi]
    ! theta [0, 2*pi]
    upper(1:2) = (/pi, 2.0d0*pi/)
    lower(1:2) = (/0.0d0, 0.0d0/)

    ! Range
    range(1:2) = upper(1:2) - lower(1:2)

    ! rho, phi and theta coordinates of the point on the surface
    par_pos_org(1) = rho ! Fixed distance from the center of the torus
    par_pos_org(2:3) = lower(1:2) + xx(1:2)*range(1:2)

    ! Convert to cartesian coordinates
    par_pos = Convert_to_xyz(par_pos_org)
    !par_pos(1) = rho*sin(par_pos_org(3))
    !par_pos(2) = (R + rho*cos(par_pos_org(3)))*cos(par_pos_org(2))
    !par_pos(3) = (R + rho*cos(par_pos_org(3)))*sin(par_pos_org(2))

    ! Jacobian of the transformation from the hypercube (Check this!)
    h_theta = rho
    h_phi = sqrt((R_y + rho*cos(par_pos_org(3)))**2*sin(par_pos_org(2))**2 &
          & + (R_z + rho*cos(par_pos_org(3)))**2*cos(par_pos_org(2))**2)
    jacobian = h_theta * h_phi * (range(1)*range(2)) ! Jacobian of the transformation from the hypercube to the surface

    !ff(1) = 1.0d0 ! Integrand results debug

    integrand_cuba_torus_simple = 0 ! Return value to Cuba, 0 = success

    ! Calculate the electric field on the surface
    field = Calc_Field_at(par_pos)
    field_norm = Field_normal(par_pos, par_pos_org, field)

    ! Add to the average field
    F_avg = F_avg + field

    ! Check if the field is favorable for emission
    if (field_norm < 0.0d0) then
      ! The field is favorable for emission
      ! Calculate the electron supply at this point
      ff(1) = Elec_Supply_V2(field_norm, par_pos)*Escape_Prob(field_norm, par_pos)
      !print *, 'Electron supply = ', Elec_Supply_V2(field_norm, par_pos)
      !print *, 'Escape probability = ', Escape_Prob(field_norm, par_pos)
      !print *, 'ff(1) = ', ff(1)
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
    !ff(1) = 1.0d0*jacobian ! Area test
    !if (userdata == sec_corner) then
    !  print *, 'ff(1) = ', ff(1)
    !end if

    !write (ud_field, "(*(E16.8, tr2))", iostat=IFAIL) &
    !                                  par_pos(1), par_pos(2), par_pos(3), &
    !                                  field(1), field(2), field(3), &
    !                                  field_norm
  end function integrand_cuba_torus_simple

end module mod_torus