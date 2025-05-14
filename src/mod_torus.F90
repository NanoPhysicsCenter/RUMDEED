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

    ! Parameters for the torus
    double precision :: R = 10.0d0
    double precision :: rho = 1.0d0

    ! ----------------------------------------------------------------------------
    ! Variables for the Metropolis-Hastings algorithm
    integer, parameter                 :: N_MH_step = 55 ! Number of steps to do in the MH algorithm

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
    !call Create_KD_Tree()
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
    character (len=*), parameter :: filename_meshdata = "Cyl_mesh_data.txt"
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
        phi = atan2(pos(2), pos(1)) ! See notes
        x_c = R*cos(phi)
        y_c = 0.0d0
        z_c = R*sin(phi)

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
    !field_E = field_E * 0.125d1 ! Debug to make field larger
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
    if (mod(step, steps/10) == 0) then
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
    end if
  
    call Do_Surface_Integration_simple(N_sup)
    N_round = Rand_Poisson(N_sup)
  
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
          print *, 'Warning: F > 0.0d0 (Do_Field_Emission_Cylinder)'
          print *, 'F_norm = ', F_norm
          print *, 'F = ', F
          print *, 'par_pos = ', par_pos
          print *, 'sec = ', sec
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

    ! Generate random angles
    phi = 1.0d0 * pi * rand()
    theta = 2.0d0 * pi * rand()

    pos_torus(1) = rho
    pos_torus(2) = phi
    pos_torus(3) = theta

    ! Calculate the position on the surface of the torus
    pos(1) = (R + rho * cos(theta)) * cos(phi)
    pos(2) = rho * sin(theta)
    pos(3) = (R + rho * cos(theta)) * sin(phi)
end subroutine Get_Random_Surface_Pos

function Field_Normal(pos, pos_torus, F)
    double precision                             :: Field_Normal
    double precision, dimension(1:3), intent(in) :: pos ! Position to calculate the field normal at
    double precision, dimension(1:3), intent(in) :: pos_torus ! Position of in toroidal coordinates
    double precision, dimension(1:3), intent(in) :: F ! Electric field at the point
    double precision, dimension(1:3)             :: unit_vector ! Unit vector from the surface at pos

    double precision :: phi, theta

    ! Calculate the angles from the position
    phi = pos_torus(2)
    theta = pos_torus(3)

    ! Get a unit vector from the surface at pos
    ! Check this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    unit_vector(1) = (R + rho * cos(theta)) * cos(phi)
    unit_vector(2) = rho * sin(theta)
    unit_vector(3) = (R + rho * cos(theta)) * sin(phi)

    ! Normalize the unit vector
    unit_vector = unit_vector / sqrt(sum(unit_vector**2))

    ! Calculate the normal vector to the surface
    Field_Normal = dot_product(unit_vector, F)

    ! Normalize the normal vector
    !Field_Normal = Field_Normal / sqrt(sum(Field_Normal**2))
end function Field_Normal

!-------------------------------------------!
! Metropolis-Hastings algorithm for the torus
subroutine Metropolis_Hastings_Torus(N_MH_step, F_out, F_norm_out, pos_xyz_out, n_vec)
    integer, intent(in) :: N_MH_step
    double precision, dimension(1:3), intent(out) :: F_out
    double precision, intent(out) :: F_norm_out
    double precision, dimension(1:3), intent(out) :: pos_xyz_out
    double precision, dimension(1:3), intent(out) :: n_vec

    double precision, dimension(1:3) :: F_cur, F_new
    double precision, dimension(1:3) :: pos_cur, pos_new, pos_torus
    double precision, dimension(1:3) :: n_vec_cur, n_vec_new
    double precision :: F_norm, F_norm_new, F_norm_cur

    integer :: IFAIL, k

    ! Local variables
    integer :: i

    ! Initialize the output variables
    F_out = 0.0d0
    F_norm_out = 0.0d0
    pos_xyz_out = 0.0d0

    ! Get an initial position for the particle
    k = 0
    do
        call Get_Random_Surface_Pos(pos_cur, pos_torus)

        F_cur = field_E_Torus(pos_cur)
        F_norm_cur = Field_Normal(pos_cur, pos_torus, F_cur)

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


end subroutine Metropolis_Hastings_Torus


subroutine Do_Surface_Integration_simple(N_sup)
    double precision, intent(out) :: N_sup

    N_sup = 0.0d0
end subroutine Do_Surface_Integration_simple

end module mod_torus