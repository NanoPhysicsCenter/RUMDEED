!-------------------------------------------!
! Module for photo emission                 !
! Kristinn Torfason                         !
! 05.04.13                                  !
!-------------------------------------------!

Module mod_photo_emission
  use mod_global
  use mod_verlet
  use mod_pair
  implicit none

  ! ----------------------------------------------------------------------------
  ! Variables
  integer, dimension(:), allocatable          :: nrEmitted_emitters
  integer                                     :: posInit
  logical                                     :: EmitGauss = .false.
  integer                                     :: maxElecEmit = -1
  integer                                     :: nrEmitted

  integer                                     :: ud_gauss

  ! ----------------------------------------------------------------------------
  ! Parameters
  integer, parameter                          :: MAX_EMISSION_TRY = 1000

contains
  subroutine Init_Photo_Emission()
    ! Allocate the number of emitters
    allocate(nrEmitted_emitters(1:nrEmit))

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_ElecHole_Planar

    ! Function for the electric field in the system
    ptr_field_E => field_E_planar

    ! The function that does the emission
    ptr_Do_Emission => Do_Photo_Emission
  end subroutine Init_Photo_Emission

  subroutine Clean_Up_Photo_Emission()
    deallocate(nrEmitted_emitters)
  end subroutine Clean_Up_Photo_Emission

  subroutine Do_Photo_Emission(step)
    integer, intent(in) :: step
    integer             :: i, IFAIL
    double precision    :: cur_time

    posInit = 0
    nrEmitted_emitters = 0

    !if (step > 2) then
    !  return
    !end if

    ! Check if we are doing a Gaussian distributed emission
    ! and set the max number of electrons allowed to be emitted if we are
    if (EmitGauss .eqv. .TRUE.) then
      maxElecEmit = Gauss_Emission(step)
    end if

    !i = 1

    ! Loop through all of the emitters
    do i = 1, nrEmit

      ! Check the type of the emitter CIRCLE / RECTANGLE
      if (emitters_delay(i) < step) then
        select case (emitters_type(i))
        case (EMIT_CIRCLE)
          !print *, 'Doing Circle'
          call Do_Photo_Emission_Circle(step, i)
        case (EMIT_RECTANGLE)
          !print *, 'Doing Rectangle'
          call Do_Photo_Emission_Rectangle(step, i)
        case default
          print *, 'Vacuum: ERROR unknown emitter type!!'
          stop
          print *, emitters_type(i)
        end select
      end if
    end do

    cur_time = time_step * step / time_scale ! Scaled in units of time_scale
    write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, posInit, &
    & nrEmitted, nrElec, (nrEmitted_emitters(i), i = 1, nrEmit), maxElecEmit

  end subroutine Do_Photo_Emission

  subroutine Do_Photo_Emission_Circle(step, emit)
    integer, intent(in)              :: step, emit
    integer                          :: nrElecEmit, nrTry
    double precision, dimension(1:3) :: par_pos, par_vel, field
    double precision                 :: r_e, theta_e

    par_pos = 0.0d0
    nrTry = 0
    nrElecEmit = 0
    nrEmitted_emitters(emit) = 0

    do while (nrTry <= MAX_EMISSION_TRY)
      ! Check if we have reached the max number
      ! of electrons to be emitted
      if ((nrElecEmit >= maxElecEmit) .and. (maxElecEmit /= -1)) then
        exit
      end if

      if (nrElec == MAX_PARTICLES-1) then
        print *, 'WARNING: Reached maximum number of electrons!!!'
        exit
      end if

      CALL RANDOM_NUMBER(par_pos(1:2)) ! Gives a random number between [0,1]

      r_e = emitters_dim(1, emit) * par_pos(1)
      theta_e = 2.0d0*pi * par_pos(2)
      par_pos(1) = r_e * cos(theta_e) + emitters_pos(1, emit)
      par_pos(2) = r_e * sin(theta_e) + emitters_pos(2, emit)

      nrTry = nrTry + 1 ! Add one more try
      par_pos(3) = 0.0d0 * length_scale ! Check in plane
      field = Calc_Field_at(par_pos)

      ! Check if the field is favourable at this point
      if (field(3) < 0.0d0) then
        par_pos(3) = 1.0d0 * length_scale ! Place above plane
        par_vel = 0.0d0 ! Set the velocity

        ! Escape velocity from image charge partner
        par_vel(3) = ( q_0 / sqrt(8*pi*epsilon_0*epsilon_r*m_0*par_pos(3))) )*0.80d0
        
        ! Speed needed to reach over the gap spacing. Includes image charge partners behind the cathode and annode but not the electric field.
        !par_vel(3) = q_0/(4.0d0*sqrt(pi*epsilon_0*epsilon_r*m_0))*(d-2.0d0*par_pos(3))/sqrt(d*par_pos(3)*(d-par_pos(3)))
        
        call Add_Particle(par_pos, par_vel, species_elec, step, emit)

        nrElecEmit = nrElecEmit + 1
        nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
        nrTry = 0
      end if
    end do

    posInit = posInit + nrElecEmit
    nrEmitted = nrEmitted + nrElecEmit
  end subroutine Do_Photo_Emission_Circle

  subroutine Do_Photo_Emission_Rectangle(step, emit)
    integer, intent(in)              :: step, emit
    integer                          :: nrElecEmit, nrTry
    double precision, dimension(1:3) :: par_pos, par_vel, field

    par_pos = 0.0d0
    nrTry = 0
    nrElecEmit = 0
    nrEmitted_emitters(emit) = 0

    do while (nrTry <= MAX_EMISSION_TRY)

      ! Check if we have reached the max number
      ! of electrons to be emitted
      if ((nrElecEmit >= maxElecEmit) .and. (maxElecEmit /= -1)) then
        exit
      end if


      if (nrElec == MAX_PARTICLES-1) then
        print *, 'WARNING: Reached maximum number of electrons!!!'
        exit
      end if

      CALL RANDOM_NUMBER(par_pos(1:2)) ! Gives a random number [0,1]

      par_pos(1:2) = emitters_pos(1:2, emit) + emitters_dim(1:2, emit)*par_pos(1:2)
      !par_pos(2) = emitters_pos(2, emit) + emitters_dim(2, emit)*par_pos(2)

      nrTry = nrTry + 1
      par_pos(3) = 0.0d0 * length_scale !Check in plane
      field = Calc_Field_at(par_pos)

      if (field(3) < 0.0d0) then
        par_pos(3) = 1.0d0 * length_scale !Above plane
        field = Calc_Field_at(par_pos)

        if (field(3) < 0.0d0) then

          par_pos(3) = 1.0d0 * length_scale ! Place above plane
          par_vel = 0.0d0
          call Add_Particle(par_pos, par_vel, species_elec, step, emit)

          !print *, 'field = ', field
          !pause
          !call Add_Plane_Graph_emitt(par_pos, par_vel)

          nrElecEmit = nrElecEmit + 1
          nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
          nrTry = 0

        end if
      end if
    end do

    posInit = posInit + nrElecEmit
    nrEmitted = nrEmitted + nrElecEmit
  end subroutine Do_Photo_Emission_Rectangle

  ! ! Add the absorption event to the absorption graph.
  ! ! This graph is a histogram.
  ! subroutine Add_Plane_Graph_emitt(pos, vel)
  !   double precision, dimension(1:3), intent(in) :: pos, vel
  !   integer :: int_x, int_y
  !
  !   if ((pos(1) < int_end_x) .and. (pos(1) > int_start_x) .and. (pos(2) < int_end_y) .and. (pos(2) > int_start_y)) then
  !     int_x = ceiling((pos(1)+abs(int_start_x))/dx)
  !     int_y = ceiling((pos(2)+abs(int_start_y))/dy)
  !
  !     plane_graph_mesh_emitt(int_x, int_y) = plane_graph_mesh_emitt(int_x, int_y) + 1
  !   else
  !     print *, 'Vaccum: Particle outside plane range, abs'
  !     print *, 'x = ', pos(1)
  !     print *, 'y = ', pos(2)
  !     print *, 'z = ', pos(3)
  !
  !     print *, 'v_x = ', vel(1)
  !     print *, 'v_y = ', vel(2)
  !     print *, 'v_z = ', vel(3)
  !
  !     print *, ''
  !   end if
  ! end subroutine Add_Plane_Graph_emitt
  !
  ! ! Write out the absorption graph data for both the histogram
  ! ! plot and the kernel density estimation.
  ! subroutine Write_Plane_Graph_emitt()
  !   integer :: i, j, IFAIL
  !   double precision :: x, y
  !
  !   !do j = 1, res_y
  !   !    write (ud_plane_graph_matlab, "(*(E16.8, tr8))", iostat=IFAIL) (plane_graph_mesh(i, j), i = 1, res_x)
  !   !end do
  !
  !   !do j = 1, (res_y_kdens+1)
  !   !    write (ud_kdens_graph_matlab, "(*(E16.8, tr8))", iostat=IFAIL) (plane_graph_kdens(i, j), i = 1, (res_x_kdens+1))
  !   !end do
  !
  !   !rewind(ud_plane_graph)
  !
  !   do i = 1, res_x
  !     x = ((i - res_x*0.5d0) + 0.5d0)*dx
  !     x = x / length_scale
  !     do j = 1, res_y
  !       y = ((j - res_y*0.5d0) + 0.5d0)*dy
  !       y = y / length_scale
  !       !if (plane_graph_mesh(i, j) /= 0) then
  !         write (ud_plane_graph_emitt, "(i6, tr2, i6, tr2, E16.8, tr2, E16.8, tr2, i6)", iostat=IFAIL) &
  !         i, j, x, y, plane_graph_mesh_emitt(i, j)
  !       !end if
  !     end do
  !
  !     write (ud_plane_graph_emitt, *) ''
  !   end do
  ! end subroutine Write_Plane_Graph_emitt

  ! Gives a gaussian emission curve
  ! where step is the current time step
  ! returns the number of electrons allowed to be emitted in that time step
  integer function Gauss_Emission(step)
    integer, intent(in)         :: step ! Current time step
    integer                     :: IFAIL
    double precision, parameter :: sigma = 1000.0d0 ! Width / standard deviation
    double precision, parameter :: mu = 0.0d0 ! Center
    double precision, parameter :: A = 6.0d0 ! Height
    double precision, parameter :: b = 1.0d0/(2.0d0*pi*sigma**2)

    Gauss_Emission = IDNINT(  A * exp( -1.0d0*b*(step - mu)**2 )  )

    write (ud_gauss, "(i6, tr2, i6)", iostat=IFAIL) step, Gauss_Emission
  end function Gauss_Emission
end module mod_photo_emission
