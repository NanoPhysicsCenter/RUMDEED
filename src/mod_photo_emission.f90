!-------------------------------------------!
! Module for photo emission                 !
! Kristinn Torfason                         !
! 05.04.13                                  !
! Modifications                             !
! Hákon Örn Árnason                         !
! 04.12.20                                  !
!-------------------------------------------!

module mod_photo_emission
  use mod_global
  use mod_verlet
  !use mod_velocity
  use mod_pair
  use mod_work_function
  use mod_velocity
  implicit none

  ! ----------------------------------------------------------------------------
  ! Variables
  integer, dimension(:), allocatable          :: nrEmitted_emitters
  integer                                     :: posInit
  logical                                     :: EmitGauss                  = .FALSE.
  integer                                     :: maxElecEmit                = -1
  integer                                     :: nrEmitted
  integer                                     :: LASER_MODE, PHOTON_MODE, GAUSS_MODE
  double precision                            :: laser_energy               ! Photon energy in eV
  double precision                            :: laser_variation            ! Standard deviation of photon energy in eV
  double precision                            :: Gauss_pulse_center         ! mu
  double precision                            :: Gauss_pulse_width          ! sigma
  double precision                            :: Gauss_pulse_amplitude      ! A
  !integer                                     :: ud_gauss

  ! ----------------------------------------------------------------------------
  ! Parameters
  integer, parameter                          :: MAX_EMISSION_TRY = 100
  ! Type of laser input pulse
  integer, parameter                          :: LASER_FIXED      = 1  ! Fixed laser energy
  integer, parameter                          :: LASER_POISSON    = 2  ! Poisson distributed laser energy
  ! Type of photon initial velocity distribution
  integer, parameter                          :: PHOTON_ZERO      = 1  ! Zero inital velocity
  integer, parameter                          :: PHOTON_WFD       = 2  ! Workfunction dependant velocity
  ! Gauss pulse 
  integer, parameter                          :: GAUSS_ON         = 1  ! Gauss pulse on
  integer, parameter                          :: GAUSS_OFF        = 2  ! Gauss pulse off

  ! ----------------------------------------------------------------------------
  ! Spots
  integer                                        :: num_spots
  double precision, dimension(:, :), allocatable :: spot_pos       ! Position of spots
  double precision, dimension(:),    allocatable :: spot_radius_sq ! Radius squared of spots

contains

  subroutine Read_laser_parameters()
    integer :: ud_laser, IFAIL
    character(256) :: iomsg

    ! Open the file that contains information about the work function
    open(newunit=ud_laser, iostat=IFAIL, iomsg=iomsg, file='laser', &
      & status='OLD', form='FORMATTED', access='SEQUENTIAL', action='READ')
      
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file laser. ABORTING'
      print *, IFAIL
      print *, iomsg
      stop

    end if

    ! Read the type of laser pulse and photon velocity model to use
    read(unit=ud_laser, FMT=*) GAUSS_MODE, LASER_MODE, PHOTON_MODE

    select case (GAUSS_MODE)
      case (GAUSS_ON)
        ! Gaussian pulse
        print '(a)', 'RUMDEED: Using Gaussian pulse model'
        EmitGauss = .TRUE.

        select case (LASER_MODE)
          case (LASER_FIXED)    
            ptr_Get_Photo_Emission_Energy => Get_Fixed_Laser_Energy
            print '(a)', 'RUMDEED: Using fixed photon energy' 
            ! Read mean and std photon energy from the file
            read(unit=ud_laser, FMT=*) laser_energy, laser_variation
            ! Read Gaussian parameters from the file (mu, std and Amplitude) 
            read(unit=ud_laser, FMT=*) Gauss_pulse_center, Gauss_pulse_width, Gauss_pulse_amplitude

          case (LASER_POISSON)
            ptr_Get_Photo_Emission_Energy => Get_Laser_Energy
            print '(a)', 'RUMDEED: Using poisson distributed photon energy'
            ! Read mean and std photon energy from the file
            read(unit=ud_laser, FMT=*) laser_energy, laser_variation
            ! Read Gaussian parameters from the file (mu, std and Amplitude) 
            read(unit=ud_laser, FMT=*) Gauss_pulse_center, Gauss_pulse_width, Gauss_pulse_amplitude
            ! print *, Gauss_pulse_center, Gauss_pulse_width, Gauss_pulse_amplitude, laser_energy, laser_variation

          case DEFAULT
            print '(a)', 'RUMDEED: ERROR UNKNOWN LASER MODE'
            print *, LASER_MODE
            stop
        end select  


      case (GAUSS_OFF)
        ! Gauss Mode off
        print '(a)', 'RUMDEED: Using standard emission model'
        !EmitGauss = .FALSE.

        select case (LASER_MODE)
          case (LASER_FIXED)
            ptr_Get_Photo_Emission_Energy => Get_Fixed_Laser_Energy
            print '(a)', 'RUMDEED: Using fixed laser energy'    
            ! Read mean and std photon energy from the file
            read(unit=ud_laser, FMT=*) laser_energy, laser_variation

          case (LASER_POISSON)
            ptr_Get_Photo_Emission_Energy => Get_Laser_Energy
            print '(a)', 'RUMDEED: Using poisson distributed photon energy'
            ! Read mean and std photon energy from the file
            read(unit=ud_laser, FMT=*) laser_energy, laser_variation

          case DEFAULT
            print '(a)', 'RUMDEED: ERROR UNKNOWN LASER MODE'
            print *, LASER_MODE
            stop
        end select

      case DEFAULT
        print '(a)', 'RUMDEED: ERROR UNKNOWN LASER TYPE'
        print *, LASER_MODE
        stop
    end select

    select case (PHOTON_MODE)
      case (PHOTON_ZERO)
        print '(a)', 'RUMDEED: Using zero initial velocity'
      case (PHOTON_WFD)
        print '(a)', 'RUMDEED: Using work function velocity'
      case DEFAULT
        print '(a)', 'RUMDEED: ERROR UNKNOWN LASER MODE'
        print *, PHOTON_MODE
        stop
    end select
    close(unit=ud_laser)
  end subroutine Read_laser_parameters

  subroutine Init_Photo_Emission()
    integer :: i

    ! Allocate the number of emitters
    allocate(nrEmitted_emitters(1:nrEmit))

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_Planar
    
    ! Function for the electric field in the system
    ptr_field_E => field_E_planar

    ! The function that does the emission
    ptr_Do_Emission => Do_Photo_Emission

    ! The function to do image charge effects
    ptr_Image_Charge_effect => Force_Image_charges_v2

    call Read_work_function()

    call Read_laser_parameters()
    
    ! Initial velocity
    !call Init_Emission_Velocity(VELOCITY_ZERO)

    do i = 1, nrEmit
      if (emitters_type(i) == EMIT_RECTANGLE_SPOTS) then
        if (nrEmit > 1) then
          print *, 'ERROR: nrEmit must be 1 when using EMIT_RECTANGLE_SPOTS'
          stop
        end if

        print '(a)', 'RUMDEED: Doing spot emission'
        call Read_Spots_File()
      end if
    end do
  end subroutine Init_Photo_Emission

  subroutine Clean_Up_Photo_Emission()
    deallocate(nrEmitted_emitters)

    call Work_fun_cleanup()
  end subroutine Clean_Up_Photo_Emission

  ! Read in the file for the spot emission
  ! Reads x and y coords in nm and radius in nm
  subroutine Read_Spots_File()
    integer        :: i, IFAIL
    integer        :: ud_spot
    character(256) :: iomsg

    open(newunit=ud_spot, iostat=IFAIL, iomsg=iomsg, file='w_spots', &
       & status='OLD', form='FORMATTED', access='SEQUENTIAL', action='READ')
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file w_spots. ABORTING'
      print *, IFAIL
      print *, iomsg
      stop
    end if

    read(unit=ud_spot, FMT=*) num_spots
    print '(a)', 'RUMDEED: Reading spot file'
    
    !num_spots = 3

    allocate(spot_pos(1:2, 1:num_spots))
    allocate(spot_radius_sq(1:num_spots))

    do i = 1, num_spots
      read(unit=ud_spot, FMT=*) spot_pos(1:2, i), spot_radius_sq(i)
      spot_pos(1:2, i) = spot_pos(1:2, i) * length_scale
      spot_radius_sq(i) = (spot_radius_sq(i) * length_scale)**2
    end do

  end subroutine Read_Spots_File

  subroutine Do_Photo_Emission_Spot(step, emit)
    integer, intent(in)              :: step, emit
    integer                          :: nrElecEmit, nrTry
    double precision, dimension(1:3) :: par_pos, par_vel, field
    double precision                 :: p_eV

    ! Get Energy distribution of the photons
    p_eV = ptr_Get_Photo_Emission_Energy()
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
      par_pos(3) = 0.0d0 * length_scale ! Check in plane

      if (Emit_From_Spot(par_pos) .eqv. .false.) then
        cycle ! We do not emit from this position
      end if

      ! Check if photon energy exceeds work function
      if (w_theta_xy(par_pos, emit) <= p_eV) then

        field = Calc_Field_at(par_pos)
        if (field(3) < 0.0d0) then
          par_pos(3) = 1.0d0 * length_scale ! Above plane
          field = Calc_Field_at(par_pos)

          if (field(3) < 0.0d0) then
            par_pos(3) = 1.0d0 * length_scale ! Place above plane
            
            if (PHOTON_MODE == 1) then
              par_vel = 0.0d0
            else if (PHOTON_MODE == 2) then
              !par_vel = 0.0d0 ! Need to modify this to K = h (v - v_0) for excess energy
              par_vel(1:2) = 0.0d0 ! Set x and y components to zero
              par_vel(3) = sqrt((2.0d0 * ((p_eV - w_theta_xy(par_pos, emit))*q_0))/m_0) ! <-- Newtonian
            else
              print *, "WARNING: Unknown photon velocity mode!"
            end if

            call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1)

            nrElecEmit = nrElecEmit + 1
            nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
            nrTry = 0
          end if
        end if
      end if
    end do

    posInit   = posInit   + nrElecEmit
    nrEmitted = nrEmitted + nrElecEmit
  end subroutine Do_Photo_Emission_Spot

  ! Return true if we are allowed to emit from this position
  logical function Emit_From_Spot(par_pos)
    double precision, dimension(1:3), intent(in) :: par_pos
    integer                                      :: i
    double precision                             :: r_sq, k1, k2, x, y

    Emit_From_Spot = .false. ! Default to false

    x = par_pos(1)
    y = par_pos(2)

    ! Loop through all spots
    ! Equation for circle
    ! (x-k1)^2 + (y-k2)^2 = r^2
    do i = 1, num_spots
      k1 = spot_pos(1, i)
      k2 = spot_pos(2, i)

      r_sq = (x - k1)**2 + (y - k2)**2

      if (r_sq <= spot_radius_sq(i)) then
        Emit_From_Spot = .true. ! We emit from this position
        exit ! We can exit the loop. No need to check other spots
      end if
    end do
  end function Emit_From_Spot

  subroutine Do_Photo_Emission(step)
    integer, intent(in) :: step
    integer             :: i, IFAIL
    ! double precision    :: cur_time

    posInit            = 0
    nrEmitted_emitters = 0

    ! Check if we are doing a Gaussian distributed emission
    ! and set the max number of electrons allowed to be emitted if we are
    if (EmitGauss .eqv. .TRUE.) then
      maxElecEmit = Rand_Poisson( Gauss_Emission(step) )
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

        case (EMIT_RECTANGLE_SPOTS)
          !print *, 'Doing spots'
          !call Do_Photo_Emission_Rectangle_Spot(step, i)
          call Do_Photo_Emission_Spot(step, i)
          
        case default
          print *, 'RUMDEED: ERROR unknown emitter type!!'
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
    integer, intent(in)              :: step,       emit
    integer                          :: nrElecEmit, nrTry
    double precision, dimension(1:3) :: par_pos, par_vel, field
    !double precision                :: r_e, theta_e
    double precision                 :: r_e, r_e2, r2, p_eV, schk_wf, schk_red
    
    ! Get Energy distribution of the photons
    p_eV = ptr_Get_Photo_Emission_Energy()
    par_pos = 0.0d0
    nrTry = 0
    nrElecEmit = 0
    nrEmitted_emitters(emit) = 0

    do while (nrTry <= MAX_EMISSION_TRY)
      ! Check if we have reached the max number of electrons to be emitted,
      ! if we are using Gaussian limited emission.
      if ((nrElecEmit >= maxElecEmit) .and. (EmitGauss .eqv. .True.)) then
       exit
      end if

      ! Check if we have reached the maximum number of electrons.
      ! If this happends then the parameter MAX_PARTICLES (in mod_gobal)
      ! needs to be increased and the code recompiled.
      if (nrElec == MAX_PARTICLES-1) then
       print *, 'WARNING: Reached maximum number of electrons!!!'
       exit
      end if

      ! Find a random spot on the cathode.
      r_e = emitters_dim(1, emit) ! Radius of emitter
      r_e2 = r_e**2 ! Radius of emitter squared
      r2 = r_e2 + 1.0d0 ! Must be larger then r_e2 for the do while loop to run
      do while (r2 > r_e2)
        call random_number(par_pos(1:2)) ! Gives a random number [0,1]

        par_pos(1:2) = 2.0d0*(par_pos(1:2) - 0.5d0)*r_e ! Range is -r_e to +r_e
        r2 = par_pos(1)**2 + par_pos(2)**2 ! Radius squared of our random point
      end do
      ! Shifting position to comply with checkerboard workfunction location function
      par_pos(1:2) = emitters_pos(1:2, emit) + par_pos(1:2) + emitters_dim(1:2, emit)

      nrTry = nrTry + 1 ! Add one more try
      par_pos(3) = 0.0d0 * length_scale ! Check in plane
      
      ! Schottkey Effect, reduction of workfunction due to electric field
      ! Workfunction_reduction = sqrt ( (electron_charge^3 * electric_fielc)/(4 * pi * epsilon_0) )
      
      field = Calc_Field_at(par_pos)  
      if (field(3) < 0.0d0) then
        schk_red = sqrt( ( -field(3) *  (q_0**3)) / (4 * pi * epsilon_0)) / q_0
        schk_wf = w_theta_xy(par_pos, emit) - schk_red
        if (schk_wf <= p_eV) then
          par_pos(3) = 1.0d0 * length_scale ! Above plane
          field = Calc_Field_at(par_pos)

          if (field(3) < 0.0d0) then
            par_pos(3) = 1.0d0 * length_scale ! Place above plane
                        
            if (PHOTON_MODE == 1) then
              par_vel = 0.0d0
            else if (PHOTON_MODE == 2) then
              par_vel(1:2) = 0.0d0 ! Set x and y components to zero
              par_vel(3) = sqrt((2.0d0 * ((p_eV - schk_wf)*q_0))/m_0) ! <-- Newtonian
            else
              print *, "WARNING: Unknown photon velocity mode!"
            end if

            ! Escape velocity from image charge partner
            !par_vel(3) = q_0 / sqrt(8.0d0*pi*epsilon_0*epsilon_r*m_0*par_pos(3)) 
            
            ! Speed needed to reach over the gap spacing. Includes image charge partners behind the cathode and annode but not the electric field.
            !par_vel(3) = q_0/(4.0d0*sqrt(pi*epsilon_0*epsilon_r*m_0))*(d-2.0d0*par_pos(3))/sqrt(d*par_pos(3)*(d-par_pos(3)))
          
            call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1)

            nrElecEmit = nrElecEmit + 1
            nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
            nrTry = 0
          end if
        end if
      end if
    
    end do

    posInit   = posInit   + nrElecEmit
    nrEmitted = nrEmitted + nrElecEmit
  end subroutine Do_Photo_Emission_Circle

  subroutine Do_Photo_Emission_Rectangle(step, emit)
    integer, intent(in)              :: step,       emit
    integer                          :: nrElecEmit, nrTry
    double precision, dimension(1:3) :: par_pos,    par_vel, field
    double precision                 :: p_eV

    ! Get Energy distribution of the photons
    p_eV = ptr_Get_Photo_Emission_Energy()
    par_pos = 0.0d0
    nrTry = 0
    nrElecEmit = 0
    nrEmitted_emitters(emit) = 0

    do while (nrTry <= MAX_EMISSION_TRY)

      ! Check if we have reached the max number of electrons to be emitted
      if ((nrElecEmit >= maxElecEmit) .and. (maxElecEmit /= -1)) then
        exit
      end if


      if (nrElec == MAX_PARTICLES-1) then
        print *, 'WARNING: Reached maximum number of electrons!!!'
        exit
      end if

      call random_number(par_pos(1:2)) ! Gives a random number [0,1]

      par_pos(1:2) = emitters_pos(1:2, emit) + emitters_dim(1:2, emit)*par_pos(1:2)

      nrTry = nrTry + 1
      par_pos(3) = 0.0d0 * length_scale ! Check in plane
      
      ! Check if work function at position is lower then photon energy
      if (w_theta_xy(par_pos, emit) <= p_eV) then

        field = Calc_Field_at(par_pos)
        if (field(3) < 0.0d0) then
          par_pos(3) = 1.0d0 * length_scale ! Above plane
          field = Calc_Field_at(par_pos)

          if (field(3) < 0.0d0) then
            par_pos(3) = 1.0d0 * length_scale ! Place above plane
            
            if (PHOTON_MODE == 1) then
              par_vel = 0.0d0
            else if (PHOTON_MODE == 2) then
              par_vel(1:2) = 0.0d0 ! Set x and y components to zero
              par_vel(3) = sqrt((2.0d0 * ((p_eV - w_theta_xy(par_pos, emit))*q_0))/m_0) ! <-- Newtonian
            else
              print *, "WARNING: Unknown photon velocity mode!"
            end if
            
            call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1)

            nrElecEmit = nrElecEmit + 1
            nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
            nrTry = 0
          end if
        end if
      end if
    end do

    posInit = posInit + nrElecEmit
    nrEmitted = nrEmitted + nrElecEmit
  end subroutine Do_Photo_Emission_Rectangle

  ! Emission for a rectangular area with circular spots that are active.
  subroutine Do_Photo_Emission_Rectangle_Spot(step, emit)
    integer, intent(in)              :: step, emit
    integer                          :: nrElecEmit, nrTry
    double precision, dimension(1:3) :: par_pos, par_vel, field
    double precision                 :: p_eV
    
    ! Get Energy distribution of the photons
    p_eV = ptr_Get_Photo_Emission_Energy()
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

      call random_number(par_pos(1:2)) ! Gives a random number [0,1]

      par_pos(1:2) = emitters_pos(1:2, emit) + emitters_dim(1:2, emit)*par_pos(1:2)
      !par_pos(2) = emitters_pos(2, emit) + emitters_dim(2, emit)*par_pos(2)

      nrTry = nrTry + 1
      par_pos(3) = 0.0d0 * length_scale ! Check in plane
      
      if (w_theta_xy(par_pos, emit) <= p_eV) then

        field = Calc_Field_at(par_pos)
        if (field(3) < 0.0d0) then
          par_pos(3) = 1.0d0 * length_scale ! Above plane
          field = Calc_Field_at(par_pos)

          if (field(3) < 0.0d0) then
            par_pos(3) = 1.0d0 * length_scale ! Place above plane

            if (PHOTON_MODE == 1) then
              par_vel = 0.0d0
            else if (PHOTON_MODE == 2) then
              par_vel(1:2) = 0.0d0 ! Set x and y components to zero
              par_vel(3) = sqrt((2.0d0 * ((p_eV - w_theta_xy(par_pos, emit))*q_0))/m_0) ! <-- Newtonian
            else
              print *, "WARNING: Unknown photon velocity mode!"
            end if
            
            call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1)

            !print *, 'field = ', field
            !pause
            !call Add_Plane_Graph_emitt(par_pos, par_vel)

            nrElecEmit = nrElecEmit + 1
            nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
            nrTry = 0

          end if
        end if
      end if
    end do

    posInit   = posInit   + nrElecEmit
    nrEmitted = nrEmitted + nrElecEmit
  end subroutine Do_Photo_Emission_Rectangle_Spot

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

  ! ----------------------------------------------------------------------------    
  ! Gives a Gaussian emission curve
  ! where step is the current time step
  ! returns the number of electrons allowed to be emitted in that time step
  ! Needs rework from if to case for pulse series input - Hákon 15.12.21
  double precision function Gauss_Emission(step)
    integer, intent(in)          :: step ! Current time step
    integer                      :: IFAIL
    double precision             :: b, b1
    b = 1.0d0 / ( 2.0d0 * pi * Gauss_pulse_width**2 )
    b1 =  -1.0d0 * b * (step - Gauss_pulse_center)**2 
    ! For mutiple pulses - WIP
    !if (step <= 15000) then
    if (b1 < -500) then
      Gauss_Emission = 0.0d0
    else
      Gauss_Emission = Gauss_pulse_amplitude * exp(b1)
    end if
    !else
    !  Gauss_Emission = A * exp( -1.0d0 * b * (step - mu2)**2 )
    !end if
    write (ud_gauss, "(i6, tr2, i6)", iostat=IFAIL) step, Gauss_Emission
  end function Gauss_Emission

  ! ----------------------------------------------------------------------------
  ! Repurposed Maxwell-Boltzmann distribution for velocity generation
  ! https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution#Distribution_for_the_velocity_vector
  ! 
  ! Box-Muller method for normal distribution of photon energies
  ! Input is read from laser file with mean and standard deviation (std) in electronVolts (eV)
  function Get_Laser_Energy()
    double precision                 :: Get_Laser_Energy
    double precision, dimension(1:2) :: std
    double precision, dimension(1:2) :: mean
    double precision, dimension(1:3) :: Set_Laser_Energy

    mean = laser_energy
    std  = laser_variation ! Standard deviation of the Maxwell-Boltzmann distribution
    
    ! Get normal distributed numbers.
    ! The Box Muller method gives two numbers.
    ! We overwrite the second element in the array.
    Set_Laser_Energy(1:2) = box_muller(mean, std)
    Set_Laser_Energy(2:3) = box_muller(mean, std)

    Get_Laser_Energy = abs(Set_Laser_Energy(3)) ! Positive velocity in the z-direction
  end function Get_Laser_Energy

  function Get_Fixed_Laser_Energy()
    double precision                 :: Get_Fixed_Laser_Energy

    Get_Fixed_Laser_Energy = laser_energy
  end function Get_Fixed_Laser_Energy

  integer function Gauss_Emission_int(step)
    integer, intent(in)         :: step ! Current time step
    integer                     :: IFAIL
    double precision, parameter :: sigma = 1000.0d0 ! Width / standard deviation
    double precision, parameter :: mu = 3000.0d0 ! Center
    double precision, parameter :: A = 6.0d0 ! Height
    double precision, parameter :: b = 1.0d0/(2.0d0*pi*sigma**2)

    Gauss_Emission_int = IDNINT(  A * exp( -1.0d0*b*(step - mu)**2 )  )

    write (ud_gauss, "(i6, tr2, i6)", iostat=IFAIL) step, Gauss_Emission_int
  end function Gauss_Emission_int

end module mod_photo_emission
