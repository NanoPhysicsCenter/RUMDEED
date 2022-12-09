!
! Kristinn Torfason
! Reykjavík University Molecular Dynamics code for Electron Emission and Dynamics (RUMDEED)
!
program RUMDEED
! Vacuum
  use iso_fortran_env
#if defined(_OPENMP)
  use omp_lib
#endif
  use mod_global
  use mod_verlet
  use mod_photo_emission
  use mod_field_emission
  use mod_field_emission_v2
  use mod_emission_tip
  use mod_field_emission_2D
  use mod_field_thermo_emission
  !use mod_therminoic_emission
  use mod_pair
  use mod_unit_tests
  use mod_manual_emission
#if defined(__INTEL_COMPILER)
  use IFPORT ! Needed for getpid()
#endif
  implicit none
#if defined(__PGI)
  interface
    integer function getpid()
    end function getpid
  end interface
#endif


  integer                 :: i, nthreads
  integer, dimension(1:9) :: progress
  integer, dimension(8)   :: values ! Date and time

  call date_and_time(VALUES=values)

  print '(a)', 'Reykjavík University Molecular Dynamics code for Electron Emission and Dynamics'
  print '(a, i0.2, a, i0.2, a, i4, a, i0.2, a, i0.2, a, i0.2)', &
        'RUMDEED: Starting (d/m/y) ', values(3), '/', values(2), '/', values(1), &
        ' - (h:m:s) ', values(5), ':', values(6), ':', values(7)

  print '(4a)', 'This program was compiled by ', &
                compiler_version(), ' using the options ', &
                compiler_options()

#if defined(_GIT_VERSION_)
  print '(2a)', 'GIT version: ', _GIT_VERSION_
#endif

#if defined(_OPENMP)
  print '(a)', 'RUMDEED: Using OpenMP'
  nthreads = omp_get_max_threads()
  print '(tr1, a, i0)', 'Number of threads ', nthreads
#else
  print '(a)', 'RUMDEED: Single threaded'
  nthreads = 1
#endif


  print '(a)', 'RUMDEED: Reading input values'
  call Read_Input_Variables()

  print '(a)', 'RUMDEED: Initialzing'
  call Init()

  SELECT CASE (EMISSION_MODE)
  case(EMISSION_UNIT_TEST)
    print '(a)', 'RUMDEED: **** UNIT TESTING IS ACTIVE ****'
    call Init_Unit_Test()
  case(EMISSION_PHOTO)
    print '(a)', 'RUMDEED: Doing Photo emission'
    call Init_Photo_Emission()
  case(EMISSION_FIELD)
    print '(a)', 'RUMDEED: Doing Field emission'
    call Init_Field_Emission()
  case(EMISSION_TIP)
    print '(a)', 'RUMDEED: Doing Field emission from a tip'
    call Init_Emission_Tip()
  case(EMISSION_FIELD_2D_2DEG_C, EMISSION_FIELD_2D_2DEG_NC, EMISSION_FIELD_2D_DIRAC_C, EMISSION_FIELD_2D_DIRAC_NC)
    print '(a)', 'RUMDEED: Doing Field emission from 2D material'
    call Init_Field_Emission_2D()
  case(EMISSION_THERMIONIC)
    print '(a)', 'RUMDEED: Doing Thermionic emission'
    !call Init_Thermionic_Emission()
  case(EMISSION_FIELD_THERMO)
    print '(a)', 'RUMDEED: Doing General Field+Thermionic emission'
    call Init_Field_Thermo_Emission()
  case(EMISSION_FIELD_V2)
    print '(a)', 'RUMDEED: Doing Field emission V2'
    call Init_Field_Emission_v2()
  case(EMISSION_MANUAL)
    print '(a)', 'RUMDEED: Doing manual emission'
    call Init_Manual()
  case DEFAULT
    print '(a)', 'RUMDEED: ERROR UNKNOWN EMISSION MODEL'
    print *, EMISSION_MODE
    stop
  END SELECT

  print '(tr1, a, i0)', 'Number of emitters ', nrEmit

  print '(a)', 'RUMDEED: Writing out variables'
  call Write_Initial_Variables()

  !print '(a)', 'RUMDEED: Setting up dipoles'
  !call Setup_Dipoles(5, 5, 5)
  !call Init_Dipoles(0, 0, 0)
  !call Write_Dipole_data()

  print '(a)', 'RUMDEED: Starting main loop'
  print '(tr1, a, i0, a, ES12.4, a)', 'Doing ', steps, ' time steps of size ', time_step/1.0E-12, ' ps'

  print *, ''

  cur_time = 0
  call Set_Voltage(0) ! Set voltage for time step 0

  do i = 1, steps

    ! Do Emission
    !print *, 'Emission'
    call ptr_Do_Emission(i)

    ! Update the position of all particles
    !print *, 'Update position'
    call Update_Position(i)
    call Write_Position(i)
    !call Write_Position_XYZ_Step(i)

    ! Remove particles from the system
    !print *, 'Remove particles'
    call Remove_Particles(i)

    ! Do Collisions
    call Do_Collisions(i)

    ! Flush data
    call Flush_Data()

    !do j = 1, 9
    !  if (i == progress(j)) then
    !    call PrintProgress(j)
    !    call Flush_Data()
    !    exit
    !  end if
    !end do

    ! Check the progress, if we are at 10%, 20%, 30%, ...
    ! then print the progress and flush data
    if (any(i .eq. progress) .eqv. .true.) then
      call PrintProgress(nint(i*10.0d0/steps))
      !call Flush_Data()
    end if

    if (cought_stop_signal .eqv. .true.) then
      print '(a)', 'RUMDEED: Got signal SIGINT stopping main loop.'
      exit ! We cought the signal to stop. Exit the main loop.
    end if

    ! if (i == progress(1)) then
    !   call PrintProgress(1)
    ! else if (i == progress(2)) then
    !   call PrintProgress(2)
    ! else if (i == progress(3)) then
    !   call PrintProgress(3)
    ! else if (i == progress(4)) then
    !   call PrintProgress(4)
    ! else if (i == progress(5)) then
    !   call PrintProgress(5)
    ! else if (i == progress(6)) then
    !   call PrintProgress(6)
    ! else if (i == progress(7)) then
    !   call PrintProgress(7)
    ! else if (i == progress(8)) then
    !   call PrintProgress(8)
    ! else if (i == progress(9)) then
    !   call PrintProgress(9)
    ! end if

  end do

  call PrintProgress(10)
  print '(a)', 'RUMDEED: Main loop finished'


  print '(a)', 'RUMDEED: Writing data'
  call Write_Life_Time()

  print '(a)', 'RUMDEED: Emission clean up'
  SELECT CASE (EMISSION_MODE)
  case(EMISSION_UNIT_TEST)
    call Clean_Up_Unit_Test()
  case(EMISSION_PHOTO)
    call Clean_Up_Photo_Emission()
  case(EMISSION_FIELD)
    call Clean_Up_Field_Emission()
  case(EMISSION_TIP)
    call Clean_Up_Emission_Tip()
  case(EMISSION_FIELD_2D_2DEG_C, EMISSION_FIELD_2D_2DEG_NC, EMISSION_FIELD_2D_DIRAC_C, EMISSION_FIELD_2D_DIRAC_NC)
    call Clean_Up_Field_Emission_2D()
  !case(EMISSION_THERMIONIC)
    !call Clean_Up_Thermionic_Emission()
  case(EMISSION_FIELD_THERMO)
    call Clean_Up_Field_Thermo_Emission()
  case(EMISSION_FIELD_V2)
    call Clean_Up_Field_Emission_v2()
  case(EMISSION_MANUAL)
    call Clean_Up_Manual()
  END SELECT

  print '(a)', 'RUMDEED: Final clean up'
  call Clean_up()
  call date_and_time(VALUES=values)
  print '(a, i0.2, a, i0.2, a, i4, a, i0.2, a, i0.2, a, i0.2)', &
        'RUMDEED: Ended (d/m/y) ', values(3), '/', values(2), '/', values(1), &
        ' - (h:m:s) ', values(5), ':', values(6), ':', values(7)
  print '(a)', 'RUMDEED: Program finished'
contains

  subroutine PrintProgress(i)
    integer, intent(in) :: i

    print '(a, i3, a)', 'RUMDEED: ', (i*10), '% done'
    print '(tr1, i0, a)', nrPart, ' particles in the system'
    print '(tr2, i0, a, i0, a)', nrElec, ' electrons and ', nrHole, ' holes'
    print *, ''
  end subroutine PrintProgress

  subroutine Signal_Handler()
    cought_stop_signal = .true.
  end subroutine Signal_Handler

  ! ----------------------------------------------------------------------------
  ! A subroutine to read the input variable from the file 'input'
  subroutine Read_Input_Variables()
    integer :: ud_input, IFAIL


    !Open the file 'input' for reading and check for errors
    open(newunit=ud_input, iostat=IFAIL, file='input', status='OLD', action='read')
    if (IFAIL /= 0) then
      print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file input'
      stop
    end if

    allocate(emitters_pos(1:3, 1:MAX_EMITTERS))
    allocate(emitters_dim(1:3, 1:MAX_EMITTERS))
    allocate(emitters_type(1:MAX_EMITTERS))
    allocate(emitters_delay(1:MAX_EMITTERS))
    allocate(nrElec_remove_top_emit(1:MAX_EMITTERS))
    allocate(ramo_current_emit(1:MAX_SECTIONS, 1:MAX_EMITTERS))

    emitters_dim = 0.0d0
    emitters_pos = 0.0d0
    emitters_type = 0
    emitters_delay = 0
    nrElec_remove_top_emit(1:MAX_EMITTERS) = 0

    ! Read the input file
    ! We read in
    ! V_s: Voltage over the system
    ! box_dim: Size of the system
    ! time_step: Size of the time step, given in ps
    ! steps: Number of time_steps to do
    ! nrEmit: Number of emitters
    ! emitters_pos: Positions of emitters (x and y, z is ignored)
    ! emitters_dim: Dimensions of emitters (x and y, z is ignored)
    ! emitters_type: Type of emitter (Circle or square)
    ! emitters_delay: When the emitter should start emitting ()
    ! EMISSION_MODE: What emission to do (Photo/CL, Field emission, e.t.c)
    read(ud_input, NML=input)

    ! Close the 'input' file
    close(unit=ud_input, iostat=IFAIL, status='keep')

    ! box_dim: Dimensions of the system
    ! d: Gap spacingelse
    box_dim = box_dim * length_scale
    d = box_dim(3)

    ! V_s: Voltage over the gap, d
    V_d = V_s
    ! Electric field in the system for a planar case
    E_z = -1.0d0*V_d/d

    ! Electric field for unit voltage (See ramo current)
    E_zunit = -1.0d0/d

    ! Emitters position and dimensions are given in length_scale (nm)
    emitters_dim = emitters_dim * length_scale
    emitters_pos = emitters_pos * length_scale
  
    !dens_x_d = box_dim(1) / (N_x_densmap-1)
    !dens_y_d = box_dim(2) / (N_y_densmap-1)

    !time_step: The size of the time step
    time_step = time_step * time_scale
    time_step2 = time_step**2 ! Time_step squared

    ! Temperature, pressure and density
    P_abs = P_abs * P_ntp 
    n_d = P_abs/(k_b*T_temp)

    ! steps: Number of time steps that the program will simulate
    if (steps <= 0) then
      print '(a)', 'ERROR: steps <= 0'
      stop
    end if

    call Read_Cross_Section_Data()

    open(newunit=ud_debug, iostat=IFAIL, file='debug.dt', status='replace', action='write')
    if (IFAIL /= 0) then
      print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file debug.dt'
      stop
    end if

    !write(ud_debug, NML=input_test)

    !close(unit=ud_debug, iostat=IFAIL, status='keep')

  end subroutine Read_Input_Variables

  ! ----------------------------------------------------------------------------
  ! Initialize the program
  ! allocate and initilize variables, open data files for writing, etc.
  subroutine Init()
    integer :: IFAIL
    character(len=128)  :: filename
    !integer, dimension(:), allocatable :: my_seed

    ! Allocate arrays
    allocate(particles_cur_pos(1:3, 1:MAX_PARTICLES))
    allocate(particles_prev_pos(1:3, 1:MAX_PARTICLES))
    allocate(particles_last_col_pos(1:3, 1:MAX_PARTICLES))
    allocate(particles_cur_vel(1:3, 1:MAX_PARTICLES))
    !allocate(particles_prev_vel(1:3, 1:MAX_PARTICLES))
    allocate(particles_cur_accel(1:3, 1:MAX_PARTICLES))
    allocate(particles_prev_accel(1:3, 1:MAX_PARTICLES))
    allocate(particles_prev2_accel(1:3, 1:MAX_PARTICLES))
    allocate(particles_charge(1:MAX_PARTICLES))
    allocate(particles_species(1:MAX_PARTICLES))
    allocate(particles_mass(1:MAX_PARTICLES))
    allocate(particles_step(1:MAX_PARTICLES))
    allocate(particles_mask(1:MAX_PARTICLES))
    allocate(particles_emitter(1:MAX_PARTICLES))
    allocate(particles_section(1:MAX_PARTICLES))
    allocate(particles_life(1:MAX_PARTICLES))
    allocate(particles_id(1:MAX_PARTICLES))

    allocate(life_time(1:MAX_LIFE_TIME, 1:2))
    allocate(ramo_current(1:nrSpecies))

    allocate(density_map_elec(1:N_x_densmap, 1:N_y_densmap))
    allocate(density_map_hole(1:N_x_densmap, 1:N_y_densmap))

    particles_cur_pos     = 0.0d0
    particles_prev_pos    = 0.0d0
    particles_cur_vel     = 0.0d0
    particles_cur_accel   = 0.0d0
    particles_prev_accel  = 0.0d0
    particles_prev2_accel = 0.0d0
    particles_charge      = 0.0d0
    particles_species     = 0
    particles_mass        = 0.0d0
    particles_step        = 0
    particles_mask        = .true.
    particles_id          = 0

    ramo_current = 0.0d0
    ramo_current_emit(1:MAX_SECTIONS, 1:MAX_EMITTERS) = 0.0d0
    life_time = 0

    density_map_elec = 0
    density_map_hole = 0

    !V_cur = 0.0d0
    !V_prev = 0.0d0

    ! Start with an empty system
    nrPart      = 0
    nrElec      = 0
    nrHole      = 0
    nrElecHole  = 0

    !startElecHoles = 1
    !endElecHoles   = 0


    nrPart_remove = 0
    nrElec_remove = 0
    nrHole_remove = 0

    nrPart_remove_top = 0
    nrPart_remove_bot = 0
    nrElec_remove_top = 0
    nrElec_remove_bot = 0
    nrHole_remove_top = 0
    nrHole_remove_bot = 0

    ! ID starts with 0
    nrID = 0

    !progress(10) = 1.0*steps
    progress(9) = nint(0.9d0*steps)
    progress(8) = nint(0.8d0*steps)
    progress(7) = nint(0.7d0*steps)
    progress(6) = nint(0.6d0*steps)
    progress(5) = nint(0.5d0*steps)
    progress(4) = nint(0.4d0*steps)
    progress(3) = nint(0.3d0*steps)
    progress(2) = nint(0.2d0*steps)
    progress(1) = nint(0.1d0*steps)

    !call RANDOM_SEED(size = n)
    !allocate(my_seed(n))
    !my_seed = SEED
    !call RANDOM_SEED(PUT = my_seed)
    !deallocate(my_seed)

    call init_random_seed()

    ! Register a subroutine to catch the SIGINT signal
#if defined(__GNUC__)
    IFAIL = SIGNAL(SIGINT, Signal_Handler)
#endif

    ! Create folder for output files
#if defined(__PGI)
    call system('mkdir -p out/')
#else
    call execute_command_line ('mkdir -p out/')
#endif

    ! Set the number of cores Cuba should use to 0, i.e. no parallelization in cuba.
    call cubacores(0, 1000)


    ! Open data file for writing
    open(newunit=ud_pos, iostat=IFAIL, file='out/position.bin', status='REPLACE', action='write', access='STREAM')
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file position.bin. ABORTING'
      stop
    end if

    ! open(newunit=ud_vel, iostat=IFAIL, file='velocity.dt', status='REPLACE', action='write')
    ! if (IFAIL /= 0) then
    !   print *, 'RUMDEED: Failed to open file velocity.dt. ABORTING'
    !   stop
    ! end if

    open(newunit=ud_emit, iostat=IFAIL, file='out/emitted.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file emitted.dt. ABORTING'
      stop
    end if

    open(newunit=ud_absorb, iostat=IFAIL, file='out/absorbed.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file absorbed.dt. ABORTING'
      stop
    end if

    open(newunit=ud_absorb_top, iostat=IFAIL, file='out/absorbed_top.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file absorbed_top.dt. ABORTING'
      stop
    end if

    open(newunit=ud_absorb_bot, iostat=IFAIL, file='out/absorbed_bot.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file absorbed_bot.dt. ABORTING'
      stop
    end if

    open(newunit=ud_ramo, iostat=IFAIL, file='out/ramo_current.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file ramo_current.dt. ABORTING'
      stop
    end if

    open(newunit=ud_ramo_sec, iostat=IFAIL, file='out/ramo_current.bin', status='REPLACE', action='WRITE', &
       & access='STREAM', asynchronous='YES')
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file ramo_current.bin. ABORTING'
      stop
    end if

    open(newunit=ud_volt, iostat=IFAIL, file='out/volt.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file volt.dt. ABORTING'
      stop
    end if

    open(newunit=ud_field, iostat=IFAIL, file='out/field.dt', status='replace', action='write')
    if (IFAIL /= 0) then
      print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file field.dt'
      stop
    end if

    open(newunit=ud_coll, iostat=IFAIL, file='out/collisions.dt', status='replace', action='write')
    if (IFAIL /= 0) then
      print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file collisions.dt'
      stop
    end if

    open(newunit=ud_integrand, iostat=IFAIL, file='out/integration.dt', status='replace', action='write')
    if (IFAIL /= 0) then
      print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file integration.dt'
      stop
    end if

    open(newunit=ud_gauss, iostat=IFAIL, file='out/gauss.dt', status='replace', action='write')
    if (IFAIL /= 0) then
      print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file gauss.dt'
      stop
    end if

    open(newunit=ud_mh, iostat=IFAIL, file='out/mh.bin', &
    status='REPLACE', action='WRITE', access='STREAM')
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file mh.bin. ABORTING'
      stop
    end if
    
    !------------------------------------------------------------------------------------
    ! Emission density
    open(newunit=ud_density_emit, iostat=IFAIL, file='out/density_emit.bin', &
    status='REPLACE', action='WRITE', access='STREAM')
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file density_emit.bin. ABORTING'
      stop
    end if

    open(newunit=ud_density_ion, iostat=IFAIL, file='out/density_ion.bin', &
    status='REPLACE', action='WRITE', access='STREAM')
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file density_ion.bin. ABORTING'
      stop
    end if

    !-------------------------------------------------------------------------------------
    ! Absorbsion density
    open(newunit=ud_density_absorb_top, iostat=IFAIL, file='out/density_absorb_top.bin', &
         status='REPLACE', action='WRITE', access='STREAM')
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file density_absorb_top.bin. ABORTING'
      stop
    end if

    open(newunit=ud_density_absorb_bot, iostat=IFAIL, file='out/density_absorb_bot.bin', &
         status='REPLACE', action='WRITE', access='STREAM')
    if (IFAIL /= 0) then
      print *, 'RUMDEED: Failed to open file density_absorb_bin.dt. ABORTING'
      stop
    end if

    !-------------------------------------------------------------------------------------
    ! Files for planes
    do i = 1, planes_N
      if (planes_z(i) > 0.0d0) then
        write(filename, '(a11, i0, a3)') 'out/planes-', i, '.dt'
        open(newunit=planes_ud(i), iostat=IFAIL, file=filename, status='REPLACE', action='WRITE', access='STREAM')
        if (IFAIL /= 0) then
          print *, 'RUMDEED: Failed to open file for planes. ABORTING'
          stop
        end if
      else
        planes_ud(i) = -1
      end if
    end do

  end subroutine Init

  !---------------------------------------------------------------------------------------
  ! Set a random seed to use
  ! Taken from
  ! https://gcc.gnu.org/onlinedocs/gcc-4.7.4/gfortran/RANDOM_005fSEED.html
  subroutine init_random_seed()
    implicit none
    !integer, allocatable :: my_seed(:)
    integer :: i, n, un, istat, dt(8), pid, t(2), s
    integer(8) :: count, tms
  
    call random_seed(size = n)
    allocate(my_seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) my_seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(count)
       if (count /= 0) then
          t = transfer(count, t)
       else
          call date_and_time(values=dt)
          tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24 * 60 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
          t = transfer(tms, t)
       end if
       s = ieor(t(1), t(2))
       pid = getpid() + 1099313 ! Add a prime
       s = ieor(s, pid)
       if (n >= 3) then
          my_seed(1) = t(1) + 36343
          my_seed(2) = t(2) + 72679
          my_seed(3) = pid
          if (n > 3) then
             my_seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
          end if
       else
          my_seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
       end if
    end if
    call random_seed(put=my_seed)
  end subroutine init_random_seed

  ! ----------------------------------------------------------------------------
  ! Flush data written to files such that it can be read
  ! subroutine Flush_Data()

  !   flush(ud_emit)
  !   flush(ud_absorb)
  !   flush(ud_volt)
  !   flush(ud_field)
  !   flush(ud_debug)

  !   flush(ud_density_emit)
  !   flush(ud_density_ion)

  !   flush(ud_density_absorb_top)
  !   flush(ud_density_absorb_bot)

  ! end subroutine Flush_Data


  ! ----------------------------------------------------------------------------
  ! Write out parameters, constants and initial variables used when the program starts
  subroutine Write_Initial_Variables()
    implicit none
    integer                     :: ud_init, IFAIL
    character(len=*), parameter :: fmt_int = "(a, tr8, i8, tr6, a)"
    character(len=*), parameter :: fmt_rel = "(a, tr2, ES16.6, tr2, a)"
    !character(len=*), parameter :: fmt_cmp = "(a, tr2, '(', G16.6, ',', G16.6, ')', tr2, a)"

    ! Open file 'init.dt' for writing
    open(newunit=ud_init, iostat=IFAIL, file='out/init.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Failed to open file init.dt. ABORTING'
      stop
    end if


    write(ud_init, *) 'Vacuum Electronics Molecular Dynamics'
    write(ud_init, *) 'Values used for constans and variables'
    write(ud_init, *)

    write(ud_init, *) '---------------------------------------------------------'
    write(ud_init, '(a)') material
    write(ud_init, fmt_rel) 'epsilon_r  = ', epsilon_r, 'Relative permittivity'
    write(ud_init, fmt_rel) 'm_eeff     = ', m_eeff,    'Effective mass for electrons'
    write(ud_init, fmt_rel) 'm_heff     = ', m_heff,    'Effective mass for holes'
    write(ud_init, *) '---------------------------------------------------------'
    write(ud_init, *)
    write(ud_init, *)

    write(ud_init, *) '---------------------------------------------------------'
    write(ud_init, fmt_rel) 'length_scale    = ', length_scale, 'The length scale used for output (m)'
    write(ud_init, fmt_rel) 'time_scale      = ', time_scale,   'The time scale used for output (s)'
    write(ud_init, fmt_rel) 'vel_scale       = ', vel_scale,    'The velocity scale used for output (m/s)'
    write(ud_init, fmt_rel) 'cur_scale       = ', cur_scale,    'The current scale used for output (A)'
    write(ud_init, fmt_rel) 'time_step       = ', time_step,    'Size of the time_step used in the Verlet integration (s)'
    write(ud_init, fmt_int) 'steps           = ', steps,        'Number of time steps in the Verlet integration'
    write(ud_init, *) '---------------------------------------------------------'
    write(ud_init, *)
    write(ud_init, *)

    !write(ud_init, fmt_rel) 'pre_fac_c     = ', pre_fac_c,     'Coulomb prefactor e_q^2/(4 \pi \epsilon_0 m_e)'
    !write(ud_init, fmt_rel) 'pre_fac_E     = ', pre_fac_E,     'Electric field prefactor q_e/m_e'
    !write(ud_init, *)

    write(ud_init, *) '---------------------------------------------------------'
    write(ud_init, fmt_int) 'MAX_PARTICLES       = ', MAX_PARTICLES, 'Maximum number of electrons in the system'
    !write(ud_init, fmt_int) 'SEED                = ', SEED,          'Seed value used in the random number generator'
    !write(ud_init, fmt_int) 'MAX_EMISSION_TRY    = ', MAX_EMISSION_TRY,    'Maximum number of failed emission attempts'
    !write(ud_init, fmt_int) 'MAX_TIME_STEP_WRITE = ', MAX_TIME_STEP_WRITE, 'Maximum number of times steps to output data'
    write(ud_init, *) '---------------------------------------------------------'
    write(ud_init, *)
    write(ud_init, *)

    write(ud_init, *) '---------------------------------------------------------'
    write(ud_init, fmt_rel) 'V_s           = ', V_s,        'The Voltage from the source'
    write(ud_init, fmt_rel) 'box_dim_x     = ', box_dim(1), 'System x-length'
    write(ud_init, fmt_rel) 'box_dim_y     = ', box_dim(2), 'System y-length'
    write(ud_init, fmt_rel) 'box_dim_z     = ', box_dim(3), 'System z-length'
    write(ud_init, fmt_rel) 'd             = ', d,          'Gap spacing in the system'
    write(ud_init, fmt_rel) 'E_z           = ', E_z,        'Y-component of the electric field in the system'
    write(ud_init, *) '---------------------------------------------------------'
    !write(ud_init, *)
    !write(ud_init, *)
    write(ud_init, *) '---------------------------------------------------------'
    write(ud_init, fmt_int) 'NrEmit        = ', NrEmit,        'Number of emitters'
    write(ud_init, fmt_int) 'EMISSION_MODE = ', EMISSION_MODE, 'The emission mechanism'
    write(ud_init, *) '---------------------------------------------------------'
    write(ud_init, *) 'GIT VERSION'
#if defined(_GIT_VERSION_)
    write(ud_init, *) _GIT_VERSION_
#endif

    ! Close file 'init.dt'
    close(unit=ud_init, iostat=IFAIL, status='keep')

    open(newunit=ud_init, iostat=IFAIL, file='out/init.bin', status='REPLACE', action='write', access='STREAM')
    if (IFAIL /= 0) then
      print *, 'Failed to open file init.dt. ABORTING'
      stop
    end if

    write(unit=ud_init) epsilon_r, m_eeff, m_heff, length_scale, time_scale, vel_scale, cur_scale, &
                        MAX_PARTICLES, MAX_EMITTERS, MAX_SECTIONS, MAX_LIFE_TIME

    close(unit=ud_init, iostat=IFAIL, status='keep')
  end subroutine Write_Initial_Variables


  ! ----------------------------------------------------------------------------
  ! Clean up after the program
  ! Deallocate variables, close files, etc.
  subroutine Clean_up()
    integer :: i

    ! Close file descriptors
    close(unit=ud_pos, status='keep')
    !close(unit=ud_vel, iostat=IFAIL, status='keep')
    close(unit=ud_emit,  status='keep')
    close(unit=ud_absorb,  status='keep')
    close(unit=ud_absorb_top, status='keep')
    close(unit=ud_absorb_bot, status='keep')
    close(unit=ud_ramo, status='keep')
    close(unit=ud_ramo_sec, status='keep')
    close(unit=ud_volt, status='keep')
    close(unit=ud_field,status='keep')
    close(unit=ud_coll,status='keep')
    close(unit=ud_integrand, status='keep')
    close(unit=ud_gauss, status='keep')
    close(unit=ud_debug, status='keep')

    close(unit=ud_density_emit, status='keep')
    close(unit=ud_density_ion, status='keep')
    close(unit=ud_density_absorb_top, status='keep')
    close(unit=ud_density_absorb_bot, status='keep')

    do i = 1, planes_N
      close(unit=planes_ud(i), status='keep')
    end do

    ! Deallocate arrays
    deallocate(particles_cur_pos)
    deallocate(particles_prev_pos)
    deallocate(particles_last_col_pos)
    deallocate(particles_cur_vel)
    !deallocate(particles_prev_vel)
    deallocate(particles_cur_accel)
    deallocate(particles_prev_accel)
    deallocate(particles_prev2_accel)
    deallocate(particles_charge)
    deallocate(particles_species)
    deallocate(particles_mass)
    deallocate(particles_mask)
    deallocate(particles_emitter)
    deallocate(particles_section)
    deallocate(particles_life)
    deallocate(particles_id)
    deallocate(particles_step)

    deallocate(emitters_pos)
    deallocate(emitters_dim)
    deallocate(emitters_type)
    deallocate(emitters_delay)
    deallocate(nrElec_remove_top_emit)
    deallocate(ramo_current_emit)

    deallocate(ramo_current)
    deallocate(life_time)

    deallocate(density_map_elec)
    deallocate(density_map_hole)

    deallocate(my_seed)

  end subroutine Clean_up
end program RUMDEED