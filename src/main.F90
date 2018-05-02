program VacuumMD
#if defined(_OPENMP)
  use omp_lib
#endif
  use mod_global
  use mod_verlet
  use mod_photo_emission
  use mod_field_emission
  !use mod_therminoic_emission
  use mod_pair
  use mod_unit_tests
  implicit none

  integer :: i, nthreads
  integer, dimension(1:9) :: progress

#if defined(_OPENMP)
  print '(a)', 'Vacuum: Using OpenMP'
  nthreads = omp_get_max_threads()
  print '(tr1, a, i0)', 'Number of threads ', nthreads
#else
  print '(a)', 'Vacuum: Single threaded'
  nthreads = 1
#endif

#if defined(_UNIT_TEST_)
  print '(a)', 'Vacuum: **** UNIT TESTING IS ACTIVE ****'
#endif

  print '(a)', 'Vacuum: Reading input values'
  call Read_Input_Variables()

  print '(a)', 'Vacuum: Initialzing'
  call Init()

  SELECT CASE (EMISSION_MODE)
  case(EMISSION_PHOTO)
    print '(a)', 'Vacuum: Doing Photo emission'
    call Init_Photo_Emission()
  case(EMISSION_FIELD)
    print '(a)', 'Vacuum: Doing Field emission'
    call Init_Field_Emission()
  case DEFAULT
    print '(a)', 'Vaccum: ERROR UNKNOWN EMISSION MODEL'
    stop
  END SELECT

  print '(a)', 'Vacuum: Writing out variables'
  call Write_Initial_Variables()

  !print '(a)', 'Vacuum: Setting up dipoles'
  !call Setup_Dipoles(5, 5, 5)
  !call Init_Dipoles(0, 0, 0)
  !call Write_Dipole_data()

#if defined(_UNIT_TEST_)
  print '(a)', 'Vacuum: Starting Unit Tests'
  print *, ''
#else
  print '(a)', 'Vacuum: Starting main loop'
  print '(tr1, a, i0, a, ES12.4, a)', 'Doing ', steps, ' time steps of size ', time_step, ' seconds'
#endif

#if defined(_UNIT_TEST_)
  call Run_Unit_Tests()
#else

  print *, ''

  cur_time = 0
  call Set_Voltage(0) ! Set voltage for time step 0

  do i = 1, steps

    ! Do Emission
    call ptr_Do_Emission(i)

    ! Update the position of all particles
    !print *, 'Update position'
    call Update_Position(i)

    ! Remove particles from the system
    call Remove_Particles(i)

    ! Flush data
    call Flush_Data()

    if (i == progress(1)) then
      call PrintProgress(1)
    else if (i == progress(2)) then
      call PrintProgress(2)
    else if (i == progress(3)) then
      call PrintProgress(3)
    else if (i == progress(4)) then
      call PrintProgress(4)
    else if (i == progress(5)) then
      call PrintProgress(5)
    else if (i == progress(6)) then
      call PrintProgress(6)
    else if (i == progress(7)) then
      call PrintProgress(7)
    else if (i == progress(8)) then
      call PrintProgress(8)
    else if (i == progress(9)) then
      call PrintProgress(9)
    end if

  end do

  call PrintProgress(10)
  print '(a)', 'Vacuum: Main loop finished'

! End of else for unit test
#endif


#if defined(_UNIT_TEST_)
  print '(a)', 'Vacuum: Unit tests finished'
#else
  print '(a)', 'Vacuum: Writing data'
  call Write_Life_Time()
#endif

  print '(a)', 'Vacuum: Program finished'
  call Clean_Up_Photo_Emission()
  call Clean_up()
contains

  subroutine PrintProgress(i)
    integer, intent(in) :: i

    print '(a, i3, a)', 'Vacuum: ', (i*10), '% done'
    print '(tr1, i0, a)', nrPart, ' particles in the system'
    print '(tr2, i0, a, i0, a)', nrElec, ' electrons and ', nrHole, ' holes'
    print *, ''
  end subroutine PrintProgress

  ! ----------------------------------------------------------------------------
  ! A subroutine to read the input variable from the file 'input'
  subroutine Read_Input_Variables()
    integer :: ud_input, IFAIL


    !Open the file 'input' for reading and check for errors
    open(newunit=ud_input, iostat=IFAIL, file='input', status='OLD', action='read')
    if (IFAIL /= 0) then
      print '(a)', 'Vacuum: ERROR UNABLE TO OPEN file input'
      stop
    end if

    allocate(emitters_pos(1:3, 1:MAX_EMITTERS))
    allocate(emitters_dim(1:3, 1:MAX_EMITTERS))
    allocate(emitters_type(1:MAX_EMITTERS))
    allocate(emitters_delay(1:MAX_EMITTERS))

    emitters_dim = 0.0d0
    emitters_pos = 0.0d0
    emitters_type = 0
    emitters_delay = 0

    ! Read the input file
    read(ud_input, NML=input)

    ! Close the 'input' file
    close(unit=ud_input, iostat=IFAIL, status='keep')

    ! V: Voltage over the gap
    ! box_dim: Dimensions of the system
    ! d: Gap spacing
    box_dim = box_dim * length_scale
    d = box_dim(3)

    V_d = V_s
    E_z = -1.0d0*V_d/d
    !print *, E_z

    E_zunit = -1.0d0/d

    emitters_dim = emitters_dim * length_scale
    emitters_pos = emitters_pos * length_scale

    ! Set the dimensions of the box used
    !box_dim(1) = d
    !box_dim(2) = d
    !box_dim(3) = d

    dens_x_d = box_dim(1) / (N_x_densmap-1)
    dens_y_d = box_dim(2) / (N_y_densmap-1)

    ! Set the current scale to mA / micrometer^2
    !cur_scale = 1.0d-3 * (box_dim(1)*1.0d6) * (box_dim(2)*1.0d6) ! mA / micrometer^2
    cur_scale = 1.0d-6 ! microAmper

    !time_step: The size of the time step
    time_step = time_step * time_scale
    time_step2 = time_step**2 ! Time_step squared

    ! steps: Number of time steps that the program will simulate
    if (steps <= 0) then
      print '(a)', 'ERROR: steps <= 0'
      stop
    end if

    open(newunit=ud_debug, iostat=IFAIL, file='debug.dt', status='replace', action='write')
    if (IFAIL /= 0) then
      print '(a)', 'Vacuum: ERROR UNABLE TO OPEN file debug.dt'
      stop
    end if

    !write(ud_debug, NML=input_test)

    !close(unit=ud_debug, iostat=IFAIL, status='keep')

  end subroutine Read_Input_Variables

  ! ----------------------------------------------------------------------------
  ! Initialize the program
  ! allocate and initilize variables, open data files for writing, etc.
  subroutine Init()
    integer :: IFAIL, n
    integer, dimension(:), allocatable :: my_seed

    ! Allocate arrays
    allocate(particles_cur_pos(1:3, 1:MAX_PARTICLES))
    allocate(particles_prev_pos(1:3, 1:MAX_PARTICLES))
    allocate(particles_cur_vel(1:3, 1:MAX_PARTICLES))
    !allocate(particles_prev_vel(1:3, 1:MAX_PARTICLES))
    allocate(particles_cur_accel(1:3, 1:MAX_PARTICLES))
    allocate(particles_prev_accel(1:3, 1:MAX_PARTICLES))
    allocate(particles_charge(1:MAX_PARTICLES))
    allocate(particles_species(1:MAX_PARTICLES))
    allocate(particles_mass(1:MAX_PARTICLES))
    allocate(particles_step(1:MAX_PARTICLES))
    allocate(particles_mask(1:MAX_PARTICLES))

    allocate(life_time(1:MAX_LIFE_TIME, 1:2))
    allocate(ramo_current(1:nrSpecies))

    allocate(density_map_elec(1:N_x_densmap, 1:N_y_densmap))
    allocate(density_map_hole(1:N_x_densmap, 1:N_y_densmap))

    particles_cur_pos = 0.0d0
    particles_prev_pos = 0.0d0
    particles_cur_vel = 0.0d0
    particles_cur_accel = 0.0d0
    particles_prev_accel = 0.0d0
    particles_charge = 0.0d0
    particles_species = 0
    particles_mass = 0.0d0
    particles_step = 0
    particles_mask = .true.

    ramo_current = 0.0d0
    ramo_cur_prev = 0.0d0
    ramo_integral = 0.0d0
    life_time = 0

    density_map_elec = 0
    density_map_hole = 0

    V_cur = 0.0d0
    V_prev = 0.0d0

    ! Start with an empty system
    nrPart      = 0
    nrElec      = 0
    nrHole      = 0
    nrElecHole  = 0

    startElecHoles = 1
    endElecHoles   = 0


    nrPart_remove = 0
    nrElec_remove = 0
    nrHole_remove = 0

    nrPart_remove_top = 0
    nrPart_remove_bot = 0
    nrElec_remove_top = 0
    nrElec_remove_bot = 0
    nrHole_remove_top = 0
    nrHole_remove_bot = 0

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

    call RANDOM_SEED(size = n)
    allocate(my_seed(n))
    my_seed = SEED
    call RANDOM_SEED(PUT = my_seed)
    deallocate(my_seed)

    ! Create folder for output files
    call execute_command_line ('mkdir -p out/')


    ! Open data file for writing
    !open(newunit=ud_pos, iostat=IFAIL, file='position.dt', status='REPLACE', action='write')
    !if (IFAIL /= 0) then
    !  print *, 'Vacuum: Failed to open file position.dt. ABORTING'
    !  stop
    !end if

    ! open(newunit=ud_vel, iostat=IFAIL, file='velocity.dt', status='REPLACE', action='write')
    ! if (IFAIL /= 0) then
    !   print *, 'Vacuum: Failed to open file velocity.dt. ABORTING'
    !   stop
    ! end if

    open(newunit=ud_emit, iostat=IFAIL, file='out/emitted.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file emitted.dt. ABORTING'
      stop
    end if

    open(newunit=ud_absorb, iostat=IFAIL, file='out/absorbed.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file absorbed.dt. ABORTING'
      stop
    end if

    open(newunit=ud_absorb_top, iostat=IFAIL, file='out/absorbed_top.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file absorbed_top.dt. ABORTING'
      stop
    end if

    open(newunit=ud_absorb_bot, iostat=IFAIL, file='out/absorbed_bot.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file absorbed_bot.dt. ABORTING'
      stop
    end if

    open(newunit=ud_ramo, iostat=IFAIL, file='out/ramo_current.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file ramo_current.dt. ABORTING'
      stop
    end if

    open(newunit=ud_volt, iostat=IFAIL, file='out/volt.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file volt.dt. ABORTING'
      stop
    end if

    open(newunit=ud_dipole_pos, iostat=IFAIL, file='out/dipole_pos.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file dipole_pos.dt. ABORTING'
      stop
    end if

    open(newunit=ud_dipole_vec, iostat=IFAIL, file='out/dipole_vec.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file dipole_vec.dt. ABORTING'
      stop
    end if

    open(newunit=ud_field, iostat=IFAIL, file='out/long_field.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file long_field.dt. ABORTING'
      stop
    end if

    open(newunit=ud_density_map_elec, iostat=IFAIL, file='out/density_elec.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file density_elec.dt. ABORTING'
      stop
    end if

    open(newunit=ud_density_map_hole, iostat=IFAIL, file='out/density_hole.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file density_hole.dt. ABORTING'
      stop
    end if

    open(newunit=ud_density_map_total, iostat=IFAIL, file='out/density_total.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file density_total.dt. ABORTING'
      stop
    end if
  end subroutine Init

  ! ----------------------------------------------------------------------------
  ! Flush data written to files such that it can read
  subroutine Flush_Data()

    flush(ud_emit)
    flush(ud_absorb)
    flush(ud_volt)

  end subroutine Flush_Data


  ! ----------------------------------------------------------------------------
  ! Write out parameters, constants and initial variables used when the program starts
  subroutine Write_Initial_Variables()
    implicit none
    integer :: ud_init, IFAIL
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
    write(ud_init, fmt_int) 'SEED                = ', SEED,          'Seed value used in the random number generator'
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

    ! Close file 'init.dt'
    close(unit=ud_init, iostat=IFAIL, status='keep')
  end subroutine Write_Initial_Variables


  ! ----------------------------------------------------------------------------
  ! Clean up after the program
  ! Deallocate variables, close files, etc.
  subroutine Clean_up()
    integer :: IFAIL

    ! Close file descriptors
    close(unit=ud_pos, iostat=IFAIL, status='keep')
    !close(unit=ud_vel, iostat=IFAIL, status='keep')
    close(unit=ud_emit, iostat=IFAIL, status='keep')
    close(unit=ud_absorb, iostat=IFAIL, status='keep')
    close(unit=ud_absorb_top, iostat=IFAIL, status='keep')
    close(unit=ud_absorb_bot, iostat=IFAIL, status='keep')
    close(unit=ud_ramo, iostat=IFAIL, status='keep')
    close(unit=ud_volt, iostat=IFAIL, status='keep')
    close(unit=ud_dipole_pos, iostat=IFAIL, status='keep')
    close(unit=ud_dipole_vec, iostat=IFAIL, status='keep')
    close(unit=ud_field, iostat=IFAIL, status='keep')
    close(unit=ud_density_map_elec, iostat=IFAIL, status='keep')
    close(unit=ud_density_map_hole, iostat=IFAIL, status='keep')
    close(unit=ud_density_map_total, iostat=IFAIL, status='keep')

    ! Deallocate arrays
    deallocate(particles_cur_pos)
    deallocate(particles_prev_pos)
    deallocate(particles_cur_vel)
    !deallocate(particles_prev_vel)
    deallocate(particles_cur_accel)
    deallocate(particles_prev_accel)
    deallocate(particles_charge)
    deallocate(particles_species)
    deallocate(particles_mass)
    deallocate(particles_mask)

    deallocate(emitters_pos)
    deallocate(emitters_dim)
    deallocate(emitters_type)
    deallocate(emitters_delay)

    deallocate(ramo_current)
    deallocate(life_time)

    deallocate(density_map_elec)
    deallocate(density_map_hole)

    !Dipole arrays
    !deallocate(dipoles_theta)
    !deallocate(dipoles_omega)
    !deallocate(dipoles_alpha)
    !deallocate(dipoles_center)

  end subroutine Clean_up
end program VacuumMD
