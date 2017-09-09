program VacuumMD
#if defined(_OPENMP)
  use omp_lib
#endif
#if defined(__INTEL_COMPILER)
  USE IFPORT, only: GETPID
#endif
  use mod_global
  use mod_verlet
  use mod_pair
  use mod_unit_tests
  implicit none

  integer :: i, nthreads
  integer, dimension(1:9) :: progress

#if defined(_OPENMP)
  nthreads = omp_get_max_threads()
#else
  nthreads = 1
#endif

#if defined(_UNIT_TEST_)
  print '(a)', 'Vacuum: **** UNIT TESTING IS ACTIVE ****'
#endif

  print '(a)', 'Vacuum: Reading input values'
  call Read_Input_Variables()

  print '(a)', 'Vacuum: Initialzing'
  call Init()

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

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, tid)

#if defined(_UNIT_TEST_)
  call Run_Unit_Tests()
  !$OMP BARRIER
  stop
#endif

#if defined(_OPENMP)
  tid = omp_get_thread_num()
#else
  tid = 0
#endif
  !$OMP SINGLE
  print '(tr1, a, i0)', 'Number of threads ', nthreads
  print *, ''

  at_step = 0
  cur_time = 0
  !$OMP END SINGLE

  do i = 1, steps

    ! Do Emission
    !call emission code

    ! Update the position of all particles
    !print *, 'Update position'
    call Update_Position(i)

    ! Remove particles from the system
    call Remove_Particles(i)

    !$OMP MASTER
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

    !$OMP END MASTER

    !$OMP BARRIER
  end do

  !$OMP MASTER
  call PrintProgress(10)
  print '(a)', 'Vacuum: Main loop finished'
  !$OMP END MASTER

  !$OMP END PARALLEL

  print '(a)', 'Vacuum: Writing data'
  call Write_Life_Time()

  print '(a)', 'Vacuum: Program finished'
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

    ! Read the input file
    read(ud_input, NML=input)

    ! Close the 'input' file
    close(unit=ud_input, iostat=IFAIL, status='keep')

    ! V: Voltage over the gap
    ! box_dim: Dimensions of the system
    ! d: Gap spacing
    box_dim = box_dim * length_scale
    d = box_dim(3)
    E_z = -1.0d0*V/d
    V_a = V

    E_zunit = -1.0d0/d

    ! Set the dimensions of the box used
    !box_dim(1) = d
    !box_dim(2) = d
    !box_dim(3) = d

    dens_x_d = box_dim(1) / (N_x_densmap-1)
    dens_y_d = box_dim(2) / (N_y_densmap-1)

    ! Set the current scale to mA / cm^2
    cur_scale = 1.0d-3 * (box_dim(1)*1.0d2) * (box_dim(3)*1.0d2) ! mA / cm^2

    !time_step: The size of the time step
    time_step = time_step * time_scale
    time_step2 = time_step**2 ! Time_step squared

    ! steps: Number of time steps that the program will simulate
    if (steps <= 0) then
      print '(a)', 'ERROR: steps <= 0'
      stop
    end if

    !open(newunit=ud_debug, iostat=IFAIL, file='debug.dt', status='replace', action='write')
    !if (IFAIL /= 0) then
    !  print '(a)', 'Vacuum: ERROR UNABLE TO OPEN file debug.dt'
    !  stop
    !end if

    !write(ud_debug, NML=input_test)

    !close(unit=ud_debug, iostat=IFAIL, status='keep')

  end subroutine Read_Input_Variables

  ! ----------------------------------------------------------------------------
  ! Initialize the program
  ! allocate and initilize variables, open data files for writing, etc.
  subroutine Init()
    integer :: IFAIL

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
    life_time = 0

    density_map_elec = 0
    density_map_hole = 0

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


    ! Open data file for writing
    open(newunit=ud_pos, iostat=IFAIL, file='position.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file position.dt. ABORTING'
      stop
    end if

    ! open(newunit=ud_vel, iostat=IFAIL, file='velocity.dt', status='REPLACE', action='write')
    ! if (IFAIL /= 0) then
    !   print *, 'Vacuum: Failed to open file velocity.dt. ABORTING'
    !   stop
    ! end if

    open(newunit=ud_emit, iostat=IFAIL, file='emitted.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file emitted.dt. ABORTING'
      stop
    end if

    open(newunit=ud_absorb, iostat=IFAIL, file='absorbed.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file absorbed.dt. ABORTING'
      stop
    end if

    open(newunit=ud_absorb_top, iostat=IFAIL, file='absorbed_top.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file absorbed_top.dt. ABORTING'
      stop
    end if

    open(newunit=ud_absorb_bot, iostat=IFAIL, file='absorbed_bot.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file absorbed_bot.dt. ABORTING'
      stop
    end if

    open(newunit=ud_ramo, iostat=IFAIL, file='ramo_current.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file ramo_current.dt. ABORTING'
      stop
    end if

    open(newunit=ud_volt, iostat=IFAIL, file='volt.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file volt.dt. ABORTING'
      stop
    end if

    open(newunit=ud_dipole_pos, iostat=IFAIL, file='dipole_pos.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file dipole_pos.dt. ABORTING'
      stop
    end if

    open(newunit=ud_dipole_vec, iostat=IFAIL, file='dipole_vec.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file dipole_vec.dt. ABORTING'
      stop
    end if

    open(newunit=ud_field, iostat=IFAIL, file='long_field.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file long_field.dt. ABORTING'
      stop
    end if

    open(newunit=ud_density_map_elec, iostat=IFAIL, file='density_elec.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file density_elec.dt. ABORTING'
      stop
    end if

    open(newunit=ud_density_map_hole, iostat=IFAIL, file='density_hole.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file density_hole.dt. ABORTING'
      stop
    end if

    open(newunit=ud_density_map_total, iostat=IFAIL, file='density_total.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file density_total.dt. ABORTING'
      stop
    end if
  end subroutine Init


  ! ----------------------------------------------------------------------------
  ! Write out parameters, constants and initial variables used when the program starts
  subroutine Write_Initial_Variables()
    implicit none
    integer :: ud_init, IFAIL
    character(len=*), parameter :: fmt_int = "(a, tr8, i8, tr6, a)"
    character(len=*), parameter :: fmt_rel = "(a, tr2, ES16.6, tr2, a)"
    !character(len=*), parameter :: fmt_cmp = "(a, tr2, '(', G16.6, ',', G16.6, ')', tr2, a)"

    ! Open file 'init.dt' for writing
    open(newunit=ud_init, iostat=IFAIL, file='init.dt', status='REPLACE', action='write')
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
    write(ud_init, fmt_rel) 'lambda_pos = ', lambda_pos,'Mean value in the Poisson distribution'
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
    write(ud_init, fmt_rel) 'V             = ', V,          'The Voltage in the system'
    write(ud_init, fmt_rel) 'box_dim_x     = ', box_dim(1), 'System x-length'
    write(ud_init, fmt_rel) 'box_dim_y     = ', box_dim(2), 'System y-length'
    write(ud_init, fmt_rel) 'box_dim_z     = ', box_dim(3), 'System z-length'
    write(ud_init, fmt_rel) 'd             = ', d,          'Gap spacing in the system'
    write(ud_init, fmt_rel) 'E_z           = ', E_z,        'Y-component of the electric field in the system'
    write(ud_init, *) '---------------------------------------------------------'
    !write(ud_init, *)
    !write(ud_init, *)

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
