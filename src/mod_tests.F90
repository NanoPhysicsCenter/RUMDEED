!-------------------------------------------!
! Module for unit and integration tests     !
! Kristinn Torfason                         !
! 10.09.17                                  !
!-------------------------------------------!

module mod_tests
  use mod_global
  use mod_pair
  use mod_verlet
  use mod_field_emission_v2
  implicit none

  PRIVATE
  PUBLIC :: Init_Test, Clean_Up_Test, Run_Tests, tests_passed, tests_failed

  ! ----------------------------------------------------------------------------
  ! Variables
  integer, dimension(:), allocatable :: nrEmitted_emitters
  integer                            :: nrElecEmitAll

  ! Test result counters. These are made public so that the main program can
  ! query them after the tests have run and set a non-zero exit code on failure.
  integer :: tests_passed = 0
  integer :: tests_failed = 0

  ! Constant used in MC integration
  double precision :: time_step_div_q0

  ! Passed / failed color strings
  character(len=16), parameter :: passed_str = achar(27)//'[32m PASSED'//achar(27)//'[0m'
  character(len=16), parameter :: failed_str = achar(27)//'[31m FAILED'//achar(27)//'[0m'

contains
    !-----------------------------------------------------------------------------
    ! Initialize the Field Emission
  subroutine Init_Test()
    ! Allocate the number of emitters
    nrEmit = 1
    allocate(nrEmitted_emitters(1:nrEmit))

    ! Initialize variables
    nrEmitted_emitters = 0 ! Number of electrons emitted from emitter number i in this time step

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_Planar

    ! Function for the electric field in the system
    ptr_field_E => field_E_planar

    ! The function that does the emission
    ptr_Do_Emission => Do_Test_Emission

    ! The function to do image charge effects
    ptr_Image_Charge_effect => Force_Image_charges_v2

    ! The unit electric field used by the Verlet integrator (Ramo current).
    ! Without this the integrator dereferences a null procedure pointer in the
    ! transit time test.
    ptr_E_zunit => E_zunit_planar

    ! Parameters used in the module
    time_step_div_q0 = time_step / q_0

    ! Initialize the Ziggurat algorithm
    !call zigset(my_seed(1))

  end subroutine Init_Test

  subroutine Clean_Up_Test()
    deallocate(nrEmitted_emitters)
  end subroutine Clean_Up_Test

  !-----------------------------------------------------------------------------
  ! Shared assertion helpers
  !
  ! These keep track of how many checks have passed or failed in the module
  ! level counters tests_passed / tests_failed. The main program inspects
  ! tests_failed after the tests have run and sets a non-zero exit code if any
  ! check failed. This is what allows the test suite to gate automation (CI).
  !-----------------------------------------------------------------------------

  ! Assert that a logical condition is true.
  subroutine Assert_True(cond, label)
    logical, intent(in)          :: cond
    character(len=*), intent(in) :: label

    if (cond .eqv. .true.) then
      tests_passed = tests_passed + 1
      print *, trim(label)//passed_str
    else
      tests_failed = tests_failed + 1
      print *, trim(label)//failed_str
    end if
  end subroutine Assert_True

  ! Assert that a scalar value is close to an expected value.
  ! Uses a relative tolerance (tolerance_rel), falling back to an absolute
  ! tolerance (tolerance_abs) when the expected value is (near) zero.
  subroutine Assert_Close(actual, expected, label)
    double precision, intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    logical                      :: ok

    if (abs(expected) < tolerance_abs) then
      ok = (abs(actual - expected) < tolerance_abs)
    else
      ok = (abs(actual - expected)/abs(expected) < tolerance_rel)
    end if

    if (ok .eqv. .true.) then
      tests_passed = tests_passed + 1
      print *, trim(label)//passed_str
    else
      tests_failed = tests_failed + 1
      print *, trim(label)//failed_str
      print *, '  actual   = ', actual
      print *, '  expected = ', expected
    end if
  end subroutine Assert_Close

  ! Assert that a 3-vector is close to an expected 3-vector.
  ! Each component is compared with a relative tolerance (tolerance_rel),
  ! falling back to an absolute tolerance (tolerance_abs) for (near) zero
  ! components. This matches the component-wise idiom used throughout the file.
  subroutine Assert_Close_Vec(actual, expected, label)
    double precision, dimension(1:3), intent(in) :: actual, expected
    character(len=*), intent(in)                 :: label
    logical                                      :: ok
    integer                                      :: k

    ok = .true.
    do k = 1, 3
      if (abs(expected(k)) < tolerance_abs) then
        if (abs(actual(k) - expected(k)) >= tolerance_abs) ok = .false.
      else
        if (abs(actual(k) - expected(k))/abs(expected(k)) >= tolerance_rel) ok = .false.
      end if
    end do

    if (ok .eqv. .true.) then
      tests_passed = tests_passed + 1
      print *, trim(label)//passed_str
    else
      tests_failed = tests_failed + 1
      print *, trim(label)//failed_str
      print *, '  actual   = ', actual
      print *, '  expected = ', expected
    end if
  end subroutine Assert_Close_Vec

    !-----------------------------------------------------------------------------
    ! This subroutine gets called from main when the emitters should emit the electrons
  subroutine Do_Test_Emission(step)
    integer, intent(in)              :: step
    integer                          :: i, IFAIL = 0, emit_step, species
    double precision                 :: cur_time, x, y, z
    integer                          :: ud_file
    double precision, dimension(1:3) :: par_pos, par_vel
    double precision                 :: x_0, y_0, z_0
    double precision                 :: x_1, y_1, z_1

    nrElecEmitAll = 0
    nrEmitted_emitters = 0

    ! open(newunit=ud_file, iostat=IFAIL, file='Test_1.bin', status='OLD', action='READ', access='STREAM')

    ! do
    !   read(unit=ud_file, iostat=IFAIL) x, y, z, emit_step, species
    !   if (IFAIL == 0) then
    !     if (emit_step == step) then
    !       par_pos(1) = x
    !       par_pos(2) = y
    !       par_pos(3) = z
    !       par_vel = 0.0d0

    !       ! Add a particle to the system
    !       call Add_Particle(par_pos, par_vel, species, step, 1, -1, 1)
    !     end if
    !   else
    !     exit
    !   end if
    ! end do

    ! close(unit=ud_file, status='keep')

    ! P1
    x_0 = 5.5d0*length_scale
    y_0 = -70.4d0*length_scale
    z_0 = 36.2d0*length_scale

    par_pos(1) = x_0
    par_pos(2) = y_0
    par_pos(3) = z_0
    par_vel = 0.0d0
    species = species_elec

    call Add_Particle(par_pos, par_vel, species, step, 1, -1, 1)

    ! P2
    x_1 = -33.4d0*length_scale
    y_1 = 4.0d0*length_scale
    z_1 = 45.3d0*length_scale

    par_pos(1) = x_1
    par_pos(2) = y_1
    par_pos(3) = z_1
    par_vel = 0.0d0
    species = species_elec

    call Add_Particle(par_pos, par_vel, species, step, 1, -1, 1)

    cur_time = time_step * step / time_scale ! Scaled in units of time_scale
    write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, nrElecEmitAll, nrElec, &
                                                      & (nrEmitted_emitters(i), i = 1, nrEmit)
  end subroutine Do_Test_Emission

  subroutine Run_Tests()
    logical :: fe_v_y_ok

    ! Reset the counters in case the tests are run more than once.
    tests_passed = 0
    tests_failed = 0

    print *, ''
    call Test_Field_Emission_Module(fe_v_y_ok)
    call Assert_True(fe_v_y_ok, 'Field emission v_y')

    print *, ''
    call Test_Acceleration_Without_Image_Charge()

    print *, ''
    call Test_Image_Charge()

    !call Test_Acceleration_With_Image_Charge()

    print *, ''
    call Test_Particle_Bookkeeping()

    print *, ''
    call Test_Planar_Field()

    print *, ''
    call Test_Transit_Time()

    !call Test_Many_Particles()

    print *, ''
    call Test_Tip_Acceleration()

#if defined(_OPENMP)
    print '(a)', 'Do OpenMP test'
#else
    print '(a)', 'Do Serial test'
#endif

    ! Print a summary of the results.
    print *, ''
    print '(a)', '----------------------------------------'
    print '(a, i0, a, i0, a)', 'RUMDEED: Test summary: ', tests_passed, ' passed, ', &
                             & tests_failed, ' failed'
    print '(a, i0)', 'RUMDEED: Total checks: ', (tests_passed + tests_failed)
    if (tests_failed > 0) then
      print *, achar(27)//'[31mRUMDEED: ONE OR MORE TESTS FAILED'//achar(27)//'[0m'
    else
      print *, achar(27)//'[32mRUMDEED: ALL TESTS PASSED'//achar(27)//'[0m'
    end if
    print '(a)', '----------------------------------------'
  end subroutine Run_Tests

  !-----------------------------------------------------------------------------
  ! Initialize an empty system with the given parameters for testing
  subroutine Setup_Test_System(d_test, delta_t_test, steps_test, V_test, use_ic, N_ic)
    double precision, intent(in) :: d_test, delta_t_test, V_test
    integer, intent(in)          :: steps_test, N_ic
    logical, intent(in)          :: use_ic

    ! Clear variables
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
    life_time = 0

    ! Start with an empty system
    nrPart     = 0
    nrElec     = 0
    nrIon      = 0
    nrElecIon  = 0

    nrPart_remove = 0
    nrElec_remove = 0
    nrIon_remove  = 0

    nrPart_remove_top = 0
    nrPart_remove_bot = 0
    nrElec_remove_top = 0
    nrElec_remove_bot = 0
    nrIon_remove_top = 0
    nrIon_remove_bot = 0

    ! ID starts with 0
    nrID = 0

    ! Set input variables
    box_dim = (/ 100.0d0*length_scale, 100.0d0*length_scale, d_test /)
    time_step = delta_t_test
    steps = steps_test
    V_s = V_test

    ! Init
    d = box_dim(3)
    V_d = V_s
    E_z = -1.0d0*V_d/d
    !E_zunit = -1.0d0/d ! Needs to be updated

    ! Temperature, pressure and density
    P_abs = P_abs * P_ntp 
    n_d = P_abs/(k_b*T_temp)

    !time_step: The size of the time step
    !time_step = time_step * time_scale
    time_step2 = time_step**2 ! Time_step squared

    image_charge = use_ic
    N_ic_max = N_ic
  end subroutine Setup_Test_System

  subroutine Test_Acceleration_Without_Image_Charge()
    double precision, dimension(1:3) :: R_1, R_1a, R_1b
    double precision, dimension(1:3) :: R_2, R_2a, R_2b
    double precision, dimension(1:3) :: R_3, R_3a, R_3b
    double precision, dimension(1:3) :: E_test, par_vel, E_pos, E_python, E_F

    double precision, dimension(1:3) :: a_1, a_12, a_13, a_1E, a_1_res
    double precision, dimension(1:3) :: a_2, a_21, a_23, a_2E, a_2_res
    double precision, dimension(1:3) :: a_3, a_31, a_32, a_3E, a_3_res

    double precision, parameter      :: d_test = 100.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d0 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step

    double precision, parameter      :: pre_fac_c_test = q_0**2/(4.0d0*pi*epsilon_0*epsilon_r)
    !double precision, parameter      :: pre_fac_E_test = q_0/m_0

    print *, 'Running test for acceleration without image charge'

    ! Electric field in the system (Planar case)
    E_test = (/ 0.0d0, 0.0d0, -1.0d0*V_test/d_test /)

    !--------------------------------------------------------------------------
    ! Location of particles

    ! Particle 1 with negative charge
    R_1 = (/  3.0d0, -10.0d0, 2.0d0 /) * length_scale 

    ! Particle 2 with negative charge
    R_2 = (/ -9.0d0,  26.0d0, 80.0d0 /) * length_scale

    ! Particle 3 with positive charge
    R_3 = (/  6.0d0, -24.0d0, 56.53d0 /) * length_scale

    !--------------------------------------------------------------------------
    ! Calculate the acceleration
    ! Particle 1
    ! Acceleration from other particles
    a_12  = +1.0d0*pre_fac_c_test/m_0 * (R_1 - R_2) / norm2(R_1 - R_2)**3
    a_13  = -1.0d0*pre_fac_c_test/m_0 * (R_1 - R_3) / norm2(R_1 - R_3)**3

    ! Acceleration from the electric field in the diode
    a_1E = -1.0d0*q_0/m_0 * E_test

    ! Total acceleration 1
    a_1 = a_12 + a_13 + a_1E

    ! Particle 2
    ! Acceleration from other particles
    a_21  = +1.0d0*pre_fac_c_test/m_0 * (R_2 - R_1) / norm2(R_2 - R_1)**3
    a_23  = -1.0d0*pre_fac_c_test/m_0 * (R_2 - R_3) / norm2(R_2 - R_3)**3

    ! Acceleration from the electric field in the diode
    a_2E = -1.0d0*q_0/m_0 * E_test

    ! Total aceleration on particle 2
    a_2 = a_21 + a_23 + a_2E

    ! Particle 3
    ! Acceleration from other particles
    a_31  = -1.0d0*pre_fac_c_test/m_N2p * (R_3 - R_1) / norm2(R_3 - R_1)**3
    a_32  = -1.0d0*pre_fac_c_test/m_N2p * (R_3 - R_2) / norm2(R_3 - R_2)**3

    ! Acceleration from the electric field in the diode
    a_3E = +1.0d0*q_0/m_N2p * E_test

    ! Total acceleration on particle 3
    a_3 = a_31 + a_32 + a_3E

    ! Set input variables
    call Setup_Test_System(d_test, delta_t_test, 100, V_test, .false., 0)

    ! Test the field at this location
    par_vel = (/ -4.55, -2.34, 96.44 /) * length_scale
    ! Field without any electrons in the system
    E_F = Calc_Field_at(par_vel)

    ! Add particles
    par_vel = 0.0d0
    call Add_Particle(R_1, par_vel, species_elec, 1, 0, -1)
    call Add_Particle(R_2, par_vel, species_elec, 1, 0, -1)
    call Add_Particle(R_3, par_vel, species_ion, 1, 0, -1)

    call Calculate_Acceleration_Particles()

    ! Test the field at this location
    par_vel = (/ -4.55, -2.34, 96.44 /) * length_scale
    E_pos = Calc_Field_at(par_vel)

    ! Results from the acceleration
    a_1_res = particles_cur_accel(:, 1)
    a_2_res = particles_cur_accel(:, 2)
    a_3_res = particles_cur_accel(:, 3)

    call Assert_Close_Vec(a_1_res, a_1, 'Acceleration (no IC) particle 1')
    call Assert_Close_Vec(a_2_res, a_2, 'Acceleration (no IC) particle 2')
    call Assert_Close_Vec(a_3_res, a_3, 'Acceleration (no IC) particle 3')

    call Assert_Close_Vec(E_F, E_test, 'Electric field without particles in planar case')

    E_python = (/ -314559.29097098, 1423979.07058996, -20246038.87978313 /) ! Results from Python script (no image charge)
    call Assert_Close_Vec(E_pos, E_python, 'Electric field with particles')

    print *, 'Acceleration test without image charge finished'
    print *, ''
  end subroutine Test_Acceleration_Without_Image_Charge

  !-----------------------------------------------------------------------------
  ! Test the function Force_Image_charges_v2(pos_1, pos_2)
  subroutine Test_Image_Charge()
    ! Parameters for the test system
    double precision, parameter      :: d_test = 100.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d0 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step

    ! Particles
    double precision, parameter      :: q_1 = -1.0d0*q_0, q_2 = +1.0d0*q_0
    double precision, parameter      :: q_1_ic = -1.0d0*q_1, q_2_ic = -1.0d0*q_2
    double precision, dimension(1:3) :: R_1, R_2, R_1_IC, R_2_IC
    double precision, dimension(1:3) :: R_1_IC_1, R_1_IC_2, R_1_IC_3, R_1_IC_4, R_1_IC_5
    double precision, dimension(1:3) :: R_2_IC_1, R_2_IC_2, R_2_IC_3, R_2_IC_4, R_2_IC_5
    double precision                 :: q_1_IC_1, q_1_IC_2, q_1_IC_3, q_1_IC_4, q_1_IC_5
    double precision                 :: q_2_IC_1, q_2_IC_2, q_2_IC_3, q_2_IC_4, q_2_IC_5

    ! Calculation results
    double precision, dimension(1:3) :: F_R1, F_R2, F_IC_11, F_IC_12, F_IC_21, F_IC_13, F_IC_14, F_IC_15
    double precision, dimension(1:3) :: force_ic_11, force_ic_22, force_ic_12, force_ic_21

    print *, 'Testing image charge function "Force_Image_charges_v2(pos_1, pos_2)"'

    ! Set input variables
    ! The zero at the end means take into account only 1 image charge partner
    call Setup_Test_System(d_test, delta_t_test, 100, V_test, .True., 0)

    ! We have two particles, an electron at,
    R_1 = (/  3.0d0, -10.0d0, 2.0d0 /) * length_scale
    ! and ion at,
    R_2 = (/  6.0d0, -24.0d0, 98.53d0 /) * length_scale

    ! Let's test the self interaction first
    force_ic_11 = q_1**2 * div_fac_c * Force_Image_charges_v2(R_1, R_1)

    ! The force should be given by Coulombs law
    ! F = 1/(4*pi*epsilon_0) * (r - r')/|r - r'|^3
    ! Here r is the location of the particle and r' is the location of the image charge
    ! The image charge is located below emitter. The force should only be in the z-direction,
    ! F_R1 = 1/(4*pi*epsilon_0) * (-q_1)*(+q_1)*(z - z')/|z - z'|^3 .
    ! We know that z' = -z so z - z' = z - (-z) = 2z, or,
    ! F_R1 = -q_1**2/(4*pi*epsilon_0) * 2z/(2z)^3 = -q_1**2/(4*pi*epsilon_0) * 1/(2z)^2

    F_R1 = 0.0d0
    F_R1(3) = -1.0d0*q_1**2/(4.0d0*pi*epsilon_0)*1.0d0/(2.0d0*R_1(3))**2

    call Assert_Close_Vec(force_ic_11, F_R1, 'Self interaction on particle 1')

    force_ic_22 = q_2**2 * div_fac_c * Force_Image_charges_v2(R_2, R_2)

    ! F_R2 = -q_2**2/(4*pi*epsilon_0) * 2z/(2z)^3 = -q_2**2/(4*pi*epsilon_0) * 1/(2z)^2

    F_R2 = 0.0d0
    F_R2(3) = -1.0d0*q_2**2/(4.0d0*pi*epsilon_0)*1.0d0/(2.0d0*R_2(3))**2

    call Assert_Close_Vec(force_ic_22, F_R2, 'Self interaction on particle 2')

    ! Test image charge force on particle 1 from particle 2
    ! Image charge of particle 2 is located at -z
    R_2_IC = R_2
    R_2_IC(3) = -1.0d0*R_2_Ic(3)

    F_IC_12 = q_1*q_2_ic/(4.0d0*pi*epsilon_0) * (R_1 - R_2_IC) / norm2(R_1 - R_2_IC)**3

    force_ic_12 = q_1*q_2 * div_fac_c * Force_Image_charges_v2(R_1, R_2)

    call Assert_Close_Vec(force_ic_12, F_IC_12, 'Image charge interaction on particle 1 from particle 2')

    ! Test image charge force on particle 2 from particle 1
    ! Image charge of particle 1 is located at -z
    R_1_IC = R_1
    R_1_IC(3) = -1.0d0*R_1_Ic(3)

    F_IC_21 = q_2*q_1_ic/(4.0d0*pi*epsilon_0) * (R_2 - R_1_IC) / norm2(R_2 - R_1_IC)**3

    force_ic_21 = q_2*q_1 * div_fac_c * Force_Image_charges_v2(R_2, R_1)

    call Assert_Close_Vec(force_ic_21, F_IC_21, 'Image charge interaction on particle 2 from particle 1')

    !-------------------------------------------------------------------------------------------
    ! Test 5 image charge partners of each particle
    call Setup_Test_System(d_test, delta_t_test, 100, V_test, .True., 1)

    ! Test self interaction of particle 1 (-)
    ! n =  0: z = -z_0 (+)
    ! n = -1: z = -2d-z_0 (+), z = -2d+z_0 (-)
    ! n = +1: z =  2d-z_0 (+), z =  2d+z_0 (-)

    R_1_IC_1 = R_1
    R_1_IC_1(3) = -1.0d0*R_1_IC_1(3)
    q_1_IC_1 = +1.0d0*q_0

    R_1_IC_2 = R_1
    R_1_IC_2(3) = -2.0d0*d_test - R_1_IC_2(3)
    q_1_IC_2 = +1.0d0*q_0

    R_1_IC_3 = R_1
    R_1_IC_3(3) = -2.0d0*d_test + R_1_IC_3(3)
    q_1_IC_3 = -1.0d0*q_0

    R_1_IC_4 = R_1
    R_1_IC_4(3) = 2.0d0*d_test - R_1_IC_4(3)
    q_1_IC_4 = +1.0d0*q_0
    
    R_1_IC_5 = R_1
    R_1_IC_5(3) = 2.0d0*d_test + R_1_IC_5(3)
    q_1_IC_5 = -1.0d0*q_0

    F_IC_11 = q_1*q_1_IC_1/(4.0d0*pi*epsilon_0) * (R_1 - R_1_IC_1) / norm2(R_1 - R_1_IC_1)**3
    F_IC_12 = q_1*q_1_IC_2/(4.0d0*pi*epsilon_0) * (R_1 - R_1_IC_2) / norm2(R_1 - R_1_IC_2)**3
    F_IC_13 = q_1*q_1_IC_3/(4.0d0*pi*epsilon_0) * (R_1 - R_1_IC_3) / norm2(R_1 - R_1_IC_3)**3
    F_IC_14 = q_1*q_1_IC_4/(4.0d0*pi*epsilon_0) * (R_1 - R_1_IC_4) / norm2(R_1 - R_1_IC_4)**3
    F_IC_15 = q_1*q_1_IC_5/(4.0d0*pi*epsilon_0) * (R_1 - R_1_IC_5) / norm2(R_1 - R_1_IC_5)**3

    F_R1 = F_IC_11 + F_IC_12 + F_IC_13 + F_IC_14 + F_IC_15

    force_ic_11 = q_1**2 * div_fac_c * Force_Image_charges_v2(R_1, R_1)

    call Assert_Close_Vec(force_ic_11, F_R1, 'Self interaction on particle 1 (N_IC_MAX = 1)')

    ! Test self interaction of particle 2 (+)
    ! n =  0: z = -z_0 (-)
    ! n = -1: z = -2d-z_0 (-), z = -2d+z_0 (+)
    ! n = +1: z =  2d-z_0 (-), z =  2d+z_0 (+)

    R_2_IC_1 = R_2
    R_2_IC_1(3) = -1.0d0*R_2_IC_1(3)
    q_2_IC_1 = -1.0d0*q_0

    R_2_IC_2 = R_2
    R_2_IC_2(3) = -2.0d0*d_test - R_2_IC_2(3)
    q_2_IC_2 = -1.0d0*q_0

    R_2_IC_3 = R_2
    R_2_IC_3(3) = -2.0d0*d_test + R_2_IC_3(3)
    q_2_IC_3 = +1.0d0*q_0

    R_2_IC_4 = R_2
    R_2_IC_4(3) = 2.0d0*d_test - R_2_IC_4(3)
    q_2_IC_4 = -1.0d0*q_0
    
    R_2_IC_5 = R_2
    R_2_IC_5(3) = 2.0d0*d_test + R_2_IC_5(3)
    q_2_IC_5 = +1.0d0*q_0

    F_IC_11 = q_2*q_2_IC_1/(4.0d0*pi*epsilon_0) * (R_2 - R_2_IC_1) / norm2(R_2 - R_2_IC_1)**3
    F_IC_12 = q_2*q_2_IC_2/(4.0d0*pi*epsilon_0) * (R_2 - R_2_IC_2) / norm2(R_2 - R_2_IC_2)**3
    F_IC_13 = q_2*q_2_IC_3/(4.0d0*pi*epsilon_0) * (R_2 - R_2_IC_3) / norm2(R_2 - R_2_IC_3)**3
    F_IC_14 = q_2*q_2_IC_4/(4.0d0*pi*epsilon_0) * (R_2 - R_2_IC_4) / norm2(R_2 - R_2_IC_4)**3
    F_IC_15 = q_2*q_2_IC_5/(4.0d0*pi*epsilon_0) * (R_2 - R_2_IC_5) / norm2(R_2 - R_2_IC_5)**3

    F_R2 = F_IC_11 + F_IC_12 + F_IC_13 + F_IC_14 + F_IC_15

    force_ic_22 = q_2**2 * div_fac_c * Force_Image_charges_v2(R_2, R_2)

    call Assert_Close_Vec(force_ic_22, F_R2, 'Self interaction on particle 2 (N_IC_MAX = 1)')

    !--------------------------------------------------------------------------------------
    ! Testing interaction on particle 1 from particle 2 with N_IC_MAX = 1

    F_IC_11 = q_1*q_2_IC_1/(4.0d0*pi*epsilon_0) * (R_1 - R_2_IC_1) / norm2(R_1 - R_2_IC_1)**3
    F_IC_12 = q_1*q_2_IC_2/(4.0d0*pi*epsilon_0) * (R_1 - R_2_IC_2) / norm2(R_1 - R_2_IC_2)**3
    F_IC_13 = q_1*q_2_IC_3/(4.0d0*pi*epsilon_0) * (R_1 - R_2_IC_3) / norm2(R_1 - R_2_IC_3)**3
    F_IC_14 = q_1*q_2_IC_4/(4.0d0*pi*epsilon_0) * (R_1 - R_2_IC_4) / norm2(R_1 - R_2_IC_4)**3
    F_IC_15 = q_1*q_2_IC_5/(4.0d0*pi*epsilon_0) * (R_1 - R_2_IC_5) / norm2(R_1 - R_2_IC_5)**3

    F_R1 = F_IC_11 + F_IC_12 + F_IC_13 + F_IC_14 + F_IC_15

    force_ic_12 = q_1*q_2 * div_fac_c * Force_Image_charges_v2(R_1, R_2)

    call Assert_Close_Vec(force_ic_12, F_R1, 'Image charge interaction on particle 1 from particle 2 (N_IC_MAX = 1)')

    print *, 'Image charge function "Force_Image_charges_v2(pos_1, pos_2)" test finished'
    print *, ''
  end subroutine Test_Image_Charge

  !-----------------------------------------------------------------------------
  ! Test the Acceleration calculations
  ! Two electrons and one ion are placed in the system
  ! The acceleration is calculated using Coulomb's law and also with the
  ! Calculate Acceleration subroutine.
  subroutine Test_Acceleration_With_Image_Charge()
    double precision, dimension(1:3) :: R_1, R_1a, R_1b
    double precision, dimension(1:3) :: R_2, R_2a, R_2b
    double precision, dimension(1:3) :: R_3, R_3a, R_3b
    double precision, dimension(1:3) :: E_test, par_vel, E_pos, E_python

    double precision, dimension(1:3) :: a_1, a_12, a_13, a_1E, a_1_res
    double precision, dimension(1:3) :: a_2, a_21, a_23, a_2E, a_2_res
    double precision, dimension(1:3) :: a_3, a_31, a_32, a_3E, a_3_res
    double precision, dimension(1:3) :: a_11a, a_11b, a_12a, a_12b, a_13a, a_13b
    double precision, dimension(1:3) :: a_21a, a_21b, a_22a, a_22b, a_23a, a_23b
    double precision, dimension(1:3) :: a_31a, a_31b, a_32a, a_32b, a_33a, a_33b

    double precision, parameter      :: d_test = 100.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d0 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step

    double precision, parameter      :: pre_fac_c_test = q_0**2/(4.0d0*pi*epsilon_0*epsilon_r)
    double precision, parameter      :: pre_fac_E_test = q_0/m_0

    ! We must use OMP MASTER here because the variables declared above
    ! are private for each thread and are NOT shared.
    print *, 'Running test for acceleration with image charge'

    ! Electric field in the system
    E_test = (/ 0.0d0, 0.0d0, -1.0d0*V_test/d_test /)

    !--------------------------------------------------------------------------
    ! Location of particles
    R_1 = (/  3.0d0, -10.0d0, 2.0d0 /) * length_scale
    R_2 = (/ -9.0d0,  26.0d0, 80.0d0 /) * length_scale
    R_3 = (/  6.0d0, -24.0d0, 56.53d0 /) * length_scale

    ! Location of image charge partners
    ! Particle 1
    ! Above (Charge: Positive)
    R_1a = R_1
    R_1a(3) = 2.0d0*d_test - R_1a(3)

    ! Below (Charge: Positive)
    R_1b = R_1
    R_1b(3) = -1.0d0*R_1b(3)

    ! Particle 2
    ! Above (Charge: Positive)
    R_2a = R_2
    R_2a(3) = 2.0d0*d_test - R_2a(3)

    ! Below (Charge: Positive)
    R_2b = R_2
    R_2b(3) = -1.0d0*R_2b(3)

    ! Particle 3
    ! Above (Charge: Negative)
    R_3a = R_3
    R_3a(3) = 2.0d0*d_test - R_3a(3)

    ! Below (Charge: Negative)
    R_3b = R_3
    R_3b(3) = -1.0d0*R_3b(3)

    !--------------------------------------------------------------------------
    ! Calculate the acceleration
    ! Particle 1
    ! Acceleration from other particles
    a_12  = +1.0d0*pre_fac_c_test * (R_1 - R_2) / norm2(R_1 - R_2)**3
    a_13  = -1.0d0*pre_fac_c_test * (R_1 - R_3) / norm2(R_1 - R_3)**3

    ! Acceleration from image charge partners
    ! Self (-q*+q = -q**2)
    a_11a = -1.0d0*pre_fac_c_test * (R_1 - R_1a) / norm2(R_1 - R_1a)**3
    a_11b = -1.0d0*pre_fac_c_test * (R_1 - R_1b) / norm2(R_1 - R_1b)**3
    ! Other electron (-q*+q = -q**2)
    a_12a = -1.0d0*pre_fac_c_test * (R_1 - R_2a) / norm2(R_1 - R_2a)**3
    a_12b = -1.0d0*pre_fac_c_test * (R_1 - R_2b) / norm2(R_1 - R_2b)**3
    ! Other ion (-q*-q = +q**2)
    a_13a = +1.0d0*pre_fac_c_test * (R_1 - R_3a) / norm2(R_1 - R_3a)**3
    a_13b = +1.0d0*pre_fac_c_test * (R_1 - R_3b) / norm2(R_1 - R_3b)**3

    ! Acceleration from the electric field in the diode
    a_1E = -1.0d0*pre_fac_E_test * E_test

    a_1 = a_12 + a_13 + a_11a + a_11b + a_12a + a_12b + a_13a + a_13b + a_1E

    ! Particle 2
    ! Acceleration from other particles
    a_21  = +1.0d0*pre_fac_c_test * (R_2 - R_1) / norm2(R_2 - R_1)**3
    a_23  = -1.0d0*pre_fac_c_test * (R_2 - R_3) / norm2(R_2 - R_3)**3

    ! Acceleration from image charge partners
    ! Other electron (-q*+q = -q**2)
    a_21a = -1.0d0*pre_fac_c_test * (R_2 - R_1a) / norm2(R_2 - R_1a)**3
    a_21b = -1.0d0*pre_fac_c_test * (R_2 - R_1b) / norm2(R_2 - R_1b)**3
    ! Self (-q*+q = -q**2)
    a_22a = -1.0d0*pre_fac_c_test * (R_2 - R_2a) / norm2(R_2 - R_2a)**3
    a_22b = -1.0d0*pre_fac_c_test * (R_2 - R_2b) / norm2(R_2 - R_2b)**3
    ! Other ion (-q*-q = +q**2)
    a_23a = +1.0d0*pre_fac_c_test * (R_2 - R_3a) / norm2(R_2 - R_3a)**3
    a_23b = +1.0d0*pre_fac_c_test * (R_2 - R_3b) / norm2(R_2 - R_3b)**3

    ! Acceleration from the electric field in the diode
    a_2E = -1.0d0*pre_fac_E_test * E_test

    a_2 = a_21 + a_23 + a_21a + a_21b + a_22a + a_22b + a_23a + a_23b + a_2E

    ! Particle 3
    ! Acceleration from other particles
    a_31  = -1.0d0*pre_fac_c_test * (R_3 - R_1) / norm2(R_3 - R_1)**3
    a_32  = -1.0d0*pre_fac_c_test * (R_3 - R_2) / norm2(R_3 - R_2)**3

    ! Acceleration from image charge partners
    ! Other electron (+q*+q = +q**2)
    a_31a = +1.0*pre_fac_c_test * (R_3 - R_1a) / norm2(R_3 - R_1a)**3
    a_31b = +1.0*pre_fac_c_test * (R_3 - R_1b) / norm2(R_3 - R_1b)**3
    ! Other electron (+q*+q = +q**2) 
    a_32a = +1.0*pre_fac_c_test * (R_3 - R_2a) / norm2(R_3 - R_2a)**3
    a_32b = +1.0*pre_fac_c_test * (R_3 - R_2b) / norm2(R_3 - R_2b)**3
    ! Self (+q*-q = -q**2) 
    a_33a = -1.0*pre_fac_c_test * (R_3 - R_3a) / norm2(R_3 - R_3a)**3
    a_33b = -1.0*pre_fac_c_test * (R_3 - R_3b) / norm2(R_3 - R_3b)**3

    ! Acceleration from the electric field in the diode
    a_3E = +1.0d0*pre_fac_E_test * E_test

    a_3 = a_31 + a_32 + a_31a + a_31b + a_32a + a_32b + a_33a + a_33b + a_3E

    ! Set input variables
    call Setup_Test_System(d_test, delta_t_test, 100, V_test, .true., 0)

    ! Add particles
    par_vel = 0.0d0
    call Add_Particle(R_1, par_vel, species_elec, 1, 0, -1)
    call Add_Particle(R_2, par_vel, species_elec, 1, 0, -1)
    call Add_Particle(R_3, par_vel, species_ion, 1, 0, -1)

    call Calculate_Acceleration_Particles()

    ! Test the field at this location
    par_vel = (/ -4.55, -2.34, 96.44 /) * length_scale
    E_pos = Calc_Field_at(par_vel)

    ! Results from the acceleration
    a_1_res = particles_cur_accel(:, 1)
    a_2_res = particles_cur_accel(:, 2)
    a_3_res = particles_cur_accel(:, 3)

    if (all(abs(a_1 - a_1_res)/a_1 < tolerance_rel)) then
      print *, 'Particle 1 PASSED'
    else
      print *, 'Particle 1 FAILED'
      print *, a_1_res
      print *, a_1
      print *, ''
    end if

    if (all(abs(a_2 - a_2_res)/a_2 < tolerance_rel)) then
      print *, 'Particle 2 PASSED'
    else
      print *, 'Particle 2 FAILED'
      print *, 'Results'
      print *, a_2_res
      print *, 'Expected'
      print *, a_2
      print *, 'Difference'
      print *, abs(a_2 - a_2_res)
      print *, all(abs(a_2 - a_2_res) < tolerance_rel)
      print *, ''
    end if

    if (all(abs(a_3 - a_3_res)/a_3 < tolerance_rel)) then
      print *, 'Particle 3 PASSED'
    else
      print *, 'Particle 3 FAILED'
      print *, a_3_res
      print *, a_3
      print *, ''
    end if

    !E_python = (/ 6192429.94450208d0,  -4866704.16373579d0, -17454923.69763095d0 /) ! Results from Python script (no image charge)
    E_python = (/ -102526.05673208d0, 421022.84663293d0, -20456407.74487634d0 /)
    if (all(abs(E_python - E_pos)/E_python < tolerance_rel)) then
      print *, 'Electric field PASSED'
    else
      print *, 'Electric field FAILED'
      print *, 'E_python = ', E_python
      print *, 'E_pos = ', E_pos
      print *, 'abs(E_python - E_pos)/E_python = ', abs(E_python - E_pos)/E_python
      print *, ''
    end if

    print *, 'Acceleration test with image charge finished'
    print *, ''
  end subroutine Test_Acceleration_With_Image_Charge

  !-----------------------------------------------------------------------------
  ! We can compute the transit time of a single electron in the system using
  ! Z = Z_0 + V_0*t + 1/2*a*t^2, where a = qV/md, V_0 = 0, Z_0 = 1.0 nm and
  ! Z = d. This gives t = sqrt(2*m*d*(d-z_0)/(q*V)). If we dived this number
  ! with the time step and round up. We should know how many time steps should
  ! pass before the electron is absorbed. This test is to see if the
  ! Verlet algorithm is working properly.
  subroutine Test_Transit_Time()
    double precision, parameter      :: d_test = 1000.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d3 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step

    double precision                 :: time_exp ! Expected time
    integer                          :: steps_exp ! Exptected number of times steps

    integer                          :: i, steps_res

    logical                          :: found_res = .false. ! Set to true if
                                                            ! the electron has
                                                            ! been absorbed

    double precision, dimension(1:3) :: R_1, par_vel

    print *, 'Starting transit time test'
    R_1 = (/ 0.0d0, 0.0d0, 1.0d0*length_scale /)

    time_exp = sqrt(2.0d0*d_test*(d_test - R_1(3))*m_0/(q_0*V_test))
    steps_exp = ceiling(time_exp / delta_t_test)

    ! Set input variables
    call Setup_Test_System(d_test, delta_t_test, steps_exp+1000, V_test, .true., 0)

    ! Add particle
    par_vel = 0.0d0
    call Add_Particle(R_1, par_vel, species_elec, 1, 1, -1)

    do i = 1, steps

      ! Do Emission
      !call emission code

      ! Update the position of all particles
      !print *, 'Update position'
      call Update_Position(i)

      ! Remove particles from the system
      call Remove_Particles(i)

      !print *, 'i = ', i
      !print *, 'nrPart = ', nrPart
      if ((nrPart == 0) .and. (found_res .eqv. .false.)) then
        steps_res = i
        found_res = .true.
      end if

    end do

      do i = 1, MAX_LIFE_TIME
        if (life_time(i, species_elec) /= 0) then
          print *, i
          print *, life_time(i, species_elec)
          print *, ''
        end if
      end do

      print *, 'Transit time = ', steps_res
      print *, 'Exptected time = ', steps_exp
      ! The transit time should be within 1% of the expected number of steps.
      call Assert_True((abs(steps_exp - steps_res) < 0.01d0*steps_exp), 'Transit time test')

      print *, 'Transit time test finished'
      print *, ''
  end subroutine Test_Transit_Time

  !-----------------------------------------------------------------------------
  ! Test the acceleration and the batched field evaluation for the
  ! hyperboloid tip geometry (field_E_Hyperboloid + Sphere_IC_field).
  ! The expected values are computed here on the host with the same pair
  ! loop as the generic OpenMP version of Calculate_Acceleration_Particles,
  ! using direct calls to the geometry functions. In OpenACC builds
  ! Calculate_Acceleration_Particles and Calc_Field_at_Batch dispatch to the
  ! GPU kernels for this geometry, so this test verifies the device code
  ! against the host functions.
  subroutine Test_Tip_Acceleration()
    use mod_emission_tip, only: Init_Emission_Tip

    double precision, dimension(1:3)      :: par_vel
    double precision, dimension(1:3, 1:3) :: R, a_expected
    double precision, dimension(1:3)      :: pos_1, pos_2, diff, force_c, force_ic, force_ic_N
    double precision                      :: r_dis, inv_r3, pre_fac_c
    double precision, dimension(1:3, 1:4) :: field_pos, field_batch
    double precision, dimension(1:3)      :: field_single
    integer                               :: i, j, k
    character(len=64)                     :: label

    print '(a)', 'Starting hyperboloid tip test'

    ! d = 1000 nm gap with a 100 nm high tip: d_tip = 900 nm to the anode,
    ! base radius 100 nm.
    call Setup_Test_System(1000.0d0*length_scale, 0.25d-3*time_scale, 100, 2.0d3, .true., 0)
    emitters_dim(1:3, 1) = (/ 900.0d0, 100.0d0, 100.0d0 /) * length_scale
    emitters_pos(1:3, 1) = 0.0d0

    ! Sets the tip geometry pointers and computes the tip parameters
    ! (a_foci, eta_1, r_tip, shift_z, pre_fac_E_tip, ...).
    call Init_Emission_Tip()

    ! Three particles just above the tip apex (apex at z = h_tip = 100 nm).
    ! They are placed close to the tip and to each other on purpose: there
    ! the sphere image charge force is a large fraction of the total force,
    ! so the test is sensitive to errors in the image charge parameters
    ! (sphere center and radius) used by the OpenACC kernels.
    R(:, 1) = (/   2.0d0,   1.0d0, 103.0d0 /) * length_scale
    R(:, 2) = (/  -2.0d0,   2.5d0, 106.0d0 /) * length_scale
    R(:, 3) = (/   1.5d0,  -3.0d0, 110.0d0 /) * length_scale

    par_vel = 0.0d0
    call Add_Particle(R(:, 1), par_vel, species_elec, 1, 0, -1)
    call Add_Particle(R(:, 2), par_vel, species_elec, 1, 0, -1)
    call Add_Particle(R(:, 3), par_vel, species_ion, 1, 0, -1)

    ! Expected accelerations: the same pair loop as the OpenMP version of
    ! Calculate_Acceleration_Particles, using the geometry pointers directly.
    a_expected = 0.0d0
    do i = 1, 3
      pos_1 = particles_cur_pos(:, i)

      do j = i+1, 3
        pos_2 = particles_cur_pos(:, j)
        pre_fac_c = particles_charge(i) * particles_charge(j) * div_fac_c

        diff = pos_1 - pos_2
        r_dis = sqrt( sum(diff**2) ) + length_scale**2
        inv_r3 = 1.0d0 / (r_dis*r_dis*r_dis)
        force_c = (pre_fac_c*inv_r3) * diff

        force_ic = pre_fac_c * ptr_Image_Charge_effect(pos_1, pos_2)
        force_ic_N(1:2) = -1.0d0*force_ic(1:2)
        force_ic_N(3)   = force_ic(3)

        a_expected(:, j) = a_expected(:, j) + (force_ic_N - force_c) / particles_mass(j)
        a_expected(:, i) = a_expected(:, i) + (force_c + force_ic) / particles_mass(i)
      end do

      a_expected(:, i) = a_expected(:, i) &
                     & + particles_charge(i) * ptr_field_E(pos_1) / particles_mass(i)
    end do

    ! The acceleration from the kernel (GPU in OpenACC builds)
    call Calculate_Acceleration_Particles()

    do i = 1, 3
      write(label, '(a, i0)') 'Tip acceleration particle ', i
      call Assert_Close_Vec(particles_cur_accel(:, i), a_expected(:, i), trim(label))
    end do

    ! Batched field evaluation against the point-by-point Calc_Field_at
    field_pos(:, 1) = (/   0.0d0,   0.0d0, 102.0d0 /) * length_scale
    field_pos(:, 2) = (/  30.0d0, -10.0d0, 130.0d0 /) * length_scale
    field_pos(:, 3) = (/ -50.0d0,  40.0d0, 300.0d0 /) * length_scale
    field_pos(:, 4) = (/  15.0d0,   8.0d0, 500.0d0 /) * length_scale

    call Calc_Field_at_Batch(4, field_pos, field_batch)

    do k = 1, 4
      field_single = Calc_Field_at(field_pos(:, k))
      write(label, '(a, i0)') 'Tip batched field point ', k
      call Assert_Close_Vec(field_batch(:, k), field_single, trim(label))
    end do

    print '(a)', 'Hyperboloid tip test finished'
  end subroutine Test_Tip_Acceleration


  !-----------------------------------------------------------------------------
  ! Test the particle bookkeeping in mod_pair.
  ! Add a known mix of electrons and ions to an empty system and check that the
  ! particle counters (nrElec, nrIon, nrElecIon, nrPart) are updated correctly.
  subroutine Test_Particle_Bookkeeping()
    double precision, parameter      :: d_test = 100.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d0 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step
    double precision, dimension(1:3) :: R, par_vel

    print *, 'Starting particle bookkeeping test'

    ! Start with an empty system
    call Setup_Test_System(d_test, delta_t_test, 100, V_test, .false., 0)
    call Assert_True((nrPart == 0), 'Bookkeeping: empty system has no particles')

    par_vel = 0.0d0

    ! Add two electrons
    R = (/  3.0d0, -10.0d0,  2.0d0 /) * length_scale
    call Add_Particle(R, par_vel, species_elec, 1, 0, -1)
    R = (/ -9.0d0,  26.0d0, 80.0d0 /) * length_scale
    call Add_Particle(R, par_vel, species_elec, 1, 0, -1)

    ! Add one ion
    R = (/  6.0d0, -24.0d0, 56.53d0 /) * length_scale
    call Add_Particle(R, par_vel, species_ion, 1, 0, -1)

    call Assert_True((nrElec == 2),    'Bookkeeping: electron count')
    call Assert_True((nrIon == 1),     'Bookkeeping: ion count')
    call Assert_True((nrElecIon == 3), 'Bookkeeping: electron + ion count')
    call Assert_True((nrPart == 3),    'Bookkeeping: total particle count')

    print *, 'Particle bookkeeping test finished'
    print *, ''
  end subroutine Test_Particle_Bookkeeping

  !-----------------------------------------------------------------------------
  ! Test the planar electric field with an empty system.
  ! For a planar diode the field is uniform and given by E = -V/d in the
  ! z-direction with no x or y component.
  subroutine Test_Planar_Field()
    double precision, parameter      :: d_test = 100.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d0 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step
    double precision, dimension(1:3) :: pos, E_res, E_expected

    print *, 'Starting planar field test'

    ! Empty system, no image charge
    call Setup_Test_System(d_test, delta_t_test, 100, V_test, .false., 0)

    ! Evaluate the field at an arbitrary point inside the diode
    pos = (/ -4.55d0, -2.34d0, 96.44d0 /) * length_scale
    E_res = Calc_Field_at(pos)

    E_expected = (/ 0.0d0, 0.0d0, -1.0d0*V_test/d_test /)
    call Assert_Close_Vec(E_res, E_expected, 'Planar field with empty system')

    print *, 'Planar field test finished'
    print *, ''
  end subroutine Test_Planar_Field

  !-----------------------------------------------------------------------------
  ! A test with 1000 particles
  !
  subroutine Test_Many_Particles()
    double precision, parameter      :: d_test = 100.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d0 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step

    integer, parameter               :: n_par = 1000 ! Number of particles to read from the input file
    double precision, dimension(1:3) :: par_pos, par_vel

    integer                          :: ud_data, i, IFAIL

    print *, 'Starting many particles test'

    ! Set up the test system for 1000 time steps
    call Setup_Test_System(d_test, delta_t_test, 1000, V_test, .true., 0)

    ! Open the file with initial positions of particles
    open(newunit=ud_data, iostat=IFAIL, file='rand_pos_init.dt', status='OLD', action='read')
    if (IFAIL /= 0) then
      print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file input'
      return
    end if

    do i = 1, n_par
      ! Set the intial position and velocity
      read(unit=ud_data, fmt=*) par_pos ! Read the next line
      par_vel = 0.0d0

      ! Add the particle to the system
      call Add_Particle(par_pos, par_vel, species_elec, 1, 0, -1)
    end do

    ! Start the simulation
    do i = 1, steps

      ! Do Emission
      !call emission code

      ! Update the position of all particles
      !print *, 'Update position'
      call Update_Position(i)

      ! Remove particles from the system
      call Remove_Particles(i)
    end do

    print *, 'Many particles test finished'
  end subroutine Test_Many_Particles
end module mod_tests
