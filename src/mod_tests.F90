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
  use mod_work_function
  use mod_collisions, only: normal_dist, folded_normal_dist, folded_normal_max, &
                          & Calculate_Kramers_Cross_Section
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
    call Test_FN_Functions()

    print *, ''
    call Test_Collision_Math()

    print *, ''
    call Test_Work_Function_Checkerboard()

    print *, ''
    call Test_Acceleration_Without_Image_Charge()

    print *, ''
    call Test_Image_Charge()

    !call Test_Acceleration_With_Image_Charge()

    print *, ''
    call Test_Particle_Bookkeeping()

    print *, ''
    call Test_Particle_Removal()

    print *, ''
    call Test_Remove_All_Particles()

    print *, ''
    call Test_Pointer_Arrays()

    print *, ''
    call Test_Planar_Field()

    print *, ''
    call Test_Planar_Batch_Field()

    print *, ''
    call Test_Beeman_Kinematics()

    print *, ''
    call Test_TTS_Equivalence()

    print *, ''
    call Test_Transit_Time()

    !call Test_Many_Particles()

    print *, ''
    call Test_Tip_Acceleration()

    ! The two whole-system tests run last: they replace the emission,
    ! boundary and field pointers with the production ones for their
    ! geometry (Init_Field_Emission_v2 / Init_Emission_Tip).
    print *, ''
    call Test_Planar_System()

    print *, ''
    call Test_Tip_System()

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
    nrAtom     = 0
    nrElecIon  = 0

    nrPart_remove = 0
    nrElec_remove = 0
    nrIon_remove  = 0
    nrAtom_remove = 0

    nrPart_remove_top = 0
    nrPart_remove_bot = 0
    nrPart_remove_recom = 0
    nrPart_remove_ion = 0
    nrElec_remove_top = 0
    nrElec_remove_bot = 0
    nrElec_remove_recom = 0
    nrIon_remove_top = 0
    nrIon_remove_bot = 0
    nrIon_remove_recom = 0
    nrAtom_remove_ion = 0

    nrElec_remove_top_emit = 0

    ! Tests that exercise the two-time-step path set this themselves and are
    ! expected to restore it, but reset it here so a failing test cannot leak
    ! the mode into the tests that follow it.
    two_time_step = .false.

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

    call Update_Particle_Acceleration(1)
    ! call Calculate_Acceleration_Particles()

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

    call Update_Particle_Acceleration(1)
    ! call Calculate_Acceleration_Particles()

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
  ! Constant work function used by the tests below. Matches the Work_fun
  ! abstract interface in mod_work_function.
  double precision function Const_Work_Fun_Test(pos, emit, sec)
    double precision, dimension(1:3), intent(in) :: pos
    integer, intent(in)                          :: emit
    integer, intent(out), optional               :: sec

    if (present(sec) .eqv. .true.) sec = 1
    Const_Work_Fun_Test = 4.7d0
  end function Const_Work_Fun_Test

  !-----------------------------------------------------------------------------
  ! Count the number of particles whose species pointer arrays are
  ! inconsistent. For every particle i the inverse pointer k must map back to
  ! i through the pointer array of its own species (see Check_Pointer_Arrays
  ! in mod_pair). Returns 0 when the bookkeeping is consistent.
  integer function Count_Pointer_Mismatches()
    integer :: i, k, j

    Count_Pointer_Mismatches = 0
    do i = 1, nrPart
      k = particles_inverse_pointer(i)
      j = -1
      select case (particles_species(i))
        case (species_elec)
          if ((k >= 1) .and. (k <= nrElec)) j = particles_elec_pointer(k)
        case (species_ion)
          if ((k >= 1) .and. (k <= nrIon)) j = particles_ion_pointer(k)
        case (species_atom)
          if ((k >= 1) .and. (k <= nrAtom)) j = particles_atom_pointer(k)
      end select
      if (j /= i) Count_Pointer_Mismatches = Count_Pointer_Mismatches + 1
    end do
  end function Count_Pointer_Mismatches

  !-----------------------------------------------------------------------------
  ! Test Mark_Particles_Remove and Remove_Particles in mod_pair.
  ! A known mix of electrons and ions is added, three particles are marked for
  ! removal (top, bottom and top) and removed. The test checks the removal
  ! counters, that the compaction (compact_array_*) keeps the data of the
  ! survivors intact and in order, that the life time histogram is recorded,
  ! and that the masks are reset afterwards.
  subroutine Test_Particle_Removal()
    double precision, parameter      :: d_test = 100.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d0 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step

    integer, parameter                    :: n_add = 7
    double precision, dimension(1:3, 1:n_add) :: R, R_prev, Vel
    integer, dimension(1:n_add)           :: spc
    integer, dimension(1:4)               :: keep, expect_spc, expect_id
    character(len=64)                     :: label
    integer                               :: i, k

    print *, 'Starting particle removal test'

    call Setup_Test_System(d_test, delta_t_test, 100, V_test, .false., 0)

    ! Electrons everywhere except ions at 3 and 6
    spc = species_elec
    spc(3) = species_ion
    spc(6) = species_ion

    do i = 1, n_add
      ! Distinct data for every particle so the compaction can be checked
      R(:, i)   = (/ 1.0d0*i, -2.0d0*i, 10.0d0*i /) * length_scale
      Vel(:, i) = (/ 100.0d0*i, -50.0d0*i, 25.0d0*i /)
      call Add_Particle(R(:, i), Vel(:, i), spc(i), 1, 1, -1)

      ! Overwrite the previous position with a distinct value as well, so a
      ! second 2D array with different content goes through the compaction.
      R_prev(:, i) = R(:, i) + 0.5d0*length_scale
      particles_prev_pos(:, i) = R_prev(:, i)
    end do

    call Assert_True((nrPart == 7) .and. (nrElec == 5) .and. (nrIon == 2), &
                   & 'Removal: system before removal')

    ! Mark an electron at the top, an electron at the bottom and an ion at the
    ! top. Particle 2 is marked twice: the second call must be a no-op.
    call Mark_Particles_Remove(2, remove_top)
    call Mark_Particles_Remove(2, remove_top)
    call Mark_Particles_Remove(4, remove_bot)
    call Mark_Particles_Remove(6, remove_top)

    call Assert_True((nrPart_remove == 3), 'Removal: marked count (double mark is a no-op)')
    call Assert_True((nrElec_remove == 2) .and. (nrIon_remove == 1), &
                   & 'Removal: marked count per species')
    call Assert_True((nrElec_remove_top == 1) .and. (nrElec_remove_bot == 1) &
                   & .and. (nrIon_remove_top == 1), 'Removal: marked count per boundary')
    call Assert_True((particles_charge(2) == 0.0d0), 'Removal: marked particle is neutralized')

    ! Particles were added at step 1 and are removed at step 11:
    ! life time = 10 steps.
    call Remove_Particles(11)

    call Assert_True((nrPart == 4), 'Removal: total count after removal')
    call Assert_True((nrElec == 3) .and. (nrIon == 1) .and. (nrElecIon == 4), &
                   & 'Removal: species counts after removal')

    ! The survivors are the old particles 1, 3, 5 and 7 in that order
    keep       = (/ 1, 3, 5, 7 /)
    expect_spc = (/ species_elec, species_ion, species_elec, species_elec /)
    expect_id  = (/ 0, 2, 4, 6 /) ! IDs are assigned 0, 1, 2, ... in Add_Particle

    do k = 1, 4
      i = keep(k)
      write(label, '(a, i0)') 'Removal: position of survivor ', k
      call Assert_Close_Vec(particles_cur_pos(:, k)/length_scale, R(:, i)/length_scale, trim(label))
      write(label, '(a, i0)') 'Removal: prev position of survivor ', k
      call Assert_Close_Vec(particles_prev_pos(:, k)/length_scale, R_prev(:, i)/length_scale, trim(label))
      write(label, '(a, i0)') 'Removal: velocity of survivor ', k
      call Assert_Close_Vec(particles_cur_vel(:, k), Vel(:, i), trim(label))
    end do

    call Assert_True(all(particles_species(1:4) == expect_spc), 'Removal: species of survivors')
    call Assert_True(all(particles_id(1:4) == expect_id), 'Removal: IDs of survivors')
    call Assert_True(all(particles_emitter(1:4) == 1) .and. all(particles_section(1:4) == 1), &
                   & 'Removal: emitter and section of survivors')

    ! Charge and mass of the survivors (scaled to O(1) numbers so that the
    ! relative tolerance in Assert_Close applies)
    call Assert_Close(particles_charge(1)/q_0, -1.0d0, 'Removal: charge of surviving electron')
    call Assert_Close(particles_charge(2)/q_0, +1.0d0, 'Removal: charge of surviving ion')
    call Assert_Close(particles_mass(1)/m_0, m_eeff, 'Removal: mass of surviving electron')
    call Assert_Close(particles_mass(2)/m_N2p, 1.0d0, 'Removal: mass of surviving ion')

    ! Life times: two electrons and one ion lived 10 steps
    call Assert_True((life_time(10, species_elec) == 2), 'Removal: electron life time recorded')
    call Assert_True((life_time(10, species_ion) == 1), 'Removal: ion life time recorded')

    ! The masks over the old live range must be reset and the counters cleared
    call Assert_True(all(particles_mask(1:7)), 'Removal: masks reset after removal')
    call Assert_True((nrPart_remove == 0) .and. (nrElec_remove == 0) .and. (nrIon_remove == 0), &
                   & 'Removal: remove counters cleared')

    print *, 'Particle removal test finished'
    print *, ''
  end subroutine Test_Particle_Removal

  !-----------------------------------------------------------------------------
  ! Remove every particle in the system. This takes the branch in
  ! Remove_Particles that skips the compaction entirely, which the test above
  ! does not reach. Afterwards the system must be empty and usable again.
  subroutine Test_Remove_All_Particles()
    double precision, parameter      :: d_test = 100.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d0 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step
    double precision, dimension(1:3) :: R, par_vel

    print *, 'Starting remove-all-particles test'

    call Setup_Test_System(d_test, delta_t_test, 100, V_test, .false., 0)

    par_vel = 0.0d0
    R = (/  3.0d0, -10.0d0,  2.0d0 /) * length_scale
    call Add_Particle(R, par_vel, species_elec, 1, 1, -1)
    R = (/ -9.0d0,  26.0d0, 80.0d0 /) * length_scale
    call Add_Particle(R, par_vel, species_elec, 1, 1, -1)

    call Mark_Particles_Remove(1, remove_top)
    call Mark_Particles_Remove(2, remove_top)
    call Remove_Particles(21)

    call Assert_True((nrPart == 0) .and. (nrElec == 0), 'Remove all: system is empty')
    call Assert_True((nrPart_remove == 0), 'Remove all: remove counters cleared')
    call Assert_True(all(particles_mask(1:2)), 'Remove all: masks reset')

    ! The system must accept new particles again
    R = (/ 5.0d0, 6.0d0, 50.0d0 /) * length_scale
    call Add_Particle(R, par_vel, species_elec, 2, 1, -1)
    call Assert_True((nrPart == 1) .and. (nrElec == 1), 'Remove all: particle added to empty system')
    call Assert_Close_Vec(particles_cur_pos(:, 1)/length_scale, R/length_scale, &
                        & 'Remove all: position of new particle')

    print *, 'Remove-all-particles test finished'
    print *, ''
  end subroutine Test_Remove_All_Particles

  !-----------------------------------------------------------------------------
  ! Test the species pointer bookkeeping through a removal in two-time-step
  ! mode. In that mode Remove_Particles rebuilds the pointer arrays
  ! (Update_Pointer_Arrays_parallel) before compacting them, and the
  ! two-time-step integrator and the POLARSO solver rely on the result. The
  ! test adds electrons, ions and an atom, removes one of each species and
  ! checks that every pointer round trip (particle -> inverse pointer ->
  ! species pointer -> particle) is consistent afterwards.
  subroutine Test_Pointer_Arrays()
    double precision, parameter      :: d_test = 100.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d0 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step

    integer, parameter               :: n_add = 7
    integer, dimension(1:n_add)      :: spc
    integer, dimension(1:4)          :: expect_spc, expect_id
    double precision, dimension(1:3) :: R, par_vel
    integer                          :: i

    print *, 'Starting pointer array test'

    call Setup_Test_System(d_test, delta_t_test, 100, V_test, .false., 0)
    two_time_step = .true. ! Take the pointer rebuild path in Remove_Particles

    spc = (/ species_elec, species_ion, species_elec, species_atom, &
           & species_elec, species_ion, species_atom /)

    par_vel = 0.0d0
    do i = 1, n_add
      R = (/ 2.0d0*i, -1.0d0*i, 5.0d0*i /) * length_scale
      call Add_Particle(R, par_vel, spc(i), 1, 1, -1)
    end do

    call Assert_True((nrElec == 3) .and. (nrIon == 2) .and. (nrAtom == 2), &
                   & 'Pointers: counts after adding')
    call Assert_True((Count_Pointer_Mismatches() == 0), 'Pointers: consistent after adding')

    ! Remove one particle of each species
    call Mark_Particles_Remove(3, remove_top) ! electron
    call Mark_Particles_Remove(2, remove_bot) ! ion
    call Mark_Particles_Remove(7, remove_ion) ! atom
    call Remove_Particles(31)

    call Assert_True((nrPart == 4) .and. (nrElec == 2) .and. (nrIon == 1) .and. (nrAtom == 1), &
                   & 'Pointers: counts after removal')

    ! Survivors are the old particles 1 (e), 4 (a), 5 (e) and 6 (i)
    expect_spc = (/ species_elec, species_atom, species_elec, species_ion /)
    expect_id  = (/ 0, 3, 4, 5 /)
    call Assert_True(all(particles_species(1:4) == expect_spc), 'Pointers: species of survivors')
    call Assert_True(all(particles_id(1:4) == expect_id), 'Pointers: IDs of survivors')

    call Assert_True((Count_Pointer_Mismatches() == 0), 'Pointers: consistent after removal')

    ! With no marked particles the pointer rebuild must be a no-op
    call Update_Pointer_Arrays_sequential()
    call Assert_True((Count_Pointer_Mismatches() == 0), 'Pointers: sequential rebuild is a no-op')

    two_time_step = .false.

    print *, 'Pointer array test finished'
    print *, ''
  end subroutine Test_Pointer_Arrays

  !-----------------------------------------------------------------------------
  ! Test the Beeman integrator and the Ramo current with a single electron in
  ! a uniform field, where every step can be written down exactly.
  ! The electron starts at rest in z with a transverse drift velocity v_x.
  ! With a constant acceleration a = q V / (m d) and zeroed acceleration
  ! history the Beeman scheme gives after three steps:
  !   z3 - z0 = 3 a dt^2,   v_z3 = 5/2 a dt,   x3 - x0 = 3 v_x dt
  ! and the Ramo current of the last step is I = q_0 v_z3 / d.
  subroutine Test_Beeman_Kinematics()
    double precision, parameter      :: d_test = 1000.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d3 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step
    double precision, parameter      :: v_x0 = 1.0d3 ! Transverse drift velocity (m/s)

    double precision, dimension(1:3) :: R_0, par_vel
    double precision                 :: a_z, v_z_exp, dz_exp, dx_exp, ramo_exp
    integer                          :: i

    print *, 'Starting Beeman kinematics test'

    call Setup_Test_System(d_test, delta_t_test, 10, V_test, .false., 0)

    ! Make sure the planar geometry is active
    ptr_field_E => field_E_planar
    ptr_Image_Charge_effect => Force_Image_charges_v2
    ptr_Check_Boundary => Check_Boundary_Planar
    ptr_E_zunit => E_zunit_planar

    R_0 = (/ 0.0d0, 0.0d0, 500.0d0*length_scale /)
    par_vel = (/ v_x0, 0.0d0, 0.0d0 /)
    call Add_Particle(R_0, par_vel, species_elec, 1, 1, -1)

    do i = 1, 3
      call Update_Position(i)
    end do

    ! Expected values, see the header of this subroutine
    a_z      = q_0*V_test/(m_0*d_test)
    dz_exp   = 3.0d0*a_z*delta_t_test**2
    v_z_exp  = 2.5d0*a_z*delta_t_test
    dx_exp   = 3.0d0*v_x0*delta_t_test
    ramo_exp = q_0*v_z_exp/d_test

    call Assert_Close((particles_cur_pos(3, 1) - R_0(3))/length_scale, dz_exp/length_scale, &
                    & 'Beeman: z-displacement after 3 steps')
    call Assert_Close(particles_cur_vel(3, 1), v_z_exp, 'Beeman: z-velocity after 3 steps')
    call Assert_Close((particles_cur_pos(1, 1) - R_0(1))/length_scale, dx_exp/length_scale, &
                    & 'Beeman: transverse drift after 3 steps')
    call Assert_Close(particles_cur_vel(1, 1), v_x0, 'Beeman: transverse velocity is unchanged')
    call Assert_Close((particles_cur_pos(2, 1) - R_0(2))/length_scale, 0.0d0, &
                    & 'Beeman: no drift in y')
    ! Scaled to nA so that the relative tolerance in Assert_Close applies
    call Assert_Close(ramo_current(species_elec)*1.0d9, ramo_exp*1.0d9, &
                    & 'Beeman: Ramo current of a single electron')

    print *, 'Beeman kinematics test finished'
    print *, ''
  end subroutine Test_Beeman_Kinematics

  !-----------------------------------------------------------------------------
  ! The two-time-step integrator with atom_time_interval = 1 must reproduce
  ! the standard integrator: every species is then updated on every step with
  ! the same Beeman formulas, only through the species pointer lists instead
  ! of one loop over all particles. The same three-particle system is run for
  ! three steps in both modes and the trajectories, the Ramo current and
  ! Calc_Field_at (which dispatches to Calc_Field_at_tts) are compared.
  subroutine Test_TTS_Equivalence()
    double precision, parameter      :: d_test = 1000.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d3 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step

    double precision, dimension(1:3, 1:3) :: R_0, V_0, pos_ref, vel_ref
    double precision, dimension(1:3)      :: probe, E_ots, E_tts
    double precision, dimension(1:2)      :: ramo_ref
    integer, dimension(1:3)               :: spc
    character(len=64)                     :: label
    integer                               :: i, run

    print *, 'Starting two-time-step equivalence test'

    ! Two electrons and one ion a few nm apart: the pair forces are a large
    ! part of the acceleration, so the test is sensitive to the pair loops of
    ! the two-time-step routines, not only to the vacuum field.
    R_0(:, 1) = (/  1.0d0, -2.0d0, 500.0d0 /) * length_scale
    R_0(:, 2) = (/  3.0d0,  4.0d0, 503.0d0 /) * length_scale
    R_0(:, 3) = (/ -2.0d0,  1.0d0, 505.0d0 /) * length_scale
    V_0(:, 1) = (/  1.0d3,  0.0d0, 0.0d0 /)
    V_0(:, 2) = (/  0.0d0, -1.0d3, 0.0d0 /)
    V_0(:, 3) = 0.0d0
    spc = (/ species_elec, species_elec, species_ion /)

    probe = (/ 0.0d0, 1.0d0, 502.0d0 /) * length_scale

    do run = 1, 2
      call Setup_Test_System(d_test, delta_t_test, 10, V_test, .true., 0)

      ! Make sure the planar geometry is active
      ptr_field_E => field_E_planar
      ptr_Image_Charge_effect => Force_Image_charges_v2
      ptr_Check_Boundary => Check_Boundary_Planar
      ptr_E_zunit => E_zunit_planar

      if (run == 2) then
        two_time_step = .true.
        atom_time_interval = 1 ! Update ions on every step
      end if

      do i = 1, 3
        call Add_Particle(R_0(:, i), V_0(:, i), spc(i), 1, 1, -1)
      end do

      if (run == 1) then
        ! The field summed over the species pointer lists must equal the
        ! field summed over all particles.
        E_ots = Calc_Field_at(probe)
        two_time_step = .true.
        E_tts = Calc_Field_at(probe)
        two_time_step = .false.
        call Assert_Close_Vec(E_tts, E_ots, 'TTS: field equals the standard field')
      end if

      do i = 1, 3
        call Update_Position(i)
      end do

      if (run == 1) then
        pos_ref = particles_cur_pos(:, 1:3)
        vel_ref = particles_cur_vel(:, 1:3)
        ramo_ref(1) = ramo_current(species_elec)
        ramo_ref(2) = ramo_current(species_ion)
      else
        do i = 1, 3
          write(label, '(a, i0)') 'TTS: trajectory of particle ', i
          call Assert_Close_Vec(particles_cur_pos(:, i)/length_scale, &
                              & pos_ref(:, i)/length_scale, trim(label))
          write(label, '(a, i0)') 'TTS: velocity of particle ', i
          call Assert_Close_Vec(particles_cur_vel(:, i), vel_ref(:, i), trim(label))
        end do
        ! Scaled to nA so that the relative tolerance in Assert_Close applies
        call Assert_Close(ramo_current(species_elec)*1.0d9, ramo_ref(1)*1.0d9, &
                        & 'TTS: electron Ramo current')
        call Assert_Close(ramo_current(species_ion)*1.0d9, ramo_ref(2)*1.0d9, &
                        & 'TTS: ion Ramo current')
      end if
    end do

    two_time_step = .false.

    print *, 'Two-time-step equivalence test finished'
    print *, ''
  end subroutine Test_TTS_Equivalence

  !-----------------------------------------------------------------------------
  ! Test the batched field evaluation for the planar geometry against the
  ! point-by-point Calc_Field_at, with image charge partners enabled
  ! (N_IC_MAX = 1 exercises the image charge series). In OpenACC builds
  ! Calc_Field_at_Batch dispatches to the planar GPU kernel, so this verifies
  ! the device code against the host functions; the hyperboloid tip geometry
  ! has the same check in Test_Tip_Acceleration.
  subroutine Test_Planar_Batch_Field()
    double precision, parameter      :: d_test = 100.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d0 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step

    double precision, dimension(1:3, 1:4) :: field_pos, field_batch
    double precision, dimension(1:3)      :: R, par_vel, field_single, E_vac
    character(len=64)                     :: label
    integer                               :: k

    print *, 'Starting planar batched field test'

    call Setup_Test_System(d_test, delta_t_test, 100, V_test, .true., 1)

    ! Make sure the planar geometry is active
    ptr_field_E => field_E_planar
    ptr_Image_Charge_effect => Force_Image_charges_v2

    field_pos(:, 1) = (/   0.0d0,   0.0d0, 10.0d0 /) * length_scale
    field_pos(:, 2) = (/  12.0d0,  -4.0d0, 21.0d0 /) * length_scale
    field_pos(:, 3) = (/ -20.0d0,  30.0d0, 80.0d0 /) * length_scale
    field_pos(:, 4) = (/  40.0d0,  40.0d0, 95.0d0 /) * length_scale

    ! Empty system: the batch must return the vacuum field at every point
    E_vac = (/ 0.0d0, 0.0d0, -1.0d0*V_test/d_test /)
    call Calc_Field_at_Batch(4, field_pos, field_batch)
    do k = 1, 4
      write(label, '(a, i0)') 'Planar batched vacuum field point ', k
      call Assert_Close_Vec(field_batch(:, k), E_vac, trim(label))
    end do

    ! Add particles, point 2 is close to the first electron on purpose
    par_vel = 0.0d0
    R = (/  10.0d0, -5.0d0, 20.0d0 /) * length_scale
    call Add_Particle(R, par_vel, species_elec, 1, 1, -1)
    R = (/ -15.0d0,  8.0d0, 60.0d0 /) * length_scale
    call Add_Particle(R, par_vel, species_elec, 1, 1, -1)
    R = (/   5.0d0, 25.0d0, 40.0d0 /) * length_scale
    call Add_Particle(R, par_vel, species_ion, 1, 1, -1)

    call Calc_Field_at_Batch(4, field_pos, field_batch)
    do k = 1, 4
      field_single = Calc_Field_at(field_pos(:, k))
      write(label, '(a, i0)') 'Planar batched field point ', k
      call Assert_Close_Vec(field_batch(:, k), field_single, trim(label))
    end do

    print *, 'Planar batched field test finished'
    print *, ''
  end subroutine Test_Planar_Batch_Field

  !-----------------------------------------------------------------------------
  ! Test the Fowler-Nordheim functions v_y, t_y and Escape_Prob in
  ! mod_field_emission_v2 against reference values computed independently
  ! (Python) with the same physical constants, for a constant work function
  ! of 4.7 eV. Also checks the identities v_y = t_y = 1 without image charge
  ! and the clamping of l = y^2 to 1 for very strong fields.
  subroutine Test_FN_Functions()
    double precision, parameter      :: d_test = 1000.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d3 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step

    double precision, dimension(1:3) :: pos
    double precision                 :: F

    print *, 'Starting Fowler-Nordheim function test'

    call Setup_Test_System(d_test, delta_t_test, 100, V_test, .false., 0)

    ! Constant work function of 4.7 eV
    ptr_Work_fun => Const_Work_Fun_Test

    pos = 0.0d0

    ! Without the image charge effect both correction factors are exactly 1
    F = -2.0d9 ! V/m
    call Assert_Close(v_y(F, pos, 1), 1.0d0, 'FN: v_y = 1 without image charge')
    call Assert_Close(t_y(F, pos, 1), 1.0d0, 'FN: t_y = 1 without image charge')

    ! With the image charge effect, compare against the reference values
    image_charge = .true.

    F = -2.0d9 ! V/m, l = 0.13037...
    call Assert_Close(v_y(F, pos, 1), 0.8253581935658024d0, 'FN: v_y at F = -2 GV/m')
    call Assert_Close(t_y(F, pos, 1), 1.0292422630703880d0, 'FN: t_y at F = -2 GV/m')
    ! The escape probability itself (3.354e-13) is below the absolute
    ! tolerance of Assert_Close, so compare its logarithm instead.
    call Assert_Close(log(Escape_Prob(F, pos, 1)), -28.723444978828507d0, &
                    & 'FN: log escape probability at F = -2 GV/m')

    F = -4.0d9 ! V/m, l = 0.26074...
    call Assert_Close(v_y(F, pos, 1), 0.6808388366998485d0, 'FN: v_y at F = -4 GV/m')
    call Assert_Close(t_y(F, pos, 1), 1.0484437096180323d0, 'FN: t_y at F = -4 GV/m')
    call Assert_Close(log(Escape_Prob(F, pos, 1)), -11.846999895226958d0, &
                    & 'FN: log escape probability at F = -4 GV/m')

    ! For fields so strong that l = y^2 > 1 (barrier fully suppressed) v_y
    ! clamps l to 1, and v_y = 1 - 1 + 1/6 * 1 * log(1) = 0
    F = -1.0d12 ! V/m
    call Assert_Close(v_y(F, pos, 1), 0.0d0, 'FN: v_y clamps to 0 above the barrier maximum')

    image_charge = .false.

    print *, 'Fowler-Nordheim function test finished'
    print *, ''
  end subroutine Test_FN_Functions

  !-----------------------------------------------------------------------------
  ! Test the statistical helper functions and the Kramers recombination cross
  ! section in mod_collisions against reference values computed independently
  ! (Python) with the same physical constants.
  subroutine Test_Collision_Math()
    double precision :: fmax, f, x
    logical          :: envelope_ok
    integer          :: k

    print *, 'Starting collision math test'

    ! Normal distribution
    call Assert_Close(normal_dist(0.0d0, 1.0d0, 0.0d0), 0.3989422804014327d0, &
                    & 'Collisions: standard normal at 0')
    call Assert_Close(normal_dist(5.0d0, 25.0d0, 12.5d0), 0.015255512618420964d0, &
                    & 'Collisions: normal distribution')

    ! Folded normal distribution: reference value and the identity
    ! f(x; mu, sigma) = N(x; mu, sigma) + N(x; -mu, sigma)
    call Assert_Close(folded_normal_dist(5.0d0, 25.0d0, 30.0d0), 0.015667927606195526d0, &
                    & 'Collisions: folded normal distribution')
    call Assert_Close(folded_normal_dist(12.0d0, 40.0d0, 77.0d0), &
                    & normal_dist(12.0d0, 40.0d0, 77.0d0) + normal_dist(-12.0d0, 40.0d0, 77.0d0), &
                    & 'Collisions: folded normal sum identity')

    ! folded_normal_max is used as the envelope constant in the
    ! acceptance-rejection sampling of the scattering angle: it must bound the
    ! distribution on the angular domain [0, 180].
    fmax = folded_normal_max(5.0d0, 25.0d0)
    envelope_ok = .true.
    do k = 0, 180
      x = 1.0d0*k
      f = folded_normal_dist(5.0d0, 25.0d0, x)
      if (f > fmax*(1.0d0 + 1.0d-12)) envelope_ok = .false.
    end do
    call Assert_True(envelope_ok, 'Collisions: folded normal max bounds the distribution')

    ! Kramers cross section. The values (~1e-26 m^2) are far below the
    ! absolute tolerance of Assert_Close, so compare them scaled to O(1).
    call Assert_Close(Calculate_Kramers_Cross_Section(10.0d0)*1.0d26, 3.995353707087291d0, &
                    & 'Collisions: Kramers cross section at 10 eV')
    call Assert_Close(Calculate_Kramers_Cross_Section(100.0d0)*1.0d28, 8.842728751351863d0, &
                    & 'Collisions: Kramers cross section at 100 eV')
    call Assert_True(Calculate_Kramers_Cross_Section(10.0d0) > Calculate_Kramers_Cross_Section(100.0d0), &
                   & 'Collisions: Kramers cross section decreases with energy')

    print *, 'Collision math test finished'
    print *, ''
  end subroutine Test_Collision_Math

  !-----------------------------------------------------------------------------
  ! Test the checkerboard work function model in mod_work_function with a
  ! 2x2 pattern on a 100x100 nm emitter: the value lookup (with the reversed
  ! y-direction in the matrix), the section numbering and the clamping of
  ! positions outside the emitter.
  subroutine Test_Work_Function_Checkerboard()
    double precision, dimension(1:3) :: pos
    integer                          :: sec

    print *, 'Starting checkerboard work function test'

    ! A single 100 x 100 nm rectangular emitter centered at the origin of the
    ! unit square mapping (emitters_pos is the lower left corner).
    nrEmit = 1
    emitters_Type(1) = EMIT_RECTANGLE
    emitters_pos(1:3, 1) = 0.0d0
    emitters_dim(1:3, 1) = (/ 100.0d0, 100.0d0, 0.0d0 /) * length_scale

    ! 2x2 checkerboard. Row 1 of the matrix is the TOP row of the pattern
    ! (the y-direction is reversed in the array).
    x_num = 2
    y_num = 2
    x_len = 1.0d0/x_num
    y_len = 1.0d0/y_num
    if (allocated(w_theta_arr)) deallocate(w_theta_arr)
    allocate(w_theta_arr(1:y_num, 1:x_num))
    w_theta_arr(1, 1:2) = (/ 4.10d0, 4.20d0 /) ! Top row
    w_theta_arr(2, 1:2) = (/ 4.30d0, 4.40d0 /) ! Bottom row

    ptr_Work_fun => w_theta_checkerboard

    ! Centers of the four sections. The section numbering starts at the
    ! bottom left and runs row-wise:  | 3 4 |
    !                                 | 1 2 |
    pos = (/ 25.0d0, 25.0d0, 0.0d0 /) * length_scale ! Bottom left
    call Assert_Close(w_theta_xy(pos, 1, sec), 4.30d0, 'Work function: bottom left value')
    call Assert_True((sec == 1), 'Work function: bottom left section')

    pos = (/ 75.0d0, 25.0d0, 0.0d0 /) * length_scale ! Bottom right
    call Assert_Close(w_theta_xy(pos, 1, sec), 4.40d0, 'Work function: bottom right value')
    call Assert_True((sec == 2), 'Work function: bottom right section')

    pos = (/ 25.0d0, 75.0d0, 0.0d0 /) * length_scale ! Top left
    call Assert_Close(w_theta_xy(pos, 1, sec), 4.10d0, 'Work function: top left value')
    call Assert_True((sec == 3), 'Work function: top left section')

    pos = (/ 75.0d0, 75.0d0, 0.0d0 /) * length_scale ! Top right
    call Assert_Close(w_theta_xy(pos, 1, sec), 4.20d0, 'Work function: top right value')
    call Assert_True((sec == 4), 'Work function: top right section')

    ! Positions outside the emitter clamp to the nearest section
    pos = (/ 150.0d0, 25.0d0, 0.0d0 /) * length_scale ! Right of the emitter
    call Assert_Close(w_theta_xy(pos, 1, sec), 4.40d0, 'Work function: clamped in x')
    call Assert_True((sec == 2), 'Work function: clamped section in x')

    pos = (/ -10.0d0, 130.0d0, 0.0d0 /) * length_scale ! Left of and above the emitter
    call Assert_Close(w_theta_xy(pos, 1, sec), 4.10d0, 'Work function: clamped in x and y')
    call Assert_True((sec == 3), 'Work function: clamped section in x and y')

    ! Leave a valid work function behind for any test that runs after this one
    deallocate(w_theta_arr)
    ptr_Work_fun => Const_Work_Fun_Test

    print *, 'Checkerboard work function test finished'
    print *, ''
  end subroutine Test_Work_Function_Checkerboard

  !-----------------------------------------------------------------------------
  ! Seed the random number generator with a fixed seed so that the system
  ! tests below are reproducible from run to run (they are still not bitwise
  ! identical across compilers, whose RANDOM_NUMBER generators differ; the
  ! assertions use invariants and generous bands for that reason).
  subroutine Seed_RNG_Deterministic()
    integer                            :: n, i
    integer, dimension(:), allocatable :: seed

    call random_seed(size=n)
    allocate(seed(1:n))
    do i = 1, n
      seed(i) = 12345 + 37*(i - 1)
    end do
    call random_seed(put=seed)
    deallocate(seed)
  end subroutine Seed_RNG_Deterministic

  !-----------------------------------------------------------------------------
  ! Whole system test for the planar diode: run the production field emission
  ! (Init_Field_Emission_v2 / Do_Field_Emission with the Cuba supply integral
  ! and Metropolis-Hastings sampling, using the work function read from the
  ! file 'work') together with the integrator and the absorption for a few
  ! hundred steps, exactly like the main loop does. The checks are physical
  ! invariants:
  !   * charge bookkeeping: emitted = absorbed + still in flight (exact)
  !   * particles never end up outside the gap and positions stay finite
  !   * Ramo charge accounting: the time-integrated Ramo current must equal
  !     q_0/d times the total z-displacement of all electrons (Shockley-Ramo
  !     theorem), which is known from the absorption counts and the positions
  !     of the electrons still in flight
  !   * the steady-state current lies in a band around the value this diode
  !     (1 kV over 500 nm, 100x100 nm emitter, 2.0 eV) is known to give
  subroutine Test_Planar_System()
    double precision, parameter :: d_test = 500.0d0*length_scale
    double precision, parameter :: V_test = 1.0d3 ! V
    double precision, parameter :: delta_t_test = 0.25d-15 ! Time step
    integer, parameter          :: n_steps = 250
    double precision, parameter :: z_emit = 1.0d0*length_scale ! Emission height used by Do_Field_Emission

    ! Band for the steady-state current (second half of the run). The value
    ! is stochastic (Metropolis-Hastings) and the RNG differs between
    ! compilers, so the band is generous.
    double precision, parameter :: I_low  = 0.5d-3 ! A
    double precision, parameter :: I_high = 4.0d-3 ! A

    integer          :: i, n_before, n_emit_tot, n_abs_top, n_abs_bot, n_peak
    double precision :: Q_ramo, Q_ramo_ss, Q_exp, I_avg_ss, sum_z
    logical          :: bounds_ok, finite_ok

    print *, 'Starting planar field emission system test'

    call Setup_Test_System(d_test, delta_t_test, n_steps, V_test, .true., 0)

    ! One 100 x 100 nm emitter centered in the box
    nrEmit = 1
    emitters_pos(1:3, 1) = (/ -50.0d0, -50.0d0, 0.0d0 /) * length_scale
    emitters_dim(1:3, 1) = (/ 100.0d0, 100.0d0, 0.0d0 /) * length_scale
    emitters_type(1)     = EMIT_RECTANGLE
    emitters_delay(1)    = 0

    ! The production init: sets the planar pointers, allocates the emission
    ! bookkeeping and reads the work function from the file 'work' (2.0 eV).
    call Init_Field_Emission_v2()

    call Seed_RNG_Deterministic()

    n_emit_tot = 0
    n_abs_top  = 0
    n_abs_bot  = 0
    n_peak     = 0
    Q_ramo     = 0.0d0
    Q_ramo_ss  = 0.0d0
    bounds_ok  = .true.
    finite_ok  = .true.

    ! The same sequence as the main loop
    do i = 1, n_steps
      n_before = nrElec
      call ptr_Do_Emission(i)
      n_emit_tot = n_emit_tot + (nrElec - n_before)

      call Update_Position(i)

      ! The Ramo current of this step, and the absorption counts marked
      ! during the position update (Remove_Particles clears them)
      Q_ramo = Q_ramo + ramo_current(species_elec)*time_step
      if (i > n_steps/2) then
        Q_ramo_ss = Q_ramo_ss + ramo_current(species_elec)*time_step
      end if
      n_abs_top = n_abs_top + nrElec_remove_top
      n_abs_bot = n_abs_bot + nrElec_remove_bot

      call Remove_Particles(i)

      if (nrPart > 0) then
        if (.not. (all(particles_cur_pos(3, 1:nrPart) >= 0.0d0) .and. &
                 & all(particles_cur_pos(3, 1:nrPart) <= box_dim(3)))) bounds_ok = .false.
        ! A NaN also fails this comparison
        if (.not. all(abs(particles_cur_pos(:, 1:nrPart)) < 1.0d0)) finite_ok = .false.
      end if
      n_peak = max(n_peak, nrElec)
    end do

    ! Expected Ramo charge: every electron contributes q_0 * dz/d, absorbed
    ! electrons traveled from z_emit to the anode (or back to the cathode),
    ! electrons still in flight from z_emit to their current position.
    sum_z = 0.0d0
    do i = 1, nrPart
      sum_z = sum_z + (particles_cur_pos(3, i) - z_emit)
    end do
    Q_exp = q_0*( n_abs_top*(d_test - z_emit) - n_abs_bot*z_emit + sum_z )/d_test

    I_avg_ss = Q_ramo_ss / (time_step*(n_steps - n_steps/2))

    print '(a, i0, a, i0, a, i0, a, i0)', '  emitted = ', n_emit_tot, ', absorbed top = ', n_abs_top, &
                                        & ', absorbed bottom = ', n_abs_bot, ', in flight = ', nrElec
    print '(a, i0)', '  peak number of electrons = ', n_peak
    print '(a, f8.4, a)', '  steady-state current = ', I_avg_ss*1.0d3, ' mA'
    print '(a, f8.4, a, f8.4, a)', '  Ramo charge = ', Q_ramo*1.0d15, ' fC, expected = ', Q_exp*1.0d15, ' fC'

    call Assert_True((n_emit_tot > 100), 'Planar system: emission is active')
    call Assert_True((n_abs_top > 50), 'Planar system: electrons reach the anode')
    call Assert_True((nrElec == n_emit_tot - n_abs_top - n_abs_bot), &
                   & 'Planar system: charge bookkeeping is conserved')
    call Assert_True(bounds_ok, 'Planar system: no particle outside the gap')
    call Assert_True(finite_ok, 'Planar system: positions stay finite')
    call Assert_True((abs(Q_ramo - Q_exp) < 0.10d0*abs(Q_exp)), &
                   & 'Planar system: Ramo charge accounting (within 10%)')
    call Assert_True((I_avg_ss > I_low) .and. (I_avg_ss < I_high), &
                   & 'Planar system: steady-state current in the expected band')

    ! Clean up the emission module state allocated by Init_Field_Emission_v2
    ! and leave a valid work function pointer behind.
    call Clean_Up_Field_Emission_v2()
    ptr_Work_fun => Const_Work_Fun_Test

    print *, 'Planar field emission system test finished'
    print *, ''
  end subroutine Test_Planar_System

  !-----------------------------------------------------------------------------
  ! Whole system test for the hyperboloid tip: run the production tip field
  ! emission (Init_Emission_Tip / Do_Emission_Tip, surface supply integral
  ! plus Metropolis-Hastings on the tip surface) together with the integrator
  ! and the absorption, exactly like the main loop does. The checks mirror
  ! the planar system test. For the Ramo accounting the tip weighting field
  ! (E_zunit_tip) is not uniform, so instead of an exact identity the
  ! Shockley-Ramo bounds are used: every electron absorbed at the anode
  ! contributes q_0, every electron still in flight between 0 and q_0, and
  ! electrons that returned to the tip contribute nothing.
  subroutine Test_Tip_System()
    use mod_emission_tip, only: Init_Emission_Tip

    double precision, parameter :: d_test = 1000.0d0*length_scale ! Total gap (tip + vacuum)
    ! The apex field of this tip is about 9.8 MV/m per volt (49 nm apex
    ! radius). 800 V gives about 7.9 GV/m at the apex: strong emission, but
    ! still below the field where the image charge correction saturates
    ! (l > 1) and the barrier is fully suppressed.
    double precision, parameter :: V_test = 800.0d0 ! V
    double precision, parameter :: delta_t_test = 0.25d-15 ! Time step
    integer, parameter          :: n_steps = 350

    integer          :: i, n_before, n_emit_tot, n_abs_top, n_abs_bot, n_peak
    double precision :: Q_ramo, Q_transits
    logical          :: bounds_ok, finite_ok

    print *, 'Starting tip field emission system test'

    call Setup_Test_System(d_test, delta_t_test, n_steps, V_test, .true., 0)

    ! The same tip as in Test_Tip_Acceleration: 1000 nm gap with a 100 nm
    ! high tip, 900 nm from apex to anode, 100 nm base radius.
    nrEmit = 1
    emitters_dim(1:3, 1) = (/ 900.0d0, 100.0d0, 100.0d0 /) * length_scale
    emitters_pos(1:3, 1) = 0.0d0
    ! NOTE: the tip module has its own emitter type convention, 1 means
    ! field emission (it is not EMIT_CIRCLE).
    emitters_type(1)     = 1
    emitters_delay(1)    = 0

    ! The production init: tip geometry parameters and pointers
    call Init_Emission_Tip()

    call Seed_RNG_Deterministic()

    n_emit_tot = 0
    n_abs_top  = 0
    n_abs_bot  = 0
    n_peak     = 0
    Q_ramo     = 0.0d0
    bounds_ok  = .true.
    finite_ok  = .true.

    ! The same sequence as the main loop
    do i = 1, n_steps
      n_before = nrElec
      call ptr_Do_Emission(i)
      n_emit_tot = n_emit_tot + (nrElec - n_before)

      call Update_Position(i)

      Q_ramo = Q_ramo + ramo_current(species_elec)*time_step
      n_abs_top = n_abs_top + nrElec_remove_top
      n_abs_bot = n_abs_bot + nrElec_remove_bot

      call Remove_Particles(i)

      if (nrPart > 0) then
        if (.not. (all(particles_cur_pos(3, 1:nrPart) >= 0.0d0) .and. &
                 & all(particles_cur_pos(3, 1:nrPart) <= box_dim(3)))) bounds_ok = .false.
        ! A NaN also fails this comparison
        if (.not. all(abs(particles_cur_pos(:, 1:nrPart)) < 1.0d0)) finite_ok = .false.
      end if
      n_peak = max(n_peak, nrElec)
    end do

    ! Ramo charge in units of full cathode-anode transits
    Q_transits = Q_ramo/q_0

    print '(a, i0, a, i0, a, i0, a, i0)', '  emitted = ', n_emit_tot, ', absorbed top = ', n_abs_top, &
                                        & ', absorbed bottom = ', n_abs_bot, ', in flight = ', nrElec
    print '(a, i0)', '  peak number of electrons = ', n_peak
    print '(a, f10.3, a, i0, a, i0)', '  Ramo charge = ', Q_transits, &
                                    & ' transits, bounds: > 0.8 x ', n_abs_top, ', < in flight + ', n_abs_top

    call Assert_True((n_emit_tot > 20), 'Tip system: emission is active')
    call Assert_True((n_abs_top > 5), 'Tip system: electrons reach the anode')
    call Assert_True((nrElec == n_emit_tot - n_abs_top - n_abs_bot), &
                   & 'Tip system: charge bookkeeping is conserved')
    call Assert_True(bounds_ok, 'Tip system: no particle outside the gap')
    call Assert_True(finite_ok, 'Tip system: positions stay finite')
    call Assert_True((Q_transits > 0.8d0*n_abs_top - 2.0d0), &
                   & 'Tip system: Ramo charge above the Shockley-Ramo lower bound')
    call Assert_True((Q_transits < 1.0d0*(n_abs_top + nrElec) + 2.0d0), &
                   & 'Tip system: Ramo charge below the Shockley-Ramo upper bound')

    print *, 'Tip field emission system test finished'
    print *, ''
  end subroutine Test_Tip_System

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
