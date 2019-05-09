!-------------------------------------------!
! Module for unit testing                   !
! Kristinn Torfason                         !
! 10.09.17                                  !
! Just some simple tests to check if        !
! some things are working correctly         !
!-------------------------------------------!

module mod_unit_tests
  use mod_global
  use mod_pair
  use mod_verlet
  implicit none

  double precision, parameter :: tolerance_rel = 0.02d0 ! 2% relative error tolerance
  double precision, parameter :: tolerance_abs = 1.0d-6 ! absolute error tolerance 
contains
  subroutine Run_Unit_Tests()
    call Test_Acceleration_Without_Image_Charge()

    call Test_Image_Charge()

    !call Test_Acceleration_With_Image_Charge()

    !call Test_Transit_Time()

    !call Test_Many_Particles()
  end subroutine Run_Unit_Tests

  !-----------------------------------------------------------------------------
  ! Initialize an empty system with the given parameters for testing
  subroutine Setup_Test_System(d_test, delta_t_test, steps_test, V_test, use_ic, N_ic)
    double precision, intent(in) :: d_test, delta_t_test, V_test
    integer, intent(in)          :: steps_test, N_ic
    logical, intent(in)          :: use_ic

    ! Clear variables
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

    ! Set input variables
    box_dim = (/ 100.0d0*length_scale, 100.0d0*length_scale, d_test /)
    time_step = delta_t_test
    steps = steps_test
    V_s = V_test

    ! Init
    d = box_dim(3)
    V_d = V_s
    E_z = -1.0d0*V_d/d
    E_zunit = -1.0d0/d

    !dens_x_d = box_dim(1) / (N_x_densmap-1)
    !dens_y_d = box_dim(2) / (N_y_densmap-1)

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

    double precision, parameter      :: pre_fac_c_test = q_0**2/(4.0d0*pi*m_0*epsilon_0)
    double precision, parameter      :: pre_fac_E_test = q_0/m_0

    ! We must use OMP MASTER here because the variables declared above
    ! are private for each thread and are NOT shared.
    print *, 'Running test for acceleration without image charge'

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

    ! Acceleration from the electric field in the diode
    a_1E = -1.0d0*pre_fac_E_test * E_test

    ! Total acceleration 1
    a_1 = a_12 + a_13 + a_1E

    ! Particle 2
    ! Acceleration from other particles
    a_21  = +1.0d0*pre_fac_c_test * (R_2 - R_1) / norm2(R_2 - R_1)**3
    a_23  = -1.0d0*pre_fac_c_test * (R_2 - R_3) / norm2(R_2 - R_3)**3

    ! Acceleration from the electric field in the diode
    a_2E = -1.0d0*pre_fac_E_test * E_test

    ! Total aceleration on particle 2
    a_2 = a_21 + a_23 + a_2E

    ! Particle 3
    ! Acceleration from other particles
    a_31  = -1.0d0*pre_fac_c_test * (R_3 - R_1) / norm2(R_3 - R_1)**3
    a_32  = -1.0d0*pre_fac_c_test * (R_3 - R_2) / norm2(R_3 - R_2)**3

    ! Acceleration from the electric field in the diode
    a_3E = +1.0d0*pre_fac_E_test * E_test

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
    call Add_Particle(R_3, par_vel, species_hole, 1, 0, -1)

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

    if (abs(E_test(3) - E_F(3))/E_test(3) < tolerance_rel) then
      if ( (abs(E_test(1) - 0.0d0) < tolerance_abs) .and. (abs(E_test(2) - 0.0d0) < tolerance_abs) ) then
        print *, 'Electric field without particles PASSED'
      else
        print *, 'Electric field without particles FAILED'
        print *, 'x and y component not zero'
      end if
    else
      print *, 'Electric field without particles FAILED'
      print *, 'z component not correct'
      print *, 'E_test = ', E_test
      print *, 'E_F = ', E_F
      print *, 'abs(E_test - E_F)/E_test = ', abs(E_test - E_F)/E_test
      print *, ''
    end if

    E_python = (/ -314559.29097098, 1423979.07058996, -20246038.87978313 /) ! Results from Python script (no image charge)
    if (all(abs(E_python - E_pos)/E_python < tolerance_rel)) then
      print *, 'Electric field with particles PASSED'
    else
      print *, 'Electric field with particles FAILED'
      print *, 'E_python = ', E_python
      print *, 'E_pos = ', E_pos
      print *, 'abs(E_python - E_pos)/E_python = ', abs(E_python - E_pos)/E_python
      print *, ''
    end if

    print *, 'Acceleration test without image charge finished'
    print *, ''
  end subroutine

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
    ! and hole at,
    R_2 = (/  6.0d0, -24.0d0, 98.53d0 /) * length_scale

    ! Let's test the self interaction first
    force_ic_11 = q_1**2 * div_fac_c * Force_Image_charges_v2(R_1, R_1)

    ! The force should be given by Coulombs law
    ! F = 1/(4*pi*epsilon_0) * (r - r')/|r - r'|^3
    ! Here r is the location of the particle and r' is the location of the image charge
    ! The image charge is located below emitter. The force should only be in the z-direction,
    ! F_R1 = 1/(4*pi*epsilon_0) * (-q_1)*(+q_1)*(z - z')/|z - z'|^3 .
    ! We know that z' = -z so z - z' = z - (-z) = 2z, or,
    ! F_R1 = -q_1**2/(4*pi*epsilon_0) * 2z/(2z)^(3/2) 1/(4*pi*epsilon_0) * 1/sqrt(2z)

    F_R1 = 0.0d0
    F_R1(3) = -1.0d0*q_1**2/(4.0d0*pi*epsilon_0)*1.0d0/sqrt(2.0d0*R_1(3))

    if (abs(force_ic_11(3) - F_R1(3))/F_R1(3) < tolerance_rel) then
      if ( (abs(force_ic_11(1) - 0.0d0) < tolerance_abs) .and. (abs(force_ic_11(2) - 0.0d0) < tolerance_abs) ) then
        print *, 'Self interaction on particle 1 PASSED'
      else
        print *, 'Self interaction on particle 1 FAILED'
        print *, 'x and y component are not zero'
        print *, 'force_ic_11(1) = ', force_ic_11(1)
        print *, 'force_ic_11(2) = ', force_ic_11(2)
      endif
    else
      print *, 'Self interaction on image charge particle 1 FAILED'
      print *, 'z component not correct'
      print *, 'F_R1 = ', F_R1(3)
      print *, 'force_ic_11 = ', force_ic_11(3)
      print *, 'abs(force_ic_11 - F_R1)/F_R1 = ', abs(F_R1(3) - force_ic_11(3))/F_R1(3)
      print *, ''
    end if

    force_ic_22 = q_2**2 * div_fac_c * Force_Image_charges_v2(R_1, R_1)

    ! F_R2 = -q_1**2/(4*pi*epsilon_0) * 2z/(2z)^(3/2) 1/(4*pi*epsilon_0) * 1/sqrt(2z)

    F_R2 = 0.0d0
    F_R2(3) = -1.0d0*q_2**2/(4.0d0*pi*epsilon_0)*1.0d0/sqrt(2.0d0*R_2(3))

    if (abs(force_ic_22(3) - F_R2(3))/F_R2(3) < tolerance_rel) then
      if ( (abs(force_ic_22(1) - 0.0d0) < tolerance_abs) .and. (abs(force_ic_22(2) - 0.0d0) < tolerance_abs) ) then
        print *, 'Self interaction on particle 2 PASSED'
      else
        print *, 'Self interaction on particle 2 FAILED'
        print *, 'x and y component are not zero'
        print *, 'force_ic_22(1) = ', force_ic_22(1)
        print *, 'force_ic_22(2) = ', force_ic_22(2)
      endif
    else
      print *, 'Self interaction on image charge particle 2 FAILED'
      print *, 'z component not correct'
      print *, 'F_R2 = ', F_R2(3)
      print *, 'force_ic_22 = ', force_ic_22(3)
      print *, 'abs(force_ic_22 - F_R1)/F_R1 = ', abs(F_R2(3) - force_ic_22(3))/F_R2(3)
      print *, ''
    end if

    ! Test image charge force on particle 1 from particle 2
    ! Image charge of particle 2 is located at -z
    R_2_IC = R_2
    R_2_IC(3) = -1.0d0*R_2_Ic(3)

    F_IC_12 = q_1*q_2_ic/(4.0d0*pi*epsilon_0) * (R_1 - R_2_IC) / norm2(R_1 - R_2_IC)**3

    force_ic_12 = q_1*q_2 * div_fac_c * Force_Image_charges_v2(R_1, R_2)

    if (all(abs(force_ic_12 - F_IC_12)/F_IC_12 < tolerance_rel)) then
      print *, 'Image charge interaction on particle 1 from particle 2 PASSED'
    else
      print *, 'Image charge interaction on particle 1 from particle 2 FAILED'
      print *, 'F_IC_12 = ', F_IC_12(3)
      print *, 'force_ic_12 = ', force_ic_12(3)
      print *, 'abs(force_ic_12 - F_IC_12)/F_IC_12 = ', abs(force_ic_12 - F_IC_12)/F_IC_12
      print *, ''
    end if

    ! Test image charge force on particle 2 from particle 1
    ! Image charge of particle 1 is located at -z
    R_1_IC = R_1
    R_1_IC(3) = -1.0d0*R_1_Ic(3)

    F_IC_21 = q_2*q_1_ic/(4.0d0*pi*epsilon_0) * (R_2 - R_1_IC) / norm2(R_2 - R_1_IC)**3

    force_ic_21 = q_2*q_1 * div_fac_c * Force_Image_charges_v2(R_2, R_1)

    if (all(abs(force_ic_21 - F_IC_21)/F_IC_21 < tolerance_rel)) then
      print *, 'Image charge interaction on particle 2 from particle 1 PASSED'
    else
      print *, 'Image charge interaction on particle 2 from particle 1 FAILED'
      print *, 'F_IC_21 = ', F_IC_21
      print *, 'force_ic_21 = ', force_ic_21
      print *, 'abs(force_ic_21 - F_IC_21)/F_IC_21 = ', abs(force_ic_21 - F_IC_21)/F_IC_21
      print *, ''
    end if

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

    if (abs(force_ic_11(3) - F_R1(3))/F_R1(3) < tolerance_rel) then
      if ( (abs(force_ic_11(1) - 0.0d0) < tolerance_abs) .and. (abs(force_ic_11(2) - 0.0d0) < tolerance_abs) ) then
        print *, 'Self interaction on particle 1 (N_IC_MAX = 1) PASSED'
      else
        print *, 'x and y component are not zero'
        print *, 'force_ic_11(1) = ', force_ic_11(1)
        print *, 'force_ic_11(2) = ', force_ic_11(2)
      end if
    else
      print *, 'Self interaction on particle 1 (N_IC_MAX = 1) FAILED'
      print *, 'z component if not correct'
      print *, 'F_R1(3) = ', F_R1(3)
      print *, 'force_ic_11(3) = ', force_ic_11(3)
      print *, 'abs(force_ic_11(3) - F_R1(3))/F_R1(3) = ', abs(force_ic_11(3) - F_R1(3))/F_R1(3)
      print *, ''
    end if

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

    if (abs(force_ic_22(3) - F_R2(3))/F_R2(3) < tolerance_rel) then
      if ( (abs(force_ic_22(1) - 0.0d0) < tolerance_abs) .and. (abs(force_ic_22(2) - 0.0d0) < tolerance_abs) ) then
        print *, 'Self interaction on particle 2 (N_IC_MAX = 1) PASSED'
      else
        print *, 'x and y component are not zero'
        print *, 'force_ic_11(1) = ', force_ic_22(1)
        print *, 'force_ic_11(2) = ', force_ic_22(2)
      end if
    else
      print *, 'Self interaction on particle 2 (N_IC_MAX = 1) FAILED'
      print *, 'z component if not correct'
      print *, 'F_R2 = ', F_R2
      print *, 'force_ic_22 = ', force_ic_22
      print *, 'abs(force_ic_22(3) - F_R2(3))/F_R2(3) = ', abs(force_ic_22(3) - F_R2(3))/F_R2(3)
      print *, ''
    end if

    !--------------------------------------------------------------------------------------
    ! Testing interaction on particle 1 from particle 2 with N_IC_MAX = 1

    F_IC_11 = q_1*q_2_IC_1/(4.0d0*pi*epsilon_0) * (R_1 - R_2_IC_1) / norm2(R_1 - R_2_IC_1)**3
    F_IC_12 = q_1*q_2_IC_2/(4.0d0*pi*epsilon_0) * (R_1 - R_2_IC_2) / norm2(R_1 - R_2_IC_2)**3
    F_IC_13 = q_1*q_2_IC_3/(4.0d0*pi*epsilon_0) * (R_1 - R_2_IC_3) / norm2(R_1 - R_2_IC_3)**3
    F_IC_14 = q_1*q_2_IC_4/(4.0d0*pi*epsilon_0) * (R_1 - R_2_IC_4) / norm2(R_1 - R_2_IC_4)**3
    F_IC_15 = q_1*q_2_IC_5/(4.0d0*pi*epsilon_0) * (R_1 - R_2_IC_5) / norm2(R_1 - R_2_IC_5)**3

    F_R1 = F_IC_11 + F_IC_12 + F_IC_13 + F_IC_14 + F_IC_15

    force_ic_12 = q_1*q_2 * div_fac_c * Force_Image_charges_v2(R_1, R_2)

    if (all(abs(force_ic_12 - F_R1)/F_R1 < tolerance_rel)) then
      print *, 'Image charge interaction on particle 1 from particle 2 (N_IC_MAX = 1) PASSED'
    else
      print *, 'Image charge interaction on particle 1 from particle 2 (N_IC_MAX = 1) FAILED'
      print *, 'F_IC_21 = ', F_R1
      print *, 'force_ic_12 = ', force_ic_12
      print *, 'abs(force_ic_12 - F_R1)/F_R1 = ', abs(force_ic_12 - F_R1)/F_R1
      print *, ''
    end if

    print *, 'Image charge function "Force_Image_charges_v2(pos_1, pos_2)" test finished'
    print *, ''
  end subroutine Test_Image_Charge

  !-----------------------------------------------------------------------------
  ! Test the Acceleration calculations
  ! Two electrons and one hole are placed in the system
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

    double precision, parameter      :: pre_fac_c_test = q_0**2/(4.0d0*pi*m_0*epsilon_0)
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
    ! Other hole (-q*-q = +q**2)
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
    ! Other hole (-q*-q = +q**2)
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
    call Add_Particle(R_3, par_vel, species_hole, 1, 0, -1)

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

    print *, 'Starting transite time test'
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

      print *, 'Transite time = ', steps_res
      print *, 'Exptected time = ', steps_exp
      if (abs(steps_exp - steps_res)/steps_exp < 0.01) then
      !if (steps_res == steps_exp) then
        print *, 'Transite time test PASSED'
        print *, 'Difference was less than 1% or ', abs(steps_res - steps_exp)
      else
        print *, 'Transite time test FAILED'
        print *, 'Expected = ', steps_exp
        print *, 'Results = ', steps_res
      end if

      print *, 'Transit time test finished'
      print *, ''
  end subroutine Test_Transit_Time

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
      print '(a)', 'Vacuum: ERROR UNABLE TO OPEN file input'
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
end module mod_unit_tests
