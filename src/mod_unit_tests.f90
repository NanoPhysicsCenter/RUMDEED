!-------------------------------------------!
! Module for unit testing                   !
! Kristinn Torfason                         !
! 10.09.17                                  !
!-------------------------------------------!

module mod_unit_tests
  use mod_global
  use mod_pair
  use mod_verlet
  implicit none

  double precision, parameter :: tolerance = 1.0d-10
contains
  subroutine Run_Unit_Tests()
    call Test_Acceleration()

    call Test_Transit_Time()
  end subroutine Run_Unit_Tests

  !-----------------------------------------------------------------------------
  ! Initialize an empty system with the given parameters for testing
  subroutine Setup_Test_System(d_test, delta_t_test, steps_test, V_test)
    double precision, intent(in) :: d_test, delta_t_test, V_test
    integer, intent(in)          :: steps_test

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
    V = V_test

    ! Init
    d = box_dim(3)
    E_z = -1.0d0*V/d
    V_a = V
    E_zunit = -1.0d0/d

    dens_x_d = box_dim(1) / (N_x_densmap-1)
    dens_y_d = box_dim(2) / (N_y_densmap-1)

    ! Set the current scale to mA / cm^2
    cur_scale = 1.0d-3 * (box_dim(1)*1.0d2) * (box_dim(3)*1.0d2) ! mA / cm^2

    !time_step: The size of the time step
    !time_step = time_step * time_scale
    time_step2 = time_step**2 ! Time_step squared
  end subroutine Setup_Test_System

  !-----------------------------------------------------------------------------
  ! Test the Acceleration calculations
  ! Two electrons and one hole are placed in the system
  ! The acceleration is calculated using Coulomb's law and also with the
  ! Calculate Acceleration subroutine. The two are then compared to see if the
  ! optimizations and the OpenMP are working correctly.
  subroutine Test_Acceleration()
    double precision, dimension(1:3) :: R_1, R_2, R_3
    double precision, dimension(1:3) :: E_test, par_vel, E_pos

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
    !$OMP MASTER
    print *, 'Running test for acceleration'

    R_1 = (/  3.0d0, -10.0d0, 101.0d0 /) * length_scale
    R_3 = (/  6.0d0, -24.0d0, 118.0d0 /) * length_scale
    R_2 = (/ -9.0d0,  26.0d0,  80.0d0 /) * length_scale

    E_test = (/ 0.0d0, 0.0d0, -1.0d0*V_test/d_test /)

    ! Particle 1
    a_12 = +1.0d0*pre_fac_c_test * (R_1 - R_2) / norm2(R_1 - R_2)**3
    a_13 = -1.0d0*pre_fac_c_test * (R_1 - R_3) / norm2(R_1 - R_3)**3
    a_1E = -1.0d0*pre_fac_E_test * E_test

    a_1 = a_12 + a_13 + a_1E

    ! Particle 2
    a_21 = +1.0d0*pre_fac_c_test * (R_2 - R_1) / norm2(R_2 - R_1)**3
    a_23 = -1.0d0*pre_fac_c_test * (R_2 - R_3) / norm2(R_2 - R_3)**3
    a_2E = -1.0d0*pre_fac_E_test * E_test

    a_2 = a_21 + a_23 + a_2E

    ! Particle 3
    a_31 = -1.0d0*pre_fac_c_test * (R_3 - R_1) / norm2(R_3 - R_1)**3
    a_32 = -1.0d0*pre_fac_c_test * (R_3 - R_2) / norm2(R_3 - R_2)**3
    a_3E = +1.0d0*pre_fac_E_test * E_test

    a_3 = a_31 + a_32 + a_3E

    ! Set input variables
    call Setup_Test_System(d_test, delta_t_test, 100, V_test)

    ! Add particles
    par_vel = 0.0d0
    call Add_Particle(R_1, par_vel, species_elec, 1)
    call Add_Particle(R_2, par_vel, species_elec, 1)
    call Add_Particle(R_3, par_vel, species_hole, 1)

    !$OMP END MASTER

    !$OMP BARRIER
    call Calculate_Acceleration_Particles()

    par_vel = (/ -4.55, -2.34, 96.44 /) * length_scale
    E_pos = Calc_Field_at(par_vel)

    !$OMP MASTER
    a_1_res = particles_cur_accel(:, 1)
    a_2_res = particles_cur_accel(:, 2)
    a_3_res = particles_cur_accel(:, 3)

    if (all(abs(a_1 - a_1_res)/a_1 < tolerance)) then
      print *, 'Particle 1 PASSED'
    else
      print *, 'Particle 1 FAILED'
      print *, a_1_res
      print *, a_1
      print *, ''
    end if

    if (all(abs(a_2 - a_2_res)/a_2 < tolerance)) then
      print *, 'Particle 2 PASSED'
    else
      print *, 'Particle 2 FAILED'
      print *, 'Results'
      print *, a_2_res
      print *, 'Expected'
      print *, a_2
      print *, 'Difference'
      print *, abs(a_2 - a_2_res)
      print *, all(abs(a_2 - a_2_res) < tolerance)
      print *, ''
    end if

    if (all(abs(a_3 - a_3_res)/a_3 < tolerance)) then
      print *, 'Particle 3 PASSED'
    else
      print *, 'Particle 3 FAILED'
      print *, a_2_res
      print *, a_2
      print *, ''
    end if

    print *, 'Electric field at'
    print *, 'E_pos = ', E_pos

    print *, 'Acceleration test finished'
    print *, ''

    !$OMP END MASTER
  end subroutine Test_Acceleration

  !-----------------------------------------------------------------------------
  ! We can compute the transit time of a single electron in the system using
  ! Z = Z_0 + V_0*t + 1/2*a*t^2, where a = qV/md, V_0 = 0, Z_0 = 1.0 nm and
  ! Z = d. This gives t = sqrt(2*m*d*(d-z_0)/(q*V)). If we dived this number
  ! with the time step and round up. We should know how many time steps should
  ! pass before the electron is absorbed. This test is to see if the
  ! Verlet algorithm is working properly.
  subroutine Test_Transit_Time()
    double precision, parameter      :: d_test = 100.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d0 ! V
    double precision, parameter      :: delta_t_test = 0.25d-15 ! Time step

    double precision                 :: time_exp ! Expected time
    integer                          :: steps_exp ! Exptected number of times steps

    integer                          :: i, steps_res

    logical                          :: found_res = .false. ! Set to true if
                                                            ! the electron has
                                                            ! been absorbed

    double precision, dimension(1:3) :: R_1, par_vel

    !$OMP MASTER
    print *, 'Starting transite time test'
    R_1 = (/ 0.0d0, 0.0d0, 1.0d0*length_scale /)

    time_exp = sqrt(2.0d0*d_test*(d_test - R_1(3))*m_0/(q_0*V_test))
    steps_exp = ceiling(time_exp / delta_t_test)

    ! Set input variables
    call Setup_Test_System(d_test, delta_t_test, 1000, V_test)

    ! Add particle
    par_vel = 0.0d0
    call Add_Particle(R_1, par_vel, species_elec, 1)
    !$OMP END MASTER
    !$OMP BARRIER

    do i = 1, steps

      ! Do Emission
      !call emission code

      ! Update the position of all particles
      !print *, 'Update position'
      call Update_Position(i)

      ! Remove particles from the system
      call Remove_Particles(i)
      !$OMP MASTER
      !print *, 'i = ', i
      !print *, 'nrPart = ', nrPart
      if ((nrPart == 0) .and. (found_res .eqv. .false.)) then
        steps_res = i
        found_res = .true.
      end if
      !$OMP END MASTER
      !$OMP BARRIER
    end do

    !$OMP MASTER
      do i = 1, MAX_LIFE_TIME
        if (life_time(i, species_elec) /= 0) then
          print *, i
          print *, life_time(i, species_elec)
          print *, ''
        end if
      end do

      print *, 'Transite time = ', steps_res
      if (steps_res == steps_exp) then
        print *, 'Transite time test PASSED'
      else
        print *, 'Transite time test FAILED'
      end if

      print *, 'Transit time test finished'
      print *, ''
    !$OMP END MASTER
  end subroutine Test_Transit_Time
end module mod_unit_tests
