!-------------------------------------------!
! Module for unit tests                     !
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
  end subroutine Run_Unit_Tests

  subroutine Test_Acceleration()
    double precision, dimension(1:3) :: R_1, R_2, R_3
    double precision, dimension(1:3) :: E_test, par_vel

    double precision, dimension(1:3) :: a_1, a_12, a_13, a_1E, a_1_res
    double precision, dimension(1:3) :: a_2, a_21, a_23, a_2E, a_2_res
    double precision, dimension(1:3) :: a_3, a_31, a_32, a_3E, a_3_res

    double precision, parameter      :: d_test = 100.0d0*length_scale
    double precision, parameter      :: V_test = 2.0d0 ! V
    double precision, parameter      :: delta_t_test = 0.25d-12 ! Time step

    double precision, parameter      :: pre_fac_c_test = q_0**2/(4.0d0*pi*m_0*epsilon_0)
    double precision, parameter      :: pre_fac_E_test = q_0/m_0

    !$OMP FLUSH
    !$OMP MASTER

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
    box_dim = (/ 100.0d0, 100.0d0, d_test /)
    time_step = delta_t_test
    steps = 100
    T = 5.0d0
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
    time_step = time_step * time_scale
    time_step2 = time_step**2 ! Time_step squared

    ! Add particles
    par_vel = 0.0d0
    call Add_Particle(R_1, par_vel, species_elec, 0)
    call Add_Particle(R_2, par_vel, species_elec, 0)
    call Add_Particle(R_3, par_vel, species_hole, 0)

    !$OMP END MASTER

    print *, 'Hi'
    !$OMP BARRIER

    !!$OMP SINGLE

    call Calculate_Acceleration_Particles()

    !!$OMP END SINGLE

    print *, 'Hi2'
    !$OMP BARRIER
    print *, 'Hi3'

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

    !$OMP END MASTER
    print *, 'Hi4'
  end subroutine Test_Acceleration
end module mod_unit_tests
