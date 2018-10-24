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

  double precision, parameter :: tolerance = 1.0d-6
contains
  subroutine Run_Unit_Tests()
    call Test_Acceleration()

    call Test_Transit_Time()

    call Test_Many_Particles()
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
    V_s = V_test

    ! Init
    d = box_dim(3)
    V_d = V_s
    E_z = -1.0d0*V_d/d
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
    print *, 'Running test for acceleration'

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
    call Setup_Test_System(d_test, delta_t_test, 100, V_test)

    ! Add particles
    par_vel = 0.0d0
    call Add_Particle(R_1, par_vel, species_elec, 1, 0)
    call Add_Particle(R_2, par_vel, species_elec, 1, 0)
    call Add_Particle(R_3, par_vel, species_hole, 1, 0)

    call Calculate_Acceleration_Particles()

    ! Test the field at this location
    par_vel = (/ -4.55, -2.34, 96.44 /) * length_scale
    E_pos = Calc_Field_at(par_vel)

    ! Results from the acceleration
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
      print *, a_3_res
      print *, a_3
      print *, ''
    end if

    !E_python = (/ 6192429.94450208d0,  -4866704.16373579d0, -17454923.69763095d0 /) ! Results from Python script (no image charge)
    E_python = (/ -102526.05673208d0, 421022.84663293d0, -20456407.74487634d0 /)
    if (all(abs(E_python - E_pos)/E_python < tolerance)) then
      print *, 'Electric field PASSED'
    else
      print *, 'Electric field FAILED'
      print *, 'E_python = ', E_python
      print *, 'E_pos = ', E_pos
      print *, 'abs(E_python - E_pos)/E_python = ', abs(E_python - E_pos)/E_python
      print *, ''
    end if

    print *, 'Acceleration test finished'
    print *, ''
  end subroutine Test_Acceleration

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
    call Setup_Test_System(d_test, delta_t_test, steps_exp+1000, V_test)

    ! Add particle
    par_vel = 0.0d0
    call Add_Particle(R_1, par_vel, species_elec, 1, 1)

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
    call Setup_Test_System(d_test, delta_t_test, 1000, V_test)

    ! Open the file with initial positions of particles
    open(newunit=ud_data, iostat=IFAIL, file='rand_pos_init.dt', status='OLD', action='read')
    if (IFAIL /= 0) then
      print '(a)', 'Vacuum: ERROR UNABLE TO OPEN file input'
      stop
    end if

    do i = 1, n_par
      ! Set the intial position and velocity
      read(unit=ud_data, fmt=*) par_pos ! Read the next line
      par_vel = 0.0d0

      ! Add the particle to the system
      call Add_Particle(par_pos, par_vel, species_elec, 1, 0)
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
