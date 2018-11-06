!
! Kristinn Torfason
! 05.11.2018
! A program that does unit tests on the Verlet module
!
program Verlet_Tests
    use mod_global
    use mod_verlet
    use mod_pair

    implicit none

    call Init()
    call Test_Transit_Time()
    call Clean_up()
contains
    subroutine Init()
        integer :: istat, ud_bin, ud_text

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

        ! Allocate arrays
        allocate(particles_cur_pos(1:3, 1:MAX_PARTICLES))
        allocate(particles_prev_pos(1:3, 1:MAX_PARTICLES))
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

        allocate(life_time(1:MAX_LIFE_TIME, 1:2))
        allocate(ramo_current(1:nrSpecies))

        open(newunit=ud_text, file='/dev/null', status='REPLACE', action='WRITE', iostat=istat) ! Linux / Unix
        if (istat /= 0) then
            open(newunit=ud_text, file='NUL', status='REPLACE', action='WRITE', iostat=istat) ! Windows
            if (istat /= 0) then
                print *, 'ERROR: Unable to redirect file output'
            end if
        end if

        open(newunit=ud_bin, file='/dev/null', status='REPLACE', action='WRITE', iostat=istat, access='STREAM') ! Linux / Unix
        if (istat /= 0) then
            open(newunit=ud_bin, file='NUL', status='REPLACE', action='WRITE', iostat=istat, access='STREAM') ! Windows
            if (istat /= 0) then
                print *, 'ERROR: Unable to redirect file output'
            end if
        end if

        ud_pos = ud_bin
        ud_emit = ud_text
        ud_absorb = ud_text
        ud_absorb_top = ud_text
        ud_absorb_bot = ud_text
        ud_ramo = ud_text
        ud_volt = ud_text
        ud_debug = ud_text
        ud_field = ud_text
        ud_ramo_sec = ud_text

        ud_density_emit = ud_bin

        ud_density_absorb_top = ud_bin
        ud_density_absorb_bot = ud_bin
    end subroutine Init

    subroutine Clean_up()
      ! Deallocate arrays
      deallocate(particles_cur_pos)
      deallocate(particles_prev_pos)
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

    !-----------------------------------------------------------------------------
    ! Initialize an empty system with the given parameters for testing
    subroutine Setup_Test_System(d_test, delta_t_test, steps_test, V_test, ic, ic_max)
        double precision, intent(in) :: d_test, delta_t_test, V_test
        logical, intent(in)          :: ic
        integer, intent(in)          :: steps_test, ic_max

        ! Clear variables
        particles_cur_pos = 0.0d0
        particles_prev_pos = 0.0d0
        particles_cur_vel = 0.0d0
        particles_cur_accel = 0.0d0
        particles_prev_accel = 0.0d0
        particles_prev2_accel = 0.0d0
        particles_charge = 0.0d0
        particles_species = 0
        particles_mass = 0.0d0
        particles_step = 0
        particles_mask = .true.

        ramo_current = 0.0d0
        life_time = 0

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
        ! System size is x: 100nm, y: 100nm and z: X d_test
        box_dim = (/ 100.0d0*length_scale, 100.0d0*length_scale, d_test /)

        ! The time step
        time_step = delta_t_test
        time_step2 = time_step**2 ! Time_step squared

        ! Number of time steps
        steps = steps_test
        ! Gap spaceing
        d = box_dim(3)

        ! Voltage
        V_s = V_test
        V_d = V_s

        ! Electric field (planar)
        E_z = -1.0d0*V_d/d
        E_zunit = -1.0d0/d

        ! Image charge
        image_charge = ic
        N_ic_max = ic_max
    end subroutine Setup_Test_System

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
        call Setup_Test_System(d_test, delta_t_test, steps_exp+1000, V_test, .true., 2)
    
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
end program Verlet_Tests