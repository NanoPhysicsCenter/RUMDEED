!
! Kristinn Torfason
! 05.11.2018
! A program that does unit tests on the Verlet module
!
program Verlet_Tests
    implicit none
    use mod_global
    use mod_verlet
    use mod_pair

    call Init()
    call Transit_time_test()
    call Clean_up()
contains
    subroutine Init()
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
    end subroutine Setup_Test_System

    subroutine Transit_time_test()
    end subroutine Transit_time_test
end program Verlet_Tests