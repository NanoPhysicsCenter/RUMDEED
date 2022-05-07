!-------------------------------------------!
! Module for Electron hole pair creation    !
! and removal                               !
! Kristinn Torfason                         !
! 08.05.15                                  !
!-------------------------------------------!

module mod_pair
  use mod_global
  implicit none

  ! ----------------------------------------------------------------------------
  ! Interface to combine the functions into one
#if defined(__INTEL_COMPILER)
  interface compact_array
    module procedure compact_array_2D_double_alt, compact_array_1D_double_alt, compact_array_1D_int_alt
  end interface
#else
interface compact_array
  module procedure compact_array_2D_double, compact_array_1D_double, compact_array_1D_int
end interface
#endif
contains
  ! ----------------------------------------------------------------------------
  ! Subroutine to add a particle to the system
  ! note: This subroutine should be called inside a OpenMP critical section
  ! Keyword arguments:
  ! par_pos: Particle postion (x, y, z)
  ! par_vel: Particle initial velocity (v_x, v_y, v_z)
  ! par_species: The type of particle, search for species_elec in mod_global to see a list
  ! step: The current time step, i.e. when the particle is emitted
  ! emit: The number of the emitter that the particle came from
  subroutine Add_Particle(par_pos, par_vel, par_species, step, emit, life, opt_sec)
    double precision, dimension(1:3), intent(in) :: par_pos,     par_vel
    integer, intent(in)                          :: par_species, step, emit, life
    integer, intent(in), optional                :: opt_sec
    integer                                      :: sec

    ! Check if we have reach the maximum number of paticles allowed
    if (nrPart+1 > MAX_PARTICLES) then
      print '(a)', 'Vacuum: WARNING MAX_PARTICLES REACHED. INCREASE MAX_PARTICLES'
      ! Add code to reallocate arrays to larger size?
    else

      if (present(opt_sec) .eqv. .true.) then
        if (opt_sec > MAX_SECTIONS) then
          sec = MAX_SECTIONS
          print '(a)', 'Vacuum: WARNING MAX_SECTIONS REACHED. INCREASE MAX_SECTIONS'
          print *, opt_sec
        else
          sec = opt_sec
        end if
      else
        ! To do: Have some default section rules, like 10x10 for square emitters?
        sec = 1
      end if

      ! Add the particle
      particles_cur_pos(:, nrPart+1)      = par_pos
      particles_prev_pos(:, nrPart+1)     = -1.0d0*length_scale
      particles_last_col_pos(:, nrPart+1) = par_pos
      particles_cur_accel(:, nrPart+1)    = 0.0d0
      particles_prev_accel(:, nrPart+1)   = 0.0d0
      particles_prev2_accel(:, nrPart+1)  = 0.0d0
      particles_cur_vel(:, nrPart+1)      = par_vel
      particles_step(nrPart+1)            = step
      particles_mask(nrPart+1)            = .true.
      particles_species(nrPart+1)         = par_species
      particles_emitter(nrPart+1)         = emit
      particles_section(nrPart+1)         = sec
      particles_life(nrPart+1)            = life
      particles_id(nrPart+1)              = nrID ! ID is updated near the end

      if (par_species == species_elec) then ! Electron
        particles_charge(nrPart+1) = -1.0d0*q_0
        particles_mass(nrPart+1) = m_eeff*m_0
        nrElec = nrElec + 1
        ! Write out the x and y position of the emitted particle
        ! along with which emitter and section it came from.
        !if (abs(par_pos(3) - 1.0d0*length_scale) < 1.0E-3) then
          write(unit=ud_density_emit) (par_pos(1) / length_scale), (par_pos(2) / length_scale), emit, sec, nrID
        !end if
      else ! N_2^+ Ion
        particles_charge(nrPart+1) = +1.0d0*q_0
        particles_mass(nrPart+1) = m_N2p
        nrHole = nrHole + 1
        write(unit=ud_density_ion) (par_pos(1) / length_scale), (par_pos(2) / length_scale), (par_pos(3) / length_scale), nrID
      end if

      ! Update the number of particles in the system
      nrElecHole = nrElec + nrHole
      nrPart = nrElecHole
      !endElecHoles = nrPart

      nrID = nrID + 1
    end if
  end subroutine Add_Particle

  ! ----------------------------------------------------------------------------
  ! Mark a particle for removal
  ! Keyword arguments:
  ! i -- The particle to be removed
  ! m -- Why the particle is being removed
  ! See mod_global for list of removal flags
  ! if particles_mask(i) is .true. then the particle is active
  ! if particles_mask(i) is .false. then the particle is inactive and should be removed
  subroutine Mark_Particles_Remove(i, m)
    integer, intent(in) :: i, m
    integer             :: emit

    ! Check if the particle has already been marked for removal
    ! if so just return
    if (particles_mask(i) .eqv. .false.) return

    ! Mark the particle for removal
    ! The particle is actually removed later when Remove_Particles is called.
    particles_mask(i) = .false.

    ! Set the charge to zero. That way it no longer has any effects on calculations.
    particles_charge(i) = 0.0d0

    ! Take care of the book keeping
    !$OMP ATOMIC UPDATE
    nrPart_remove = nrPart_remove + 1

    if (particles_species(i) == species_elec) then

      !$OMP ATOMIC UPDATE
      nrElec_remove = nrElec_remove + 1

      SELECT CASE (m)
      CASE (remove_top)
          !$OMP ATOMIC UPDATE
          nrPart_remove_top = nrPart_remove_top + 1

          !$OMP ATOMIC UPDATE
          nrElec_remove_top = nrElec_remove_top + 1

          emit = particles_emitter(i)

          !$OMP ATOMIC UPDATE
          nrElec_remove_top_emit(emit) = nrElec_remove_top_emit(emit) + 1

          ! Write out the transverse position of the particle, its velocity and which emitter/section it came from.
          !$OMP CRITICAL(DENSITY_ABSORB_TOP)
          write(unit=ud_density_absorb_top) particles_cur_pos(1, i)/length_scale, particles_cur_pos(2, i)/length_scale, &
                                          & particles_cur_vel(1, i), particles_cur_vel(2, i), particles_cur_vel(3, i), &
                                          & particles_emitter(i), particles_section(i), particles_id(i)
          !$OMP END CRITICAL(DENSITY_ABSORB_TOP)
        CASE (remove_bot)
          !$OMP ATOMIC UPDATE
          nrPart_remove_bot = nrPart_remove_bot + 1

          !$OMP ATOMIC UPDATE
          nrElec_remove_bot = nrElec_remove_bot + 1

          ! Write out the x and y position of the particle along with which emitter it came from.
          !$OMP CRITICAL(DENSITY_ABSORB_BOT)
          write(unit=ud_density_absorb_bot) particles_cur_pos(1, i)/length_scale, particles_cur_pos(2, i)/length_scale, &
                                          & particles_emitter(i), particles_section(i), particles_id(i)
          !$OMP END CRITICAL(DENSITY_ABSORB_BOT)
        CASE DEFAULT
          print *, 'Error unkown remove case ', m
      END SELECT

    else if (particles_species(i) == species_hole) then

      !$OMP ATOMIC UPDATE
      nrHole_remove = nrHole_remove + 1

      SELECT CASE (m)
      CASE (remove_top)
          !$OMP ATOMIC UPDATE
          nrPart_remove_top = nrPart_remove_top + 1

          !$OMP ATOMIC UPDATE
          nrHole_remove_top = nrHole_remove_top + 1
        CASE (remove_bot)
          !$OMP ATOMIC UPDATE
          nrPart_remove_bot = nrPart_remove_bot + 1

          !$OMP ATOMIC UPDATE
          nrHole_remove_bot = nrHole_remove_bot + 1
        CASE DEFAULT
          print *, 'Error unkown remove case ', m
      END SELECT
    else
      ! This should probably never happen
      print *, 'Error: Removing unknow particle'
      print *, 'particles_species(i) = ', particles_species(i)
    end if
  end subroutine Mark_Particles_Remove


  ! ----------------------------------------------------------------------------
  ! A subroutine to remove particles from the system.
  ! It takes the array's that hold information and about the particles and
  ! compacts them.
  !
  ! If particle 3 has been marked for removal then
  ! | 1 | 2 | - | 4 | 5 | -> | 1 | 2 | 4 | 5 | - |
  !
  ! Keyword arguments:
  ! step -- The current time step
  subroutine Remove_Particles(step)
    integer, intent(in) :: step
    integer             :: m, k


    call Write_Absorbed(step)

    if ((nrPart_remove > 0) .and. (nrPart > 0)) then ! Check if we have some thing to do

      if ((nrPart - nrPart_remove) > 0) then ! Check if we can skip this
        !$OMP PARALLEL DEFAULT(NONE) PRIVATE(m, k) &
        !$OMP& SHARED(particles_mask, particles_cur_pos, particles_prev_pos, particles_last_col_pos) &
        !$OMP& SHARED(particles_cur_vel, particles_cur_accel, particles_prev_accel, particles_prev2_accel) &
        !$OMP& SHARED(particles_step, particles_species, step, particles_mass, particles_charge) &
        !$OMP& SHARED(particles_emitter, particles_section, particles_life, particles_id, nrPart, life_time)
        !$OMP SINGLE

        k = 1
        m = nrPart

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_cur_pos, particles_mask)
        !call compact_array_2D_double_verbal(particles_cur_pos, particles_mask, k, m)
        call compact_array(particles_cur_pos, particles_mask, k, m)
        !$OMP END TASK

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_prev_pos, particles_mask)
        call compact_array(particles_prev_pos, particles_mask, k, m)
        !$OMP END TASK

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_last_col_pos, particles_mask)
        call compact_array(particles_last_col_pos, particles_mask, k, m)
        !$OMP END TASK

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_cur_vel, particles_mask)
        call compact_array(particles_cur_vel, particles_mask, k, m)
        !$OMP END TASK

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_cur_accel, particles_mask)
        call compact_array(particles_cur_accel, particles_mask, k, m)
        !$OMP END TASK

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_prev_accel, particles_mask)
        call compact_array(particles_prev_accel, particles_mask, k, m)
        !$OMP END TASK

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_prev2_accel, particles_mask)
        call compact_array(particles_prev2_accel, particles_mask, k, m)
        !$OMP END TASK

        !$OMP TASK FIRSTPRIVATE(k, m, step) SHARED(particles_species, particles_step, particles_mask, life_time)
        call record_lifetime(particles_step, particles_mask, k, m, step) ! This also compacts the array
        call compact_array(particles_species, particles_mask, k, m)
        !$OMP END TASK

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_mass, particles_mask)
        call compact_array(particles_mass, particles_mask, k, m)
        !$OMP END TASK

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_charge, particles_mask)
        call compact_array(particles_charge, particles_mask, k, m)
        !$OMP END TASK

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_emitter, particles_mask)
        call compact_array(particles_emitter, particles_mask, k, m)
        !$OMP END TASK

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_section, particles_mask)
        call compact_array(particles_section, particles_mask, k, m)
        !$OMP END TASK

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_life, particles_mask)
        call compact_array(particles_life, particles_mask, k, m)
        !$OMP END TASK

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_id, particles_mask)
        call compact_array(particles_id, particles_mask, k, m)
        !$OMP END TASK

        !$OMP END SINGLE NOWAIT

        ! Wait for all tasks to finish
        !$OMP TASKWAIT

        !$OMP END PARALLEL
      end if

      ! !Sanity check
      ! l = 0
      ! do j = 1, MAX_PARTICLES
      !   if (particles_mask(j) .eqv. .false.) then
      !     l = l + 1
      !     !print *, 'hi'
      !   end if
      ! end do
      !
      ! if (l /= nrPart_remove) then
      !   print *, 'l = ', l
      !   print *, 'nrPart_remove = ', nrPart_remove
      !   pause
      ! end if

      ! Update the number of particles in the system
      !nrPart = nrPart - nrPart_remove
      nrElec = nrElec - nrElec_remove
      nrHole = nrHole - nrHole_remove

      !print *, 'nrElec_remove = ', nrElec_remove
      !print *, 'nrHole_remove = ', nrHole_remove
      !print *, 'nrPart_remove = ', nrPart_remove
      !print *, 'nrPart - nrPart_remove = ', nrPart - nrPart_remove

      ! This should not happen, but just in case.
      if (nrElec < 0) nrElec = 0
      if (nrHole < 0) nrHole = 0

      nrElecHole = nrElec + nrHole
      nrPart = nrElecHole

      !startElecHoles = 1
      !endElecHoles = startElecHoles + nrElecHole - 1

      particles_mask = .true. ! Reset the mask. .true. means all particles are active.

      ! Reset the number of particles to remove
      nrPart_remove = 0
      nrElec_remove = 0
      nrHole_remove = 0

      nrPart_remove_top = 0
      nrPart_remove_bot = 0
      nrElec_remove_top = 0
      nrElec_remove_bot = 0
      nrHole_remove_top = 0
      nrHole_remove_bot = 0

      nrElec_remove_top_emit(1:MAX_EMITTERS) = 0
    end if

  end subroutine Remove_Particles

  ! --------------------------------------------------------------------------
  ! Write out the life time of particles
  subroutine Write_Life_Time()
    integer          :: IFAIL, i, ud_lifetime
    double precision :: time

    open(newunit=ud_lifetime, iostat=IFAIL, file='lifetime.dt', status='REPLACE', action='write')
    do i = 1, MAX_LIFE_TIME
      time = i * time_step / time_scale
      write(ud_lifetime, '(ES12.4, tr2, i6, tr2, i6, tr2, i6)', iostat=IFAIL) time, i, life_time(i, 1), life_time(i, 2)
    end do
    close(unit=ud_lifetime, iostat=IFAIL, status='keep')
  end subroutine Write_Life_Time

  ! ----------------------------------------------------------------------------
  ! Write data about number of particles absorbed. Called in Remove_Particle
  ! Keyword arguments:
  ! step -- The current time step
  subroutine Write_Absorbed(step)
    integer, intent(in) :: step
    integer             :: IFAIL, i

    write (ud_absorb, "(ES12.4, *(tr2, i8))", iostat=IFAIL) cur_time, step, &
    & nrPart_remove, nrElec_remove, nrHole_remove

    write (ud_absorb_top, "(ES12.4, *(tr2, i8))", iostat=IFAIL) cur_time, step, &
    & nrPart_remove_top, nrElec_remove_top, nrHole_remove_top, (nrElec_remove_top_emit(i), i = 1, nrEmit)

    write (ud_absorb_bot, "(ES12.4, *(tr2, i8))", iostat=IFAIL) cur_time, step, &
    & nrPart_remove_bot, nrElec_remove_bot, nrHole_remove_bot
  end subroutine Write_Absorbed

  ! ----------------------------------------------------------------------------
  ! Write out the current and information about the number of particles in
  ! the system
  ! Keyword arguments:
  ! step -- The current time step
  subroutine Write_Ramo_Current(step)
    integer, intent(in) :: step
    integer             :: IFAIL, i
    double precision    :: ramo_cur
    double precision    :: avg_speed

    avg_speed = sqrt(avg_vel(1)**2 + avg_vel(2)**2 + avg_vel(3)**2)

    ! Write out the current for each emitter and section
    if (write_ramo_sec .eqv. .true.) then
      !write(unit=ud_ramo_sec, asynchronous='YES') ramo_current_emit
      write(unit=ud_ramo_sec) ramo_current_emit
    end if

    ! Write the total current along with other data
    ramo_cur = sum(ramo_current) / cur_scale
    write(unit=ud_ramo, &
    & fmt="(ES12.4, tr2, i8, tr2, ES12.4, tr2, ES12.4, tr2, i6, tr2, i6, tr2, i6, tr2, ES12.4, tr2, ES12.4, *(tr2, ES12.4))", &
    & iostat=IFAIL) &
    & cur_time, step, ramo_cur, V_d, nrPart, nrElec, nrHole, avg_mob, avg_speed, (ramo_current(i), i = 1, nrSpecies)

  end subroutine Write_Ramo_current

  !-----------------------------------------------------------------------------
  ! Write the current position of all particles to a file
  subroutine Write_Position(step)
    integer, intent(in)              :: step
    integer                          :: i
    integer, parameter               :: N_steps = 10
    double precision, dimension(1:3) :: par_pos

    if (write_position_file .eqv. .True.) then
      if (step == 1) then
        ! Write at the start of the file the total number of time steps
        write(unit=ud_pos) steps, N_steps
      end if

      ! Write position out every N_steps
      if (mod(step, N_steps) == 0) then

        ! Write out what time step we are on and the current number of particles
        write(unit=ud_pos) step, nrPart

        do i = 1, nrPart
          par_pos(:) = particles_cur_pos(:, i) ! Position of the particle

          ! Write out x, y, z and which emitter the particle came from
          write(unit=ud_pos) par_pos(1), par_pos(2), par_pos(3), particles_emitter(i), particles_section(i), particles_id(i)
        end do
      end if
    end if
  end subroutine Write_Position

  !-----------------------------------------------------------------------------
  ! Write the current position of all particles to a files.
  ! The files will be in XYZ format an each timestep is written to a seperate file.
  ! The file name format is position.xyz.? where ? is the timestep.
  ! Datastructure in each file is x, y, z, species, vel_x, vel_y, vel_z, speed
  subroutine Write_Position_XYZ_Step(step)
    integer, intent(in)              :: step
    integer                          :: i, par_species
    integer                          :: IFAIL, ud_pos_data
    integer, parameter               :: N_steps = 100
    double precision, dimension(1:3) :: par_pos, par_vel
    double precision                 :: par_speed
    character(len=128)               :: filename

    ! Check if we are writing a position file
    if (write_position_file .eqv. .True.) then

      ! Write position out every N_steps
      if (mod(step, N_steps) == 0) then
        ! Generate the file name posotion.xyz.? and open the file for writing
        write(filename, '(a12, i0)') 'out/pos.xyz.', step
        open(newunit=ud_pos_data, iostat=IFAIL, file=filename, status='REPLACE', action='write')
        if (IFAIL /= 0) then
          print *, 'Vacuum: Failed to open the file to write position data'
          return
        end if

        ! Loop over all particles
        do i = 1, nrPart
          ! Get particle position and scale to nm
          par_pos = particles_cur_pos(:, i) / length_scale
          ! Ger particle velocity and calulate speed
          par_vel = particles_cur_vel(:, i)
          par_speed = sqrt(par_vel(1)**2 + par_vel(2)**2 + par_vel(3)**2)
          ! Get particle species, i.e. electron, ion, hole, ...
          par_species = particles_species(i)

          write(unit=ud_pos_data, fmt="(F6.2, F6.2, F6.2, i2, ES12.4, ES12.4, ES12.4, ES12.4)", iostat=IFAIL) &
            par_pos(1), par_pos(2), par_pos(3), par_species, par_vel(1), par_vel(2), par_vel(3), par_speed
        end do

        ! Close the file
        close(unit=ud_pos_data, iostat=IFAIL, status='keep')
      end if
    end if
  end subroutine

  !-----------------------------------------------------------------------------
  ! Write out the positions
  ! Keyword arguments:
  ! step -- The current time step
  subroutine Write_Particle_Data(step)
    integer, intent(in) :: step
    integer             :: i, IFAIL, ud_data
    character(len=128)  :: filename

    ! Prepare the name of the output file
    ! each file is named particles-0.dt where the number
    ! represents the current time step.
    write(filename, '(a14, i0, a3)') 'out/particles-', step, '.dt'

    ! Open the output file
    open(newunit=ud_data, iostat=IFAIL, file=filename, status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open the data file.'
      return
    end if

    ! The first line in the file is the number of particles
    write(unit=ud_data, fmt="(i8)", iostat=IFAIL) nrPart

    ! All the other lines are data about the particles, with each particle on its
    ! own line.
    ! Position, Velocity, Acceleration
    do i = 1, nrPart
      write(unit=ud_data, fmt="(*(ES12.4))", iostat=IFAIL) &
        particles_cur_pos(:, i), particles_cur_vel(:, i), particles_cur_accel(:, i)
    end do

    ! Close the file
    close(unit=ud_data, iostat=IFAIL, status='keep')
  end subroutine

  ! ----------------------------------------------------------------------------
  ! Record the life time of particles.
  ! The subroutine should be called before the particles_species array
  ! is compacted. This subroutine will also compact the particles_step array.
  subroutine record_lifetime(A, mask, k, m, step)
    integer, dimension(:), intent(inout) :: A
    logical, dimension(:), intent(in)    :: mask
    integer, intent(in)                  :: m, k, step
    integer                              :: i, j, lt, s

    j = 0
    !k = lbound(A, 1)
    !k = startElecHoles
    do i = k, m
      if (mask(i) .eqv. .true.) then
        j = j + 1
        A(j) = A(i)
      else
        lt = step - A(i)
        if (lt <= 0) lt = 1
        if (lt > MAX_LIFE_TIME) lt = MAX_LIFE_TIME

        s = particles_species(i)

        life_time(lt, s) = life_time(lt, s) + 1

      end if
    end do
  end subroutine


  ! ----------------------------------------------------------------------------
  ! A subroutine to compactify the 2D arrays used to store data
  ! The mask says which particles we should keep and which should be removed.
  ! It is .true. for particles that should be keept and .false. for particles
  ! that it should be removed. The loop jumps over particles that are marked
  ! as .false..
  !
  ! Note: The pack subroutine also does this, but returns 1D arrays!!
  !

  ! Todo: Check if this has better performance than the pack method.
  subroutine compact_array_2D_double_alt(A, mask, k, m)
    double precision, dimension(:, :), intent(inout) :: A
    logical, dimension(:), intent(in)                :: mask
    integer, intent(in)                              :: m, k ! m = nrPart
    integer                                          :: i, j

    j = 0
    do i = k, m
      if (mask(i) .eqv. .true.) then ! Check if mask(i) == .true.
        j = j + 1
        A(:, j) = A(:, i)
      !else
      !  print *, 'Removed particle ', i, j, ' at ', particles_cur_pos(:, i) / length_scale, particles_cur_pos(:, j) / length_scale
      end if
    end do
    !print *, ''
  end subroutine compact_array_2D_double_alt

  ! This subroutine uses the pack method.
  subroutine compact_array_2D_double(A, mask, k, m)
    double precision, dimension(:, :), intent(inout) :: A
    logical, dimension(:), intent(in)                :: mask
    !logical, allocatable, dimension(:, :)            :: mask_2d
    integer, intent(in)                              :: m, k ! m = nrPart

    A(1, :) = pack(A(1, :), mask, A(1, :))
    A(2, :) = pack(A(2, :), mask, A(2, :))
    A(3, :) = pack(A(3, :), mask, A(3, :))

    ! We could also do! Performance?
    ! A = reshape(pack(A, mask), (/ 3, N/))
  end subroutine compact_array_2D_double


  ! ----------------------------------------------------------------------------
  ! A subroutine to compactify the 1D arrays used to store data
  subroutine compact_array_1D_double(A, mask, k, m)
    double precision, dimension(:), intent(inout) :: A
    logical, dimension(:), intent(in)             :: mask
    integer, intent(in)                           :: m, k

    A = pack(A, mask, A)
  end subroutine compact_array_1D_double

  subroutine compact_array_1D_double_alt(A, mask, k, m)
    double precision, dimension(:), intent(inout) :: A
    logical, dimension(:), intent(in)             :: mask
    integer, intent(in)                           :: m, k ! m = nrPart
    integer                                       :: i, j

    j = 0
    do i = k, m
      if (mask(i) .eqv. .true.) then ! Check if mask(i) == .true.
        j = j + 1
        A(j) = A(i)
      !else
      !  print *, 'Removed particle ', i, j, ' at ', particles_cur_pos(:, i) / length_scale, particles_cur_pos(:, j) / length_scale
      end if
    end do
    !print *, ''
  end subroutine compact_array_1D_double_alt

  ! ----------------------------------------------------------------------------
  ! A subroutine to compactify the 1D arrays used to store data
  subroutine compact_array_1D_int(A, mask, k, m)
    integer, dimension(:), intent(inout) :: A
    logical, dimension(:), intent(in)    :: mask
    integer, intent(in)                  :: m, k

    A = pack(A, mask, A)
  end subroutine compact_array_1D_int

  subroutine compact_array_1D_int_alt(A, mask, k, m)
    integer, dimension(:), intent(inout) :: A
    logical, dimension(:), intent(in)    :: mask
    integer, intent(in)                  :: m, k ! m = nrPart
    integer                              :: i, j

    j = 0
    do i = k, m
      if (mask(i) .eqv. .true.) then ! Check if mask(i) == .true.
        j = j + 1
        A(j) = A(i)
      !else
      !  print *, 'Removed particle ', i, j, ' at ', particles_cur_pos(:, i) / length_scale, particles_cur_pos(:, j) / length_scale
      end if
    end do
    !print *, ''
  end subroutine compact_array_1D_int_alt
end module mod_pair
