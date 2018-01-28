!-------------------------------------------!
! Module for Electron hole pair creation    !
! and removal                               !
! Kristinn Torfason                         !
! 08.05.15                                  !
!-------------------------------------------!

module mod_pair
  use mod_global
  implicit none

  ! Defines how the energy is split up between the electron and hole
  ! This comes from the convervation of momentum, m_e*v_ei = m_h*v_hi,
  ! where v_ei = sqrt(2*k*\Delta E/m_e) and v_hi = sqrt(2*(1-k)*\Delta E/m_h)
  double precision, parameter :: k_p = m_heff / (m_eeff + m_heff)
  double precision, parameter :: E_gj = E_g*q_0
  ! ----------------------------------------------------------------------------
  ! Interface to combine the two functions into one
  interface compact_array
    module procedure compact_array_2D_double, compact_array_1D_double, compact_array_1D_int
  end interface
contains
  subroutine Add_Particle(par_pos, par_vel, par_species, step)
    double precision, dimension(1:3), intent(in) :: par_pos, par_vel
    integer, intent(in)                          :: par_species, step

    ! Add the electron
    particles_cur_pos(:, nrPart+1) = par_pos
    particles_prev_pos(:, nrPart+1) = -1.0d0*length_scale
    particles_cur_accel(:, nrPart+1) = 0.0d0
    particles_prev_accel(:, nrPart+1) = 0.0d0
    particles_cur_vel(:, nrPart+1) = par_vel
    particles_step(nrPart+1) = step
    particles_mask(nrPart+1) = .true.
    particles_species(nrPart+1) = par_species

    if (par_species == species_elec) then
      particles_charge(nrPart+1) = -1.0d0*q_0
      particles_mass(nrPart+1) = m_eeff*m_0
      nrElec = nrElec + 1
    else
      particles_charge(nrPart+1) = +1.0d0*q_0
      particles_mass(nrPart+1) = m_heff*m_0
      nrHole = nrHole + 1
    end if

    if (particles_mass(nrPart+1) == 0.0d0) then
      print *, 'WTF'
      pause
    end if

    ! Update the number of particles in the system
    nrElecHole = nrElec + nrHole
    nrPart = nrElecHole
    endElecHoles = nrPart
  end subroutine Add_Particle

  ! ----------------------------------------------------------------------------
  ! Mark a particle for removal
  ! Keyword arguments:
  ! i -- The particle to be removed
  ! m -- Why the particle is being removed
  !    m = 0: Unknown (remove_unknown)
  !    m = 1: Exited left boundary (remove_left)
  !    m = 2: Exited right boundary (remove_right)
  !    m = 3: Recombination (remove_recomb)
  subroutine Mark_Particles_Remove(i, m)
    integer, intent(in) :: i, m

    ! Check if the particle has already been marked for removal
    ! if so just return
    if (particles_mask(i) .eqv. .false.) return

    ! Mark the particle for removal
    particles_mask(i) = .false.

    ! Set the charge to zero. That way it no longer has any effects on calculations.
    particles_charge(i) = 0.0d0

    ! Take care of the book keeping
    !$OMP ATOMIC
    nrPart_remove = nrPart_remove + 1

    if (particles_species(i) == species_elec) then

      !$OMP ATOMIC
      nrElec_remove = nrElec_remove + 1

      SELECT CASE (m)
      CASE (remove_top)
          !$OMP ATOMIC
          nrPart_remove_top = nrPart_remove_top + 1

          !$OMP ATOMIC
          nrElec_remove_top = nrElec_remove_top + 1
        CASE (remove_bot)
          !$OMP ATOMIC
          nrPart_remove_bot = nrPart_remove_bot + 1

          !$OMP ATOMIC
          nrElec_remove_bot = nrElec_remove_bot + 1
        CASE DEFAULT
          print *, 'Error unkown remove case ', m
      END SELECT

    else if (particles_species(i) == species_hole) then

      !$OMP ATOMIC
      nrHole_remove = nrHole_remove + 1

      SELECT CASE (m)
      CASE (remove_top)
          !$OMP ATOMIC
          nrPart_remove_top = nrPart_remove_top + 1

          !$OMP ATOMIC
          nrHole_remove_top = nrHole_remove_top + 1
        CASE (remove_bot)
          !$OMP ATOMIC
          nrPart_remove_bot = nrPart_remove_bot + 1

          !$OMP ATOMIC
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

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m, k)

    !$OMP SINGLE

    call Write_Absorbed(step)

    if ((nrPart_remove > 0) .and. (nrPart > 0)) then ! Check if we have some thing to do

      if ((nrPart - nrPart_remove) > 0) then ! Check if we can skip this
        k = startElecHoles
        m = endElecHoles

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_cur_pos, particles_mask)
        !call compact_array_2D_double_verbal(particles_cur_pos, particles_mask, k, m)
        call compact_array(particles_cur_pos, particles_mask, k, m)
        !$OMP END TASK

        !$OMP TASK FIRSTPRIVATE(k, m) SHARED(particles_prev_pos, particles_mask)
        call compact_array(particles_prev_pos, particles_mask, k, m)
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

        !$OMP TASKWAIT
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

      if (nrElec < 0) nrElec = 0
      if (nrHole < 0) nrHole = 0

      nrElecHole = nrElec + nrHole
      nrPart = nrElecHole

      startElecHoles = 1
      endElecHoles = startElecHoles + nrElecHole - 1

      particles_mask = .true. ! Reset the mask

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
    end if

    !$OMP END SINGLE

    !$OMP END PARALLEL
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
    integer             :: IFAIL

    write (ud_absorb, "(ES12.4, *(tr2, i8))", iostat=IFAIL) cur_time, step, &
    & nrPart_remove, nrElec_remove, nrHole_remove

    write (ud_absorb_top, "(ES12.4, *(tr2, i8))", iostat=IFAIL) cur_time, step, &
    & nrPart_remove_top, nrElec_remove_top, nrHole_remove_top

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
    integer             :: i, IFAIL
    double precision    :: ramo_cur = 0.0d0

    !ramo_cur = 0.0d0

    !do i = 1, nrSpecies
    !  ramo_cur = ramo_cur + ramo_current(i)
    !end do
    ramo_cur = sum(ramo_current) / cur_scale

    write (ud_ramo, fmt="(ES12.4, tr2, i8, tr2, E12.4, tr2, E12.4, tr2, i6, tr2, i6, tr2, i6)", iostat=IFAIL) &
    & cur_time, step, ramo_cur, V_d, nrPart, nrElec, nrHole

  end subroutine Write_Ramo_current

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
  ! The mask says which particles we should keep and wich should be removed.
  ! It is .true. for particles that should be keept and .false. for particles
  ! that should be removed. The loop jumps over particles that are marked
  ! as .false..
  !
  ! Note: The pack subroutine also does this, but returns 1D arrays!!
  !
  subroutine compact_array_2D_double(A, mask, k, m)
    double precision, dimension(:, :), intent(inout) :: A
    logical, dimension(:), intent(in)                :: mask
    !logical, allocatable, dimension(:, :)            :: mask_2d
    integer, intent(in)                              :: m, k ! m = nrPart

    A(1, :) = pack(A(1, :), mask)
    A(2, :) = pack(A(2, :), mask)
    A(3, :) = pack(A(3, :), mask)
  end subroutine compact_array_2D_double


  subroutine compact_array_2D_double_verbal(A, mask, k, m)
    double precision, dimension(:, :), intent(inout) :: A
    logical, dimension(:), intent(in)                :: mask
    integer, intent(in)                              :: m, k ! m = nrPart
    integer                                          :: i, j

    j = 0
    do i = k, m
      if (mask(i) .eqv. .true.) then ! Check if mask(i) == .true.
        j = j + 1
        A(:, j) = A(:, i)
      else
        print *, 'Removed particle ', i, j, ' at ', particles_cur_pos(:, i) / length_scale, particles_cur_pos(:, j) / length_scale
      end if
    end do
    !print *, ''
  end subroutine compact_array_2D_double_verbal


  ! ----------------------------------------------------------------------------
  ! A subroutine to compactify the 1D arrays used to store data
  subroutine compact_array_1D_double(A, mask, k, m)
    double precision, dimension(:), intent(inout) :: A
    logical, dimension(:), intent(in)             :: mask
    integer, intent(in)                           :: m, k

    A = pack(A, mask)
  end subroutine compact_array_1D_double

  ! ----------------------------------------------------------------------------
  ! A subroutine to compactify the 1D arrays used to store data
  subroutine compact_array_1D_int(A, mask, k, m)
    integer, dimension(:), intent(inout) :: A
    logical, dimension(:), intent(in)    :: mask
    integer, intent(in)                  :: m, k

    A = pack(A, mask)
  end subroutine compact_array_1D_int
end module mod_pair
