!-------------------------------------------!
! Module for Verlet Integration             !
! Kristinn Torfason                         !
! 21.01.16                                  !
!-------------------------------------------!

module mod_verlet
  use mod_global
  use mod_pair
  implicit none
contains

  ! ----------------------------------------------------------------------------
  ! Call the method used to update the positions
  subroutine Update_Position(step)
    integer, intent(in) :: step

    call Velocity_Verlet(step)
  end subroutine Update_Position


  ! ----------------------------------------------------------------------------
  ! Velocity Verlet
  subroutine Velocity_Verlet(step)
    integer, intent(in) :: step

    !$OMP SINGLE
    ramo_current = 0.0d0

    ! Update the current time
    cur_time = time_step * at_step / time_scale ! Scaled in units of time_scale
    at_step = at_step + 1

    ! Update the voltage in the system
    call Set_Voltage(step)

    ! Add electron/holes to the system
    !call Pair_Creation(step, i)

    !$OMP END SINGLE

    ! Update the position of particles
    if (nrElecHole > 0) then
      call Update_Fast_Position(step)

      call Calculate_Acceleration_Particles()

      call Update_Fast_Velocity(step)
    end if
    !call Update_Dipoles(step, i)

    !!!$OMP SINGLE

    ! Write out the current in the system
    call Write_Ramo_Current(step, at_step)

    ! Remove particles from the system
    call Remove_Particles(step)
    !!!$OMP END SINGLE

    ! Write out the field in the system
    !call Write_Field_Longitudinal(step)
    !call Write_Density_Map(step)

  end subroutine Velocity_Verlet

  ! ----------------------------------------------------------------------------
  !
  subroutine Update_Fast_Position(step)
    integer, intent(in) :: step
    integer             :: i

    !$OMP DO PRIVATE(i) SCHEDULE(STATIC)
    do i = startElecHoles, endElecHoles
      particles_prev_pos(:, i)   = particles_cur_pos(:, i) ! Store the previous position
      particles_cur_pos(:, i)    = particles_cur_pos(:, i) + particles_cur_vel(:, i)*time_step &
                               & + 0.5d0*particles_cur_accel(:, i)*time_step2

      particles_prev_accel(:, i) = particles_cur_accel(:, i)
      particles_cur_accel(:, i)  = 0.0d0

      !if (isnan(particles_cur_pos(1, i)) .eqv. .true.) then
      !  print *, 'Error'
      !  !pause
      !end if


      ! Mark particles that should be removed with .false. in the mask array
      call Check_Boundary_ElecHole(i)
    end do
    !$OMP END DO
  end subroutine Update_Fast_Position

  ! ----------------------------------------------------------------------------
  ! Checks the boundary conditions of the box. Enforce periodic boundary conditions
  ! and check which particles to remove.
  subroutine Check_Boundary_ElecHole(i)
    integer, intent(in) :: i
    double precision    :: x, y, z


    x = particles_cur_pos(1, i)
    y = particles_cur_pos(2, i)
    z = particles_cur_pos(3, i)

    if (x < 0.0d0) then
      !particles_cur_pos(1, i) = z + d
      particles_cur_pos(1, i) = box_dim(1) - x
    else if (x > box_dim(1)) then
      !particles_cur_pos(1, i) = z - d
      particles_cur_pos(1, i) = x - box_dim(1)
    end if

    if (z < 0.0d0) then
      !particles_cur_pos(3, i) = z + d
      particles_cur_pos(3, i) = box_dim(3) - z
    else if (z > box_dim(3)) then
      !particles_cur_pos(3, i) = z - d
      particles_cur_pos(3, i) = z - box_dim(3)
    end if

    ! Check if the particle should be removed from the system
    if (y < 0.0d0) then
      ! Only remove holes on the left side
      !if (particles_species(i) == species_hole) then
        call Mark_Particles_Remove(i, remove_bot)
      !else
      !  particles_cur_pos(2, i) = -1.0d0*y
      !  particles_cur_vel(2, i) = -1.0d0*particles_cur_vel(2, i)
      !end if
    else if (z > box_dim(3)) then
      ! Only remove electrons on the right side
      !if (particles_species(i) == species_elec) then
        call Mark_Particles_Remove(i, remove_top)
      !else
      !  particles_cur_pos(2, i) = 2.0d0*d - y
      !  particles_cur_vel(2, i) = -1.0d0*particles_cur_vel(2, i)
      !end if
    end if

    !nrPart_remove = nrElec_remove + nrHole_remove
    !nrPart_remove_left = nrElec_remove_left + nrHole_remove_left
    !nrPart_remove_right = nrElec_remove_right  + nrHole_remove_right
  end subroutine Check_Boundary_ElecHole


  ! ----------------------------------------------------------------------------
  !
  subroutine Update_Fast_Velocity(step)
    integer, intent(in) :: step
    integer             :: i, k
    double precision    :: q

    !$OMP DO PRIVATE(i, q, k)
    do i = startElecHoles, endElecHoles
      particles_cur_vel(:, i) = particles_cur_vel(:, i) &
                            & + 0.5d0*( particles_prev_accel(:, i) &
                            & + particles_cur_accel(:, i) )*time_step

      ! if (isnan(particles_cur_vel(2, i)) == .true.) then
      !   print *, 'i = ', i
      !   print *, 'particles_prev_accel(2, i) = ', particles_prev_accel(2, i)
      !   print *, 'particles_cur_accel(2, i) = ', particles_cur_accel(2, i)
      !   pause
      ! end if

      ! if (isnan(particles_cur_vel(1, i)) .eqv. .true.) then
      !   print *, 'Error'
      ! end if

      q = particles_charge(i)
      k = particles_species(i)

      !$OMP ATOMIC
      ramo_current(k) = ramo_current(k) + q * E_zunit * particles_cur_vel(2, i)
    end do
    !$OMP END DO
  end subroutine Update_Fast_Velocity


  ! ----------------------------------------------------------------------------
  ! Acceleration
  subroutine Calculate_Acceleration_Particles()
    double precision, dimension(1:3) :: force_E, force_c
    double precision, dimension(1:3) :: pos_1, pos_2, diff
    double precision                 :: r
    double precision                 :: q_1, q_2
    double precision                 :: im_1, im_2
    double precision                 :: pre_fac_c
    integer                          :: i, j, k_1, k_2

    !$OMP DO PRIVATE(i, j, k_1, k_2, pos_1, pos_2, diff, r, force_E, force_c, im_1, q_1, im_2, q_2, pre_fac_c) &
    !$OMP& REDUCTION(+:particles_cur_accel) SCHEDULE(GUIDED)
    do i = 1, nrPart
      ! Information about the particle we are calculating the force/acceleration on
      pos_1 = particles_cur_pos(:, i)
      im_1 = 1.0d0 / particles_mass(i)
      q_1 = particles_charge(i)
      k_1 = particles_species(i)

      ! Acceleration due to electric field
      force_E = q_1 * field_E(pos_1)

      ! Loop over particles from i+1 to nrElec.
      ! There is no need to loop over all particles since
      ! The forces are equal but in opposite directions
      do j = i+1, nrPart

        ! Information about the particle that is acting on the particle at pos_1
        pos_2 = particles_cur_pos(:, j)
        im_2 = 1.0d0 / particles_mass(j)
        q_2 = particles_charge(j)
        k_2 = particles_species(j)


        pre_fac_c = q_1*q_2 * div_fac_c ! q_1*q_2 / (4*pi*epsilon)

        ! Calculate the distance between the two particles
        diff = pos_1 - pos_2
        r = NORM2(diff) + length_scale**3 ! Prevent singularity

        ! Calculate the Coulomb force
        ! F = (r_1 - r_2) / |r_1 - r_2|^3
        ! F = (diff / r) * 1/r^2
        ! (diff / r) is a unit vector
        force_c = diff * r**(-3)

        ! if (isnan(force_c(2)) .eqv. .true.) then
        !   print *, 'force_c(2) = ', force_c(2)
        !   print *, 'r = ', r / length_scale
        !   print *, 'diff = ', diff / length_scale
        !   print *, 'i = ', i
        !   print *, 'j = ', j
        !   print *, 'pos_1 = ', pos_1 / length_scale
        !   print *, 'pos_2 = ', pos_2 / length_scale
        !   print *, 'particles_mask(i) = ', particles_mask(i)
        !   print *, 'particles_mask(j) = ', particles_mask(j)
        !   print *, 'q_1 = ', q_1
        !   print *, 'q_2 = ', q_2
        !   print *, 'pos_1(2) == pos_2(2)', (pos_1(2) == pos_2(2))
        !   print *, 'particles_species(i) = ', particles_species(i)
        !   print *, 'particles_species(j) = ', particles_species(j)
        !   print *, 'nrElec = ', nrElec
        !   print *, 'nrHole = ', nrHole
        !   print *, 'nrElecHole = ', nrElecHole
        !   pause
        ! else
          ! Record the acceleration using the fact that they are equal but
          ! in opposite directions.
          particles_cur_accel(:, j) = particles_cur_accel(:, j) - pre_fac_c * force_c * im_2
          particles_cur_accel(:, i) = particles_cur_accel(:, i) + pre_fac_c * force_c * im_1
        ! end if
      end do

      particles_cur_accel(:, i) = particles_cur_accel(:, i) + force_E * im_1
    end do
    !$OMP END DO
  end subroutine Calculate_Acceleration_Particles


  !-----------------------------------------------------------------------------
  ! Find field in point
  function Calc_Field_at(pos)
    double precision, dimension(1:3)             :: Calc_Field_at
    double precision, dimension(1:3), intent(in) :: pos

    double precision, dimension(1:3) :: force_c
    double precision, dimension(1:3) :: pos_1, pos_2, diff
    double precision                 :: r
    double precision                 :: q_2
    double precision                 :: pre_fac_c
    integer                          :: j


    !$OMP SINGLE
    ! Position of the particle we are calculating the force/acceleration on
    pos_1 = pos

    ! Acceleration due to electric field
    force_tot = field_E(pos_1)
    !$OMP END SINGLE

    !$OMP DO PRIVATE(j, pos_2, diff, r, force_c, q_2, pre_fac_c) &
    !$OMP& REDUCTION(+:force_tot) SCHEDULE(GUIDED)
    do j = 1, nrPart

      ! Position of the particle that is acting on the particle at pos_1
      pos_2 = particles_cur_pos(:, j)
      q_2 = particles_charge(j)

      pre_fac_c = q_2 * div_fac_c ! q_1*q_2 / (4*pi*epsilon)

      ! Calculate the distance between the two particles
      diff = pos_1 - pos_2
      r = NORM2(diff) + length_scale**3 ! distance + Prevent singularity

      ! Calculate the Coulomb force
      ! F = (r_1 - r_2) / |r_1 - r_2|^3
      ! F = (diff / r) * 1/r^2
      ! (diff / r) is a unit vector
      force_c = diff * r**(-3)

      if (isnan(force_c(2)) .eqv. .false.) then
        force_tot = force_tot + pre_fac_c * force_c
      end if
    end do
    !$OMP END DO

    Calc_Field_at = force_tot

  end function Calc_Field_at


  ! ----------------------------------------------------------------------------
  ! The vacuum electric field
  pure function field_E(pos)
    double precision, dimension(1:3), intent(in) :: pos
    double precision, dimension(1:3)             :: field_E

    ! Electric field
    ! x
    field_E(1) = 0.0d0

    ! y
    field_E(2) = 0.0d0

    ! z
    field_E(3) = E_z

  end function field_E


  !-----------------------------------------------------------------------------
  ! The voltage
  subroutine Set_Voltage(step)
    integer, intent(in) :: step
    integer             :: IFAIL

    V = V_a
    E_z = -1.0d0*V/d
    !E_zunit = -1.0d0*sign(1.0d0, V)/d

    write (ud_volt, "(ES12.4, tr2, i8, tr2, i8, tr2, ES12.4)", iostat=IFAIL) cur_time, step, at_step, V
  end subroutine Set_Voltage

end module mod_verlet
