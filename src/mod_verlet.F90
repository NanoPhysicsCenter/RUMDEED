!-------------------------------------------!
! Module for Verlet Integration             !
! Kristinn Torfason                         !
! 21.01.16                                  !
!-------------------------------------------!

Module mod_verlet
  use mod_global
  use mod_pair
  use mod_collisions
  implicit none

  ! RLC
  double precision :: V_rf = 0.0d0, V_rf_prev = 0.0d0, V_rf_next = 0.0d0
  double precision :: I_cur = 0.0d0, I_prev = 0.0d0

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

    ! Update the current time
    cur_time = time_step * step / time_scale ! Scaled in units of time_scale

    ! Update the voltage in the system
    call Set_Voltage(step)

    ! Update the position of particles (Electrons / Holes)
    !if (nrElecHole > 0) then
      call Update_ElecHole_Position(step)

      if (EMISSION_MODE == EMISSION_UNIT_TEST) then
        call Write_Position_Test(step)
      end if

      call Calculate_Acceleration_Particles()

      if ((EMISSION_MODE == EMISSION_UNIT_TEST) .or. (EMISSION_MODE == EMISSION_MANUAL)) then
        call Write_Acceleration_Test(step)
      end if

      ! Reset the ramo current. Note, this should be done after set_voltage,
      ! since the ramo current may be used there.
      wait(ud_ramo_sec) ! ud_ramo_sec is done asynchronously. We must make sure it is finished.
      ramo_current = 0.0d0
      ramo_current_emit = 0.0d0
      call Update_Velocity(step)
    !end if

    ! Write out the current in the system
    call Write_Ramo_Current(step)

  end subroutine Velocity_Verlet

  subroutine Do_Collisions(step)
    integer, intent(in) :: step

    if (collisions .eqv. .true.) then
      call Do_Ion_Collisions(step)
    end if
  end subroutine Do_Collisions

  subroutine Read_Cross_Section_Data()
    if (collisions .eqv. .true.) then
      print '(a)', 'Vacuum: Doing collisions reading in data'
      call Read_Cross_Section()
    end if
   end subroutine Read_Cross_Section_Data

  ! ----------------------------------------------------------------------------
  ! Update the position of particles in the verlet integration
  subroutine Update_ElecHole_Position(step)
    integer, intent(in) :: step
    integer             :: i

    !$OMP PARALLEL DO PRIVATE(i) &
    !$OMP& SHARED(particles_prev_pos, particles_cur_pos, particles_cur_vel, particles_cur_accel) &
    !$OMP& SHARED(time_step, time_step2, particles_prev2_accel, particles_prev_accel) &
    !$OMP& SCHEDULE(GUIDED, CHUNK_SIZE)
    do i = 1, nrPart
      ! Verlet
      !particles_prev_pos(:, i) = particles_cur_pos(:, i) ! Store the previous position
      !particles_cur_pos(:, i)  = particles_cur_pos(:, i) + particles_cur_vel(:, i)*time_step &
      !                       & + 0.5d0*particles_cur_accel(:, i)*time_step2

      ! Beeman
      particles_prev_pos(:, i) = particles_cur_pos(:, i) ! Store the previous position
      particles_cur_pos(:, i)  = particles_cur_pos(:, i) + particles_cur_vel(:, i)*time_step &
                             & + 1.0d0/6.0d0*( 4.0d0*particles_cur_accel(:, i) - particles_prev_accel(:, i) )*time_step2

      particles_prev2_accel(:, i) = particles_prev_accel(:, i) ! Beeman
      particles_prev_accel(:, i) = particles_cur_accel(:, i)
      particles_cur_accel(:, i)  = 0.0d0

      ! Mark particles that should be removed with .false. in the mask array
      call ptr_Check_Boundary(i)
    end do
    !$OMP END PARALLEL DO
  end subroutine Update_ElecHole_Position


  subroutine Write_Position_Test(step)
    integer, intent(in)              :: step
    integer                          :: ud_pos_test, IFAIL, i
    character(len=1024)              :: filename
    double precision, dimension(1:3) :: par_pos

    ! Prepare the name of the output file
    ! each file is named accel-0.dt where the number
    ! represents the current time step.
    write(filename, '(a8, i0, a4)') 'pos/pos-', step, '.bin'

    open(newunit=ud_pos_test, iostat=IFAIL, file=filename, status='replace', action='write', access='STREAM')
    if (IFAIL /= 0) then
      print '(a)', 'Vacuum: ERROR UNABLE TO OPEN file for position'
      print *, filename
      print *, step
    end if

    do i = 1, nrPart
      par_pos(:) = particles_cur_pos(:, i) ! Position of the particle

      ! Write out x, y, z
      write(unit=ud_pos_test) par_pos(1), par_pos(2), par_pos(3)
    end do

    close(unit=ud_pos_test, iostat=IFAIL, status='keep')
  end subroutine Write_Position_Test

  ! ----------------------------------------------------------------------------
  ! Checks the boundary conditions of the box.
  ! Check which particles to remove
  subroutine Check_Boundary_ElecHole_Planar(i)
    integer, intent(in) :: i
    double precision    :: z

    z = particles_cur_pos(3, i)

    ! Check if the particle should be removed from the system
    if (z < 0.0d0) then
        call Mark_Particles_Remove(i, remove_bot)
    else if (z > box_dim(3)) then
        call Mark_Particles_Remove(i, remove_top)
    end if

  end subroutine Check_Boundary_ElecHole_Planar


    ! ----------------------------------------------------------------------------
  ! Checks the boundary conditions of the box.
  ! Check which particles to remove
  ! Enforce periodic boundary conditions
  subroutine Check_Boundary_ElecHole_Periodic(i)
    integer, intent(in) :: i
    double precision    :: x, y, z

    x = particles_cur_pos(1, i)
    y = particles_cur_pos(2, i)
    z = particles_cur_pos(3, i)

    ! Check if the particle should be removed from the system
    if (z < 0.0d0) then
        call Mark_Particles_Remove(i, remove_bot)
    else if (z > box_dim(3)) then
        call Mark_Particles_Remove(i, remove_top)
    end if

    ! Do periodic boundary conditions
    particles_cur_pos(1, i) = modulo(x, box_dim(1))
    particles_cur_pos(1, i) = modulo(y, box_dim(2))

  end subroutine Check_Boundary_ElecHole_Periodic

  ! ----------------------------------------------------------------------------
  ! Update the velocity in the verlet integration
  subroutine Update_Velocity(step)
    integer, intent(in)              :: step
    integer                          :: i, k, emit, sec
    double precision                 :: q

    avg_vel(:) = 0.0d0

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, q, k, emit, sec) &
    !$OMP& SHARED(nrPart, particles_cur_vel, particles_prev_accel, particles_cur_accel, time_step) &
    !$OMP& SHARED(ramo_current, ramo_current_emit, E_zunit, particles_section, particles_charge) &
    !$OMP& SHARED(particles_species, particles_emitter, particles_prev2_accel) &
    !$OMP& REDUCTION(+:avg_vel) SCHEDULE(GUIDED, CHUNK_SIZE)
    do i = 1, nrPart
      ! Verlet
      !particles_cur_vel(:, i) = particles_cur_vel(:, i) &
      !                      & + 0.5d0*( particles_prev_accel(:, i) &
      !                      & + particles_cur_accel(:, i) )*time_step

      !! Beeman
      particles_cur_vel(:, i) = particles_cur_vel(:, i) &
                            & + 1.0d0/6.0d0*( 2.0d0*particles_cur_accel(:, i) &
                            & + 5.0d0*particles_prev_accel(:, i) & 
                            & - particles_prev2_accel(:, i) )*time_step

      q = particles_charge(i)
      k = particles_species(i)
      emit = particles_emitter(i)
      sec  = particles_section(i)

      ! We use OMP ATOMIC here because the index k is not a loop index
      !$OMP ATOMIC UPDATE
      ramo_current(k) = ramo_current(k) + q * E_zunit * particles_cur_vel(3, i)

      ! We use OMP ATOMIC here because the indexes sec and emit are not loop indexes
      !$OMP ATOMIC UPDATE
      ramo_current_emit(sec, emit) = ramo_current_emit(sec, emit) + q * E_zunit * particles_cur_vel(3, i)

      avg_vel(:) = avg_vel(:) + particles_cur_vel(:, i)
    end do
    !$OMP END PARALLEL DO

    avg_mob = sqrt(avg_vel(1)**2 + avg_vel(2)**2 + avg_vel(3)**2) / (-1.0d0*E_z)
  end subroutine Update_Velocity


  ! ----------------------------------------------------------------------------
  ! Acceleration
  ! Update the acceleration for all the particles
  subroutine Calculate_Acceleration_Particles()
    double precision, dimension(1:3) :: force_E, force_c, force_ic, force_ic_N, force_ic_self
    double precision, dimension(1:3) :: pos_1, pos_2, diff, pos_ic, pos_2_per
    double precision                 :: r
    double precision                 :: q_1, q_2
    double precision                 :: im_1, im_2
    double precision                 :: pre_fac_c
    integer                          :: i, j, k_1, k_2, u, v

    ! We do not use GUIDED scheduling in OpenMP here because the inner loop changes size.
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k_1, k_2, pos_1, pos_2, pos_2_per, diff, r, pos_ic) &
    !$OMP& PRIVATE(force_E, force_c, force_ic, force_ic_N, force_ic_self, im_1, q_1, im_2, q_2, pre_fac_c) &
    !$OMP SHARED(nrPart, particles_cur_pos, particles_mass, particles_species, ptr_field_E, Num_per, box_dim) &
    !$OMP SHARED(ptr_Image_Charge_effect, particles_charge, d) &
    !$OMP& SCHEDULE(DYNAMIC, 1) &
    !$OMP& REDUCTION(+:particles_cur_accel)
    do i = 1, nrPart
      ! Information about the particle we are calculating the force/acceleration on
      pos_1 = particles_cur_pos(:, i)
      im_1 = 1.0d0 / particles_mass(i)
      q_1 = particles_charge(i)
      k_1 = particles_species(i)

      ! Acceleration due to electric field
      force_E = q_1 * ptr_field_E(pos_1)

      ! Do image charge, self interaction
      ! It is questionable if the self interaction of the image charge is valid at these
      ! short distances. The escape velocity for an electron at z = 1nm from its image charge partner
      ! is v_esc = e/sqrt(8*pi*epsilon_0*m_e*z) = 355853 m/s. This gives a de Broglie wavelength of
      ! w_l = h/(m_e*v_esc) = 2.04 nm. This equal to the distance between the electron and its
      ! image charge partner.
      !force_ic_self = q_1**2 * div_fac_c * Force_Image_charges_v2(pos_1, pos_1)
      force_ic_self = 0.0d0

      ! Loop over particles from i+1 to nrElec.
      ! There is no need to loop over all particles since
      ! The forces are equal but in opposite directions
      do j = 1, nrPart
        ! Information about the particle that is acting on the particle at pos_1
        pos_2 = particles_cur_pos(:, j)
        !if (particles_mass(j) == 0.0d0) then
        !  print *, 'Hi'
        !end if
        im_2 = 1.0d0 / particles_mass(j)
        q_2 = particles_charge(j)
        k_2 = particles_species(j)

        ! Prefactor for Coloumb's law
        pre_fac_c = q_1*q_2 * div_fac_c ! q_1*q_2 / (4*pi*epsilon)

        do v = -1*Num_per, Num_per, 1 ! x
          do u = -1*Num_per, Num_per, 1 ! y

            !print *, 'u = ', u
            !print *, 'v = ', v

            ! The inner coulomb loop usally goes from j = 1+i, nrPart
            ! We have to skip the self interaction. But we want the periodic part of it.
            !if ((i == j) .and. (u == 0) .and. (v == 0)) then
            !  cycle
            !end if

            ! Shift the position
            pos_2_per(1) = pos_2(1) + v*(box_dim(1) + per_padding)
            pos_2_per(2) = pos_2(2) + u*(box_dim(2) + per_padding)
            pos_2_per(3) = pos_2(3)

            ! Calculate the distance between the two particles
            diff = pos_1 - pos_2_per
            ! There are fours ways to calculate the distance
            ! Number 1: Use the intrinsic function NORM2(v)
            ! Number 2: Use the equation for it sqrt( v(1)**2 + v(2)**2 + v(3)**2 )
            ! Number 3: Or do sqrt( dot_product(v, v) )
            ! Number 4: Or use sqrt( sum(v**2) )
            ! It turns out number 1 is the slowest by far. Number 2, 3 and 4 are
            ! often similar in speed. The difference is small and they fluctuate a lot,
            ! with no clear winner.
            !
            ! We add a small number (length_scale**3) to the results to
            ! prevent a singularity when calulating 1/r**3
            !
            r = sqrt( sum(diff**2) ) + length_scale**3
            !print *, r/length_scale
            !r = sqrt( dot_product(diff, diff) ) + length_scale**3
            !r = NORM2(diff) + length_scale**3

            ! Calculate the Coulomb force
            ! F = (r_1 - r_2) / |r_1 - r_2|^3
            ! F = (diff / r) * 1/r^2
            ! (diff / r) is a unit vector
            force_c = pre_fac_c * diff / r**3
            !print *, force_c
            !print *, im_1
            !print *, force_c*im_1

            ! Do image charge
            force_ic = pre_fac_c * ptr_Image_Charge_effect(pos_1, pos_2_per)

            ! The image charge force of particle i on particle j is the same in the z-direction
            ! but we reverse the x and y directions of the force due to symmetry.
            force_ic_N(1:2) = -1.0d0*force_ic(1:2)
            force_ic_N(3)   = +1.0d0*force_ic(3)

            !print *, force_ic

            ! ! Below plane
            ! pos_ic(1:2) = pos_2(1:2)
            ! pos_ic(3) = -1.0d0*pos_2(3)
            ! diff = pos_1 - pos_ic
            ! r = sqrt( sum(diff**2) ) + length_scale**3
            ! force_ic = (-1.0d0)*pre_fac_c * diff / r**3

            ! ! Above plane
            ! pos_ic(1:2) = pos_2(1:2)
            ! pos_ic(3) = 2*d - pos_2(3)
            ! diff = pos_1 - pos_ic
            ! r = sqrt( sum(diff**2) ) + length_scale**3
            ! force_ic = force_ic + (-1.0d0)*pre_fac_c * diff / r**3


            !!!$OMP CRITICAL(ACCEL_UPDATE)
            !particles_cur_accel(:, i) = particles_cur_accel(:, i) + force_c*im_1 + force_ic*im_1
            if (j /= i) then ! Do not double count!!!
              particles_cur_accel(:, j) = particles_cur_accel(:, j) - force_c * im_2 + force_ic_N * im_2
            end if
            !particles_cur_accel(:, i) = particles_cur_accel(:, i) + force_c * im_1 !+ force_ic     * im_1
            !!!$OMP END CRITICAL(ACCEL_UPDATE)

            !print *, ''

            ! if ((isnan(force_c(1)*im_2) .eqv. .true.) &
            !    & .or. (isnan(force_c(2)*im_2) .eqv. .true.) &
            !    & .or. (isnan(force_c(3)*im_2) .eqv. .true.)) then
            !   print *, 'u = ', u, ', v = ', v
            !   print *, force_c*im_2
            !   print *, pos_1/length_scale
            !   print *, pos_2_per/length_scale
            !   print *, ''
            !   pause
            ! end if
          end do
        end do
      end do

      !!!$OMP CRITICAL(ACCEL_UPDATE)
      particles_cur_accel(:, i) = particles_cur_accel(:, i) + force_E * im_1 + force_ic_self * im_1
      !!!$OMP END CRITICAL(ACCEL_UPDATE)
    end do
    !$OMP END PARALLEL DO

    !print *, 'Accel'
    !print *, particles_cur_accel(:, 1:nrPart)
    !print *, ''
    !print *, 'Pos'
    !print *, particles_cur_pos(:, 1:nrPart)/length_scale
    !print *, 'Done'
    !print *, nrPart
    !stop
  end subroutine Calculate_Acceleration_Particles

  subroutine Write_Acceleration_Test(step)
    integer, intent(in)              :: step
    integer                          :: ud_accel, IFAIL, i
    character(len=1024)              :: filename
    double precision, dimension(1:3) :: par_accel

    print *, 'ACCEL DATA'

    ! Prepare the name of the output file
    ! each file is named accel-0.dt where the number
    ! represents the current time step.
    write(filename, '(a12, i0, a4)') 'accel/accel-', step, '.bin'

    open(newunit=ud_accel, iostat=IFAIL, file=filename, status='replace', action='write', access='STREAM')
    if (IFAIL /= 0) then
      print '(a)', 'Vacuum: ERROR UNABLE TO OPEN file for acceleration'
      print *, filename
      print *, step
    end if

    do i = 1, nrPart
      par_accel(:) = particles_cur_accel(:, i) ! Position of the particle

      ! Write out x, y, z and which emitter the particle came from
      write(unit=ud_accel) par_accel(1), par_accel(2), par_accel(3)
      
      print *, i
      print *, par_accel
      print *, ''
    end do

    close(unit=ud_accel, iostat=IFAIL, status='keep')
  end subroutine Write_Acceleration_Test


  !-----------------------------------------------------------------------------
  ! Find field in a point
  ! Returns the electric field in V/m
  function Calc_Field_at(pos)
    double precision, dimension(1:3)             :: Calc_Field_at
    double precision, dimension(1:3), intent(in) :: pos

    double precision, dimension(1:3) :: force_c, force_tot, force_ic
    double precision, dimension(1:3) :: pos_1, pos_2, diff, pos_2_per
    double precision                 :: r
    double precision                 :: q_1, q_2
    double precision                 :: pre_fac_c
    integer                          :: j, v, u

    ! Position of the particle we are calculating the force/acceleration on
    pos_1 = pos

    ! Electric field in the system
    force_tot = ptr_field_E(pos_1)
    !print *, force_tot

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(j, pos_2, pos_2_per, diff, r, force_c, force_ic, q_2, pre_fac_c) &
    !$OMP& SHARED(nrPart, particles_cur_pos, particles_charge, ptr_Image_Charge_effect, pos_1, Num_per, box_dim) &
    !$OMP& REDUCTION(+:force_tot) SCHEDULE(GUIDED, CHUNK_SIZE)
    do j = 1, nrPart

      ! Position of the particle that is acting on the particle at pos_1
      pos_2 = particles_cur_pos(:, j)
      q_2 = particles_charge(j)

      pre_fac_c = q_2 * div_fac_c ! q_2 / (4*pi*epsilon)

        ! Loop over the periodic systems
        !
        ! ---------------------------------------------------
        ! |         |         |         |         |         |
        ! | -2x,+2y | -1x,+2y | +0x,+2y | +1x,+2y | +2x,+2y |
        ! |         |         |         |         |         |
        ! |--------------------------------------------------
        ! |         |         |         |         |         |
        ! | -2x,+1y | -1x,+1y | +0x,+1y | +1x,+1y | +2x,+1y |
        ! |         |         |         |         |         |
        ! |--------------------------------------------------
        ! |         |         |         |         |         |
        ! | -2x,+0y | -1x,+0y | +0x,+0y | +1x,+0y | +2x,+0y |
        ! |         |         |         |         |         |
        ! |--------------------------------------------------
        ! |         |         |         |         |         |
        ! | -2x,-1y | -1x,-1y | +0x,-1y | +1x,-1y | +2x,-1y |
        ! |         |         |         |         |         |
        ! |--------------------------------------------------
        ! |         |         |         |         |         |
        ! | -2x,-2y | -1x,-2y | +0x,-2y | +1x,-2y | +2x,-2y |
        ! |         |         |         |         |         |
        ! |--------------------------------------------------
      do v = -1*Num_per, Num_per, 1 ! x
        do u = -1*Num_per, Num_per, 1 ! y

          ! Shift the position
          pos_2_per(1) = pos_2(1) + v*(box_dim(1) + per_padding)
          pos_2_per(2) = pos_2(2) + u*(box_dim(2) + per_padding)
          pos_2_per(3) = pos_2(3)

          ! Calculate the distance between the two particles
          diff = pos_1 - pos_2_per
          r = sqrt( sum(diff**2) ) + length_scale**3
          !r = sqrt( dot_product(diff, diff) ) + length_scale**3
          !r = NORM2(diff) + length_scale**3 ! distance + Prevent singularity

          ! Calculate the Coulomb force
          ! F = (r_1 - r_2) / |r_1 - r_2|^3
          ! F = (diff / r) * 1/r^2
          ! (diff / r) is a unit vector
          force_c = diff / r**3
          !force_c = diff / (r*r*r)

          ! Image charge effect
          force_ic = ptr_Image_Charge_effect(pos_1, pos_2_per)

          ! The total force
          force_tot = force_tot + pre_fac_c * force_c + pre_fac_c * force_ic
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    Calc_Field_at = force_tot

  end function Calc_Field_at

  !-----------------------------------------------------------------------------
  ! We have a particle at location pos_1 and we want to calculate the effects of
  ! the image charge partners of the particle at location pos_2.
  ! For particle 2, its image charge partners are located at the same x and y coordinates
  ! but the z coordinates change. If the particle is at z_0 then,
  ! for partners with the opposite charge we have z_n = 2*n*d - z_0, n = 0,±1,±2,±3,... ,
  ! for partners with the same charge we have z_n = 2*n*d + z_0, n = ±1, ±2, ±3, ... .
  !
  ! This function return
  ! F = (pos_1 - pos_ic)/r**3
  ! i.e. without q_1*q_2/(4\pi\epsilon_0)
  function Force_Image_charges_v2(pos_1, pos_2)
    double precision, intent(in), dimension(1:3) :: pos_1, pos_2
    double precision, dimension(1:3)             :: Force_Image_charges_v2
    integer                                      :: n
    double precision, dimension(1:3)             :: pos_ic, diff
    double precision                             :: r

    ! Check if we are doing image charge or not
    if (image_charge .eqv. .false.) then
      ! Return 0 if we are not using image charge
      Force_Image_charges_v2 = 0.0d0
    else
      ! Start with n = 0
      n = 0
      pos_ic(1:2) = pos_2(1:2) ! Same x and y
      pos_ic(3) = -1.0d0*pos_2(3) ! Change z

      diff = pos_1 - pos_ic
      r = sqrt( sum(diff**2) ) + length_scale**2
      Force_Image_charges_v2 = (-1.0d0)*diff/r**3 ! -1.0d0 because of the opposite charge

      do n = 1, N_ic_max
        ! The charges with the opposite charges first
        ! Plus n
        pos_ic(3) = 2.0d0*n*d - pos_2(3) ! Change z
        diff = pos_1 - pos_ic
        r = sqrt( sum(diff**2) ) + length_scale**2
        Force_Image_charges_v2 = Force_Image_charges_v2 + (-1.0d0)*diff/r**3 ! -1.0d0 because of the opposite charge

        ! Negative n
        pos_ic(3) = -2.0d0*n*d - pos_2(3) ! Change z
        diff = pos_1 - pos_ic
        r = sqrt( sum(diff**2) ) + length_scale**2
        Force_Image_charges_v2 = Force_Image_charges_v2 + (-1.0d0)*diff/r**3 ! -1.0d0 because of the opposite charge

        ! Now do the charges with the same charge
        ! Plus n
        pos_ic(3) = 2.0d0*n*d + pos_2(3) ! Change z
        diff = pos_1 - pos_ic
        r = sqrt( sum(diff**2) ) + length_scale**2
        Force_Image_charges_v2 = Force_Image_charges_v2 + (+1.0d0)*diff/r**3 ! +1.0d0 because of the same charge

        ! Negative n
        pos_ic(3) = -2.0d0*n*d + pos_2(3) ! Change z
        diff = pos_1 - pos_ic
        r = sqrt( sum(diff**2) ) + length_scale**2
        Force_Image_charges_v2 = Force_Image_charges_v2 + (+1.0d0)*diff/r**3 ! +1.0d0 because of the same charge
      end do
    end if
  end function Force_Image_charges_v2

  ! ----------------------------------------------------------------------------
  ! The vacuum electric field
  pure function field_E_planar(pos) result(field_E)
    double precision, dimension(1:3), intent(in) :: pos
    double precision, dimension(1:3)             :: field_E

    ! Electric field
    ! x
    field_E(1) = 0.0d0

    ! y
    field_E(2) = 0.0d0

    ! z
    field_E(3) = E_z

  end function field_E_planar


  !-----------------------------------------------------------------------------
  ! The voltage
  subroutine Set_Voltage(step)
    integer, intent(in) :: step
    integer             :: IFAIL
    ! !double precision    :: V_R, V_C
    ! !double precision    :: ramo_cur

    ! ! Calculate the total ramo current
    ! !ramo_cur = sum(ramo_current)
    ! !I_prev = I_cur
    ! !I_cur = ramo_cur

    ! !I = Get_Current_Const(step)
    ! I_prev = I_cur
    ! I_cur = sum(ramo_current)

    ! !V_R = Voltage_Resistor()
    ! !V_C = Voltage_Capacitor()
    ! !V_C = 0.0d0
    ! !V_d = V_s + V_R + V_C
    
    ! ! if (step == 1) then
    ! !   print *, step
    ! !   print *, V_rf
    ! ! end if

    ! if (step == 0) then
    !   I_prev = 0.0d0
    !   V_rf_prev = 0.0d0
    !   V_rf = V_s
    ! end if
    
    ! V_rf_next = Parallel_RLC_FD(step, I_cur, I_prev, V_rf, V_rf_prev)
    ! V_rf_prev = V_rf
    ! V_rf = V_rf_next

    ! V_d = V_s + V_rf ! DC voltage + RF voltage
    ! !V_d = V_s

    ! !V_C = Voltage_Parallel_Capacitor(step)
    ! !V_d = V_C

    ! !V_d = V_s

    ! E_z = -1.0d0*V_d/d
    ! !E_zunit = -1.0d0*sign(1.0d0, V)/d

    V_d = V_s
    V_rf = 0.0d0
    E_z = -1.0d0*V_d/d

    write (ud_volt, "(ES12.4, tr2, i8, tr2, ES18.8, tr2, ES18.8)", iostat=IFAIL) &
          & cur_time, step, V_d, V_rf
  end subroutine Set_Voltage

  double precision pure function Voltage_Resistor()
    double precision            :: ramo_cur

    ! Calculate the total ramo current
    ramo_cur = sum(ramo_current)

    ! Calculate the voltage drop over the resistor
    Voltage_Resistor = -1.0d0*R_s*ramo_cur
  end function Voltage_Resistor

  ! double precision function Get_Current_Const(step)
  !   integer, intent(in) :: step
  !   Get_Current_Const = 10.0E-3 ! Amper
  ! end function Get_Current_Const

  double precision function Parallel_RLC_FD(step, I_cur, I_prev, V_cur, V_prev)
    integer, intent(in) :: step
    double precision, intent(in) :: V_cur, V_prev ! Current V(t) and previous V(t-Δt) values of the voltage
    double precision, intent(in) :: I_cur, I_prev ! Current I(t) and previous I(t-Δt) values of the current

    !double precision, parameter  :: R = 2.0d0   ! Ohm
    !double precision, parameter  :: L = 1.0d-7  ! Henry
    !double precision, parameter  :: C = 1.0d-22 ! Farad

    !double precision, parameter  :: R_p = 7500.0d0 ! Ohm
    !double precision, parameter  :: L_p = 1.04d-9  ! Henry
    !double precision, parameter  :: C_p = 0.53d-17 ! Farad

    double precision :: RC
    double precision :: LC

    !double precision, parameter :: R = 13.5d0
    !double precision, parameter :: omega_0 = 2.0d0*pi*3.1d9 ! rad/s [omega_0 = 1/sqrt(LC)]
    !double precision, parameter :: Q = 150.0d0 ! Dimensionless [Q = R*sqrt(C/L)]
    
    double precision :: V_next ! Next value of the voltage V(t+Δt)

    RC = R_p*C_p
    LC = L_p*C_p

    V_next = time_step/C_p * I_cur + V_cur * (2.0d0 + time_step/RC - time_step2/LC) &
         & - V_prev - time_step/C_p * I_prev
    V_next = V_next / (1.0d0 + time_step/RC)

    !V_next = time_step*R*omega_0/Q * I_cur + V_cur * (2.0d0 + time_step*omega_0/Q - time_step2*omega_0**2) &
    !     & - V_prev - time_step*R*omega_0/Q * I_prev
    !V_next = V_next / (1.0d0 + time_step*omega_0/Q)

    Parallel_RLC_FD = V_next

    ! if (step < 2) then
    !   print *, 'Parallel_RLC_FD'
    !   print *, I_cur
    !   print *, I_prev
    !   print *, V_cur
    !   print *, V_prev
    !   print *, time_step

    !   print *, time_step/C
    !   print *, (2.0d0 + time_step/(R*C) - time_step2/(L*C))
    !   print *, (1.0d0 + time_step/(R*C))
    ! end if
  end function Parallel_RLC_FD

end module mod_verlet
