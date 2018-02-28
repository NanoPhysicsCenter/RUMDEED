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

    ! Update the current time
    cur_time = time_step * step / time_scale ! Scaled in units of time_scale

    ! Update the voltage in the system
    call Set_Voltage(step)

    ! Reset the ramo current. Note, this should be done after set_voltage,
    ! since the ramo current may be used there.
    ramo_current = 0.0d0


    ! Update the position of particles (Electrons / Holes)
    if (nrElecHole > 0) then
      call Update_ElecHole_Position(step)

      call Calculate_Acceleration_Particles()

      call Update_Velocity(step)
    end if

    ! Write out the current in the system
    call Write_Ramo_Current(step)

    ! Write out the field in the system
    !call Write_Field_Longitudinal(step)
    !call Write_Density_Map(step)

  end subroutine Velocity_Verlet

  ! ----------------------------------------------------------------------------
  !
  subroutine Update_ElecHole_Position(step)
    integer, intent(in) :: step
    integer             :: i

    !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC)
    do i = startElecHoles, endElecHoles
      particles_prev_pos(:, i) = particles_cur_pos(:, i) ! Store the previous position
      particles_cur_pos(:, i)  = particles_cur_pos(:, i) + particles_cur_vel(:, i)*time_step &
                             & + 0.5d0*particles_cur_accel(:, i)*time_step2

      particles_prev_accel(:, i) = particles_cur_accel(:, i)
      particles_cur_accel(:, i)  = 0.0d0

      ! Mark particles that should be removed with .false. in the mask array
      !call Check_Boundary_ElecHole(i)
      call ptr_Check_Boundary(i)
    end do
    !$OMP END PARALLEL DO
  end subroutine Update_ElecHole_Position

  ! ----------------------------------------------------------------------------
  ! Checks the boundary conditions of the box.
  ! Check which particles to remove
  ! Enforce periodic boundary conditions (ToDo)
  subroutine Check_Boundary_ElecHole(i)
    integer, intent(in) :: i
    double precision    :: z

    z = particles_cur_pos(3, i)

    ! Check if the particle should be removed from the system
    if (z < 0.0d0) then
        call Mark_Particles_Remove(i, remove_bot)
    else if (z > box_dim(3)) then
        call Mark_Particles_Remove(i, remove_top)
    end if

  end subroutine Check_Boundary_ElecHole


  ! ----------------------------------------------------------------------------
  !
  subroutine Update_Velocity(step)
    integer, intent(in) :: step
    integer             :: i, k
    double precision    :: q

    !$OMP PARALLEL DO PRIVATE(i, q, k)
    do i = startElecHoles, endElecHoles
      particles_cur_vel(:, i) = particles_cur_vel(:, i) &
                            & + 0.5d0*( particles_prev_accel(:, i) &
                            & + particles_cur_accel(:, i) )*time_step

      q = particles_charge(i)
      k = particles_species(i)

      ! We use OMP ATOMIC here because the index k is not a loop index
      !$OMP ATOMIC
      ramo_current(k) = ramo_current(k) + q * E_zunit * particles_cur_vel(3, i)
    end do
    !$OMP END PARALLEL DO
  end subroutine Update_Velocity


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

    !$OMP PARALLEL DO PRIVATE(i, j, k_1, k_2, pos_1, pos_2, diff, r, force_E, force_c, im_1, q_1, im_2, q_2, pre_fac_c) &
    !$OMP& REDUCTION(+:particles_cur_accel) SCHEDULE(GUIDED)
    do i = 1, nrPart
      ! Information about the particle we are calculating the force/acceleration on
      pos_1 = particles_cur_pos(:, i)
      im_1 = 1.0d0 / particles_mass(i)
      q_1 = particles_charge(i)
      k_1 = particles_species(i)

      ! Acceleration due to electric field
      force_E = q_1 * ptr_field_E(pos_1)

      ! Loop over particles from i+1 to nrElec.
      ! There is no need to loop over all particles since
      ! The forces are equal but in opposite directions
      do j = i+1, nrPart

        ! Information about the particle that is acting on the particle at pos_1
        pos_2 = particles_cur_pos(:, j)
        if (particles_mass(j) == 0.0d0) then
          print *, 'Hi'
        end if
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

        particles_cur_accel(:, j) = particles_cur_accel(:, j) - pre_fac_c * force_c * im_2
        particles_cur_accel(:, i) = particles_cur_accel(:, i) + pre_fac_c * force_c * im_1
      end do

      particles_cur_accel(:, i) = particles_cur_accel(:, i) + force_E * im_1
    end do
    !$OMP END PARALLEL DO
  end subroutine Calculate_Acceleration_Particles


  !-----------------------------------------------------------------------------
  ! Find field in point
  ! Return the electric field in V/m
  function Calc_Field_at(pos)
    double precision, dimension(1:3)             :: Calc_Field_at
    double precision, dimension(1:3), intent(in) :: pos

    double precision, dimension(1:3) :: force_c, force_tot
    double precision, dimension(1:3) :: pos_1, pos_2, diff
    double precision                 :: r
    double precision                 :: q_2
    double precision                 :: pre_fac_c
    integer                          :: j

    ! Position of the particle we are calculating the force/acceleration on
    pos_1 = pos

    ! Electric field in the system
    force_tot = ptr_field_E(pos_1)

    !$OMP PARALLEL DO PRIVATE(j, pos_2, diff, r, force_c, q_2, pre_fac_c) &
    !$OMP& REDUCTION(+:force_tot) SCHEDULE(GUIDED)
    do j = 1, nrPart

      ! Position of the particle that is acting on the particle at pos_1
      pos_2 = particles_cur_pos(:, j)
      q_2 = particles_charge(j)

      pre_fac_c = q_2 * div_fac_c ! q_2 / (4*pi*epsilon)

      ! Calculate the distance between the two particles
      diff = pos_1 - pos_2
      r = NORM2(diff) + length_scale**3 ! distance + Prevent singularity

      ! Calculate the Coulomb force
      ! F = (r_1 - r_2) / |r_1 - r_2|^3
      ! F = (diff / r) * 1/r^2
      ! (diff / r) is a unit vector
      force_c = diff * r**(-3)

      force_tot = force_tot + pre_fac_c * force_c
    end do
    !$OMP END PARALLEL DO

    Calc_Field_at = force_tot

  end function Calc_Field_at


  ! ----------------------------------------------------------------------------
  ! The vacuum electric field
  pure function field_E_planar(pos)
    double precision, dimension(1:3), intent(in) :: pos
    double precision, dimension(1:3)             :: field_E_planar

    ! Electric field
    ! x
    field_E_planar(1) = 0.0d0

    ! y
    field_E_planar(2) = 0.0d0

    ! z
    field_E_planar(3) = E_z

  end function field_E_planar


  !-----------------------------------------------------------------------------
  ! The voltage
  subroutine Set_Voltage(step)
    integer, intent(in) :: step
    integer             :: IFAIL
    double precision    :: V_R = 0.0d0, V_C = 0.0d0

    !V_R = Voltage_Resistor()
    !V_C = Voltage_Capacitor()
    !V_d = V_s + V_R + V_C

    !V_C = Voltage_Parallel_Capacitor(step)
    !V_d = V_C

    V_d = V_s

    E_z = -1.0d0*V_d/d
    !E_zunit = -1.0d0*sign(1.0d0, V)/d

    write (ud_volt, "(ES12.4, tr2, i8, tr2, ES12.4, tr2, ES12.4, tr2, ES12.4)", iostat=IFAIL) cur_time, step, V_d, V_R, V_C
  end subroutine Set_Voltage

  double precision pure function Voltage_Resistor()
    double precision, parameter :: R = 5.0d5 ! Ohm
    double precision            :: ramo_cur

    ! Calculate the total ramo current
    ramo_cur = sum(ramo_current)

    ! Calculate the voltage drop over the resistor
    Voltage_Resistor = -1.0d0*R*ramo_cur
  end function Voltage_Resistor

  double precision function Voltage_Capacitor()
    double precision, parameter :: C = 10.0E-18 ! Farad (Atto Farads?)
    double precision            :: ramo_cur

    ! Calculate the total ramo current
    ramo_cur = sum(ramo_current)

    ! Colculate the voltage drop over the capacitor
    ramo_integral = ramo_integral + time_step * (ramo_cur_prev + ramo_cur) * 0.5d0
    ramo_cur_prev = ramo_cur
    Voltage_Capacitor = -1.0d0/C * ramo_integral
  end function Voltage_Capacitor

  double precision function Voltage_Parallel_Capacitor_Resistor()
    double precision, parameter :: C   = 10.0d-18 ! Farad
    double precision, parameter :: R_C = 0.5d6 ! Ohm
    double precision, parameter :: R   = 0.5d6 ! Ohm
    double precision, parameter :: ib   = 1.0d0/(C*(R+R_C))
    double precision            :: ramo_cur

    ! Calculate the total ramo current
    ramo_cur = sum(ramo_current)

    ! Caclulate the voltage over the diode
    ramo_integral = ramo_integral + time_step * (ramo_cur_prev*exp(-1.0d0*time_step*ib) + ramo_cur) * 0.5d0
    ramo_cur_prev = ramo_cur

    Voltage_Parallel_Capacitor_Resistor = (R**2)/(C*(R+R_C)**2) * ramo_integral - R*R_C/(R+R_C)*ramo_cur
  end function Voltage_Parallel_Capacitor_Resistor

  double precision function Voltage_Parallel_Capacitor(step)
    integer, intent(in)         :: step
    double precision, parameter :: C = 1.0E-15 ! Farad
    double precision, parameter :: R = 1.0d3 ! Ohm
    double precision, parameter :: RC = R*C
    double precision, parameter :: iRC = 1.0d0/RC
    double precision            :: ramo_cur, cur_time_s

    ! Calculate the total ramo current
    ramo_cur = sum(ramo_current)

    ! Current time in seconds
    cur_time_s = time_step * step

    ! Caclulate the voltage over the diode
    ramo_integral = ramo_integral + time_step * (ramo_cur_prev*exp(-1.0d0*time_step*iRC) + ramo_cur) * 0.5d0
    ramo_cur_prev = ramo_cur

    Voltage_Parallel_Capacitor = V_s*exp(-1.0d0*cur_time_s*iRC)*( RC*(exp(cur_time_s*iRC) - 1.0d0) + 0.0d0 ) &
                               - 1.0d0/C * ramo_integral
  end function Voltage_Parallel_Capacitor

  double precision function Parallel_Capacitor_MNA(step, I_D)
    integer, intent(in)          :: step
    double precision, intent(in) :: I_D
    double precision, parameter  :: C = 10.0d-18 ! Farad
    double precision, parameter  :: R_D = 1.0d6  ! Ohm
    double precision, parameter  :: R_S = 1.0d6  ! Ohm
    double precision, parameter  :: R_C = 1.0d6  ! Ohm

    double precision, dimension(1:5, 1:5) :: A, A_inv
    double precision, dimension(1:5)      :: b
    logical                               :: OK_FLAG

    A = 0.0d0
    A_inv = 0.0d0
    b = 0.0d0

    A(1, 1) = 1.0d0/R_D + 1.0d0/R_C
    A(1, 2) = -1.0d0/R_D
    A(1, 3) = -1.0d0/R_C
    A(1, 4) = 1.0d0/R_S
    b(1)    = 0.0

    A(2, 1) = 1.0d0
    A(2, 4) = -1.0d0
    b(2)    = V_S

    A(3, 1) = -1.0d0
    A(3, 2) = 1.0d0
    b(3)    = -1.0d0*I_D*R_D

    A(4, 3) = 2.0d0*C/time_step
    A(4, 5) = -1.0d0
    b(4)    = V_prev(4) + 2.0d0*C/time_step*V_prev(2)

    A(5, 4) = 1.0d0/R_S
    A(5, 5) = 1.0d0
    b(5)    = -I_D

    ! Store previous values of the voltages
    V_prev = V_cur

    ! Solve the system of equations Ax=b or x=A^-1*b
    ! Find the inverse of A
    call M55INV(A, A_inv, OK_FLAG)

    ! x = A^-1*b
    V_cur = matmul(A_inv, b)

    ! Calculate voltages
    !V_D(step) = V_cur(1)             ! Voltage over the diode
    !V_C(step) = V_cur(2)             ! Voltage of the capacitor
    !V_SC(step) = V_cur(0) - V_cur(3) ! Source voltage (Should be equal to V_S)

    ! Calculate currents from voltages over resistors
    !I_T(step)  = ( -1.0d0*V_cur(3) / R_S ) / 1.0d-6          ! Total current
    !!I_C(step)  = ( (V_cur(0) - V_cur(2)) / R_C ) / 1.0d-6  ! Capacitor current
    !I_C(step)  = V_cur(4) / 1.0d-6                         ! Capacitor current
    !I_DC(step) = ( (V_cur(0) - V_cur(1)) / R_D ) / 1.0d-6  ! Diode current

    !I_D = Current_Diode_Child(V_D(step), step)
    !I_D = Current_Diode(V_D[step], step)

    ! Store the current time
    !time = step*time_step / time_scale
  end function Parallel_Capacitor_MNA

end module mod_verlet
