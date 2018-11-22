!-------------------------------------------!
! Module for emission from tip              !
! Kristinn Torfason                         !
! 04.06.18                                  !
!-------------------------------------------!

Module mod_field_emission_tip
  use mod_global
  use mod_hyperboloid_tip
  use mod_verlet
  use mod_pair
  use ieee_arithmetic
  implicit none

  ! ----------------------------------------------------------------------------
  ! Variables
  integer                            :: posInit
  integer                            :: nrEmitted
  double precision                   :: res_s ! residual
  double precision                   :: d_tip

  ! ----------------------------------------------------------------------------
  ! Constants for field emission
  ! Fyrst Fowler-Nordheim constant in units [ A eV V^{-2} ]
  double precision, parameter :: a_FN = q_02/(16.0d0*pi**2*h_bar)

  ! Second Fowler-Nodheim constant in units [ eV^{-3/2} V m^{-1} ]
  double precision, parameter :: b_FN = -4.0d0/(3.0d0*h_bar) * sqrt(2.0d0*m_0*q_0)

  ! Constant used for calculation of the l in v_y and t_y.
  ! The units are [ eV^{2} V^{-1} m ]
  ! See Forbes, R. G., & Deane, J. H. (2007, November).
  ! "Reformulation of the standard theory of Fowler–Nordheim tunnelling and cold field electron emission."
  ! In Proceedings of the Royal Society of London A: Mathematical,
  ! Physical and Engineering Sciences (Vol. 463, No. 2087, pp. 2907-2927). The Royal Society.
  double precision, parameter :: l_const = q_0 / (4.0d0*pi*epsilon_0)

  ! The work function. Unit [ eV ]
  double precision, parameter :: w_theta = 4.7d0

  ! Use image Charge or not
  !logical, parameter          :: image_charge = .true.

contains

  !-----------------------------------------------------------------------------
  ! Initialize the Field Emission
  subroutine Init_Field_Emission_Tip()

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_ElecHole_Tip

    ! Function for the electric field in the system
    ptr_field_E => field_E_Hyperboloid

    ! The function that does the emission
    ptr_Do_Emission => Do_Field_Emission_Tip


    ! Tip Parameters

    d_tip = emitters_dim(1, 1)
    R_base = emitters_dim(2, 1)
    h_tip = emitters_dim(3, 1)

    d = d_tip + h_tip
    E_z = -1.0d0*V_s/d
    E_zunit = -1.0d0/d

    ! Other Parameters
    max_xi = h_tip/d_tip + 1.0d0
    a_foci = sqrt(d_tip**2*R_base**2 / (h_tip**2 + 2*d_tip*h_tip) + d_tip**2)
    eta_1 = -1.0d0 * d_tip / a_foci

    theta_tip = acos(d_tip/a_foci)
    r_tip = a_foci * sin(theta_tip) * tan(theta_tip)
    shift_z = abs(a_foci * eta_1 * max_xi)

    pre_fac_E_tip = log( (1.0d0 + eta_1)/(1.0d0 - eta_1) * (1.0d0 - eta_2)/(1.0d0 + eta_2) )
    pre_fac_E_tip = 2.0d0 * V_s / (a_foci * pre_fac_E_tip)

  end subroutine Init_Field_Emission_Tip

  subroutine Clean_Up_Field_Emission_Tip()
    ! Nothing to do here
  end subroutine Clean_Up_Field_Emission_Tip

  subroutine Do_Field_Emission_Tip(step)
    integer, intent(in) :: step
    integer             :: IFAIL

    posInit = 0
    nrEmitted = 0

    if (emitters_Type(1) == 1) then
      call Do_Field_Emission_Tip_1(step)
    else
      ! Reverse the voltage
      call Do_Field_Emission_Tip_2(step)
    end if

    cur_time = time_step * step / time_scale ! Scaled in units of time_scale
    write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, posInit, &
    & nrEmitted, nrElec
  end subroutine Do_Field_Emission_Tip

!----------------------------------------------------------------------------------------
subroutine Do_Field_Emission_Tip_2(step)
  integer, intent(in)              :: step
  double precision, dimension(1:3) :: par_pos, par_vel
  integer                          :: nrElecEmit, emit
  double precision                 :: x, y

  nrElecEmit = 0

  par_pos(1:2) = box_muller((/0.0d0, 0.0d0/), (/10.0d0*length_scale, 10.0d0*length_scale/)) ! Random x and y
  par_pos(3) = d - 1.0d0*length_scale ! 1 nm below the plane in z

  x = par_pos(1)
  y = par_pos(2)
  if ((x > 0.0d0) .and. (y > 0.0d0)) then
    emit = 1
  else if ((x > 0.0d0) .and. (y < 0.0d0)) then
    emit = 2
  else if ((x < 0.0d0) .and. (y > 0.0d0)) then
    emit = 3
  else
    emit = 4
  end if

  par_vel = 0.0d0
  call Add_Particle(par_pos, par_vel, species_elec, step, emit)
  nrElecEmit = nrElecEmit + 1

  posInit = posInit + nrElecEmit
  nrEmitted = nrEmitted + nrElecEmit
end subroutine Do_Field_Emission_Tip_2

!----------------------------------------------------------------------------------------

  subroutine Do_Field_Emission_Tip_1(step)
    integer, intent(in)              :: step
    double precision                 :: F
    double precision, dimension(1:3) :: par_pos, surf_norm, par_vel
    !double precision, dimension(1)   :: rnd
    double precision                 :: rnd
    integer                          :: i, j, s, IFAIL, nrElecEmit, n_r
    double precision                 :: A_f, D_f, n_s, F_avg, n_add
    double precision                 :: len_phi, len_xi
    integer                          :: nr_phi, nr_xi, ndim, emit
    double precision                 :: xi_1, phi_1, xi_2, phi_2, xi_c, phi_c
    double precision, dimension(1:3) :: field

    !!!$OMP SINGLE
    emit = 0

    nr_phi = 100
    nr_xi = nr_phi
    len_phi = 2.0d0*pi / nr_phi
    len_xi = (max_xi - 1.0d0) / nr_xi

    !print *, len_x, len_x / length
    !print *, len_y, len_y / length

    nrElecEmit = 0

    !par_pos = 0.0d0
    n_s = 0.0d0
    !F_avg = 0.0d0

    !!!$OMP END SINGLE

    !!!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, par_pos, par_elec, xi_c, phi_c, xi_1, phi_1, xi_2, phi_2, A_f, F, n_add, s, D_f, surf_norm)

    !$OMP PARALLEL DO PRIVATE(i, j, xi_c, phi_c, par_pos, field, F, n_add, xi_1, xi_2, phi_1, phi_2, A_f) &
    !$OMP& REDUCTION(+:n_s,F_avg)
    do i = 1, nr_xi
      do j = 1, nr_phi
        xi_c = 1.0d0 + (i - 0.5d0)*len_xi
        phi_c = (j - 0.5d0)*len_phi

        par_pos(1) = x_coor(xi_c, eta_1, phi_c)
        par_pos(2) = y_coor(xi_c, eta_1, phi_c)
        par_pos(3) = z_coor(xi_c, eta_1, phi_c)

        field = Calc_Field_at(par_pos)
        F = Field_normal(par_pos, field)

        F_avg = F_avg + F


        if (F >= 0.0d0) then
          n_add = 0.0d0
        else
          xi_1 = 1.0d0 + (i - 1.0d0)*len_xi
          phi_1 = (j - 1.0d0)*len_phi
          xi_2 = 1.0d0 + (i + 0.0d0)*len_xi
          phi_2 = (j + 0.0d0)*len_phi
          A_f = Tip_Area(xi_1, xi_2, phi_1, phi_2)
          n_add = Elec_supply(A_f, F, par_pos)
          !if (isnan(n_add) == .true.) then
          !if (n_add < 0.0d0) then
          !  print *, 'A_f = ', A_f
          !  print *, 'F = ', F
          !  print *, 'n_add = ', n_add
            !stop
          !end if
        end if

        n_s = n_s + n_add
        !!$OMP CRITICAL
        !print *, n_s
        !!$OMP END CRITICAL
      end do
    end do
    !$OMP END PARALLEL DO

    !!$OMP SINGLE
    !print *, 'F_avg = ', F_avg
    F_avg = F_avg / (nr_phi*nr_xi)
    write (ud_debug, "(i8, tr2, E16.8)", iostat=IFAIL) step, F_avg

    !n_s = n_s - res_s
    n_r = nint(n_s)
    !res_s = n_r - n_s

    !print *, 'n_r = ', n_r
    if (n_r < 0) then
      print *, 'n_r < 0'
      print *, 'n_s = ', n_s
      print *, 'n_r = ', n_r
      stop
    end if

    !!!$OMP END SINGLE

    print *, 'Doing emission'
    print *, n_r

    !$OMP PARALLEL DO PRIVATE(s, ndim, par_pos, field, F, D_f, surf_norm, xi_1, phi_1)
    do s = 1, n_r

      !!!$OMP FLUSH (particles, nrElec)
      ndim = 25
      par_pos = Metro_algo_tip(ndim, xi_1, phi_1)
      !par_pos = Metro_algo_rec(ndim)
      field = Calc_Field_at(par_pos)
      F = Field_normal(par_pos, field)

      if (F >= 0.0d0) then
        !Try again
        !par_pos = Metro_algo_rec(ndim)
        par_pos = Metro_algo_tip(ndim, xi_1, phi_1)
        field = Calc_Field_at(par_pos)
        F = Field_normal(par_pos, field)

        if (F >= 0.0d0) then
          D_f = 0.0d0
          print *, 'Warning: F > 0.0d0'
        end if
      else
        D_f = Escape_Prob(F, par_pos)
      end if

      CALL RANDOM_NUMBER(rnd)
      print *, 'rnd ', rnd
      print *, 'D_f ', D_f
      print *, ''
      if (rnd <= D_f) then
        !par_pos(3) = par_pos(3) + 1.0d0*length
        surf_norm = surface_normal(par_pos)
        par_pos = par_pos + surf_norm*length_scale
        !$OMP CRITICAL
          !call Add_kdens_Graph_emit(par_pos, step)
          !call Add_Plane_Graph_emit(par_pos, step)
          !call Add_J_tip(xi_1, eta_1, phi_1, step)
          !call Add_Plane_Graph_xi(xi_1, eta_1, phi_1, step)
          !call Add_Plane_Graph_arc(xi_1, eta_1, phi_1, step)
          !write (ud_debug2, "(E16.8, tr2, E16.8, tr2, E16.8)", iostat=IFAIL) par_pos(1), par_pos(2), par_pos(3)
          !particles(nrElec+1)%cur_pos = par_pos

          ! Add a particle to the system
          par_vel = 0.0d0
          call Add_Particle(par_pos, par_vel, species_elec, step, emit)

          !nrElec = nrElec + 1
          nrElecEmit = nrElecEmit + 1
          print *, 'Particle emitted'
        !$OMP END CRITICAL
      end if
    end do
    !$OMP END PARALLEL DO

    !!!$OMP MASTER

    posInit = posInit + nrElecEmit
    nrEmitted = nrEmitted + nrElecEmit
    print *, 'Emission done'
    print *, posInit
    print *, ''
    !!!$OMP END MASTER

    !!!$OMP END PARALLEL
  end subroutine Do_Field_Emission_Tip_1

  function Metro_algo_tip(ndim, xi, phi)
    integer, intent(in)                          :: ndim
    double precision, intent(out)                :: xi, phi
    double precision, dimension(1:3)             :: Metro_algo_tip
    double precision, dimension(ndim, 1:2)       :: r
    double precision                             :: rnd
    double precision, dimension(1:3)             :: par_pos, new_pos
    double precision, dimension(1:2, 1:2)        :: T
    double precision                             :: old_val, new_val, alpha, eta_f
    integer                                      :: old_val_s, new_val_s
    double precision                             :: d_xi
    integer                                      :: i, IFAIL = 0
    integer                                      :: ud_metro
    double precision, dimension(1:3)             :: field

    !open(newunit=ud_metro, iostat=IFAIL, file='metro.dt', status='REPLACE', action='write')

    T = 0.0d0
    !T(1, 1) = 0.5d0 ! xi
    !T(2, 2) = 0.5d0 ! phi
    T(1, 1) = (max_xi - 1.0d0) * 0.05d0
    T(2, 2) = 2.0d0*pi * 0.05d0

    par_pos(1) = 0.0d0
    par_pos(2) = 0.0d0

    !IFAIL = vdRngGaussianMV(VSL_RNG_METHOD_GAUSSIANMV_ICDF, stream, ndim, r(:, 1:2), 2, VSL_MATRIX_STORAGE_FULL, par_pos(1:2), T)

    !print *, 'Random start'

    ! IFAIL = vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, ndim, r(:, 1), par_pos(1), T(1, 1))
    ! if (IFAIL /= 0) then
    !   print *, 'IFAIL'
    !   print *, 'ndim = ', ndim
    ! end if
    !
    ! IFAIL = vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, ndim, r(:, 2), par_pos(2), T(2, 2))
    ! if (IFAIL /= 0) then
    !   print *, 'IFAIL'
    !   print *, 'ndim = ', ndim
    ! end if
    !
    ! IFAIL = vdrnguniform(VSL_RNG_METHOD_UNIFORM_STD, stream, ndim, alpha_r(:), 0.0d0, 1.0d0)
    ! if (IFAIL /= 0) then
    !   print *, 'IFAIL'
    !   print *, 'ndim = ', ndim
    ! end if
    !
    ! IFAIL = vdrnguniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 2, par_pos(1:2), 0.0d0, 1.0d0)
    ! if (IFAIL /= 0) then
    !   print *, 'IFAIL'
    !   print *, 'ndim = ', ndim
    ! end if

    do i = 1, ndim
      par_pos(1:2) = box_muller((/0.0d0, 0.0d0/), (/T(1,1), T(1, 1)/))
      r(i, 1) = par_pos(1)
      par_pos(1:2) = box_muller((/0.0d0, 0.0d0/), (/T(2, 2), T(2, 2)/))
      r(i, 2) = par_pos(2)
    end do

    !print *, 'Random end'
    !print *, ''
    par_pos(1:2) = box_muller((/0.0d0, 0.0d0/), (/T(1, 1), T(1, 1)/))
    xi = (max_xi - 1.0d0)*par_pos(1) + 1.0d0

    par_pos(1:2) = box_muller((/0.0d0, 0.0d0/), (/T(2, 2), T(2, 2)/))
    phi = 2.0d0*pi*par_pos(2)
    !xi = 1.0d0
    !phi = 0.0d0

    par_pos(1) = x_coor(xi, eta_1, phi)
    par_pos(2) = y_coor(xi, eta_1, phi)
    par_pos(3) = z_coor(xi, eta_1, phi)

    field = Calc_Field_at(par_pos)
    eta_f = Field_normal(par_pos, field)

    !old_val = norm2(par_elec%accel)
    old_val = eta_f
    if ((eta_f < 0.0d0) .and. (field(3) < 0.0d0)) then
      old_val_s = +1
    else
      old_val_s = -1
    end if

    new_pos = 0.0d0
    do i = 1, ndim

      xi = xi + r(i, 1)
      phi = phi + r(i, 2)

      call Check_xi(xi)

      ! Keep \phi between 0 and 2\pi
      if (phi > 2.0d0*pi) then
        phi = phi - 2.0d0*pi
      end if
      if (phi < 0.0d0) then
        phi = phi + 2.0d0*pi
      end if

      new_pos(1) = x_coor(xi, eta_1, phi)
      new_pos(2) = y_coor(xi, eta_1, phi)
      new_pos(3) = z_coor(xi, eta_1, phi)

      field = Calc_Field_at(par_pos)
      eta_f = Field_normal(par_pos, field)
      !new_val = norm2(par_elec%accel)
      new_val = eta_f

      ! Check if this point is favourable for emission
      if ((eta_f < 0.0d0) .and. (field(3) < 0.0d0)) then
        new_val_s = +1
      else
        new_val_s = -1
      end if

      if ((old_val_s < 0) .and. (new_val_s < 0)) then
        alpha = old_val / new_val
      else if ((old_val_s > 0) .and. (new_val_s < 0)) then
        alpha = 0.0d0 ! Do not jump to this point
      else if ((old_val_s < 0) .and. (new_val_s > 0)) then
        alpha = 1.0d0 ! Do jump this point
      else
        alpha = new_val / old_val
      end if

!       if (par_elec%accel(3) < 0.0d0) then
!         new_val = 1.0d0-12
!       end if

!       alpha = new_val / old_val
      !print *, alpha
      !if (alpha < 1.0d0) then
      !  alpha = alpha * 0.05d0
      !end if
      CALL RANDOM_NUMBER(rnd)
      if (rnd < alpha) then
        old_val = new_val
        old_val_s = new_val_s
        par_pos = new_pos
        !write(ud_metro, '(i5, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8)', iostat=IFAIL) i, par_pos(1)/length_scale, par_pos(2)/length_scale, par_pos(3)/length_scale, alpha
      end if

      !write(ud_metro, '(i5, tr2, E16.8, tr2, E16.8, tr2, E16.8)', iostat=IFAIL) i, par_pos(1)/length_scale, par_pos(2)/length_scale, par_pos(3)/length_scale
    end do

    !close(unit=ud_metro, iostat=IFAIL, status='keep')

    !stop

    !ndim = new_val_s
    Metro_algo_tip = par_pos
  end function Metro_algo_tip

  subroutine Check_xi(xi)
    double precision, intent(inout) :: xi
    double precision             :: d_xi, max_d

    max_d = max_xi - 1.0d0

    if (xi > max_xi) then
      d_xi = xi - max_xi
      if (d_xi > max_d) then
        xi = max_xi
        !print *, 'max_d'
      else
        xi = max_xi - d_xi
      end if
    end if

    if (xi < 0.0d0) then
      xi = 1.0d0
    end if

    if (xi < 1.0d0) then
      d_xi = 1.0 - xi
      if (d_xi > max_d) then
        xi = 1.0d0
      else
        xi = 1.0d0 + d_xi
      end if
    end if
  end subroutine Check_xi

  ! ----------------------------------------------------------------------------
  ! Checks the boundary conditions of the box.
  ! Check which particles to remove
  ! Enforce periodic boundary conditions (ToDo)
  subroutine Check_Boundary_ElecHole_Tip(i)
    integer, intent(in) :: i
    double precision    :: x, y, z, eta

    x = particles_cur_pos(1, i)
    y = particles_cur_pos(2, i)
    z = particles_cur_pos(3, i)
    eta = eta_coor(x, y, z)

    ! Check if the particle should be removed from the system
    if (z < 0.0d0) then
        call Mark_Particles_Remove(i, remove_bot)
    else if (z > box_dim(3)) then
        call Mark_Particles_Remove(i, remove_top)
    end if

    if (eta < eta_1) then
      call Mark_Particles_Remove(i, remove_bot)
    end if

  end subroutine Check_Boundary_ElecHole_Tip

!----------------------------------------------------------------------------------------
! The functions v_y and t_y are because of the image charge effect in the FN equation.
! The approximation for v_y and t_y are taken from
! Forbes, R. G., & Deane, J. H. (2007, November).
! "Reformulation of the standard theory of Fowler–Nordheim tunnelling and cold field electron emission."
! In Proceedings of the Royal Society of London A: Mathematical,
! Physical and Engineering Sciences (Vol. 463, No. 2087, pp. 2907-2927). The Royal Society.
!
  double precision function v_y(F, pos)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos
    double precision                             :: l

    if (image_charge .eqv. .true.) then
      l = l_const * (-1.0d0*F) / w_theta_xy(pos)**2 ! l = y^2, y = 3.79E-4 * sqrt(F_old) / w_theta
      if (l > 1.0d0) then
        l = 1.0d0
      end if
      v_y = 1.0d0 - l + 1.0d0/6.0d0 * l * log(l)
    else
      v_y = 1.0d0
    end if
  end function v_y

  double precision function t_y(F, pos)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos
    double precision                             :: l

    if (image_charge .eqv. .true.) then
      l = l_const * (-1.0d0*F) / w_theta_xy(pos)**2 ! l = y^2, y = 3.79E-4 * sqrt(F_old) / w_theta
      if (l > 1.0d0) then
        print *, 'Error: l > 1.0'
        print *, 'l = ', l, ', F = ', F, ', t_y = ', t_y
        print *, 'x = ', pos(1), 'y = ', pos(2)
        l = 1.0d0
        !call Write_Current_Position()
        !stop
      end if

      t_y = 1.0d0 + l*(1.0d0/9.0d0 - 1.0d0/18.0d0*log(l))
    else
      t_y = 1.0d0
    end if
  end function t_y

  !-----------------------------------------------------------------------------
  ! The Fowler-Nordheim equation is
  ! J = a_FN*F^2/(t_y^2*w_theta)*exp(-b_FN*w_theta^(3/2)*v_y/F)
  ! We break it into two parts
  ! J = Elec_Supply*Escape_Prob .
  ! Elec_supply is the the part before the exponental
  ! and Escape_prob is the exponental.

  ! The electron supply part
  ! This functions returns the number of electrons
  ! Elec_supply = a_FN*F^2/(t_y^2*w_theta) * (A*time_step/q_0),
  ! To get the number of electrons we multiply the first part of the FN
  ! equation with the area (A) and the time step (time_step). This gives
  ! the current. The divide that with the charge of the electron (q_0) to get
  ! the number of electrons.
  double precision function Elec_Supply(A, F, pos)
    double precision, intent(in)                 :: A, F
    double precision, dimension(1:3), intent(in) :: pos
    double precision                             :: n

    n = A * a_FN * F**2 * time_step / (q_0 * w_theta_xy(pos) * (t_y(F, pos))**2)

    Elec_supply = n
  end function Elec_supply

  ! This function returns the escape probability of the Electrons.
  ! Escape_prob = exp(-b_FN*w_theta^(3/2)*v_y/F) .
  double precision function Escape_Prob(F, pos)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos

    Escape_Prob = exp(b_FN * (sqrt(w_theta_xy(pos)))**3 * v_y(F, pos) / (-1.0d0*F))

    if (Escape_Prob > 1.0d0) then
      print *, 'Escape_prob is larger than 1.0'
      print *, 'Escape_prob = ', Escape_prob
      print *, ''
    end if
    if (isnan(Escape_Prob)) then
      print *, 'b_FN ', b_FN
      print *, 'w_theta_xy ', w_theta_xy(pos)
      print *, 'pos ', pos
      print *, 'v_y ', v_y(F, pos)
      print *, 'F ', F
    end if
    print *, 'Escape_Prob ', Escape_Prob
    print *, isnan(Escape_Prob)
    print *, ieee_is_nan(Escape_Prob)
    print *, '-'
  end function Escape_Prob

  double precision function w_theta_xy(pos)
    double precision, dimension(1:3), intent(in) :: pos

    w_theta_xy = w_theta

    !if (pos(1) > pos(2)) then
    !  w_theta_xy = 2.0d0
    !else
    !  w_theta_xy = 2.4d0
    !end if
  end function w_theta_xy

end module mod_field_emission_tip
