!-------------------------------------------!
! Module for emission from tip              !
! Kristinn Torfason                         !
! 04.06.18                                  !
!-------------------------------------------!

Module mod_emission_tip
  use mod_global
  use mod_hyperboloid_tip
  use mod_verlet
  use mod_pair
  use mod_ic
  !use ieee_arithmetic
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

  double precision :: time_step_div_q0

  ! Use image Charge or not
  !logical, parameter          :: image_charge = .true.

  ! Gauss photo emission
  logical                                     :: EmitGauss = .True.
  integer                                     :: maxElecEmit = -1

  ! Photo emission
  integer, parameter                          :: MAX_EMISSION_TRY_PHOTO = 100

contains

  !-----------------------------------------------------------------------------
  ! Initialize the Field Emission
  subroutine Init_Emission_Tip()

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_ElecHole_Tip

    ! Function for the electric field in the system
    ptr_field_E => field_E_Hyperboloid

    ! The function that does the emission
    ptr_Do_Emission => Do_Emission_Tip
    !ptr_Do_Emission => Do_Simple_Field_Emission_tip
    !ptr_Do_Emission => Do_Field_Emission_Tip_Test

    ! The function to do image charge effects
    ptr_Image_Charge_effect => Sphere_IC_field
    !ptr_Image_Charge_effect => NULL()

    ! Parameters used in the module
    time_step_div_q0 = time_step / q_0
    res_s = 0.0d0


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

  end subroutine Init_Emission_Tip

  subroutine Clean_Up_Emission_Tip()
    ! Nothing to do here
  end subroutine Clean_Up_Emission_Tip

  subroutine Do_Emission_Tip(step)
    integer, intent(in) :: step
    integer             :: IFAIL

    posInit = 0
    nrEmitted = 0

    print *, 'Emission from tip'

    select case (emitters_type(1))
    case (1)
      ! Field emission
      call Do_Field_Emission_Tip_Test(step)

    case (2)
      ! Reverse the voltage
      call Do_Field_Emission_Tip_2(step)

    case (3)
      ! Photo emission
      call Do_Photo_Emission_Tip(step)
      print *, 'Photo'

    case default
      print *, 'Vacuum: ERROR unknown emitter type!!'
      stop
      print *, emitters_type(1)
    end select

    cur_time = time_step * step / time_scale ! Scaled in units of time_scale
    write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, posInit, &
    & nrEmitted, nrElec
  end subroutine Do_Emission_Tip

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
  call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1)
  nrElecEmit = nrElecEmit + 1

  posInit = posInit + nrElecEmit
  nrEmitted = nrEmitted + nrElecEmit
end subroutine Do_Field_Emission_Tip_2

subroutine Do_Field_Emission_Tip_Test(step)
  integer, intent(in)               :: step
  integer                           :: ndim, N_round, nrElecEmit
  integer                           :: i, s
  double precision, dimension(1:3)  :: par_pos, Field, par_vel, surf_norm
  double precision                  :: F, N_sup, xi, phi, D_f
  integer, parameter                :: MAX_EMISSION_TRY = 10000

  ! print *, 'Hi Test'
  ! print *, ''
  ! print *, 'Field at top'
  ! par_pos(1) = x_coor(1.0d0, eta_1, 0.0d0)
  ! par_pos(2) = y_coor(1.0d0, eta_1, 0.0d0)
  ! par_pos(3) = z_coor(1.0d0, eta_1, 0.0d0)
  ! Field = field_E_Hyperboloid(par_pos)
  ! F = Field_normal(par_pos, Field)
  ! print *, Field
  ! print *, F
  ! print *, 'Supply'
  ! print *, Elec_Supply_tip(F, par_pos)
  ! print *, 'Escape Prob'
  ! print *, Escape_Prob_Tip(F, par_pos)
  ! print *, ''
  ! print *, 'Field at bottom'
  ! par_pos(1) = x_coor(max_xi, eta_1, 0.0d0)
  ! par_pos(2) = y_coor(max_xi, eta_1, 0.0d0)
  ! par_pos(3) = z_coor(max_xi, eta_1, 0.0d0)
  ! Field = field_E_Hyperboloid(par_pos)
  ! F = Field_normal(par_pos, Field)
  ! print *, Field
  ! print *, F
  ! print *, 'Supply'
  ! print *, Elec_Supply_tip(F, par_pos)
  ! print *, 'Escape Prob'
  ! print *, Escape_Prob_Tip(F, par_pos)
  ! pause

  nrElecEmit = 0

  call Do_Cuba_Suave_FE_Tip_Test(1, N_sup)

  N_round = nint(N_sup + res_s)
  res_s = N_sup - N_round

  !print *, 'Number of electrons = ', N_round, N_sup
  !pause
  !!$OMP END SINGLE

  !!$OMP DO PRIVATE(i, s, ndim, par_pos, xi, phi, surf_norm) REDUCTION(nrElecEmit)
  do i = 1, N_round
    !print *, 'i = ', i
    do s = 1, MAX_EMISSION_TRY
      ndim = 30
      call Metro_algo_tip_v2(ndim, xi, phi, F, D_f, par_pos)

      if (F >= 0.0d0) then
        D_f = 0.0d0
      end if

      if (F < 0.0d0) then

        surf_norm = surface_normal(par_pos)
        par_pos = par_pos + surf_norm*length_scale

        !call Add_particle(par_pos, step)
        par_vel = 0.0d0
        !!$OMP CRITICAL
        call Add_Particle(par_pos, par_vel, species_elec, step, 1, -1)
        nrElecEmit = nrElecEmit + 1

        !!$OMP END CRITICAL
        exit
      end if
    end do
  end do
  !!$OMP END DO
  !!$OMP END PARALLEL

  !print *, 'Loop done'

  posInit = posInit + nrElecEmit
  nrEmitted = nrEmitted + nrElecEmit
end subroutine Do_Field_Emission_Tip_Test

!----------------------------------------------------------------------------------------

subroutine Do_Photo_Emission_Tip(step)
  integer, intent(in)              :: step
  integer                          :: emit, nrElecEmit, nrTry
  double precision                 :: r_pos
  double precision                 :: xi, phi
  !double precision                 :: x, y, z
  double precision, dimension(1:3) :: par_pos, par_vel
  double precision, dimension(1:3) :: field, surf_norm

  emit = 1
  nrElecEmit = 0

  ! Check if we are doing a Gaussian distributed emission
  ! and set the max number of electrons allowed to be emitted if we are
  if (EmitGauss .eqv. .TRUE.) then
    maxElecEmit = Gauss_Emission(step)
  end if

  do while (nrTry <= MAX_EMISSION_TRY_PHOTO)
    ! Check if we have reached the max number
    ! of electrons to be emitted
    if ((nrElecEmit >= maxElecEmit) .and. (maxElecEmit /= -1)) then
      exit
    end if


    if (nrElec == MAX_PARTICLES-1) then
      print *, 'WARNING: Reached maximum number of electrons!!!'
      exit
    end if

    CALL RANDOM_NUMBER(par_pos(1:2)) ! Gives a random number [0,1]
    ! Get random x and y coordinates
    par_pos(1:2) = par_pos(1:2)*R_base

    ! Check if position is inside tip
    r_pos = sqrt(par_pos(1)**2 + par_pos(2)**2)
    if (r_pos > R_base) then
      cycle ! Not within the emitter area. Try again
    end if

    ! Calculate the z-coordinate (see eq. 38 in sec. 3.3 in doc)
    par_pos(3) = eta_1/sqrt(1.0d0 - eta_1**2) * sqrt(par_pos(1)**2 + par_pos(2)**2 + a_foci**2*(1.0d0 - eta_1**2))

    nrTry = nrTry + 1
    !Check in plane
    field = Calc_Field_at(par_pos)

    if (field(3) < 0.0d0) then
      !1 nm Above plane
      surf_norm = surface_normal(par_pos)
      par_pos = par_pos + surf_norm*length_scale
      field = Calc_Field_at(par_pos)

      if (field(3) < 0.0d0) then

        ! Place 1 nm above plane
        par_vel = 0.0d0
        call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1)

        !print *, 'field = ', field
        !pause
        !call Add_Plane_Graph_emitt(par_pos, par_vel)

        nrElecEmit = nrElecEmit + 1
        nrTry = 0

      end if
    end if

  end do

  posInit = posInit + nrElecEmit
  nrEmitted = nrEmitted + nrElecEmit
end subroutine Do_Photo_Emission_Tip

!----------------------------------------------------------------------------------------

  subroutine Do_Field_Emission_Tip_1(step)
    integer, intent(in)              :: step
    double precision                 :: F
    double precision, dimension(1:3) :: par_pos, surf_norm, par_vel
    !double precision, dimension(1)   :: rnd
    double precision                 :: rnd, N_sup
    integer                          :: i, j, s, IFAIL, nrElecEmit, n_r
    double precision                 :: A_f, D_f, n_s, F_avg, n_add
    double precision                 :: len_phi, len_xi
    integer                          :: nr_phi, nr_xi, ndim, emit
    double precision                 :: xi_1, phi_1, xi_2, phi_2, xi_c, phi_c
    double precision, dimension(1:3) :: field

    !!!$OMP SINGLE
    emit = 1


    !print *, len_x, len_x / length
    !print *, len_y, len_y / length

    nrElecEmit = 0


    !!$OMP SINGLE
    !print *, 'F_avg = ', F_avg
    !F_avg = 0.0d0
    !write (ud_debug, "(i8, tr2, E16.8)", iostat=IFAIL) step, F_avg

    call Do_Cuba_Suave_FE_Tip(emit, N_sup)

    n_r = nint(N_sup + res_s)
    res_s = N_sup - N_r

    !print *, 'n_r = ', n_r
    if (n_r < 0) then
      print *, 'n_r < 0'
      print *, 'n_s = ', n_s
      print *, 'n_r = ', n_r
      stop
    end if

    !!!$OMP END SINGLE

    !print *, 'Doing emission'
    !print *, n_r

    !!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(s, ndim, par_pos, field, F, D_f, surf_norm, xi_1, phi_1, rnd, par_vel) SCHEDULE(GUIDED, 2500)
    do s = 1, n_r

      !!!$OMP FLUSH (particles, nrElec)
      ndim = 25
      call Metro_algo_tip_v2(ndim, xi_1, phi_1, F, D_f, par_pos)

      if (F >= 0.0d0) then
        D_f = 0.0d0
      end if

      CALL RANDOM_NUMBER(rnd)
      !print *, 'rnd ', rnd
      !print *, 'D_f ', D_f
      !print *, ''
      if (rnd <= D_f) then
        surf_norm = surface_normal(par_pos)
        par_pos = par_pos + surf_norm*length_scale
        !$OMP CRITICAL

        ! Add a particle to the system
        par_vel = 0.0d0
        call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1, 1)

        !nrElec = nrElec + 1
        nrElecEmit = nrElecEmit + 1
        !print *, 'Particle emitted'
        !$OMP END CRITICAL
      end if
    end do
    !!$OMP END PARALLEL DO

    !!!$OMP MASTER

    posInit = posInit + nrElecEmit
    nrEmitted = nrEmitted + nrElecEmit
    !print *, 'Emission done'
    !print *, posInit
    !print *, ''
    !!!$OMP END MASTER

    !!!$OMP END PARALLEL
  end subroutine Do_Field_Emission_Tip_1

  subroutine Do_Simple_Field_Emission_tip(step)
    integer, intent(in) :: step
    integer             :: i, j, s
    integer             :: nrElecEmit, ndim, n_addi, IFAIL
    integer             :: nr_phi, nr_xi
    double precision    :: len_phi, len_xi
    double precision    :: n_s, F_avg, A_f, F, n_add, A_tot
    double precision    :: xi_start, xi_end, xi_c
    double precision    :: phi_start, phi_end, phi_c
    integer, parameter  :: MAX_EMISSION_TRY = 10000
    double precision, parameter :: pre_fac_a = m_0 / (-1.0d0*q_0) ! Acceleration prefactor
    double precision, dimension(1:3) :: par_pos, surf_norm, par_accel, par_vel
    !type(electron)                   :: par_elec

    nr_phi = 50
    nr_xi = 50
    !len_phi = 2.0d0*pi / nr_phi
    !len_xi = (max_xi - 1.0d0) / nr_xi
    len_phi = 2.0d0*pi
    len_xi = max_xi - 1.0d0

    !print *, len_x, len_x / length
    !print *, len_y, len_y / length

    nrElecEmit = 0

    par_pos = 0.0d0
    n_s = 0.0d0
    F_avg = 0.0d0
    n_add = 0.0d0
    A_tot = 0.0d0

    !print *, 'len_xi = '

    !!$OMP PARALLEL DEFAULT(SHARED)

    !!$OMP DO PRIVATE(i, j, xi_1, phi_1, xi_2, phi_2)&
    !!$OMP& PRIVATE(xi_c, phi_c, par_pos, par_elec, F, A_f) REDUCTION(+:n_add)
    do i = 0, (nr_xi-1)
      xi_start = 1.0 + i*len_xi/nr_xi
      xi_end = 1.0 + (i+1)*len_xi/nr_xi
      xi_c = (xi_end + xi_start) * 0.5

      do j = 0, (nr_phi-1)

        phi_start = j*len_phi/nr_phi
        phi_end = (j+1)*len_phi/nr_phi
        phi_c = (phi_end + phi_start) * 0.5

        par_pos = xyz_corr(xi_c, eta_1, phi_c)


        par_accel = Calc_Field_at(par_pos)
        F = Field_normal(par_pos, par_accel)


        if (F < 0.0d0) then
          A_f = Tip_Area(xi_start, xi_end, phi_start, phi_end)
          A_tot = A_tot + A_f
          n_add = n_add + Elec_Supply(A_f, F, par_pos) * Escape_Prob_Tip(F, par_pos)
        end if

      end do
    end do
    !!$OMP END DO

    !!OMP SINGLE
    n_addi = nint(n_add)
    !print *, 'Number of electrons = ', n_addi, n_add
    !print *, 'Area = ', A_tot
    !pause
    !!$OMP END SINGLE

    !!$OMP DO PRIVATE(i, s, ndim, par_pos, xi_c, phi_c, surf_norm)
    do i = 1, n_addi
      do s = 1, MAX_EMISSION_TRY
        ndim = 30
        par_pos = Metro_algo_tip(ndim, xi_c, phi_c)

        if (ndim > 0) then
          !!$OMP CRITICAL

          surf_norm = surface_normal(par_pos)
          par_pos = par_pos + surf_norm*length_scale

          !call Add_particle(par_pos, step)
          par_vel = 0.0d0
          call Add_Particle(par_pos, par_vel, species_elec, step, 1, -1)
          nrElecEmit = nrElecEmit + 1

          !!$OMP END CRITICAL
          exit
        end if
      end do
    end do
    !!$OMP END DO
    !!$OMP END PARALLEL

    posInit = posInit + nrElecEmit
    nrEmitted = nrEmitted + nrElecEmit
  end subroutine Do_Simple_Field_Emission_tip

  subroutine Metro_algo_tip_v2(ndim, xi, phi, eta_f, df_cur, par_pos)
    integer, intent(in)              :: ndim
    double precision, intent(out)    :: xi, phi, eta_f, df_cur
    double precision                 :: new_xi, new_phi, new_eta_f, df_new, rnd, alpha
    double precision, dimension(1:3) :: std, new_pos, cur_pos, par_pos, field
    integer                          :: i, count

    std(1) = max_xi*0.075d0/100.d0 ! Standard deviation for the normal distribution is 0.075% of the emitter length.
    std(2) = 2.0d0*pi*0.075d0/100.d0
    ! This means that 68% of jumps are less than this value.
    ! The expected value of the absolute value of the normal distribution is std*sqrt(2/pi).

    ! Get a random initial position on the surface of the tip.
    ! We pick this location from a uniform distribution.
    count = 0
    do ! Infinite loop, we try to find a favourable position to start from
      CALL RANDOM_NUMBER(cur_pos(1:2))

      xi = (max_xi - 1.0d0)*cur_pos(1) + 1.0d0
      phi = 2.0d0*pi*cur_pos(2)

      cur_pos(1) = x_coor(xi, eta_1, phi)
      cur_pos(2) = y_coor(xi, eta_1, phi)
      cur_pos(3) = z_coor(xi, eta_1, phi)

      ! Calculate the electric field at this position
      field = Calc_Field_at(cur_pos)
      eta_f = Field_normal(cur_pos, field) ! Component normal to the surface

      if (eta_f < 0.0d0) then
        exit ! We found a nice spot so we exit the loop
      else
        count = count + 1
        if (count > 10000) exit ! The loop is infnite, must stop it at some point.
        ! In field emission it is rare the we reach the CL limit.
      end if
    end do

    ! Calculate the escape probability at this location
    if (eta_f < 0.0d0) then
      df_cur = Escape_Prob_tip(field(3), cur_pos)
    else
      df_cur = 0.0d0 ! Zero escape probabilty if field is not favourable
    end if

    !---------------------------------------------------------------------------
    ! We now pick a random distance and direction to jump to from our
    ! current location. We do this ndim times.
    do i = 1, ndim
      ! Find a new position using a normal distribution.
      new_pos(1:2) = box_muller((/0.0d0, 0.0d0/), std)
      new_xi = new_pos(1) + xi
      new_phi = new_pos(2) + phi

      ! Make sure that the new position is within the limits of the emitter area.
      call check_limits_metro_rec_tip(new_xi, new_phi)

      new_pos(1) = x_coor(new_xi, eta_1, new_phi)
      new_pos(2) = y_coor(new_xi, eta_1, new_phi)
      new_pos(3) = z_coor(new_xi, eta_1, new_phi)

      ! Calculate the field at the new position
      field = Calc_Field_at(new_pos)
      new_eta_f = Field_normal(new_pos, field)

      ! Check if the field is favourable for emission at the new position.
      ! If it is not then cycle, i.e. we reject this location and
      ! pick another one.
      if (new_eta_f > 0.0d0) cycle ! Do the next loop iteration, i.e. find a new position.

      ! Calculate the escape probability at the new position, to compair with
      ! the current position.
      df_new = Escape_Prob_Tip(new_eta_f, new_pos)

      ! If the escape probability is higher in the new location,
      ! then we jump to that location. If it is not then we jump to that
      ! location with the probabilty df_new / df_cur.
      if (df_new > df_cur) then
        cur_pos = new_pos ! New position becomes the current position
        df_cur = df_new
        eta_f = new_eta_f
        xi = new_xi
        phi = new_phi
      else
        alpha = df_new / df_cur

        CALL RANDOM_NUMBER(rnd)
        ! Jump to this position with probability alpha, i.e. if rnd is less than alpha
        if (rnd < alpha) then
          cur_pos = new_pos ! New position becomes the current position
          df_cur = df_new
          eta_f = new_eta_f
          xi = new_xi
          phi = new_phi
        end if
      end if
    end do

    par_pos = cur_pos
  end subroutine Metro_algo_tip_v2

  function Sphere_IC_field(pos_1, pos_2)
    double precision, dimension(1:3)             :: Sphere_IC_field
    double precision, dimension(1:3), intent(in) :: pos_1, pos_2
    double precision                 :: x_a, y_a, z_a
    double precision                 :: x_b, y_b, z_b
    double precision                 :: x, y, z
    double precision                 :: z_0, Sphere_R
    double precision                 :: dis_a, tmp_dis_a, tmp_dis_b

    if (image_charge .eqv. .true.) then

      z_0 = h_tip - r_tip
      Sphere_R = r_tip

      !print *, 'Sphere_R = ', Sphere_R
      !print *, 'z_0 = ', z_0

      x_a = pos_1(1)
      y_a = pos_1(2)
      z_a = pos_1(3)

      x = pos_2(1)
      y = pos_2(2)
      z = pos_2(3)

      dis_a = sqrt(x_a**2 + y_a**2 + (z_a - z_0)**2)

      z_b = z_0 + Sphere_R**2 / ( sqrt(1 + x_a**2/(z_a - z_0)**2 + y_a**2/(z_a - z_0)**2) * dis_a )
      x_b = (z_b - z_0) * x_a / (z_a - z_0)
      y_b = (z_b - z_0) * y_a / (z_a - z_0)

      tmp_dis_a = ( (x - x_a)**2 + (y - y_a)**2 + (z - z_a)**2 )**(3.0d0/2.0d0)
      tmp_dis_b = ( (x - x_b)**2 + (y - y_b)**2 + (z - z_b)**2 )**(3.0d0/2.0d0)

      Sphere_IC_field(1) = 1.0d0*q_0/(4.0d0*pi*epsilon_0) * ( (x_a - x)/tmp_dis_a - (Sphere_R*(x_b - x))/(dis_a*tmp_dis_b) )
      Sphere_IC_field(2) = 1.0d0*q_0/(4.0d0*pi*epsilon_0) * ( (y_a - y)/tmp_dis_a - (Sphere_R*(y_b - y))/(dis_a*tmp_dis_b) )
      Sphere_IC_field(3) = 1.0d0*q_0/(4.0d0*pi*epsilon_0) * ( (z_a - z)/tmp_dis_a - (Sphere_R*(z_b - z))/(dis_a*tmp_dis_b) )

    else
      Sphere_IC_field = 0.0d0
    end if

  end function Sphere_IC_field

  subroutine check_limits_metro_rec_tip(xi, phi)
    double precision, intent(inout) :: xi, phi
    double precision                :: d_xi, max_d

    ! Keep \phi between 0 and 2\pi
    if (phi > 2.0d0*pi) then
      phi = phi - 2.0d0*pi
    end if
    if (phi < 0.0d0) then
      phi = phi + 2.0d0*pi
    end if

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
  end subroutine check_limits_metro_rec_tip

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
      l = l_const * (-1.0d0*F) / w_theta_pos_tip(pos)**2 ! l = y^2, y = 3.79E-4 * sqrt(F_old) / w_theta
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
      l = l_const * (-1.0d0*F) / w_theta_pos_tip(pos)**2 ! l = y^2, y = 3.79E-4 * sqrt(F_old) / w_theta
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

    n = A * a_FN * F**2 * time_step / (q_0 * w_theta_pos_tip(pos) * (t_y(F, pos))**2)

    Elec_supply = n
  end function Elec_supply

  double precision function Elec_Supply_tip(F, pos)
  double precision, intent(in)                 :: F
  double precision, dimension(1:3), intent(in) :: pos

  !n = a_FN * F**2 * time_step / (q_0 * w_theta_pos_tip(pos) * (t_y(F, pos))**2)
  Elec_supply_tip = time_step_div_q0 * a_FN/(t_y(F, pos)**2*w_theta_pos_tip(pos)) * F**2
  !Elec_supply_tip = time_step_div_q0 * a_FN/(1.0d0*w_theta_pos_tip(pos)) * F**2

  !Elec_supply_tip = a_FN/(t_y(F, pos)**2*w_theta_pos_tip(pos)) * F**2 * Escape_Prob_Tip(F, pos)

end function Elec_supply_tip

  ! This function returns the escape probability of the Electrons.
  ! Escape_prob = exp(-b_FN*w_theta^(3/2)*v_y/F) .
  double precision function Escape_Prob_Tip(F, pos)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos

    Escape_Prob_Tip = exp(b_FN * (sqrt(w_theta_pos_tip(pos)))**3 * v_y(F, pos) / abs(F))

    if (Escape_Prob_Tip > 1.0d0) then
      print *, 'Escape_prob is larger than 1.0'
      print *, 'Escape_prob = ', Escape_prob_Tip
      print *, ''
    end if
    if (isnan(Escape_Prob_tip) .eqv. .True.) then
    print *, 'Escape_Prob++ ', Escape_Prob_Tip
      print *, 'b_FN ', b_FN
      print *, 'w_theta_pos_tip ', w_theta_pos_tip(pos)
      print *, 'pos ', pos
      print *, 'v_y ', v_y(F, pos)
      print *, 'F ', F
      print *, 'WTF'
      !pause
    end if
    !print *, 'Escape_Prob-- ', Escape_Prob_Tip
    !print *, isnan(Escape_Prob_Tip)
    !print *, ieee_is_nan(Escape_Prob_Tip)
    !print *, '--'
  end function Escape_Prob_Tip

  double precision function w_theta_pos_tip(pos)
    double precision, dimension(1:3), intent(in) :: pos

    w_theta_pos_tip = w_theta
  end function w_theta_pos_tip


  ! ----------------------------------------------------------------------------
  ! The integration function for the Cuba library
  !
  integer function integrand_cuba_fe_tip(ndim, xx, ncomp, ff, userdata)
    ! Input / output variables
  integer, intent(in) :: ndim ! Number of dimensions (Should be 2)
  integer, intent(in) :: ncomp ! Number of vector-components in the integrand (Always 1 here)
  integer, intent(in) :: userdata ! Additional data passed to the integral function (In our case the number of the emitter)
  double precision, intent(in), dimension(1:ndim)   :: xx ! Integration points, between 0 and 1
  double precision, intent(out), dimension(1:ncomp) :: ff ! Results of the integrand function

  ! Variables used for calculations
  double precision, dimension(1:3) :: par_pos, field
  double precision                 :: A ! Emitter area
  double precision                 :: eta_f ! Component normal to the surface
  double precision                 :: xi, phi ! Prolate coordinates
  double precision                 :: h_xi, h_phi ! Scale factors

  !! Emitter area
  !A = Tip_Area(1.0d0, max_xi, 0.0d0, 2.0d0*pi)

  ! Surface position
  ! Cuba does the intergration over the hybercube.
  ! It gives us coordinates between 0 and 1.
  xi = (max_xi - 1.0d0)*xx(1) + 1.0d0
  phi = 2.0d0*pi*xx(2)

  !par_pos(1) = x_coor(xi, eta_1, phi)
  !par_pos(2) = y_coor(xi, eta_1, phi)
  !par_pos(3) = z_coor(xi, eta_1, phi)
  par_pos = xyz_corr(xi, eta_1, phi)

  ! Calculate the electric field on the surface
  field = Calc_Field_at(par_pos)
  eta_f = Field_normal(par_pos, field)

  ! Add to the average field
  !F_avg = F_avg + field

  ! Check if the field is favourable for emission
  if (eta_f < 0.0d0) then
    ! The field is favourable for emission

    ! Calculate the scale factors
    h_xi = a_foci*sqrt((xi**2 - eta_1**2)/(xi**2 - 1.0d0))
    h_phi = a_foci*sqrt((xi**2 - 1.0d0)*(1 - eta_1**2))

    ! Calculate the current density at this point
    ff(1) = Elec_Supply_tip(eta_f, par_pos) * h_xi * h_phi
  else
    ! The field is NOT favourable for emission
    ! This point does not contribute
    ff(1) = 0.0d0
  end if

  ! We mutiply with 2.0*pi (max_xi - 1.0) because Cuba does the 
  ! integration over the hybercube, i.e. from 0 to 1.
  ff(1) = 2.0d0*pi*(max_xi - 1.0d0)*ff(1)
  
  integrand_cuba_fe_tip = 0 ! Return value to Cuba, 0 = success
  end function integrand_cuba_fe_tip

  ! ----------------------------------------------------------------------------
  ! The integration function for the Cuba library
  !
  integer function integrand_cuba_fe_tip_test(ndim, xx, ncomp, ff, userdata)
    ! Input / output variables
    integer, intent(in) :: ndim ! Number of dimensions (Should be 2)
    integer, intent(in) :: ncomp ! Number of vector-components in the integrand (Always 1 here)
    integer, intent(in) :: userdata ! Additional data passed to the integral function (In our case the number of the emitter)
    double precision, intent(in), dimension(1:ndim)   :: xx ! Integration points, between 0 and 1
    double precision, intent(out), dimension(1:ncomp) :: ff ! Results of the integrand function

    ! Variables used for calculations
    double precision, dimension(1:3) :: par_pos, field
    double precision                 :: A ! Emitter area
    double precision                 :: eta_f ! Component normal to the surface
    double precision                 :: xi, phi ! Prolate coordinates
    double precision                 :: h_xi, h_phi ! Scale factors

    !! Emitter area
    !A = Tip_Area(1.0d0, max_xi, 0.0d0, 2.0d0*pi)

    ! Surface position
    ! Cuba does the intergration over the hybercube.
    ! It gives us coordinates between 0 and 1.
    xi = (max_xi - 1.0d0)*xx(1) + 1.0d0
    phi = 2.0d0*pi*xx(2)

    !par_pos(1) = x_coor(xi, eta_1, phi)
    !par_pos(2) = y_coor(xi, eta_1, phi)
    !par_pos(3) = z_coor(xi, eta_1, phi)
    par_pos = xyz_corr(xi, eta_1, phi)

    ! Calculate the electric field on the surface
    field = Calc_Field_at(par_pos)
    eta_f = Field_normal(par_pos, field)

    ! Add to the average field
    !F_avg = F_avg + field

    ! Check if the field is favourable for emission
    if (eta_f < 0.0d0) then
      ! The field is favourable for emission

      ! Calculate the scale factors
      h_xi = a_foci*sqrt((xi**2 - eta_1**2)/(xi**2 - 1.0d0))
      h_phi = a_foci*sqrt((xi**2 - 1.0d0)*(1 - eta_1**2))

      ! Calculate the current density at this point
      ff(1) = Elec_Supply_tip(eta_f, par_pos)*Escape_Prob_Tip(eta_f, par_pos) * h_xi * h_phi
    else
      ! The field is NOT favourable for emission
      ! This point does not contribute
      ff(1) = 0.0d0
    end if

    ! We mutiply with 2.0*pi (max_xi - 1.0) because Cuba does the 
    ! integration over the hybercube, i.e. from 0 to 1.
    ff(1) = 2.0d0*pi*(max_xi - 1.0d0)*ff(1)
    
    integrand_cuba_fe_tip_test = 0 ! Return value to Cuba, 0 = success
  end function integrand_cuba_fe_tip_test

  subroutine Do_Cuba_Suave_FE_Tip(emit, N_sup)
    ! Input / output variables
    integer, intent(in)           :: emit
    double precision, intent(out) :: N_sup

    ! Cuba integration variables
    integer, parameter :: ndim = 2 ! Number of dimensions
    integer, parameter :: ncomp = 1 ! Number of components in the integrand
    integer            :: userdata = 0 ! User data passed to the integrand
    integer, parameter :: nvec = 1 ! Number of points given to the integrand function
    double precision   :: epsrel = 1.0d-2 ! Requested relative error
    double precision   :: epsabs = 1.0d-4 ! Requested absolute error
    integer            :: flags = 0+4 ! Flags
    integer            :: seed = 0 ! Seed for the rng. Zero will use Sobol.
    integer            :: mineval = 10000 ! Minimum number of integrand evaluations
    integer            :: maxeval = 5000000 ! Maximum number of integrand evaluations
    integer            :: nnew = 2500 ! Number of integrand evaluations in each subdivision
    integer            :: nmin = 1000 ! Minimum number of samples a former pass must contribute to a subregion to be considered in the region's compound integral value.
    double precision   :: flatness = 5.0d0 ! Determine how prominently out-liers, i.e. samples with a large fluctuation, 
                                           ! figure in the total fluctuation, which in turn determines how a region is split up.
                                           ! As suggested by its name, flatness should be chosen large for 'flat" integrand and small for 'volatile' integrands
                                           ! with high peaks.
    character          :: statefile = "" ! File to save the state in. Empty string means don't do it.
    integer            :: spin = -1 ! Spinning cores
    integer            :: nregions ! <out> The actual number of subregions nedded
    integer            :: neval ! <out> The actual number of integrand evaluations needed
    integer            :: fail ! <out> Error flag (0 = Success, -1 = Dimension out of range, >0 = Accuracy goal was not met)
    double precision, dimension(1:ncomp) :: integral ! <out> The integral of the integrand over the unit hybercube
    double precision, dimension(1:ncomp) :: error ! <out> The presumed absolute error
    double precision, dimension(1:ncomp) :: prob ! <out> The chi-square probability


    ! Initialize the average field to zero
    !F_avg = 0.0d0

    ! Pass the number of the emitter being integraded over to the integrand as userdata
    userdata = emit

    call suave(ndim, ncomp, integrand_cuba_fe_tip, userdata, nvec, &
     & epsrel, epsabs, flags, seed, &
     & mineval, maxeval, nnew, nmin, flatness, &
     & statefile, spin, &
     & nregions, neval, fail, integral, error, prob)

     if (fail /= 0) then
      print '(a)', 'Vacuum: WARNING Cuba did not return 0'
      print *, fail
      print *, error
      print *, prob
      print *, integral(1)
     end if

     ! Round the results to the nearest integer
     !N_sup = nint( integral(1) )
     N_sup = integral(1)

    !  print *, 'Integral results'
    !  print *, integral(1)
    !  print *, N_sup
    !  print *, Tip_Area(1.0d0, max_xi, 0.0d0, 2.0d0*pi)
    !  print *, ''
    !  print *, a_foci
    !  pause

     ! Finish calculating the average field
     !F_avg = F_avg / neval
  end subroutine Do_Cuba_Suave_FE_Tip

  subroutine Do_Cuba_Suave_FE_Tip_Test(emit, N_sup)
    ! Input / output variables
    integer, intent(in)           :: emit
    double precision, intent(out) :: N_sup

    ! Cuba integration variables
    integer, parameter :: ndim = 2 ! Number of dimensions
    integer, parameter :: ncomp = 1 ! Number of components in the integrand
    integer            :: userdata = 0 ! User data passed to the integrand
    integer, parameter :: nvec = 1 ! Number of points given to the integrand function
    double precision   :: epsrel = 1.0d-2 ! Requested relative error
    double precision   :: epsabs = 1.0d-4 ! Requested absolute error
    integer            :: flags = 0+4 ! Flags
    integer            :: seed = 0 ! Seed for the rng. Zero will use Sobol.
    integer            :: mineval = 1000 ! Minimum number of integrand evaluations
    integer            :: maxeval = 50000000 ! Maximum number of integrand evaluations
    integer            :: nnew = 250 ! Number of integrand evaluations in each subdivision
    integer            :: nmin = 100 ! Minimum number of samples a former pass must contribute to a subregion to be considered in the region's compound integral value.
    double precision   :: flatness = 5.0d0 ! Determine how prominently out-liers, i.e. samples with a large fluctuation, 
                                           ! figure in the total fluctuation, which in turn determines how a region is split up.
                                           ! As suggested by its name, flatness should be chosen large for 'flat" integrand and small for 'volatile' integrands
                                           ! with high peaks.
    character          :: statefile = "" ! File to save the state in. Empty string means don't do it.
    integer            :: spin = -1 ! Spinning cores
    integer            :: nregions ! <out> The actual number of subregions nedded
    integer            :: neval ! <out> The actual number of integrand evaluations needed
    integer            :: fail ! <out> Error flag (0 = Success, -1 = Dimension out of range, >0 = Accuracy goal was not met)
    double precision, dimension(1:ncomp) :: integral ! <out> The integral of the integrand over the unit hybercube
    double precision, dimension(1:ncomp) :: error ! <out> The presumed absolute error
    double precision, dimension(1:ncomp) :: prob ! <out> The chi-square probability

    !integer            :: verbose = 0
    !integer, parameter :: last = 4
    integer, parameter :: key = 0


    ! Initialize the average field to zero
    !F_avg = 0.0d0

    ! Pass the number of the emitter being integraded over to the integrand as userdata
    userdata = emit

    !call suave(ndim, ncomp, integrand_cuba_fe_tip_test, userdata, nvec, &
    ! & epsrel, epsabs, flags, seed, &
    ! & mineval, maxeval, nnew, nmin, flatness, &
    ! & statefile, spin, &
    ! & nregions, neval, fail, integral, error, prob)

     call cuhre(ndim, ncomp, integrand_cuba_fe_tip_test, userdata, nvec, &
     & epsrel, epsabs, flags, &
     & mineval, maxeval, key, &
     & statefile, spin, &
     & nregions, neval, fail, integral, error, prob)


     if (fail /= 0) then
      print '(a)', 'Vacuum: WARNING Cuba did not return 0'
      print *, fail
      print *, error
      print *, prob
      print *, integral(1)
     end if

     ! Round the results to the nearest integer
     !N_sup = nint( integral(1) )
     N_sup = integral(1)

    !  print *, 'Integral results'
    !  print *, integral(1)
    !  print *, N_sup
    !  print *, Tip_Area(1.0d0, max_xi, 0.0d0, 2.0d0*pi)
    !  print *, ''
    !  print *, a_foci
    !  pause

     ! Finish calculating the average field
     !F_avg = F_avg / neval
  end subroutine Do_Cuba_Suave_FE_Tip_Test

  ! Gives a gaussian emission curve
  ! where step is the current time step
  ! returns the number of electrons allowed to be emitted in that time step
  integer function Gauss_Emission(step)
    integer, intent(in)         :: step ! Current time step
    integer                     :: IFAIL
    double precision, parameter :: sigma = 1000.0d0 ! Width / standard deviation
    double precision, parameter :: mu = 3000.0d0 ! Center
    double precision, parameter :: A = 6.0d0 ! Height
    double precision, parameter :: b = 1.0d0/(2.0d0*pi*sigma**2)

    Gauss_Emission = IDNINT(  A * exp( -1.0d0*b*(step - mu)**2 )  )

    write (ud_gauss, "(i6, tr2, i6)", iostat=IFAIL) step, Gauss_Emission
  end function Gauss_Emission
end module mod_emission_tip
