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
  use mod_kevin_rjgtf_v2
  use mod_cuba_integration
  !use ieee_arithmetic
  implicit none

  PRIVATE
  PUBLIC :: Init_Emission_Tip, Clean_Up_Emission_Tip

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

  ! MH algorithm parameters
  double precision :: a_rate = 0.5d0, MH_std = 1.0d0 ! Acceptance rate and standard deviation for the MH algorithm
  integer          :: jump_a = 0, jump_r = 0 ! Number of jumps accepted and rejected

  ! Tuning constants of the v3 GTF Metropolis-Hastings sampler
  ! (Metro_algo_tip_gtf_v3). The shared step MH_std is fixed within a chain
  ! and gets one clamped multiplicative update per chain, towards the
  ! acceptance rate MH_target_tip. Same scheme as the planar samplers in
  ! mod_field_emission_v2.
  double precision, parameter :: MH_warmup_tip  = 0.25d0   ! Fraction of the jumps done with the fixed initial step
  double precision, parameter :: MH_init_tip    = 0.10d0   ! Initial step as a fraction of the (xi, phi) ranges
  double precision, parameter :: MH_target_tip  = 0.35d0   ! Acceptance rate the MH_std update aims for (optimum for a 2D random walk)
  double precision, parameter :: MH_gain_tip    = 0.025d0  ! Gain of the per-chain MH_std update
  double precision, parameter :: MH_std_max_tip = 0.125d0  ! Upper limit on MH_std
  double precision, parameter :: MH_std_min_tip = 0.0005d0 ! Lower limit on MH_std

  double precision, dimension(1:3)   :: F_avg = 0.0d0

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
    ptr_Check_Boundary => Check_Boundary_Tip

    ! Function for the electric field in the system
    ptr_field_E => field_E_Hyperboloid

    ! The function that does the emission
    ptr_Do_Emission => Do_Emission_Tip
    !ptr_Do_Emission => Do_Simple_Field_Emission_tip
    !ptr_Do_Emission => Do_Field_Emission_Tip_Test

    ! The function to do image charge effects
    ptr_Image_Charge_effect => Sphere_IC_field
    !ptr_Image_Charge_effect => NULL()

    ptr_E_zunit => E_zunit_tip

    ! Parameters used in the module
    time_step_div_q0 = time_step / q_0
    res_s = 0.0d0


    ! Tip Parameters

    d_tip = emitters_dim(1, 1)
    R_base = emitters_dim(2, 1)
    h_tip = emitters_dim(3, 1)

    d = d_tip + h_tip
    !E_z = -1.0d0*V_s/d

    ! Other Parameters
    max_xi = h_tip/d_tip + 1.0d0
    a_foci = sqrt(d_tip**2*R_base**2 / (h_tip**2 + 2*d_tip*h_tip) + d_tip**2)
    eta_1 = -1.0d0 * d_tip / a_foci

    theta_tip = acos(d_tip/a_foci)
    r_tip = a_foci * sin(theta_tip) * tan(theta_tip)
    shift_z = abs(a_foci * eta_1 * max_xi)

    pre_fac_E_tip = log( (1.0d0 + eta_1)/(1.0d0 - eta_1) * (1.0d0 - eta_2)/(1.0d0 + eta_2) )
    pre_fac_E_tip_unit_voltage = 2.0d0 * 1.0d0 / (a_foci * pre_fac_E_tip)
    pre_fac_E_tip = 2.0d0 * V_s / (a_foci * pre_fac_E_tip)

  end subroutine Init_Emission_Tip

  subroutine Clean_Up_Emission_Tip()
    ! Nothing to do here
  end subroutine Clean_Up_Emission_Tip

  function E_zunit_tip(pos)
    double precision, dimension(1:3), intent(in) :: pos
    double precision, dimension(1:3)             :: E_zunit_tip

    ! Calculate the electric field at the tip
    E_zunit_tip = field_E_Hyperboloid(pos)
    E_zunit_tip = E_zunit_tip * pre_fac_E_tip_unit_voltage / pre_fac_E_tip 

  end function E_zunit_tip

  subroutine Do_Emission_Tip(step)
    integer, intent(in) :: step
    integer             :: IFAIL

    posInit = 0
    nrEmitted = 0

    select case (emitters_type(1))
    case (1)
      ! Field emission
      !call Do_Field_Emission_Tip(step)
      call Do_Field_Emission_Tip_OLDCODE(step)

    case (2)
      ! Reverse the voltage
      call Do_Field_Emission_Tip_2(step)

    case (3)
      ! Photo emission
      call Do_Photo_Emission_Tip(step)

    case (4)
      ! GTF emission
      call Do_GTF_Emission_Tip(step)

    case default
      print *, 'RUMDEED: ERROR unknown emitter type!!'
      stop
      print *, emitters_type(1)
    end select

    cur_time = time_step * step / time_scale ! Scaled in units of time_scale
    write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, nrEmitted, nrElec
  end subroutine Do_Emission_Tip

!----------------------------------------------------------------------------------------
! GTF emission
subroutine Do_GTF_Emission_Tip(step)
  integer, intent(in)               :: step
  integer                           :: nrElecEmit, N_round, ndim, i, IFAIL
  double precision                  :: N_sup, xi, phi, F, D_f
  double precision, dimension(1:3)  :: par_pos, par_vel, surf_norm

  nrElecEmit = 0

  call Do_Surface_Integration_GTF_Tip(1, N_sup)

  N_round = Rand_Poisson(N_sup)
  !res_s = N_sup - N_round

  do i = 1, N_round
    ndim = 80 ! Number of jumps in the Metropolis-Hastings chain

    call Metro_algo_tip_gtf_v3(ndim, xi, phi, F, D_f, par_pos)

    if (F < 0.0d0) then
      surf_norm = surface_normal(par_pos)
      par_pos = par_pos + surf_norm*length_scale ! Move 1 nm above the surface
      par_vel = 0.0d0 ! No initial velocity

      ! Add the particle to the system
      call Add_Particle(par_pos, par_vel, species_elec, step, 1, -1)
      nrElecEmit = nrElecEmit + 1
    end if
  end do
  
  posInit = posInit + nrElecEmit
  nrEmitted = nrEmitted + nrElecEmit

  write (ud_field, "(i8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8)", iostat=IFAIL) &
                                      step, F_avg(1), F_avg(2), F_avg(3), N_sup, a_rate, MH_std
end subroutine Do_GTF_Emission_Tip

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

subroutine Do_Field_Emission_Tip(step)
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

  call Do_Surface_Integration_FE_Tip(1, N_sup)

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
end subroutine Do_Field_Emission_Tip

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
  nrTry = 0

  ! Check if we are doing a Gaussian distributed emission
  ! and set the max number of electrons allowed to be emitted if we are.
  ! We treat the number from the gauss distribution as a mean number of a
  ! poisson random number generator.
  if (EmitGauss .eqv. .TRUE.) then
    maxElecEmit = Rand_Poisson( Gauss_Emission(step) )
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
    par_pos(3) = eta_1/sqrt(1.0d0 - eta_1**2) * sqrt(par_pos(1)**2 + par_pos(2)**2 + a_foci**2*(1.0d0 - eta_1**2)) + shift_z

    nrTry = nrTry + 1
    !Check in plane
    field = Calc_Field_at(par_pos)
    !print *, par_pos/length_scale

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


  subroutine Do_Field_Emission_Tip_OLDCODE(step)
    integer, intent(in)              :: step
    double precision                 :: F
    double precision, dimension(1:3) :: par_pos, surf_norm, par_vel
    !double precision, dimension(1)   :: rnd
    double precision, allocatable, dimension(:) :: rnd
    integer                          :: i, j, s, IFAIL, nrElecEmit, n_r
    double precision                 :: A_f, D_f, n_s, F_avg, n_add
    double precision                 :: len_phi, len_xi
    integer                          :: nr_phi, nr_xi, ndim
    double precision                 :: xi_1, phi_1, xi_2, phi_2, xi_c, phi_c

    !!$OMP SINGLE

    nr_phi = 100
    nr_xi = nr_phi
    len_phi = 2.0d0*pi / nr_phi
    len_xi = (max_xi - 1.0d0) / nr_xi

    !print *, len_x, len_x / length
    !print *, len_y, len_y / length

    nrElecEmit = 0
    n_s = 0.0d0
    F_avg = 0.0d0

    !print *, 'Integration loop'
    do i = 1, nr_xi
      do j = 1, nr_phi
        xi_c = 1.0d0 + (i - 0.5d0)*len_xi
        phi_c = (j - 0.5d0)*len_phi

        par_pos(1) = x_coor(xi_c, eta_1, phi_c)
        par_pos(2) = y_coor(xi_c, eta_1, phi_c)
        par_pos(3) = z_coor(xi_c, eta_1, phi_c)

        F = Field_normal(par_pos, Calc_Field_at(par_pos))
        F_avg = F_avg + F
        !if (abs(F_avg) > 1.0d11) then
        !  print *, 'F = ', F
        !  print *, 'nrElec = ', nrElec
          !stop
        !end if

        if (F >= 0.0d0) then
          n_add = 0.0d0
        else
          xi_1 = 1.0d0 + (i - 1.0d0)*len_xi
          phi_1 = (j - 1.0d0)*len_phi
          xi_2 = 1.0d0 + (i + 0.0d0)*len_xi
          phi_2 = (j + 0.0d0)*len_phi
          A_f = Tip_Area(xi_1, xi_2, phi_1, phi_2)
          n_add = Elec_Supply(A_f, F, par_pos)
          !if (isnan(n_add) == .true.) then
          !if (n_add < 0.0d0) then
          !  print *, 'A_f = ', A_f
          !  print *, 'F = ', F
          !  print *, 'n_add = ', n_add
            !stop
          !end if
        end if

        n_s = n_s + n_add
      end do
    end do

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

    allocate(rnd(1:n_r))
    !IFAIL = vdrnguniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_r, rnd(1:n_r), 0.0d0, 1.0d0)
    call random_number(rnd)

    !print *, 'Emission loop'
    do s = 1, n_r

      ndim = 80 ! Number of jumps in the Metropolis-Hastings chain
      call Metro_algo_tip_v3(ndim, xi_1, phi_1, F, D_f, par_pos)
      ! A failed chain comes back with F = 1.0 and D_f = 0, so no electron
      ! is emitted from it.

      ! Emit the electron with the escape probability at the sampled position
      if ((F < 0.0d0) .and. (rnd(s) <= D_f)) then
        !par_pos(3) = par_pos(3) + 1.0d0*length
        surf_norm = surface_normal(par_pos)
        par_pos = par_pos + surf_norm*length_scale
        !!$OMP CRITICAL
          par_vel = 0.0d0
          call Add_Particle(par_pos, par_vel, species_elec, step, 1, -1)
          !nrElec = nrElec + 1
          nrElecEmit = nrElecEmit + 1
        !!$OMP END CRITICAL
      end if
    end do
    !!$OMP END DO

    !!$OMP MASTER
    deallocate(rnd)

    posInit = posInit + nrElecEmit
    nrEmitted = nrEmitted + nrElecEmit
    !!$OMP END MASTER

    !!!$OMP END PARALLEL        
  end subroutine Do_Field_Emission_Tip_OLDCODE

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
    do ! Infinite loop, we try to find a favorable position to start from
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
        if (count > 10000) exit ! The loop is infinite, must stop it at some point.
        ! In field emission it is rare the we reach the CL limit.
      end if
    end do

    ! Calculate the escape probability at this location
    if (eta_f < 0.0d0) then
      df_cur = Escape_Prob_tip(eta_f, cur_pos)
    else
      df_cur = 0.0d0 ! Zero escape probability if field is not favorable
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

      ! Check if the field is favorable for emission at the new position.
      ! If it is not then cycle, i.e. we reject this location and
      ! pick another one.
      if (new_eta_f > 0.0d0) cycle ! Do the next loop iteration, i.e. find a new position.

      ! Calculate the escape probability at the new position, to compar with
      ! the current position.
      df_new = Escape_Prob_Tip(new_eta_f, new_pos)

      ! If the escape probability is higher in the new location,
      ! then we jump to that location. If it is not then we jump to that
      ! location with the probability df_new / df_cur.
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

  !--------------------------------------------------------------------------------
  ! MH algorithm with Robbins-Monro tuning
  ! Not used at the moment, but can be used in the future.
  subroutine Metro_algo_tip_gtf_v2(ndim, xi, phi, eta_f, df_cur, par_pos)
    integer, intent(in)              :: ndim
    double precision, intent(out)    :: xi, phi, eta_f, df_cur
    double precision                 :: xi_log ! Logarithmic xi value
    double precision                 :: new_xi, new_phi, new_eta_f, df_new, rnd, alpha
    double precision                 :: new_xi_log ! Logarithmic xi value
    double precision                 :: sigma_xiphi = 1.0d0 ! Standard deviations for the normal distribution
    double precision                 :: gamma = 0.01 ! Tuning parameter for the Robbins-Monro algorithm
    double precision, dimension(1:3) :: std, new_pos, cur_pos, par_pos, field
    integer                          :: count, max_tune = 1000
    integer                          :: accepted = 0 ! Number of accepted jumps
    integer                          :: rejected = 0 ! Number of rejected jumps
    double precision                 :: accepted_rate
    logical                          :: done_tune = .false. ! Flag to indicate if we are tuning the parameters
    integer                          :: iter_after_tuning = 0 ! Number of iterations after tuning is done
    integer                          :: iter_before_tuning = 0 ! Number of iterations before tuning is done
    integer                          :: total_iter = 0 ! Total number of iterations

    ! Find a starting position on the tip
    count = 0
    do ! Infinite loop, we try to find a favourable position to start from
      CALL RANDOM_NUMBER(cur_pos(1:2))
      
      ! xi = 1 + (M - 1)*u^k, the k value skews the distribution towards the tip.
      xi = 1.0d0 + (max_xi - 1.0d0)*cur_pos(1)**2 ! xi in [1, max_xi] skewed towards the tip
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

    ! Using logarithmic xi values keeps xi between 1 and max_xi
    xi_log = log(xi) ! Logarithmic xi value

    ! Calculate the escape probability at this location
    if (eta_f < 0.0d0) then
      df_cur = Escape_Prob_tip(eta_f, cur_pos)
    else
      df_cur = 0.0d0 ! Zero escape probabilty if field is not favourable
    end if

    do
      ! New position using a normal distribution
      ! mean = 0 and std = 1
      new_pos(1:2) = box_muller((/0.0d0, 0.0d0/), (/1.0d0, 1.0d0/))

      new_xi_log = xi_log + sigma_xiphi*new_pos(1) ! Logarithmic xi value
      new_xi = exp(new_xi_log) ! Convert back to linear xi value
      ! Make sure that the new xi value is not greater than max_xi
      if (new_xi > max_xi) then
        print *, 'Warning: new_xi > max_xi, new_xi = ', new_xi
        sigma_xiphi = sigma_xiphi * 0.5d0 ! Reduce the standard deviation
        cycle ! Reject this jump, i.e. do not use this position
      end if
      

      new_phi = phi + sigma_xiphi*new_pos(2) ! New phi value
      ! Make sure phi is between 0 and 2*pi
      new_phi = mod(new_phi, 2.0d0*pi)

      new_pos(1) = x_coor(new_xi, eta_1, new_phi)
      new_pos(2) = y_coor(new_xi, eta_1, new_phi)
      new_pos(3) = z_coor(new_xi, eta_1, new_phi)
      ! Calculate the field at the new position
      field = Calc_Field_at(new_pos)
      new_eta_f = Field_normal(new_pos, field)
      ! Check if the field is favourable for emission at the new position.
      ! If it is not then cycle, i.e. we reject this location and
      ! pick another one.
      if (new_eta_f > 0.0d0) then
        print *, 'Warning: new_eta_f > 0.0d0, new_eta_f = ', new_eta_f
        cycle ! Do the next loop iteration, i.e. find a new position.
      end if
      ! Calculate the escape probability at the new position, to compair with
      df_new = Escape_Prob_tip(new_eta_f, new_pos)

      ! If the escape probability is higher in the new location,
      ! then we jump to that location. If it is not then we jump to that
      ! location with the probabilty df_new / df_cur.
      if (df_new > df_cur) then
        cur_pos = new_pos ! New position becomes the current position
        df_cur = df_new
        eta_f = new_eta_f
        xi_log = new_xi_log ! Update logarithmic xi value
        xi = new_xi ! Update linear xi value
        phi = new_phi ! Update phi value
        accepted = accepted + 1 ! Count this as an accepted jump
      else
        alpha = df_new / df_cur

        CALL RANDOM_NUMBER(rnd)
        ! Jump to this position with probability alpha, i.e. if rnd is less than alpha
        if (rnd < alpha) then
          cur_pos = new_pos ! New position becomes the current position
          df_cur = df_new
          eta_f = new_eta_f
          xi_log = new_xi_log ! Update logarithmic xi value
          xi = new_xi ! Update linear xi value
          phi = new_phi ! Update phi value
          accepted = accepted + 1 ! Count this as an accepted jump
        else
          rejected = rejected + 1 ! Count this as a rejected jump
        end if
      end if

      if (done_tune .eqv. .true.) then
        iter_after_tuning = iter_after_tuning + 1
        if (iter_after_tuning >= ndim) then
          exit
        end if
      end if

      ! Every 100 iterations we adjust the standard deviation
      if (mod(total_iter, 100) == 0) then
        accepted_rate = DBLE(accepted) / DBLE(accepted + rejected)
        ! Adjust the standard deviation based on the acceptance rate
        sigma_xiphi = sigma_xiphi * exp(gamma*((accepted_rate - 0.234d0)))
        print *, 'Iteration: ', total_iter, ' Accepted rate: ', accepted_rate, ' Alpha xiphi: ', sigma_xiphi


        if (abs(accepted_rate - 0.234d0) < 0.05d0) then
          done_tune = .true. ! If the acceptance rate is close to 0.234, we stop tuning
          print *, 'Tuning done after ', total_iter, ' iterations.'
          !iter_before_tuning = 0 ! Reset the counter for iterations
        end if
      end if

      iter_before_tuning = iter_before_tuning + 1
      if (iter_before_tuning >= max_tune) then
        print *, 'Warning: Maximum tuning iterations reached without convergence.'
        exit
      end if

      total_iter = total_iter + 1
    end do

  end subroutine Metro_algo_tip_gtf_v2

  subroutine Metro_algo_tip_gtf(ndim, xi, phi, eta_f, df_cur, par_pos)
    integer, intent(in)              :: ndim
    double precision, intent(out)    :: xi, phi, eta_f, df_cur
    double precision                 :: new_xi, new_phi, new_eta_f, df_new, rnd, alpha
    double precision, dimension(1:3) :: std, new_pos, cur_pos, par_pos, field
    integer                          :: i, count
    double precision                 :: gamma = 0.01 ! Tuning parameter for the Robbins-Monro algorithm


    ! Try to keep the acceptance ratio around 50% by
    MH_std = MH_std * exp(gamma*(0.50d0 - a_rate)) ! Adjust the standard deviation based on the acceptance rate
    ! Limits on how big or low the standard deviation can be.
    if (MH_std > 0.1250d0) then
      MH_std = 0.1250d0
    else if (MH_std < 0.0005d0) then
      MH_std = 0.0005d0
    end if

    ! Standard deviation for the normal distribution
    std(1) = max_xi*2.5d0/100.d0 * MH_std
    std(2) = 2.0d0*pi*2.5d0/100.d0 * MH_std
    ! This means that 68% of jumps are less than this value.
    ! The expected value of the absolute value of the normal distribution is std*sqrt(2/pi).

    ! Get a random initial position on the surface of the tip.
    ! We pick this location from a uniform distribution.
    count = 0
    do ! Infinite loop, we try to find a favourable position to start from
      CALL RANDOM_NUMBER(cur_pos(1:2))
      
      xi = 1.0d0 + (max_xi - 1.0d0)*cur_pos(1)**2.0d0 ! xi in [1, max_xi] skewed towards the tip
      phi = 2.0d0*pi*cur_pos(2)

      cur_pos = xyz_corr(xi, eta_1, phi) ! xyz coordinates of the position on the tip surface

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
      df_cur = Escape_Prob_tip(eta_f, cur_pos)
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

      new_pos = xyz_corr(new_xi, eta_1, new_phi) ! xyz coordinates of the new position on the tip surface

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
        jump_a = jump_a + 1 ! Count the accepted jumps
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

          jump_a = jump_a + 1 ! Count the accepted jumps
        else
          jump_r = jump_r + 1 ! Count the rejected jumps
        end if
      end if
    end do

    par_pos = cur_pos

    ! Acceptance rate
    a_rate = DBLE(jump_a) / DBLE(jump_r + jump_a)
  end subroutine Metro_algo_tip_gtf

  !--------------------------------------------------------------------------------
  ! Log of the target distribution of the v3 GTF Metropolis-Hastings sampler:
  ! the GTF current density times the h_xi*h_phi area element of the
  ! (xi, phi) parametrization of the tip surface, i.e. the same expression
  ! that integrand_cuba_gtf_tip integrates for the number of electrons to
  ! emit. h_xi*h_phi = a_foci**2 * sqrt((xi**2 - eta_1**2)*(1 - eta_1**2));
  ! the constant factors cancel in the Metropolis ratio and are left out,
  ! as is the time_step/q_0 factor of the integrand. The tiny() guard
  ! catches a current density that underflows to zero (Get_Kevin_Jgtf_v2
  ! also returns exactly zero for a negligible field). eta_f must be < 0.
  double precision function Tip_gtf_target_log(eta_f, xi)
    double precision, intent(in) :: eta_f, xi

    Tip_gtf_target_log = log(max(Get_Kevin_Jgtf_v2(eta_f, T_temp, w_theta), tiny(1.0d0))) &
                     & + 0.5d0*log(xi**2 - eta_1**2)
  end function Tip_gtf_target_log

  !--------------------------------------------------------------------------------
  ! Metropolis-Hastings sampling of the emission position for the GTF
  ! emission from the tip (third version, used by Do_GTF_Emission_Tip).
  !
  ! The chain moves in the (xi, phi) parametrization of the tip surface and
  ! samples positions proportional to J_gtf * h_xi * h_phi, the integrand
  ! that Do_Cuba_Suave_GTF_Tip integrates for the number of electrons to
  ! emit. Do_GTF_Emission_Tip emits every candidate, so the emitted
  ! positions follow the current density over the tip surface and the total
  ! matches the surface integral.
  ! (Metro_algo_tip_gtf above targets the Fowler-Nordheim escape
  ! probability instead of the GTF current density, and the sign of its
  ! step size feedback is inverted: a high acceptance rate shrinks the
  ! step, which raises the acceptance rate further, until the chain is
  ! frozen at its starting position. That is why its starting point had to
  ! be skewed towards the apex. It is kept for reference, like the
  ! Robbins-Monro attempt Metro_algo_tip_gtf_v2.)
  !
  ! The proposal is a Gaussian step in (xi, phi), reflected at the xi
  ! limits, with phi wrapping around the tip. The first MH_warmup_tip of
  ! the jumps use the large fixed step MH_init_tip so the chain finds the
  ! high current region from a uniform start (no skew needed). After that
  ! the adapted step MH_std is used, fixed within the chain and updated
  ! once at the end from the acceptance rate of this chain.
  subroutine Metro_algo_tip_gtf_v3(ndim, xi, phi, eta_f, df_cur, par_pos)
    integer, intent(in)                           :: ndim
    double precision, intent(out)                 :: xi, phi, eta_f, df_cur
    double precision, dimension(1:3), intent(out) :: par_pos

    double precision                 :: new_xi, new_phi, new_eta_f
    double precision                 :: sup_cur, sup_new ! Log of the chain target
    double precision                 :: rnd, alpha
    double precision, dimension(1:2) :: std, step
    double precision, dimension(1:3) :: cur_pos, new_pos, field
    integer                          :: i, count, ndim_first
    integer                          :: acc, rej ! Accepted/rejected jumps after the warmup

    ! Bring the shared step size into the limits before use. It starts at
    ! 1.0 and the older samplers above use other limits.
    if (MH_std > MH_std_max_tip) then
      MH_std = MH_std_max_tip
    else if (MH_std < MH_std_min_tip) then
      MH_std = MH_std_min_tip
    end if

    acc = 0
    rej = 0
    ndim_first = nint(ndim*MH_warmup_tip)

    ! Warmup step size, MH_init_tip of the (xi, phi) ranges
    std(1) = (max_xi - 1.0d0)*MH_init_tip
    std(2) = 2.0d0*pi*MH_init_tip

    ! Get a random initial position on the tip surface, uniform in (xi, phi)
    count = 0
    do ! Infinite loop, we try to find a favourable position to start from
      CALL RANDOM_NUMBER(step)
      xi = 1.0d0 + (max_xi - 1.0d0)*step(1)
      phi = 2.0d0*pi*step(2)

      cur_pos = xyz_corr(xi, eta_1, phi)

      ! Calculate the electric field at this position
      field = Calc_Field_at(cur_pos)
      eta_f = Field_normal(cur_pos, field) ! Component normal to the surface

      if (eta_f < 0.0d0) then
        exit ! We found a nice spot so we exit the loop
      else
        count = count + 1
        if (count > 10000) then ! The loop is infinite, must stop it at some point.
          ! Mark the chain as failed with defined outputs. The caller
          ! checks the sign of eta_f and emits no electron from this chain.
          xi = 1.0d0
          phi = 0.0d0
          par_pos = xyz_corr(xi, eta_1, phi)
          eta_f = 1.0d0
          df_cur = 0.0d0
          print *, 'Failed to find spot for emission on the tip'
          return
        end if
      end if
    end do

    sup_cur = Tip_gtf_target_log(eta_f, xi)

    !---------------------------------------------------------------------------
    ! We now pick a random distance and direction to jump to from our
    ! current location. We do this ndim times.
    do i = 1, ndim

      ! After the warmup the jumps use the adapted step size. MH_std is
      ! only updated between chains, so the kernel is fixed within a chain.
      if (i > ndim_first) then
        std(1) = (max_xi - 1.0d0)*MH_std
        std(2) = 2.0d0*pi*MH_std
      end if

      ! Find a new position using a normal distribution
      step = box_muller((/0.0d0, 0.0d0/), std)
      new_xi = xi + step(1)
      new_phi = modulo(phi + step(2), 2.0d0*pi) ! phi wraps around the tip

      ! Reflect xi at the limits of the surface. With the step limited to
      ! MH_std_max_tip a jump needs at most one reflection; a longer one
      ! (only possible in the far tail of the warmup step) is rejected.
      if (new_xi > max_xi) new_xi = 2.0d0*max_xi - new_xi
      if (new_xi < 1.0d0) new_xi = 2.0d0 - new_xi
      if ((new_xi < 1.0d0) .or. (new_xi > max_xi)) then
        if (i > ndim_first) rej = rej + 1
        cycle ! Do the next loop iteration, i.e. find a new position.
      end if

      new_pos = xyz_corr(new_xi, eta_1, new_phi)

      ! Calculate the field at the new position
      field = Calc_Field_at(new_pos)
      new_eta_f = Field_normal(new_pos, field)

      ! Check if the field is favourable for emission at the new position.
      ! If it is not the current density there is zero, i.e. we reject this
      ! location and pick another one.
      if (new_eta_f >= 0.0d0) then
        if (i > ndim_first) rej = rej + 1
        cycle ! Do the next loop iteration, i.e. find a new position.
      end if

      ! Calculate the target at the new position, to compare with the
      ! current position
      sup_new = Tip_gtf_target_log(new_eta_f, new_xi)

      alpha = sup_new - sup_cur ! Log of the Metropolis ratio

      if (sup_new >= sup_cur) then
        cur_pos = new_pos
        xi = new_xi
        phi = new_phi
        eta_f = new_eta_f
        sup_cur = sup_new
        if (i > ndim_first) acc = acc + 1
      else
        CALL RANDOM_NUMBER(rnd)
        if (log(rnd) <= alpha) then
          cur_pos = new_pos
          xi = new_xi
          phi = new_phi
          eta_f = new_eta_f
          sup_cur = sup_new
          if (i > ndim_first) acc = acc + 1
        else
          if (i > ndim_first) rej = rej + 1
        end if
      end if
    end do

    ! The acceptance rate of this chain drives one update of the shared
    ! step size. Note the sign: a high acceptance rate means the steps are
    ! too small, so the step grows.
    if (acc + rej > 0) then
      a_rate = DBLE(acc) / DBLE(acc + rej)
      MH_std = MH_std * exp(MH_gain_tip*(a_rate - MH_target_tip))
      ! Limits on how big or low the standard deviation can be
      if (MH_std > MH_std_max_tip) then
        MH_std = MH_std_max_tip
      else if (MH_std < MH_std_min_tip) then
        MH_std = MH_std_min_tip
      end if
    end if

    ! Return the current position. The escape probability is returned for
    ! diagnostics only, the caller emits every candidate.
    par_pos = cur_pos
    df_cur = Escape_Prob_Tip(eta_f, cur_pos)
  end subroutine Metro_algo_tip_gtf_v3

  !--------------------------------------------------------------------------------
  ! Log of the target distribution of the v3 field emission sampler: the
  ! Fowler-Nordheim electron supply times the h_xi*h_phi area element of the
  ! (xi, phi) parametrization, i.e. the same supply that
  ! Do_Field_Emission_Tip_OLDCODE integrates over the tip surface for the
  ! number of candidate electrons. The escape probability is NOT part of
  ! the target: the caller applies it as the emission test, so the emitted
  ! positions follow the current density supply * D (see the corresponding
  ! pairing in mod_field_emission_v2). Constant factors cancel in the
  ! Metropolis ratio. The tiny() guard catches a supply that underflows to
  ! zero. eta_f must be < 0.
  double precision function Tip_fe_target_log(eta_f, xi, pos)
    double precision, intent(in)                 :: eta_f, xi
    double precision, dimension(1:3), intent(in) :: pos

    Tip_fe_target_log = log(max(Elec_Supply_tip(eta_f, pos), tiny(1.0d0))) &
                    & + 0.5d0*log(xi**2 - eta_1**2)
  end function Tip_fe_target_log

  !--------------------------------------------------------------------------------
  ! Metropolis-Hastings sampling of the emission position for the field
  ! emission from the tip (third version, used by
  ! Do_Field_Emission_Tip_OLDCODE).
  !
  ! The chain moves in the (xi, phi) parametrization of the tip surface and
  ! samples positions proportional to the electron supply times the area
  ! element (Tip_fe_target_log above). The caller then emits each candidate
  ! with the escape probability D returned in df_cur, so the emitted
  ! positions follow the current density supply * D over the tip surface
  ! and the expected total matches the surface integral of the current.
  ! (Metro_algo_tip_v2 above targets the escape probability with a step of
  ! 0.075% of the surface, so its chains stay at their uniform starting
  ! point: the emitted positions followed D alone, without the supply
  ! weight and the area element. It is kept for reference.)
  !
  ! The mechanics are identical to Metro_algo_tip_gtf_v3: Gaussian steps in
  ! (xi, phi) reflected at the xi limits with phi wrapping, a warmup with
  ! the large fixed step MH_init_tip, then the adapted step MH_std, fixed
  ! within a chain and updated once per chain from its acceptance rate.
  subroutine Metro_algo_tip_v3(ndim, xi, phi, eta_f, df_cur, par_pos)
    integer, intent(in)                           :: ndim
    double precision, intent(out)                 :: xi, phi, eta_f, df_cur
    double precision, dimension(1:3), intent(out) :: par_pos

    double precision                 :: new_xi, new_phi, new_eta_f
    double precision                 :: sup_cur, sup_new ! Log of the chain target
    double precision                 :: rnd, alpha
    double precision, dimension(1:2) :: std, step
    double precision, dimension(1:3) :: cur_pos, new_pos, field
    integer                          :: i, count, ndim_first
    integer                          :: acc, rej ! Accepted/rejected jumps after the warmup

    ! Bring the shared step size into the limits before use. It starts at
    ! 1.0 and the older samplers above use other limits.
    if (MH_std > MH_std_max_tip) then
      MH_std = MH_std_max_tip
    else if (MH_std < MH_std_min_tip) then
      MH_std = MH_std_min_tip
    end if

    acc = 0
    rej = 0
    ndim_first = nint(ndim*MH_warmup_tip)

    ! Warmup step size, MH_init_tip of the (xi, phi) ranges
    std(1) = (max_xi - 1.0d0)*MH_init_tip
    std(2) = 2.0d0*pi*MH_init_tip

    ! Get a random initial position on the tip surface, uniform in (xi, phi)
    count = 0
    do ! Infinite loop, we try to find a favourable position to start from
      CALL RANDOM_NUMBER(step)
      xi = 1.0d0 + (max_xi - 1.0d0)*step(1)
      phi = 2.0d0*pi*step(2)

      cur_pos = xyz_corr(xi, eta_1, phi)

      ! Calculate the electric field at this position
      field = Calc_Field_at(cur_pos)
      eta_f = Field_normal(cur_pos, field) ! Component normal to the surface

      if (eta_f < 0.0d0) then
        exit ! We found a nice spot so we exit the loop
      else
        count = count + 1
        if (count > 10000) then ! The loop is infinite, must stop it at some point.
          ! Mark the chain as failed with defined outputs. Zero escape
          ! probability means the caller emits no electron from this chain.
          xi = 1.0d0
          phi = 0.0d0
          par_pos = xyz_corr(xi, eta_1, phi)
          eta_f = 1.0d0
          df_cur = 0.0d0
          print *, 'Failed to find spot for emission on the tip'
          return
        end if
      end if
    end do

    sup_cur = Tip_fe_target_log(eta_f, xi, cur_pos)

    !---------------------------------------------------------------------------
    ! We now pick a random distance and direction to jump to from our
    ! current location. We do this ndim times.
    do i = 1, ndim

      ! After the warmup the jumps use the adapted step size. MH_std is
      ! only updated between chains, so the kernel is fixed within a chain.
      if (i > ndim_first) then
        std(1) = (max_xi - 1.0d0)*MH_std
        std(2) = 2.0d0*pi*MH_std
      end if

      ! Find a new position using a normal distribution
      step = box_muller((/0.0d0, 0.0d0/), std)
      new_xi = xi + step(1)
      new_phi = modulo(phi + step(2), 2.0d0*pi) ! phi wraps around the tip

      ! Reflect xi at the limits of the surface. With the step limited to
      ! MH_std_max_tip a jump needs at most one reflection; a longer one
      ! (only possible in the far tail of the warmup step) is rejected.
      if (new_xi > max_xi) new_xi = 2.0d0*max_xi - new_xi
      if (new_xi < 1.0d0) new_xi = 2.0d0 - new_xi
      if ((new_xi < 1.0d0) .or. (new_xi > max_xi)) then
        if (i > ndim_first) rej = rej + 1
        cycle ! Do the next loop iteration, i.e. find a new position.
      end if

      new_pos = xyz_corr(new_xi, eta_1, new_phi)

      ! Calculate the field at the new position
      field = Calc_Field_at(new_pos)
      new_eta_f = Field_normal(new_pos, field)

      ! Check if the field is favourable for emission at the new position.
      ! If it is not the supply there is zero, i.e. we reject this location
      ! and pick another one.
      if (new_eta_f >= 0.0d0) then
        if (i > ndim_first) rej = rej + 1
        cycle ! Do the next loop iteration, i.e. find a new position.
      end if

      ! Calculate the target at the new position, to compare with the
      ! current position
      sup_new = Tip_fe_target_log(new_eta_f, new_xi, new_pos)

      alpha = sup_new - sup_cur ! Log of the Metropolis ratio

      if (sup_new >= sup_cur) then
        cur_pos = new_pos
        xi = new_xi
        phi = new_phi
        eta_f = new_eta_f
        sup_cur = sup_new
        if (i > ndim_first) acc = acc + 1
      else
        CALL RANDOM_NUMBER(rnd)
        if (log(rnd) <= alpha) then
          cur_pos = new_pos
          xi = new_xi
          phi = new_phi
          eta_f = new_eta_f
          sup_cur = sup_new
          if (i > ndim_first) acc = acc + 1
        else
          if (i > ndim_first) rej = rej + 1
        end if
      end if
    end do

    ! The acceptance rate of this chain drives one update of the shared
    ! step size. Note the sign: a high acceptance rate means the steps are
    ! too small, so the step grows.
    if (acc + rej > 0) then
      a_rate = DBLE(acc) / DBLE(acc + rej)
      MH_std = MH_std * exp(MH_gain_tip*(a_rate - MH_target_tip))
      ! Limits on how big or low the standard deviation can be
      if (MH_std > MH_std_max_tip) then
        MH_std = MH_std_max_tip
      else if (MH_std < MH_std_min_tip) then
        MH_std = MH_std_min_tip
      end if
    end if

    ! Return the current position and the escape probability there, which
    ! the caller uses for the emission test.
    par_pos = cur_pos
    df_cur = Escape_Prob_Tip(eta_f, cur_pos)
  end subroutine Metro_algo_tip_v3

  ! Sphere_IC_field has moved to mod_hyperboloid_tip, so that mod_verlet can
  ! reference it for the OpenACC geometry dispatch (this module uses
  ! mod_verlet, so it cannot be referenced the other way). It is still
  ! available here through use association.

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

    ! Keep \xi between 1 and max_xi
    max_d = max_xi - 1.0d0

    ! Reflect the position if it is outside the limits
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
      d_xi = 1.0d0 - xi
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

      field = Calc_Field_at(new_pos)
      eta_f = Field_normal(new_pos, field)
      !new_val = norm2(par_elec%accel)
      new_val = eta_f

      ! Check if this point is favorable for emission
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
      d_xi = 1.0d0 - xi
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
  subroutine Check_Boundary_Tip(i)
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

  end subroutine Check_Boundary_Tip

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
  ! Elec_supply is the the part before the exponential
  ! and Escape_prob is the exponential.

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
    !if (isnan(Escape_Prob_tip) .eqv. .True.) then
    !print *, 'Escape_Prob++ ', Escape_Prob_Tip
    !  print *, 'b_FN ', b_FN
    !  print *, 'w_theta_pos_tip ', w_theta_pos_tip(pos)
    !  print *, 'pos ', pos
    !  print *, 'v_y ', v_y(F, pos)
    !  print *, 'F ', F
    !  !pause
    !end if
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

  ! ----------------------------------------------------------------------------
  ! The integration function for the Cuba library
  !
  integer(c_int) function integrand_cuba_fe_tip_test(ndim, xx, ncomp, ff, userdata)
    ! Called from the C Cuba library: the integer arguments and the result
    ! must be integer(c_int) (see mod_cuba_integration)
    use, intrinsic :: iso_c_binding, only: c_int
    ! Input / output variables
    integer(c_int), intent(in) :: ndim ! Number of dimensions (Should be 2)
    integer(c_int), intent(in) :: ncomp ! Number of vector-components in the integrand (Always 1 here)
    integer(c_int), intent(in) :: userdata ! Additional data passed to the integral function (In our case the number of the emitter)
    double precision, intent(in), dimension(1:ndim)   :: xx ! Integration points, between 0 and 1
    double precision, intent(out), dimension(1:ncomp) :: ff ! Results of the integrand function

    ! Variables used for calculations
    double precision, dimension(1:3) :: par_pos, field
    double precision                 :: A ! Emitter area
    double precision                 :: eta_f ! Component normal to the surface
    double precision                 :: xi, phi ! Prolate coordinates
    double precision                 :: h_area ! Area element of the surface

    !! Emitter area
    !A = Tip_Area(1.0d0, max_xi, 0.0d0, 2.0d0*pi)

    ! Surface position
    ! Cuba does the integration over the hybercube.
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
    F_avg = F_avg + field

    ! Check if the field is favorable for emission
    if (eta_f < 0.0d0) then
      ! The field is favorable for emission

      ! Area element of the (xi, phi) parametrization. The scale factors are
      ! each singular at the apex (xi = 1): h_xi diverges and h_phi vanishes,
      ! so forming them separately gives Inf*0 = NaN there. Their product is
      ! finite and has to be formed analytically, because Divonne evaluates
      ! the boundary of the hypercube (border = 0), i.e. xi = 1 exactly.
      h_area = a_foci**2 * sqrt((xi**2 - eta_1**2)*(1 - eta_1**2))

      ! Calculate the current density at this point
      ff(1) = Elec_Supply_tip(eta_f, par_pos)*Escape_Prob_Tip(eta_f, par_pos) * h_area
    else
      ! The field is NOT favorable for emission
      ! This point does not contribute
      ff(1) = 0.0d0
    end if

    ! We multiply with 2.0*pi (max_xi - 1.0) because Cuba does the 
    ! integration over the hybercube, i.e. from 0 to 1.
    ff(1) = 2.0d0*pi*(max_xi - 1.0d0)*ff(1)
    
    integrand_cuba_fe_tip_test = 0 ! Return value to Cuba, 0 = success
  end function integrand_cuba_fe_tip_test

  ! ----------------------------------------------------------------------------
  ! This subroutine does the Cuba integration for the GTF emission tip.
  integer(c_int) function integrand_cuba_gtf_tip(ndim, xx, ncomp, ff, userdata)
    ! Called from the C Cuba library: the integer arguments and the result
    ! must be integer(c_int) (see mod_cuba_integration)
    use, intrinsic :: iso_c_binding, only: c_int
    integer(c_int), intent(in) :: ndim ! Number of dimensions (Should be 2)
    integer(c_int), intent(in) :: ncomp ! Number of vector-components in the integrand (Always 1 here)
    integer(c_int), intent(in) :: userdata ! Additional data passed to the integral function (In our case the number of the emitter)
    double precision, intent(in), dimension(1:ndim)   :: xx ! Integration points, between 0 and 1
    double precision, intent(out), dimension(1:ncomp) :: ff ! Results of the integrand function

    ! Variables used for calculations
    double precision, dimension(1:3) :: par_pos, field
    !double precision                 :: w_theta ! Emitter area
    double precision                 :: eta_f ! Component normal to the surface
    double precision                 :: xi, phi ! Prolate coordinates
    double precision                 :: h_area ! Area element of the surface

    ! Surface position
    ! Cuba does the intergration over the hybercube.
    ! It gives us coordinates between 0 and 1.
    xi = (max_xi - 1.0d0)*xx(1) + 1.0d0
    phi = 2.0d0*pi*xx(2)

    ! xyz position on the surface
    par_pos = xyz_corr(xi, eta_1, phi)

    ! Calculate the electric field on the surface
    field = Calc_Field_at(par_pos)
    eta_f = Field_normal(par_pos, field)

    ! Add to the average field
    F_avg = F_avg + field

    ! Check if the field is favourable for emission
    if (eta_f < 0.0d0) then
      ! The field is favourable for emission

      ! Area element of the (xi, phi) parametrization, the same expression
      ! Tip_gtf_target_log samples. The scale factors are each singular at the
      ! apex (xi = 1): h_xi diverges and h_phi vanishes, so forming them
      ! separately gives Inf*0 = NaN there. Their product is finite and has to
      ! be formed analytically, because Divonne evaluates the boundary of the
      ! hypercube (border = 0), i.e. xi = 1 exactly.
      h_area = a_foci**2 * sqrt((xi**2 - eta_1**2)*(1 - eta_1**2))

      ! Calculate the current density at this point and convert it to
      ! electrons supplied per time step
      !w_theta = w_theta_xy(par_pos, userdata)
      ff(1) = Get_Kevin_Jgtf_v2(eta_f, T_temp, w_theta) * h_area * time_step_div_q0
    else
      ! The field is NOT favourable for emission
      ! This point does not contribute
      !print *, 'Warning: Field is not favourable for emission'
      ff(1) = 0.0d0
    end if

    ! We mutiply with the range of the coordinates because
    ! Cuba does the  integration over the hybercube, i.e. from 0 to 1.
    ff(1) = ff(1) * 2.0d0*pi*(max_xi - 1.0d0)
    
    integrand_cuba_gtf_tip = 0 ! Return value to Cuba, 0 = success
  end function integrand_cuba_gtf_tip


  ! ----------------------------------------------------------------------------
  ! This subroutine does the Cuba integration for the GTF emission tip.
  subroutine Do_Surface_Integration_GTF_Tip(emit, N_sup)
    ! Input / output variables
    integer, intent(in)           :: emit
    double precision, intent(out) :: N_sup

    integer                          :: nregions, neval, fail
    double precision, dimension(1:1) :: integral, error, prob

    ! Initialize the average field to zero
    F_avg = 0.0d0

    call Cuba_Integrate(integrand_cuba_gtf_tip, 1, emit, integral, error, prob, nregions, neval, fail)

    ! The integrand is in electrons supplied per time step
    N_sup = integral(1)

    ! Finish calculating the average field
    F_avg = F_avg / neval
  end subroutine Do_Surface_Integration_GTF_Tip

  ! ----------------------------------------------------------------------------
  ! This subroutine does the Cuba integration for the field emission tip.
  subroutine Do_Surface_Integration_FE_Tip(emit, N_sup)
    ! Input / output variables
    integer, intent(in)           :: emit
    double precision, intent(out) :: N_sup

    integer                          :: nregions, neval, fail
    double precision, dimension(1:1) :: integral, error, prob

    ! Initialize the average field to zero
    F_avg = 0.0d0

    call Cuba_Integrate(integrand_cuba_fe_tip_test, 1, emit, integral, error, prob, nregions, neval, fail)

    N_sup = integral(1)

    ! Finish calculating the average field
    F_avg = F_avg / neval
  end subroutine Do_Surface_Integration_FE_Tip

  ! Gives a gaussian emission curve
  ! where step is the current time step
  ! returns the number of electrons allowed to be emitted in that time step
  double precision function Gauss_Emission(step)
    integer, intent(in)         :: step ! Current time step
    integer                     :: IFAIL
    double precision, parameter :: sigma = 3000.0d0 ! Width / standard deviation
    double precision, parameter :: mu = 20000.0d0 ! Center
    double precision, parameter :: A = 1.0d0 ! Height
    double precision, parameter :: b = 1.0d0/(2.0d0*pi*sigma**2)

    Gauss_Emission = A * exp( -1.0d0*b*(step - mu)**2 )

    write (ud_gauss, "(i6, tr2, E16.8)", iostat=IFAIL) step, Gauss_Emission
  end function Gauss_Emission

  double precision function Photon_Emission(photon_energy, freq_var)
    double precision, intent(in):: photon_energy, freq_var
    double precision :: rand_photon, photon_rand, Photon_power, Freq, Photo_pois
    
    Photon_power = photon_energy*10000
    Freq = freq_var*10000
    call random_number(rand_photon)
    photon_rand = (Photon_power-Freq) + FLOOR(((Photon_power+Freq)-(Photon_power-Freq))*rand_photon)
    Photo_pois = Rand_Poisson(photon_rand)
    Photon_Emission = Photo_pois/10000

  end function Photon_Emission
end module mod_emission_tip
