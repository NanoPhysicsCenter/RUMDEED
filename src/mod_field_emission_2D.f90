!-------------------------------------------!
! Module for field emission                 !
! Work function can vary with position      !
! Kristinn Torfason                         !
! 20.06.18                                  !
!-------------------------------------------!

Module mod_field_emission_2D
  use mod_global
  use mod_verlet
  use mod_pair
  implicit none

  PRIVATE
  PUBLIC :: Init_Field_Emission_2D, Clean_Up_Field_Emission_2D

  ! ----------------------------------------------------------------------------
  ! Variables
  integer, dimension(:), allocatable :: nrEmitted_emitters
  integer                            :: posInit
  integer                            :: nrEmitted
  double precision                   :: res_s ! residual

  ! ----------------------------------------------------------------------------
  ! Constants for field emission

  ! V_⊥, L_⊥ and e_⊥
  double precision, parameter :: e_o = 38.3d0*q_0 !Bound state energy [eV -> J]
  double precision, parameter :: v_o = sqrt(2.0d0*e_o/m_0) ! Cross plane velocity [m/s]
  double precision, parameter :: L_o = 0.335d-9 ! 2D material thickness [m] [0.335nm, 0.670nm]

  double precision, parameter :: g_sv = 4.0d0 ! Spin-valley degeneracy []
  double precision, parameter :: m_eff = 0.03d0*m_0 ! Electron effective mass [kg]

  double precision, parameter :: e_f = 0.075d0*q_0 ! Fermi-energy [eV -> J]
  double precision, parameter :: v_f = 10.0d6 ! Fermi-velocity [m/s]

  double precision, parameter :: phi_B0 = 4.5d0*q_0 ! [eV -> J]
  double precision, parameter :: phi_B = phi_B0 - e_f ! [eV -> J]

  double precision, parameter :: lambda = 10.0d-4 ! Strength of the NCLM scattering processes

  ! Use image Charge or not
  !logical, parameter          :: image_charge = .true.

  interface
    double precision function Elec_Supply(A)
      double precision, intent(in) :: A ! Area of the emitter
    end function Elec_Supply

    double precision function Escape_Prob(F)
      double precision, intent(in) :: F ! Electric field
    end function Escape_prob
  end interface

  procedure(Elec_Supply), pointer    :: ptr_Elec_Supply => null()
  procedure(Escape_Prob), pointer    :: ptr_Escape_Prob => null()

contains
  !-----------------------------------------------------------------------------
  ! Initialize the Field Emission
  subroutine Init_Field_Emission_2D()
    ! Allocate the number of emitters
    allocate(nrEmitted_emitters(1:nrEmit))

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_Planar

    ! Function for the electric field in the system
    ptr_field_E => field_E_planar

    ! The function that does the emission
    ptr_Do_Emission => Do_Field_Emission

    ! The function to do image charge effects
    ptr_Image_Charge_effect => Force_Image_charges_v2

    SELECT CASE (EMISSION_MODE)
    case(EMISSION_FIELD_2D_2DEG_C)
      ptr_Elec_Supply => Elec_Supply_2DEG_C
      ptr_Escape_Prob => Escape_Prob_2DEG_C
    case(EMISSION_FIELD_2D_2DEG_NC)
      ptr_Elec_Supply => Elec_Supply_2DEG_NC
      ptr_Escape_Prob => Escape_Prob_2DEG_NC
    case(EMISSION_FIELD_2D_DIRAC_C)
      ptr_Elec_Supply => Elec_Supply_DIRAC_C
      ptr_Escape_Prob => Escape_Prob_DIRAC_C
    case(EMISSION_FIELD_2D_DIRAC_NC)
      ptr_Elec_Supply => Elec_Supply_DIRAC_NC
      ptr_Escape_Prob => Escape_Prob_DIRAC_NC
    END SELECT
  end subroutine Init_Field_Emission_2D

  subroutine Clean_Up_Field_Emission_2D()
    deallocate(nrEmitted_emitters)
  end subroutine Clean_Up_Field_Emission_2D

  !-----------------------------------------------------------------------------
  ! This subroutine gets called from main when the emitters should emit the electrons
  subroutine Do_Field_Emission(step)
    integer, intent(in) :: step
    integer             :: i, IFAIL
    double precision    :: cur_time

    posInit = 0
    nrEmitted_emitters = 0

    ! Loop through all of the emitters
    do i = 1, nrEmit

      ! Check the type of the emitter CIRCLE / RECTANGLE
      if (emitters_delay(i) < step) then
        !select case (emitters_type(i))
        !case (EMIT_CIRCLE)
          !print *, 'Doing Circle'
        !  call Do_Photo_Emission_Circle(step, i)
        !case (EMIT_RECTANGLE)
          !print *, 'Doing Rectangle'
          call Do_Field_Emission_Planar_rectangle(step, i)
        !case default
        !  print *, 'RUMDEED: WARNING unknown emitter type!!'
        !  print *, emitters_type(i)
        !end select
      end if
    end do

    cur_time = time_step * step / time_scale ! Scaled in units of time_scale
    write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, posInit, &
    & nrEmitted, nrElec, (nrEmitted_emitters(i), i = 1, nrEmit)
  end subroutine Do_Field_Emission

  !-----------------------------------------------------------------------------
  ! Do the field emission from a planar rectangular emitter
  ! step: The current time step
  ! emit: Number of the emitter to use
  ! This version allows the work function to vary with position
  subroutine Do_Field_Emission_Planar_rectangle(step, emit)
    integer, intent(in)              :: step, emit

    ! MC integration variables
    double precision                 :: mc_err ! Error in the Monte Carlo integration
    integer                          :: N_mc, Nmc_try
    double precision                 :: A ! Area of the emitter
    double precision, dimension(1:3) :: par_pos, field, F_avg
    double precision                 :: e_sup, e_sup_avg, e_sup_res
    double precision                 :: e_sup2, e_sup_avg2
    double precision                 :: N_sup_db
    integer                          :: N_sup, N_res

    ! Emission variables
    double precision                 :: D_f, Df_avg, F, rnd
    integer                          :: s, IFAIL
    integer                          :: nrElecEmit
    double precision, dimension(1:3) :: par_vel

    ! Calculate the area of the emitter
    A = emitters_dim(1, emit) * emitters_dim(2, emit)

    ! Calculate the electron supply
    N_sup_db = ptr_Elec_Supply(A)
    N_sup = nint(N_sup_db) ! Round the number to integer

    !print *, 'V_2'
    !print *, N_sup_db
    !print *, N_sup
    !print *, mc_err
    !print *, N_mc
    !pause

    !---------------------------------------------------------------------------
    ! Loop over the supply of electrons and place them on the emitter.
    ! Then test the emission probability if we emit them.

    nrElecEmit = 0
    nrEmitted_emitters(emit) = 0

    ! Set the average escape probability to zero.
    df_avg = 0.0d0

    ! Loop over the electrons to be emitted.
    !$OMP PARALLEL DO PRIVATE(s, par_pos, field, F, D_f, rnd, par_vel) REDUCTION(+:df_avg)
    do s = 1, N_sup

      par_pos = Metropolis_Hastings_rectangle_v2(30, emit, D_f, F)
      !print *, 'D_f = ', D_f
      !print *, 'F = ', F
      !print *, ''
      !pause

      ! Check if the field is favourable for emission or not
      if (F >= 0.0d0) then
        D_f = 0.0d0
        !print *, 'Warning: F > 0.0d0'
      else
        if (D_f > 1.0d0) then
          print *, 'Warning D_f > 1.0d0'
          print *, 'D_f = ', D_f
        end if
      end if
      df_avg = df_avg + D_f

      CALL RANDOM_NUMBER(rnd)
      if (rnd <= D_f) then
        par_pos(3) = 1.0d0*length_scale
        !$OMP CRITICAL(EMIT_PAR)

          ! Add a particle to the system
          par_vel = 0.0d0
          call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1)

          nrElecEmit = nrElecEmit + 1
          nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
          !call Add_Plane_Graph_emit(par_pos, step)
          !call Add_Plane_Graph_emitt_xy(par_pos)
        !$OMP END CRITICAL(EMIT_PAR)
      end if
    end do
    !$OMP END PARALLEL DO

    df_avg = df_avg / N_sup

    write (ud_debug, "(i8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, i8, tr2, E16.8)", iostat=IFAIL) &
                                      step, F_avg(1), F_avg(2), F_avg(3), N_sup, df_avg

    posInit = posInit + nrElecEmit
    nrEmitted = nrEmitted + nrElecEmit

    !print *, ''
    !print *, 'posInit'
    !print *, posInit
    !print *, ''
    !pause
  end subroutine Do_Field_Emission_Planar_rectangle

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
  double precision function Elec_Supply_2DEG_C(A)
    double precision, intent(in)                 :: A

    Elec_supply_2DEG_C = A*time_step/q_0 * v_o/L_o * g_sv*q_0*m_eff/(2.0d0*pi*h_bar**2) * e_f

  end function Elec_supply_2DEG_C

  ! This function returns the escape probability of the Electrons.
  ! Escape_prob = exp(-b_FN*w_theta^(3/2)*v_y/F) .
  double precision function Escape_Prob_2DEG_C(F)
    double precision, intent(in)                 :: F
    double precision                             :: d_f, Df

    d_f = q_0*h_bar*(-1.0d0*F)/(2.0d0*sqrt(2.0d0*m_0*phi_B))
    Df = exp(-2.0d0*phi_B/(3.0d0*d_f))

    !Escape_Prob = exp(b_FN * (sqrt(w_theta_xy(pos)))**3 * v_y(F, pos) / (-1.0d0*F))
    Escape_Prob_2DEG_C = Df * exp(-e_f/d_f)

    if (Escape_Prob_2DEG_C > 1.0d0) then
      print *, 'Escape_prob is larger than 1.0'
      print *, 'Escape_prob = ', Escape_Prob_2DEG_C
      print *, ''
      !print *, 'F = ', F
      !print *, 'd_f = ', d_f
      !print *, 'Df = ', Df
      !pause
    end if
  end function Escape_Prob_2DEG_C

  double precision function Elec_Supply_DIRAC_C(A)
    double precision, intent(in)                 :: A

    Elec_supply_DIRAC_C = A*time_step/q_0 * v_o/L_o * g_sv*q_0/(2.0d0*pi*h_bar**2) * e_f**2/(2.0d0*v_f**2)

  end function Elec_supply_DIRAC_C

  ! This function returns the escape probability of the Electrons.
  ! Escape_prob = exp(-b_FN*w_theta^(3/2)*v_y/F) .
  double precision function Escape_Prob_DIRAC_C(F)
    double precision, intent(in)                 :: F
    double precision                             :: d_f, Df

    d_f = q_0*h_bar*(-1.0d0*F)/(2.0d0*sqrt(2.0d0*m_0*phi_B))
    Df = exp(-2.0d0*phi_B/(3.0d0*d_f))

    !Escape_Prob = exp(b_FN * (sqrt(w_theta_xy(pos)))**3 * v_y(F, pos) / (-1.0d0*F))
    Escape_Prob_DIRAC_C = Df * exp(-e_f/d_f)

    if (Escape_Prob_DIRAC_C > 1.0d0) then
      print *, 'Escape_prob is larger than 1.0'
      print *, 'Escape_prob = ', Escape_Prob_DIRAC_C
      print *, ''
      !print *, 'F = ', F
      !print *, 'd_f = ', d_f
      !print *, 'Df = ', Df
      !pause
    end if
  end function Escape_Prob_DIRAC_C


  double precision function Elec_Supply_2DEG_NC(A)
    double precision, intent(in)                 :: A

    Elec_supply_2DEG_NC = lambda * Elec_Supply_2DEG_C(A)

  end function Elec_supply_2DEG_NC

  ! This function returns the escape probability of the Electrons.
  ! Escape_prob = exp(-b_FN*w_theta^(3/2)*v_y/F) .
  double precision function Escape_Prob_2DEG_NC(F)
    double precision, intent(in)                 :: F
    double precision                             :: d_f, Df

    d_f = q_0*h_bar*(-1.0d0*F)/(2.0d0*sqrt(2.0d0*m_0*phi_B))
    Df = exp(-2.0d0*phi_B/(3.0d0*d_f))

    !Escape_Prob = exp(b_FN * (sqrt(w_theta_xy(pos)))**3 * v_y(F, pos) / (-1.0d0*F))
    Escape_Prob_2DEG_NC = Df * (1.0d0 - exp(-e_f/d_f)) * d_f

    if (Escape_Prob_2DEG_NC > 1.0d0) then
      print *, 'Escape_prob is larger than 1.0'
      print *, 'Escape_prob = ', Escape_Prob_2DEG_NC
      print *, ''
      !print *, 'F = ', F
      !print *, 'd_f = ', d_f
      !print *, 'Df = ', Df
      !pause
    end if
  end function Escape_Prob_2DEG_NC

  double precision function Elec_Supply_DIRAC_NC(A)
    double precision, intent(in)                 :: A

    Elec_supply_DIRAC_NC = lambda * Elec_Supply_DIRAC_C(A)

  end function Elec_supply_DIRAC_NC

  ! This function returns the escape probability of the Electrons.
  ! Escape_prob = exp(-b_FN*w_theta^(3/2)*v_y/F) .
  double precision function Escape_Prob_DIRAC_NC(F)
    double precision, intent(in)                 :: F
    double precision                             :: d_f, Df

    d_f = q_0*h_bar*(-1.0d0*F)/(2.0d0*sqrt(2.0d0*m_0*phi_B))
    Df = exp(-2.0d0*phi_B/(3.0d0*d_f))

    !Escape_Prob = exp(b_FN * (sqrt(w_theta_xy(pos)))**3 * v_y(F, pos) / (-1.0d0*F))
    Escape_Prob_DIRAC_NC = Df * (e_f/d_f + exp(-e_f/d_f) - 1.0d0)*d_f**2

    if (Escape_Prob_DIRAC_NC > 1.0d0) then
      print *, 'Escape_prob is larger than 1.0'
      print *, 'Escape_prob = ', Escape_Prob_DIRAC_NC
      print *, ''
      !print *, 'F = ', F
      !print *, 'd_f = ', d_f
      !print *, 'Df = ', Df
      !pause
    end if
  end function Escape_Prob_DIRAC_NC


  !-----------------------------------------------------------------------------
  ! Metropolis-Hastings algorithm
  ! Includes that the work function can vary with position
  function Metropolis_Hastings_rectangle_v2(ndim, emit, df_out, F_out)
    integer, intent(in)              :: ndim, emit
    double precision, intent(out)    :: df_out, F_out
    double precision, dimension(1:3) :: Metropolis_Hastings_rectangle_v2
    integer                          :: count, i
    double precision                 :: std, rnd, alpha
    double precision, dimension(1:3) :: cur_pos, new_pos, field
    double precision                 :: df_cur, df_new

    std = (emitters_dim(1, emit)*0.05d0 + emitters_dim(2, emit)*0.05d0) / (2.0d0)

    ! Get a random position on the surface
    count = 0
    do ! Infinite loop
      CALL RANDOM_NUMBER(cur_pos(1:2))
      cur_pos(1:2) = cur_pos(1:2)*emitters_dim(1:2, emit) + emitters_pos(1:2, emit)
      cur_pos(3) = 0.0d0 ! On the surface

      ! Calculate the electric field at this position
      field = Calc_Field_at(cur_pos)
      if (field(3) < 0.0d0) then
        exit ! The loop is infinite
      else
        count = count + 1
        if (count > 1000) exit ! The loop is infnite, must stop it at some point
      end if
    end do

    F_out = field(3)

    ! Calculate the escape probability at this location
    if (field(3) < 0.0d0) then
      df_cur = ptr_Escape_Prob(field(3))
    else
      df_cur = 0.0d0 ! Zero escape probabilty if field is not favourable
    end if

    !---------------------------------------------------------------------------
    ! We now pick a random distance and direction to jump to from our
    ! current location. We do this ndim times.
    do i = 1, ndim
      ! Find a new position
      new_pos(1:2) = box_muller(cur_pos(1:2), (/std, std/))

      ! Make sure that the new position is within the limits of the emitter area.
      call check_limits_metro_rec(new_pos, emit)

      ! Calculate the field at the new position
      field = Calc_Field_at(new_pos)
      F_out = field(3)

      ! Check if the field is favourable for emission at the new position.
      ! If it is not then cycle, i.e. we reject this location and
      ! pick another one.
      if (field(3) > 0.0d0) cycle ! Do the next loop iteration

      ! Calculate the escape probability at the new position, to compair with
      ! the current position.
      df_new = ptr_Escape_Prob(field(3))

      ! If the escape probability is higher in the new location,
      ! then we jump to that location. If it is not then we jump to that
      ! location with the probabilty df_new / df_cur.
      if (df_new > df_cur) then
        cur_pos = new_pos ! New position becomes the current position
        df_cur = df_new
        F_out = field(3)
      else
        alpha = df_new / df_cur

        CALL RANDOM_NUMBER(rnd)
        if (rnd < alpha) then
          cur_pos = new_pos ! New position becomes the current position
          df_cur = df_new
          F_out = field(3)
        end if
      end if
    end do

    ! Return the current position
    Metropolis_Hastings_rectangle_v2 = cur_pos
    df_out = df_cur
  end function Metropolis_Hastings_rectangle_v2

  subroutine check_limits_metro_rec(par_pos, emit)
    double precision, dimension(1:3), intent(inout) :: par_pos
    integer, intent(in)                             :: emit
    double precision                                :: x_max, x_min, y_max, y_min
    double precision                                :: d_x, d_y


    x_max = emitters_pos(1, emit) + emitters_dim(1, emit)
    x_min = emitters_pos(1, emit)

    y_max = emitters_pos(2, emit) + emitters_dim(2, emit)
    y_min = emitters_pos(2, emit)

    !Check x ----------------------------------------
    if (par_pos(1) > x_max) then
      d_x = par_pos(1) - x_max
      par_pos(1) = x_max - d_x

      if(d_x > emitters_dim(1, emit)) then
        print *, 'Warning: d_x to large >'
        print *, d_x
      end if
    end if

    if (par_pos(1) < x_min) then
      d_x = x_min - par_pos(1)
      par_pos(1) = d_x + x_min

      if(d_x > emitters_dim(1, emit)) then
        print *, 'Warning: d_x to large <'
        print *, d_x
      end if
    end if

    !Check y ----------------------------------------
    if (par_pos(2) > y_max) then
      d_y = par_pos(2) - y_max
      par_pos(2) = y_max - d_y

      if(d_y > emitters_dim(2, emit)) then
        print *, 'Warning: d_y to large >'
        print *, d_y
      end if
    end if

    if (par_pos(2) < y_min) then
      d_y = y_min - par_pos(2)
      par_pos(2) = d_y + y_min

      if(d_y > emitters_dim(2, emit)) then
        print *, 'Warning: d_x to large <'
        print *, d_y
      end if
    end if
  end subroutine check_limits_metro_rec

end Module mod_field_emission_2D
