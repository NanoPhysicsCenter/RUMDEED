!-------------------------------------------!
! Module for field emission                 !
! Work function can vary with position      !
! Kristinn Torfason                         !
! 20.06.18                                  !
!-------------------------------------------!

Module mod_field_emission_v2
  use mod_global
  use mod_verlet
  use mod_pair
  implicit none

  PRIVATE
  PUBLIC :: Init_Field_Emission_v2, Clean_Up_Field_Emission_v2, t_y, v_y

  ! ----------------------------------------------------------------------------
  ! Variables
  integer, dimension(:), allocatable :: nrEmitted_emitters
  integer                            :: posInit
  integer                            :: nrEmitted
  double precision                   :: res_s ! residual
  double precision, dimension(1:3)   :: F_avg = 0.0d0

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
  double precision, parameter :: w_theta = 2.0d0

  ! Constant used in MC integration
  double precision :: time_step_div_q0

  ! Use image Charge or not
  logical, parameter          :: image_charge = .true.

  interface
      ! Interface for the work function submodule
      module double precision function w_theta_xy(pos, sec)
        double precision, intent(in), dimension(1:3) :: pos ! Position on the surface
        integer, intent(out), optional               :: sec ! Return the section
      end function w_theta_xy

      ! Interface for the MC integration submodule
      module subroutine Do_Surface_Integration(emit, N_sup)
        integer, intent(in)  :: emit ! The emitter to do the integration on
        integer, intent(out) :: N_sup ! Number of electrons
      end subroutine Do_Surface_Integration
  end interface
contains
  !-----------------------------------------------------------------------------
  ! Initialize the Field Emission
  subroutine Init_Field_Emission_v2()
    ! Allocate the number of emitters
    allocate(nrEmitted_emitters(1:nrEmit))

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_ElecHole_Planar

    ! Function for the electric field in the system
    ptr_field_E => field_E_planar

    ! The function that does the emission
    ptr_Do_Emission => Do_Field_Emission

    ! Parameters used in the module
    time_step_div_q0 = time_step / q_0
  end subroutine Init_Field_Emission_v2

  subroutine Clean_Up_Field_Emission_v2()
    deallocate(nrEmitted_emitters)
  end subroutine Clean_Up_Field_Emission_v2

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
        !  print *, 'Vacuum: WARNING unknown emitter type!!'
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

    ! Integration
    integer                          :: N_sup

    ! Emission variables
    double precision                 :: D_f, Df_avg, F, rnd
    integer                          :: s, IFAIL
    integer                          :: nrElecEmit
    double precision, dimension(1:3) :: par_pos, par_vel

    call Do_Surface_Integration(emit, N_sup)

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
    !$OMP PARALLEL DO PRIVATE(s, par_pos, F, D_f, rnd, par_vel) REDUCTION(+:df_avg)
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
          call Add_Particle(par_pos, par_vel, species_elec, step, emit)

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
  end function Escape_Prob

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
      df_cur = Escape_Prob(field(3), cur_pos)
    else
      df_cur = 0.0d0 ! Zero escape probabilty if field is not favourable
    end if

    !---------------------------------------------------------------------------
    ! We now pick a random distance and direction to jump to from our
    ! current location. We do this ndim times.
    do i = 1, ndim
      ! Find a new position
      new_pos(1:2) = cur_pos(1:2) + box_muller(0.0d0, std)

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
      df_new = Escape_Prob(field(3), new_pos)

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

  ! ----------------------------------------------------------------------------
  !
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

end Module mod_field_emission_v2
