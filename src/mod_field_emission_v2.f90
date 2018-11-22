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
  use ziggurat
  implicit none

  PRIVATE
  PUBLIC :: Init_Field_Emission_v2, Clean_Up_Field_Emission_v2, t_y, v_y, Escape_Prob, F_avg, Elec_Supply_V2

  ! ----------------------------------------------------------------------------
  ! Variables
  integer, dimension(:), allocatable :: nrEmitted_emitters
  integer                            :: nrElecEmitAll
  !integer                            :: nrEmitted
  double precision, dimension(1:3)   :: F_avg = 0.0d0
  integer, parameter                 :: N_MH_step = 10 ! Number of steps to do in the MH algorithm

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
  ! This is now declared in the work function submodule
  !double precision, parameter :: w_theta = 2.0d0

  ! Constant used in MC integration (function Elec_Supply_V2)
  double precision :: time_step_div_q0

  interface 
      ! Interface for the work function submodule
      double precision module function w_theta_xy(pos, sec)
        double precision, intent(in), dimension(1:3) :: pos ! Position on the surface
        integer, intent(out), optional               :: sec ! Return the section
      end function w_theta_xy

      module subroutine Read_work_function()
      end subroutine Read_work_function

      module subroutine Work_fun_cleanup()
      end subroutine Work_fun_cleanup

      ! Interface for the MC integration submodule
      module subroutine Do_Surface_Integration_FE(emit, N_sup)
        integer, intent(in)  :: emit ! The emitter to do the integration on
        integer, intent(out) :: N_sup ! Number of electrons
      end subroutine Do_Surface_Integration_FE

      ! Interface for the Metropolis-Hastings submodule
      module subroutine Metropolis_Hastings_rectangle_v2(ndim, emit, df_out, F_out, pos_out)
        integer, intent(in)                           :: ndim, emit
        double precision, intent(out)                 :: df_out, F_out
        double precision, intent(out), dimension(1:3) :: pos_out
      end subroutine Metropolis_Hastings_rectangle_v2

      module subroutine Metropolis_Hastings_rectangle_v2_field(ndim, emit, df_out, F_out, pos_out)
        integer, intent(in)                           :: ndim, emit
        double precision, intent(out)                 :: df_out, F_out
        double precision, intent(out), dimension(1:3) :: pos_out
      end subroutine Metropolis_Hastings_rectangle_v2_field
  end interface
contains
  !-----------------------------------------------------------------------------
  ! Initialize the Field Emission
  subroutine Init_Field_Emission_v2()
    ! Allocate the number of emitters
    allocate(nrEmitted_emitters(1:nrEmit))

    ! Initialize variables
    nrEmitted_emitters = 0 ! Number of electrons emitted from emitter number i in this time step

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_ElecHole_Planar

    ! Function for the electric field in the system
    ptr_field_E => field_E_planar

    ! The function that does the emission
    ptr_Do_Emission => Do_Field_Emission

    call Read_work_function()

    ! Parameters used in the module
    time_step_div_q0 = time_step / q_0

    ! Initialize the Ziggurat algorithm
    call zigset(my_seed(1))

  end subroutine Init_Field_Emission_v2

  subroutine Clean_Up_Field_Emission_v2()
    deallocate(nrEmitted_emitters)

    call Work_fun_cleanup()
  end subroutine Clean_Up_Field_Emission_v2

  !-----------------------------------------------------------------------------
  ! This subroutine gets called from main when the emitters should emit the electrons
  subroutine Do_Field_Emission(step)
    integer, intent(in) :: step
    integer             :: i, IFAIL
    double precision    :: cur_time

    nrElecEmitAll = 0
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
    write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, nrElecEmitAll, nrElec, &
                                                       & (nrEmitted_emitters(i), i = 1, nrEmit)
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
    integer                          :: s, sec, IFAIL
    integer                          :: nrElecEmit
    double precision, dimension(1:3) :: par_pos, par_vel

    call Do_Surface_Integration_FE(emit, N_sup)
    !call Calc_Field_old_method(step, emit)

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
    !$OMP PARALLEL DO PRIVATE(s, par_pos, F, D_f, rnd, par_vel) REDUCTION(+:df_avg) SCHEDULE(AUTO)
    do s = 1, N_sup

      !call Metropolis_Hastings_rectangle_v2(N_MH_step, emit, D_f, F, par_pos)
      call Metropolis_Hastings_rectangle_v2_field(N_MH_step, emit, D_f, F, par_pos)
      !print *, 'D_f = ', D_f
      !print *, 'F = ', F
      !print *, ''
      !pause

      ! Check if the field is favourable for emission or not
      if (F >= 0.0d0) then
        D_f = 0.0d0
        !print *, 'Warning: F > 0.0d0'
      !else
      !  if (D_f > 1.0d0) then
      !    print *, 'Warning D_f > 1.0d0'
      !    print *, 'D_f = ', D_f
      !  end if
      end if
      df_avg = df_avg + D_f

      CALL RANDOM_NUMBER(rnd)
      if (rnd <= D_f) then
        par_pos(3) = 1.0d0*length_scale
        !$OMP CRITICAL(EMIT_PAR)

          ! Add a particle to the system
          par_vel = 0.0d0
          rnd = w_theta_xy(par_pos, sec) ! Get the section
          call Add_Particle(par_pos, par_vel, species_elec, step, emit, sec)

          nrElecEmit = nrElecEmit + 1
          nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
          !call Add_Plane_Graph_emit(par_pos, step)
          !call Add_Plane_Graph_emitt_xy(par_pos)
        !$OMP END CRITICAL(EMIT_PAR)
      end if
    end do
    !$OMP END PARALLEL DO

    df_avg = df_avg / N_sup

    write (ud_field, "(i8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, i8, tr2, E16.8)", iostat=IFAIL) &
                                      step, F_avg(1), F_avg(2), F_avg(3), N_sup, df_avg

    nrElecEmitAll = nrElecEmitAll + nrElecEmit
    !nrEmitted = nrEmitted + nrElecEmit
  end subroutine Do_Field_Emission_Planar_rectangle

  !----------------------------------------------------------------------------------------
  ! The old ways, for comparison
  subroutine Calc_Field_old_method(step, emit)
    integer, intent(in)              :: step, emit
    integer                          :: nr_x, nr_y, i, j, IFAIL
    double precision                 :: len_x, len_y
    double precision, dimension(1:3) :: par_pos, F, F_avg

    nr_x = 1000
    nr_y = nr_x
    len_x = emitters_dim(1, emit) / nr_x
    len_y = emitters_dim(2, emit) / nr_y

    F_avg = 0.0d0

    !$OMP PARALLEL DO PRIVATE(i, j, par_pos, F) REDUCTION(+:F_avg)
    do i = 1, nr_x
      do j = 1, nr_y
        par_pos(1) = (i - 0.5d0)*len_x + emitters_pos(1, emit)
        par_pos(2) = (j - 0.5d0)*len_y + emitters_pos(2, emit)
        par_pos(3) = 0.0d0
        
        F = Calc_Field_at(par_pos)
        F_avg(1:3) = F_avg(1:3) + F(1:3)
      end do
    end do
    !$OMP END PARALLEL DO

    F_avg(1:3) = F_avg(1:3) / (nr_x*nr_y)

    write (ud_debug, "(i8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, i8, tr2, E16.8)", iostat=IFAIL) &
                                      step, F_avg(1), F_avg(2), F_avg(3)
  end subroutine Calc_Field_old_method

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
        print *, 'x = ', pos(1)/length_scale, 'y = ', pos(2)/length_scale, ' z = ,', pos(3)/length_scale
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
  ! A simple function that calculates
  ! A_FN/(t**2(l)*w_theta(x,y)) F**2(x,y)
  ! pos: Position to calculate the function
  ! F: The z-component of the field at par_pos, it should be F < 0.0d0.
  double precision function Elec_Supply_V2(F, pos)
    double precision, dimension(1:3), intent(in) :: pos
    double precision,                 intent(in) :: F

    Elec_Supply_V2 = time_step_div_q0 * a_FN/(t_y(F, pos)**2*w_theta_xy(pos)) * F**2
  end function Elec_Supply_V2

end Module mod_field_emission_v2
