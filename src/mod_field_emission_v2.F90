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
  use mod_work_function
  use mod_polarso
  use mod_cuba_integration
  implicit none

  PRIVATE
  PUBLIC :: Init_Field_Emission_v2, Clean_Up_Field_Emission_v2, t_y, v_y, Escape_Prob, F_avg, Elec_Supply_V2, &
            Test_Field_Emission_Module, E_zunit_planar

  ! ----------------------------------------------------------------------------
  ! Variables
  integer, dimension(:), allocatable :: nrEmitted_emitters
  integer                            :: nrElecEmitAll
  !integer                            :: nrEmitted
  double precision, dimension(1:3)   :: F_avg = 0.0d0
  integer, parameter                 :: N_MH_step = 10*3 ! Number of steps to do in the MH algorithm
  double precision, dimension(:), allocatable :: residual ! Rounding residual, one per emitter

  ! Cuba
  !double precision, allocatable :: xgiven(:,:) ! xgiven(ldxgiven,ngiven) <in>, a list of points where the integrand

  ! ----------------------------------------------------------------------------
  ! Constants for field emission
  ! First Fowler-Nordheim constant in units [ A eV V^{-2} ]
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

  ! MH Acceptance rate
  double precision :: a_rate = 1.0d0
  double precision :: MH_std = 0.0125d0

  integer          :: jump_a = 0, jump_r = 0 ! Number of jumps accepted and rejected

  ! Tuning constants of the Metropolis-Hastings sampling, shared by the serial
  ! (Metropolis_Hastings_rectangle_J / Metropolis_Hastings_ring_J) and the
  ! batched (Metropolis_Hastings_rectangle_J_batch) version, so all variants
  ! use the same schedule.
  ! The proposal step MH_std is adapted towards MH_target_rate with one
  ! multiplicative update per acceptance-rate window (one chain in the serial
  ! samplers, one jump iteration in the batched sampler), see MH_std_update.
  ! Within a chain the proposal kernel is fixed.
  integer,          parameter :: MH_ndim        = 25*8      ! Number of jumps per chain
  double precision, parameter :: MH_warmup_frac = 0.25d0    ! Fraction of the jumps done with the fixed initial std
  double precision, parameter :: MH_init_std    = 0.10d0    ! Initial std as a fraction of the emitter size
  double precision, parameter :: MH_target_rate = 0.2525d0  ! Acceptance rate the MH_std adaptation aims for
  double precision, parameter :: MH_std_gain    = 0.025d0   ! Gain of the multiplicative MH_std update per window
  double precision, parameter :: MH_std_max     = 0.1250d0  ! Upper limit on MH_std
  double precision, parameter :: MH_std_min     = 0.00005d0 ! Lower limit on MH_std

  ! interface 
  !     ! ! Interface for the work function submodule
  !     ! double precision module function w_theta_xy(pos, sec)
  !     !   double precision, intent(in), dimension(1:3) :: pos ! Position on the surface
  !     !   integer, intent(out), optional               :: sec ! Return the section
  !     ! end function w_theta_xy

  !     ! module subroutine Read_work_function()
  !     ! end subroutine Read_work_function

  !     ! module subroutine Work_fun_cleanup()
  !     ! end subroutine Work_fun_cleanup

  !     ! Interface for the MC integration submodule
  !     !module subroutine Do_Surface_Integration_FE(emit, N_sup)
  !     !  integer, intent(in)           :: emit ! The emitter to do the integration on
  !     !  double precision, intent(out) :: N_sup ! Number of electrons
  !     !end subroutine Do_Surface_Integration_FE

  !     ! Interface for the Metropolis-Hastings submodule
  !     module subroutine Metropolis_Hastings_rectangle_v2(ndim, emit, df_out, F_out, pos_out)
  !       integer, intent(in)                           :: ndim, emit
  !       double precision, intent(out)                 :: df_out, F_out
  !       double precision, intent(out), dimension(1:3) :: pos_out
  !     end subroutine Metropolis_Hastings_rectangle_v2

  !     module subroutine Metropolis_Hastings_rectangle_v2_field(ndim, emit, df_out, F_out, pos_out)
  !       integer, intent(in)                           :: ndim, emit
  !       double precision, intent(out)                 :: df_out, F_out
  !       double precision, intent(out), dimension(1:3) :: pos_out
  !     end subroutine Metropolis_Hastings_rectangle_v2_field
  ! end interface
contains
  !-----------------------------------------------------------------------------
  ! Initialize the Field Emission
  subroutine Init_Field_Emission_v2()
    ! Allocate the number of emitters
    allocate(nrEmitted_emitters(1:nrEmit))
    allocate(residual(1:nrEmit))

    ! Initialize variables
    nrEmitted_emitters = 0 ! Number of electrons emitted from emitter number i in this time step
    residual = 0.0d0 ! Rounding residual carried over between time steps, per emitter

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_Planar

    ! Function for the electric field in the system
    ptr_field_E => field_E_planar

    ! The function that does the emission
    ptr_Do_Emission => Do_Field_Emission

    ! The function to do image charge effects
    ptr_Image_Charge_effect => Force_Image_charges_v2

    ! Ramo E_z
    ptr_E_zunit => E_zunit_planar

    call Read_work_function()

    ! Parameters used in the module
    time_step_div_q0 = time_step / q_0

    ! Initialize the Ziggurat algorithm
    !call zigset(my_seed(1))

    !allocate(xgiven(1:2, 1:MAX_PARTICLES))

  end subroutine Init_Field_Emission_v2

  subroutine Clean_Up_Field_Emission_v2()
    deallocate(nrEmitted_emitters)
    deallocate(residual)

    call Work_fun_cleanup()

    !deallocate(xgiven)
  end subroutine Clean_Up_Field_Emission_v2

  function E_zunit_planar(pos)
    double precision, dimension(1:3), intent(in) :: pos ! Position to calculate the E_z unit vector at
    double precision, dimension(1:3)               :: E_zunit_planar

    ! The z-component of the electric field is always -1.0d0/d in planar geometry
    E_zunit_planar(1:2) = 0.0d0
    E_zunit_planar(3) = -1.0d0/d

  end function E_zunit_planar

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
        select case (emitters_type(i))
          case (EMIT_RING)
            ! print *, 'Doing Ring'
            call Do_Field_Emission_Planar_ring(step, i)
          case (EMIT_CIRCLE, EMIT_RECTANGLE)
            ! print *, 'Doing Rectangle'
            ! The circle emitter is sampled over its bounding box, the same way
            ! it was before the ring emitter added this dispatch. Without this
            ! case a circle emitter would fall through to the default below and
            ! silently stop emitting.
            call Do_Field_Emission_Planar_rectangle(step, i)
            ! call Do_Field_Emission_Planar_simple(step, i)
          case default
            print *, 'RUMDEED: WARNING unknown emitter type!!'
            print *, emitters_type(i)
        end select
      end if
    end do

    cur_time = time_step * step / time_scale ! Scaled in units of time_scale
    write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, nrElecEmitAll, nrElec, &
                                                       & (nrEmitted_emitters(i), i = 1, nrEmit)
  end subroutine Do_Field_Emission

  subroutine Do_Field_Emission_Planar_simple(step, emit)
    integer, intent(in)              :: step, emit

    ! Integration
    double precision                 :: N_sup
    integer                          :: N_round

    ! Emission variables
    integer                          :: i, sec, nrElecEmit, IFAIL
    double precision                 :: rnd, D_f, F, Df_avg
    double precision, dimension(1:3) :: par_pos, par_vel


    ! Do integration
    call Do_Surface_Integration_Simple(emit, N_sup)

    N_round = nint(N_sup + residual(emit)) ! Round to whole number
    residual(emit) = N_sup - N_round

    nrElecEmit = 0

    ! Loop over all electrons and place them
    do i = 1, N_round
      call Metropolis_Hastings_rectangle_v2(N_MH_step, emit, D_f, F, par_pos)
    

      par_pos(3) = 1.0d0*length_scale
      par_vel = 0.0d0
      rnd = w_theta_xy(par_pos, emit, sec) ! Get the section

      ! Add a particle to the system
      call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1, sec)
      ! Recompute field after adding particle
      if (use_polarso .eqv. .true.) then
        call PL_Update_Field(nrPart)
      end if

      nrElecEmit = nrElecEmit + 1
      nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
    end do

    Df_avg = 0.0d0

    write (ud_field, "(i8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8)", iostat=IFAIL) &
                                      step, F_avg(1), F_avg(2), F_avg(3), N_sup, Df_avg, a_rate, MH_std

    nrElecEmitAll = nrElecEmitAll + nrElecEmit
  end subroutine

  !-----------------------------------------------------------------------------
  ! Do the field emission from a planar rectangular emitter
  ! step: The current time step
  ! emit: Number of the emitter to use
  ! This version allows the work function to vary with position
  subroutine Do_Field_Emission_Planar_rectangle(step, emit)
    integer, intent(in)              :: step, emit

    ! Integration
    double precision                 :: N_sup
    integer                          :: N_round, i

    ! Emission variables
    double precision                 :: D_f, Df_avg, F, rnd
    integer                          :: s, sec, IFAIL
    integer                          :: nrElecEmit
    double precision, dimension(1:3) :: par_pos, par_vel

    ! Batched Metropolis-Hastings results (used when mh_batch is .true.)
    double precision, allocatable    :: mh_df(:), mh_F(:), mh_pos(:, :)

    ! The particle configuration does not change during the surface
    ! integration and the Metropolis-Hastings sampling, so upload it to the
    ! device once and let all the batched field evaluations reuse it
    ! (a no-op in builds without OpenACC).
    call Particles_To_Device()

    call Do_Surface_Integration_FE(emit, N_sup)
    N_round = nint(N_sup + residual(emit))
    residual(emit) = N_sup - N_round

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

    ! Run the Metropolis-Hastings chains for all electrons in the supply.
    ! With mh_batch the chains advance in lockstep and the surface field of
    ! all chains is evaluated in batches (one GPU kernel per jump iteration
    ! in OpenACC builds). Otherwise each chain runs to completion in turn,
    ! doing one field evaluation per jump.
    if ((mh_batch .eqv. .true.) .and. (N_round > 0)) then
      allocate(mh_df(1:N_round), mh_F(1:N_round), mh_pos(1:3, 1:N_round))
      call Metropolis_Hastings_rectangle_J_batch(N_round, emit, mh_df, mh_F, mh_pos)
    end if

    ! Close the device residency window: from here on electrons are added to
    ! the system, so the particle configuration changes.
    call Release_Device_Particles()

    ! Loop over the electrons to be emitted.
    !!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(s, par_pos, F, D_f, rnd, par_vel) REDUCTION(+:df_avg) SCHEDULE(GUIDED, CHUNK_SIZE)
    do s = 1, N_round

      if (mh_batch .eqv. .true.) then
        D_f = mh_df(s)
        F = mh_F(s)
        par_pos = mh_pos(:, s)
      else
        !call Metropolis_Hastings_rectangle_v2(N_MH_step, emit, D_f, F, par_pos)
        i = N_MH_step
        call Metropolis_Hastings_rectangle_J(i, emit, D_f, F, par_pos)
        !call Metropolis_Hastings_rectangle_v2_field(N_MH_step, emit, D_f, F, par_pos)
        !print *, 'D_f = ', D_f
        !print *, 'F = ', F
        !print *, ''
        !pause
      end if

      ! Check if the field is favorable for emission or not
      if (F >= 0.0d0) then
        D_f = -huge(1.0d0)
        print *, 'Warning: F > 0.0d0'
        print *, 'done'
      !else
      !  if (D_f > 1.0d0) then
      !    print *, 'Warning D_f > 1.0d0'
      !    print *, 'D_f = ', D_f
      !  end if
      end if
      df_avg = df_avg + exp(D_f)

      CALL RANDOM_NUMBER(rnd)
      if (log(rnd) <= D_f) then
        !par_vel(1:2) = box_muller((/1.0d0, 1.0d0/), (/0.25d0, 0.25d0/))
        !if (par_vel(1) < 0.0d0) then
        !  if (par_vel(2) < 0.0d0) then
        !    par_pos(3) = 0.0d-6 * length_scale
        !  else
        !    par_pos(3) = par_vel(2) * length_scale
        !  end if
        !else
        !  par_pos(3) = par_vel(1) * length_scale
        !end if
        par_pos(3) = 1.0d0*length_scale
        par_vel = 0.0d0
        rnd = w_theta_xy(par_pos, emit, sec) ! Get the section
        !!$OMP CRITICAL(EMIT_PAR)

          ! Add a particle to the system
          call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1, sec)
          ! Recompute field after adding particle
          if (use_polarso .eqv. .true.) then
            call PL_Update_Field(nrPart)
          end if

          nrElecEmit = nrElecEmit + 1
          nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
        !!$OMP END CRITICAL(EMIT_PAR)
      end if
    end do
    !!$OMP END PARALLEL DO

    if (allocated(mh_df)) then
      deallocate(mh_df, mh_F, mh_pos)
    end if

    df_avg = df_avg / N_sup

    write (ud_field, "(i8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8)", iostat=IFAIL) &
                                      step, F_avg(1), F_avg(2), F_avg(3), N_sup, df_avg, a_rate, MH_std

    ! write (ud_field, "(i8, tr2, E16.8)", iostat=IFAIL) step, df_avg
    nrElecEmitAll = nrElecEmitAll + nrElecEmit
    !nrEmitted = nrEmitted + nrElecEmit
  end subroutine Do_Field_Emission_Planar_rectangle

!-----------------------------------------------------------------------------

  subroutine Do_Field_Emission_Planar_ring(step, emit)
    integer, intent(in)              :: step, emit

    ! Integration
    double precision                 :: N_sup
    integer                          :: N_round, i

    ! Emission variables
    double precision                 :: D_f, Df_avg, F, rnd
    integer                          :: s, sec, IFAIL
    integer                          :: nrElecEmit
    double precision, dimension(1:3) :: par_pos, par_vel

    ! print *, 'Do_Field_Emission_Planar_ring'

    call Do_Surface_Integration_FE(emit, N_sup)
    N_round = nint(N_sup + residual(emit))
    residual(emit) = N_sup - N_round

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
    !!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(s, par_pos, F, D_f, rnd, par_vel) REDUCTION(+:df_avg) SCHEDULE(GUIDED, CHUNK_SIZE)
    do s = 1, N_round

      !call Metropolis_Hastings_rectangle_v2(N_MH_step, emit, D_f, F, par_pos)
      i = N_MH_step
      call Metropolis_Hastings_ring_J(i, emit, D_f, F, par_pos)
      !call Metropolis_Hastings_rectangle_v2_field(N_MH_step, emit, D_f, F, par_pos)
      !print *, 'D_f = ', D_f
      !print *, 'F = ', F
      !print *, ''
      !pause

      ! Check if the field is favourable for emission or not
      if (F >= 0.0d0) then
        D_f = -huge(1.0d0)
        print *, 'Warning: F > 0.0d0'
        print *, 'F = ', F
        print *, 'par_pos = ', par_pos
        print *, 'done'
        exit
      !else
      !  if (D_f > 1.0d0) then
      !    print *, 'Warning D_f > 1.0d0'
      !    print *, 'D_f = ', D_f
      !  end if
      end if
      df_avg = df_avg + exp(D_f)

      CALL RANDOM_NUMBER(rnd)
      if (log(rnd) <= D_f) then
        !par_vel(1:2) = box_muller((/1.0d0, 1.0d0/), (/0.25d0, 0.25d0/))
        !if (par_vel(1) < 0.0d0) then
        !  if (par_vel(2) < 0.0d0) then
        !    par_pos(3) = 0.0d-6 * length_scale
        !  else
        !    par_pos(3) = par_vel(2) * length_scale
        !  end if
        !else
        !  par_pos(3) = par_vel(1) * length_scale
        !end if
        par_pos(3) = 1.0d0*length_scale
        par_vel = 0.0d0
        rnd = w_theta_xy(par_pos, emit, sec) ! Get the section
        !!$OMP CRITICAL(EMIT_PAR)

          ! Add a particle to the system
          call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1, sec)

          ! Recompute field after adding particle
          if (use_polarso .eqv. .true.) then
            call PL_Update_Field(nrPart)
          end if

          nrElecEmit = nrElecEmit + 1
          nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
        !!$OMP END CRITICAL(EMIT_PAR)
      end if
    end do
    !!$OMP END PARALLEL DO

    df_avg = df_avg / N_sup

    write (ud_field, "(i8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8)", iostat=IFAIL) &
                                      step, F_avg(1), F_avg(2), F_avg(3), N_sup, df_avg, a_rate, MH_std

    nrElecEmitAll = nrElecEmitAll + nrElecEmit
    !nrEmitted = nrEmitted + nrElecEmit
  end subroutine Do_Field_Emission_Planar_ring

!----------------------------------------------------------------------------------------
! The functions v_y and t_y are because of the image charge effect in the FN equation.
! The approximation for v_y and t_y are taken from
! Forbes, R. G., & Deane, J. H. (2007, November).
! "Reformulation of the standard theory of Fowler–Nordheim tunnelling and cold field electron emission."
! In Proceedings of the Royal Society of London A: Mathematical,
! Physical and Engineering Sciences (Vol. 463, No. 2087, pp. 2907-2927). The Royal Society.
!
  double precision function v_y(F, pos, emit)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos
    integer, intent(in)                          :: emit
    double precision                             :: l

    if (image_charge .eqv. .true.) then
      l = l_const * (-1.0d0*F) / w_theta_xy(pos, emit)**2 ! l = y^2, y = 3.79E-4 * sqrt(F_old) / w_theta
      if (l > 1.0d0) then
        l = 1.0d0
      end if
      v_y = 1.0d0 - l + 1.0d0/6.0d0 * l * log(l)
    else
      v_y = 1.0d0
    end if
  end function v_y

  double precision function t_y(F, pos, emit)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos
    integer, intent(in)                          :: emit
    double precision                             :: l

    if (image_charge .eqv. .true.) then
      l = l_const * (-1.0d0*F) / w_theta_xy(pos, emit)**2 ! l = y^2, y = 3.79E-4 * sqrt(F_old) / w_theta
      if (l > 1.0d0) then
        print *, 'Error: l > 1.0'
        print *, 'l = ', l, ', F = ', F
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
  double precision function Escape_Prob(F, pos, emit)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos
    integer, intent(in)                          :: emit

    Escape_Prob = exp(b_FN * (sqrt(w_theta_xy(pos, emit)))**3 * v_y(F, pos, emit) / (-1.0d0*F))

  end function Escape_Prob

  ! Returns the log of the escape probability.
  double precision function Escape_Prob_log(F, pos, emit)
  double precision, intent(in)                 :: F
  double precision, dimension(1:3), intent(in) :: pos
  integer, intent(in)                          :: emit

  Escape_Prob_log = b_FN * (sqrt(w_theta_xy(pos, emit)))**3 * v_y(F, pos, emit) / (-1.0d0*F)

end function Escape_Prob_log

  !-----------------------------------------------------------------------------
  ! Log of the position dependent part of the electron supply,
  ! log( F^2 / (t_y^2 * w_theta) ).
  ! This is the target distribution of the Metropolis-Hastings chains below:
  ! the number of electrons to try comes from the surface integral of the
  ! supply, the chains place them proportional to the supply, and the caller
  ! then emits each one with the escape probability, so the positions of the
  ! emitted electrons follow the current density J = supply * D.
  ! (Sampling proportional to D instead would count D twice, once in the
  ! placement and once in the emission test.)
  ! The constant prefactor a_FN * time_step / q_0 of the supply cancels in
  ! the Metropolis ratio and is left out. F must be < 0.
  double precision function Elec_Supply_log(F, pos, emit)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos
    integer, intent(in)                          :: emit

    Elec_Supply_log = 2.0d0*log(-1.0d0*F) - 2.0d0*log(t_y(F, pos, emit)) - log(w_theta_xy(pos, emit))
  end function Elec_Supply_log

  !-----------------------------------------------------------------------------
  ! One multiplicative update of the shared proposal step MH_std, driven by
  ! the acceptance rate measured over the last window (one chain in the
  ! serial samplers, one jump iteration in the batched sampler).
  ! Low acceptance shrinks the step, high acceptance grows it. The
  ! exponential form keeps MH_std positive and the clamps bound the excursion.
  subroutine MH_std_update(rate)
    double precision, intent(in) :: rate

    MH_std = MH_std * exp(MH_std_gain*(rate - MH_target_rate))
    if (MH_std > MH_std_max) then
      MH_std = MH_std_max
    else if (MH_std < MH_std_min) then
      MH_std = MH_std_min
    end if
  end subroutine MH_std_update

  !-----------------------------------------------------------------------------
  ! A simple function that calculates
  ! A_FN/(t**2(l)*w_theta(x,y)) F**2(x,y)
  ! pos: Position to calculate the function
  ! F: The z-component of the field at par_pos, it should be F < 0.0d0.
  double precision function Elec_Supply_V2(F, pos, emit)
    double precision, dimension(1:3), intent(in) :: pos
    double precision,                 intent(in) :: F
    integer, intent(in)                          :: emit

    Elec_Supply_V2 = time_step_div_q0 * a_FN/(t_y(F, pos, emit)**2*w_theta_xy(pos, emit)) * F**2
  end function Elec_Supply_V2


  ! ----------------------------------------------------------------------------
  ! This function is called to do the surface integration.
  ! The method and tolerances are selected with the cuba_* input variables
  ! (see mod_global and mod_cuba_integration). All methods use the vectorized
  ! integrand: Cuba hands it blocks of up to 256 points, which the batched
  ! field evaluation turns into one GPU kernel per block in OpenACC builds.
  !
  subroutine Do_Surface_Integration_FE(emit, N_sup)
    integer, intent(in)           :: emit ! The emitter to do the integration on
    double precision, intent(out) :: N_sup ! Number of electrons

    integer, parameter               :: nvec = 256 ! Number of points given to the integrand per call
    integer                          :: nregions, neval, fail, IFAIL
    double precision, dimension(1:1) :: integral, error, prob

    ! Initialize the average field to zero
    F_avg = 0.0d0

    call Cuba_Integrate(integrand_cuba_fe_v, nvec, emit, integral, error, prob, nregions, neval, fail)

    N_sup = integral(1)

    ! Finish calculating the average field
    F_avg = F_avg / neval

    ! Write the output variables of the integration to a file
    write(ud_integrand, '(i3, tr2, i8, tr2, i8, tr2, i4, tr2, ES12.4, tr2, ES12.4, tr2, ES12.4)', iostat=IFAIL) &
                         & emit, nregions, neval, fail, integral(1), error(1), prob(1)
  end subroutine Do_Surface_Integration_FE

  ! ----------------------------------------------------------------------------
  ! Vectorized integrand for the Cuba library.
  ! When Cuba is called with nvec > 1 it passes blocks of up to nvec
  ! integration points to the integrand in one call, with the actual number of
  ! points in the block as an extra argument (Cuba also appends further
  ! arguments, core and phase, which we do not need and therefore do not
  ! declare).
  ! The integration points and results do not depend on the block size, but
  ! the electric field for the whole block is computed in one batch, which
  ! runs as a single kernel on the GPU in OpenACC builds.
  integer function integrand_cuba_fe_v(ndim, xx, ncomp, ff, userdata, nvec)
    ! Input / output variables
    integer, intent(in) :: ndim ! Number of dimensions (Should be 2)
    integer, intent(in) :: ncomp ! Number of vector-components in the integrand (Always 1 here)
    integer, intent(in) :: userdata ! Additional data passed to the integral function (In our case the number of the emitter)
    integer, intent(in) :: nvec ! Number of integration points in this block
    double precision, intent(in), dimension(1:ndim, 1:nvec)   :: xx ! Integration points, between 0 and 1
    double precision, intent(out), dimension(1:ncomp, 1:nvec) :: ff ! Results of the integrand function

    ! Variables used for calculations
    double precision, dimension(1:3, 1:nvec) :: par_pos, field
    double precision, dimension(1:nvec)      :: A ! Area factor, per point (the ring area depends on the radius)
    double precision                         :: angle, radius
    integer                                  :: k

    ! Surface positions
    ! Cuba does the intergration over the hybercube.
    ! It gives us coordinates between 0 and 1.
    do k = 1, nvec
      select case (emitters_Type(userdata))
        case (EMIT_CIRCLE, EMIT_RECTANGLE)
          ! The circle emitter integrates over its bounding box, like the
          ! rectangle. Leaving it out would hit the default below and stop the run.
          par_pos(1:2, k) = emitters_pos(1:2, userdata) + xx(1:2, k)*emitters_dim(1:2, userdata) ! x and y position on the surface
          par_pos(3, k) = 0.0d0 ! Height, i.e. on the surface
          A(k) = emitters_dim(1, userdata)*emitters_dim(2, userdata)
        case (EMIT_RING)
          angle = xx(1, k)*2.0d0*pi
          radius = emitters_ring(2, userdata) + xx(2, k)*emitters_ring(3, userdata)
          par_pos(1, k) = radius*cos(angle)
          par_pos(2, k) = radius*sin(angle)
          par_pos(3, k) = 0.0d0 ! Height, i.e. on the surface
          A(k) = 2.0d0*pi*emitters_ring(3, userdata)*radius
        case default
          print *, 'RUMDEED: ERROR UNKNOWN EMITTER TYPE'
          stop
      end select
    end do

    ! Calculate the electric field on the surface for the whole block
    if (use_polarso .eqv. .true.) then
      ! The POLARSO solver evaluates the field per point from its grid
      do k = 1, nvec
        field(:, k) = PL_Calculate_Field_At(par_pos(:, k))
      end do
    else
      call Calc_Field_at_Batch(nvec, par_pos, field)
    end if

    do k = 1, nvec
      ! Add to the average field
      F_avg = F_avg + field(:, k)

      ! Check if the field is favorable for emission
      if (field(3, k) < 0.0d0) then
        ! The field is favorable for emission
        ! Calculate the electron supply at this point
        ff(1, k) = Elec_Supply_V2(field(3, k), par_pos(:, k), userdata)
      else
        ! The field is NOT favorable for emission
        ! This point does not contribute
        ff(1, k) = 0.0d0
      end if

      ! We multiply with the area of the emitter because Cuba does the
      ! integration over the hybercube, i.e. from 0 to 1.
      ff(1, k) = A(k)*ff(1, k)
    end do

    integrand_cuba_fe_v = 0 ! Return value to Cuba, 0 = success
  end function integrand_cuba_fe_v



  ! ----------------------------------------------------------------------------
  ! The integration function for the Cuba library
  !
  integer function integrand_cuba_simple(ndim, xx, ncomp, ff, userdata)
    ! Input / output variables
    integer, intent(in) :: ndim ! Number of dimensions (Should be 2)
    integer, intent(in) :: ncomp ! Number of vector-components in the integrand (Always 1 here)
    integer, intent(in) :: userdata ! Additional data passed to the integral function (In our case the number of the emitter)
    double precision, intent(in), dimension(1:ndim)   :: xx ! Integration points, between 0 and 1
    double precision, intent(out), dimension(1:ncomp) :: ff ! Results of the integrand function

    ! Variables used for calculations
    double precision, dimension(1:3) :: par_pos, field
    double precision                 :: A ! Emitter area

    ! Emitter area
    A = emitters_dim(1, userdata)*emitters_dim(2, userdata)

    ! Surface position
    ! Cuba does the intergration over the hybercube.
    ! It gives us coordinates between 0 and 1.
    par_pos(1:2) = emitters_pos(1:2, userdata) + xx(1:2)*emitters_dim(1:2, userdata) ! x and y position on the surface
    par_pos(3) = 0.0d0 ! Height, i.e. on the surface

    ! Calculate the electric field on the surface
    if (use_polarso .eqv. .true.) then
      field = PL_Calculate_Field_At(par_pos)
    else
      field = Calc_Field_at(par_pos)
    end if

    ! Add to the average field
    F_avg = F_avg + field

    ! Check if the field is favorable for emission
    if (field(3) < 0.0d0) then
      ! The field is favorable for emission
      ! Calculate the current density at this point
      ff(1) = Elec_Supply_V2(field(3), par_pos, userdata) * Escape_Prob(field(3), par_pos, userdata)
    else
      ! The field is NOT favorable for emission
      ! This point does not contribute
      ff(1) = 0.0d0
    end if

    ! We mutiply with the area of the emitter because Cuba does the 
    ! integration over the hybercube, i.e. from 0 to 1.
    ff(1) = A*ff(1)
    
    integrand_cuba_simple = 0 ! Return value to Cuba, 0 = success
  end function integrand_cuba_simple

  subroutine Do_Surface_Integration_Simple(emit, N_sup)
    ! Input / output variables
    integer, intent(in)           :: emit
    double precision, intent(out) :: N_sup

    integer                          :: nregions, neval, fail
    double precision, dimension(1:1) :: integral, error, prob

    ! Initialize the average field to zero
    F_avg = 0.0d0

    call Cuba_Integrate(integrand_cuba_simple, 1, emit, integral, error, prob, nregions, neval, fail)

    N_sup = integral(1)

    ! Finish calculating the average field
    F_avg = F_avg / neval
  end subroutine Do_Surface_Integration_Simple

  ! ----------------------------------------------------------------------------
  ! Plain 2D Monte Carlo integration
  !
  subroutine Do_2D_MC_plain_FE(emit, N_sup)
    integer, intent(in)  :: emit ! The emitter to do the integration on
    integer, intent(out) :: N_sup ! Number of electrons

    ! MC integration variables
    double precision                 :: mc_err ! Error in the Monte Carlo integration
    integer                          :: N_mc, Nmc_try
    double precision                 :: A ! Area of the emitter
    double precision, dimension(1:3) :: par_pos, field
    double precision                 :: e_sup, e_sup_avg, e_sup_res
    double precision                 :: e_sup2, e_sup_avg2
    double precision                 :: N_sup_db

    ! Calculate the area of the emitter
    A = emitters_dim(1, emit) * emitters_dim(2, emit)

    !---------------------------------------------------------------------------
    ! Use 2D Monte Carlo integration to calculate the electron supply.
    ! We continue until the error in the integration is less than half an electron.
    !
    ! Monte Carlo Integration
    ! N_sup = \int \Delta t / e * a_FN / (t_y**2*\phi) * F^2 dA
    ! f = \Delta t / e * a_FN / (t_y**2\phi) * F^2
    ! N_sup \approx = A <f> \pm A*sqrt((<f^2> - <f>^2)/N_mc)
    ! <f> = 1/N_mc * sum f(x), <f^2> = 1/N_mc * sum f**2(x)

    mc_err = 1.0d0 ! Set the error in the MC integration to some thing higher than 0.5
    N_mc = 0 ! Number of points in the MC integration
    e_sup = 0.0d0
    e_sup2 = 0.0d0
    F_avg = 0.0d0 ! The average field on the surface
    Nmc_try = 0

    do
      ! Get a random position on the emitter
      CALL RANDOM_NUMBER(par_pos(1:2))
      par_pos(1:2) = emitters_pos(1:2, emit) + par_pos(1:2)*emitters_dim(1:2, emit)
      par_pos(3) = 0.0d0 ! z = 0, emitter surface

      ! Calculate the field on the emitter surface
      if (use_polarso .eqv. .true.) then
        field = PL_Calculate_Field_At(par_pos)
      else
        field = Calc_Field_at(par_pos)
      end if

      ! Check if the field is favorable for emission
      if (field(3) < 0.0d0) then
        Nmc_try = 0
        F_avg(1:3) = F_avg(1:3) + field(1:3)
        N_mc = N_mc + 1

        !Calculate <f> and <f^2>
        e_sup_res = Elec_Supply_V2(field(3), par_pos, emit)
        e_sup = e_sup + e_sup_res
        e_sup2 = e_sup2 + e_sup_res**2

        e_sup_avg = e_sup / N_mc
        e_sup_avg2 = e_sup2 / N_mc

        ! Always do at least 10000 points
        if (N_mc > 1000) then
          ! Calculate the error, A*\sqrt( (<f^2> - <f>^2) / N_mc )
          mc_err = A*sqrt( (e_sup_avg2 - e_sup_avg**2) / N_mc )
          if (mc_err < 1.0d0) exit ! Stop if less than one electron in error
        end if

        ! Stop the integration if it is taking to long.
        if (N_mc > 10000000) then
          print *, 'RUMDEED: Warning MC integration taking to long, stoping it'
          print *, 'mc_err = ', mc_err
          !print *, 'step = ', step
          exit
        end if

      else ! field(3) < 0.0d0
        Nmc_try = Nmc_try + 1

        ! Stop if we are taking to long to find a favorable point.
        ! This should be rare in field emission.
        if (Nmc_try > 1000) then
          print *, 'RUMDEED: Warning to many field attempts at finding a favorable location in MC integration'
          print *, 'mc_err = ', mc_err
          !print *, 'step = ', step
          exit
        end if

      end if ! field(3) < 0.0d0
    end do

    ! Finish calculating the average field on the surface
    F_avg = F_avg / N_mc

    ! Calculate the electron supply
    N_sup_db = A*e_sup_avg
    N_sup = nint(N_sup_db) ! Round the number to integer
  end subroutine Do_2D_MC_plain_FE


!-----------------------------------------------------------------------------
  ! Metropolis-Hastings algorithm
  ! Includes that the work function can vary with position
  subroutine Metropolis_Hastings_rectangle_v2(ndim, emit, df_out, F_out, pos_out)
    ! The interface is declared in the parent module
    integer, intent(in)                           :: ndim, emit
    double precision, intent(out)                 :: df_out, F_out
    double precision, intent(out), dimension(1:3) :: pos_out
    integer                                       :: count, i
    double precision                              :: rnd, alpha
    double precision, dimension(1:2)              :: std
    double precision, dimension(1:3)              :: cur_pos, new_pos, field
    double precision                              :: df_cur, df_new

    std(1:2) = emitters_dim(1:2, emit)*0.025d0/100.d0 ! Standard deviation for the normal distribution is 0.025% of the emitter length.
    ! This means that 68% of jumps are less than this value.
    ! The expected value of the absolute value of the normal distribution is std*sqrt(2/pi).

    ! Get a random initial position on the surface.
    ! We pick this location from a uniform distribution.
    count = 0
    do ! Infinite loop, we try to find a favorable position to start from
      CALL RANDOM_NUMBER(cur_pos(1:2))
      cur_pos(1:2) = cur_pos(1:2)*emitters_dim(1:2, emit) + emitters_pos(1:2, emit)
      cur_pos(3) = 0.0d0 ! On the surface

      ! Calculate the electric field at this position
      if (use_polarso .eqv. .true.) then
        field = PL_Calculate_Field_At(cur_pos)
      else
        field = Calc_Field_at(cur_pos)
      end if
      if (field(3) < 0.0d0) then
        exit ! We found a nice spot so we exit the loop
      else
        count = count + 1
        if (count > 10000) then ! The loop is infnite, must stop it at some point.
          print *, 'WARNING: MH was unable to find a favorable spot for emission!'
          exit
        end if
        ! In field emission it is rare the we reach the CL limit.
      end if
    end do

    F_out = field(3)

    ! Calculate the escape probability at this location
    if (field(3) < 0.0d0) then
      df_cur = Escape_Prob(field(3), cur_pos, emit)
    else
      df_cur = 0.0d0 ! Zero escape probabilty if field is not favorable
    end if

    !---------------------------------------------------------------------------
    ! We now pick a random distance and direction to jump to from our
    ! current location. We do this ndim times.
    do i = 1, ndim
      ! Find a new position using a normal distribution.
      !new_pos(1:2) = ziggurat_normal(cur_pos(1:2), std)
      new_pos(1:2) = box_muller(cur_pos(1:2), std)
      new_pos(3) = 0.0d0 ! At the surface

      ! Make sure that the new position is within the limits of the emitter area.
      call check_limits_metro_rec(new_pos, emit)

      ! Calculate the field at the new position
      if (use_polarso .eqv. .true.) then
        field = PL_Calculate_Field_At(new_pos)
      else
        field = Calc_Field_at(new_pos)
      end if

      ! Check if the field is favorable for emission at the new position.
      ! If it is not then cycle, i.e. we reject this location and
      ! pick another one.
      if (field(3) > 0.0d0) cycle ! Do the next loop iteration, i.e. find a new position.

      ! Calculate the escape probability at the new position, to compair with
      ! the current position.
      df_new = Escape_Prob(field(3), new_pos, emit)

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
        ! Jump to this position with probability alpha, i.e. if rnd is less than alpha
        if (rnd < alpha) then
          cur_pos = new_pos ! New position becomes the current position
          df_cur = df_new
          F_out = field(3)
        end if
      end if
    end do

    ! Return the current position

    pos_out = cur_pos
    df_out = df_cur
  end subroutine Metropolis_Hastings_rectangle_v2

  subroutine Metropolis_Hastings_rectangle_v2_field(ndim, emit, df_out, F_out, pos_out)
    ! The interface is declared in the parent module
    integer, intent(in)                           :: ndim, emit
    double precision, intent(out)                 :: df_out, F_out
    double precision, intent(out), dimension(1:3) :: pos_out

    double precision, dimension(1:3)              :: cur_field, new_field
    double precision, dimension(1:3)              :: cur_pos, new_pos
    double precision, dimension(1:2)              :: std
    double precision                              :: rnd, alpha
    integer                                       :: i, count

    std(1:2) = emitters_dim(1:2, emit)*0.075d0/100.d0 ! Standard deviation for the normal distribution is 0.075% of the emitter length.
    ! This means that 68% of jumps are less than this value.
    ! The expected value of the absolute value of the normal distribution is std*sqrt(2/pi).

    ! Get a random initial position on the surface.
    ! We pick this location from a uniform distribution.
    count = 0
    do ! Infinite loop, we try to find a favorable position to start from
      CALL RANDOM_NUMBER(cur_pos(1:2))
      cur_pos(1:2) = cur_pos(1:2)*emitters_dim(1:2, emit) + emitters_pos(1:2, emit)
      cur_pos(3) = 0.0d0 ! On the surface

      ! Calculate the electric field at this position
      if (use_polarso .eqv. .true.) then
        cur_field = PL_Calculate_Field_At(cur_pos)
      else
        cur_field = Calc_Field_at(cur_pos)
      end if
      if (cur_field(3) < 0.0d0) then
        exit ! We found a nice spot so we exit the loop
      else
        count = count + 1
        if (count > 10000) exit ! The loop is infnite, must stop it at some point.
        ! In field emission it is rare the we reach the CL limit.
      end if
    end do

    do i = 1, ndim
      ! Find a new position using a normal distribution.
      !new_pos(1:2) = ziggurat_normal(cur_pos(1:2), std)
      new_pos(1:2) = box_muller(cur_pos(1:2), std)
      new_pos(3) = 0.0d0 ! At the surface

      ! Make sure that the new position is within the limits of the emitter area.
      call check_limits_metro_rec(new_pos, emit)

      ! Calculate the field at the new position
      if (use_polarso .eqv. .true.) then
        new_field = PL_Calculate_Field_At(new_pos)
      else
        new_field = Calc_Field_at(new_pos)
      end if

      ! Check if the field is favorable for emission at the new position.
      ! If it is not then cycle, i.e. we reject this location and
      ! pick another one.
      if (new_field(3) > 0.0d0) cycle ! Do the next loop iteration, i.e. find a new position.

      ! Keep in mind that the field is negative
      ! -2 < -1 = True (More negative field is more favorable for emission)
      if (new_field(3) < cur_field(3)) then
        cur_pos = new_pos ! New position becomes the current position
        cur_field = new_field
      else
        ! Here we have some thing like -2 < -3
        ! so alpha = -2/-3 = 2/3 = 0.67
        alpha = new_field(3) / cur_field(3)
        CALL RANDOM_NUMBER(rnd)
        ! Jump to this position with probability alpha, i.e. if rnd is less than alpha
        if (rnd < alpha) then
          cur_pos = new_pos
          cur_field = new_field
        end if
      end if

    end do

    F_out = cur_field(3)
    df_out = Escape_Prob(F_out, cur_pos, emit)
    pos_out = cur_pos
  end subroutine Metropolis_Hastings_rectangle_v2_field

    !-----------------------------------------------------------------------------
  ! Metropolis-Hastings algorithm
  ! Includes that the work function can vary with position
  subroutine Metropolis_Hastings_rectangle_J(ndim_in, emit, df_out, F_out, pos_out)
    integer, intent(inout)                        :: ndim_in
    integer, intent(in)                           :: emit
    double precision, intent(out)                 :: df_out, F_out
    integer                                       :: ndim, ndim_first
    double precision, intent(out), dimension(1:3) :: pos_out
    integer                                       :: count, i
    double precision                              :: rnd, alpha
    double precision, dimension(1:2)              :: std
    double precision, dimension(1:3)              :: cur_pos, new_pos, field
    double precision                              :: sup_cur, sup_new ! Log of the electron supply (the chain target)

    ! Reset the accepted/rejected jump counters so the acceptance rate
    ! (and the adaptive MH_std it drives) is computed per chain rather than
    ! accumulating over the whole simulation.
    jump_a = 0
    jump_r = 0
    ndim_in = 0

    ndim = MH_ndim
    ndim_first = nint(ndim*MH_warmup_frac)
    std(1:2) = emitters_dim(1:2, emit)*MH_init_std ! Warmup std, 10% of the emitter size

    ! Get a random initial position on the surface.
    ! We pick this location from a uniform distribution.
    count = 0
    do ! Infinite loop, we try to find a favorable position to start from
      CALL RANDOM_NUMBER(cur_pos(1:2))
      cur_pos(1:2) = cur_pos(1:2)*emitters_dim(1:2, emit) + emitters_pos(1:2, emit)
      cur_pos(3) = 0.0d0 ! On the surface

      ! Calculate the electric field at this position
      if (use_polarso .eqv. .true.) then
        field = PL_Calculate_Field_At(cur_pos)
      else
        field = Calc_Field_at(cur_pos)
      end if

      if (field(3) < 0.0d0) then
        exit ! We found a nice spot so we exit the loop
      else
        count = count + 1
        if (count > 10000) then ! The loop is infnite, must stop it at some point.
        ! In field emission it is rare the we reach the CL limit.
          ! Mark the chain as failed: unfavorable field and zero escape
          ! probability, so the caller does not emit an electron from it.
          ndim_in = -1
          F_out = 1.0d0
          df_out = -huge(1.0d0)
          pos_out(1:2) = emitters_pos(1:2, emit)
          pos_out(3) = 0.0d0
          print *, 'Failed to find spot for emission'
          return ! Exit the function
        end if
      end if
    end do

    F_out = field(3)

    ! The chain samples positions proportional to the electron supply,
    ! see Elec_Supply_log.
    sup_cur = Elec_Supply_log(field(3), cur_pos, emit)

    !---------------------------------------------------------------------------
    ! We now pick a random distance and direction to jump to from our
    ! current location. We do this ndim times.
    do i = 1, ndim

      ! The first MH_warmup_frac of the jumps use the fixed initial std,
      ! after that the adapted step size. MH_std is only updated between
      ! chains (see below), so the proposal kernel is fixed within a chain.
      if (i > ndim_first) then
        std(1:2) = emitters_dim(1:2, emit)*MH_std ! Standard deviation for the normal distribution is MH_std of the emitter length.
        ! This means that 68% of jumps are less than this value.
        ! The expected value of the absolute value of the normal distribution is std*sqrt(2/pi).
      end if

      ! Find a new position using a normal distribution.
      new_pos(1:2) = box_muller(cur_pos(1:2), std)
      new_pos(3) = 0.0d0 ! At the surface

      ! Make sure that the new position is within the limits of the emitter area.
      call check_limits_metro_rec(new_pos, emit)

      ! Calculate the field at the new position
      if (use_polarso .eqv. .true.) then
        field = PL_Calculate_Field_At(new_pos)
      else
        field = Calc_Field_at(new_pos)
      end if

      ! Check if the field is favorable for emission at the new position.
      ! If it is not the supply there is zero, i.e. we reject this location
      ! and pick another one.
      if (field(3) >= 0.0d0) then
        if (i > ndim_first) then
          jump_r = jump_r + 1
        end if
        cycle ! Do the next loop iteration, i.e. find a new position.
      end if

      ! Calculate the electron supply at the new position, to compare with
      ! the current position.
      sup_new = Elec_Supply_log(field(3), new_pos, emit)

      alpha = sup_new - sup_cur

      if (sup_new >= sup_cur) then
        cur_pos = new_pos
        sup_cur = sup_new
        F_out = field(3)
        if (i > ndim_first) then
          jump_a = jump_a + 1
        end if
      else
        CALL RANDOM_NUMBER(rnd)
        if (log(rnd) <= alpha) then
          cur_pos = new_pos
          sup_cur = sup_new
          F_out = field(3)
          if (i > ndim_first) then
            jump_a = jump_a + 1
          end if
        else
          if (i > ndim_first) then
            jump_r = jump_r + 1
          end if
        end if
      end if
    end do

    ! The acceptance rate of this chain drives one update of the shared
    ! proposal step size.
    if (jump_a + jump_r > 0) then
      a_rate = DBLE(jump_a) / DBLE(jump_r + jump_a)
      call MH_std_update(a_rate)
    end if

    ! Return the current position.
    ! The caller decides the emission with the escape probability at this
    ! position; the chain target above deliberately does not include it.
    pos_out = cur_pos
    df_out = Escape_Prob_log(F_out, cur_pos, emit)
  end subroutine Metropolis_Hastings_rectangle_J

  ! ----------------------------------------------------------------------------
  ! Batched (lockstep) version of Metropolis_Hastings_rectangle_J.
  ! All M chains advance together: in every jump iteration the proposals of
  ! all chains are collected and the surface field is evaluated for the whole
  ! set in one call to Calc_Field_at_Batch (a single GPU kernel in OpenACC
  ! builds), instead of one O(N) field sum per chain per jump.
  !
  ! Differences from running Metropolis_Hastings_rectangle_J M times
  ! (enabled with mh_batch = .true. in the input file, off by default):
  !  * All chains sample the particle configuration from the start of the
  !    time step. In the serial version each chain sees the electrons that
  !    earlier chains emitted within the same step.
  !  * The adaptive step size MH_std gets one update per jump iteration,
  !    driven by the acceptance rate of that iteration across the whole
  !    batch, instead of one update per chain as in the serial version.
  !  * The random number stream is consumed in a different order, so
  !    individual positions differ; the sampled distribution is the same.
  subroutine Metropolis_Hastings_rectangle_J_batch(M, emit, df_out, F_out, pos_out)
    integer, intent(in)                              :: M ! Number of chains
    integer, intent(in)                              :: emit
    double precision, dimension(1:M), intent(out)    :: df_out, F_out
    double precision, dimension(1:3, 1:M), intent(out) :: pos_out

    integer                                       :: ndim, ndim_first
    integer                                       :: i, k, mc, count, n_act
    integer                                       :: it_a, it_r ! Accepted/rejected jumps of the current iteration
    double precision                              :: rnd, alpha, sup_new
    double precision, dimension(1:2)              :: std
    integer,          dimension(:),    allocatable :: act ! Chains active in the current batch
    double precision, dimension(:, :), allocatable :: cur_pos, w_pos
    double precision, dimension(:, :), allocatable :: w_field
    double precision, dimension(:),    allocatable :: sup_cur ! Log of the electron supply (the chain target)
    logical,          dimension(:),    allocatable :: ok ! Chain has a favorable position

    ! The working arrays go on the heap: M is the full electron supply of the
    ! time step, which has no upper bound that would make automatic (stack)
    ! arrays safe. They are deallocated automatically on return.
    allocate(act(1:M), cur_pos(1:3, 1:M), w_pos(1:3, 1:M), w_field(1:3, 1:M), sup_cur(1:M), ok(1:M))

    ndim = MH_ndim
    ndim_first = nint(ndim*MH_warmup_frac)
    std(1:2) = emitters_dim(1:2, emit)*MH_init_std ! 10% of emitter size

    !---------------------------------------------------------------------------
    ! Get random initial positions on the surface for all chains.
    ! We keep proposing uniformly distributed positions for the chains that
    ! have not yet found a spot with a favorable field.
    ok = .false.
    n_act = M
    do k = 1, M
      act(k) = k
    end do

    count = 0
    do while (n_act > 0)
      do k = 1, n_act
        CALL RANDOM_NUMBER(w_pos(1:2, k))
        w_pos(1:2, k) = w_pos(1:2, k)*emitters_dim(1:2, emit) + emitters_pos(1:2, emit)
        w_pos(3, k) = 0.0d0 ! On the surface
      end do

      call Calc_Field_at_Batch(n_act, w_pos(:, 1:n_act), w_field(:, 1:n_act))

      i = n_act ! Old count, reuse i as a temporary
      n_act = 0
      do k = 1, i
        mc = act(k)
        if (w_field(3, k) < 0.0d0) then
          ! We found a nice spot for this chain.
          ! The chain samples positions proportional to the electron
          ! supply, see Elec_Supply_log.
          cur_pos(:, mc) = w_pos(:, k)
          F_out(mc) = w_field(3, k)
          sup_cur(mc) = Elec_Supply_log(w_field(3, k), w_pos(:, k), emit)
          ok(mc) = .true.
        else
          n_act = n_act + 1
          act(n_act) = mc
        end if
      end do

      count = count + 1
      if ((count > 10000) .and. (n_act > 0)) then
        ! The loop is infinite, must stop it at some point.
        ! In field emission it is rare the we reach the CL limit.
        print *, 'Failed to find spot for emission'
        do k = 1, n_act
          mc = act(k)
          ! Mark the chain as failed: unfavorable field and zero
          ! escape probability, so no electron is emitted from it.
          F_out(mc) = 1.0d0
          sup_cur(mc) = -huge(1.0d0)
          cur_pos(1:2, mc) = emitters_pos(1:2, emit)
          cur_pos(3, mc) = 0.0d0
        end do
        exit
      end if
    end do

    !---------------------------------------------------------------------------
    ! We now pick a random distance and direction to jump to from the
    ! current location of each chain. We do this ndim times, advancing all
    ! chains together.
    do i = 1, ndim

      ! The first MH_warmup_frac of the jumps use the fixed initial std,
      ! after that the adapted step size. MH_std is updated once per jump
      ! iteration (below), so within an iteration all chains share the same
      ! fixed proposal kernel.
      if (i > ndim_first) then
        std(1:2) = emitters_dim(1:2, emit)*MH_std
      end if

      ! Collect the proposed positions of all active chains.
      n_act = 0
      do mc = 1, M
        if (ok(mc) .eqv. .false.) cycle
        n_act = n_act + 1
        act(n_act) = mc

        ! Find a new position using a normal distribution.
        w_pos(1:2, n_act) = box_muller(cur_pos(1:2, mc), std)
        w_pos(3, n_act) = 0.0d0 ! At the surface

        ! Make sure that the new position is within the limits of the emitter area.
        call check_limits_metro_rec(w_pos(:, n_act), emit)
      end do

      if (n_act == 0) exit ! All chains failed in the initialization

      ! Calculate the field at the new positions of all chains in one batch.
      call Calc_Field_at_Batch(n_act, w_pos(:, 1:n_act), w_field(:, 1:n_act))

      ! Do the Metropolis accept/reject step for each chain.
      it_a = 0
      it_r = 0
      do k = 1, n_act
        mc = act(k)

        ! Check if the field is favorable for emission at the new position.
        ! If it is not the supply there is zero, i.e. we reject this location.
        if (w_field(3, k) >= 0.0d0) then
          it_r = it_r + 1
          cycle
        end if

        ! Calculate the electron supply at the new position, to compare
        ! with the current position.
        sup_new = Elec_Supply_log(w_field(3, k), w_pos(:, k), emit)

        alpha = sup_new - sup_cur(mc)

        if (sup_new >= sup_cur(mc)) then
          cur_pos(:, mc) = w_pos(:, k)
          sup_cur(mc) = sup_new
          F_out(mc) = w_field(3, k)
          it_a = it_a + 1
        else
          CALL RANDOM_NUMBER(rnd)
          if (log(rnd) <= alpha) then
            cur_pos(:, mc) = w_pos(:, k)
            sup_cur(mc) = sup_new
            F_out(mc) = w_field(3, k)
            it_a = it_a + 1
          else
            it_r = it_r + 1
          end if
        end if
      end do

      ! The acceptance rate of this jump iteration across the whole batch
      ! drives one update of the shared proposal step size, mirroring the
      ! one update per chain of the serial version.
      if ((i > ndim_first) .and. (it_a + it_r > 0)) then
        a_rate = DBLE(it_a) / DBLE(it_a + it_r)
        call MH_std_update(a_rate)
      end if
    end do

    ! Return the current positions.
    ! The caller decides the emission with the escape probability at the
    ! final position of each chain; the chain target above deliberately
    ! does not include it. Failed chains keep their sentinels.
    pos_out = cur_pos
    do mc = 1, M
      if (ok(mc) .eqv. .true.) then
        df_out(mc) = Escape_Prob_log(F_out(mc), cur_pos(:, mc), emit)
      else
        df_out(mc) = -huge(1.0d0)
      end if
    end do
  end subroutine Metropolis_Hastings_rectangle_J_batch

  ! Adaptive MH in log space
  !subroutine Adpative_MH_log()
  !end subroutine Adpative_MH_log

  ! ----------------------------------------------------------------------------
  ! Checks the limits of the rectangular region of the emitter
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
    else if (par_pos(1) < x_min) then
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
    else if (par_pos(2) < y_min) then
      d_y = y_min - par_pos(2)
      par_pos(2) = d_y + y_min

      if(d_y > emitters_dim(2, emit)) then
        print *, 'Warning: d_y to large <'
        print *, d_y
      end if
    end if
  end subroutine check_limits_metro_rec

! ----------------------------------------------------------------------------
! --------------------- Ring Emitter -----------------------------------------
! ----------------------------------------------------------------------------
  subroutine Metropolis_Hastings_ring_J(ndim_in, emit, df_out, F_out, pos_out)
    integer, intent(inout)                        :: ndim_in
    integer, intent(in)                           :: emit
    double precision, intent(out)                 :: df_out, F_out
    integer                                       :: ndim, ndim_first
    double precision, intent(out), dimension(1:3) :: pos_out
    integer                                       :: count, i
    double precision                              :: rnd, alpha
    double precision, dimension(1:2)              :: std, prop
    double precision, dimension(1:3)              :: cur_pos, new_pos, field
    double precision                              :: sup_cur, sup_new ! Log of the electron supply (the chain target)
    double precision                              :: cur_angle, cur_radius, new_angle, new_radius
    double precision                              :: r_min, r_max

    ! Reset the accepted/rejected jump counters so the acceptance rate
    ! (and the adaptive MH_std it drives) is computed per chain rather than
    ! accumulating over the whole simulation.
    jump_a = 0
    jump_r = 0
    ndim_in = 0

    ndim = MH_ndim
    ndim_first = nint(ndim*MH_warmup_frac)

    r_min = emitters_ring(2, emit)                          ! Inner radius of the ring
    r_max = emitters_ring(2, emit) + emitters_ring(3, emit) ! Outer radius of the ring

    ! The chain moves in the (angle, radius) parametrization of the ring.
    ! Warmup std: MH_init_std of the full angle range and of the ring width.
    std(1) = 2.0d0*pi*MH_init_std
    std(2) = emitters_ring(3, emit)*MH_init_std

    ! Get a random initial position on the surface.
    ! We pick this location from a uniform distribution in (angle, radius).
    count = 0
    do ! Infinite loop, we try to find a favourable position to start from
      call random_number(prop)
      cur_angle = prop(1)*2.0d0*pi
      cur_radius = r_min + prop(2)*emitters_ring(3, emit)
      cur_pos(1) = cur_radius*cos(cur_angle)
      cur_pos(2) = cur_radius*sin(cur_angle)
      cur_pos(3) = 0.0d0 ! On the surface

      ! Calculate the electric field at this position
      if (use_polarso .eqv. .true.) then
        field = PL_Calculate_Field_At(cur_pos)
      else
        field = Calc_Field_at(cur_pos)
      end if

      if (field(3) < 0.0d0) then
        exit ! We found a nice spot so we exit the loop
      else
        count = count + 1
        if (count > 10000) then ! The loop is infnite, must stop it at some point.
        ! In field emission it is rare the we reach the CL limit.
          ! Mark the chain as failed: unfavorable field and zero escape
          ! probability, so the caller does not emit an electron from it.
          ndim_in = -1
          F_out = 1.0d0
          df_out = -huge(1.0d0)
          pos_out(1) = r_min
          pos_out(2) = 0.0d0
          pos_out(3) = 0.0d0
          print *, 'Failed to find spot for emission'
          return ! Exit the function
        end if
      end if
    end do

    F_out = field(3)

    ! The chain samples positions proportional to the electron supply.
    ! The log(radius) term is the Jacobian of the (angle, radius)
    ! parametrization of the ring, so the stationary distribution on the
    ! ring surface is proportional to the supply per unit area.
    sup_cur = Elec_Supply_log(field(3), cur_pos, emit) + log(cur_radius)

    !---------------------------------------------------------------------------
    ! We now pick a random distance and direction to jump to from our
    ! current location. We do this ndim times.
    do i = 1, ndim

      ! The first MH_warmup_frac of the jumps use the fixed initial std,
      ! after that the adapted step size. MH_std is only updated between
      ! chains (see below), so the proposal kernel is fixed within a chain.
      if (i > ndim_first) then
        std(1) = 2.0d0*pi*MH_std
        std(2) = emitters_ring(3, emit)*MH_std
      end if

      ! Find a new position using a normal distribution in (angle, radius).
      ! The angle wraps around the ring. Proposals outside the radial limits
      ! of the ring are rejected, the supply is zero off the emitter.
      prop = box_muller((/cur_angle, cur_radius/), std)
      new_angle = modulo(prop(1), 2.0d0*pi)
      new_radius = prop(2)

      if ((new_radius < r_min) .or. (new_radius > r_max)) then
        if (i > ndim_first) then
          jump_r = jump_r + 1
        end if
        cycle ! Do the next loop iteration, i.e. find a new position.
      end if

      new_pos(1) = new_radius*cos(new_angle)
      new_pos(2) = new_radius*sin(new_angle)
      new_pos(3) = 0.0d0 ! On the surface

      ! Calculate the field at the new position
      if (use_polarso .eqv. .true.) then
        field = PL_Calculate_Field_At(new_pos)
      else
        field = Calc_Field_at(new_pos)
      end if

      ! Check if the field is favourable for emission at the new position.
      ! If it is not the supply there is zero, i.e. we reject this location
      ! and pick another one.
      if (field(3) >= 0.0d0) then
        if (i > ndim_first) then
          jump_r = jump_r + 1
        end if
        cycle ! Do the next loop iteration, i.e. find a new position.
      end if

      ! Calculate the electron supply at the new position, to compare with
      ! the current position.
      sup_new = Elec_Supply_log(field(3), new_pos, emit) + log(new_radius)

      alpha = sup_new - sup_cur

      if (sup_new >= sup_cur) then
        cur_pos = new_pos
        cur_angle = new_angle
        cur_radius = new_radius
        sup_cur = sup_new
        F_out = field(3)
        if (i > ndim_first) then
          jump_a = jump_a + 1
        end if
      else
        CALL RANDOM_NUMBER(rnd)
        if (log(rnd) <= alpha) then
          cur_pos = new_pos
          cur_angle = new_angle
          cur_radius = new_radius
          sup_cur = sup_new
          F_out = field(3)
          if (i > ndim_first) then
            jump_a = jump_a + 1
          end if
        else
          if (i > ndim_first) then
            jump_r = jump_r + 1
          end if
        end if
      end if
    end do

    ! The acceptance rate of this chain drives one update of the shared
    ! proposal step size.
    if (jump_a + jump_r > 0) then
      a_rate = DBLE(jump_a) / DBLE(jump_r + jump_a)
      call MH_std_update(a_rate)
    end if

    ! Return the current position.
    ! The caller decides the emission with the escape probability at this
    ! position; the chain target above deliberately does not include it.
    pos_out = cur_pos
    df_out = Escape_Prob_log(F_out, cur_pos, emit)
  end subroutine Metropolis_Hastings_ring_J

  ! ----------------------------------------------------------------------------
  ! Checks the limits of the rectangular region of the emitter
  subroutine check_limits_metro_ring(par_pos, emit)
    double precision, dimension(1:3), intent(inout) :: par_pos
    integer, intent(in)                             :: emit
    double precision                                :: max_rad,min_rad,max_x,min_x,max_y,min_y,cur_rad,rnd


    max_rad = emitters_ring(1,emit)
    min_rad = emitters_ring(2,emit)

    max_x = max_rad
    min_x = -max_rad

    max_y = max_rad
    min_y = -max_rad

    cur_rad = sqrt(par_pos(1)**2 + par_pos(2)**2)


    if ((cur_rad > max_rad) .or. (cur_rad < min_rad)) then
      call random_number(rnd) ! Random number between 0 and 1 will tell us which component goes to max

      if (rnd <= 0.5) then ! x
        if (par_pos(1)>max_x) then
          par_pos(1) = max_x
        end if
        par_pos(2) = sqrt(max_rad**2 - par_pos(1)**2)
      else ! y
        if (par_pos(2)>max_y) then
          par_pos(2) = max_y
        end if
        par_pos(1) = sqrt(max_rad**2 - par_pos(2)**2)
      end if
    end if
  end subroutine check_limits_metro_ring

  ! Adaptive Metropolis algorithm is log space
  ! See: Haario, Heikki, Eero Saksman, and Johanna Tamminen. "An adaptive Metropolis algorithm." Bernoulli 7.2 (2001): 223-242.
  ! https://projecteuclid.org/download/pdf_1/euclid.bj/1080222083
  !subroutine Adaptive_MH_log()
    
  !end subroutine

  ! Test the field emission module.
  ! The optional argument v_y_ok reports whether the v_y check passed, so that
  ! the caller (mod_tests) can fold the result into the global pass/fail count.
  subroutine Test_Field_Emission_Module(v_y_ok)
    logical, intent(out), optional :: v_y_ok
    double precision :: results_scalar, F
    double precision, dimension(1:3) :: pos
    integer :: emit, sec
    logical :: ok
    ! Test v_y
    ! Test t_y
    ! Test Escape_Prob
    ! Test Elec_Supply_V2
    ! Test check_limits_metro_rec

    print *, 'Running tests for the field emission module'

    ! Set the work function
    ptr_Work_fun => FE_Test_Work_fun

    ! v_y(F, pos, emit)

    ! Test 1 for v_y. Normal test
    ! The field in a planar diode is E_z = -V/d, i.e. it points towards the
    ! cathode and is therefore negative. v_y must be evaluated with this
    ! (negative) field, matching how it is called from the emission routines.
    F = -2.0d3 / (1000.0d-9) ! E_z = -V/d, 2 kV across 1000 nm
    pos = (/ 0.0d0, 0.0d0, 0.0d0 /)
    emit = 1
    results_scalar = v_y(F, pos, emit)
    ok = (abs(results_scalar - 0.8253581935658024) < tolerance_abs)
    if (ok .eqv. .true.) then
      print *, 'v_y PASSED'
    else
      print *, 'v_y FAILED'
    end if

    if (present(v_y_ok)) v_y_ok = ok

  end subroutine Test_Field_Emission_Module

  ! Set a constant work function for tests
  double precision function FE_Test_Work_fun(pos, emit, sec)
    double precision, dimension(1:3), intent(in) :: pos
    integer, intent(in)                          :: emit
    integer, intent(out), optional               :: sec

    FE_Test_Work_fun = 4.7d0
  end function FE_Test_Work_fun

end Module mod_field_emission_v2
