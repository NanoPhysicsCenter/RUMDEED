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
  implicit none

  PRIVATE
  PUBLIC :: Init_Field_Emission_v2, Clean_Up_Field_Emission_v2, t_y, v_y, Escape_Prob, F_avg, Elec_Supply_V2

  ! ----------------------------------------------------------------------------
  ! Variables
  integer, dimension(:), allocatable :: nrEmitted_emitters
  integer                            :: nrElecEmitAll
  !integer                            :: nrEmitted
  double precision, dimension(1:3)   :: F_avg = 0.0d0
  integer, parameter                 :: N_MH_step = 10*2 ! Number of steps to do in the MH algorithm
  double precision                   :: residual = 0.0d0 ! Should be a array the size of the number of emitters

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

    ! Initialize variables
    nrEmitted_emitters = 0 ! Number of electrons emitted from emitter number i in this time step

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_ElecHole_Planar

    ! Function for the electric field in the system
    ptr_field_E => field_E_planar

    ! The function that does the emission
    ptr_Do_Emission => Do_Field_Emission

    ! The function to do image charge effects
    ptr_Image_Charge_effect => Force_Image_charges_v2

    call Read_work_function()

    ! Parameters used in the module
    time_step_div_q0 = time_step / q_0

    ! Initialize the Ziggurat algorithm
    !call zigset(my_seed(1))

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
          !call Do_Field_Emission_Planar_simple(step, i)
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
    call Do_Cuba_Suave_Simple(emit, N_sup)

    N_round = nint(N_sup + residual) ! Round to whole number
    residual = N_sup - N_round

    ! Loop over all electrons and place them
    do i = 1, N_round
      call Metropolis_Hastings_rectangle_v2(N_MH_step, emit, D_f, F, par_pos)
    

      par_pos(3) = 1.0d0*length_scale
      par_vel = 0.0d0
      rnd = w_theta_xy(par_pos, emit, sec) ! Get the section

      ! Add a particle to the system
      call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1, sec)

      nrElecEmit = nrElecEmit + 1
      nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
    end do

    Df_avg = 0.0d0

    write (ud_field, "(i8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8)", iostat=IFAIL) &
                                      step, F_avg(1), F_avg(2), F_avg(3), N_sup, df_avg

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
    integer                          :: N_round

    ! Emission variables
    double precision                 :: D_f, Df_avg, F, rnd
    integer                          :: s, sec, IFAIL
    integer                          :: nrElecEmit
    double precision, dimension(1:3) :: par_pos, par_vel

    call Do_Surface_Integration_FE(emit, N_sup)
    N_round = nint(N_sup + residual)
    residual = N_sup - N_round

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

      call Metropolis_Hastings_rectangle_v2(N_MH_step, emit, D_f, F, par_pos)
      !call Metropolis_Hastings_rectangle_v2_field(N_MH_step, emit, D_f, F, par_pos)
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

          nrElecEmit = nrElecEmit + 1
          nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
        !!$OMP END CRITICAL(EMIT_PAR)
      end if
    end do
    !!$OMP END PARALLEL DO

    df_avg = df_avg / N_sup

    write (ud_field, "(i8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, E16.8)", iostat=IFAIL) &
                                      step, F_avg(1), F_avg(2), F_avg(3), N_sup, df_avg

    nrElecEmitAll = nrElecEmitAll + nrElecEmit
    !nrEmitted = nrEmitted + nrElecEmit
  end subroutine Do_Field_Emission_Planar_rectangle

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
  double precision function Escape_Prob(F, pos, emit)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos
    integer, intent(in)                          :: emit

    Escape_Prob = exp(b_FN * (sqrt(w_theta_xy(pos, emit)))**3 * v_y(F, pos, emit) / (-1.0d0*F))

  end function Escape_Prob

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
  ! Here we select the method to do it.
  !
  subroutine Do_Surface_Integration_FE(emit, N_sup)
    integer, intent(in)           :: emit ! The emitter to do the integration on
    double precision, intent(out) :: N_sup ! Number of electrons

    !call Do_2D_MC_plain(emit, N_sup)
    call Do_Cuba_Suave_FE(emit, N_sup)
  end subroutine Do_Surface_Integration_FE

  ! ----------------------------------------------------------------------------
  ! The integration function for the Cuba library
  !
  integer function integrand_cuba_fe(ndim, xx, ncomp, ff, userdata)
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
    field = Calc_Field_at(par_pos)

    ! Add to the average field
    F_avg = F_avg + field

    ! Check if the field is favourable for emission
    if (field(3) < 0.0d0) then
      ! The field is favourable for emission
      ! Calculate the electron supply at this point
      ff(1) = Elec_Supply_V2(field(3), par_pos, userdata)
    else
      ! The field is NOT favourable for emission
      ! This point does not contribute
      ff(1) = 0.0d0
    end if

    ! We mutiply with the area of the emitter because Cuba does the 
    ! integration over the hybercube, i.e. from 0 to 1.
    ff(1) = A*ff(1)
    
    integrand_cuba_fe = 0 ! Return value to Cuba, 0 = success
  end function integrand_cuba_fe

  ! ----------------------------------------------------------------------------
  ! Use the Cuba library to do the surface integration
  ! http://www.feynarts.de/cuba/
  !
  subroutine Do_Cuba_Suave_FE(emit, N_sup)
    ! Input / output variables
    integer, intent(in)           :: emit
    double precision, intent(out) :: N_sup
    integer                       :: IFAIL

    ! Cuba integration variables
    integer, parameter :: ndim = 2 ! Number of dimensions
    integer, parameter :: ncomp = 1 ! Number of components in the integrand
    integer            :: userdata = 0 ! User data passed to the integrand
    integer, parameter :: nvec = 1 ! Number of points given to the integrand function
    double precision   :: epsrel = 1.0d-10 ! Requested relative error
    double precision   :: epsabs = 0.5d-1 ! Requested absolute error
    integer            :: flags = 0+4 ! Flags
    integer            :: seed = 0 ! Seed for the rng. Zero will use Sobol.
    integer            :: mineval = 100000 ! Minimum number of integrand evaluations
    integer            :: maxeval = 10000000 ! Maximum number of integrand evaluations
    integer            :: nnew = 1250 ! Number of integrand evaluations in each subdivision
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
    F_avg = 0.0d0

    ! Pass the number of the emitter being integraded over to the integrand as userdata
    userdata = emit

    call suave(ndim, ncomp, integrand_cuba_fe, userdata, nvec, &
     & epsrel, epsabs, flags, seed, &
     & mineval, maxeval, nnew, nmin, flatness, &
     & statefile, spin, &
     & nregions, neval, fail, integral, error, prob)

     if (fail /= 0) then
      print '(a)', 'Vacuum: WARNING Cuba did not return 0'
      print *, fail
      print *, error
      print *, prob
      call Flush_all_files()
     end if


     !! Round the results to the nearest integer
     !N_sup = nint( integral(1) )
     N_sup = integral(1)

     ! Finish calculating the average field
     F_avg = F_avg / neval

     ! Write the output variables of the integration to a file
     write(ud_integrand, '(i3, tr2, i8, tr2, i8, tr2, i4, tr2, ES12.4, tr2, ES12.4, tr2, ES12.4)', iostat=IFAIL) &
                          & emit, nregions, neval, fail, integral(1), error(1), prob(1)
  end subroutine Do_Cuba_Suave_FE

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
    field = Calc_Field_at(par_pos)

    ! Add to the average field
    F_avg = F_avg + field

    ! Check if the field is favourable for emission
    if (field(3) < 0.0d0) then
      ! The field is favourable for emission
      ! Calculate the current density at this point
      ff(1) = Elec_Supply_V2(field(3), par_pos, userdata) * Escape_Prob(field(3), par_pos, userdata)
    else
      ! The field is NOT favourable for emission
      ! This point does not contribute
      ff(1) = 0.0d0
    end if

    ! We mutiply with the area of the emitter because Cuba does the 
    ! integration over the hybercube, i.e. from 0 to 1.
    ff(1) = A*ff(1)
    
    integrand_cuba_simple = 0 ! Return value to Cuba, 0 = success
  end function integrand_cuba_simple

  subroutine Do_Cuba_Suave_Simple(emit, N_sup)
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
    F_avg = 0.0d0

    ! Pass the number of the emitter being integraded over to the integrand as userdata
    userdata = emit

    call suave(ndim, ncomp, integrand_cuba_simple, userdata, nvec, &
     & epsrel, epsabs, flags, seed, &
     & mineval, maxeval, nnew, nmin, flatness, &
     & statefile, spin, &
     & nregions, neval, fail, integral, error, prob)

     if (fail /= 0) then
      print '(a)', 'Vacuum: WARNING Cuba did not return 0'
      print *, fail
      print *, error
      print *, prob
     end if

     !! Round the results to the nearest integer
     !N_sup = nint( integral(1) )
     N_sup = integral(1)

     ! Finish calculating the average field
     F_avg = F_avg / neval
  end subroutine

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
    F_avg = 0.0d0 ! The average field on the surface
    Nmc_try = 0

    do
      ! Get a random position on the emitter
      CALL RANDOM_NUMBER(par_pos(1:2))
      par_pos(1:2) = emitters_pos(1:2, emit) + par_pos(1:2)*emitters_dim(1:2, emit)
      par_pos(3) = 0.0d0 ! z = 0, emitter surface

      ! Calculate the field on the emitter surface
      field = Calc_Field_at(par_pos)

      ! Check if the field is favourable for emission
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
          print *, 'Vacuum: Warning MC integration taking to long, stoping it'
          print *, 'mc_err = ', mc_err
          !print *, 'step = ', step
          exit
        end if

      else ! field(3) < 0.0d0
        Nmc_try = Nmc_try + 1

        ! Stop if we are taking to long to find a favourable point.
        ! This should be rare in field emission.
        if (Nmc_try > 1000) then
          print *, 'Vacuum: Warning to many field attempts at finding a favourable location in MC integration'
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
    do ! Infinite loop, we try to find a favourable position to start from
      CALL RANDOM_NUMBER(cur_pos(1:2))
      cur_pos(1:2) = cur_pos(1:2)*emitters_dim(1:2, emit) + emitters_pos(1:2, emit)
      cur_pos(3) = 0.0d0 ! On the surface

      ! Calculate the electric field at this position
      field = Calc_Field_at(cur_pos)
      if (field(3) < 0.0d0) then
        exit ! We found a nice spot so we exit the loop
      else
        count = count + 1
        if (count > 10000) then ! The loop is infnite, must stop it at some point.
          print *, 'WARNING: MH was unable to find a favourable spot for emission!'
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
      df_cur = 0.0d0 ! Zero escape probabilty if field is not favourable
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
      field = Calc_Field_at(new_pos)

      ! Check if the field is favourable for emission at the new position.
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
    do ! Infinite loop, we try to find a favourable position to start from
      CALL RANDOM_NUMBER(cur_pos(1:2))
      cur_pos(1:2) = cur_pos(1:2)*emitters_dim(1:2, emit) + emitters_pos(1:2, emit)
      cur_pos(3) = 0.0d0 ! On the surface

      ! Calculate the electric field at this position
      cur_field = Calc_Field_at(cur_pos)
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
      new_field = Calc_Field_at(new_pos)

      ! Check if the field is favourable for emission at the new position.
      ! If it is not then cycle, i.e. we reject this location and
      ! pick another one.
      if (new_field(3) > 0.0d0) cycle ! Do the next loop iteration, i.e. find a new position.

      ! Keep in mind that the field is negative
      ! -2 < -1 = True (More negative field is more favourable for emission)
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

      !if(d_x > emitters_dim(1, emit)) then
      !  print *, 'Warning: d_x to large >'
      !  print *, d_x
      !end if
    else if (par_pos(1) < x_min) then
      d_x = x_min - par_pos(1)
      par_pos(1) = d_x + x_min

      !if(d_x > emitters_dim(1, emit)) then
      !  print *, 'Warning: d_x to large <'
      !  print *, d_x
      !end if
    end if

    !Check y ----------------------------------------
    if (par_pos(2) > y_max) then
      d_y = par_pos(2) - y_max
      par_pos(2) = y_max - d_y

      !if(d_y > emitters_dim(2, emit)) then
      !  print *, 'Warning: d_y to large >'
      !  print *, d_y
      !end if
    else if (par_pos(2) < y_min) then
      d_y = y_min - par_pos(2)
      par_pos(2) = d_y + y_min

      !if(d_y > emitters_dim(2, emit)) then
      !  print *, 'Warning: d_x to large <'
      !  print *, d_y
      !end if
    end if
  end subroutine check_limits_metro_rec

end Module mod_field_emission_v2
