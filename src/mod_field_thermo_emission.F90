!-------------------------------------------!
! Test Module for field+thermo emission     !
!                                           !
! Kristinn Torfason                         !
! 08.11.19                                  !
!-------------------------------------------!

Module mod_field_thermo_emission
  use mod_global
  use mod_verlet
  use mod_pair
  use ziggurat
  use mod_work_function
  use mod_kevin_rjgtf
  !use mod_mc_integration
  implicit none

  PRIVATE
  PUBLIC :: Init_Field_Thermo_Emission, Clean_Up_Field_Thermo_Emission

  ! ----------------------------------------------------------------------------
  ! Variables
  integer, dimension(:), allocatable :: nrEmitted_emitters
  integer                            :: nrElecEmitAll
  !integer                            :: nrEmitted
  double precision, dimension(1:3)   :: F_avg = 0.0d0
  integer, parameter                 :: N_MH_step = 10 ! Number of steps to do in the MH algorithm
  double precision                   :: residual = 0.0d0 ! Should be a array the size of the number of emitters

  ! Constant used in MC integration (function Elec_Supply_V2)
  double precision :: time_step_div_q0

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

contains

  !-----------------------------------------------------------------------------
  ! Initialize the Field Emission
subroutine Init_Field_Thermo_Emission()
    ! Allocate the number of emitters
    allocate(nrEmitted_emitters(1:nrEmit))

    ! Initialize variables
    nrEmitted_emitters = 0 ! Number of electrons emitted from emitter number i in this time step

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_ElecHole_Planar

    ! Function for the electric field in the system
    ptr_field_E => field_E_planar

    ! The function that does the emission
    ptr_Do_Emission => Do_Field_Thermo_Emission

    ! The function to do image charge effects
    ptr_Image_Charge_effect => Force_Image_charges_v2

    call Read_work_function()

    ! Parameters used in the module
    time_step_div_q0 = time_step / q_0

    ! Initialize the Ziggurat algorithm
    call zigset(my_seed(1))

  end subroutine Init_Field_Thermo_Emission

  subroutine Clean_Up_Field_Thermo_Emission()
    deallocate(nrEmitted_emitters)

    call Work_fun_cleanup()
  end subroutine Clean_Up_Field_Thermo_Emission

  !-----------------------------------------------------------------------------
  ! This subroutine gets called from main when the emitters should emit the electrons
  subroutine Do_Field_Thermo_Emission(step)
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
          !!!!!!!call Do_Field_Emission_Planar_rectangle(step, i)
          call Do_Field_Thermo_Emission_Planar_simple(step, i)
        !case default
        !  print *, 'Vacuum: WARNING unknown emitter type!!'
        !  print *, emitters_type(i)
        !end select
      end if
    end do

    cur_time = time_step * step / time_scale ! Scaled in units of time_scale
    write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, nrElecEmitAll, nrElec, &
                                                       & (nrEmitted_emitters(i), i = 1, nrEmit)
  end subroutine Do_Field_Thermo_Emission


  subroutine Do_Field_Thermo_Emission_Planar_simple(step, emit)
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

    ! Use the number electrons as an average for a Poission distribution to get the number of electrons to emitt 
    N_round = Rand_Poission(N_sup)
    !N_round = nint(N_sup + residual) ! Round to whole number
    !residual = N_sup - N_round


    nrElecEmit = 0
    ! Loop over all electrons and place them
    do i = 1, N_round
      !call Metropolis_Hastings_rectangle_v2_field(N_MH_step, emit, D_f, F, par_pos)
      call Metropolis_Hastings_rectangle_v2(N_MH_step, emit, D_f, F, par_pos)

      par_pos(3) = 1.0d0*length_scale
      par_vel = 0.0d0 ! Get velocity from MB dist!!!!
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
  end subroutine Do_Field_Thermo_Emission_Planar_simple

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
      field = Calc_Field_at(cur_pos)
      if (field(3) < 0.0d0) then
        exit ! We found a nice spot so we exit the loop
      else
        count = count + 1
        if (count > 10000) exit ! The loop is infnite, must stop it at some point.
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
    !df_out = Escape_Prob(F_out, cur_pos)
    df_out = 0.0d0
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
    double precision                 :: A, w_theta ! Emitter area

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
      w_theta = w_theta_xy(par_pos, userdata)
      ff(1) = Get_Kevin_Jgtf(field(3), T_temp, w_theta)
      !ff(1) = Elec_Supply_V2(field(3), par_pos, userdata)*Escape_Prob(field(3), par_pos, userdata)
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
     N_sup = integral(1) * time_step_div_q0

     ! Finish calculating the average field
     F_avg = F_avg / neval
  end subroutine Do_Cuba_Suave_Simple



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
  double precision function Elec_Supply_V2(F, pos, emit)
    double precision, dimension(1:3), intent(in) :: pos
    double precision,                 intent(in) :: F
    integer, intent(in)                          :: emit

    Elec_Supply_V2 = time_step_div_q0 * a_FN/(t_y(F, pos, emit)**2*w_theta_xy(pos, emit)) * F**2
  end function Elec_Supply_V2

end module mod_field_thermo_emission