!-------------------------------------------!
! Submodule for Monte Carlo Integration     !
! routines for field emission               !
! Kristinn Torfason                         !
! 21.10.18                                  !
!-------------------------------------------!
submodule (mod_field_emission_v2) smod_mc_integration
  use mod_global
  use mod_verlet

  implicit none
contains
  ! ----------------------------------------------------------------------------
  ! This function is called to do the surface integration.
  ! Here we select the method to do it.
  !
  module subroutine Do_Surface_Integration_FE(emit, N_sup)
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
      ff(1) = Elec_Supply_V2(field(3), par_pos)
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
     end if

     !! Round the results to the nearest integer
     !N_sup = nint( integral(1) )
     N_sup = integral(1)

     ! Finish calculating the average field
     F_avg = F_avg / neval
  end subroutine Do_Cuba_Suave_FE

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
        e_sup_res = Elec_Supply_V2(field(3), par_pos)
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

end submodule smod_mc_integration