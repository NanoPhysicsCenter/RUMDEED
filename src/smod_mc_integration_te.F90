!-------------------------------------------!
! Submodule for Monte Carlo Integration     !
! routines for thermionic emission          !
! Kristinn Torfason                         !
! 21.11.18                                  !
!-------------------------------------------!
submodule (mod_therminoic_emission) smod_mc_integration_te
  use mod_global
  use mod_verlet

  implicit none
contains
  ! ----------------------------------------------------------------------------
  ! This function is called to do the surface integration.
  ! Here we select the method to do it.
  !
  module subroutine Do_Surface_Integration_TE(emit, N_sup)
    integer, intent(in)  :: emit ! The emitter to do the integration on
    integer, intent(out) :: N_sup ! Number of electrons

    !call Do_2D_MC_plain(emit, N_sup)
    call Do_Cuba_Suave_TE(emit, N_sup)
  end subroutine Do_Surface_Integration_TE

  ! ----------------------------------------------------------------------------
  ! The integration function for the Cuba library
  !
  integer function integrand_cuba_te(ndim, xx, ncomp, ff, userdata)
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
    par_pos(1) = emitters_pos(1, userdata) + xx(1)*emitters_dim(1, userdata) ! x and y position on the surface
    par_pos(2) = emitters_pos(2, userdata) + xx(2)*emitters_dim(2, userdata) ! x and y position on the surface
    par_pos(3) = 0.0d0 ! Height, i.e. on the surface

    ! Calculate the electric field on the surface
    field = Calc_Field_at(par_pos)

    ! Add to the average field
    F_avg = F_avg + field

    ! Check if the field is favourable for emission
    if (field(3) < 0.0d0) then
      ! The field is favourable for emission
      ff(1) = Richardson_Dushman(field(3), par_pos)
    else
      ! The field is NOT favourable for emission
      ! This point does not contribute
      ff(1) = 0.0d0
    end if

    ! We mutiply with the area of the emitter because Cuba does the 
    ! integration over the hybercube, i.e. from 0 to 1.
    ff(1) = A*ff(1)
    
    integrand_cuba_te = 0 ! Return value to Cuba, 0 = success
  end function integrand_cuba_te

  ! ----------------------------------------------------------------------------
  ! Use the Cuba library to do the surface integration
  ! http://www.feynarts.de/cuba/
  !
  subroutine Do_Cuba_Suave_TE(emit, N_sup)
    ! Input / output variables
    integer, intent(in)  :: emit
    integer, intent(out) :: N_sup

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

    call suave(ndim, ncomp, integrand_cuba_te, userdata, nvec, &
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

     ! Round the results to the nearest integer
     N_sup = nint( integral(1) )

     ! Finish calculating the average field
     F_avg = F_avg / neval
  end subroutine Do_Cuba_Suave_TE

end submodule smod_mc_integration_te