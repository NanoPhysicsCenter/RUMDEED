!-------------------------------------------!
! Module for surface integration with CUBA  !
! http://www.feynarts.de/cuba/              !
! Kristinn Torfason                         !
!-------------------------------------------!
!
! Single home for the calls into the CUBA library. Every emission module
! integrates its electron supply over the unit square through Cuba_Integrate
! below. The integration method and the convergence criteria are the cuba_*
! variables in mod_global (settable in the input file), while the
! method-specific tuning knobs live here so they are defined once.
!
! The integrands must be in units of electrons supplied per time step, so
! that cuba_epsabs (default 0.5) means half an electron.

Module mod_cuba_integration
  use mod_global
  ! libcuba is C code compiled with 32-bit ints, so every integer that
  ! crosses into it must be integer(c_int) — including in ILP64 builds
  ! (POLARSO=yes compiles with 8-byte default integers). The integrand
  ! callbacks receive their integer arguments from C and must declare them
  ! integer(c_int) as well. In default builds c_int equals the default
  ! integer kind, so this changes nothing there.
  use, intrinsic :: iso_c_binding, only: c_int
  implicit none

  PRIVATE
  PUBLIC :: Cuba_Integrate

  ! Common Cuba arguments
  integer(c_int), parameter :: ndim = 2  ! Number of dimensions (emitter surfaces)
  integer(c_int), parameter :: ncomp = 1 ! Number of components in the integrand
  integer(c_int), parameter :: flags = 4 ! Bits 0-1: verbosity level 0.
                                  ! Bit 2 (Suave only): use only the last (largest)
                                  ! set of samples for the final result.
  integer(c_int), parameter :: seed = 0  ! 0 = Sobol quasi-random points (deterministic)
  character, parameter :: statefile = "" ! Empty string means don't save the state to a file
  integer(c_int), parameter   :: spin = -1 ! No spinning cores

  ! Suave specific
  integer(c_int), parameter   :: nnew = 1000 ! Number of integrand evaluations in each subdivision
  integer(c_int), parameter   :: nmin = 2    ! Minimum number of samples a former pass must contribute
                                             ! to a subregion to be considered in the region's compound
                                             ! integral value.
  double precision, parameter :: flatness = 25.0d0 ! Determines how prominently outliers, i.e. samples
                                             ! with a large fluctuation, figure in the total fluctuation,
                                             ! which in turn determines how a region is split up.
                                             ! Choose large for flat integrands, small for volatile ones.
                                             ! NOTE: flatness = 5 combined with a small nmin makes Suave
                                             ! converge to a biased result (~10% low even on smooth
                                             ! integrands, verified against known integrals); 25 is the
                                             ! Cuba demo value and passes the unit tests in mod_tests.

  ! Divonne specific
  integer(c_int), parameter :: key1 = 47 ! Sampling in the partitioning phase:
                                  ! 7, 9, 11, 13 selects the cubature rule of that degree
                                  ! (degree 11 only in 3D, degree 13 only in 2D). Other values
                                  ! use a quasi-random sample of |key1| points, Korobov if
                                  ! key1 > 0, Sobol/pseudo-random (by seed) if key1 < 0.
  integer(c_int), parameter :: key2 = 1  ! Sampling in the final integration phase. Same convention as
                                  ! key1; |key2| < 40 samples as many points as Divonne estimates
                                  ! it needs to reach the prescribed accuracy.
  integer(c_int), parameter :: key3 = -1 ! Strategy for the refinement phase: 0 = no further treatment,
                                  ! 1 = split the subregion once more, otherwise sample a third
                                  ! time with key3 as the sampling parameter (same as key2).
  integer(c_int), parameter :: maxpass = 2 ! Thoroughness of the partitioning phase: number of safety
                                  ! iterations without improvement before the partition is
                                  ! accepted as final.
  double precision, parameter :: border = 0.0d0 ! Width of the border extrapolation region; zero
                                  ! because the integrands can be evaluated on the boundary.
  double precision, parameter :: maxchisq = 10.0d0 ! Maximum chi-square a subregion may have in the
                                  ! final integration phase before it moves on to refinement.
  double precision, parameter :: mindeviation = 0.25d0 ! Fraction of the requested error below which
                                  ! a region that failed the chi-square test is not refined further.
  ! No list of suspected peak locations and no peak-finder subroutine are
  ! given. (If peak hints are ever added, they must be points in the unit
  ! square, NOT raw surface coordinates.)
  integer(c_int), parameter :: ngiven = 0, ldxgiven = ndim, nextra = 0

  ! Cuhre specific
  integer(c_int), parameter :: key = 0 ! Cubature rule: 0 uses the default, the degree-13 rule in 2D.

contains
  ! ----------------------------------------------------------------------------
  ! Integrate an emission integrand over the unit square with the method given
  ! by cuba_method, stopping at an estimated absolute error of
  ! max(cuba_epsabs, cuba_epsrel*|integral|).
  !
  ! integrand: Function with the Cuba integrand signature
  !            (ndim, xx, ncomp, ff, userdata [, nvec]), see the Cuba manual.
  !            When nvec > 1 it must accept blocks of points (the vectorized
  !            form with the extra nvec argument).
  ! nvec:      Maximum number of points handed to the integrand per call.
  ! userdata:  Passed through to the integrand (the emitter or section number).
  subroutine Cuba_Integrate(integrand, nvec, userdata, integral, error_out, prob, nregions, neval, fail)
    integer(c_int), external :: integrand
    integer, intent(in)  :: nvec
    integer, intent(in)  :: userdata
    double precision, dimension(1:ncomp), intent(out) :: integral  ! The integral over the unit square
    double precision, dimension(1:ncomp), intent(out) :: error_out ! The presumed absolute error
    double precision, dimension(1:ncomp), intent(out) :: prob      ! The chi-square probability
    integer, intent(out) :: nregions ! The actual number of subregions needed
    integer, intent(out) :: neval    ! The actual number of integrand evaluations needed
    integer, intent(out) :: fail     ! 0 = success, -1 = dimension out of range,
                                     ! > 0 = accuracy goal not met

    ! c_int copies of everything passed by reference into libcuba (the
    ! dummies above keep the default kind so callers need no changes)
    integer(c_int) :: nvec_c, userdata_c, mineval_c, maxeval_c
    integer(c_int) :: nregions_c, neval_c, fail_c

    nvec_c     = int(nvec,         kind=c_int)
    userdata_c = int(userdata,     kind=c_int)
    mineval_c  = int(cuba_mineval, kind=c_int)
    maxeval_c  = int(cuba_maxeval, kind=c_int)

    select case (cuba_method)
    case (cuba_method_suave)
      call suave(ndim, ncomp, integrand, userdata_c, nvec_c, &
                 cuba_epsrel, cuba_epsabs, flags, seed, &
                 mineval_c, maxeval_c, nnew, nmin, flatness, &
                 statefile, spin, &
                 nregions_c, neval_c, fail_c, integral, error_out, prob)
    case (cuba_method_divonne)
      call divonne(ndim, ncomp, integrand, userdata_c, nvec_c, &
                   cuba_epsrel, cuba_epsabs, flags, seed, &
                   mineval_c, maxeval_c, &
                   key1, key2, key3, maxpass, &
                   border, maxchisq, mindeviation, &
                   ngiven, ldxgiven, 0, nextra, 0, &
                   statefile, spin, &
                   nregions_c, neval_c, fail_c, integral, error_out, prob)
    case (cuba_method_cuhre)
      call cuhre(ndim, ncomp, integrand, userdata_c, nvec_c, &
                 cuba_epsrel, cuba_epsabs, flags, &
                 mineval_c, maxeval_c, key, &
                 statefile, spin, &
                 nregions_c, neval_c, fail_c, integral, error_out, prob)
    case default
      print '(a)', 'RUMDEED: ERROR UNKNOWN INTEGRATION METHOD'
      print *, cuba_method
      stop
    end select

    nregions = nregions_c
    neval    = neval_c
    fail     = fail_c

    ! A positive fail just means the accuracy goal was not met within
    ! cuba_maxeval evaluations. Only warn if the result missed the requested
    ! tolerance by more than 5%; a negative fail is always an error.
    if (fail /= 0) then
      if ((fail < 0) .or. &
          (error_out(1) > 1.05d0*max(cuba_epsabs, cuba_epsrel*abs(integral(1))))) then
        print '(a)', 'RUMDEED: WARNING Cuba did not return 0'
        print *, 'method = ', cuba_method
        print *, 'fail = ', fail
        print *, 'userdata = ', userdata
        print *, 'nregions = ', nregions
        print *, 'neval = ', neval, ' max is ', cuba_maxeval
        print *, 'integral(1) = ', integral(1)
        print *, 'error(1) = ', error_out(1)
        print *, 'epsabs = ', cuba_epsabs
        print *, 'integral(1)*epsrel = ', integral(1)*cuba_epsrel
        print *, 'prob(1) = ', prob(1)
        call Flush_Data()
      end if
    end if
  end subroutine Cuba_Integrate
end module mod_cuba_integration
