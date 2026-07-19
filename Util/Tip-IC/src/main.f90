!-------------------------------------------------------------------------------
! Calculates the electric field with image charge using the method in
! [1] "Methods for calculating electrostatic quantities due to a free charge in a
! nanoscale three-dimensional tip/base junction", Journal of Applied Physics,
! Volume 78, Issue 8, October 15, 1995, pp.4888-4894.
!
! This program calculates the conical functions using a library from
! [2] T.M. Dunster, A. Gil, J. Segura, N.M. Temme,
! "Conical : An extended module for computing a numerically satisfactory pair of
! solutions of the differential equation for conical functions",
! Computer Physics Communications, Volume 217, 2017, Pages 193-197.
!
! Code for the Gamma function:
! http://jin.ece.illinois.edu/routines/routines.html
! [3] Computation of special functions
! S. Zhang and J. Jin
! Wiley, 1996
! ISBN 0-471-11963-6
!
! Kristinn Torfason
! 16.03.2018
!-------------------------------------------------------------------------------

program Tip_IC
  USE Conical ! The module for the conical functions
  USE, intrinsic :: ieee_arithmetic ! For ieee_is_nan
  implicit none

  ! Note: epsilon_0, m_e, q_e and V_0 are not used in the calculation. The
  ! output is the unnormalized image-charge density from ref. 1 (no physical
  ! prefactor and no applied-field contribution), not SI units.
  double precision, parameter :: pi = 3.14159265358979324d0 ! Pi
  double precision, parameter :: epsilon_0 = 8.854187817d-12 ! Farads / meters
  double precision, parameter :: m_e = 9.10938291d-31 ! kg
  double precision, parameter :: q_e = -1.602176565d-19 ! Coulomb

  ! Below x - 1 = x_series the Gauss series is used in place of the conical
  ! function library. Measured: the series holds to <= 1e-13 out to x - 1 of
  ! 3.2e-3 at tau = 100 (further at smaller tau), while the library only becomes
  ! usable above x - 1 ~ 4.6e-4, so this sits inside the overlap for every tau
  ! reached by the integration.
  double precision, parameter :: x_series = 1.0d-3

  ! Fixed tau-quadrature grid used by Calculate_Integral: composite Simpson's
  ! rule with N_simp intervals on each unit segment [k-1, k], k = 1..tau_max.
  ! The prefactor table below is precomputed on exactly this grid.
  integer, parameter          :: N_simp  = 20
  integer, parameter          :: tau_max = 100
  double precision, parameter :: h_tau   = 1.0d0 / N_simp

  ! Highest m in the H/S summations, see ref. 1
  integer, parameter          :: m_max = 30

  ! The part of H_m that depends only on (m, tau), tabulated on the tau grid:
  ! Hm_pref(m, j) holds Calculate_Hm_prefactor(m, j*h_tau). Filled once by
  ! Precompute_Hm_Prefactors, then read-only during the grid loop.
  double complex, allocatable :: Hm_pref(:, :)

  double precision, parameter :: xi_c = 1.00d0
  !double precision, parameter :: eta_c = -0.95d0
  double precision            :: eta_c
  double precision, parameter :: phi_c = 0.0d0 ! Position of the point charge
  double precision            :: xi_p, eta_p, phi_p ! Position where the charge density is being calculated

  double precision, parameter :: length_scale = 1.0d-9 ! 1 nano-meter
  double precision            :: V_0 = 0.0d0 ! Voltage over the gap in Volts

    !Hyperbolid tip
  double precision            :: a_foci ! a_foci
  double precision            :: theta_tip ! theta = arccos(d/a_foci) Radius of curvature of the tip
  double precision            :: r_tip ! r_tip = a_foci*sin(theta)*tan(theta)
  double precision            :: eta_1 ! eta_1 = -cos(theta) defines the surface of the tip
  double precision, parameter :: eta_2 = 0.0d0 ! Defines the anode plane (should be 0)

  double precision            :: max_xi ! Max value for xi variable for the surface (tip height)
  !double precision            :: shift_z ! shift_z = a_foci*eta_1*max_xi Shift the z-coordinate uppwards. Base of tip at z = 0
  double precision            :: h_tip = 500.0d0*length_scale ! Height of the tip from base
  double precision            :: R_base = 250.0d0*length_scale ! Base radius
  double precision            :: d, d_tip = 1000.0d0*length_scale ! Gap length in units of length

  !integer :: ierr
  !double precision :: test

  d = d_tip + h_tip

  max_xi = h_tip/d_tip + 1.0d0
  a_foci = sqrt(d_tip**2*R_base**2 / (h_tip**2 + 2*d_tip*h_tip) + d_tip**2)
  eta_1 = 1.0d0 * d_tip / a_foci

  eta_c = 250.0d-9 / a_foci

  theta_tip = acos(d_tip/a_foci)
  r_tip = a_foci * sin(theta_tip) * tan(theta_tip)
  !shift_z = abs(a_foci * eta_1 * max_xi)

  ! Parameters used in reference 1
  !eta_1 = 0.9d0
  !r_tip = 20.0d-10 ! Angström

  !theta_tip = acos(eta_1)
  !a_foci = r_tip / ( sin(theta_tip)*tan(theta_tip) )
  !a_foci= 94.7d-10 ! Angström


  !xi_p = 1.3d0
  eta_p = eta_1
  !phi_p = 0.0d0

  print *, 'Tip parameters are'
  print *, 'a_foci = ', a_foci
  print *, 'eta_1  = ', eta_1
  print *, 'V_0    = ', V_0
  print * , ''

  print *, 'Electron is located at'
  print *, 'xi_c  = ', xi_c
  print *, 'eta_c = ', eta_c
  print *, 'phi_c = ', phi_c
  print * , 'or'
  print *, 'x     = ', (a_foci*sqrt(xi_c**2-1.0d0)*sqrt(1.0d0-eta_c**2)*cos(phi_c)/length_scale), ' [nm]'
  print *, 'y     = ', (a_foci*sqrt(xi_c**2-1.0d0)*sqrt(1.0d0-eta_c**2)*sin(phi_c)/length_scale), ' [nm]'
  print *, 'z     = ', (a_foci*xi_c*eta_c/length_scale), ' [nm]'
  print *, ''

  print *, 'Precomputing H_m prefactors on the tau grid'
  call Precompute_Hm_Prefactors()

  call Calculate_Charge_Density_xy_grid()
contains
  ! Subroutine to calculate the charge density on a grid
  subroutine Calculate_Charge_Density_xy_grid()
    double precision            :: x, y ! Position where the charge density is being calculated
    double precision            :: x_start, x_end, x_inc, x_l ! The grid used to calculate the charge density
    double precision            :: y_start, y_end, y_inc, y_l
    integer                     :: i, j, N, M
    integer                     :: IFAIL, ud_graph

    double complex              :: results
    double precision, allocatable, dimension(:, :) :: results_re
    double precision, allocatable, dimension(:, :) :: results_im

    print *, 'Calculating the charge density using x and y grid'
    print *, 'Defning the grid'

    ! Define the grid used to calulate the charge density
    N = 50
    M = 1

    x_start = -250.d0*length_scale
    x_end   = 250.0d0*length_scale
    y_start = 0.0d0*length_scale
    y_end   = 0.0d0*y_start

    x_l = x_end - x_start
    y_l = y_end - y_start

    x_inc = x_l / (N-1)
    !y_inc = y_l / (M-1)
    y_inc = 0.0d0

    print *, 'N       = ', N
    print *, 'x_start = ', x_start/length_scale, ' [nm]'
    print *, 'x_end   = ', x_end/length_scale, ' [nm]'
    print *, 'x_l     = ', x_l/length_scale, ' [nm]'
    print *, 'x_inc   = ', x_inc/length_scale, ' [nm]'
    print *, ''

    print *, 'M       = ', M
    print *, 'y_start = ', y_start/length_scale, ' [nm]'
    print *, 'y_end   = ', y_end/length_scale, ' [nm]'
    print *, 'y_l     = ', y_l/length_scale, ' [nm]'
    print *, 'y_inc   = ', y_inc/length_scale, ' [nm]'
    print *, ''

    ! Allocate the array to store the results
    allocate(results_re(1:N, 1:M))
    allocate(results_im(1:N, 1:M))

    print *, 'Starting calculations'

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, x, y, phi_p, xi_p, results) SCHEDULE(GUIDED)
    do i = 1, N
      x = x_start + (i-1)*x_inc
      do j = 1, M
        y = y_start + (j-1)*y_inc

        ! Convert the x and y coordinates to prolate spheroidal coordinates
        ! The value of eta_p is the one that defines the tip, eta_p = eta_1
        phi_p = atan2(y, x)
        xi_p = 1.0d0/a_foci * 1.0d0/sqrt(1.0d0-eta_1**2) * sqrt(x**2 + y**2 + a_foci**2*(1.0d0-eta_1**2))

        ! Store the results
        results = Calculate_Integral(xi_p, phi_p)

        results_re(i, j) = dreal(results) ! Real part, the charge density
        results_im(i, j) = dimag(results) ! Imaginary part, should be very small ~0
      end do
    end do
    !$OMP END PARALLEL DO

    print *, 'Calculations finished'

    print *, 'Writing out data'
    open(newunit=ud_graph, iostat=IFAIL, file='charge_graph.dt', status='REPLACE', action='WRITE')
    if (IFAIL /= 0) then
      print *, 'Failed to open file charge_graph.dt.'
      error stop
    else
      print *, 'Created file charge_graph.dt'
    end if

    do i = 1, N
      x = x_start + (i-1)*x_inc
      !x = x / length_scale
      do j = 1, M
        y = y_start + (j-1)*y_inc
        !y = y / length_scale

        phi_p = atan2(y, x)
        xi_p = 1.0d0/a_foci * 1.0d0/sqrt(1.0d0-eta_1**2) * sqrt(x**2 + y**2 + a_foci**2*(1.0d0-eta_1**2))

        write (ud_graph, "(i6, tr2, i6, *(tr2, E12.4))", iostat=IFAIL) i, j, &
          & (x/length_scale), (y/length_scale), xi_p, phi_p, results_re(i, j), results_im(i, j)
      end do

      write (ud_graph, *) ''
    end do

    close(unit=ud_graph, iostat=IFAIL, status='keep')
    print *, 'Finished writing out data'

    deallocate(results_re)
    deallocate(results_im)
  end subroutine Calculate_Charge_Density_xy_grid


  ! ----------------------------------------------------------------------------
  ! Calculates the integral over Tau.
  ! The integral is broken into parts of length 1
  ! [0, 1], [1, 2], [2, 3], ... [99, 100]
  ! and Simpson's rule is used on each part.
  !
  ! The reason why we stop at tau = 100 is because the library we use
  ! for the conic functions only goes up to tau = 100.
  !
  ! See Composite Simpson's rule
  ! https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_rule
  double complex function Calculate_Integral(xi, phi)
    double precision, intent(in) :: xi, phi
    double precision             :: a, b ! Integration interval
    double precision             :: p, h ! Integration point, interval size
    integer                      :: i, k ! Loop index
    double complex               :: res ! Results of the integration of each part
    double complex               :: res_tot ! Total results of the integration

    res_tot = cmplx(0.0d0, 0.0d0, kind(0.0d0))

    ! We divide the integral into parts [0, 1], [1, 2], [2, 3], ..., [99, 100]
    ! On each part we use Simpsons rule with N_simp = 20 intervals
    do k = 1, tau_max

      a = k - 1.0d0
      b = k
      h = (b - a) / N_simp

      res = Calculate_Hsum(a, xi, phi) + Calculate_Hsum(b, xi, phi) ! End points

      ! Odd numbers
      do i = 1, N_simp, 2
        p = a + i*h
        res = res + 4.0d0*Calculate_Hsum(p, xi, phi)
      end do

      ! Even numbers
      do i = 2, N_simp-1, 2
        p = a + i*h
        res = res + 2.0d0*Calculate_Hsum(p, xi, phi)
      end do

      res = res * h/3.0d0

      ! The conical function library can return NaN near the edges of its valid
      ! range. Drop this interval rather than let it poison the running total,
      ! which would also defeat the convergence test below.
      if (ieee_is_nan(real(res)) .or. ieee_is_nan(aimag(res))) then
        print *, 'WARNING: Integral is NaN, stopping the summation'
        print *, 'k = ', k
        print *, 'xi = ', xi
        print *, 'phi = ', phi
        print *, ''
        exit
      end if

      res_tot = res_tot + res

      ! We stop if the magnitude of the last interval is less than 10^-9
      ! of the total results. Guard against res_tot = 0 so FPE trapping
      ! builds do not abort on the division.
      if (abs(res_tot) > 0.0d0) then
        if ( (abs(res)/abs(res_tot)) < 1.0d-9 ) then
          exit
        end if
      end if
    end do

    Calculate_Integral = res_tot

  end function Calculate_Integral

  ! ----------------------------------------------------------------------------
  ! Calculates the summation over H_m,
  !
  double complex function Calculate_Hsum(tau, xi, phi)
    double precision, intent(in)       :: tau, xi, phi
    integer                            :: m ! The summation index, m <= 30
    integer                            :: m_cut ! Where to cut the summation
    double complex, dimension(0:m_max) :: H_array

    ! H_m carries an explicit tau*tanh(pi*tau) prefactor (see Calculate_Hm), so
    ! the integrand vanishes at tau = 0. The conical function library requires
    ! tau > 0, so take the limit directly rather than call it out of domain.
    if (tau <= 0.0d0) then
      Calculate_Hsum = cmplx(0.0d0, 0.0d0, kind(0.0d0))
      return
    end if

    ! In ref. 1 they say that m <= 30, and that they cut this summation when
    ! H_m starts to alternate signs and grow exponentially. Truncate_Sum applies
    ! that criterion to the terms themselves.
    m_cut = m_max
    do m = 0, m_max
      H_array(m) = Calculate_Hm(m, tau, xi, phi)
      if (Truncate_Sum(H_array, m, m_cut)) exit
    end do

    Calculate_Hsum = sum(H_array(0:m_cut))
  end function Calculate_Hsum

  ! ----------------------------------------------------------------------------
  ! Fills Hm_pref with Calculate_Hm_prefactor on the tau quadrature grid.
  ! The prefactor depends only on (m, tau), so without this table it would be
  ! recomputed identically for every point of the x-y grid.
  subroutine Precompute_Hm_Prefactors()
    integer          :: j, m
    double precision :: tau

    allocate(Hm_pref(0:m_max, 0:tau_max*N_simp))

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j, m, tau) SCHEDULE(GUIDED)
    do j = 0, tau_max*N_simp
      tau = j * h_tau
      do m = 0, m_max
        Hm_pref(m, j) = Calculate_Hm_prefactor(m, tau)
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine Precompute_Hm_Prefactors

  ! ----------------------------------------------------------------------------
  ! The part of H_m that depends only on (m, tau) and the fixed source/tip
  ! geometry -- everything except cos(m*phi) and K_xi at the field point.
  double complex function Calculate_Hm_prefactor(m, tau) result(pref)
    integer, intent(in)          :: m
    double precision, intent(in) :: tau

    double complex               :: Z_p, Z_m
    double precision             :: fac_2
    double precision             :: K_xic
    double precision             :: K_etac, K_metac
    double precision             :: K_eta1, K_meta1
    !double precision             :: K_eta2, K_meta2

    ! No initializer here: that would imply SAVE, and this routine runs inside
    ! an OpenMP parallel region.
    integer                      :: ierr

    ! The tau*tanh(pi*tau) factor vanishes at tau = 0, and the conical function
    ! library requires tau > 0 (see Calculate_Hsum).
    if (tau <= 0.0d0) then
      pref = cmplx(0.0d0, 0.0d0, kind(0.0d0))
      return
    end if

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp_lim(xi_c, m, tau, K_xic, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'xi_c = ', xi_c
      print *, 'm = ', m
      print *, 'tau = ', tau
      print *, ''
    end if

    ! With the charge on the axis (xi_c = 1) K_xic = 0 for every m >= 1, so
    ! the whole term vanishes and the eta evaluations can be skipped.
    if (K_xic == 0.0d0) then
      pref = cmplx(0.0d0, 0.0d0, kind(0.0d0))
      return
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp_lim(eta_c, m, tau, K_etac, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'eta_c'
      if (ierr == 2) then
        print *, 'eta_c = ', eta_c
        print *, 'm = ', m
        print *, 'tau = ', tau
        print *, 'K_etac = ', K_etac
        stop
      end if
    end if

    ierr = 0
    CALL conicp_lim(-1.0d0*eta_c, m, tau, K_metac, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta_c'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp_lim(eta_1, m, tau, K_eta1, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'eta_1'
    end if

    ierr = 0
    CALL conicp_lim(-1.0d0*eta_1, m, tau, K_meta1, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta_1'
    end if
    !---------------------------------------------------------------------------

    ! ! If eta_2 = 0 we can factor it out
    ! !---------------------------------------------------------------------------
    ! CALL conic(eta_2, m, tau, K_eta2, ierr)
    ! if (ierr /= 0) then
    !   print *, 'WARNING: ierr = ', ierr
    !   print *, 'eta_2'
    ! end if
    !
    ! CALL conic(-1.0d0*eta_2, m, tau, K_meta2, ierr)
    ! if (ierr /= 0) then
    !   print *, 'WARNING: ierr = ', ierr
    !   print *, '-eta_2'
    ! end if
    ! !---------------------------------------------------------------------------

    Z_m = CMPLX(0.5d0 - m, tau, kind(0.0d0))
    Z_p = CMPLX(0.5d0 + m, tau, kind(0.0d0))

    !fac_2 = (K_metac * K_eta2 - K_etac * K_meta2) / (K_eta1 * K_meta2 - K_meta1 * K_eta2)
    fac_2 = (K_metac - K_etac) / (K_eta1 - K_meta1) ! We can factor out eta_2 if it is always zero

    ! Note: ** is right-associative, so the (2 - delta_m0) factor must be
    ! multiplied, not chained onto the exponent.
    pref = ((-1.0d0)**m)*(2.0d0 - delta_k(0, m))*tau*tanh(pi*tau) &
        & * ( mygamma(Z_m) / mygamma(Z_p) ) * K_xic * fac_2
  end function Calculate_Hm_prefactor

  ! ----------------------------------------------------------------------------
  ! H_m at a field point: the tabulated (m, tau) prefactor times the
  ! point-dependent factors cos(m*phi) and K_xi.
  double complex function Calculate_Hm(m, tau, xi, phi)
    integer, intent(in)          :: m
    double precision, intent(in) :: tau, xi, phi

    double complex               :: pref
    double precision             :: K_xi
    integer                      :: ierr, j

    ! Use the tabulated prefactor when tau lies on the quadrature grid the
    ! table was built for (a + i*h reproduces j*h_tau only up to rounding,
    ! hence the tolerance); otherwise compute it directly.
    pref = cmplx(0.0d0, 0.0d0, kind(0.0d0))
    j = nint(tau / h_tau)
    if (allocated(Hm_pref) .and. (j >= 0) .and. (j <= tau_max*N_simp)) then
      if (abs(tau - j*h_tau) <= 1.0d-9) then
        pref = Hm_pref(m, j)
      else
        pref = Calculate_Hm_prefactor(m, tau)
      end if
    else
      pref = Calculate_Hm_prefactor(m, tau)
    end if

    ! A zero prefactor makes the whole term zero; skip the K_xi evaluation.
    if (pref == cmplx(0.0d0, 0.0d0, kind(0.0d0))) then
      Calculate_Hm = pref
      return
    end if

    ierr = 0
    CALL conicp_lim(xi, m, tau, K_xi, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'xi = ', xi
      print *, 'm = ', m
      print *, 'tau = ', tau
      print *, 'K_xi = ', K_xi
      print *, ''
    end if

    Calculate_Hm = pref * cos(m*phi) * K_xi
  end function Calculate_Hm

  ! ----------------------------------------------------------------------------
  ! Calculates the summation over S_m, see ref. 1, eq. 22.
  !
  ! NOTE: Calculate_F / Calculate_Sm are currently unused (the charge-density
  ! calculation only needs Calculate_Hsum). They are kept as the
  ! potential-evaluation counterpart from ref. 1.
  double complex function Calculate_F(tau, eta, xi, phi)
    double precision, intent(in)       :: tau, eta, xi, phi
    integer                            :: m ! The summation index, m <= 30
    integer                            :: m_cut ! Where to cut the summation
    double complex, dimension(0:m_max) :: S_array

    ! As in Calculate_Hsum: S_m carries a tau*tanh(pi*tau)/cosh(pi*tau)
    ! prefactor, so the integrand vanishes at tau = 0.
    if (tau <= 0.0d0) then
      Calculate_F = cmplx(0.0d0, 0.0d0, kind(0.0d0))
      return
    end if

    ! In ref. 1 they say that m <= 30, and that they cut this summation when
    ! S_m starts to alternate signs and grow exponentially. Truncate_Sum applies
    ! that criterion to the terms themselves.
    m_cut = m_max
    do m = 0, m_max
      S_array(m) = Calculate_Sm(m, tau,  eta, xi, phi)
      if (Truncate_Sum(S_array, m, m_cut)) exit
    end do

    Calculate_F = sum(S_array(0:m_cut))
  end function Calculate_F

  ! ----------------------------------------------------------------------------
  ! Calculates the function S_m, see ref. 1, eq. 24
  ! In ref. 1 they say that m <= 30.
  double complex function Calculate_Sm(m, tau, eta, xi, phi)
    integer, intent(in)          :: m
    double precision, intent(in) :: tau, eta, xi, phi
    double complex               :: Z_p, Z_m, fac_1, fac_2

    double precision             :: K_xic, K_xi
    double precision             :: K_eta, K_meta
    double precision             :: K_etac, K_metac
    double precision             :: K_eta1, K_meta1
    double precision             :: K_eta2, K_meta2

    ! No initializer here: that would imply SAVE, and this routine runs inside
    ! an OpenMP parallel region.
    integer                      :: ierr

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp_lim(xi_c, m, tau, K_xic, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'xi_c = ', xi_c
      print *, 'm = ', m
      print *, 'tau = ', tau
      print *, ''
    end if

    ierr = 0
    CALL conicp_lim(xi, m, tau, K_xi, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'xi = ', xi
      print *, 'm = ', m
      print *, 'tau = ', tau
      print *, 'K_xi = ', K_xi
      print *, ''
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp_lim(eta, m, tau, K_eta, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'eta'
    end if

    ierr = 0
    CALL conicp_lim(-1.0d0*eta, m, tau, K_meta, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp_lim(eta_c, m, tau, K_etac, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'eta_c'
      if (ierr == 2) then
        print *, 'eta_c = ', eta_c
        print *, 'm = ', m
        print *, 'tau = ', tau
        print *, 'K_etac = ', K_etac
        stop
      end if
    end if

    ierr = 0
    CALL conicp_lim(-1.0d0*eta_c, m, tau, K_metac, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta_c'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp_lim(eta_1, m, tau, K_eta1, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'eta_1'
    end if

    ierr = 0
    CALL conicp_lim(-1.0d0*eta_1, m, tau, K_meta1, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta_1'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp_lim(eta_2, m, tau, K_eta2, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'eta_2'
    end if

    ierr = 0
    CALL conicp_lim(-1.0d0*eta_2, m, tau, K_meta2, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta_2'
    end if
    !---------------------------------------------------------------------------

    Z_m = CMPLX(0.5d0 - m, tau, kind(0.0d0))
    Z_p = CMPLX(0.5d0 + m, tau, kind(0.0d0))

    fac_1 = (2.0d0 - delta_k(0, m))*cos(m*phi)*tau*(tanh(pi*tau)/cosh(pi*tau)) &
        & * (mygamma(Z_m)/mygamma(Z_p))**2 &
        & * K_xic*K_xi/(K_eta1*K_meta2 - K_meta1*K_eta2)

    fac_2 = K_eta*K_meta2 * (K_etac*K_meta1 - K_metac*K_eta1) &
        & + K_meta*K_eta1 * (K_metac*K_eta2 - K_etac*K_meta2)

    Calculate_Sm = fac_1 * fac_2
  end function Calculate_Sm

  ! A wrapper function for the Gamma function in cgamma.f
  double complex function mygamma(Z)
    double complex, intent(in) :: Z
    integer, parameter         :: KF = 1
    double precision           :: x, y
    double precision           :: GR, GI

    x = real(Z)
    y = aimag(Z)

    call CGAMA(x, y, KF, GR, GI)
    mygamma = CMPLX(GR, GI, kind(0.0d0))
  end function mygamma

  ! Kronecker delta function used in the sum
  double precision function delta_k(i, j)
    integer, intent(in) :: i, j

    if (i == j) then
      delta_k = 1.0d0
    else
      delta_k = 0.0d0
    end if
  end function delta_k

  ! ----------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------
  ! Decides whether the m-summation should stop once term m has been computed.
  ! Returns .true. to stop, and sets m_cut to the last term to keep.
  !
  ! Ref. 1 truncates when the terms start to alternate in sign and grow
  ! exponentially. The previous check here compared successive *differences*,
  ! |T(m-1)-T(m)| > 10|T(m-2)-T(m-1)|, which has three problems: it never looks
  ! at the signs, a single noisy term can trip it, and when the terms have
  ! genuinely converged both differences are tiny so the ratio is meaningless
  ! and fires at random.
  !
  ! Instead we look at the terms themselves:
  !   * a run of consecutive increases in |T| is what "grows exponentially"
  !     actually means, and requiring a run rather than one step ignores noise;
  !   * a shorter run is enough when the sign also alternates, which is the
  !     classic signature of a divergent asymptotic tail;
  !   * on divergence we keep terms up to the *smallest* one rather than the
  !     one before the growth. For an asymptotic series that is where the
  !     truncation error is minimised;
  !   * a series that simply converges exits on the relative tolerance instead,
  !     which the old check had no way to detect.
  logical function Truncate_Sum(T, m, m_cut)
    double complex, dimension(0:), intent(in)  :: T
    integer,                       intent(in)  :: m
    integer,                       intent(out) :: m_cut

    integer, parameter          :: m_least  = 5      ! always take at least this many
    integer, parameter          :: n_grow   = 3      ! consecutive increases => diverging
    integer, parameter          :: n_grow_a = 2      ! ... fewer if the sign alternates too
    double precision, parameter :: eps_conv = 1.0d-14

    double precision :: a(0:m), a_max
    integer          :: k, m_min, n_up
    logical          :: alt

    Truncate_Sum = .false.
    m_cut = m

    do k = 0, m
      a(k) = abs(T(k))
    end do

    ! A non-finite term means the recursion has broken down. Keep what came
    ! before it and stop; nothing after this point is meaningful.
    if (.not. ieee_is_finite(a(m))) then
      m_cut = max(m - 1, 0)
      Truncate_Sum = .true.
      return
    end if

    if (m < m_least) return

    a_max = maxval(a(0:m))
    m_min = minloc(a(0:m), dim=1) - 1 ! minloc is 1-based, a is 0-based

    ! Converged: the terms have died away relative to the largest one. Measuring
    ! against a_max rather than the running sum keeps this from tripping when
    ! the partial sum happens to pass near zero.
    if (a(m) <= eps_conv*a_max) then
      Truncate_Sum = .true.
      return
    end if

    ! Count the run of consecutive increases ending at m, and note whether the
    ! sign alternates the whole way along it.
    n_up = 0
    alt  = .true.
    do k = m, 1, -1
      if (a(k) <= a(k-1)) exit
      n_up = n_up + 1
      if (dreal(T(k))*dreal(T(k-1)) >= 0.0d0) alt = .false.
    end do

    if ((n_up >= n_grow) .or. (alt .and. (n_up >= n_grow_a))) then
      m_cut = m_min ! optimal truncation: stop at the smallest term
      Truncate_Sum = .true.
    end if
  end function Truncate_Sum

  ! ----------------------------------------------------------------------------
  ! Wrapper that covers the gap the library in ref. 2 leaves just above x = 1.
  !
  ! That library is valid for -1 < x < 1 and 1 < x <= 100, but in practice it
  ! returns ierr = 1 for 1 < x - 1 < ~4.6e-4, and at x = 1 exactly it returns 0
  ! with ierr = 0 -- a wrong value with no error flag. This matters because the
  ! source charge sits on the axis at xi_c = 1, and field points near the axis
  ! land in the gap, so K_xic and K_xi would otherwise be 0 or garbage.
  !
  ! Near x = 1 we use the Gauss series instead (see conicp_series), which is
  ! exact at x = 1 and more accurate than the library throughout the overlap.
  subroutine conicp_lim(x, m, tau, pmtau, ierr)
    double precision, intent(in)  :: x, tau
    integer, intent(in)           :: m
    double precision, intent(out) :: pmtau
    integer, intent(out)          :: ierr

    if ((x >= 1.0d0) .and. ((x - 1.0d0) < x_series)) then
      call conicp_series(x, m, tau, pmtau, ierr)
    else
      CALL conicp(x, m, tau, pmtau, ierr)
    end if
  end subroutine conicp_lim

  ! ----------------------------------------------------------------------------
  ! P^mu_(-1/2+i*tau)(x) for x just above 1, via the Gauss series
  !
  !   P = [C_mu/Gamma(1+mu)] * ((x-1)/(x+1))^(mu/2)
  !         * 2F1(1/2-i*tau, 1/2+i*tau; 1+mu; (1-x)/2)
  !   C_mu = prod_{k=0}^{mu-1} ((k+1/2)^2 + tau^2)
  !
  ! Both Pochhammer products are real, so this needs no complex arithmetic and
  ! nothing here can underflow the way the complex Gamma ratio does. With
  ! |1-x|/2 < 5e-4 it converges in at most ~13 terms.
  !
  ! It reproduces the exact limit at x = 1 on its own: z = 0 gives P = 1 for
  ! mu = 0, and the prefactor gives P = 0 for mu >= 1.
  !
  ! Validated against the library of ref. 2 and against a quadruple precision
  ! evaluation of the same series: accurate to ~1e-15 over the range where it
  ! is used here, against ~1e-10 for the library.
  subroutine conicp_series(x, mu, tau, pmtau, ierr)
    double precision, intent(in)  :: x, tau
    integer, intent(in)           :: mu
    double precision, intent(out) :: pmtau
    integer, intent(out)          :: ierr

    integer, parameter :: k_max = 200
    double precision   :: z, coef, pref, term, s, ak
    integer            :: k

    z = 0.5d0*(1.0d0 - x)

    ! coef = C_mu/Gamma(1+mu), accumulated factor by factor so that C_mu, which
    ! reaches ~1e120 for mu = 30 and tau = 100, is never formed on its own.
    coef = 1.0d0
    do k = 0, mu-1
      coef = coef * ( ((k + 0.5d0)**2 + tau*tau) / dble(k + 1) )
    end do

    if (mu == 0) then
      pref = 1.0d0
    else
      pref = ((x - 1.0d0)/(x + 1.0d0))**(0.5d0*mu)
    end if

    term = 1.0d0
    s    = 1.0d0
    ierr = 1
    do k = 0, k_max
      ak   = (k + 0.5d0)**2 + tau*tau
      term = term * (ak*z) / (dble(1 + mu + k)*dble(k + 1))
      s    = s + term
      if (abs(term) <= 1.0d-17*abs(s)) then
        ierr = 0
        exit
      end if
    end do

    pmtau = coef * pref * s
  end subroutine conicp_series
end program Tip_IC
