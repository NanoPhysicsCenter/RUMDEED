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
  implicit none

  double precision, parameter :: pi = 3.14159265358979324d0 ! Pi
  double precision, parameter :: epsilon_0 = 8.854187817d-12 ! Farads / meters
  double precision, parameter :: m_e = 9.10938291d-31 ! kg
  double precision, parameter :: q_e = -1.602176565d-19 ! Coulomb

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
  print *, 'x     = ', (a_foci*sqrt(xi_c**2-1.0d0)*sqrt(1.0d0-eta_c**2)*cos(phi_p)/length_scale), ' [nm]'
  print *, 'y     = ', (a_foci*sqrt(xi_c**2-1.0d0)*sqrt(1.0d0-eta_c**2)*sin(phi_p)/length_scale), ' [nm]'
  print *, 'z     = ', (a_foci*xi_c*eta_c/length_scale), ' [nm]'
  print *, ''

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
    integer, parameter           :: N = 20 ! Number of intervals
    double complex               :: res ! Results of the integration of each part
    double complex               :: res_tot ! Total results of the integration

    res_tot = cmplx(0.0d0, 0.0d0)

    ! We divide the integral into parts [0, 1], [1, 2], [2, 3], ..., [99, 100]
    ! On each part we use Simpsons rule with N = 20 intervals
    do k = 1, 100

      a = k - 1.0d0
      b = k
      h = (b - a) / N

      res = Calculate_Hsum(a, xi, phi) + Calculate_Hsum(b, xi, phi) ! End points

      ! Odd numbers
      do i = 1, N, 2
        p = a + i*h
        res = res + 4.0d0*Calculate_Hsum(p, xi, phi)
      end do

      ! Even numbers
      do i = 2, N-1, 2
        p = a + i*h
        res = res + 2.0d0*Calculate_Hsum(p, xi, phi)
      end do

      res = res * h/3.0d0
      res_tot = res_tot + res

      ! We stop if the magnitude of the last interval is less than 10^-9
      ! of the total results
      if ( (abs(res)/abs(res_tot)) < 1.0d-9 ) then
        exit
      end if
    end do

    Calculate_Integral = res_tot

  end function Calculate_Integral

  ! ----------------------------------------------------------------------------
  ! Calculates the summation over H_m,
  !
  double complex function Calculate_Hsum(tau, xi, phi)
    double precision, intent(in)       :: tau, xi, phi
    integer, parameter                 :: m_max = 30
    integer                            :: m ! The summation index, m <= 30
    integer                            :: m_cut ! Where to cut the summation
    double complex, dimension(1:m_max) :: H_array
    double precision                   :: diff, last_diff

    ! In ref. 1 they say that m <= 30. And that they cut this summation when
    ! S_m starts to alternate signs and grow exponentially.
    ! We do the same here.
    m_cut = m_max
    do m = 1, m_max

      ! Calculate H_m(tau)
      H_array(m) = Calculate_Hm(m, tau, xi, phi)

      if (m > 5) then ! Always do at least 5
        diff = abs(H_array(m-1) - H_array(m))
        last_diff = abs(H_array(m-2) - H_array(m-1))

        ! Check for exponental growth
        if (diff > 10.0d0*last_diff) then
          m_cut = m - 1
          exit ! Ok cut the summation and exit the loop
        end if

      end if
    end do

    Calculate_Hsum = sum(H_array(1:m_cut))
  end function Calculate_Hsum

  double complex function Calculate_Hm(m, tau, xi, phi)
    integer, intent(in)          :: m
    double precision, intent(in) :: tau, xi, phi

    double complex               :: Z_p, Z_m, fac_1

    double precision             :: fac_2
    double precision             :: K_xic, K_xi
    double precision             :: K_etac, K_metac
    double precision             :: K_eta1, K_meta1
    !double precision             :: K_eta2, K_meta2

    integer                      :: ierr = 0

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp(xi_c, m, tau, K_xic, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'xi_c = ', xi_c
      print *, 'm = ', m
      print *, 'tau = ', tau
      print *, ''
    end if

    ierr = 0
    CALL conicp(xi, m, tau, K_xi, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'xi'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp(eta_c, m, tau, K_etac, ierr)
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
    CALL conicp(-1.0d0*eta_c, m, tau, K_metac, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta_c'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp(eta_1, m, tau, K_eta1, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'eta_1'
    end if

    ierr = 0
    CALL conicp(-1.0d0*eta_1, m, tau, K_meta1, ierr)
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

    Z_m = CMPLX(0.5d0 - m, tau)
    Z_p = CMPLX(0.5d0 + m, tau)

    fac_1 = cos(m*phi)*tau*tanh(pi*tau) * ( mygamma(Z_m) / mygamma(Z_p) ) * (K_xic * K_xi)
    !fac_2 = (K_metac * K_eta2 - K_etac * K_meta2) / (K_eta1 * K_meta2 - K_meta1 * K_eta2)
    fac_2 = (K_metac - K_etac) / (K_eta1 - K_meta1) ! We can factor out eta_2 if it is always zero

    Calculate_Hm = (-1.0d0)**m**(2.0d0-delta_k(0, m))*fac_1 * fac_2
  end function Calculate_Hm

  ! ----------------------------------------------------------------------------
  ! Calculates the summation over S_m, see ref. 1, eq. 22.
  !
  double complex function Calculate_F(tau, eta, xi, phi)
    double precision, intent(in)       :: tau, eta, xi, phi
    integer, parameter                 :: m_max = 30
    integer                            :: m ! The summation index, m <= 30
    integer                            :: m_cut ! Where to cut the summation
    double complex, dimension(1:m_max) :: S_array
    double precision                   :: diff, last_diff

    ! In ref. 1 they say that m <= 30. And that they cut this summation when
    ! S_m starts to alternate signs and grow exponentially.
    do m = 1, m_max

      ! Calculate S_m(tau)
      S_array(m) = Calculate_Sm(m, tau,  eta, xi, phi)

      if (m > 5) then ! Always do at least 5
        diff = abs(S_array(m-1) - S_array(m))
        last_diff = abs(S_array(m-2) - S_array(m-1))

        ! Check for exponental growth
        if (diff > 10.0d0*last_diff) then
          m_cut = m - 1
          exit ! Ok cut the summation and exit the loop
        end if

      end if
    end do

    Calculate_F = sum(S_array(1:m_cut))
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

    integer                      :: ierr = 0

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp(xi_c, m, tau, K_xic, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'xi_c = ', xi_c
      print *, 'm = ', m
      print *, 'tau = ', tau
      print *, ''
    end if

    ierr = 0
    CALL conicp(xi, m, tau, K_xi, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'xi'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp(eta, m, tau, K_eta, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'eta'
    end if

    ierr = 0
    CALL conicp(-1.0d0*eta, m, tau, K_meta, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp(eta_c, m, tau, K_etac, ierr)
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
    CALL conicp(-1.0d0*eta_c, m, tau, K_metac, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta_c'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp(eta_1, m, tau, K_eta1, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'eta_1'
    end if

    ierr = 0
    CALL conicp(-1.0d0*eta_1, m, tau, K_meta1, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta_1'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conicp(eta_2, m, tau, K_eta2, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'eta_2'
    end if

    ierr = 0
    CALL conicp(-1.0d0*eta_2, m, tau, K_meta2, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta_2'
    end if
    !---------------------------------------------------------------------------

    Z_m = CMPLX(0.5d0 - m, tau)
    Z_p = CMPLX(0.5d0 + m, tau)

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
    mygamma = CMPLX(GR, GI)
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
end program Tip_IC
