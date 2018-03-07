!-------------------------------------------------------------------------------
! Kristinn Torfason
! 04.03.15
! Code to calculate image charge for the prolate spheroidal tip
!
! Based on:
! [1] Methods for calculating electrostatic quantities due to a free charge
! in a nanoscale three‐dimensional tip/base junction
! Published in:
! Journal of applied physics
! 78, no. 8 (1995): 4888-4894.
! http://dx.doi.org/10.1063/1.359777
!
! Code for the conical functions:
! [2] An improved algorithm and a Fortran 90 module for computing the
! conical function P^m_{-1/2+i\tau} (\chi)
! Published in:
! Computer Physics Communications
! Volume 183, Issue 3, March 2012, Pages 794–799
! http://dx.doi.org/10.1016/j.cpc.2011.11.025
!
! Code for the Gamma function:
! http://jin.ece.illinois.edu/routines/routines.html
! [3] Computation of special functions
! S. Zhang and J. Jin
! Wiley, 1996
! ISBN 0-471-11963-6
!-------------------------------------------------------------------------------
program IC
  USE Conical ! The module for the conical functions
  USE someconstants ! Constants used in the conical module

  implicit none
  intrinsic cmplx, real, aimag, tanh, atan2, sin, cos, tan, log, sqrt, abs, nint

  !double precision, parameter :: pi = 3.14159265358979324d0 ! Pi
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
  !call Write_Sm()

  !CALL conic(-eta_1, 0, 0.0d0, test, ierr)

  !print *, 'test = ', test
  !print *, 'ierr = ', ierr
contains
  subroutine Calculate_Charge_Density_xiphi_grid()
    double precision            :: x, y
    double precision            :: xi, phi
    double precision            :: xi_start, xi_end, xi_inc, xi_l ! The grid used to calculate the charge density
    double precision            :: phi_start, phi_end, phi_inc, phi_l
    integer                     :: i, j, N, M
    integer                     :: IFAIL, ud_graph

    double complex              :: results
    double precision, allocatable, dimension(:, :) :: results_re
    double precision, allocatable, dimension(:, :) :: results_im

    print *, 'Calculating the charge density using xi and phi grid'
    print *, 'Defning the grid'

    N = 10
    M = 10

    xi_start = 1.0d0
    xi_end = 3.0d0
    phi_start = pi
    phi_end = 2.0d0*pi

    xi_l = xi_end - xi_start
    phi_l = phi_end - phi_start

    xi_inc = xi_l / (N-1)
    phi_inc = phi_l / (M-1)

    print *, 'N = ', N
    print *, 'xi_start = ', xi_start
    print *, 'xi_end  = ', xi_end
    print *, 'xi_l  = ', xi_l
    print *, 'xi_inc  = ', xi_inc
    print *, ''

    print *, 'M = ', M
    print *, 'phi_start = ', phi_start
    print *, 'phi_end  = ', phi_end
    print *, 'phi_l  = ', phi_l
    print *, 'phi_inc  = ', phi_inc
    print *, ''

    ! Allocate the array to store the results
    allocate(results_re(1:N, 1:M))
    allocate(results_im(1:N, 1:M))

    print *, 'Starting calculations'
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, phi, xi, results) SCHEDULE(GUIDED)
    do i = 1, N
      xi = xi_start + (i-1)*xi_inc
      do j = 1, M
        phi = phi_start + (j-1)*phi_inc

        ! Store the results
        results = SumM(xi, phi)
        results_re(i, j) = real(results) ! Real part, the charge density
        results_im(i, j) = aimag(results) ! Imaginary part, should be very small ~0
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
      xi = xi_start + (i-1)*xi_inc
      x = a_foci*sqrt(xi**2-1.0d0)*sqrt(1.0d0-eta_1**2)*cos(phi) / length_scale

      do j = 1, M
        phi = phi_start + (j-1)*phi_inc
        y = a_foci*sqrt(xi**2-1.0d0)*sqrt(1.0d0-eta_1**2)*sin(phi) / length_scale

        write (ud_graph, "(i6, tr2, i6, *(tr2, E12.4))", iostat=IFAIL) i, j, &
          & x, y, xi, phi, results_re(i, j), results_im(i, j)
      end do

      write (ud_graph, *) ''
    end do

    close(unit=ud_graph, iostat=IFAIL, status='keep')
    print *, 'Finished writing out data'

    deallocate(results_re)
    deallocate(results_im)
  end subroutine Calculate_Charge_Density_xiphi_grid

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
        results = TauInt(xi_p, phi_p)
        if (isnan(real(results)) .eqv. .true.) then
          print *, 'Warning: TauInt returned NaN trying SumM'
          results = SumM(xi_p, phi_p)
          if (isnan(real(results)) .eqv. .true.) then
            print *, 'Warning: SumM also returned NaN'
          end if
        end if
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

  subroutine Write_Sm()
    integer                            :: m, ud_sm, IFAIL
    integer, parameter                 :: m_max = 30
    double complex                     :: results
    double complex, dimension(0:m_max) :: results_re, results_im
    double precision                   :: tau, eta, xi, phi

    tau = 1.0d0
    xi = 1.2d0
    eta = 0.75d0
    phi = 0.0d0

    print *, 'Calculating S_m for'
    print *, 'tau = ', tau
    print *, 'xi = ', xi
    print *, 'eta = ', eta
    print *, 'phi = ', phi
    print *, ''

    do m = 0, m_max
      results = Calculate_Sm(m, tau, eta, xi, phi)
      results_re(m) = Real( results )
      results_im(m) = aimag( results )
    end do

    open(newunit=ud_sm, iostat=IFAIL, file='sm_graph.dt', status='REPLACE', action='write')
    if (IFAIL /= 0) then
      print *, 'Failed to open file sm_graph.dt.'
    else
      print *, 'Created file sm_graph.dt'
    end if

    do m = 0, m_max
      write (ud_sm, "(i6, tr2, E12.4, tr2, E12.4)", iostat=IFAIL) m, &
        & results_re(m), results_im(m)
    end do

    close(unit=ud_sm, iostat=IFAIL, status='keep')
  end subroutine Write_Sm


  ! Calculates the function S_m from ref. 1
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
    CALL conic(xi_c, m, tau, K_xic, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'xi_c = ', xi_c
      print *, 'm = ', m
      print *, 'tau = ', tau
      print *, ''
    end if

    ierr = 0
    CALL conic(xi, m, tau, K_xi, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'xi'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conic(eta, m, tau, K_eta, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'eta'
    end if

    ierr = 0
    CALL conic(-1.0d0*eta, m, tau, K_meta, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conic(eta_c, m, tau, K_etac, ierr)
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
    CALL conic(-1.0d0*eta_c, m, tau, K_metac, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta_c'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conic(eta_1, m, tau, K_eta1, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'eta_1'
    end if

    ierr = 0
    CALL conic(-1.0d0*eta_1, m, tau, K_meta1, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta_1'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conic(eta_2, m, tau, K_eta2, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'eta_2'
    end if

    ierr = 0
    CALL conic(-1.0d0*eta_2, m, tau, K_meta2, ierr)
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
        & * K_meta*K_eta1 * (K_metac*K_eta2 - K_etac*K_meta2)

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

  ! ! A wrapper function for the Gamma function in cgamma.f90
  ! ! This one gives strange results. The integral doesn't seem to converge.
  ! double complex function mygamma(Z)
  !   double complex, intent(in) :: Z
  !   double complex             :: W
  !   integer, parameter         ::mo = 0
  !
  !   call cgamma(mo, Z, W)
  !   mygamma = W
  ! end function mygamma


  ! ! A wrapper function for the Gamma function in cdgamma.f
  ! ! This one seems to work also.
  ! double complex function mygamma(Z)
  !   double complex, intent(in) :: Z
  !   double complex, external :: cdgamma
  !
  !   mygamma = cdgamma(z)
  ! end function mygamma


  ! The function that calculates the product of the Gamma and Legender
  ! functions inside the integral.
  double complex function LegendFun(tau, m, xi)
    double precision, intent(in) :: tau, xi
    integer, intent(in)          :: m
    double complex               :: Z_p, Z_m, fac_1

    double precision             :: fac_2
    double precision             :: K_xic, K_xi
    double precision             :: K_etac, K_metac
    double precision             :: K_eta1, K_meta1
    !double precision             :: K_eta2, K_meta2

    integer                      :: ierr = 0

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conic(xi_c, m, tau, K_xic, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'xi_c = ', xi_c
      print *, 'm = ', m
      print *, 'tau = ', tau
      print *, ''
    end if

    ierr = 0
    CALL conic(xi, m, tau, K_xi, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'xi'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conic(eta_c, m, tau, K_etac, ierr)
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
    CALL conic(-1.0d0*eta_c, m, tau, K_metac, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, '-eta_c'
    end if
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ierr = 0
    CALL conic(eta_1, m, tau, K_eta1, ierr)
    if (ierr /= 0) then
      print *, 'WARNING: ierr = ', ierr
      print *, 'eta_1'
    end if

    ierr = 0
    CALL conic(-1.0d0*eta_1, m, tau, K_meta1, ierr)
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

    fac_1 = tau*tanh(pi*tau) * ( mygamma(Z_m) / mygamma(Z_p) ) * (K_xic * K_xi)
    !fac_2 = (K_metac * K_eta2 - K_etac * K_meta2) / (K_eta1 * K_meta2 - K_meta1 * K_eta2)
    fac_2 = (K_metac - K_etac) / (K_eta1 - K_meta1) ! We can factor out eta_2 if it is always zero

    LegendFun = fac_1 * fac_2
  end function LegendFun

  double complex function SumM(xi, phi)
    double precision, intent(in) :: xi, phi
    integer                      :: m, i, k
    double precision, parameter  :: h = 0.05d0 !\Delta\tau = 0.05
    integer, parameter           :: m_max = 30
    double precision             :: tau_0, tau_1, tau_2
    double complex               :: int_res, sum_res
    double complex               :: int_part
    !double complex               :: int_res_last, sum_res_last
    double precision             :: sum_fac, int_fac, E_vac

    ! Vacuum field
    E_vac = -2.0d0*epsilon_0*V_0/a_foci * 1.0d0/(sqrt((xi**2-eta_1**2)*(1-eta_1**2))) &
        & * 1.0d0 / log((1.0d0 + eta_2)/(1.0d0 - eta_2) * (1.0d0 - eta_1)/(1.0d0 + eta_1))

    ! Prefactor before the sum
    sum_fac = q_e/(twopi*a_foci**2) * 1.0d0/( sqrt( (xi**2 - eta_1**2)*(1.0d0 - eta_1**2) ) )

    sum_res = (0.0d0, 0.0d0)
    !sum_res_last = (0.0d0, 0.0d0)
    !!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(m, int_fac, int_res, int_res_last, k, i, tau_0, tau_1, tau_2) REDUCTION(+:sum_res)
    do m = 0, m_max
      !print *, 'm = ', m
      !print *, ''

      ! Prefactor in the sum before the integral
      int_fac = (-1.0d0)**m * (2.0d0 - delta_k(0, m)) * cos(m*phi)

      ! Initialize the varibles used in the tau integration loop
      int_res = (0.0d0, 0.0d0)
      !int_res_last = (0.0d0, 0.0d0)

      ! The loop that calculates the tau integral
      k = 0
      do

        int_part = (0.0d0, 0.0d0)
        !Split the integral into parts [0, 1], [1, 2], ..., etc
        ! Each part is calculated using Simpson's rule with \Delta\tau = 0.05
        do i = 0, 19
          ! The function that calculates the conic function doesn't seem to
          ! like tau = 0, therefor, we add 1.0d-35
          tau_0 = k + i*h + 1.0d-35 !Start point
          tau_2 = k + (i+1)*h + 1.0d-35 !End point
          tau_1 = (tau_0 + tau_2) * 0.5d0 !Mid Point

          int_part = int_part + h*onesix * ( LegendFun(tau_0, m, xi) &
                & + 4.0d0 * LegendFun(tau_1, m, xi) + LegendFun(tau_2, m, xi) )
        end do

        int_res = int_res + int_part
        ! Stop the integration when the magnitude of the last interval
        ! is less than 10^-9 of the magnitude of the total results
        if ((abs(int_part)/abs(int_res)) < 1.0d-9) then
          !print *, 'Interval [', k, ',', k+1, ']'
          !print *, 'Error = ', (abs(int_res - int_res_last)/abs(int_res))
          !print *, 'Value integral = ', int_res
          !print *, 'Value sum = ', (sum_res+int_fac*int_res)
          !print *, 'Stopping'
          !print *, ''
          exit
        else
          !print *, 'Interval [', k, ',', k+1, ']'
          !print *, 'Error = ', (abs(int_res - int_res_last)/abs(int_res))
          !print *, 'Value integral = ', int_res
          !print *, 'Value sum = ', (sum_res+int_fac*int_res)
          !print *, ''

          !int_res_last = int_res
          !int_res = int_res + int_part
          k = k + 1

          if (k > 100) then
            print *, 'Warning tau to large'
            print *, 'k = ', k
            exit
          end if
        end if
      end do

      sum_res = sum_res + int_fac * int_res
      if ((abs(int_fac * int_res)/abs(sum_res)) < 1.0d-9) then
        !print *, 'sum_res error less than 1.0E-9'
        !print *, 'm = ', m
        exit
      !else
        !sum_res_last = sum_res
        !sum_res = sum_res + int_fac * int_res
      end if

      ! if ((abs(sum_res - sum_res_last)/abs(sum_res)) > 1.0d9) then
      !   print *, 'sum_res error larger than 1.0E9'
      !   print *, 'm = ', m
      !   sum_res = sum_res_last
      !   exit
      ! end if

      if (isnan(real(sum_res)) .eqv. .true.) then
        print *, 'Error Real part of sum_res is NaN'
        print *, 'm = ', m
      else if (isnan(aimag(sum_res)) .eqv. .true.) then
        print *, 'Error Imaginary part of sum_res is NaN'
        print *, 'm = ', m
      end if
    end do
    !!!$OMP END PARALLEL DO

    SumM = E_vac + sum_fac*sum_res
  end function SumM

  ! Kronecker delta function used in the sum
  double precision function delta_k(i, j)
    integer, intent(in) :: i, j

    if (i == j) then
      delta_k = 1.0d0
    else
      delta_k = 0.0d0
    end if
  end function delta_k

  ! The function that calculates the integral over tau
  double complex function TauInt(xi, phi)
    double precision, intent(in) :: xi, phi
    integer                      :: i, k
    double complex               :: int_res, int_part
    double precision, parameter  :: h = 0.025d0 !\Delta\tau = 0.05
    integer, parameter           :: n = nint((1.0d0/h-1))
    !double precision             :: tau_0 ,tau_1, tau_2
    double precision             :: a, b
    double precision, parameter  :: gtau_1 = -sqrt(3.0d0/7.0d0 - 2.0d0/7.0d0*sqrt(6.0d0/5.0d0))
    double precision, parameter  :: gtau_2 = -gtau_1
    double precision, parameter  :: gtau_3 = -sqrt(3.0d0/7.0d0 + 2.0d0/7.0d0*sqrt(6.0d0/5.0d0))
    double precision, parameter  :: gtau_4 = -gtau_3
    double precision, parameter  :: w_1 = (18.0d0 + sqrt(30.0d0))/36.0d0
    double precision, parameter  :: w_2 = w_1
    double precision, parameter  :: w_3 = (18.0d0 - sqrt(30.0d0))/36.0d0
    double precision, parameter  :: w_4 = w_3
    double precision             :: E_vac, fac

    ! Initialize the varibles used in the tau integration loop
    int_res = (0.0d0, 0.0d0) !Integral results
    ! Vacuum field
    E_vac = -2.0d0*epsilon_0*V_0/a_foci * 1.0d0/(sqrt((xi**2-eta_1**2)*(1-eta_1**2))) &
        & * 1.0d0 / log((1.0d0 + eta_2)/(1.0d0 - eta_2) * (1.0d0 - eta_1)/(1.0d0 + eta_1))

    ! Prefactor before the sum
    fac = q_e/(twopi*a_foci**2) * 1.0d0/( sqrt( (xi**2 - eta_1**2)*(1.0d0 - eta_1**2) ) )

    ! The loop that calculates the tau integral
    k = 0
    do
      int_part = (0.0d0, 0.0d0)
      ! Split the integral into parts [0, 1], [1, 2], ..., etc
      ! Each part is calculated using Simpson's rule with \Delta\tau = 0.05
      ! do i = 0, 39
      !   tau_0 = k + i*h + 1.0d-35 !Start point
      !   tau_2 = k + (i+1)*h + 1.0d-35 !End point
      !   tau_1 = (tau_0 + tau_2) * 0.5d0 !Mid Point
      !
      !   int_part = int_part + h*onesix * (DoSum(tau_0, xi, phi) &
      !         & + 4.0d0 * DoSum(tau_1, xi, phi) + DoSum(tau_2, xi, phi))
      ! end do

      !Gaussian Quadrature
      do i = 0, n
        a = k + i*h + 1.0d-35
        b = k + (i+1)*h + 1.0d-35

        int_part = int_part + h*0.5d0*( w_1*DoSum(0.5d0*h*gtau_1 + 0.5d0*(a+b), xi, phi) &
                                    & + w_2*DoSum(0.5d0*h*gtau_2 + 0.5d0*(a+b), xi, phi) &
                                    & + w_3*DoSum(0.5d0*h*gtau_3 + 0.5d0*(a+b), xi, phi) &
                                    & + w_4*DoSum(0.5d0*h*gtau_4 + 0.5d0*(a+b), xi, phi) )
      end do

      int_res = int_res + int_part

      ! Stop the integration when the magnitude of the last interval
      ! is less than 10^-9 of the magnitude of the total results
      if ((abs(int_part)/abs(int_res)) < 1.0d-9) then
        exit
      else
        ! print *, 'Interval [', k, ',', k+1, ']'
        ! print *, 'Error = ', (abs(int_res - int_res_prev)/abs(int_res))
        ! print *, 'Value = ', int_res
        ! print *, ''

        k = k + 1

        if (k > 100) then
          print *, 'Warning tau to large'
          exit
        end if
      end if
    end do

    TauInt = E_vac + fac*int_res
  end function TauInt


  ! The function that computes the sum inside the integral
  double complex function DoSum(tau, xi, phi)
    double precision, intent(in) :: tau, xi, phi !Tau value being calculated in the integral
    integer                      :: m !Sum index
    integer, parameter           :: m_max = 40 ! Max m in the conic function is 40
    double complex               :: sum_res, sum_fac, sum_part !Sum results and prefactor

    ! The sum
    sum_res = (0.0d0, 0.0d0)
    do m = 0, m_max
      ! Caclulate the prefactor
      sum_fac = (-1.0d0)**m*(2.0d0 - delta_k(0, m)) * cos(m*phi)

      ! Add the sum results
      sum_part = sum_fac*LegendFun(tau, m, xi)
      sum_res = sum_res + sum_part
      if ((abs(sum_part)/abs(sum_res)) < 1.0d-9) then
        exit
      end if
      if (isnan(real(sum_part)) .eqv. .true.) then
        print *, 'Warning: DoSum returns NaN'
        !stop
        exit
      end if
    end do

    DoSum = sum_res
  end function DoSum
end program IC
