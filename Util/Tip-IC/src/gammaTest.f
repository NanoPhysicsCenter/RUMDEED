subroutine cgama ( x, y, kf, gr, gi )

!*****************************************************************************80
!
!! CGAMA computes the Gamma function for complex argument.
!
!  Discussion:
!
!    This procedcure computes the gamma function �(z) or ln[�(z)]
!    for a complex argument
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    26 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the real and imaginary parts of 
!    the argument Z.
!
!    Input, integer ( kind = 4 ) KF, the function code.
!    0 for ln[�(z)]
!    1 for �(z)
!
!    Output, real ( kind = 8 ) GR, GI, the real and imaginary parts of
!    the selected function.
!
  implicit none

  real ( kind = 8 ), save, dimension ( 10 ) :: a = (/ &
    8.333333333333333D-02, -2.777777777777778D-03, &
    7.936507936507937D-04, -5.952380952380952D-04, &
    8.417508417508418D-04, -1.917526917526918D-03, &
    6.410256410256410D-03, -2.955065359477124D-02, &
    1.796443723688307D-01, -1.39243221690590D+00 /)
  real ( kind = 8 ) g0
  real ( kind = 8 ) gi
  real ( kind = 8 ) gi1
  real ( kind = 8 ) gr
  real ( kind = 8 ) gr1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kf
  integer ( kind = 4 ) na
  real ( kind = 8 ) pi
  real ( kind = 8 ) si
  real ( kind = 8 ) sr
  real ( kind = 8 ) t
  real ( kind = 8 ) th
  real ( kind = 8 ) th1
  real ( kind = 8 ) th2
  real ( kind = 8 ) x
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2

  pi = 3.141592653589793D+00

  if ( y == 0.0D+00 .and. x == int ( x ) .and. x <= 0.0D+00 ) then
    gr = 1.0D+300
    gi = 0.0D+00
    return
  else if ( x < 0.0D+00 ) then
    x1 = x
    y1 = y
    x = -x
    y = -y
  end if

  x0 = x

  if ( x <= 7.0D+00 ) then
    na = int ( 7 - x )
    x0 = x + na
  end if

  z1 = sqrt ( x0 * x0 + y * y )
  th = atan ( y / x0 )
  gr = ( x0 - 0.5D+00 ) * log ( z1 ) - th * y - x0 &
    + 0.5D+00 * log ( 2.0D+00 * pi )
  gi = th * ( x0 - 0.5D+00 ) + y * log ( z1 ) - y

  do k = 1, 10
    t = z1 ** ( 1 - 2 * k )
    gr = gr + a(k) * t * cos ( ( 2.0D+00 * k - 1.0D+00 ) * th )
    gi = gi - a(k) * t * sin ( ( 2.0D+00 * k - 1.0D+00 ) * th )
  end do

  if ( x <= 7.0D+00 ) then
    gr1 = 0.0D+00
    gi1 = 0.0D+00
    do j = 0, na - 1
      gr1 = gr1 + 0.5D+00 * log ( ( x + j ) ** 2 + y * y )
      gi1 = gi1 + atan ( y / ( x + j ) )
    end do
    gr = gr - gr1
    gi = gi - gi1
  end if

  if ( x1 < 0.0D+00 ) then
    z1 = sqrt ( x * x + y * y )
    th1 = atan ( y / x )
    sr = - sin ( pi * x ) * cosh ( pi * y )
    si = - cos ( pi * x ) * sinh ( pi * y )
    z2 = sqrt ( sr * sr + si * si )
    th2 = atan ( si / sr )
    if ( sr < 0.0D+00 ) then
      th2 = pi + th2
    end if
    gr = log ( pi / ( z1 * z2 ) ) - gr
    gi = - th1 - th2 - gi
    x = x1
    y = y1
  end if

  if ( kf == 1 ) then
    g0 = exp ( gr )
    gr = g0 * cos ( gi )
    gi = g0 * sin ( gi )
  end if

  return
end
