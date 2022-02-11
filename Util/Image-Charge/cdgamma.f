! complex Gamma function in double precision
! http://www.kurims.kyoto-u.ac.jp/~ooura/gamerf.html
!
      function cdgamma(x)
      implicit real*8 (a - h, o - z)
      complex*16 cdgamma, x
      parameter (
     &    pi = 3.14159265358979324d+00,
     &    pv = 7.31790632447016203d+00,
     &    pu = 3.48064577727581257d+00,
     &    pr = 3.27673720261526849d-02,
     &    p1 = 1.05400280458730808d+01,
     &    p2 = 4.73821439163096063d+01,
     &    p3 = 9.11395751189899762d+01,
     &    p4 = 6.62756400966213521d+01,
     &    p5 = 1.32280130755055088d+01,
     &    p6 = 2.93729529320536228d-01)
      parameter (
     &    q1 = 9.99999999999975753d-01,
     &    q2 = 2.00000000000603851d+00,
     &    q3 = 2.99999999944915534d+00,
     &    q4 = 4.00000003016801681d+00,
     &    q5 = 4.99999857982434025d+00,
     &    q6 = 6.00009857740312429d+00)
      xr = dreal(x)
      xi = dimag(x)
      if (xr .lt. 0) then
          wr = 1 - xr
          wi = -xi
      else
          wr = xr
          wi = xi
      end if
      ur = wr + q6
      vr = ur * (wr + q5) - wi * wi
      vi = wi * (wr + q5) + ur * wi
      yr = p6 + (p5 * ur + p4 * vr)
      yi = p5 * wi + p4 * vi
      ur = vr * (wr + q4) - vi * wi
      ui = vi * (wr + q4) + vr * wi
      vr = ur * (wr + q3) - ui * wi
      vi = ui * (wr + q3) + ur * wi
      yr = yr + (p3 * ur + p2 * vr)
      yi = yi + (p3 * ui + p2 * vi)
      ur = vr * (wr + q2) - vi * wi
      ui = vi * (wr + q2) + vr * wi
      vr = ur * (wr + q1) - ui * wi
      vi = ui * (wr + q1) + ur * wi
      yr = yr + (p1 * ur + vr)
      yi = yi + (p1 * ui + vi)
      ur = vr * wr - vi * wi
      ui = vi * wr + vr * wi
      t = ur * ur + ui * ui
      vr = (yr * ur + yi * ui) + pr * t
      vi = yi * ur - yr * ui
      yr = wr + pv
      ur = 0.5d0 * log(yr * yr + wi * wi) - 1
      ui = atan2(wi, yr)
      yr = exp(ur * (wr - 0.5d0) - ui * wi - pu) / t
      yi = ui * (wr - 0.5d0) + ur * wi
      ur = yr * cos(yi)
      ui = yr * sin(yi)
      yr = ur * vr - ui * vi
      yi = ui * vr + ur * vi
      if (xr .lt. 0) then
          wr = pi * xr
          wi = exp(pi * xi)
          vi = 1 / wi
          ur = (vi + wi) * sin(wr)
          ui = (vi - wi) * cos(wr)
          vr = ur * yr + ui * yi
          vi = ui * yr - ur * yi
          ur = 2 * pi / (vr * vr + vi * vi)
          yr = ur * vr
          yi = ur * vi
      end if
      cdgamma = cmplx(yr, yi)
      end
!
