! Code from Kevin Jensen's books
! RJgtf subroutine pulled out

module mod_kevin_rjgtf

implicit none
PRIVATE
PUBLIC :: Get_Kevin_Jgtf

!Source: http://physics.nist.gov/constants
! (converted from MKSA)
!Units:
! electron volts, nanometers,
! femtoseconds, charge q = 1
! terms with different units are combinations of Units
! eg: [rmo] = [eV]/c ̂ 2, [c] = nm/fs
!––––––––––––––––––––––––––––––––––––––––––––––––-
    double precision, parameter :: pi = 3.14159265359d0, &
    & rmo = 5.685629853d0, rkb = 1.0d0/11604.5192d0, &
    & rhbar = 0.658211899d0, &
    & alpha = 7.29735257d-3, c = 299.7924580d0, &
    & epso = 1.0d0/(4*pi*alpha*rhbar*c), &
    & Qo = alpha*rhbar*c/4, Ampcm2 = 1.602176565d10, &
    & Arld = (rmo*rkb**2)/(2*rhbar*(pi*rhbar)**2)
    !––––––––––––––––––––––––––––––––––––––––––––––––-

contains

! Call this function to convert to and from the units Kevin uses
! Field: V/m -> eV/nm
! Current density: A/cm² -> A/m²
! Note that in Kevins units the electron charge is e = 1.
! This means the 1 V = 1 eV and that 1 GV/m = 1 eV/nm.
! See Chapter 2 page 10 in Kevins book Introduction to the Physics of Electrom Emission.
double precision function Get_Kevin_Jgtf(F, T, w_theta)
  double precision, intent(in) :: F, T, w_theta
  double precision             :: Chem = 7.0d0 ! Chemical potential
  double precision             :: F_evnm

  F_evnm = F * 1.0d-9 ! Kevin wants the field in eV/nm see note above (Basically we convert to GV/m).
  print *, 'Enter Kevin'
  Get_Kevin_Jgtf = RJgtf(w_theta, Chem, F, T) * (1.0d0/1.0d-2)**2 ! Kevin returns in A/cm² convert to A/m²
  print *, Get_Kevin_Jgtf
  print *, 'Exit Kevin'
  print *, ''
end function Get_Kevin_Jgtf

!==============================================================

!––––––––––––––––––––––––––––––––––––––––––––––––-
double precision function RNns(n,s)
! (7*pi^4 - 720)/720 = -0.10593434
! ( pi^2 - 12)/6 = -0.35506593
!––––––––––––––––––––––––––––––––––––––––––––––––-
  implicit none
  double precision :: n, s, sn, sd, x, y, z, sng

  x = n**2
  y = 1.0/x
  sn = -x*(0.10593434d0*x + 0.35506593d0)
  sd = -y*(0.10593434d0*y + 0.35506593d0)
  z = (n-1.0d0)*s
  if(abs(z) .gt. 0.001d0) then
      sng = (n**2 + 1.0d0)*((n**2)*exp(-s) - &
      &     exp(-n*s))/(n**2 - 1.0d0)
  else
      sng = (0.5d0*(n**2 + 1.0d0)*exp(-s)/(n+1.0d0))* &
          & ((1.0d0-n)*(s**2) + 2*(1.0d0+n) + 2*s)
  end if

  RNns = sng + sn*exp(-n*s) + x*sd*exp(-s)

  return
!––––––––––––––––––––––––––––––––––––––––––––––––-
end function RNns

double precision function ThetaE(Phi,Chem,Field,E)
! Cubic in p = (E - mu)/phix
!––––––––––––––––––––––––––––––––––––––––––––––––-
  implicit none
  double precision :: Qmu, Bfmu, Bfmphi, p, Phi, Chem, &
                    & Field, E, vy, ty, y, phix

  y = sqrt(4*Qo*Field)/Phi
  phix = Phi - sqrt(4*Qo*Field)
  vy = 1.0 - (3.0d0 - log(y))*(y**2)/3
  ty = 1.0 + (1.0d0 - log(y))*(y**2)/9
  Qmu = 4*sqrt(2*rmo*Phi**3)*vy/(3*rhbar*Field)
  Bfmu = 2*sqrt(2*rmo*Phi)*ty/( rhbar*Field)
  Bfmphi = pi*sqrt(y*rmo*Phi) /( rhbar*Field)
  p = (E - Chem)/phix
  ThetaE = Qmu*(2*p+1.0d0)*(1.0d0 - p)**2 - &
  & phix*Bfmu*p*(1.0d0 - p)**2 + &
  & phix*Bfmphi*(1.0d0 - p)*p**2

  return
!––––––––––––––––––––––––––––––––––––––––––––––––-
end function ThetaE

double precision function BetaE(Phi,Chem,Field,E)
  ! Quadratic in p = (E - mu)/phix
  !––––––––––––––––––––––––––––––––––––––––––––––––-
  implicit none
  double precision :: Qmu, Bfmu, Bfmphi, p, Phi, Chem, &
        & Field, E, vy, ty, y, phix

  y = sqrt(4*Qo*Field)/Phi
  phix = Phi - sqrt(4*Qo*Field)
  vy = 1.0d0 - (3.0d0 - log(y))*(y**2)/3
  ty = 1.0d0 + (1.0d0 - log(y))*(y**2)/9
  Qmu = 4*sqrt(2*rmo*Phi**3)*vy/(3*rhbar*Field)
  Bfmu = 2*sqrt(2*rmo*Phi)*ty/( rhbar*Field)
  Bfmphi = pi*sqrt(y*rmo*Phi) /( rhbar*Field)
  p = (E - Chem)/phix
  BetaE = (6*Qmu/phix)*p*(1.0d0 - p) - &
  & Bfmu*(1.0d0 - p)*(3*p - 1.0d0) + &
  & Bfmphi*p*(3*p - 2.0d0)

  return
!––––––––––––––––––––––––––––––––––––––––––––––––-
end

double precision function RJgtf(Phi,Chem,Fld,Tmp)
!Internal Units: eV-nm-fs-q
!Output Units: Amps / cm ̂ 2
!––––––––––––––––––––––––––––––––––––––––––––––––-
  implicit none
  double precision :: y, phix, vy, ty, Qmu, Bfmu, Bfmphi, Tmin, &
                    & Tmax, Tmp, Em, s, rn, Phi, Chem, &
                    & Fld, BetaT, &
                    & Ap, Bp, Cp, p

  ! Create terms needed for the evaluations
  betaT = 1.0d0/(rkb*Tmp)
  y = sqrt(4*Qo*Fld)/Phi
  phix = Phi - sqrt(4*Qo*Fld)
  vy = 1.0d0 - (3.0d0 - log(y))*(y**2)/3
  ty = 1.0d0 + (1.0d0 - log(y))*(y**2)/9
  Qmu = 4*sqrt(2*rmo*Phi**3)*vy/(3*rhbar*Fld)
  Bfmu = 2*sqrt(2*rmo*Phi)*ty/( rhbar*Fld)
  Bfmphi = pi*sqrt(y*rmo*Phi) /( rhbar*Fld)
  Tmin = (rhbar*Fld/(4*ty*rkb))*sqrt(2.0d0/(rmo*Phi))
  Tmax = (rhbar*Fld/(pi*rkb)) *sqrt(1.0d0/(y*rmo*Phi))

  ! Three regimes:
  if(Tmp.le.Tmin) then !Field Emission Regime
    s = Qmu
    rn = betaT/Bfmu
  else if(Tmp.gt.Tmax) then !Thermal Emission Regime
    s = Bfmphi*Phix
    rn = betaT/BFmphi
  else !Intermediate Regime
    rn = 1.0d0
    Ap = 3*(phix*(Bfmu + Bfmphi) - 2*Qmu)
    Bp = 2*(3*Qmu - phix*(2*Bfmu + Bfmphi))
    Cp = phix*(Bfmu - BetaT)
    p = (-Bp-sqrt(Bp**2 - 4*Ap*Cp))/(2*Ap)
    Em = Chem + p*phix
    s = betaT*(Em - Chem + &
      & (ThetaE(Phi,Chem,Fld,Em)/BetaE(Phi,Chem,Fld,Em)))
  end if

  RJgtf = (Arld/(rkb*betaT)**2)*RNns(rn,s)*Ampcm2

  return
!––––––––––––––––––––––––––––––––––––––––––––––––-
end function RJgtf
!==============================================================

end module mod_kevin_rjgtf