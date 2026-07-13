! Translation of Kevin Jensen's newer oGTF Matlab code (see the Kevin/ folder:
! FUNdamental.m, FUNgtfns.m, FUNnsu.m) for the general thermal-field emission
! current density J_gtf(F, T, phi).
!
! Differences from the older rGTF code in mod_kevin_rjgtf.F90:
!  - N(n,s) (FUNnsu.m) now has an exact n = 1 branch, N = (s+1)*exp(-s), which
!    is the one used in the intermediate regime, and a high-field floor
!    N >= n^2*exp(-s) (the infinite-n, positive-s limit) that prevents the
!    roll-over to negative values the rGTF method suffered from.
!  - In the intermediate regime theta(E_m) is evaluated directly from the
!    quadratic root p_o instead of through the ThetaE/BetaE helper functions
!    (analytically equivalent, fewer operations).
!  - Physical constants updated to CODATA 2018 values.
!  - Regime boundaries are T < Tmin / T > Tmax (rGTF used T <= Tmin), so the
!    boundary points now fall in the intermediate regime.
!
! Validity: the model assumes y = sqrt(4*Qo*F)/phi < 1, i.e. the Schottky
! lowering does not exceed the work function (F < phi^2/(4*Qo), about
! 14 eV/nm for phi = 4.5 eV).

module mod_kevin_rjgtf_v2
implicit none
PRIVATE
PUBLIC :: Get_Kevin_Jgtf_v2, oGTF_ns

!Units (FUNdamental.m):
! electron volts, nanometers,
! femtoseconds, charge q = 1
! terms with different units are combinations of Units
! eg: [mo] = [eV]/c^2, [c] = nm/fs
!--------------------------------------------------------------
  double precision, parameter :: pi   = 3.14159265358979324d0
  double precision, parameter :: kb   = 1.0d0/11604.50635d0    ! Boltzmann's const [eV/K]
  double precision, parameter :: hbar = 0.6582119571d0         ! Planck's const / 2 pi [eV fs]
  double precision, parameter :: c    = 299.7924580d0          ! speed of light [nm/fs]
  double precision, parameter :: mo   = 5.685630103d0          ! e mass [eV (fs/nm)^2]
  double precision, parameter :: afs  = 1.0d0/137.035999084d0  ! fine structure const
  double precision, parameter :: Qo   = afs*hbar*c/4.0d0       ! image charge factor [eV nm]

  ! Unit conversions: cm in nm, Amp in q/fs (= Coulomb/q / fs)
  double precision, parameter :: cm  = 1.0d7
  double precision, parameter :: Amp = 6.241509074d3

  ! Richardson-Laue-Dushman constant [A cm^-2 K^-2]
  double precision, parameter :: Arld = (mo*kb**2/(2.0d0*pi**2*hbar**3))*cm**2/Amp
!--------------------------------------------------------------

contains

! Wrapper with the same conventions as Get_Kevin_Jgtf in mod_kevin_rjgtf.F90:
! Field: V/m -> eV/nm
! Current density: A/cm^2 -> A/m^2
! Note that in Kevin's units the electron charge is e = 1.
! This means that 1 V = 1 eV and that 1 GV/m = 1 eV/nm.
! See Chapter 2 page 10 in Kevin's book Introduction to the Physics of Electron Emission.
double precision function Get_Kevin_Jgtf_v2(F, T, w_theta)
  double precision, intent(in) :: F       ! Field [V/m]
  double precision, intent(in) :: T       ! Temperature [K]
  double precision, intent(in) :: w_theta ! Work function [eV]
  double precision, parameter  :: chem = 7.0d0 ! Chemical potential [eV], cancels out of s
  double precision             :: F_evnm, Jgtf, nft, sft

  F_evnm = abs(F) * 1.0d-9 ! Kevin wants the field in eV/nm see note above (Basically we convert to GV/m).

  ! Guard against a negligible field: oGTF divides by the field and takes log(y) with
  ! y ~ sqrt(field), both of which blow up to NaN/Inf as field -> 0. No field means no emission.
  if (F_evnm < 1.0d-9) then
    Get_Kevin_Jgtf_v2 = 0.0d0
    return
  end if

  call oGTF_ns(F_evnm, T, w_theta, chem, Jgtf, nft, sft)
  Get_Kevin_Jgtf_v2 = Jgtf * 1.0d4 ! Kevin returns in A/cm^2 convert to A/m^2
end function Get_Kevin_Jgtf_v2

!==============================================================

!--------------------------------------------------------------
! oGTF formalism (FUNgtfns.m)
! for calculating n, s and Jgtf, and optionally Jrld and Jfn
!Internal Units: eV-nm-fs-q
!Output Units: Amps / cm^2
!--------------------------------------------------------------
subroutine oGTF_ns(Fo, To, Phi, chem, Jgtf, nft, sft, Jrld, Jfn)
  double precision, intent(in)            :: Fo   ! Field [eV/nm]
  double precision, intent(in)            :: To   ! Temperature [K]
  double precision, intent(in)            :: Phi  ! Work function [eV]
  double precision, intent(in)            :: chem ! Chemical potential [eV]
  double precision, intent(out)           :: Jgtf ! GTF current density [A/cm^2]
  double precision, intent(out)           :: nft, sft ! n and s in N(n,s)
  double precision, intent(out), optional :: Jrld ! Richardson-Laue-Dushman [A/cm^2]
  double precision, intent(out), optional :: Jfn  ! Fowler-Nordheim [A/cm^2]

  double precision :: yo, phix, ty, vy, Tmin, Tmax
  double precision :: betaT, betau, betap, theto
  double precision :: Ap, Bp, Cp, po, Em, theta

  yo   = sqrt(4.0d0*Qo*Fo)/Phi
  phix = Phi - sqrt(4.0d0*Qo*Fo) ! Barrier reduced by the Schottky lowering
  if (yo >= 1.0d0) then
    print *, 'mod_kevin_rjgtf_v2: Error: y >= 1.0 (field too high for the oGTF model)'
  end if

  ty = 1.0d0 + (yo**2)*(1.0d0 - log(yo))/9.0d0
  vy = 1.0d0 - (yo**2)*(3.0d0 - log(yo))/3.0d0

  Tmin = (hbar*Fo/(4.0d0*kb*ty))*sqrt(2.0d0/(mo*Phi))
  Tmax = hbar*Fo/(kb*pi*sqrt(mo*Phi*yo))

  betaT = 1.0d0/(kb*To)
  betau = (2.0d0/(hbar*Fo))*sqrt(2.0d0*mo*Phi)*ty
  betap = (pi/(hbar*Fo))*sqrt(mo*Phi*yo)
  theto = (4.0d0*sqrt(2.0d0*mo*Phi**3)/(3.0d0*hbar*Fo))*vy

  ! Three regimes:
  if (To < Tmin) then ! Field Emission Regime
    nft = betaT/betau
    sft = theto
  else if (To > Tmax) then ! Thermal Emission Regime
    nft = betaT/betap
    sft = betap*phix
  else ! Intermediate Regime
    ! Root of quadratic Cp + Bp*p + Ap*p^2 = 0
    Ap = 3.0d0*(betap + betau) - 6.0d0*theto/phix
    Bp = -2.0d0*(betap + 2.0d0*betau) + 6.0d0*theto/phix
    Cp = betau - betaT
    po = (-Bp - sqrt(Bp**2 - 4.0d0*Ap*Cp))/(2.0d0*Ap)

    Em = chem + po*phix
    theta = ((1.0d0 - po)**2)*(2.0d0*po + 1.0d0)*theto - &
          & phix*po*(1.0d0 - po)*((1.0d0 - po)*betau - po*betap)
    nft = 1.0d0
    sft = theta + betaT*(Em - chem)
  end if

  Jgtf = Arld*Nns(nft, sft)*To**2
  if (present(Jrld)) Jrld = Arld*exp(-betaT*phix)*To**2
  if (present(Jfn))  Jfn  = (Arld/(kb*betau)**2)*exp(-theto)
end subroutine oGTF_ns

!--------------------------------------------------------------
! Function FUNnsu.m: N(n,s)
! High field correction by introduction of sfn (n >> 1) limit
! this is done to prevent roll-over to negative values associated
! with s becoming small or negative in the rGTF method
! (7*pi^4 - 720)/720  = -0.10593434 (x2)
! (  pi^2 -  12)/6    = -0.35506593
!--------------------------------------------------------------
double precision function Nns(n, s)
  double precision, intent(in) :: n, s
  double precision :: x, y, z, sn, sd, sng, sfn

  if (n == 1.0d0) then
    Nns = (s + 1.0d0)*exp(-s)
  else
    x = n**2
    y = 1.0d0/x
    ! The exact sng expression is a 0/0 form as n -> 1 (T at the Tmin/Tmax
    ! regime boundaries); switch to its series expansion in z = (n-1)*s near
    ! n = 1 to avoid the cancellation. Not in the Matlab code; RNns in
    ! mod_kevin_rjgtf.F90 does the same with |z| > 0.001, but 1d-5 is where
    ! the series truncation error meets the cancellation error (~1e-11).
    z = (n - 1.0d0)*s
    if (abs(z) > 1.0d-5) then
      sng = (x + 1.0d0)*(x*exp(-s) - exp(-n*s))/(x - 1.0d0)
    else
      sng = (0.5d0*(x + 1.0d0)*exp(-s)/(n + 1.0d0)) * &
          & ((1.0d0 - n)*s**2 + 2.0d0*(1.0d0 + n) + 2.0d0*s)
    end if
    sn = -x*(0.10593434d0*x + 0.35506593d0)
    sd = -y*(0.10593434d0*y + 0.35506593d0)
    sfn = x*exp(-s) ! infinite n, positive s limit
    Nns = max(sng + sn*exp(-n*s) + x*sd*exp(-s), sfn)
  end if
end function Nns
!==============================================================

end module mod_kevin_rjgtf_v2
