!–––––––––––––––––––––––––––––––––––––––––––––––-
module Fundamentals
!Source: http://physics.nist.gov/constants
! (converted from MKSA)
!Units:
! electron volts, nanometers,
! femtoseconds, charge q = 1
! terms with different units are combinations of Units
! eg: [rmo] = [eV]/c ̂ 2, [c] = nm/fs
!––––––––––––––––––––––––––––––––––––––––––––––––-
real, parameter :: pi = 3.14159265359, &
& rmo = 5.685629853, rkb = 1.0/11604.5192, &
& rhbar = 0.658211899, &
& alpha = 7.29735257E-3, c = 299.7924580, &
& epso = 1.0/(4*pi*alpha*rhbar*c), &
& Qo = alpha*rhbar*c/4, Ampcm2 = 1.602176565E10, &
& Arld = (rmo*rkb**2)/(2*rhbar*(pi*rhbar)**2)
!––––––––––––––––––––––––––––––––––––––––––––––––-
end module Fundamentals
!––––––––––––––––––––––––––––––––––––––––––––––––-
!––––––––––––––––––––––––––––––––––––––––––––––––-
module MaterialEnvironment
!––––––––––––––––––––––––––––––––––––––––––––––––-
real :: Chem,Phi,Fmax,Fmin,Field,Tmax,Tmin,Temp
integer :: Nx
character*1:: tb = char(9)
data Chem, Phi, Nx / 7.0, 4.5, 41 /
data Fmax, Fmin, Field, Tmax, Tmin, Temp &
& / 8.0, 1.0, 2.0, 1700.0, 300.0, 1000.0 /
!––––––––––––––––––––––––––––––––––––––––––––––––-
end module MaterialEnvironment
!––––––––––––––––––––––––––––––––––––––––––––––––-
!==============================================================
PROGRAM GenThermField
!==============================================================
!.....Calculation of Generalized Thermal
! Field Current Density
!.....Units: Field = eV/nm; Temp = Kelvin;
! J = Amp/cm ̂ 2
!.....When not specified, units are eV-nm-fs-q
!––––––––––––––––––––––––––––––––––––––––––––––––-
use MaterialEnvironment !parameters used in subroutines
real :: Tmp, Fld
!....Read in needed parameters
call NameDropping
!––––––––––––––––––––––––––––––––––––––––––––––––-
!....Create temperature variation output
!....Let ln(F(j)) be evenly spaced
!––––––––––––––––––––––––––––––––––––––––––––––––-
write(*,*) 'Creating file Jfield.txt:  field variation'
write(*,*) '  Temp taken to be T = ',Temp,' Kelvin'
open(unit=11,file='Jfield.txt')
write(11,*) 'F [eV/nm]',tb,'JGTF [A/cm ̂ 2]',tb, &
& 'JRLD',tb,'JFN'
do i = 1, Nx
Fld = Fmin*(Fmax/Fmin)**(float(i-1)/(Nx-1))
write(11,*) Fld,tb,RJgtf(Phi,Chem,Fld,Temp), tb, &
& RJrldo(Phi,Fld,Temp),tb,RJfno(Phi,Fld)
end do
close(unit=11)
!––––––––––––––––––––––––––––––––––––––––––––––––-
!....Create field variation output
!....Let T(j) be evenly spaced
!––––––––––––––––––––––––––––––––––––––––––––––––-
write(*,*) 'Creating file Jtemp.txt: thermal variation'
write(*,*) '  Field taken to be F = ',Field,' eV/nm'
open(unit=12,file='JTemp.txt')
write(12,*) 'T [K]',tb,'JGTF [A/cm ̂ 2]',tb,&
& 'JRLD',tb,'JFN'
do i = 1, Nx
Tmp = Tmin + (Tmax-Tmin)*float(i-1)/(Nx-1)
write(12,*) Tmp,tb,RJgtf(Phi,Chem,Field,Tmp), tb, &
& RJrldo(Phi,Field,Tmp),tb,RJfno(Phi,Field)
end do
close(unit=12)
!––––––––––––––––––––––––––––––––––––––––––––––––-
END PROGRAM GenThermField
!==============================================================
!==============================================================
subroutine NameDropping
!READING/CREATING OF DATA INPUT FILE
!––––––––––––––––––––––––––––––––––––––––––––––––-
use MaterialEnvironment
!––––––––––––––––––––––––––––––––––––––––––––––––-
implicit none
integer :: istat, idebug
namelist/JftParam/ Chem, Phi, Fmax, Fmin, Field, &
& Tmax, Tmin, Temp, Nx
write(*,*) 'PROGRAM GenFNRLD'
write(*,*) '  <Jfrin.nml> being opened...'
open(unit=11,file='Jfrin.nml',status='old',iostat=istat)
close(unit=11)
if(istat == 0) then !A namelist file exists; read it
open(unit=11,file='Jfrin.nml')
read(11,JftParam)
close(unit=11)
if(idebug>0) write(*,JftParam)
else !A namelist file does not exist; create it
write(*,*) '  input file Jfrin.NML not found.'
write(*,*) '  Creating input data file now...'
open(unit=11,file='Jfrin.nml')
write(11,*) 'Symbol    Definition            Unit'
write(11,*) 'Chem      Chemical Potential    [eV]'
write(11,*) 'Phi       Work Function         [eV]'
write(11,*) 'Fmax      Max field for Jfield  [eV/nm]'
write(11,*) 'Fmin      Min field for Jfield  [eV/nm]'
write(11,*) 'Field     Field for for Jtemp   [eV/nm]'
write(11,*) 'Tmax      Max temp  for Jtemp   [Kelvin]'
write(11,*) 'Tmin      Min temp  for Jtemp   [Kelvin]'
write(11,*) 'Tmin      Temp for  for Jfield  [Kelvin]'
write(11,*) 'Nx        Number of plot points'
write(11,JftParam)
close(unit=11)
stop
end if
return
!––––––––––––––––––––––––––––––––––––––––––––––––-
! END subroutine NameDropping
! READING/CREATING OF DATA INPUT
!==============================================================
end
!==============================================================
function RJfno(Phi,Field)
!––––––––––––––––––––––––––––––––––––––––––––––––-
!Internal Units: eV-nm-fs-q
!Implicit use is made of the Forbes-Deane
! Approximation to v(y)
!The Nordheim paraemeter t(y) is approximated by tyo
!Output Units: Amps / cm ̂ 2
!––––––––––––––––––––––––––––––––––––––––––––––––-
use Fundamentals
implicit none
real, parameter :: Ao = 1.0/(16*rhbar*pi**2), &
& Bo = 4*sqrt(2*rmo)/(3*rhbar), &
& tyo = 1.0613132
real :: Phi, Field, RJfno, nu
nu = 2*Bo*Qo/(3*sqrt(Phi))
RJfno = (Ao/(Phi*tyo**2))*(((Phi**2)*exp(6.0)/&
& (4*Qo))**nu)*(Field**(2.0-nu))* &
& exp(-Bo*(Phi**1.5)/Field)*Ampcm2
return
!––––––––––––––––––––––––––––––––––––––––––––––––-
end
!==============================================================
!==============================================================
function RJrldo(Phi,Field,Temp)
!––––––––––––––––––––––––––––––––––––––––––––––––-
!Internal Units: eV-nm-fs-q
!Output Units: Amps / cm ̂ 2
!––––––––––––––––––––––––––––––––––––––––––––––––-
use Fundamentals
implicit none
real :: Phi, phix, Field, Temp, RJrldo
phix = Phi - sqrt(4*Qo*Field)
RJrldo = Arld*(Temp**2)*exp(-phix/(rkb*Temp))*Ampcm2
return
!––––––––––––––––––––––––––––––––––––––––––––––––-
end
!==============================================================
!==============================================================
function RJgtf(Phi,Chem,Fld,Tmp)
!Internal Units: eV-nm-fs-q
!Output Units: Amps / cm ̂ 2
!––––––––––––––––––––––––––––––––––––––––––––––––-
use Fundamentals
implicit none
real :: y, phix, vy, ty, Qmu, Bfmu, Bfmphi, Tmin, &
& Tmax, Tmp, Em, s, rn, RJgtf, Phi, Chem, &
& Fld, BetaT, &
& RNns, Ap, Bp, Cp, p, ThetaE, BetaE
! Create terms needed for the evaluations
betaT = 1.0/(rkb*Tmp)
y = sqrt(4*Qo*Fld)/Phi
phix = Phi - sqrt(4*Qo*Fld)
vy = 1.0 - (3.0 - alog(y))*(y**2)/3
ty = 1.0 + (1.0 - alog(y))*(y**2)/9
Qmu = 4*sqrt(2*rmo*Phi**3)*vy/(3*rhbar*Fld)
Bfmu = 2*sqrt(2*rmo*Phi)*ty/( rhbar*Fld)
Bfmphi = pi*sqrt(y*rmo*Phi) /( rhbar*Fld)
Tmin = (rhbar*Fld/(4*ty*rkb))*sqrt(2.0/(rmo*Phi))
Tmax = (rhbar*Fld/(pi*rkb)) *sqrt(1.0/(y*rmo*Phi))
! Three regimes:
if(Tmp.le.Tmin) then !Field Emission Regime
s = Qmu
rn = betaT/Bfmu
else if(Tmp.gt.Tmax) then !Thermal Emission Regime
s = Bfmphi*Phix
rn = betaT/BFmphi
else !Intermediate Regime
rn = 1.0
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
end
!––––––––––––––––––––––––––––––––––––––––––––––––-
!––––––––––––––––––––––––––––––––––––––––––––––––-
function RNns(n,s)
! (7*pi ̂ 4 - 720)/720 = -0.10593434
! ( pî 2 - 12)/6 = -0.35506593
!––––––––––––––––––––––––––––––––––––––––––––––––-
implicit none
real :: n, s, RNns, sn, sd, x, y, z, sng
x = n**2
y = 1.0/x
sn = -x*(0.10593434*x + 0.35506593)
sd = -y*(0.10593434*y + 0.35506593)
z = (n-1.0)*s
if(abs(z).gt.0.001) then
sng = (n**2 + 1.0)*((n**2)*exp(-s) - &
& exp(-n*s))/(n**2 - 1.0)
else
sng = (0.5*(n**2 + 1.0)*exp(-s)/(n+1.0))* &
& ((1.0-n)*(s**2) + 2*(1.0+n) + 2*s)
end if
RNns = sng + sn*exp(-n*s) + x*sd*exp(-s)
return
!––––––––––––––––––––––––––––––––––––––––––––––––-
end
!==============================================================
!==============================================================
function ThetaE(Phi,Chem,Field,E)
! Cubic in p = (E - mu)/phix
!––––––––––––––––––––––––––––––––––––––––––––––––-
use Fundamentals
implicit none
real :: Qmu, Bfmu, Bfmphi, p, ThetaE, Phi, Chem, &
& Field, E, vy, ty, y, phix
y = sqrt(4*Qo*Field)/Phi
phix = Phi - sqrt(4*Qo*Field)
vy = 1.0 - (3.0 - alog(y))*(y**2)/3
ty = 1.0 + (1.0 - alog(y))*(y**2)/9
Qmu = 4*sqrt(2*rmo*Phi**3)*vy/(3*rhbar*Field)
Bfmu = 2*sqrt(2*rmo*Phi)*ty/( rhbar*Field)
Bfmphi = pi*sqrt(y*rmo*Phi) /( rhbar*Field)
p = (E - Chem)/phix
ThetaE = Qmu*(2*p+1.0)*(1.0 - p)**2 - &
& phix*Bfmu*p*(1.0 - p)**2 + &
& phix*Bfmphi*(1.0 - p)*p**2
return
!––––––––––––––––––––––––––––––––––––––––––––––––-
end
!==============================================================
!==============================================================
function BetaE(Phi,Chem,Field,E)
! Quadratic in p = (E - mu)/phix
!––––––––––––––––––––––––––––––––––––––––––––––––-
use Fundamentals
implicit none
real :: Qmu, Bfmu, Bfmphi, p, BetaE, Phi, Chem, &
& Field, E, vy, ty, y, phix
y = sqrt(4*Qo*Field)/Phi
phix = Phi - sqrt(4*Qo*Field)
vy = 1.0 - (3.0 - alog(y))*(y**2)/3
ty = 1.0 + (1.0 - alog(y))*(y**2)/9
Qmu = 4*sqrt(2*rmo*Phi**3)*vy/(3*rhbar*Field)
Bfmu = 2*sqrt(2*rmo*Phi)*ty/( rhbar*Field)
Bfmphi = pi*sqrt(y*rmo*Phi) /( rhbar*Field)
p = (E - Chem)/phix
BetaE = (6*Qmu/phix)*p*(1.0 - p) - &
& Bfmu*(1.0 - p)*(3*p - 1.0) + &
& Bfmphi*p*(3*p - 2.0)
return
!––––––––––––––––––––––––––––––––––––––––––––––––-
end
!==============================================================