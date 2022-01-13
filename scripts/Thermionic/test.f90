program test
  implicit none

  real :: Chem,Phi,Fmax,Fmin,Field,Tmax,Tmin,Temp
  integer :: Nx

  namelist/JftParam/ Chem, Phi, Fmax, Fmin, Field, &
                     & Tmax, Tmin, Temp, Nx

  Chem = 7.00000
  PHI = 4.50000
  FMAX = 8.00000
  FMIN = 1.00000
  FIELD = 2.00000
  TMAX = 1700.00
  TMIN = 300.000
  TEMP = 1000.00
  NX = 41

  write(*,nml=JftParam)
end program test
