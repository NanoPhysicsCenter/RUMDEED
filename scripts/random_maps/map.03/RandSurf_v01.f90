! generates a random landscape of a binary function
!=====================================================================
  MODULE Mod_Map
   implicit none
   real*8, parameter :: pi=3.14159265358979324
!  lattice of sites
   real*8, parameter :: Lx=1500, Ly=1000  ! lateral dimensions of the rectangle 
   integer, parameter:: Nx=300, Ny=200, NN=Nx*Ny  ! number of lattice points 
   real*8, parameter :: dx=Lx/(Nx-1), dy=Ly/(Ny-1) ! lattice spacing 
   real*8 site(NN,3)  ! coordinates x,y and color of each site 
!  islands 
   real*8, parameter  :: avamp=30, viamp=20 ! av. amplitude and var. of sinuses
   real*8, parameter  :: avfrq=1.5, vifrq=0.5 ! av. frequency factor and var.
   integer, parameter :: nmax=5  ! max number of sin terms (or cos)
   integer, parameter :: nisl=50 ! number of islands
!  integer, parameter :: nisl=1 ! number of islands
!  real*8 amp(0:2*nmax,nisl),frq(2*nmax,nisl),phi0(nisl),xc(nisl),yc(nisl)
   real*8 amp(0:nmax,nisl),frq(nmax,nisl),phi0(0:nmax,nisl),xc(nisl),yc(nisl)
  END MODULE Mod_Map
!=====================================================================
  PROGRAM Main
    use Mod_Map
    implicit none
    real*8 x,y,Radius,rad,phi,xi,eta
    integer i,j,nsite
!   dx=Lx/(Nx-1); dy=Ly/(Ny-1)
!   write(*,*) dx,dy
    nsite=0
    do i=1,Nx; x=-Lx/2+(i-1)*dx 
      do j=1,Ny; y=-Ly/2+(j-1)*dy
        nsite=nsite+1
        site(nsite,1)=x; site(nsite,2)=y; site(nsite,3)=0
      enddo 
    enddo
!   write(*,*) nsite,NN
    open(10,file='map.otp',status='unknown')
    do nsite=1,NN
       write(10,*) site(nsite,1),site(nsite,2)
    enddo
    close(10)
    call RANDOM_SEED
    call ISLANDS_PARAMETERS

!   do i=1,nisl
!     do phi=0,2*pi,0.01 
!       rad=Radius(i,phi)
!!     write(50+i,*) phi,rad
!       xi=xc(i)+rad*sin(phi+phi0(0,i))
!       eta=yc(i)+rad*cos(phi+phi0(0,i))
!       write(50+i,*) xi,eta
!     enddo
!   enddo

   call COLOR_MAP

    stop
  END PROGRAM Main
!=====================================================================
   SUBROUTINE IRANDOM(n1,n2,n)
!  random integer n between n1 and n2 (n1<n2)
   implicit none
   integer n,n1,n2,ntime
   real x
   call RANDOM_NUMBER(x)
   n=ifix(n1+x*(n2-n1))
!  write(*,*) x,n
   return
   END SUBROUTINE IRANDOM
!=====================================================================
   SUBROUTINE ISLANDS_PARAMETERS
!  parameters for islands with random centre and contour
   use Mod_Map
   implicit none
   integer i,n,ns,nc,k,ierr
   real z 
   real*8 rad,phi,x,y

   amp=0; frq=0

   do i=1,nisl  ! i cycle 

!  island centre
   call IRANDOM(1,NN,n)
   xc(i)=site(n,1); yc(i)=site(n,2)
   call RANDOM_NUMBER(z); xc(i)=xc(i)+z*dx
   call RANDOM_NUMBER(z); yc(i)=yc(i)+z*dy
!  write(*,*) xc,yc

!  rotation angle
   call RANDOM_NUMBER(z); phi0(0,i)=z*2*pi
!  term indep on angle
   call RANDOM_NUMBER(z); amp(0,i)=avamp+z*viamp
   
!  number of sin and cos used
   call IRANDOM(1,nmax,ns) 
!  call IRANDOM(1,nmax,nc)

!  sin amplitudes, frq, and phases
   do k=1,ns
     call RANDOM_NUMBER(z); amp(k,i)=avamp+z*viamp
     call RANDOM_NUMBER(z); frq(k,i)=avfrq+z*vifrq
     call RANDOM_NUMBER(z); phi0(k,i)=z*2*pi 
   enddo

!  cos amplitudes and frq
!  do k=1,nc
!    call RANDOM_NUMBER(z); amp(nmax+k,i)=avamp+z*viamp
!    call RANDOM_NUMBER(z); frq(nmax+k,i)=avfrq+z*vifrq
!  enddo

   enddo ! i cycle

   END SUBROUTINE ISLANDS_PARAMETERS 
!=====================================================================
   FUNCTION Radius(i,phi)
!  radius as function of angle for island i with phi in [0, 2*pi]
   use Mod_Map
   implicit none
   integer i,k
   real*8 phi,ang,rad,Radius

   ang=phi
   if(phi.gt.2*pi) ang=phi-2*pi
   if(phi.lt.0) ang=phi+2*pi

   rad=amp(0,i)
   do k=1,nmax
!  sin terms
     rad=rad+amp(k,i)*(sin(frq(k,i)*k*ang+phi0(k,i)))**2*(1-exp(ang-2*pi))
     rad=rad+amp(k,i)*(sin(phi0(k,i)))**2*exp(ang-2*pi) ! continuity at 0, 2*pi
   enddo
!  cos terms
!  do k=nmax+1,2*nmax
!    rad=rad+amp(k,i)*ang*(cos(frq(k,i)*k*ang))**2*(1-exp(ang-2*pi))
!    rad=rad+amp(k,i)*ang*(cos(frq(k,i)*ang))**2*(1-exp(ang-2*pi))
!    rad=rad+amp(k,i)*(cos(frq(k,i)*ang))**2
!  enddo

   Radius=rad

   return
   END FUNCTION Radius
!=====================================================================
   SUBROUTINE COLOR_MAP
!  mark lattice points inside islands with color parameter
   use Mod_Map
   implicit none
   integer i,k
   real*8 x,y,xi,eta,dist,phi,xx,yy
   real*8 rad,Radius

   open(12,file='map.otp',status='unknown')

   do k=1,NN

     do i=1,nisl
       x=site(k,1); y=site(k,2)
       xx=x-xc(i); yy=y-yc(i)
       if(xx.eq.0.and.yy.gt.0) phi=pi/2
       if(xx.eq.0.and.yy.lt.0) phi=3*pi/2
       if(xx.ne.0) phi=atan(yy/xx)
       if(xx.lt.0.and.yy.gt.0) phi=phi+pi
       if(xx.lt.0.and.yy.lt.0) phi=phi+pi
       if(xx.gt.0.and.yy.lt.0) phi=phi+2*pi
       if(yy.eq.0) phi=0
       dist=sqrt(xx**2+yy**2)
       rad=Radius(i,phi+phi0(0,i))
       if(dist.lt.rad) site(k,3)=1
     enddo
     write(12,*) site(k,1),site(k,2),site(k,3)
   enddo

   close(12)

!   do i=1,nisl
!     write(*,*) xc(i), yc(i)
!   enddo

   END SUBROUTINE COLOR_MAP
!=====================================================================
