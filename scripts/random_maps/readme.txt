The program RandSurf_v01.f90 creates a number of islands with random positions and shapes 
on a uniform backgound 

Each island boder is defined in polar coordinates arount a point (xc,yc)
The border is a function r(phi) = const + a sum of terms (sin(~phi+ phase))**2 
See notes.pdf

Control parameters

   real*8, parameter  :: avamp=20, viamp=2 ! av. amplitude and var. of sinuses
   real*8, parameter  :: avfrq=0.5, vifrq=0.5 ! av. frequency factor and var.
   integer, parameter :: nmax=4  ! max number of sin terms (or cos)
   integer, parameter :: nisl=80 ! number of islands


Output file map.otp: each point on a lattice and a two level function with values 0 or 1
x, y, value

map.gnu is a gnuplot script to convert the output data in an image


examples in the folders map.??
