!-------------------------------------------!
! Module for global variables and constants !
! Kristinn Torfason                         !
! 07.05.15                                  !
!-------------------------------------------!

module mod_global
  implicit none

  ! ----------------------------------------------------------------------------
  ! Define material parameters
  ! GaAs
  ! http://www.iue.tuwien.ac.at/phd/palankovski/node32.html
  ! http://en.wikipedia.org/wiki/Effective_mass_%28solid-state_physics%29#Density_of_states_effective_masses_.28lightly_doped_semiconductors.29
  !character(len=*), parameter :: material = 'GaAs' ! Material used
  !double precision, parameter :: epsilon_r = 13.1d0 ! Relative permittivity
  !double precision, parameter :: m_eeff = 0.067d0 ! Effective mass for electrons
  !double precision, parameter :: m_heff = 0.45d0 ! Effective mass for holes
  !double precision, parameter :: lambda_pos = 2.0d0

  ! ! ----------------------------------------------------------------------------
  ! ! Define material parameters
  ! ! Si at 300 K
  ! ! http://en.wikipedia.org/wiki/Relative_permittivity
  ! ! http://www.semiconductors.co.uk/propiviv5431.htm
  ! ! http://en.wikipedia.org/wiki/Effective_mass_%28solid-state_physics%29#Density_of_states_effective_masses_.28lightly_doped_semiconductors.29
  ! character(len=*), parameter :: material = 'Silicon' ! Material used
  ! !double precision, parameter :: epsilon_r = 11.68d0 ! Relative permittivity
  ! double precision, parameter :: epsilon_r = 11.9d0 ! Relative permittivity
  ! !double precision, parameter :: m_eeff = 1.09d0 ! Effective mass for electrons
  ! !double precision, parameter :: m_heff = 1.15d0 ! Effective mass for holes
  ! !double precision, parameter :: m_eeff = 0.9163d0 ! Effective mass for electrons
  ! !double precision, parameter :: m_heff = 0.49d0 ! Effective mass for holes
  ! double precision, parameter :: m_eeff = 0.34d0 ! Scaps
  ! double precision, parameter :: m_heff = 0.60d0 ! Scaps
  ! !double precision, parameter :: lambda_pos = 0.0746d0 ! Average number of electron/holes pairs to generate per time step
  ! double precision            :: lambda_pos ! Average number of electron/holes pairs to generate per time step
  ! !double precision, parameter :: alpha = 1.11E6 ! Absorption coefficient in m^{-1}
  ! double precision, parameter :: E_g = 1.12d0 ! Band gap in eV
  ! double precision, parameter :: mu_e = 0.14d0 ! Electron mobility
  ! double precision, parameter :: mu_h = 0.045d0 ! Hole mobility

  ! ----------------------------------------------------------------------------
  ! Define material parameters
  ! Methylammonium
  !
  character(len=*), parameter :: material = 'Vacuum' ! Material used

  ! See "Molecular Motion and Dynamic Crystal Structures of Hybrid Halide Perovskites"
  ! it says epsilon_r = 33
  double precision, parameter :: epsilon_r = 1.0d0! Relative permittivity

  ! See "Direct measurement of the exciton binding energy and effective masses for
  !      charge carriers in organic–inorganic tri-halide perovskites"
  !   DOI: 10.1038/NPHYS3357
  ! it says m_eff = 0.104*m_e
  ! m_eff = 0.23d0
  ! m_heff = 0.29d0
  !double precision, parameter :: m_eeff = 0.23d0
  double precision, parameter :: m_eeff = 1.0d0
  double precision, parameter :: m_heff = 1.0d0

  double precision            :: lambda_pos ! Average number of electron/holes pairs to generate per time step
  double precision, parameter :: E_g = 1.55d0 ! Band gap in eV

  ! See "Impact of work function of back contact of perovskite solar cells without hole transport material analyzed by device simulation"
  double precision, parameter :: mu_e = 2.0d0*1.0d-4 ! [cm²/Vs] -> [m²/Vs] Electron mobility
  double precision, parameter :: mu_h = 2.0d0*1.0d-4 ! [cm²/Vs] -> [m²/Vs] Hole mobility


  ! ----------------------------------------------------------------------------
  ! Define physical constants
  ! See http://physics.nist.gov/cuu/Constants/index.html
  double precision, parameter :: pi = 3.14159265358979324d0 ! Pi
  double precision, parameter :: pi_2 = 2.0d0*pi ! 2*pi
  double precision, parameter :: h = 6.626070040d-34 ! Planck's constant (Js)
  double precision, parameter :: k_b = 1.38064852d-23 ! Boltzmann constant (J/K)
  double precision, parameter :: c = 299792458.0d0 ! Speed of light (m/s)
  double precision, parameter :: mu_0 = 4.0d0*pi * 1d-7 ! Vacuum permeability (H/m)
  double precision, parameter :: epsilon_0 = 1.0d0/(mu_0 * c**2) ! Vacuum permittivity (F/m)
  !double precision, parameter :: epsilon_0 = 8.854187817d-12 ! Farads / meters
  double precision, parameter :: epsilon = epsilon_r * epsilon_0 ! Permittivity
  double precision, parameter :: m_u = 1.660539040d-27 ! Atomic mass unit (kg)


  double precision, parameter :: m_0 = 9.10938356d-31 ! Free electron mass (kg)
  !double precision, parameter :: m_e = m_eeff * m_0 ! m_e* Effective electron mass (kg)
  !double precision, parameter :: m_h = m_heff * m_0! m_h* Effective hole mass (kg)

  double precision, parameter :: q_0 = 1.6021766208d-19 ! Standard charge (C)
  double precision, parameter :: q_02 = q_0**2 ! Standard charge squared (C)
  !double precision, parameter :: q_e = -1.602176565d-19 ! Electron charge (Coulomb)
  !double precision, parameter :: q_h = 1.602176565d-19 ! Hole charge (Coulomb)

  double precision, parameter :: hc_nm = h*c/(1.0d-9)
  double precision, parameter :: hc_evnm = h*c/(q_0*1.0d-9) ! E_\lambda = h*c/lambda

  double precision, parameter :: T_room = 300.0d0 ! Room temperature in Kelvin


  ! ----------------------------------------------------------------------------
  ! Define scales used when reading and writing data
  double precision, parameter :: length_scale = 1.0d-9 ! Length scale (1 nanometer)
  double precision, parameter :: time_scale = 1.0d-12 ! Time scale (1 ps)
  double precision, parameter :: vel_scale = length_scale / time_scale ! Velocity scale (1 nm / 1 ps)
  double precision            :: cur_scale  ! Current scale (1 mA/cm^2)


  ! ----------------------------------------------------------------------------
  ! Define constans
  integer, parameter :: MAX_PARTICLES = 50000 ! Maximum number of particles of each type allowed in the system
  integer, parameter :: MAX_EMITTERS  = 10    ! Maximum number of emitters in the system


  !! ----------------------------------------------------------------------------
  !! Define the charge of the electrons and holes
  !integer, parameter :: charge_elec = -1
  !integer, parameter :: charge_hole = +1

  ! ----------------------------------------------------------------------------
  ! Define the particle species
  integer, parameter :: species_unkown   = 0 ! Unknown particle
  integer, parameter :: species_elec     = 1 ! Electron
  integer, parameter :: species_hole     = 2 ! Hole
  integer, parameter :: nrSpecies        = 2 ! 2 = Elec, Hole


  ! ----------------------------------------------------------------------------
  ! Particle removal flags
  integer, parameter :: remove_unknown = 0
  integer, parameter :: remove_top    = 1
  integer, parameter :: remove_bot   = 2

  ! ----------------------------------------------------------------------------
  ! Define storage arrays for particles
  ! Fyrst dimension is x,y,z, second one is the number of the particle
  double precision, dimension(:, :), allocatable :: particles_cur_pos    ! Current position (1:3, 1:MAX_PARTICLES)
  double precision, dimension(:, :), allocatable :: particles_prev_pos   ! Previous position
  double precision, dimension(:, :), allocatable :: particles_cur_vel    ! Current velocity
  double precision, dimension(:, :), allocatable :: particles_cur_accel  ! Current acceleration
  double precision, dimension(:, :), allocatable :: particles_prev_accel ! Previous acceleration
  double precision, dimension(:)   , allocatable :: particles_charge     ! Charge
  integer         , dimension(:)   , allocatable :: particles_species    ! Type of particle
  double precision, dimension(:)   , allocatable :: particles_mass       ! Mass
  integer         , dimension(:)   , allocatable :: particles_step       ! Time step when particle was created
  logical         , dimension(:)   , allocatable :: particles_mask       ! Mask array used to indicate which particles should be removed

  ! ----------------------------------------------------------------------------
  ! Define storage arrays for emitters
  ! Fyrst dimension is x,y,z, second one is the number of the emitter
  double precision, dimension(:, :), allocatable :: emitters_pos         ! Position of the emitters (1:3, 1:MAX_EMITTERS)
  double precision, dimension(:, :), allocatable :: emitters_dim         ! Dimensions of the emitters
  integer,          dimension(:),    allocatable :: emitters_Type        ! The type of emitter
  integer,          dimension(:),    allocatable :: emitters_delay       ! The time step the emitters become active

  ! Density map
  integer, dimension(:, :), allocatable :: density_map_elec
  integer, dimension(:, :), allocatable :: density_map_hole
  integer, parameter                    :: N_x_densmap = 100, N_y_densmap = 100
  double precision                      :: dens_x_d, dens_y_d


  ! ----------------------------------------------------------------------------
  ! Define input parameters
  double precision :: V_d     ! Voltage over the gap
  double precision :: V_s     ! Voltage from the source
  double precision :: d       ! Gap spacing
  double precision :: E_z     ! Electric field in the y-direction (E_z = -V/d)
  double precision :: E_zunit ! Unit electric field (E_zunit = -1/d) (See Ramo Current)

  double precision, dimension(1:3) :: box_dim ! Dimensions of the cell

  double precision :: time_step  ! Size of the time_step
  double precision :: time_step2 ! time_step squared

  integer          :: steps      ! Number of time steps in the simulation


  ! ----------------------------------------------------------------------------
  ! Define run time variables
  integer :: nrPart ! Number of particles in the system (nrPart = nrElec + nrHole + nrFixedPart)
  integer :: nrElec ! Number of electrons in the system
  integer :: nrHole ! Number of holes in the system
  integer :: nrElecHole
  integer :: nrEmit ! Number of emitters in the system

  integer :: nrPart_remove_top
  integer :: nrPart_remove_bot

  integer :: nrElec_remove_top
  integer :: nrElec_remove_bot

  integer :: nrHole_remove_top
  integer :: nrHole_remove_bot

  integer :: nrPart_remove ! Number of particles to be removed
  integer :: nrElec_remove ! Number of electrons to be removed
  integer :: nrHole_remove ! Number of holes to be removed

  integer :: startElecHoles
  integer :: endElecHoles


  ! ----------------------------------------------------------------------------
  ! Emitter types
  integer, parameter :: EMIT_UNKNOWN   = 0
  integer, parameter :: EMIT_CIRCLE    = 1
  integer, parameter :: EMIT_RECTANGLE = 2


  double precision, dimension(:), allocatable :: ramo_current
  double precision :: ramo_cur_prev
  double precision :: ramo_integral

  double precision :: cur_time ! This is updated in the main loop, given in ps
  integer, parameter :: MAX_LIFE_TIME = 1000
  integer, dimension(:, :), allocatable :: life_time



  ! ----------------------------------------------------------------------------
  ! Define run time constants
  double precision, parameter :: div_fac_c = 1.0d-7 * c**2 / epsilon_r ! 1/(4*pi*epsilon_0*epsilon_r)


  ! ----------------------------------------------------------------------------
  ! Parameters for random number generators
  ! http://en.wikipedia.org/wiki/Mersenne_Twister#SFMT
  integer                                            :: SEED = 2815


  ! ----------------------------------------------------------------------------
  ! unit descriptors for data files
  integer :: ud_pos ! Position file
  integer :: ud_emit ! File for emitted electrons and holes
  integer :: ud_absorb ! File for absorbed electrons and holes
  integer :: ud_absorb_top ! File for absorbed electrons and holes
  integer :: ud_absorb_bot ! File for absorbed electrons and holes
  integer :: ud_ramo ! File for the Ramo current
  integer :: ud_volt ! Voltage in the system
  integer :: ud_debug ! File for debuging and testing
  integer :: ud_dipole_pos ! File for dipole positions
  integer :: ud_dipole_vec ! File for dipole orientation
  integer :: ud_field ! File for longitudinal field
  integer :: ud_density_map_elec ! Density maps for electrons
  integer :: ud_density_map_hole ! Density maps for holes
  integer :: ud_density_map_total ! Density maps for holes - electrons

  !--
  ! Used in Calc_Field_at in the Verlet module. Needs to be shared between threads.
  ! Can we move it into the function? ToDo...
  double precision, dimension(1:3) :: force_tot

  !-----------------------------------------------------------------------------
  ! Nodal Analysis
  double precision, dimension(1:5) :: V_cur, V_prev ! Voltage and branch currents for the nodal analysis

  ! ----------------------------------------------------------------------------
  ! Define namelist
  namelist /input/ V_s, box_dim, time_step, steps, &
                   nrEmit, emitters_pos, emitters_dim, &
                   emitters_type, emitters_delay

  ! ----------------------------------------------------------------------------
  ! Prodecure interfaces and pointers
  ! These are subroutines/functions that change depending on the type of
  ! emission / geometry used.
  interface
    subroutine Check_Boundary(i)
      integer, intent(in) :: i
    end subroutine Check_Boundary

    function Electric_Field(pos)
      double precision, dimension(1:3), intent(in) :: pos
      double precision, dimension(1:3)             :: Electric_Field
    end function Electric_Field

    subroutine Do_Emission(step)
      integer, intent(in) :: step
    end subroutine Do_Emission
  end interface

  procedure(Check_Boundary), pointer :: ptr_Check_Boundary => null()
  procedure(Electric_Field), pointer :: ptr_field_E => null()
  procedure(Do_Emission), pointer    :: ptr_Do_Emission => null()
contains

  logical function isinf(a)
    double precision, intent(in) :: a

    ! if the number is infinity then the results is always infinity
    if ((a-1.0d0) == a) then
      isinf = .true.
    else
      isinf = .false.
    end if
  end

! PGI compiler does not have is isnan function
#if defined(__PGI)
  logical function isnan(x)
    use ieee_arithmetic
    double precision, intent(in) :: x

    isnan = ieee_is_nan(x)
  end function isnan
#endif

! So far the PGI compiler (v. 17.4) has not implemented the NORM2 function
! from the Fortran 2008 Standard
#if defined(__PGI)
  double precision function norm2(a)
    double precision, dimension(:), intent(in) :: a
    ! integer                                    :: i

    norm2 = sqrt(sum(a**2))

    ! norm2 = 0.0d0
    ! do i = lbound(a, dim=1), ubound(a, dim=1)
    !   norm2 = norm2 + a(i)**2
    ! end do
    ! norm2 = sqrt(norm2)
  end function norm2
#endif

!***********************************************************************************************************************************
!  M55INV  -  Compute the inverse of a 5x5 matrix.
!
!  A       = input 5x5 matrix to be inverted
!  AINV    = output 5x5 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
! http://www.davidgsimpson.com/software.html
!***********************************************************************************************************************************

      SUBROUTINE M55INV (A, AINV, OK_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(5,5), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(5,5), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET, A11, A12, A13, A14, A15, A21, A22, A23, A24, &
         A25, A31, A32, A33, A34, A35, A41, A42, A43, A44, A45,   &
         A51, A52, A53, A54, A55
      DOUBLE PRECISION, DIMENSION(5,5) :: COFACTOR


      A11=A(1,1); A12=A(1,2); A13=A(1,3); A14=A(1,4); A15=A(1,5)
      A21=A(2,1); A22=A(2,2); A23=A(2,3); A24=A(2,4); A25=A(2,5)
      A31=A(3,1); A32=A(3,2); A33=A(3,3); A34=A(3,4); A35=A(3,5)
      A41=A(4,1); A42=A(4,2); A43=A(4,3); A44=A(4,4); A45=A(4,5)
      A51=A(5,1); A52=A(5,2); A53=A(5,3); A54=A(5,4); A55=A(5,5)

      DET = A15*A24*A33*A42*A51-A14*A25*A33*A42*A51-A15*A23*A34*A42*A51+    &
         A13*A25*A34*A42*A51+A14*A23*A35*A42*A51-A13*A24*A35*A42*A51-       &
         A15*A24*A32*A43*A51+A14*A25*A32*A43*A51+A15*A22*A34*A43*A51-       &
         A12*A25*A34*A43*A51-A14*A22*A35*A43*A51+A12*A24*A35*A43*A51+       &
         A15*A23*A32*A44*A51-A13*A25*A32*A44*A51-A15*A22*A33*A44*A51+       &
         A12*A25*A33*A44*A51+A13*A22*A35*A44*A51-A12*A23*A35*A44*A51-       &
         A14*A23*A32*A45*A51+A13*A24*A32*A45*A51+A14*A22*A33*A45*A51-       &
         A12*A24*A33*A45*A51-A13*A22*A34*A45*A51+A12*A23*A34*A45*A51-       &
         A15*A24*A33*A41*A52+A14*A25*A33*A41*A52+A15*A23*A34*A41*A52-       &
         A13*A25*A34*A41*A52-A14*A23*A35*A41*A52+A13*A24*A35*A41*A52+       &
         A15*A24*A31*A43*A52-A14*A25*A31*A43*A52-A15*A21*A34*A43*A52+       &
         A11*A25*A34*A43*A52+A14*A21*A35*A43*A52-A11*A24*A35*A43*A52-       &
         A15*A23*A31*A44*A52+A13*A25*A31*A44*A52+A15*A21*A33*A44*A52-       &
         A11*A25*A33*A44*A52-A13*A21*A35*A44*A52+A11*A23*A35*A44*A52+       &
         A14*A23*A31*A45*A52-A13*A24*A31*A45*A52-A14*A21*A33*A45*A52+       &
         A11*A24*A33*A45*A52+A13*A21*A34*A45*A52-A11*A23*A34*A45*A52+       &
         A15*A24*A32*A41*A53-A14*A25*A32*A41*A53-A15*A22*A34*A41*A53+       &
         A12*A25*A34*A41*A53+A14*A22*A35*A41*A53-A12*A24*A35*A41*A53-       &
         A15*A24*A31*A42*A53+A14*A25*A31*A42*A53+A15*A21*A34*A42*A53-       &
         A11*A25*A34*A42*A53-A14*A21*A35*A42*A53+A11*A24*A35*A42*A53+       &
         A15*A22*A31*A44*A53-A12*A25*A31*A44*A53-A15*A21*A32*A44*A53+       &
         A11*A25*A32*A44*A53+A12*A21*A35*A44*A53-A11*A22*A35*A44*A53-       &
         A14*A22*A31*A45*A53+A12*A24*A31*A45*A53+A14*A21*A32*A45*A53-       &
         A11*A24*A32*A45*A53-A12*A21*A34*A45*A53+A11*A22*A34*A45*A53-       &
         A15*A23*A32*A41*A54+A13*A25*A32*A41*A54+A15*A22*A33*A41*A54-       &
         A12*A25*A33*A41*A54-A13*A22*A35*A41*A54+A12*A23*A35*A41*A54+       &
         A15*A23*A31*A42*A54-A13*A25*A31*A42*A54-A15*A21*A33*A42*A54+       &
         A11*A25*A33*A42*A54+A13*A21*A35*A42*A54-A11*A23*A35*A42*A54-       &
         A15*A22*A31*A43*A54+A12*A25*A31*A43*A54+A15*A21*A32*A43*A54-       &
         A11*A25*A32*A43*A54-A12*A21*A35*A43*A54+A11*A22*A35*A43*A54+       &
         A13*A22*A31*A45*A54-A12*A23*A31*A45*A54-A13*A21*A32*A45*A54+       &
         A11*A23*A32*A45*A54+A12*A21*A33*A45*A54-A11*A22*A33*A45*A54+       &
         A14*A23*A32*A41*A55-A13*A24*A32*A41*A55-A14*A22*A33*A41*A55+       &
         A12*A24*A33*A41*A55+A13*A22*A34*A41*A55-A12*A23*A34*A41*A55-       &
         A14*A23*A31*A42*A55+A13*A24*A31*A42*A55+A14*A21*A33*A42*A55-       &
         A11*A24*A33*A42*A55-A13*A21*A34*A42*A55+A11*A23*A34*A42*A55+       &
         A14*A22*A31*A43*A55-A12*A24*A31*A43*A55-A14*A21*A32*A43*A55+       &
         A11*A24*A32*A43*A55+A12*A21*A34*A43*A55-A11*A22*A34*A43*A55-       &
         A13*A22*A31*A44*A55+A12*A23*A31*A44*A55+A13*A21*A32*A44*A55-       &
         A11*A23*A32*A44*A55-A12*A21*A33*A44*A55+A11*A22*A33*A44*A55

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = A25*A34*A43*A52-A24*A35*A43*A52-A25*A33*A44*A52+      &
         A23*A35*A44*A52+A24*A33*A45*A52-A23*A34*A45*A52-A25*A34*A42*A53+   &
         A24*A35*A42*A53+A25*A32*A44*A53-A22*A35*A44*A53-A24*A32*A45*A53+   &
         A22*A34*A45*A53+A25*A33*A42*A54-A23*A35*A42*A54-A25*A32*A43*A54+   &
         A22*A35*A43*A54+A23*A32*A45*A54-A22*A33*A45*A54-A24*A33*A42*A55+   &
         A23*A34*A42*A55+A24*A32*A43*A55-A22*A34*A43*A55-A23*A32*A44*A55+   &
         A22*A33*A44*A55

      COFACTOR(2,1) = -A15*A34*A43*A52+A14*A35*A43*A52+A15*A33*A44*A52-     &
         A13*A35*A44*A52-A14*A33*A45*A52+A13*A34*A45*A52+A15*A34*A42*A53-   &
         A14*A35*A42*A53-A15*A32*A44*A53+A12*A35*A44*A53+A14*A32*A45*A53-   &
         A12*A34*A45*A53-A15*A33*A42*A54+A13*A35*A42*A54+A15*A32*A43*A54-   &
         A12*A35*A43*A54-A13*A32*A45*A54+A12*A33*A45*A54+A14*A33*A42*A55-   &
         A13*A34*A42*A55-A14*A32*A43*A55+A12*A34*A43*A55+A13*A32*A44*A55-   &
         A12*A33*A44*A55

      COFACTOR(3,1) = A15*A24*A43*A52-A14*A25*A43*A52-A15*A23*A44*A52+      &
         A13*A25*A44*A52+A14*A23*A45*A52-A13*A24*A45*A52-A15*A24*A42*A53+   &
         A14*A25*A42*A53+A15*A22*A44*A53-A12*A25*A44*A53-A14*A22*A45*A53+   &
         A12*A24*A45*A53+A15*A23*A42*A54-A13*A25*A42*A54-A15*A22*A43*A54+   &
         A12*A25*A43*A54+A13*A22*A45*A54-A12*A23*A45*A54-A14*A23*A42*A55+   &
         A13*A24*A42*A55+A14*A22*A43*A55-A12*A24*A43*A55-A13*A22*A44*A55+   &
         A12*A23*A44*A55

      COFACTOR(4,1) = -A15*A24*A33*A52+A14*A25*A33*A52+A15*A23*A34*A52-     &
         A13*A25*A34*A52-A14*A23*A35*A52+A13*A24*A35*A52+A15*A24*A32*A53-   &
         A14*A25*A32*A53-A15*A22*A34*A53+A12*A25*A34*A53+A14*A22*A35*A53-   &
         A12*A24*A35*A53-A15*A23*A32*A54+A13*A25*A32*A54+A15*A22*A33*A54-   &
         A12*A25*A33*A54-A13*A22*A35*A54+A12*A23*A35*A54+A14*A23*A32*A55-   &
         A13*A24*A32*A55-A14*A22*A33*A55+A12*A24*A33*A55+A13*A22*A34*A55-   &
         A12*A23*A34*A55

      COFACTOR(5,1) = A15*A24*A33*A42-A14*A25*A33*A42-A15*A23*A34*A42+      &
         A13*A25*A34*A42+A14*A23*A35*A42-A13*A24*A35*A42-A15*A24*A32*A43+   &
         A14*A25*A32*A43+A15*A22*A34*A43-A12*A25*A34*A43-A14*A22*A35*A43+   &
         A12*A24*A35*A43+A15*A23*A32*A44-A13*A25*A32*A44-A15*A22*A33*A44+   &
         A12*A25*A33*A44+A13*A22*A35*A44-A12*A23*A35*A44-A14*A23*A32*A45+   &
         A13*A24*A32*A45+A14*A22*A33*A45-A12*A24*A33*A45-A13*A22*A34*A45+   &
         A12*A23*A34*A45

      COFACTOR(1,2) = -A25*A34*A43*A51+A24*A35*A43*A51+A25*A33*A44*A51-     &
         A23*A35*A44*A51-A24*A33*A45*A51+A23*A34*A45*A51+A25*A34*A41*A53-   &
         A24*A35*A41*A53-A25*A31*A44*A53+A21*A35*A44*A53+A24*A31*A45*A53-   &
         A21*A34*A45*A53-A25*A33*A41*A54+A23*A35*A41*A54+A25*A31*A43*A54-   &
         A21*A35*A43*A54-A23*A31*A45*A54+A21*A33*A45*A54+A24*A33*A41*A55-   &
         A23*A34*A41*A55-A24*A31*A43*A55+A21*A34*A43*A55+A23*A31*A44*A55-   &
         A21*A33*A44*A55

      COFACTOR(2,2) = A15*A34*A43*A51-A14*A35*A43*A51-A15*A33*A44*A51+      &
         A13*A35*A44*A51+A14*A33*A45*A51-A13*A34*A45*A51-A15*A34*A41*A53+   &
         A14*A35*A41*A53+A15*A31*A44*A53-A11*A35*A44*A53-A14*A31*A45*A53+   &
         A11*A34*A45*A53+A15*A33*A41*A54-A13*A35*A41*A54-A15*A31*A43*A54+   &
         A11*A35*A43*A54+A13*A31*A45*A54-A11*A33*A45*A54-A14*A33*A41*A55+   &
         A13*A34*A41*A55+A14*A31*A43*A55-A11*A34*A43*A55-A13*A31*A44*A55+   &
         A11*A33*A44*A55

      COFACTOR(3,2) = -A15*A24*A43*A51+A14*A25*A43*A51+A15*A23*A44*A51-     &
         A13*A25*A44*A51-A14*A23*A45*A51+A13*A24*A45*A51+A15*A24*A41*A53-   &
         A14*A25*A41*A53-A15*A21*A44*A53+A11*A25*A44*A53+A14*A21*A45*A53-   &
         A11*A24*A45*A53-A15*A23*A41*A54+A13*A25*A41*A54+A15*A21*A43*A54-   &
         A11*A25*A43*A54-A13*A21*A45*A54+A11*A23*A45*A54+A14*A23*A41*A55-   &
         A13*A24*A41*A55-A14*A21*A43*A55+A11*A24*A43*A55+A13*A21*A44*A55-   &
         A11*A23*A44*A55

      COFACTOR(4,2) = A15*A24*A33*A51-A14*A25*A33*A51-A15*A23*A34*A51+      &
         A13*A25*A34*A51+A14*A23*A35*A51-A13*A24*A35*A51-A15*A24*A31*A53+   &
         A14*A25*A31*A53+A15*A21*A34*A53-A11*A25*A34*A53-A14*A21*A35*A53+   &
         A11*A24*A35*A53+A15*A23*A31*A54-A13*A25*A31*A54-A15*A21*A33*A54+   &
         A11*A25*A33*A54+A13*A21*A35*A54-A11*A23*A35*A54-A14*A23*A31*A55+   &
         A13*A24*A31*A55+A14*A21*A33*A55-A11*A24*A33*A55-A13*A21*A34*A55+   &
         A11*A23*A34*A55

      COFACTOR(5,2) = -A15*A24*A33*A41+A14*A25*A33*A41+A15*A23*A34*A41-     &
         A13*A25*A34*A41-A14*A23*A35*A41+A13*A24*A35*A41+A15*A24*A31*A43-   &
         A14*A25*A31*A43-A15*A21*A34*A43+A11*A25*A34*A43+A14*A21*A35*A43-   &
         A11*A24*A35*A43-A15*A23*A31*A44+A13*A25*A31*A44+A15*A21*A33*A44-   &
         A11*A25*A33*A44-A13*A21*A35*A44+A11*A23*A35*A44+A14*A23*A31*A45-   &
         A13*A24*A31*A45-A14*A21*A33*A45+A11*A24*A33*A45+A13*A21*A34*A45-   &
         A11*A23*A34*A45

      COFACTOR(1,3) = A25*A34*A42*A51-A24*A35*A42*A51-A25*A32*A44*A51+      &
         A22*A35*A44*A51+A24*A32*A45*A51-A22*A34*A45*A51-A25*A34*A41*A52+   &
         A24*A35*A41*A52+A25*A31*A44*A52-A21*A35*A44*A52-A24*A31*A45*A52+   &
         A21*A34*A45*A52+A25*A32*A41*A54-A22*A35*A41*A54-A25*A31*A42*A54+   &
         A21*A35*A42*A54+A22*A31*A45*A54-A21*A32*A45*A54-A24*A32*A41*A55+   &
         A22*A34*A41*A55+A24*A31*A42*A55-A21*A34*A42*A55-A22*A31*A44*A55+   &
         A21*A32*A44*A55

      COFACTOR(2,3) = -A15*A34*A42*A51+A14*A35*A42*A51+A15*A32*A44*A51-     &
         A12*A35*A44*A51-A14*A32*A45*A51+A12*A34*A45*A51+A15*A34*A41*A52-   &
         A14*A35*A41*A52-A15*A31*A44*A52+A11*A35*A44*A52+A14*A31*A45*A52-   &
         A11*A34*A45*A52-A15*A32*A41*A54+A12*A35*A41*A54+A15*A31*A42*A54-   &
         A11*A35*A42*A54-A12*A31*A45*A54+A11*A32*A45*A54+A14*A32*A41*A55-   &
         A12*A34*A41*A55-A14*A31*A42*A55+A11*A34*A42*A55+A12*A31*A44*A55-   &
         A11*A32*A44*A55

      COFACTOR(3,3) = A15*A24*A42*A51-A14*A25*A42*A51-A15*A22*A44*A51+      &
         A12*A25*A44*A51+A14*A22*A45*A51-A12*A24*A45*A51-A15*A24*A41*A52+   &
         A14*A25*A41*A52+A15*A21*A44*A52-A11*A25*A44*A52-A14*A21*A45*A52+   &
         A11*A24*A45*A52+A15*A22*A41*A54-A12*A25*A41*A54-A15*A21*A42*A54+   &
         A11*A25*A42*A54+A12*A21*A45*A54-A11*A22*A45*A54-A14*A22*A41*A55+   &
         A12*A24*A41*A55+A14*A21*A42*A55-A11*A24*A42*A55-A12*A21*A44*A55+   &
         A11*A22*A44*A55

      COFACTOR(4,3) = -A15*A24*A32*A51+A14*A25*A32*A51+A15*A22*A34*A51-     &
         A12*A25*A34*A51-A14*A22*A35*A51+A12*A24*A35*A51+A15*A24*A31*A52-   &
         A14*A25*A31*A52-A15*A21*A34*A52+A11*A25*A34*A52+A14*A21*A35*A52-   &
         A11*A24*A35*A52-A15*A22*A31*A54+A12*A25*A31*A54+A15*A21*A32*A54-   &
         A11*A25*A32*A54-A12*A21*A35*A54+A11*A22*A35*A54+A14*A22*A31*A55-   &
         A12*A24*A31*A55-A14*A21*A32*A55+A11*A24*A32*A55+A12*A21*A34*A55-   &
         A11*A22*A34*A55

      COFACTOR(5,3) = A15*A24*A32*A41-A14*A25*A32*A41-A15*A22*A34*A41+      &
         A12*A25*A34*A41+A14*A22*A35*A41-A12*A24*A35*A41-A15*A24*A31*A42+   &
         A14*A25*A31*A42+A15*A21*A34*A42-A11*A25*A34*A42-A14*A21*A35*A42+   &
         A11*A24*A35*A42+A15*A22*A31*A44-A12*A25*A31*A44-A15*A21*A32*A44+   &
         A11*A25*A32*A44+A12*A21*A35*A44-A11*A22*A35*A44-A14*A22*A31*A45+   &
         A12*A24*A31*A45+A14*A21*A32*A45-A11*A24*A32*A45-A12*A21*A34*A45+   &
         A11*A22*A34*A45

      COFACTOR(1,4) = -A25*A33*A42*A51+A23*A35*A42*A51+A25*A32*A43*A51-     &
         A22*A35*A43*A51-A23*A32*A45*A51+A22*A33*A45*A51+A25*A33*A41*A52-   &
         A23*A35*A41*A52-A25*A31*A43*A52+A21*A35*A43*A52+A23*A31*A45*A52-   &
         A21*A33*A45*A52-A25*A32*A41*A53+A22*A35*A41*A53+A25*A31*A42*A53-   &
         A21*A35*A42*A53-A22*A31*A45*A53+A21*A32*A45*A53+A23*A32*A41*A55-   &
         A22*A33*A41*A55-A23*A31*A42*A55+A21*A33*A42*A55+A22*A31*A43*A55-   &
         A21*A32*A43*A55

      COFACTOR(2,4) = A15*A33*A42*A51-A13*A35*A42*A51-A15*A32*A43*A51+      &
         A12*A35*A43*A51+A13*A32*A45*A51-A12*A33*A45*A51-A15*A33*A41*A52+   &
         A13*A35*A41*A52+A15*A31*A43*A52-A11*A35*A43*A52-A13*A31*A45*A52+   &
         A11*A33*A45*A52+A15*A32*A41*A53-A12*A35*A41*A53-A15*A31*A42*A53+   &
         A11*A35*A42*A53+A12*A31*A45*A53-A11*A32*A45*A53-A13*A32*A41*A55+   &
         A12*A33*A41*A55+A13*A31*A42*A55-A11*A33*A42*A55-A12*A31*A43*A55+   &
         A11*A32*A43*A55

      COFACTOR(3,4) = -A15*A23*A42*A51+A13*A25*A42*A51+A15*A22*A43*A51-     &
         A12*A25*A43*A51-A13*A22*A45*A51+A12*A23*A45*A51+A15*A23*A41*A52-   &
         A13*A25*A41*A52-A15*A21*A43*A52+A11*A25*A43*A52+A13*A21*A45*A52-   &
         A11*A23*A45*A52-A15*A22*A41*A53+A12*A25*A41*A53+A15*A21*A42*A53-   &
         A11*A25*A42*A53-A12*A21*A45*A53+A11*A22*A45*A53+A13*A22*A41*A55-   &
         A12*A23*A41*A55-A13*A21*A42*A55+A11*A23*A42*A55+A12*A21*A43*A55-   &
         A11*A22*A43*A55

      COFACTOR(4,4) = A15*A23*A32*A51-A13*A25*A32*A51-A15*A22*A33*A51+      &
         A12*A25*A33*A51+A13*A22*A35*A51-A12*A23*A35*A51-A15*A23*A31*A52+   &
         A13*A25*A31*A52+A15*A21*A33*A52-A11*A25*A33*A52-A13*A21*A35*A52+   &
         A11*A23*A35*A52+A15*A22*A31*A53-A12*A25*A31*A53-A15*A21*A32*A53+   &
         A11*A25*A32*A53+A12*A21*A35*A53-A11*A22*A35*A53-A13*A22*A31*A55+   &
         A12*A23*A31*A55+A13*A21*A32*A55-A11*A23*A32*A55-A12*A21*A33*A55+   &
         A11*A22*A33*A55

      COFACTOR(5,4) = -A15*A23*A32*A41+A13*A25*A32*A41+A15*A22*A33*A41-     &
         A12*A25*A33*A41-A13*A22*A35*A41+A12*A23*A35*A41+A15*A23*A31*A42-   &
         A13*A25*A31*A42-A15*A21*A33*A42+A11*A25*A33*A42+A13*A21*A35*A42-   &
         A11*A23*A35*A42-A15*A22*A31*A43+A12*A25*A31*A43+A15*A21*A32*A43-   &
         A11*A25*A32*A43-A12*A21*A35*A43+A11*A22*A35*A43+A13*A22*A31*A45-   &
         A12*A23*A31*A45-A13*A21*A32*A45+A11*A23*A32*A45+A12*A21*A33*A45-   &
         A11*A22*A33*A45

      COFACTOR(1,5) = A24*A33*A42*A51-A23*A34*A42*A51-A24*A32*A43*A51+      &
         A22*A34*A43*A51+A23*A32*A44*A51-A22*A33*A44*A51-A24*A33*A41*A52+   &
         A23*A34*A41*A52+A24*A31*A43*A52-A21*A34*A43*A52-A23*A31*A44*A52+   &
         A21*A33*A44*A52+A24*A32*A41*A53-A22*A34*A41*A53-A24*A31*A42*A53+   &
         A21*A34*A42*A53+A22*A31*A44*A53-A21*A32*A44*A53-A23*A32*A41*A54+   &
         A22*A33*A41*A54+A23*A31*A42*A54-A21*A33*A42*A54-A22*A31*A43*A54+   &
         A21*A32*A43*A54

      COFACTOR(2,5) = -A14*A33*A42*A51+A13*A34*A42*A51+A14*A32*A43*A51-     &
         A12*A34*A43*A51-A13*A32*A44*A51+A12*A33*A44*A51+A14*A33*A41*A52-   &
         A13*A34*A41*A52-A14*A31*A43*A52+A11*A34*A43*A52+A13*A31*A44*A52-   &
         A11*A33*A44*A52-A14*A32*A41*A53+A12*A34*A41*A53+A14*A31*A42*A53-   &
         A11*A34*A42*A53-A12*A31*A44*A53+A11*A32*A44*A53+A13*A32*A41*A54-   &
         A12*A33*A41*A54-A13*A31*A42*A54+A11*A33*A42*A54+A12*A31*A43*A54-   &
         A11*A32*A43*A54

      COFACTOR(3,5) = A14*A23*A42*A51-A13*A24*A42*A51-A14*A22*A43*A51+      &
         A12*A24*A43*A51+A13*A22*A44*A51-A12*A23*A44*A51-A14*A23*A41*A52+   &
         A13*A24*A41*A52+A14*A21*A43*A52-A11*A24*A43*A52-A13*A21*A44*A52+   &
         A11*A23*A44*A52+A14*A22*A41*A53-A12*A24*A41*A53-A14*A21*A42*A53+   &
         A11*A24*A42*A53+A12*A21*A44*A53-A11*A22*A44*A53-A13*A22*A41*A54+   &
         A12*A23*A41*A54+A13*A21*A42*A54-A11*A23*A42*A54-A12*A21*A43*A54+   &
         A11*A22*A43*A54

      COFACTOR(4,5) = -A14*A23*A32*A51+A13*A24*A32*A51+A14*A22*A33*A51-     &
         A12*A24*A33*A51-A13*A22*A34*A51+A12*A23*A34*A51+A14*A23*A31*A52-   &
         A13*A24*A31*A52-A14*A21*A33*A52+A11*A24*A33*A52+A13*A21*A34*A52-   &
         A11*A23*A34*A52-A14*A22*A31*A53+A12*A24*A31*A53+A14*A21*A32*A53-   &
         A11*A24*A32*A53-A12*A21*A34*A53+A11*A22*A34*A53+A13*A22*A31*A54-   &
         A12*A23*A31*A54-A13*A21*A32*A54+A11*A23*A32*A54+A12*A21*A33*A54-   &
         A11*A22*A33*A54

      COFACTOR(5,5) = A14*A23*A32*A41-A13*A24*A32*A41-A14*A22*A33*A41+      &
         A12*A24*A33*A41+A13*A22*A34*A41-A12*A23*A34*A41-A14*A23*A31*A42+   &
         A13*A24*A31*A42+A14*A21*A33*A42-A11*A24*A33*A42-A13*A21*A34*A42+   &
         A11*A23*A34*A42+A14*A22*A31*A43-A12*A24*A31*A43-A14*A21*A32*A43+   &
         A11*A24*A32*A43+A12*A21*A34*A43-A11*A22*A34*A43-A13*A22*A31*A44+   &
         A12*A23*A31*A44+A13*A21*A32*A44-A11*A23*A32*A44-A12*A21*A33*A44+   &
         A11*A22*A33*A44

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE M55INV
end module mod_global
