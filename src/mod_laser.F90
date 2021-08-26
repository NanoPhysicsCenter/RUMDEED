!-------------------------------------------!
! Module for laser input                    !
!                                           !
! Hákon Örn Árnason                         !
! 20.08.21                                  !
!-------------------------------------------!

Module mod_laser
  use mod_global
  use mod_velocity
  use mod_photo_emission

  implicit none
  integer                                                 :: LASER_TYPE

  ! Photon energy parameters
  double precision, allocatable, dimension(:, :)          :: w_theta_arr ! 1:y_num, 1:x_num
  double precision                                        :: laser_energy, laser_variation, intensity

  ! Gaussians
  integer                                     :: num_gauss ! Number of Gaussian points
  double precision                            :: w_theta_base
  double precision, allocatable, dimension(:) :: w_gaussians_A     ! Amplitude
  double precision, allocatable, dimension(:) :: w_gaussians_x     ! x - center
  double precision, allocatable, dimension(:) :: w_gaussians_y     ! y - center
  double precision, allocatable, dimension(:) :: w_gaussians_std_x ! standard deviation / spread in x
  double precision, allocatable, dimension(:) :: w_gaussians_std_y ! standard deviation / spread in y

  ! Type of laser input
  integer, parameter :: LASER_GAUSS      = 1
  integer, parameter :: LASER_SQUARE     = 2  

  ! Type of photon velocity
  integer, parameter :: PHOTON_ZERO  = 1
  integer, parameter :: PHOTON_MB    = 2
  
  PRIVATE
  PUBLIC Init_Photon_Velocity, PHOTON_ZERO, PHOTON_MB
    

  interface
    double precision function Laser_fun(pos, emit, sec)
      double precision, dimension(1:3), intent(in) :: pos
      integer, intent(in)                          :: emit
      integer, intent(out), optional               :: sec
    end function Laser_fun
  end interface
  procedure(Laser_fun), pointer :: ptr_Laser_fun => null()

contains

  subroutine Read_laser_function()
      integer :: ud_laser, IFAIL, i, j
      character(256) :: iomsg

      ! Open the file that contains information about the work function
      open(newunit=ud_laser, iostat=IFAIL, iomsg=iomsg, file='laser', &
        & status='OLD', form='FORMATTED', access='SEQUENTIAL', action='READ')
      if (IFAIL /= 0) then
        print *, 'Vacuum: Failed to open file laser. ABORTING'
        print *, IFAIL
        print *, iomsg
        stop
      end if

      ! Read the type of laser function to use
      read(unit=ud_laser, FMT=*) LASER_TYPE

      SELECT CASE (LASER_TYPE)
      case (LASER_GAUSS)
        ! Gaussian Pulse
        print '(a)', 'Vacuum: Using Gaussian pulse model'
      case (LASER_SQUARE)
        ! Square Pulse
        print '(a)', 'Vacuum: Using Square pulse model'
      case DEFAULT
        print '(a)', 'Vacuum: ERROR UNKNOWN WORK FUNCTION TYPE'
        print *, LASER_TYPE
        stop
      END SELECT

      close(unit=ud_laser)
    end subroutine Read_laser_function

    !double precision function laser_parameters(laser_energy, wavelength, intensity)
    !end function laser_parameters

  
    subroutine Init_Photon_Velocity(PHOTON_MODE)
        integer, intent(in) :: PHOTON_MODE
        SELECT CASE (PHOTON_MODE)
        case(PHOTON_ZERO)
            ptr_Get_Photon_Velocity => Get_Zero_Photon_Velocity
            print '(a)', 'Vacuum: Using zero inital velocity'
        case(PHOTON_MB)
            ptr_Get_Photon_Velocity => Get_Photon_Energy
            print '(a)', 'Vacuum: Using Maxwell-Boltzman energy distribution for Photons'
        case DEFAULT
            print '(a)', 'Vacuum: ERROR UNKNOWN'
            print *, PHOTON_MODE
            stop
        END SELECT
    end subroutine Init_Photon_Velocity
    
    ! ----------------------------------------------------------------------------
    ! Just give zero initial velocity
    function Get_Zero_Photon_Velocity()
        double precision, dimension(1:3) :: Get_Zero_Photon_Velocity
    
        Get_Zero_Photon_Velocity = 0.0d0
    end function Get_Zero_Photon_Velocity

    ! ----------------------------------------------------------------------------
    ! Generate velocity from a Maxwell-Boltzmann distribution.
    ! https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution#Distribution_for_the_velocity_vector
    ! 
    ! Maxwell-Boltzmann velocity distribution is basically a normal distribution for each velocity component with
    ! zero mean and standard deviation \sqrt{k_b T / m}.
    function Get_Photon_Energy()
        double precision, dimension(1:3) :: Get_Photon_Energy
        double precision, dimension(1:2) :: std
        double precision, dimension(1:2) :: mean
        mean = laser_energy
        std = laser_variation ! Standard deviation of the Maxwell-Boltzmann distribution
        !mean = 4.7d0
        !std = 0.1d0 ! Standard deviation of the Maxwell-Boltzmann distribution
    
        ! Get normal distributed numbers.
        ! The Box Muller method gives two numbers.
        ! We overwrite the second element in the array.
        Get_Photon_Energy(1:2) = box_muller(mean, std)
        Get_Photon_Energy(2:3) = box_muller(mean, std)
    
        Get_Photon_Energy(3) = abs(Get_Photon_Energy(3)) ! Positive velocity in the z-direction
    end function Get_Photon_Energy
    
  end Module mod_laser