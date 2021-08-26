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
  double precision                                        :: laser_energy, wavelength, intensity

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
        print '(a)', 'Vacuum: Using Gaussian emission model'

      case DEFAULT
        print '(a)', 'Vacuum: ERROR UNKNOWN WORK FUNCTION TYPE'
        print *, WORK_TYPE
        stop
      END SELECT

      close(unit=ud_laser)
    end subroutine Read_laser_function

    !double precision function laser_parameters(laser_energy, wavelength, intensity)
    !end function laser_parameters
  end Module mod_laser