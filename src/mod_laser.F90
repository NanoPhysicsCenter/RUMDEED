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
  integer            :: LASER_TYPE, PHOTON_MODE
  double precision   :: laser_energy, laser_variation ! Laser energy
  double precision   :: mu, sigma, A ! Gauss pulse parameters

  ! Photon energy parameters of photon energy

  ! Gauss emission
  logical, parameter :: Gauss_Emission = .FALSE.

  ! Type of laser input
  integer, parameter :: LASER_GAUSS      = 1
  integer, parameter :: LASER_SQUARE     = 2  

  ! Type of photon velocity
  integer, parameter :: PHOTON_ZERO  = 1
  integer, parameter :: PHOTON_MB    = 2
  
!  PRIVATE
!  PUBLIC Init_Photon_Velocity, PHOTON_ZERO, PHOTON_MB
    

  ! interface
  !   double precision function Laser_fun(pos, emit, sec)
  !     double precision, dimension(1:3), intent(in) :: pos
  !     integer, intent(in)                          :: emit
  !     integer, intent(out), optional               :: sec
  !   end function Laser_fun
  ! end interface
  ! procedure(Laser_fun), pointer :: ptr_Laser_fun => null()

contains

  subroutine Read_laser_parameters()
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

      ! Read the type of laser pulse to use
      read(unit=ud_laser, FMT=*) LASER_TYPE, PHOTON_MODE

      SELECT CASE (LASER_TYPE)
        case (LASER_GAUSS)
          ! Gaussian Pulse
          print '(a)', 'Vacuum: Using Gaussian pulse model'
          ! Add function to change Gauss_emission to .TRUE. 
          ! in mod_photo_emission.f90 check Image charge for ref.
          ! ptr_Laser_fun => gauss_pulse_parameters
          SELECT CASE (PHOTON_MODE)
            case (PHOTON_ZERO)    
              ptr_Get_Photon_Velocity => Get_Zero_Photon_Velocity
              print '(a)', 'Vacuum: Using zero inital velocity' 
              ! Read mean and std photon energy from the file
              read(unit=ud_laser, FMT=*) laser_energy, laser_variation
              ! Read Gaussian parameters from the file (mu, std and Amplitude) 
              read(unit=ud_laser, FMT=*) mu, sigma, A
              ! These parameters should go to the Gauss_Emission
              ! in mod_photo_emission.f90

            case (PHOTON_MB)
              ptr_Get_Photon_Velocity => Get_Photon_Energy
              print '(a)', 'Vacuum: Using Maxwell-Boltzman energy distribution for Photons'
              ! Read mean and std photon energy from the file
              read(unit=ud_laser, FMT=*) laser_energy, laser_variation
              ! Read Gaussian parameters from the file (mu, std and Amplitude) 
              read(unit=ud_laser, FMT=*) mu, sigma, A
              ! These parameters should go to the Gauss_Emission
              ! in mod_photo_emission.f90
            case DEFAULT
              print '(a)', 'Vacuum: ERROR UNKNOWN PHOTON MODE'
              print *, PHOTON_MODE
              stop
          END SELECT  


        case (LASER_SQUARE)
          ! Square Pulse
          print '(a)', 'Vacuum: Using Square pulse model'
          ptr_Laser_fun => laser
          SELECT CASE (PHOTON_MODE)
            case (PHOTON_ZERO)
              ptr_Get_Photon_Velocity => Get_Zero_Photon_Velocity
              print '(a)', 'Vacuum: Using zero inital velocity'    
              ! Read mean and std photon energy from the file
              read(unit=ud_laser, FMT=*) laser_energy, laser_variation

            case (PHOTON_MB)
              ptr_Get_Photon_Velocity => Get_Photon_Energy
              print '(a)', 'Vacuum: Using Maxwell-Boltzman energy distribution for Photons'
              ! Read mean and std photon energy from the file
              read(unit=ud_laser, FMT=*) laser_energy, laser_variation
            case DEFAULT
              print '(a)', 'Vacuum: ERROR UNKNOWN PHOTON MODE'
              print *, PHOTON_MODE
              stop
          END SELECT
        case DEFAULT
          print '(a)', 'Vacuum: ERROR UNKNOWN LASER TYPE AND/OR PHOTON MODE'
          print *, LASER_TYPE, PHOTON_MODE
          stop
      END SELECT

      close(unit=ud_laser)
    end subroutine Read_laser_parameters

    ! ! Clean up routine
    ! subroutine Laser_cleanup()
    !   if (LASER_TYPE == LASER_GAUSS) then
    !     deallocate()
    !     deallocate()
    !   endif


    !   if (LASER_TYPE == LASER_SQUARE) then
    !     deallocate()
    !     deallocate()
    !   endif
    ! end subroutine Laser_cleanup

    !double precision function laser_parameters(laser_energy, wavelength, intensity)
    !end function laser_parameters

  
    ! subroutine Init_Photon_Velocity(PHOTON_MODE)
    !     integer, intent(in) :: PHOTON_MODE
    !     SELECT CASE (PHOTON_MODE)
    !     case(PHOTON_ZERO)
    !         ptr_Get_Photon_Velocity => Get_Zero_Photon_Velocity
    !         print '(a)', 'Vacuum: Using zero inital velocity'
    !     case(PHOTON_MB)
    !         ptr_Get_Photon_Velocity => Get_Photon_Energy
    !         print '(a)', 'Vacuum: Using Maxwell-Boltzman energy distribution for Photons'
    !     case DEFAULT
    !         print '(a)', 'Vacuum: ERROR UNKNOWN PHOTON MODE'
    !         print *, PHOTON_MODE
    !         stop
    !     END SELECT
    ! end subroutine Init_Photon_Velocity
    
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
    function Get_Photon_Energy() ! (laser_energy, laser_variation)
      !double precision, intent(in)     :: laser_energy, laser_variation
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