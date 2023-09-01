!-------------------------------------------!
! Module for global variables and constants !
! Kristinn Torfason                         !
! 07.05.15                                  !
!-------------------------------------------!

module mod_velocity
use mod_global
implicit none

integer, parameter :: VELOCITY_ZERO   = 1
integer, parameter :: VELOCITY_MB     = 2
integer, parameter :: VELOCITY_PHOTON = 3

PRIVATE
PUBLIC Init_Emission_Velocity, VELOCITY_ZERO, VELOCITY_MB

contains

    subroutine Init_Emission_Velocity(VELOCITY_MODE)
        integer, intent(in) :: VELOCITY_MODE
        SELECT CASE (VELOCITY_MODE)
        case(VELOCITY_ZERO)
            ptr_Get_Emission_Velocity => Get_Zero_Velocity
            print '(a)', 'RUMDEED: Using zero inital velocity'
        case(VELOCITY_MB)
            ptr_Get_Emission_Velocity => Get_MB_Velocity
            print '(a)', 'RUMDEED: Using Maxwell-Boltzman velocity distribution'
        case(VELOCITY_PHOTON)
            ptr_Get_Emission_Velocity => Get_Photon_Velocity
            print '(a)', 'RUMDEED: Using photon energy'
        case DEFAULT
            print '(a)', 'RUMDEED: ERROR UNKNOWN VELOCITY MODEL'
            print *, VELOCITY_MODE
            stop
        END SELECT
    end subroutine Init_Emission_Velocity

    ! ----------------------------------------------------------------------------
    ! Just give zero initial velocity
    function Get_Zero_Velocity()
        double precision, dimension(1:3) :: Get_Zero_Velocity

        Get_Zero_Velocity = 0.0d0
    end function Get_Zero_Velocity

    ! ----------------------------------------------------------------------------
    ! Generate velocity from a Maxwell-Boltzmann distribution.
    ! https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution#Distribution_for_the_velocity_vector
    ! 
    ! Maxwell-Boltzmann velocity distribution is basically a normal distribution for each velocity component with
    ! zero mean and standard deviation \sqrt{k_b T / m}.
    function Get_MB_Velocity()
        double precision, dimension(1:3) :: Get_MB_Velocity
        double precision, dimension(1:2) :: std
        double precision, dimension(1:2) :: mean

        mean = 0.0d0
        std = sqrt(k_b * T_temp / m_0) ! Standard deviation of the Maxwell-Boltzmann distribution

        ! Get normal distributed numbers.
        ! The Box Muller method gives two numbers.
        ! We overwrite the second element in the array.
        Get_MB_Velocity(1:2) = box_muller(mean, std)
        Get_MB_Velocity(2:3) = box_muller(mean, std)

        Get_MB_Velocity(3) = abs(Get_MB_Velocity(3)) ! Positive velocity in the z-direction
    end function Get_MB_Velocity

    ! ----------------------------------------------------------------------------
    ! Generate velocity from a depentant on the photon energy and work function.
    ! *TODO* Is implemented in the photo_emission module. Move it here.
    function Get_Photon_Velocity()
        double precision, dimension(1:3) :: Get_Photon_Velocity

        Get_Photon_Velocity(1:2) = 0.0d0
        Get_Photon_Velocity(3) = 0.0d0 !sqrt((2.0d0 * ((p_eV - w_theta_xy(par_pos, emit))*q_0))/m_0)
    end function Get_Photon_Velocity

end module mod_velocity