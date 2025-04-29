!-------------------------------------------!
! Module for photo emission                 !
! Kristinn Torfason                         !
! 05.04.13                                  !
!-------------------------------------------!

Module mod_therminoic_emission
  use mod_global
  use mod_verlet
  use mod_pair
  implicit none

  PRIVATE
  PUBLIC :: Init_Thermionic_Emission, Clean_Up_Thermionic_Emission, Richardson_Dushman

  ! ----------------------------------------------------------------------------
  ! Variables
  integer, dimension(:), allocatable :: nrEmitted_emitters
  integer                            :: nrElecEmitAll
  !integer                            :: nrEmitted

  ! ----------------------------------------------------------------------------
  ! Parameters
  double precision, parameter :: lambda_R = 0.5d0 ! Material paramater []
  double precision, parameter :: A_0 = 4.0d0*pi*m_0*k_b**2*q_0 / h**3 ! Constant for the RD eq. [A m^-2 K^-2]
  double precision, parameter :: A_G = lambda_R*A_0 ! Constant for the RD eq. [A m^-2 K^-2]
  double precision, parameter :: B_Sch_eV_prefix = q_0/(4.0d0*pi*epsilon_0) ! Schottky effect prefix [m V]

  double precision :: w_theta = 2.0d0 ! Work function [eV]
  double precision :: T_k = 100.0d0 ! Temperature [K]
  double precision :: k_b_eV = 8.6173303d-5 ! Boltzman constant [eV K^-1]

  ! Constant used in MC integration (function Elec_Supply_V2)
  double precision :: time_step_div_q0

  double precision, dimension(1:3)   :: F_avg = 0.0d0

  interface 
      ! Interface for the work function submodule
      !double precision module function w_theta_xy(pos, sec)
      !  double precision, intent(in), dimension(1:3) :: pos ! Position on the surface
      !  integer, intent(out), optional               :: sec ! Return the section
      !end function w_theta_xy

      ! Interface for the MC integration submodule
      module subroutine Do_Surface_Integration_TE(emit, N_sup)
        integer, intent(in)  :: emit ! The emitter to do the integration on
        integer, intent(out) :: N_sup ! Number of electrons
      end subroutine Do_Surface_Integration_TE

      ! Interface for the Metropolis-Hastings submodule
      !module function Metropolis_Hastings_rectangle_v2(ndim, emit, df_out, F_out)
      !  integer, intent(in)              :: ndim, emit
      !  double precision, intent(out)    :: df_out, F_out
      !  double precision, dimension(1:3) :: Metropolis_Hastings_rectangle_v2
      !end function Metropolis_Hastings_rectangle_v2
  end interface

contains
  subroutine Init_Thermionic_Emission()
! Allocate the number of emitters
    allocate(nrEmitted_emitters(1:nrEmit))

    ! Initialize variables
    nrEmitted_emitters = 0 ! Number of electrons emitted from emitter number i in this time step

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_Planar

    ! Function for the electric field in the system
    ptr_field_E => field_E_planar

    ! The function that does the emission
    ptr_Do_Emission => Do_Thermionic_Emission

    !call Read_work_function()

    ! Parameters used in the module
    !time_step_div_q0 = time_step / q_0

    ! Initialize the Ziggurat algorithm
    call zigset(my_seed(1))
  end subroutine Init_Thermionic_Emission

  subroutine Clean_Up_Thermionic_Emission()
    deallocate(nrEmitted_emitters)
  end subroutine Clean_Up_Thermionic_Emission

  subroutine Do_Thermionic_Emission(step)
    integer, intent(in) :: step
    integer             :: i, IFAIL
    double precision    :: cur_time

    nrElecEmitAll = 0
    nrEmitted_emitters = 0

    ! Loop through all of the emitters
    do i = 1, nrEmit

      ! Check the type of the emitter CIRCLE / RECTANGLE
      if (emitters_delay(i) < step) then
        call Do_Thermionic_Emission_Rectangle(step, i)
      end if
    end do

    cur_time = time_step * step / time_scale ! Scaled in units of time_scale
    write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, nrElecEmitAll, &
                                                       & nrElec, (nrEmitted_emitters(i), i = 1, nrEmit)
  end subroutine Do_Thermionic_Emission

  subroutine Do_Thermionic_Emission_Rectangle(step, nrEmit)
    integer, intent(in) :: step, nrEmit
    integer             :: N_sup ! Number of electrons to emit

    call Do_Surface_Integration_TE(nrEmit, N_sup)

  end subroutine Do_Thermionic_Emission_Rectangle

  double precision pure function Escape_Prob(F, pos)
    double precision,                 intent(in) :: F
    double precision, dimension(1:3), intent(in) :: pos
    double precision                             :: delta_W ! Schottky effect

    delta_W = sqrt(B_Sch_eV_prefix*(-1.0d0*F)) ! Units: mV*V/m = V^2, sqrt(V^2) = V = J/C = eV
    Escape_Prob = exp( -1.0d0*(w_theta - delta_W)/(k_b_eV*T_k) )

  end function Escape_Prob

  double precision pure function Elec_Supply(F, pos)
    double precision,                 intent(in) :: F
    double precision, dimension(1:3), intent(in) :: pos

    Elec_Supply = A_G * T_k**2
  end function Elec_Supply

  double precision pure function Richardson_Dushman(F, pos)
    double precision, dimension(1:3), intent(in) :: pos     ! Position of particle
    double precision                , intent(in) :: F       ! Field strength at pos

    Richardson_Dushman = Elec_Supply(F, pos) * Escape_Prob(F, pos)
  end function Richardson_Dushman
end module mod_therminoic_emission
