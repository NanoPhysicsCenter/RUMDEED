program Parallel_RLC
    implicit none

    integer :: k, ud_volt, IFAIL
    integer, parameter :: N = 10000
    double precision :: I_cur = 0.0d0, I_prev = 0.0d0, V_cur = 0.0d0, V_prev = 0.0d0, I = 0.0d0, V = 0.0d0, V_d

    open(newunit=ud_volt, iostat=IFAIL, file='volt.dt', status='REPLACE', action='write')

    do k = 1, N
        I = Get_Current(k)
        I_prev = I_cur
        I_cur = I

        V = Parallel_RLC_FD(k, I_cur, I_prev, V_cur, V_prev)
        V_prev = V_cur
        V_cur = V
        V_d = 2.0d0 + V

        write (ud_volt, "(i8, tr2, ES18.8, tr2, ES18.8)", iostat=IFAIL) k, V_d, V

        if (k == 1) then
            print *, V
        end if
    end do

    close(unit=ud_volt, status='keep')
contains
    ! Current
    double precision function Get_Current(step)
        integer, intent(in) :: step
        Get_Current = 10.0E-3 ! Amper
    end function Get_Current

    ! Voltage
    double precision function Parallel_RLC_FD(step, I_cur, I_prev, V_cur, V_prev)
        integer, intent(in) :: step
        double precision, intent(in) :: V_cur, V_prev ! Current V(t) and previous V(t-Δt) values of the voltage
        double precision, intent(in) :: I_cur, I_prev ! Current I(t) and previous I(t-Δt) values of the current

        !double precision, parameter  :: R = 2.0d0  ! Ohm
        !double precision, parameter  :: L = 1.0d-7  ! Henry
        !double precision, parameter  :: C = 1.0d-22 ! Farad

        double precision, parameter  :: R = 13.5d0  ! Ohm
        double precision, parameter  :: L = 1.04d-7  ! Henry
        double precision, parameter  :: C = 2.54d-14 ! Farad

        double precision, parameter :: time_step = 0.25E-15
        double precision, parameter :: time_step2 = time_step**2


        double precision :: V_next ! Next value of the voltage V(t+Δt)

        V_next = 1.0*time_step/C * I_cur + V_cur * (2.0d0 + time_step/(R*C) - time_step2/(L*C)) &
               & - V_prev - time_step/C * I_prev
        V_next = V_next / (1.0d0 + time_step/(R*C))

        Parallel_RLC_FD = V_next
        if (step == 1) then
            print *, I_cur
            print *, I_prev
            print *, V_cur
            print *, V_prev
            print *, time_step

            print *, time_step/C
            print *, (2.0d0 + time_step/(R*C) - time_step2/(L*C))
            print *, (1.0d0 + time_step/(R*C))
        end if
    end function Parallel_RLC_FD
end program Parallel_RLC