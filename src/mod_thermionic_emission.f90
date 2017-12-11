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

  ! ----------------------------------------------------------------------------
  ! Variables
  integer, dimension(:), allocatable :: nrEmitted_emitters
  integer                                     :: posInit
  logical                                     :: EmitGauss = .false.
  integer                                     :: maxElecEmit = -1
  integer                                     :: nrEmitted

  integer                                     :: ud_gauss

  ! ----------------------------------------------------------------------------
  ! Parameters
  integer, parameter                          :: MAX_EMISSION_TRY = 100

  double precision, parameter :: A_g = 4.0d0*pi*m_e*k_b**2*(-1.0d0*q_e) / h**3
  double precision, parameter :: B_Sch_eV_prefix = sqrt(-1.0d0*q_e**3/(4.0d0*pi*epsilon_0)) / (-1.0d0*q_e)

contains
  subroutine Init_Photo_Emission()
    allocate(nrEmitted_emitters(1:nrEmit))
  end subroutine Init_Photo_Emission

  subroutine Clean_Up_Photo_Emission()
    deallocate(nrEmitted_emitters)
  end subroutine Clean_Up_Photo_Emission

  subroutine Do_Thermionic_Emission(step)
    integer, intent(in) :: step
    integer             :: i, IFAIL
    double precision    :: cur_time
    double precision, dimension(1:3) :: pos

    posInit = 0
    nrEmitted_emitters = 0

    ! Loop through all of the emitters
    do i = 1, nrEmit

      ! Check the type of the emitter CIRCLE / RECTANGLE
      if (emitters_delay(i) < step) then
        call Do_Thermionic_Emission_Rectangle(step, i)
      end if
    end do

    cur_time = time_step * step / time_scale ! Scaled in units of time_scale
    write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, posInit, &
    & nrEmitted, nrElec, (nrEmitted_emitters(i), i = 1, nrEmit), maxElecEmit
  end subroutine Do_Thermionic_Emission

  subroutine Do_Thermionic_Emission_Rectangle(step, nrEmit)
    integer, intent(in) :: step, nrEmit

    double precision                 :: F
    double precision, dimension(1:3) :: par_pos
    !double precision, allocatable, dimension(:) :: rnd
    integer                          :: i, j, s, IFAIL, nrElecEmit, n_r
    double precision                 :: A_f, n_s, n_add, res_s
    double precision                 :: len_x, len_y
    integer                          :: nr_x, nr_y

    nr_x = 1000
    nr_y = nr_x
    len_x = emitters(emit)%dim(1) / nr_x
    len_y = emitters(emit)%dim(2) / nr_y

    A_f = len_x*len_y

    nrElecEmit = 0

    par_pos = 0.0d0
    n_s = 0.0d0

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, s, par_pos, par_elec, F, n_add, D_f)

    !$OMP DO REDUCTION(+:n_s,F_avg)
    do i = 1, nr_x
      do j = 1, nr_y
        par_pos(1) = (i - 0.5d0)*len_x + emitters(emit)%pos(1)
        par_pos(2) = (j - 0.5d0)*len_y + emitters(emit)%pos(2)
        par_pos(3) = 0.0d0

        field = Calc_Field_at(par_pos)
        F = field(3)

        if (F >= 0.0d0) then
          n_add = 0.0d0
        else
          n_add = Richardson_Dushman(par_pos, F, A)
        end if

        n_s = n_s + n_add
      end do
    end do
    !$OMP END DO

    print *, 'n_s = ', n_s
    pause


    !$OMP SINGLE
    n_s = n_s - res_s
    n_r = nint(n_s)
    res_s = n_r - n_s
    !$OMP END SINGLE

    ! !$OMP DO REDUCTION(+:df_avg)
    ! do s = 1, n_r
    !
    !   !par_pos(1:2) = Metro_algo_circle(30)
    !   par_pos(1:2) = Metro_algo_rec(30, emit)
    !   !if (isnan(par_pos(1)) == .true.) then
    !   !  print *, par_pos
    !   !  pause
    !   !end if
    !   par_pos(3) = 0.0d0*length
    !   !particles(nrElec+1)%cur_pos = par_pos
    !   !call Calculate_Acceleration(nrElec+1)
    !   par_elec%cur_pos = par_pos
    !   call Calculate_Acceleration_particle(par_elec)
    !   F = par_elec%accel(3) * pre_fac_a
    !   !F = E_z
    !
    !   !print *, 'F = ', F
    !   !print *, 's = ', s
    !   !print *, 'nrElec = ', nrElec
    !   !pause
    !
    !   if (F >= 0.0d0) then
    !     D_f = 0.0d0
    !     !print *, 'Warning: F > 0.0d0'
    !   else
    !     D_f = Escape_Prob(F, par_pos)
    !     if (D_f > 1.0d0) then
    !       print *, 'Warning D_f > 1.0d0'
    !       print *, 'D_f = ', D_f
    !     end if
    !   end if
    !
    !   CALL RANDOM_NUMBER(rnd)
    !   if (rnd <= D_f) then
    !     par_pos(3) = 1.0d0*length
    !     !$OMP CRITICAL
    !     particles(nrElec+1)%cur_pos = par_pos
    !     particles(nrElec+1)%in_step = step
    !     nrElec = nrElec + 1
    !     nrElecEmit = nrElecEmit + 1
    !     call Add_Plane_Graph_emit(par_pos, step)
    !     !call Add_Plane_Graph_emitt_xy(par_pos)
    !     !$OMP END CRITICAL
    !   end if
    ! end do
    ! !$OMP END DO

    !$OMP END PARALLEL

    posInit = posInit + nrElecEmit
    nrEmitted = nrEmitted + nrElecEmit
  end subroutine Do_Thermionic_Emission_Rectangle

  double precision pure function Richardson_Dushman(pos, F, A_emit)
    double precision, dimension(1:3), intent(in) :: pos     ! Position of particle
    double precision                , intent(in) :: F       ! Field strength at pos
    double precision,                 intent(in) :: A_emit  ! Emitter area
    double precision                             :: delta_W ! Schottky effect
    double precision                             :: n_s, escape_prob

    delta_W = B_Sch_eV_prefix*sqrt(-1.0d0*F)

    n_s = A_g * T_k**2 * A_emit * time_step / (-1.0d0*q_e)
    escape_prob = exp(-(w_theta - delta_W)/(k_b_eV*T_k) )

    Richardson_Dushman = n_s * escape_prob
  end function Richardson_Dushman
end module mod_therminoic_emission
