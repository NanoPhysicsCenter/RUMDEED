!-------------------------------------------!
! Module for emission                       !
! Kristinn Torfason                         !
! 05.04.13                                  !
!-------------------------------------------!

Module mod_field_emission
  use mod_global
  !use mod_hyperboloid_tip
  use mod_verlet
  use mod_pair
  implicit none

  ! ----------------------------------------------------------------------------
  ! Variables
  integer, dimension(:), allocatable :: nrEmitted_emitters
  integer                            :: posInit
  integer                            :: nrEmitted
  double precision                   :: res_s ! residual

  ! ----------------------------------------------------------------------------
  ! Constants for field emission
  double precision, parameter :: a_FN = q_02/(16.0d0*pi**2*h_bar) ! A eV V^{-2}
  double precision, parameter :: b_FN = -4.0d0/(3.0d0*h_bar) * sqrt(2.0d0*m_0*q_0) ! eV^{-3/2} V m^{-1}
  double precision, parameter :: l_const = q_0 / (4.0d0*pi*epsilon_0) ! eV^{2} V^{-1} m
  double precision, parameter :: w_theta = 2.0d0 ! work function in eV

  ! Image Charge
  logical, parameter          :: image_charge = .true.

contains
  !-----------------------------------------------------------------------------
  ! Initialize the Field Emission
  subroutine Init_Field_Emission()
    ! Allocate the number of emitters
    allocate(nrEmitted_emitters(1:nrEmit))

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_ElecHole_Planar

    ! Function for the electric field in the system
    ptr_field_E => field_E_planar
  end subroutine Init_Field_Emission

  !-----------------------------------------------------------------------------
  ! This subroutine gets called from main when the emitters should emit the electrons
  subroutine Do_Field_Emission(step)
    integer, intent(in) :: step
    integer             :: i, IFAIL
    double precision    :: cur_time

    posInit = 0
    nrEmitted_emitters = 0

    ! Loop through all of the emitters
    do i = 1, nrEmit

      ! Check the type of the emitter CIRCLE / RECTANGLE
      if (emitters_delay(i) < step) then
        !select case (emitters_type(i))
        !case (EMIT_CIRCLE)
          !print *, 'Doing Circle'
        !  call Do_Photo_Emission_Circle(step, i)
        !case (EMIT_RECTANGLE)
          !print *, 'Doing Rectangle'
          call Do_Field_Emission_Plane_int_rec(step, i)
        !case default
        !  print *, 'Vacuum: WARNING unknown emitter type!!'
        !  print *, emitters_type(i)
        !end select
      end if
    end do

    cur_time = time_step * step / time_scale ! Scaled in units of time_scale
    write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, posInit, &
    & nrEmitted, nrElec, (nrEmitted_emitters(i), i = 1, nrEmit)
  end subroutine Do_Field_Emission

  !-----------------------------------------------------------------------------
  ! Do the field emission for a planar rectangular emitter
  ! step is the current time step
  ! emit is the number of the emitter
  subroutine Do_Field_Emission_Plane_int_rec(step, emit)
    integer, intent(in)              :: step, emit
    double precision                 :: F
    double precision, dimension(1:3) :: par_pos, par_vel ! Position and velocity of the electron to be emitted
    !double precision, allocatable, dimension(:) :: rnd
    integer                          :: i, j, s, nrElecEmit, n_r
    double precision                 :: A_f, D_f, n_s, n_add
    double precision                 :: len_x, len_y
    integer                          :: nr_x, nr_y
    double precision, dimension(1:3) :: field, F_avg
    double precision                 :: df_avg, rnd

    par_pos = 0.0d0
    nrElecEmit = 0
    nrEmitted_emitters(emit) = 0

    A_f = len_x*len_y ! Total area of the emitter

    nr_x = 1000 ! Divide the area into this many sections in the x-direction
    nr_y = nr_x ! Divide the area into this many sections in the y-direction

    len_x = emitters_dim(1, emit) / nr_x ! Size of each section in x
    len_y = emitters_dim(2, emit) / nr_y ! Size of each section in y

    n_s = 0.0d0 ! Set the number of electrons to be emitted in this time step to zero.
    F_avg = 0.0d0 ! Set the avereage field to zero before we start.

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, s, par_pos, F, n_add, D_f)

    !$OMP DO REDUCTION(+:n_s,F_avg)
    do i = 1, nr_x
      do j = 1, nr_y

        ! Find the mid point of each section
        par_pos(1) = (i - 0.5d0)*len_x + emitters_pos(1, emit)
        par_pos(2) = (j - 0.5d0)*len_y + emitters_pos(2, emit)
        par_pos(3) = 0.0d0

        ! Calculate the field in the midpoint of each section
        par_pos(3) = 0.0d0 * length_scale !Check in plane
        field = Calc_Field_at(par_pos)

        F = field(3) ! Take the z-value of the field

        ! Calculate the average field. (Division with number of sections is done after the loop).
        F_avg(1:3) = F_avg(1:3) + field(1:3)

        ! Check if the field is favourable.
        if (F >= 0.0d0) then
          n_add = 0.0d0 ! Do not emit any electrons from this section.
        else
          ! Calculate the electron supply to see how many electrons to emit from this section.
          n_add = Elec_supply(A_f, F, par_pos)
        end if

        ! Add to the total number of electrons to be emitted.
        n_s = n_s + n_add
      end do
    end do
    !$OMP END DO


    !$OMP SINGLE
    ! Finish calculating the average field.
    F_avg(1:3) = F_avg(1:3) / (nr_x*nr_y)

    ! The number of electrons we are going to emit has to be a integer (A whole number).
    ! We round to the nearest integer and store the residual.
    ! We add the resdual from the previous time step for accurate accounting.
    n_s = n_s - res_s
    n_r = nint(n_s) ! round the number
    res_s = n_r - n_s

    ! Set the average escape probability to zero.
    df_avg = 0.0d0
    !$OMP END SINGLE


    ! Loop over the electrons to be emitted.
    !$OMP DO REDUCTION(+:df_avg)
    do s = 1, n_r

      par_pos(1:2) = Metro_algo_rec(30, emit)
      par_pos(3) = 0.0d0 * length_scale !Check in plane
      field = Calc_Field_at(par_pos)

      F = field(3)

      if (F >= 0.0d0) then
        D_f = 0.0d0
        !print *, 'Warning: F > 0.0d0'
      else
        D_f = Escape_Prob(F, par_pos)
        if (D_f > 1.0d0) then
          print *, 'Warning D_f > 1.0d0'
          print *, 'D_f = ', D_f
        end if
      end if
      df_avg = df_avg + D_f

      CALL RANDOM_NUMBER(rnd)
      if (rnd <= D_f) then
        par_pos(3) = 1.0d0*length_scale
        !$OMP CRITICAL

          ! Add a particle to the system
          par_vel = 0.0d0
          call Add_Particle(par_pos, par_vel, species_elec, step)

          nrElecEmit = nrElecEmit + 1
          nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
          !call Add_Plane_Graph_emit(par_pos, step)
          !call Add_Plane_Graph_emitt_xy(par_pos)
        !$OMP END CRITICAL
      end if
    end do
    !$OMP END DO

    !$OMP END PARALLEL


    !deallocate(rnd)

    df_avg = df_avg / n_r

    !write (ud_debug, "(i8, tr2, E16.8, tr2, E16.8, tr2, E16.8, tr2, i8, tr2, E16.8)", iostat=IFAIL) &
    !                                  step, F_avg(1), F_avg(2), F_avg(3), n_r, df_avg

    posInit = posInit + nrElecEmit
    nrEmitted = nrEmitted + nrElecEmit

  end subroutine Do_Field_Emission_Plane_int_rec

!----------------------------------------------------------------------------------------
  double precision function v_y(F, pos)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos
    double precision                             :: l

    if (image_charge .eqv. .true.) then
      l = l_const * (-1.0d0*F) / w_theta_xy(pos)**2 ! l = y^2, y = 3.79E-4 * sqrt(F_old) / w_theta
      if (l > 1.0d0) then
        l = 1.0d0
      end if
      v_y = 1.0d0 - l + 1.0d0/6.0d0 * l * log(l)
    else
      v_y = 1.0d0
    end if
  end function v_y

  double precision function t_y(F, pos)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos
    double precision                             :: l

    if (image_charge .eqv. .true.) then
      l = l_const * (-1.0d0*F) / w_theta_xy(pos)**2 ! l = y^2, y = 3.79E-4 * sqrt(F_old) / w_theta
      if (l > 1.0d0) then
        print *, 'Error: l > 1.0'
        print *, 'l = ', l, ', F = ', F, ', t_y = ', t_y
        print *, 'x = ', pos(1), 'y = ', pos(2)
        l = 1.0d0
        !call Write_Current_Position()
        !stop
      end if

      t_y = 1.0d0 + l*(1.0d0/9.0d0 - 1.0d0/18.0d0*log(l))
    else
      t_y = 1.0d0
    end if
  end function t_y

  double precision function Elec_Supply(A, F, pos)
    double precision, intent(in)                 :: A, F
    double precision, dimension(1:3), intent(in) :: pos
    double precision                             :: n

    n = A * a_FN * F**2 * time_step / (1.0d0*q_0 * w_theta_xy(pos) * (t_y(F, pos))**2)

    Elec_supply = n
  end function Elec_supply

  double precision function Escape_Prob(F, pos)
    double precision, intent(in)                 :: F
    double precision, dimension(1:3), intent(in) :: pos

    Escape_Prob = exp(b_FN * (sqrt(w_theta_xy(pos)))**3 * v_y(F, pos) / (-1.0d0*F))

    if (Escape_Prob > 1.0d0) then
      print *, 'Escape_prob is larger than 1.0'
      print *, 'Escape_prob = ', Escape_prob
      print *, ''
    end if
  end function Escape_Prob

  function Metro_algo_rec(ndim, emit)
    integer, intent(in)                          :: ndim, emit
    double precision, dimension(1:2)             :: Metro_algo_rec
    double precision, dimension(1:3)             :: par_pos, new_pos, field
    double precision                             :: old_val, new_val, alpha, std, rnd
    integer                                      :: old_val_s, new_val_s
    integer                                      :: i

    std = (emitters_dim(1, emit)*0.05d0 + emitters_dim(2, emit)*0.05d0) / (2.0d0)

    CALL RANDOM_NUMBER(par_pos)

    par_pos(1:2) = par_pos(1:2)*emitters_dim(1:2, emit) + emitters_pos(1:2, emit)

    field = Calc_Field_at(par_pos)*(-1.0d0) ! This code assumes it is using the acceleration

    old_val = norm2(field)
    if (field(3) > 0.0d0) then
      old_val_s = +1
    else
      old_val_s = -1
    end if

    new_pos = 0.0d0
    do i = 1, ndim
      !IFAIL = vdRngGaussianMV(VSL_RNG_METHOD_GAUSSIANMV_ICDF, stream, 1, r(:, 1:2), 2, VSL_MATRIX_STORAGE_FULL, par_pos(1:2), T)
      new_pos(1:2) = par_pos(1:2) + box_muller(0.0d0, std)
      call check_limits_metro_rec(new_pos, emit)

      field = Calc_Field_at(par_pos)*(-1.0d0) ! This code assumes it is using the acceleration
      new_val = norm2(field)

      if (field(3) > 0.0d0) then
        new_val_s = +1
      else
        new_val_s = -1
      end if

      if ((old_val_s < 0) .and. (new_val_s < 0)) then
        alpha = old_val / new_val
      else if ((old_val_s > 0) .and. (new_val_s < 0)) then
        alpha = 0.0d0 ! Do not jump to this point
      else if ((old_val_s < 0) .and. (new_val_s > 0)) then
        alpha = 1.0d0 ! Do jump to this point
      else
        alpha = new_val / old_val
      end if

      alpha = new_val / old_val
      CALL RANDOM_NUMBER(rnd)
      if (rnd < alpha) then
        old_val = new_val
        old_val_s = new_val_s
        par_pos = new_pos
      end if

    end do

    Metro_algo_rec = par_pos(1:2)
  end function Metro_algo_rec

  subroutine check_limits_metro_rec(par_pos, emit)
    double precision, dimension(1:3), intent(inout) :: par_pos
    integer, intent(in)                             :: emit
    double precision                                :: x_max, x_min, y_max, y_min
    double precision                                :: d_x, d_y


!     if ((new_pos(1) > emitters(1)%dim(1)) .or. (new_pos(1) < 0.0d0)) then
!       new_pos(1) = par_pos(1)
!     end if
!     if ((new_pos(2) > emitters(1)%dim(2)) .or. (new_pos(2) < 0.0d0)) then
!       new_pos(2) = par_pos(2)
!     end if

    x_max = emitters_pos(1, emit) + emitters_dim(1, emit)
    x_min = emitters_pos(1, emit)

    y_max = emitters_pos(2, emit) + emitters_dim(2, emit)
    y_min = emitters_pos(2, emit)

    !Check x ----------------------------------------
    if (par_pos(1) > x_max) then
      d_x = par_pos(1) - x_max
      par_pos(1) = x_max - d_x

      if(d_x > emitters_dim(1, emit)) then
        print *, 'Warning: d_x to large >'
        print *, d_x
      end if
    end if

    if (par_pos(1) < x_min) then
      d_x = x_min - par_pos(1)
      par_pos(1) = d_x + x_min

      if(d_x > emitters_dim(1, emit)) then
        print *, 'Warning: d_x to large <'
        print *, d_x
      end if
    end if

    !Check y ----------------------------------------
    if (par_pos(2) > y_max) then
      d_y = par_pos(2) - y_max
      par_pos(2) = y_max - d_y

      if(d_y > emitters_dim(2, emit)) then
        print *, 'Warning: d_y to large >'
        print *, d_y
      end if
    end if

    if (par_pos(2) < y_min) then
      d_y = y_min - par_pos(2)
      par_pos(2) = d_y + y_min

      if(d_y > emitters_dim(2, emit)) then
        print *, 'Warning: d_x to large <'
        print *, d_y
      end if
    end if
  end subroutine check_limits_metro_rec


  ! subroutine Add_Plane_Graph_emitt_xy(pos)
  !   double precision, dimension(1:3), intent(in) :: pos
  !   double precision                             :: x, y
  !   integer                                      :: IFAIL
  !
  !   x = pos(1) / length_scale
  !   y = pos(2) / length_scale
  !
  !   write (ud_plane_graph_emit, "(F7.2, tr2, F7.2)", iostat=IFAIL) x, y
  ! end subroutine Add_Plane_Graph_emitt_xy
  !
  ! subroutine Add_Plane_Graph_emit(pos, step)
  !   double precision, dimension(1:3), intent(in) :: pos
  !   integer, intent(in) :: step
  !   integer :: int_x, int_y
  !
  !   if ((pos(1) < int_end_x_emit) .and. (pos(1) > int_start_x_emit) &
  !      .and. (pos(2) < int_end_y_emit) .and. (pos(2) > int_start_y_emit)) then
  !     int_x = ceiling((pos(1)+abs(int_start_x_emit))/dx_emit)
  !     int_y = ceiling((pos(2)+abs(int_start_y_emit))/dy_emit)
  !
  !     plane_graph_mesh_emit(int_x, int_y) = plane_graph_mesh_emit(int_x, int_y) + 1
  !   else
  !     print *, 'Vacuum: Add_Plane_Graph_emit() error particle outside range'
  !     print *, 'x = ', pos(1)
  !     print *, 'y = ', pos(2)
  !     print *, 'step = ', step
  !     print *, ''
  !   end if
  ! end subroutine Add_Plane_Graph_emit
  !
  ! subroutine Add_kdens_Graph_emit(pos, step)
  !   double precision, dimension(1:3), intent(in) :: pos
  !   integer, intent(in) :: step
  !   integer          :: i, j
  !   double precision :: x, y
  !
  !   do i = 1, (res_x_kdens+1)
  !     x = dx_kdens_emit * (i-1) + int_start_x_emit
  !     x = (x - pos(1)) / bandwidth_emit
  !
  !     do j = 1, (res_y_kdens+1)
  !       y = dy_kdens_emit * (j-1) + int_start_y_emit
  !       y = (y - pos(2)) / bandwidth_emit
  !
  !       plane_graph_kdens_emit(i, j) = plane_graph_kdens_emit(i, j) + exp(-0.5d0*(x**2 + y**2))
  !     end do
  !
  !   end do
  ! end subroutine Add_kdens_Graph_emit
  !
  ! subroutine Write_Plane_Graph_emit()
  !   integer :: i, j, IFAIL
  !   double precision :: x, y
  !
  !   !do j = 1, res_y
  !   !    write (ud_plane_graph_matlab_emit, "(*(E16.8E3, tr8))", iostat=IFAIL) (plane_graph_mesh_emit(i, j), i = 1, res_x)
  !   !end do
  !
  !   !do j = 1, (res_y_kdens+1)
  !   !    write (ud_kdens_graph_matlab_emit, "(*(E16.8E3, tr8))", iostat=IFAIL) (plane_graph_kdens_emit(i, j), i = 1, (res_x_kdens+1))
  !   !end do
  !
  !   !rewind(ud_plane_graph_emit)
  !
  !   do i = 1, res_x
  !     x = ((i - res_x*0.5d0) + 0.5d0)*dx_emit
  !     x = x / length_scale
  !     do j = 1, res_y
  !       y = ((j - res_y*0.5d0) + 0.5d0)*dy_emit
  !       y = y / length_scale
  !       !if (plane_graph_mesh(i, j) /= 0) then
  !         write (ud_plane_graph_emit, "(i6, tr2, i6, tr2, E16.8, tr2, E16.8, tr2, i6)", iostat=IFAIL) &
  !                                     i, j, x, y, plane_graph_mesh_emit(i, j)
  !       !end if
  !     end do
  !
  !     write (ud_plane_graph_emit, *) ''
  !   end do
  !
  !   !rewind(ud_kdens_graph)
  !
  !   do i = 1, (res_x_kdens+1)
  !     x = dx_kdens_emit * (i-1) + int_start_x_emit
  !     x = x / length_scale
  !
  !     do j = 1, (res_y_kdens+1)
  !       y = dy_kdens_emit * (j-1) + int_start_y_emit
  !       y = y / length_scale
  !
  !       write (ud_kdens_graph_emit, "(i6, tr2, i6, tr2, E16.8, tr2, E16.8, tr2, E16.8E3)", iostat=IFAIL) &
  !                                      i, j, x, y, plane_graph_kdens_emit(i, j)
  !     end do
  !
  !     write (ud_kdens_graph_emit, *) ''
  !   end do
  ! end subroutine Write_Plane_Graph_emit

  double precision function w_theta_xy(pos)
    double precision, dimension(1:3), intent(in) :: pos

    w_theta_xy = w_theta

    !if (pos(1) > pos(2)) then
    !  w_theta_xy = 2.0d0
    !else
    !  w_theta_xy = 2.4d0
    !end if
  end function w_theta_xy

end Module mod_field_emission
