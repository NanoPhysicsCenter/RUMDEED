!-------------------------------------------!
! Submodule for the Metropolis-Hastings     !
! algorithm                                 !
! Kristinn Torfason                         !
! 21.10.18                                  !
!-------------------------------------------!
submodule (mod_field_emission_v2) smod_metro_algo

contains
  !-----------------------------------------------------------------------------
  ! Metropolis-Hastings algorithm
  ! Includes that the work function can vary with position
  function Metropolis_Hastings_rectangle_v2(ndim, emit, df_out, F_out)
    integer, intent(in)              :: ndim, emit
    double precision, intent(out)    :: df_out, F_out
    !double precision, dimension(1:3) :: Metropolis_Hastings_rectangle_v2 ! The interface is declared in the parent module
    integer                          :: count, i
    double precision                 :: std, rnd, alpha
    double precision, dimension(1:3) :: cur_pos, new_pos, field
    double precision                 :: df_cur, df_new

    std = (emitters_dim(1, emit)*0.05d0 + emitters_dim(2, emit)*0.05d0) / (2.0d0)

    ! Get a random position on the surface
    count = 0
    do ! Infinite loop
      CALL RANDOM_NUMBER(cur_pos(1:2))
      cur_pos(1:2) = cur_pos(1:2)*emitters_dim(1:2, emit) + emitters_pos(1:2, emit)
      cur_pos(3) = 0.0d0 ! On the surface

      ! Calculate the electric field at this position
      field = Calc_Field_at(cur_pos)
      if (field(3) < 0.0d0) then
        exit ! The loop is infinite
      else
        count = count + 1
        if (count > 1000) exit ! The loop is infnite, must stop it at some point
      end if
    end do

    F_out = field(3)

    ! Calculate the escape probability at this location
    if (field(3) < 0.0d0) then
      df_cur = Escape_Prob(field(3), cur_pos)
    else
      df_cur = 0.0d0 ! Zero escape probabilty if field is not favourable
    end if

    !---------------------------------------------------------------------------
    ! We now pick a random distance and direction to jump to from our
    ! current location. We do this ndim times.
    do i = 1, ndim
      ! Find a new position
      new_pos(1:2) = cur_pos(1:2) + box_muller(0.0d0, std)

      ! Make sure that the new position is within the limits of the emitter area.
      call check_limits_metro_rec(new_pos, emit)

      ! Calculate the field at the new position
      field = Calc_Field_at(new_pos)
      F_out = field(3)

      ! Check if the field is favourable for emission at the new position.
      ! If it is not then cycle, i.e. we reject this location and
      ! pick another one.
      if (field(3) > 0.0d0) cycle ! Do the next loop iteration

      ! Calculate the escape probability at the new position, to compair with
      ! the current position.
      df_new = Escape_Prob(field(3), new_pos)

      ! If the escape probability is higher in the new location,
      ! then we jump to that location. If it is not then we jump to that
      ! location with the probabilty df_new / df_cur.
      if (df_new > df_cur) then
        cur_pos = new_pos ! New position becomes the current position
        df_cur = df_new
        F_out = field(3)
      else
        alpha = df_new / df_cur

        CALL RANDOM_NUMBER(rnd)
        if (rnd < alpha) then
          cur_pos = new_pos ! New position becomes the current position
          df_cur = df_new
          F_out = field(3)
        end if
      end if
    end do

    ! Return the current position
    Metropolis_Hastings_rectangle_v2 = cur_pos
    df_out = df_cur
  end function Metropolis_Hastings_rectangle_v2

  ! ----------------------------------------------------------------------------
  !
  subroutine check_limits_metro_rec(par_pos, emit)
    double precision, dimension(1:3), intent(inout) :: par_pos
    integer, intent(in)                             :: emit
    double precision                                :: x_max, x_min, y_max, y_min
    double precision                                :: d_x, d_y


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

end submodule smod_metro_algo