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
  module function Metropolis_Hastings_rectangle_v2(ndim, emit, df_out, F_out)
    integer, intent(in)              :: ndim, emit
    double precision, intent(out)    :: df_out, F_out
#if defined(__INTEL_COMPILER)
    double precision, dimension(1:3) :: Metropolis_Hastings_rectangle_v2 ! The interface is declared in the parent module
#endif
    integer                          :: count, i
    double precision                 :: rnd, alpha
    double precision, dimension(1:2) :: std
    double precision, dimension(1:3) :: cur_pos, new_pos, field
    double precision                 :: df_cur, df_new

    std(1:2) = emitters_dim(1:2, emit)*0.005d0 ! Standard deviation for the normal distribution is 2.5% of the emitter length.
    ! This means that 68% of jumps are less than this value.
    ! The expected value of the absolute value of the normal distribution is std*sqrt(2/pi).

    ! Get a random initial position on the surface.
    ! We pick this location from a uniform distribution.
    count = 0
    do ! Infinite loop, we try to find a favourable position to start from
      CALL RANDOM_NUMBER(cur_pos(1:2))
      cur_pos(1:2) = cur_pos(1:2)*emitters_dim(1:2, emit) + emitters_pos(1:2, emit)
      cur_pos(3) = 0.0d0 ! On the surface

      ! Calculate the electric field at this position
      field = Calc_Field_at(cur_pos)
      if (field(3) < 0.0d0) then
        exit ! We found a nice spot so we exit the loop
      else
        count = count + 1
        if (count > 10000) exit ! The loop is infnite, must stop it at some point.
        ! In field emission it is rare the we reach the CL limit.
      end if
    end do

    F_out = field(3)

    ! Calculate the escape probability at this location
    if (field(3) < 0.0d0) then
      df_cur = -1.0d0*log(Escape_Prob(field(3), cur_pos))
    else
      df_cur = 0.0d0 ! Zero escape probabilty if field is not favourable
    end if

    !---------------------------------------------------------------------------
    ! We now pick a random distance and direction to jump to from our
    ! current location. We do this ndim times.
    do i = 1, ndim
      ! Find a new position using a normal distribution.
      new_pos(1:2) = ziggurat_normal(cur_pos(1:2), std)
      !new_pos(1:2) = cur_pos(1:2) + box_muller(0.0d0, std)

      ! Make sure that the new position is within the limits of the emitter area.
      call check_limits_metro_rec(new_pos, emit)

      ! Calculate the field at the new position
      field = Calc_Field_at(new_pos)

      ! Check if the field is favourable for emission at the new position.
      ! If it is not then cycle, i.e. we reject this location and
      ! pick another one.
      if (field(3) > 0.0d0) cycle ! Do the next loop iteration, i.e. find a new position.

      ! Calculate the escape probability at the new position, to compair with
      ! the current position.
      df_new = -1.0d0*log(Escape_Prob(field(3), new_pos))

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
        ! Jump to this position with probability alpha, i.e. if rnd is less than alpha
        if (rnd < alpha) then
          cur_pos = new_pos ! New position becomes the current position
          df_cur = df_new
          F_out = field(3)
        end if
      end if
    end do

    ! Return the current position
    Metropolis_Hastings_rectangle_v2 = cur_pos
    df_out = exp(-1.0d0*df_cur)
  end function Metropolis_Hastings_rectangle_v2

  ! ----------------------------------------------------------------------------
  ! Checks the limits of the rectangular region of the emitter
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

      !if(d_x > emitters_dim(1, emit)) then
      !  print *, 'Warning: d_x to large >'
      !  print *, d_x
      !end if
    else if (par_pos(1) < x_min) then
      d_x = x_min - par_pos(1)
      par_pos(1) = d_x + x_min

      !if(d_x > emitters_dim(1, emit)) then
      !  print *, 'Warning: d_x to large <'
      !  print *, d_x
      !end if
    end if

    !Check y ----------------------------------------
    if (par_pos(2) > y_max) then
      d_y = par_pos(2) - y_max
      par_pos(2) = y_max - d_y

      !if(d_y > emitters_dim(2, emit)) then
      !  print *, 'Warning: d_y to large >'
      !  print *, d_y
      !end if
    else if (par_pos(2) < y_min) then
      d_y = y_min - par_pos(2)
      par_pos(2) = d_y + y_min

      !if(d_y > emitters_dim(2, emit)) then
      !  print *, 'Warning: d_x to large <'
      !  print *, d_y
      !end if
    end if
  end subroutine check_limits_metro_rec

  ! ----------------------------------------------------------------------------
  ! Generate random numbers using the Ziggurat method.
  ! Normal distributed random numbers.
  ! This is faster than the Box-Muller.
  function ziggurat_normal(mean, std)
    double precision, intent(in), dimension(1:2) :: mean
    double precision, intent(in), dimension(1:2) :: std
    double precision, dimension(1:2)             :: ziggurat_normal

    ! This stuff is not thread safe.
    ! The Intel compiler barfs at this while using OpenMP.
    ! It has to do with the fact that the save attribute is used in the module.
    ! Variables with the save attribute are shared.
    !$OMP CRITICAL(ZIGGURAT)
    ziggurat_normal(1) = rnor()
    ziggurat_normal(2) = rnor()
    !$OMP END CRITICAL(ZIGGURAT)

    ziggurat_normal(:) = mean(:) + std(:)*ziggurat_normal(:)
  end function ziggurat_normal

end submodule smod_metro_algo
