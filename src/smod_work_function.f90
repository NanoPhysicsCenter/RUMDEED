!-------------------------------------------!
! Submodule for position dependent          !
! work function                             !
! Kristinn Torfason                         !
! 21.10.18                                  !
!-------------------------------------------!
submodule (mod_field_emission_v2) smod_work_function

contains

  ! ----------------------------------------------------------------------------
  ! Function that returns the position dependant work function.
  ! Here we pick the method to be used.
  double precision function w_theta_xy(pos, sec)
    double precision, dimension(1:3), intent(in) :: pos
    integer, intent(out), optional               :: sec

   
    !w_theta_xy = w_theta_triangle(pos, sec)
    w_theta_xy = w_theta_checkerboard(pos, sec)
    !w_theta_xy = w_theta_checkerboard_2x2(pos, sec)
    !w_theta_xy = w_theta_constant(pos, sec)
    !w_theta_xy = w_theta_gaussian(pos, sec)

  end function w_theta_xy

  ! ----------------------------------------------------------------------------
  ! Constant work function
  !
  double precision function w_theta_constant(pos, sec)
    double precision, dimension(1:3), intent(in) :: pos
    integer, intent(out), optional               :: sec

    w_theta_constant = 2.0d0

    if (present(sec) .eqv. .true.) then
      sec = 0
    end if
  end function w_theta_constant

  ! ----------------------------------------------------------------------------
  ! A Gaussian work function that dips
  !
  double precision function w_theta_gaussian(pos, sec)
    double precision, intent(in), dimension(1:3) :: pos
    integer, intent(out), optional               :: sec
    double precision, parameter                  :: A = -0.20d0, B = 25.0d0*length_scale
    double precision, parameter                  :: x_c = 50.0d0*length_scale, y_c = 50.0d0*length_scale
    double precision                             :: x, y

    x = pos(1)
    y = pos(2)

    w_theta_gaussian = 4.70d0 + A*exp( -1.0d0*( (x-x_c)**2/(2*B**2) + (y-y_c)**2/(2*B**2) ) )
    if (present(sec) .eqv. .true.) then
      sec = 0
    end if
  end function w_theta_gaussian

  ! ----------------------------------------------------------------------------
  ! Triangle work function that splits the emitter in two
  ! |----|
  ! |4.7/|
  ! |0 / |
  ! | / 1|
  ! |/4.5|
  ! |----|
  !
  double precision function w_theta_triangle(pos, sec)
    double precision, intent(in), dimension(1:3) :: pos
    integer, intent(out), optional               :: sec
    double precision                             :: x, y

    x = pos(1)
    y = pos(2)

    if (x > y) then
      w_theta_triangle = 4.7d0 !+ A*exp( -1.0d0*( (x-x_c)**2 + (y-y_c)**2 )/(2.0d0*B**2) )
      if (present(sec) .eqv. .true.) then
        sec = 0
      end if
    else
      w_theta_triangle = 4.5d0
      if (present(sec) .eqv. .true.) then
        sec = 1
      end if
    end if
  end function w_theta_triangle

  ! ----------------------------------------------------------------------------
  ! Checkerboard work function
  ! Takes an array and maps it to the emitter area
  ! An array like this
  ! |----|----|----|
  ! | 1  | 2  | 3  |
  ! |----|----|----|
  ! | 4  | 5  | 6  |
  ! |----|----|----|
  ! would map exactly like this to the emitter area
  !
  double precision function w_theta_checkerboard(pos, sec)
    double precision, intent(in), dimension(1:3)  :: pos
    integer, intent(out), optional                :: sec
    double precision,             dimension(1:3)  :: pos_scaled
    double precision                              :: x, y
    double precision                              :: x_len, y_len
    integer                                       :: x_i, y_i
    integer, parameter                            :: emit = 1 ! Assume emitter nr. 1 for now
    integer, parameter                            :: y_num = 4, x_num = 4
    double precision, dimension(1:y_num, 1:x_num) :: w_theta_arr

    ! To do: Read this from a file
    !w_theta_arr(1, 1:4) = (/ 4.70d0, 4.70d0, 4.70d0, 4.70d0 /)
    !w_theta_arr(2, 1:4) = (/ 4.70d0, 4.65d0, 4.65d0, 4.70d0 /)
    !w_theta_arr(3, 1:4) = (/ 4.70d0, 4.65d0, 4.65d0, 4.70d0 /)
    !w_theta_arr(4, 1:4) = (/ 4.70d0, 4.70d0, 4.70d0, 4.70d0 /)

    w_theta_arr(1, 1:4) = (/ 4.65d0, 4.65d0, 4.65d0, 4.50d0 /)
    w_theta_arr(2, 1:4) = (/ 4.65d0, 4.70d0, 4.70d0, 4.65d0 /)
    w_theta_arr(3, 1:4) = (/ 4.65d0, 4.70d0, 4.70d0, 4.65d0 /)
    w_theta_arr(4, 1:4) = (/ 4.65d0, 4.65d0, 4.65d0, 4.65d0 /)

    ! Scale x, y to unit square
    pos_scaled(:) = (pos(:) - emitters_pos(1:3, emit)) / emitters_dim(:, emit)
    x = pos_scaled(1)
    y = pos_scaled(2)

    ! Length of each section
    x_len = 1.0/x_num
    y_len = 1.0/y_num

    ! Calculate the position in the matrix
    x_i = floor(x/x_len) + 1
    y_i = floor(y/y_len) + 1

    ! Check the numbers
    if (x_i > x_num) then
      x_i = x_num
    else if (x_i < 1) then
      x_i = 1
    end if
    if (y_i > y_num) then
      y_i = y_num
    else if (y_i < 1) then
      y_i = 1
    end if

    ! Return the section on the emitter
    ! The numbering scheme is,
    ! |----|----|----|
    ! | 1  | 2  | 3  |
    ! |----|----|----|
    ! | 4  | 5  | 6  |
    ! |----|----|----|
    ! and so forth.
    if (present(sec) .eqv. .true.) then
      sec = x_num*(y_i - 1) + x_i
    end if

    ! Reverse the y-direction in the array
    y_i = y_num - y_i + 1

    ! Return the results
    w_theta_checkerboard = w_theta_arr(y_i, x_i)


  end function w_theta_checkerboard

  double precision function w_theta_checkerboard_2x2(pos, sec)
    double precision, intent(in), dimension(1:3) :: pos
    integer, intent(out), optional               :: sec
    double precision,             dimension(1:3) :: pos_scaled
    double precision                             :: x, y
    integer, parameter                           :: emit = 1

    ! Scale x, y to unit square
    pos_scaled(:) = (pos(:) - emitters_pos(1:3, emit)) / emitters_dim(:, emit)
    x = pos_scaled(1)
    y = pos_scaled(2)

    !
    ! |---|---|
    ! | 1 | 2 |
    ! |---|---|
    ! | 3 | 4 |
    ! |---|---|
    !
    ! Use unit square coordinates
    if (x > 0.5d0) then
      if (y > 0.5d0) then
        ! 2
        w_theta_checkerboard_2x2 = 4.75d0
      else
        ! 4
        w_theta_checkerboard_2x2 = 4.68d0
      end if
    else
      if (y > 0.5d0) then
        ! 1
        w_theta_checkerboard_2x2 = 4.65d0
      else
        ! 3
        w_theta_checkerboard_2x2 = 4.72d0
      end if
    end if

    if (present(sec) .eqv. .true.) then
      sec = 0
    end if
  end function w_theta_checkerboard_2x2

end submodule smod_work_function