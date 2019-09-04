!-------------------------------------------!
! module for position dependent          !
! work function                             !
! Kristinn Torfason                         !
! 21.10.18                                  !
!-------------------------------------------!
module mod_work_function
  use mod_global

  implicit none
  integer                                                 :: y_num, x_num, WORK_TYPE

  ! Checkerboard
  double precision, allocatable, dimension(:, :)          :: w_theta_arr ! 1:y_num, 1:x_num
  double precision                                        :: x_len, y_len

  ! Gaussians
  integer                                     :: num_gauss ! Number of Gaussian points
  double precision                            :: w_theta_base
  double precision, allocatable, dimension(:) :: w_gaussians_A ! Amplitude
  double precision, allocatable, dimension(:) :: w_gaussians_x ! x - center
  double precision, allocatable, dimension(:) :: w_gaussians_y ! y - center
  double precision, allocatable, dimension(:) :: w_gaussians_std_x ! standard deviation / spread in x
  double precision, allocatable, dimension(:) :: w_gaussians_std_y ! standard deviation / spread in y

  ! Type of work function models
  integer, parameter :: WORK_CHECKBOARD = 1
  integer, parameter :: WORK_GAUSS      = 2
  integer, parameter :: WORK_CIRCLE     = 3

  interface
    double precision function Work_fun(pos, sec)
      double precision, dimension(1:3), intent(in) :: pos
      integer, intent(out), optional               :: sec
    end function Work_fun
  end interface
  procedure(Work_fun), pointer :: ptr_Work_fun => null()
contains

  subroutine Read_work_function()
    integer :: ud_work, IFAIL, i, j
    character(256) :: iomsg

    ! Open the file that contains information about the work function
    open(newunit=ud_work, iostat=IFAIL, iomsg=iomsg, file='w_theta', &
       & status='OLD', form='FORMATTED', access='SEQUENTIAL', action='READ')
    if (IFAIL /= 0) then
      print *, 'Vacuum: Failed to open file w_theta. ABORTING'
      print *, IFAIL
      print *, iomsg
      stop
    end if

    ! Read the type of work function to use
    read(unit=ud_work, FMT=*) WORK_TYPE

    SELECT CASE (WORK_TYPE)
    case (WORK_CHECKBOARD)
      ! Checkerboard work function
      print '(a)', 'Vacuum: Using checkerboard work function model'
      ptr_Work_fun => w_theta_checkerboard

      ! Read the size of matrix from the file
      read(unit=ud_work, FMT=*) y_num, x_num

      ! Length of each section
      x_len = 1.0d0/x_num
      y_len = 1.0d0/y_num

      ! Allocate the matrix
      allocate(w_theta_arr(1:y_num, 1:x_num))

      ! Read the matrix from the file
      do i = 1, y_num
        read(unit=ud_work, FMT=*) (w_theta_arr(i, j), j = 1, x_num)
      end do
    case (WORK_GAUSS)
      ! Gaussian work function
      print '(a)', 'Vacuum: Using Gaussian work function model'
      ptr_Work_fun => w_theta_gaussian
      ! Read base work function
      ! Read number of gaussian points
      ! Read them one by one, x, y, height, width

      ! Read the base work function
      ! See Gauss function for details
      read(unit=ud_work, FMT=*) w_theta_base

      ! Read the number of gaussian points
      read(unit=ud_work, FMT=*) num_gauss

      ! Allocate the storage array's
      allocate(w_gaussians_A(1:num_gauss))
      allocate(w_gaussians_x(1:num_gauss))
      allocate(w_gaussians_y(1:num_gauss))
      allocate(w_gaussians_std_x(1:num_gauss))
      allocate(w_gaussians_std_y(1:num_gauss))

      ! Read the points in one by one
      ! Each line should have
      ! A, x_c, y_c, std_x, std_y
      do i = 1, num_gauss
        read(unit=ud_work, FMT=*) w_gaussians_A(i), &
                                & w_gaussians_x(i), w_gaussians_y(i), &
                                & w_gaussians_std_x(i), w_gaussians_std_y(i)
      end do
    case (WORK_CIRCLE)
      print '(a)', 'Vacuum: Using circle work function model'
      ptr_Work_fun => w_theta_circle
    case DEFAULT
      print '(a)', 'Vacuum: ERROR UNKNOWN WORK FUNCTION TYPE'
      print *, WORK_TYPE
      stop
    END SELECT

    close(unit=ud_work)
  end subroutine Read_work_function

  ! Clean up stuff
  subroutine Work_fun_cleanup()
    if (WORK_TYPE == WORK_CHECKBOARD) then
      deallocate(w_theta_arr)
    end if

    if (WORK_TYPE == WORK_GAUSS) then
      deallocate(w_gaussians_A)
      deallocate(w_gaussians_x)
      deallocate(w_gaussians_y)
      deallocate(w_gaussians_std_x)
      deallocate(w_gaussians_std_y)
    end if
  end subroutine Work_fun_cleanup

  ! ----------------------------------------------------------------------------
  ! Function that returns the position dependant work function.
  ! This function simply calls the function that was set in Read_work_function.
  double precision function w_theta_xy(pos, sec)
    double precision, dimension(1:3), intent(in) :: pos
    integer, intent(out), optional               :: sec

    ! Call the function set in Read_work_function
    w_theta_xy = ptr_Work_fun(pos, sec)
   
    !w_theta_xy = w_theta_triangle(pos, sec)
    !w_theta_xy = w_theta_checkerboard(pos, sec)
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
      sec = 1
    end if
  end function w_theta_constant

  ! ----------------------------------------------------------------------------
  ! A Gaussian work function that dips
  !
  double precision function w_theta_gaussian(pos, sec)
    double precision, intent(in), dimension(1:3) :: pos
    integer, intent(out), optional               :: sec
    double precision                             :: A
    double precision                             :: std_x, std_y
    double precision                             :: x_c, y_c
    double precision                             :: x, y
    integer                                      :: i

    x = pos(1)
    y = pos(2)

    ! Set the base value
    ! The formula is
    ! w = base + \sum A*exp( -((x-x_c)**2/(2*std_x**2) + (y-y_c)**2/(2*std_y**2)) )
    w_theta_gaussian = w_theta_base

    ! Loop over all the Gaussians
    do i = 1, num_gauss
      ! A should be in eV
      A = w_gaussians_A(i)

      ! We read in x_c and y_c in nm, convert to m
      x_c = w_gaussians_x(i) * length_scale
      y_c = w_gaussians_y(i) * length_scale

      ! We read std_x and std_y in nm, convert to m
      std_x = w_gaussians_std_x(i) * length_scale
      std_y = w_gaussians_std_y(i) * length_scale

      ! Sum upp all the Gaussians
      w_theta_gaussian = w_theta_gaussian &
                     & + A*exp(-1.0d0*( (x-x_c)**2/(2.0d0*std_x**2) + (y-y_c)**2/(2.0d0*std_y**2) ) )
    end do

    ! Section?
    if (present(sec) .eqv. .true.) then
      sec = 1
    end if
  end function w_theta_gaussian

  double precision function w_theta_circle(pos, sec)
  double precision, intent(in), dimension(1:3) :: pos
  integer, intent(out), optional               :: sec

  integer                                      :: sec_
  double precision                             :: r
  double precision, parameter                  :: r_1 = 150.0d0, r_2 = 333.0d0, r_3 = 600.0d0

  r = sqrt(pos(1)**2 + pos(2)**2) / length_scale


  if (r <= r_1) then
    w_theta_circle = 2.20d0
    sec_ = 1
  else if (r <= r_2) then
    w_theta_circle = 1.60d0
    sec_ = 2
  else if (r <= r_3) then
    w_theta_circle = 2.20d0
    sec_ = 3
  else
    w_theta_circle = 2.20d0
    sec_ = 4
  end if

  if (present(sec) .eqv. .true.) then
    sec = sec_
  end if
  end function w_theta_circle

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
        sec = 1
      end if
    else
      w_theta_triangle = 4.5d0
      if (present(sec) .eqv. .true.) then
        sec = 2
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
    double precision, dimension(1:3)              :: pos_scaled
    double precision                              :: x, y
    integer                                       :: x_i, y_i
    integer, parameter                            :: emit = 1 ! Assume emitter nr. 1 for now

    ! To do: Read this from a file
    !w_theta_arr(1, 1:4) = (/ 4.70d0, 4.70d0, 4.70d0, 4.70d0 /)
    !w_theta_arr(2, 1:4) = (/ 4.70d0, 4.60d0, 4.60d0, 4.70d0 /)
    !w_theta_arr(3, 1:4) = (/ 4.70d0, 4.60d0, 4.60d0, 4.70d0 /)
    !w_theta_arr(4, 1:4) = (/ 4.70d0, 4.70d0, 4.70d0, 4.70d0 /)

    !w_theta_arr(1, 1:4) = (/ 4.65d0, 4.65d0, 4.65d0, 4.65d0 /)
    !w_theta_arr(2, 1:4) = (/ 4.65d0, 4.70d0, 4.70d0, 4.65d0 /)
    !w_theta_arr(3, 1:4) = (/ 4.65d0, 4.70d0, 4.70d0, 4.65d0 /)
    !w_theta_arr(4, 1:4) = (/ 4.65d0, 4.65d0, 4.65d0, 4.65d0 /)

    !w_theta_arr(1, 1:4) = (/ 4.65d0, 4.65d0, 4.65d0, 4.65d0 /)
    !w_theta_arr(2, 1:4) = (/ 4.65d0, 4.65d0, 4.65d0, 4.65d0 /)
    !w_theta_arr(3, 1:4) = (/ 4.65d0, 4.65d0, 4.65d0, 4.65d0 /)
    !w_theta_arr(4, 1:4) = (/ 4.65d0, 4.65d0, 4.65d0, 4.65d0 /)

    !w_theta_arr(1, 1:4) = (/ 4.70d0, 4.70d0, 4.70d0, 4.70d0 /)
    !w_theta_arr(2, 1:4) = (/ 4.70d0, 4.70d0, 4.70d0, 4.70d0 /)
    !w_theta_arr(3, 1:4) = (/ 4.70d0, 4.70d0, 4.70d0, 4.70d0 /)
    !w_theta_arr(4, 1:4) = (/ 4.70d0, 4.70d0, 4.70d0, 4.70d0 /)

    ! Scale x, y to unit square
    pos_scaled(1:2) = (pos(1:2) - emitters_pos(1:2, emit)) / emitters_dim(1:2, emit)

    x = pos_scaled(1)
    y = pos_scaled(2)

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
    pos_scaled(1:2) = (pos(1:2) - emitters_pos(1:2, emit)) / emitters_dim(1:2, emit)

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
      sec = 1
    end if
  end function w_theta_checkerboard_2x2

end module mod_work_function
