!-------------------------------------------!
! Submodule for photo emission              !
! Laser variables                           !
! Hákon Örn Árnason                         !
! 20.08.21                                  !
!-------------------------------------------!

submodule (mod_photo_emission) smod_laser
  use mod_global
  use mod_velocity
  use mod_photo_emission

  implicit none
contains
subroutine Read_laser_function()
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
    case (WORK_VORONOI)
      print '(a)', 'Vacuum: Using Voronoi work function model'
      ptr_Work_fun => w_theta_voronoi

      ! Read the number of sites
      read(unit=ud_work, FMT=*) num_vor_sites

      ! Allocate variables
      allocate(vor_sites_x(1:num_vor_sites))
      allocate(vor_sites_y(1:num_vor_sites))
      allocate(vor_w_theta(1:num_vor_sites))
      allocate(vor_sec(1:num_vor_sites))

      ! Loop over number of sites and read in each one
      do i = 1, num_vor_sites
        read(unit=ud_work, FMT=*) vor_sites_x(i), vor_sites_y(i), vor_w_theta(i), vor_sec(i)

        ! Coordinates in the file should be given on scale from [0, 1]
        vor_sites_x(i) = emitters_pos(1, 1) + vor_sites_x(i)*emitters_dim(1, 1)
        vor_sites_y(i) = emitters_pos(2, 1) + vor_sites_y(i)*emitters_dim(2, 1)
      end do
    case DEFAULT
      print '(a)', 'Vacuum: ERROR UNKNOWN WORK FUNCTION TYPE'
      print *, WORK_TYPE
      stop
    END SELECT

    close(unit=ud_work)
  end subroutine Read_work_function