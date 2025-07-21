!-------------------------------------------!
! Module for testing/debuging               !
! Kristinn Torfason                         !
! 23.10.20                                  !
! Manual placement of electrons for testing !
!-------------------------------------------!

module mod_manual_emission
    use mod_global
    use mod_pair
    use mod_verlet
    implicit none

    PRIVATE
    PUBLIC :: Init_Manual, Clean_Up_Manual

    ! ----------------------------------------------------------------------------
    ! Variables
    integer, dimension(:), allocatable :: nrEmitted_emitters
    integer                            :: nrElecEmitAll

    ! Constant used in MC integration
    double precision :: time_step_div_q0
  contains
    !-----------------------------------------------------------------------------
    ! Initialize the Field Emission
  subroutine Init_Manual()
    ! Allocate the number of emitters
    nrEmit = 1
    allocate(nrEmitted_emitters(1:nrEmit))

    ! Initialize variables
    nrEmitted_emitters = 0 ! Number of electrons emitted from emitter number i in this time step

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_Planar

    ! Function for the electric field in the system
    ptr_field_E => field_E_planar

    ! The function that does the emission
    ptr_Do_Emission => Do_Manual_Emission

    ! The function to do image charge effects
    ptr_Image_Charge_effect => Force_Image_charges_v2

    ! Parameters used in the module
    time_step_div_q0 = time_step / q_0

    ! Initialize the Ziggurat algorithm
    !call zigset(my_seed(1))

  end subroutine Init_Manual

  subroutine Clean_Up_Manual()
    deallocate(nrEmitted_emitters)
  end subroutine Clean_Up_Manual

    !-----------------------------------------------------------------------------
    ! This subroutine gets called from main when the emitters should emit the electrons
  ! subroutine Do_Manual_Emission(step)
  !   integer, intent(in)              :: step
  !   integer                          :: i, u, v, IFAIL = 0, emit_step, species
  !   double precision                 :: cur_time, x, y, z
  !   integer                          :: ud_file
  !   double precision, dimension(1:3) :: par_pos, par_vel, per_pos
  !   integer                          :: Manual_Num_per = 0

  !   nrElecEmitAll = 0
  !   nrEmitted_emitters = 0

  !   ! open(newunit=ud_file, iostat=IFAIL, file='manual_pos.bin', status='OLD', action='READ', access='STREAM')

  !   ! do
  !   !   read(unit=ud_file, iostat=IFAIL) x, y, z, emit_step, species
  !   !   if (IFAIL == 0) then
  !   !     if (emit_step == step) then
  !   !       par_pos(1) = x
  !   !       par_pos(2) = y
  !   !       par_pos(3) = z
  !   !       par_vel = 0.0d0

  !   !       ! Add a particle to the system
  !   !       call Add_Particle(par_pos, par_vel, species, step, 1, -1, 1)
  !   !     end if
  !   !   else
  !   !     exit
  !   !   end if
  !   ! end do

  !   ! close(unit=ud_file, status='keep')

  !   par_vel = 0.0d0
  !   species = species_elec

  !   par_pos(1) = 89.0d0*length_scale
  !   par_pos(2) = 77.3d0*length_scale
  !   par_pos(3) = 3.0d0*length_scale

  !   !call Add_Particle(par_pos, par_vel, species, step, 1, -1, 1)
  !   !print *, 'Adding particle at ', par_pos

  !   do v = -1*Manual_Num_per, Manual_Num_per, 1 ! x
  !       do u = -1*Manual_Num_per, Manual_Num_per, 1 ! y

  !         ! Shift the position
  !         per_pos(1) = par_pos(1) + v*(box_dim(1) + 0.0d0)
  !         per_pos(2) = par_pos(2) + u*(box_dim(2) + 0.0d0)
  !         per_pos(3) = par_pos(3)

  !         call Add_Particle(per_pos, par_vel, species, step, 1, -1, 1)
  !         print *, 'Adding particle at ', per_pos/length_scale
  !         print *, 'u = ', u
  !         print *, 'v = ', v
  !         print *, 'nrPart = ', nrPart
  !         print *, ''
  !       end do
  !   end do

  !   cur_time = time_step * step / time_scale ! Scaled in units of time_scale
  !   write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, nrElecEmitAll, nrElec, &
  !                                                      & (nrEmitted_emitters(i), i = 1, nrEmit)
  ! end subroutine Do_Manual_Emission

  subroutine Do_Manual_Emission(step)
    integer, intent(in) :: step
    integer :: i
    double precision, dimension(1:3) :: par_pos, par_vel
    double precision :: lim_z
    logical :: do_emit

    lim_z = 0.0d0
    do_emit = .true.

    do i = 1, nrPart
      if (particles_cur_pos(3,i) < lim_z) then
        do_emit = .false.
        exit
      end if
    end do

    if (do_emit .eqv. .true.) then
      par_pos = (/0.0d0, 0.0d0, 1.0d0/)
      par_vel = (/0.0d0, 0.0d0, 1.0d0/)
      call Add_Particle(par_pos, par_vel, species_elec, step, 1, -1, 1)
    end if

  end subroutine Do_Manual_Emission

end module mod_manual_emission