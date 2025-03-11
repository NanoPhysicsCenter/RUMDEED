!--------------------------------------------!
! Module for emission from a torus           !
! Kristinn Torfason                          !
! 11.03.25                                   !
!--------------------------------------------!

Module mod_torus
    use mod_global
#if defined(_OPENMP)
    use omp_lib
#endif
    use mod_verlet
    use mod_pair
    !use mod_work_function
    use kdtree2_precision_module
    use kdtree2_module
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan

    implicit none
    PRIVATE
    PUBLIC :: Init_Torus, Clean_Up_Torus

    double precision :: time_step_div_q0

    ! ----------------------------------------------------------------------------
    ! Variables
    integer, dimension(:), allocatable :: nrEmitted_emitters
contains
!-------------------------------------------!
! Initialize the cylinder emission
! TODO: Change this for the new system
subroutine Init_Torus()
    ! Local variables
    !double precision :: int_res
    !double precision, dimension(1:3) :: E_test, pos_test, vel_test
    !integer :: IFAIL

    ! Allocate the number of emitters
    allocate(nrEmitted_emitters(1:nrEmit))

    ! Function that checks the boundary conditions for the System
    ptr_Check_Boundary => Check_Boundary_Torus

    ! Function for the electric field in the system
    ptr_field_E => field_E_Torus

    ! The function that does the emission
    ptr_Do_Emission => Do_Field_Emission_Torus_simple

    ! The function to do image charge effects
    ptr_Image_Charge_effect => Image_Charge_Torus ! Force_Image_charges for the cylinder
    !ptr_Image_Charge_effect => Force_Image_charges_v2

    !call Read_work_function()

    ! Parameters used in the module
    time_step_div_q0 = time_step / q_0

    ! Create the KD-tree
    !call Create_KD_Tree()
end subroutine Init_Torus

!-------------------------------------------!
! Clean up the cylinder emission
subroutine Clean_Up_Torus()
    deallocate(nrEmitted_emitters)

    !close(unit=ud_cyl_debug)
end subroutine Clean_Up_Torus

subroutine Check_Boundary_Torus(i)
    integer, intent(in) :: i
    integer             :: sec

    ! To do
    ! Check the boundary conditions for the particle
    !if (Check_Boundary_Cylinder_pos(particles_cur_pos(:, i), sec) .eqv. .true.) then
      ! Remove the particle from the system
    !  call Mark_Particles_Remove(i, sec)
    !end if
  
end subroutine Check_Boundary_Torus

function field_E_Torus(pos) result(field_E)
    ! Input variables
    double precision, dimension(1:3), intent(in) :: pos ! Position to calculate the electric field at

    ! Output variables
    double precision, dimension(1:3)             :: field_E ! Electric field at the point

    ! ! Local variables
    ! integer, parameter                        :: nn = 4 ! number of nearest neighbors to find
    ! integer                                   :: k ! loop variable
    ! type(kdtree2_result)                      :: kd_results(1:nn) ! results from the kd-tree
    ! double precision, dimension(1:nn)         :: w ! weights for the interpolation
    ! double precision, parameter               :: eps = length_scale**3 ! Small number to avoid division by zero
    ! integer, parameter                        :: p = 1 ! Power for the inverse distance weighting

    ! !! For debugging start with planar field
    ! !field_E = field_E_planar(pos)

    ! ! Find the nearest neighbors
    ! !print *, 'pos = ', pos
    ! ! Note: The kd-tree is not thread safe, so we have to use a critical section
    ! !$omp critical (kdtree2)
    ! call kdtree2_n_nearest(tp=kd_tree, qv=pos, nn=nn, results=kd_results)
    ! !$omp end critical (kdtree2)

    ! ! Interpolate the electric field at the point using the nearest neighbors and the inverse distance to each as weights
    ! ! Calculate the weights as 1/distance^p
    ! !print *, 'kd_results(1:nn)%idx = ', kd_results(1:nn)%idx
    ! !print *, 'kd_results(1:nn)%dis = ', kd_results(1:nn)%dis
    ! w(1:nn) = 1.0d0/(kd_results(1:nn)%dis**p + eps)

    ! ! Normalize the weights
    ! w = w / sum(w)
    
    ! ! Calculate the interpolated electric field
    ! field_E = 0.0d0
    ! do k = 1, nn
    !     field_E = field_E + w(k) * kd_data(:, kd_results(k)%idx)
    ! end do
    ! !field_E = field_E * 0.125d1 ! Debug to make field larger
end function field_E_Torus

subroutine Do_Field_Emission_Torus_simple(step)
    integer, intent(in) :: step
    integer             :: i, IFAIL, per
    integer, parameter  :: emit = 1
  
    ! ! Integration
    ! double precision                 :: N_sup
    ! integer                          :: N_round
  
    ! ! Emission variables
    ! double precision                 :: D_f, Df_avg, F_norm, rnd
    ! double precision, dimension(1:3) :: F, n_vec
    ! integer                          :: s, sec, sec_b
    ! integer                          :: nrElecEmit
    ! double precision, dimension(1:3) :: par_pos, par_vel, old_pos
    ! !logical                          :: to_pause = .false.
  
    ! nrEmitted_emitters = 0
  
    ! ! Check if the step is 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90% or 100% of the total steps
    ! if (mod(step, steps/10) == 0) then
    !   !per = nint(dble(step)/dble(steps)*100.0d0) ! Calculate the percentage of the simulation
    !   !print *, 'Step = ', step
    !   !print *, 'Percentage = ', per, '%'
  
    !   ! Calculate the field at the surface
    !   !print *, 'Calling Calc_E_edge_cyl'
    !   !call Calc_E_edge_cyl(per)
    !   !print *, 'Calling Calc_E_side_cyl'
    !   !call Calc_E_top_cyl(per)
    !   !print *, 'Calling Calc_E_top_cyl'
    !   !call Calc_E_corner_cyl(per)
  
    !   call Calc_E_cyl(per)
    !   call Calc_E_circle_cyl(per)
    ! end if
  
    ! call Do_Surface_Integration_simple(N_sup)
    ! N_round = Rand_Poission(N_sup)
  
    ! nrElecEmit = 0
    ! nrEmitted_emitters(emit) = 0
  
    ! Df_avg = 0.0d0
    ! i = 0
    ! ! Loop over the electrons to be emitted.
    ! do s = 1, N_round
    !     ! Get the position of the particle
    !     ! ndim_in, F_out, F_norm_out, pos_xyz_out, sec_out
    !     call Metropolis_Hastings_cyl_tip(N_MH_step, F, F_norm, par_pos, sec, n_vec)
  
    !     ! Check if the field is favorable for emission or not
    !     if (F_norm >= 0.0d0) then
    !       D_f = -huge(1.0d0)
    !       print *, 'Warning: F > 0.0d0 (Do_Field_Emission_Cylinder)'
    !       print *, 'F_norm = ', F_norm
    !       print *, 'F = ', F
    !       print *, 'par_pos = ', par_pos
    !       print *, 'sec = ', sec
    !       stop
    !     else
  
    !       par_pos = par_pos + n_vec*length_scale
    !       call Add_Particle(par_pos, par_vel, species_elec, step, emit, -1, sec)
  
    !       nrElecEmit = nrElecEmit + 1
    !       nrEmitted_emitters(emit) = nrEmitted_emitters(emit) + 1
    !     end if
    !   end do
  
    !   ! Df_avg = Df_avg / dble(i)
    !   ! !print *, 'df_avg = ', Df_avg
    !   ! if (Df_avg > 1.0d-4) then
    !   !   print *, 'RUMDEED: Df_avg > 1.0d-4 (Do_Field_Emission_Cylinder)'
    !   !   print *, step
    !   !   write_position_file = .true.
    !   !   call Write_Position(0)
    !   !   cought_stop_signal = .true.
    !   ! end if
  
    ! write (ud_emit, "(E14.6, *(tr8, i6))", iostat=IFAIL) cur_time, step, nrElecEmit, nrElec, &
    !                                                    & (nrEmitted_emitters(i), i = 1, nrEmit)
end subroutine Do_Field_Emission_Torus_simple

!-------------------------------------------!
! Image charge effect
function Image_Charge_Torus(pos_1, pos_2)
    double precision, dimension(1:3)             :: Image_Charge_Torus ! Image charge effect
    double precision, intent(in), dimension(1:3) :: pos_1 ! Position of the particle we are calculating the force/acceleration on
    double precision, intent(in), dimension(1:3) :: pos_2 ! Position of the particle that is acting on the particle at pos_1

    Image_Charge_Torus = 0.0d0
end function Image_Charge_Torus

end module mod_torus