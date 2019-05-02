!-------------------------------------------!
! Submodule for Monte Carlo Integration     !
! routines for field emission               !
! Kristinn Torfason                         !
! 21.10.18                                  !
!-------------------------------------------!
submodule (mod_verlet) smod_collisions
    use mod_global

    implicit none

    double precision, dimension(:), pointer                :: cross_tot_energy => null()
    double precision, dimension(:), pointer                :: cross_tot_data => null()
    double precision, dimension(:), pointer                :: cross_ion_energy => null()
    double precision, dimension(:), pointer                :: cross_ion_data => null()

    !double precision :: cross_sum

contains
! ----------------------------------------------------------------------------
  ! Ion collisions
subroutine Do_Collisions_1()
    !double precision, parameter      :: N_mean_col  = 100 ! Average number of collisions per time step
    integer, parameter               :: N_max_tries = 1000 ! Maximum number of tries before we give up looking for particles
    double precision, parameter      :: v2_min      = (2.0d0*q_0*1.0d0/m_0) ! Minimum velocity squared
    double precision, parameter      :: v2_max      = (2.0d0*q_0*100.0d0/m_0) ! Maximum velocity squared
    integer                          :: N_col           ! Number of collisions to do in this time step
    integer                          :: N_try            ! Number of tries done so far
    double precision                 :: rnd
    double precision                 :: vel2             ! Squared velocity of the current particle
    integer                          :: i
    double precision, dimension(1:3) :: par_vec

  
    ! Number of collisions per time step is poisson distributed
    N_col = Rand_Poission(collisions_mean)

    ! Keep track of what particles have had collisions
    particles_collision = .false.
    
    N_try = 0
    do while ((N_try < N_max_tries) .and. (N_col > 0))
      ! Randomly pick a particle
      call random_number(rnd)
      i = floor(rnd*nrPart) + 1

      ! Calulate the squared velocity of the particle picked
      vel2 = particles_cur_vel(1, i)**2 + particles_cur_vel(2, i)**2 + particles_cur_vel(3, i)**2

      ! Check if it is above and below the minimum and maximum
      if ((vel2 > v2_min) .and. (vel2 < v2_max) .and. (particles_collision(i) .eqv. .false.)) then
        N_col = N_col - 1 ! One less collision to do
        N_try = 0 ! Reset number of failed attempts
        particles_collision(i) = .true. ! Keep track of what particles have had collisions

        ! Pick a new random direction for the particle
        ! Todo: This should not be uniform
        call random_number(par_vec)
        par_vec = par_vec / sqrt(par_vec(1)**2 + par_vec(2)**2 + par_vec(3)**2)

        ! Reduce the energy by a random amount between 0 and 10%
        ! Todo: This should probably not be uniform
        ! Todo: Check upper and lower limit
        call random_number(rnd)
        rnd = rnd*10.0d0/100.0d0
        vel2 = vel2*rnd
        
        ! Set the new velocity
        particles_cur_vel(:, i) = par_vec*sqrt(vel2)
      else
        N_try = N_try + 1 ! Try again
      end if
    end do
  end subroutine


  ! ----------------------------------------------------------------------------
  ! Ion collisions
  ! Here we create an additional electron
  subroutine Do_Collisions_2(step)
    integer, intent(in)              :: step
    !double precision, parameter      :: N_mean_col  = 100 ! Average number of collisions per time step
    integer, parameter               :: N_max_tries = 1000 ! Maximum number of tries before we give up looking for particles
    double precision, parameter      :: v2_min      = (2.0d0*q_0*1.0d0/m_0) ! Minimum velocity squared
    double precision, parameter      :: v2_max      = (2.0d0*q_0*100.0d0/m_0) ! Maximum velocity squared
    integer                          :: N_col           ! Number of collisions to do in this time step
    integer                          :: N_try            ! Number of tries done so far
    double precision                 :: rnd
    double precision                 :: vel2             ! Squared velocity of the current particle
    integer                          :: i
    double precision, dimension(1:3) :: par_vec

  
    ! Number of collisions per time step is poisson distributed
    N_col = Rand_Poission(nrPart*0.05d0)

    ! Keep track of what particles have had collisions
    particles_collision = .false.
    
    N_try = 0
    do while ((N_try < N_max_tries) .and. (N_col > 0))
      ! Randomly pick a particle
      call random_number(rnd)
      i = floor(rnd*nrPart) + 1

      ! Calulate the squared velocity of the particle picked
      !vel2 = particles_cur_vel(1, i)**2 + particles_cur_vel(2, i)**2 + particles_cur_vel(3, i)**2

      ! Check if it is above and below the minimum and maximum
      if ((vel2 > v2_min) .and. (vel2 < v2_max) .and. (particles_collision(i) .eqv. .false.)) then
        N_col = N_col - 1 ! One less collision to do
        N_try = 0 ! Reset number of failed attempts
        particles_collision(i) = .true. ! Keep track of what particles have had collisions

        par_vec = 0.90d0*particles_cur_vel(:, i)
        particles_cur_vel(:, i) = 0.10d0*particles_cur_vel(:, i)

        call Add_Particle(particles_cur_pos(:, i), par_vec, species_elec, step, 1)
      else
        N_try = N_try + 1 ! Try again
      end if
    end do
  end subroutine

  ! ----------------------------------------------------------------------------
  ! Ion collisions
  ! Mean free path approch
  subroutine Do_Collisions_3(step)
    integer, intent(in)              :: step
    double precision, parameter      :: mean_path = 1000.0d0*length_scale ! Mean free path
    double precision, dimension(1:3) :: cur_pos, prev_pos, par_vec
    double precision, parameter      :: v2_min      = (2.0d0*q_0*10.0d0/m_0) ! Minimum velocity squared
    double precision, parameter      :: v2_max      = (2.0d0*q_0*100.0d0/m_0) ! Maximum velocity squared
    double precision                 :: d ! The distance traveled
    double precision                 :: rnd, alpha
    double precision                 :: vel2             ! Squared velocity of the current particle
    double precision, parameter      :: e_max = 0.10d0    ! Max value of the coefficient of restitution
    integer                          :: i, nrColl, IFAIL

    nrColl = 0

    !$OMP PARALLEL DO PRIVATE(i, cur_pos, prev_pos, d, alpha, rnd, par_vec, vel2) &
    !$OMP& REDUCTION(+:nrColl) SCHEDULE(GUIDED, CHUNK_SIZE)
    do i = 1, nrPart
      cur_pos(:) = particles_cur_pos(:, i)
      prev_pos(:) = particles_prev_pos(:, i)

      d = sqrt( (cur_pos(1) - prev_pos(1))**2 + (cur_pos(2) - prev_pos(2))**2 + (cur_pos(3) - prev_pos(3))**2 )
      alpha = d/mean_path

      vel2 = particles_cur_vel(1, i)**2 + particles_cur_vel(2, i)**2 + particles_cur_vel(3, i)**2
    
      ! Check if we do a collision or not
      if ((vel2 > v2_min) .and. (vel2 < v2_max)) then
        ! Check if we do a collision or not
        call random_number(rnd)
        if (rnd < alpha) then
          ! Pick a new random direction for the particle
          ! Todo: This should not be uniform
          call random_number(par_vec)
          par_vec = par_vec / sqrt(par_vec(1)**2 + par_vec(2)**2 + par_vec(3)**2)

          ! Set the new velocity
          !call random_number(rnd)
          !rnd = rnd*e_max
          !vel2 = particles_cur_vel(1, i)**2 + particles_cur_vel(2, i)**2 + particles_cur_vel(3, i)**2
          !vel2 = vel2*rnd
          particles_cur_vel(:, i) = par_vec*sqrt(vel2)*0.1d0
          par_vec = par_vec*sqrt(vel2)*0.9d0
          particles_last_col_pos(:, i) = cur_pos
          
          alpha = 0.5d0
          call random_number(rnd)
          if (rnd < alpha) then
            !$OMP CRITICAL
            call Add_Particle(particles_cur_pos(:, i), par_vec, species_elec, step, 1)
            !$OMP END CRITICAL
          end if

          ! Update the number of collisions
          nrColl = nrColl + 1
        end if
      end if
    end do
    !$OMP END PARALLEL DO

    ! Write data
    write(ud_coll, '(i6, tr2, i6)', iostat=IFAIL) step, nrColl
  end subroutine

  ! ----------------------------------------------------------------------------
  ! Ion collisions
  ! Mean free path approch
  module subroutine Do_Collisions_4(step)
    integer, intent(in)              :: step
    !double precision, parameter      :: mean_path = 1000.0d0*length_scale ! Mean free path
    double precision, dimension(1:3) :: cur_pos, prev_pos, par_vec
    double precision, parameter      :: v2_min      = (2.0d0*q_0*0.1d0/m_0) ! Minimum velocity squared
    double precision, parameter      :: v2_max      = (2.0d0*q_0*5000.0d0/m_0) ! Maximum velocity squared
    double precision, parameter      :: T_temp = 293.15d0 ! Temperature in Kelvin
    double precision, parameter      :: P_abs = 101325.0d0 ! Absolute pressure in Pa
    double precision, parameter      :: n_d = P_abs/(k_b*T_temp) ! Density
    double precision                 :: mean_path, mean_path_avg, mean_actual_avg ! Mean free path
    integer                          :: count_n
    double precision                 :: cross_tot, cross_ion
    double precision                 :: d ! The distance traveled
    double precision                 :: rnd, alpha
    double precision                 :: vel2             ! Squared velocity of the current particle
    double precision, parameter      :: e_max = 0.10d0    ! Max value of the coefficient of restitution
    integer                          :: i, nrColl, IFAIL, nrIon
    double precision                 :: KE ! Kinetic energy
    double precision                 :: cur_time

    nrColl = 0
    nrIon = 0
    count_n = 0
    mean_path_avg = 0.0d0
    mean_actual_avg = 0.0d0

    if (nrElec > 0) then

    !$OMP PARALLEL DO PRIVATE(i, cur_pos, prev_pos, d, alpha, rnd, par_vec, vel2, KE, mean_path, cross_tot, cross_ion) &
    !$OMP& REDUCTION(+:nrColl, count_n, mean_path_avg, nrIon) SCHEDULE(GUIDED, CHUNK_SIZE)
    do i = 1, nrElec
      cur_pos(:) = particles_cur_pos(:, i)
      prev_pos(:) = particles_prev_pos(:, i)
      !prev_pos(:) = particles_last_col_pos(:, i)

      d = sqrt( (cur_pos(1) - prev_pos(1))**2 + (cur_pos(2) - prev_pos(2))**2 + (cur_pos(3) - prev_pos(3))**2 )
      vel2 = particles_cur_vel(1, i)**2 + particles_cur_vel(2, i)**2 + particles_cur_vel(3, i)**2
    
      ! Check if we do a collision or not
      if ((vel2 > v2_min) .and. (vel2 < v2_max)) then

        KE = 0.5d0*m_0*vel2/q_0 ! Energy in eV

        cross_tot = Find_Cross_tot_data(KE)

        mean_path = 1.0d0/(n_d*cross_tot)
        mean_path_avg = mean_path_avg + mean_path
        count_n = count_n + 1

        alpha = d/mean_path
        if (alpha > 1.0d0) then
          print *, 'WARNING: alpha > 1 in mean path'
          print *, mean_path
          print *, mean_path/1.0E-9
          print *, n_d
          print *, cross_tot
          print *, KE
          print *, sqrt(vel2)
          print *, d
          !pause
        end if

        ! Check if we do a collision or not
        call random_number(rnd)
        if (rnd < alpha) then
          ! Pick a new random direction for the particle
          call random_number(par_vec)
          par_vec(1:2) = par_vec(1:2) - 0.5d0
          par_vec(3) = par_vec(3) - 0.25d0
          !par_vec = par_vec - 0.5d0
          par_vec = par_vec / sqrt(par_vec(1)**2 + par_vec(2)**2 + par_vec(3)**2)

          ! Set the new velocity
          call random_number(rnd)

          particles_cur_vel(:, i) = par_vec*sqrt(vel2)*rnd*1.0d0

          call random_number(rnd)
          !par_vec = par_vec*sqrt(vel2)*(1.0d0 - rnd)
          par_vec = par_vec*sqrt(2.0d0*q_0*rnd*40.0d0/m_0) ! Give the electron an initial energy between 0 and 40 eV

          !prev_pos(:) = particles_last_col_pos(:, i)
          !d = sqrt( (cur_pos(1) - prev_pos(1))**2 + (cur_pos(2) - prev_pos(2))**2 + (cur_pos(3) - prev_pos(3))**2 )
          !particles_last_col_pos(:, i) = cur_pos(:)
          !mean_actual_avg = mean_actual_avg + d

          cross_ion = Find_Cross_ion_data(KE)
          
          alpha = cross_ion/cross_tot
          !alpha = 0.20d0
          if (alpha > 1.0d0) then
            print *, 'WARNING: alpha > 1 in cross section'
          end if
          call random_number(rnd)
          if (rnd < alpha) then
            !$OMP CRITICAL
            !!call Add_Particle(particles_cur_pos(:, i), par_vec, species_elec, step, 1)
            call random_number(prev_pos)
            prev_pos = prev_pos - 0.5d0
            cur_pos = cur_pos + prev_pos*length_scale
            !call Add_Particle(cur_pos, step, par_vec)
            call Add_Particle(cur_pos, par_vec, species_elec, step, 1)
            !$OMP END CRITICAL
            nrIon = nrIon + 1
          end if

          ! Update the number of collisions
          nrColl = nrColl + 1
        end if
      else
        KE = 0.5d0*m_0*vel2/q_0 ! Energy in eV
        if (KE > 0.1d0) then
          print *, 'Outside range'
          print *, KE
          print *, ''
        end if
      end if
    end do
    !$OMP END PARALLEL DO

    end if

    mean_path_avg = mean_path_avg / count_n
    mean_actual_avg = 0.0d0

    cur_time = time_step * step / time_scale ! Scaled in units of time_scale

    ! Write data
    write(ud_coll, '(i6, tr2, ES12.4, tr2, i6, tr2, i6, tr2, i6, tr2, ES12.4, tr2, ES12.4)', iostat=IFAIL) &
             step, cur_time, nrColl, nrIon, count_n, &
             (mean_path_avg/length_scale), (mean_actual_avg/length_scale)
  end subroutine Do_Collisions_4


    ! --------------------------------------------------------------------------
    ! Read the Cross Section from a file
  module subroutine Read_Cross_Section()
    integer                      :: IFAIL, ud_cross, i, n
    character (len=*), parameter :: filename_tot="N2-tot-cross.txt"
    character (len=*), parameter :: filename_ion="N2-ion-cross.txt"

    ! Open the file for reading
    open(newunit=ud_cross, iostat=IFAIL, file=filename_tot, status='OLD', action='READ')
    if (IFAIL /= 0) then
      print '(a)', 'Vacuum: ERROR UNABLE TO OPEN file ', filename_tot
      stop
    end if

    ! Count the number of lines in the file
    n = 0
    do
      read(ud_cross, *, iostat=IFAIL)
      if (IFAIL /= 0) exit
      n = n + 1
    end do

    ! Allocate the array to hold the file
    allocate(N2_tot_cross(1:n, 1:2))

    ! Rewind to the start of the file
    rewind(ud_cross)

    ! Read the file into the array
    do i = 1, n
      read(ud_cross, *) N2_tot_cross(i, 1), N2_tot_cross(i, 2)
    end do

    ! Close the file
    close(unit=ud_cross, iostat=IFAIL, status='keep')

    cross_tot_energy => N2_tot_cross(:, 1)
    cross_tot_data => N2_tot_cross(:, 2)

    !------------------------------------------------------------------------------------------

    ! Open the file for reading
    open(newunit=ud_cross, iostat=IFAIL, file=filename_ion, status='OLD', action='READ')
    if (IFAIL /= 0) then
      print '(a)', 'Vacuum: ERROR UNABLE TO OPEN file ', filename_ion
      stop
    end if

    ! Count the number of lines in the file
    n = 0
    do
      read(ud_cross, *, iostat=IFAIL)
      if (IFAIL /= 0) exit
      n = n + 1
    end do

    ! Allocate the array to hold the file
    allocate(N2_ion_cross(1:n, 1:2))

    ! Rewind to the start of the file
    rewind(ud_cross)

    ! Read the file into the array
    do i = 1, n
      read(ud_cross, *) N2_ion_cross(i, 1), N2_ion_cross(i, 2)
    end do

    ! Close the file
    close(unit=ud_cross, iostat=IFAIL, status='keep')

    cross_ion_energy => N2_ion_cross(:, 1)
    cross_ion_data => N2_ion_cross(:, 2)
    !cross_sum  = sum(cross_data)

    !print *, 'abs_coeff', abs_coeff
    !pause
  end subroutine Read_Cross_Section

  double precision function Find_Cross_tot_data(energy)
    double precision, intent(in) :: energy
    double precision, parameter  :: a = 7.98d0, b = -0.005845d0, c = 4.628d0, d = -0.0007864d0
    double precision             :: x1, x2, y1, y2, h, q
    integer                      :: i, i1, i2

    if ((energy > 70.0d0) .and. (energy <= 3000.0d0)) then
      ! Use fitted data
      Find_Cross_tot_data = (a*exp(b*energy) + c*exp(d*energy)) * 1.0d-20
    else
      i = BinarySearch(cross_tot_energy, energy, i1, i2)
      y1 = cross_tot_data(i1)
      y2 = cross_tot_data(i2)
      x1 = cross_tot_energy(i1)
      x2 = cross_tot_energy(i2)

      h = (y1 - y2)/(x1 - x2)
      q = (y2*x1 - y1*x2)/(x1 - x2)

      Find_Cross_tot_data = (h*energy + q)*1.0d-20
      !Find_Cross_tot_data = cross_tot_data(i) * 1.0d-20 ! Data is in 10^-16 cm^2 to m^2
    end if
  end function Find_Cross_tot_data

  double precision function Find_Cross_ion_data(energy)
    double precision, intent(in) :: energy
    double precision, parameter  :: a = 2.251d0, b = -0.00311d0, c = 1.04d0, d = -0.0003378d0
    double precision             :: x1, x2, y1, y2, h, q
    integer                      :: i, i1, i2

    if ((energy > 180.0d0) .and. (energy <= 3000.0d0)) then
      Find_Cross_ion_data = (a*exp(b*energy) + c*exp(d*energy)) * 1.0d-20
    else
      i = BinarySearch(cross_ion_energy, energy)
      y1 = cross_ion_data(i1)
      y2 = cross_ion_data(i2)
      x1 = cross_ion_energy(i1)
      x2 = cross_ion_energy(i2)

      h = (y1 - y2)/(x1 - x2)
      q = (y2*x1 - y1*x2)/(x1 - x2)

      Find_Cross_ion_data = (h*energy + q)*1.0d-20

      !Find_Cross_ion_data = cross_ion_data(i) * 1.0d-20 ! Data is in 10^-16 cm^2 to m^2
    end if
  end function Find_Cross_ion_data
end submodule smod_collisions