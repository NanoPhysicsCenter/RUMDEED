!-------------------------------------------!
! Submodule for Monte Carlo Integration     !
! routines for field emission               !
! Kristinn Torfason                         !
! 21.10.18                                  !
!-------------------------------------------!
module mod_collisions
    use mod_global
    use mod_pair

    implicit none

    double precision, dimension(:), pointer                :: cross_tot_energy => null()
    double precision, dimension(:), pointer                :: cross_tot_data => null()
    double precision, dimension(:), pointer                :: cross_ion_energy => null()
    double precision, dimension(:), pointer                :: cross_ion_data => null()

    !double precision :: cross_sum

contains
  ! ----------------------------------------------------------------------------
  ! Ion collisions, Mean free path approach.
  ! Simulate collisions with N2 molecules. We read the collision cross section and use it
  ! to calculate the mean free path of the electrons. The distance the electrons travel
  ! in one time step divided with the mean free path is used as a probability to
  ! determine if a collision occurs. The density of N2 molecules is calulated from
  ! the temperature and pressure.
  subroutine Do_Ion_Collisions(step)
    integer, intent(in)              :: step
    double precision, dimension(1:3) :: cur_pos, prev_pos, par_vec, old_vel
    double precision, parameter      :: v2_min      = (2.0d0*q_0*0.1d0/m_0) ! Minimum velocity squared
    double precision, parameter      :: v2_max      = (2.0d0*q_0*5000.0d0/m_0) ! Maximum velocity squared
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
    integer                          :: ion_life_time

    nrColl = 0
    nrIon = 0
    count_n = 0
    mean_path_avg = 0.0d0
    mean_actual_avg = 0.0d0

    if (nrPart > 0) then

    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP& PRIVATE(i, cur_pos, prev_pos, d, alpha, rnd, par_vec, vel2, KE, mean_path, cross_tot, cross_ion) &
    !$OMP& PRIVATE(old_vel, ion_life_time) &
    !$OMP& SHARED(nrPart, particles_species, particles_life, step, particles_cur_pos) &
    !$OMP& SHARED(particles_prev_pos, particles_cur_vel, time_step, n_d) &
    !$OMP& REDUCTION(+:nrColl, count_n, mean_path_avg, nrIon) SCHEDULE(GUIDED, CHUNK_SIZE)
    do i = 1, nrPart
      if (particles_species(i) /= species_elec) then
        if (step >= particles_life(i)) then
          call Mark_Particles_Remove(i, remove_top) ! Mark the ion to be removed, it has reached the end of its life.
        end if
      else

        cur_pos(:) = particles_cur_pos(:, i)
        prev_pos(:) = particles_prev_pos(:, i)
        !prev_pos(:) = particles_last_col_pos(:, i)
        old_vel = particles_cur_vel(:, i)

        d = sqrt( (cur_pos(1) - prev_pos(1))**2 + (cur_pos(2) - prev_pos(2))**2 + (cur_pos(3) - prev_pos(3))**2 )
        vel2 = old_vel(1)**2 + old_vel(2)**2 + old_vel(3)**2
      
        ! The velocity should be above the minimum velocity we set.
        !if ((vel2 > v2_min) .and. (vel2 < v2_max)) then
        if (vel2 > v2_min) then

          if (vel2 > v2_max) then
            KE = 0.5d0*m_0*v2_max/q_0
          else
            KE = 0.5d0*m_0*vel2/q_0 ! Energy in eV
          end if

          cross_tot = Find_Cross_tot_data(KE)

          mean_path = 1.0d0/(n_d*cross_tot)
          mean_path_avg = mean_path_avg + mean_path
          count_n = count_n + 1

          ! Calculate the ratio between the distance travled and the mean free path
          alpha = d/mean_path
          ! if (alpha > 1.0d0) then
          !   print *, 'WARNING: alpha > 1 in mean path'
          !   print *, mean_path
          !   print *, mean_path/1.0E-9
          !   print *, n_d
          !   print *, cross_tot
          !   print *, KE
          !   print *, sqrt(vel2)
          !   print *, d
          !   !pause
          ! end if

          ! Check if we do a collision or not
          call random_number(rnd)
          if (rnd < alpha) then
            ! Pick a new random direction for the particle
            par_vec = Get_Injected_Vec(KE, old_vel)
            !call random_number(par_vec)
            !par_vec(1:2) = par_vec(1:2) - 0.5d0
            !par_vec(3) = par_vec(3) - 0.25d0
            !par_vec = par_vec - 0.5d0
            par_vec = par_vec / sqrt(par_vec(1)**2 + par_vec(2)**2 + par_vec(3)**2)

            ! Set the new velocity
            call random_number(rnd)

            particles_cur_vel(:, i) = par_vec*sqrt(vel2)*rnd*1.0d0

            call random_number(rnd)
            !par_vec = par_vec*sqrt(vel2)*(1.0d0 - rnd)

            rnd = rnd*0.80d0 + (1.0d0-0.80d0)
            par_vec = par_vec*sqrt(2.0d0*q_0*rnd*40.0d0/m_0) ! Give the electron an initial energy between 8 and 40 eV

            !prev_pos(:) = particles_last_col_pos(:, i)
            !d = sqrt( (cur_pos(1) - prev_pos(1))**2 + (cur_pos(2) - prev_pos(2))**2 + (cur_pos(3) - prev_pos(3))**2 )
            !particles_last_col_pos(:, i) = cur_pos(:)
            !mean_actual_avg = mean_actual_avg + d

            !---------------------------------------------
            ! Check if the collision ionizes the N2 or not  
            cross_ion = Find_Cross_ion_data(KE)
            
            alpha = cross_ion/cross_tot
            !alpha = 0.20d0
            ! if (alpha > 1.0d0) then
            !   print *, 'WARNING: alpha > 1 in cross section'
            ! end if
            
            call random_number(rnd)
            if (rnd < alpha) then ! Check if we ionize or not

              ! We ionize
              ! Pick a position for the new electron to appear at some where close to where the collision occurred
              call random_number(prev_pos)
              prev_pos = prev_pos - 0.5d0
              cur_pos = cur_pos + prev_pos*length_scale

              ! Pick a direction for the new electron to go in
              par_vec = Get_Ejected_Vec(KE, KE, old_vel)

              !$OMP CRITICAL
              ! Add the new electron to the system
              call Add_Particle(cur_pos, par_vec, species_elec, step, 1, -1) ! Electron
              !$OMP END CRITICAL

              par_vec = 0.0d0

              ! Pick a position for the ion to appear at some where close to where the collision occurred
              call random_number(prev_pos)
              prev_pos = prev_pos - 0.5d0
              cur_pos = cur_pos + prev_pos*length_scale

              ! Calculate the life time of the Ion
              ! This number is drawn from an exponential distribution
              ! The half life is 0.13843690559395497 ps
              alpha = -0.13843690559395497E-12 / time_step
              call random_number(rnd)
              ion_life_time = NINT(alpha*log(1.0d0 - rnd))

              !$OMP CRITICAL
              ! Add the new positively charged ion to the system
              call Add_Particle(cur_pos, par_vec, species_hole, step, 1, step+ion_life_time) ! Ion
              !$OMP END CRITICAL

              nrIon = nrIon + 1
            end if

            ! Update the number of collisions
            nrColl = nrColl + 1
          end if
        else
          KE = 0.5d0*m_0*vel2/q_0 ! Energy in eV
          ! if (KE > 0.1d0) then
          !   print *, 'Outside range'
          !   print *, KE
          !   print *, ''
          ! end if
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
  end subroutine Do_Ion_Collisions

  function Get_Injected_Vec(T, par_vel)
    double precision, dimension(1:3)             :: Get_Injected_Vec
    double precision, dimension(1:3), intent(in) :: par_vel
    double precision, intent(in)                 :: T ! Energy in eV
    double precision, parameter                  :: mu = 5.0d0, sigma = 25.0d0
    double precision                             :: angle
    double precision                             :: m_factor, rnd, alpha
    double precision                             :: dot_p, len_vec, len_vel
    double precision, dimension(1:3)             :: par_vec
    integer                                      :: n_tries

    !m_factor = 1.0d0/(sqrt(2.0d0*pi)*sigma)
    m_factor = folded_normal_dist(mu, sigma, mu)
    n_tries = 0

    do
      ! Get random numbers between 0 and 1
      call random_number(par_vec)
      ! Convert this to random numbers between -0.5 and 0.5
      par_vec = par_vec - 0.5d0
      ! Dot product
      dot_p = par_vel(1)*par_vec(1) + par_vel(2)*par_vec(2) + par_vel(3)*par_vec(3)

      ! Length of the vectors
      len_vec = sqrt(par_vec(1)**2 + par_vec(2)**2 + par_vec(3)**2)
      len_vel = sqrt(par_vel(1)**2 + par_vel(2)**2 + par_vel(3)**2)

      ! Calculate the angle between the two vector using the dot product.
      ! We want the angle in degrees.
      angle = acos(dot_p/(len_vec*len_vel)) * 180.0d0/pi

      ! Get the acceptance/rejection criteria. We use a folded normal distribution
      ! because the angle between two vector is always positive. [0, pi]/[0, 180]. 
      alpha = folded_normal_dist(mu, sigma, angle) / m_factor

      call random_number(rnd)
      if (rnd < alpha) then
        exit ! Exit the loop
      else
        n_tries = n_tries + 1
        if (n_tries >= 1000000) then ! This should never take this long
          print *, 'n_tries > 100000 (Injected)'
          print *, 'angle = ', angle
          print *, 'alpha = ', alpha
          print *, 'n_tries = ', n_tries
          print *, ''
          exit
        end if
      end if
    end do

    ! Return the normalized vector
    Get_Injected_Vec = (par_vec / len_vec)
  end function Get_Injected_Vec

  ! This function returns a normalized direction vector for the ejected electron
  ! The angle between this vector and the velocity vector of the incident electron
  ! is approximatly distrubted according the experimental values.
  ! See:
  ! Dobly Differential Cross Section for Electron Scattered by Nitrogen
  ! J. C. Nogueira, M. A. Eschiapati Ferreira and Ronaldo S. Barbieri
  function Get_Ejected_Vec(W, T, par_vel)
    double precision, dimension(1:3)             :: Get_Ejected_Vec
    double precision, intent(in)                 :: T, W
    double precision, dimension(1:3), intent(in) :: par_vel
    double precision, parameter                  :: a = -430.5d0, b = -0.5445d0, c = 89.32d0
    double precision, parameter                  :: sigma = 48.0d0
    double precision                             :: angle_max, angle
    double precision                             :: m_factor, rnd, alpha
    double precision                             :: dot_p, len_vec, len_vel
    double precision, dimension(1:3)             :: par_vec
    integer                                      :: n_tries 

    if (T < 100.d0) then
      angle_max = a*100.0d0**b + c
    else
      angle_max = a*T**b + c
    end if

    !m_factor = 1.0d0/(sqrt(2.0d0*pi)*sigma)
    m_factor = folded_normal_dist(angle_max, sigma, angle_max)
    n_tries = 0

    do
      call random_number(par_vec)
      par_vec = par_vec - 0.5d0
      dot_p = par_vel(1)*par_vec(1) + par_vel(2)*par_vec(2) + par_vel(3)*par_vec(3)
      len_vec = sqrt(par_vec(1)**2 + par_vec(2)**2 + par_vec(3)**2)
      len_vel = sqrt(par_vel(1)**2 + par_vel(2)**2 + par_vel(3)**2)

      angle = acos(dot_p/(len_vec*len_vel)) * 180.0d0/pi

      alpha = folded_normal_dist(angle_max, sigma, angle) / m_factor

      call random_number(rnd)
      if (rnd < alpha) then
        exit ! Exit the loop
      else
        n_tries = n_tries + 1
        if (n_tries >= 1000000) then ! This should never take this long
          print *, 'n_tries > 100000 (Ejected)'
          exit
        end if
      end if
    end do

    Get_Ejected_Vec = (par_vec / len_vec)
  end function Get_Ejected_Vec

  ! Normal distribution
  double precision function normal_dist(mu, sigma, x)
    double precision, intent(in) :: mu, sigma, x

    normal_dist = 1.0d0/(sqrt(2.0d0*pi)*sigma)*exp(-(x-mu)**2/(2.0d0*sigma**2))
  end function normal_dist

  ! Folded normal distribution
  ! https://en.wikipedia.org/wiki/Folded_normal_distribution
  double precision function folded_normal_dist(mu, sigma, x)
    double precision, intent(in) :: mu, sigma, x
    double precision             :: sigma2

    sigma2 = sigma**2
    folded_normal_dist = sqrt(2.0d0/(pi*sigma2)) * exp(-1.0d0*(mu**2 + x**2)/(2.0d0*sigma2))*cosh(mu*x/sigma2)
  end function folded_normal_dist


    ! --------------------------------------------------------------------------
    ! Read the Cross Section from a file
  subroutine Read_Cross_Section()
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
      i = BinarySearch(cross_ion_energy, energy, i1, i2)
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
end module mod_collisions
