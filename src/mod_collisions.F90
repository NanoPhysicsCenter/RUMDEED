!-------------------------------------------!
! Submodule for Monte Carlo Integration     !
! routines for field emission               !
! Kristinn Torfason                         !
! 21.10.18                                  !
!-------------------------------------------!
module mod_collisions
  use mod_global
  use mod_pair
  use mod_polynomialroots

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

subroutine Do_Electron_Atom_Collisions(step)
  integer, intent(in)               ::  step
  integer                           ::  IFAIL, nrCollisions=0, nrRecombinations=0, nrIonizations=0

  if (step >= collision_delay) then
  if (two_time_step .eqv. .true.) then
    call Update_Collision_Data_All_tts()

    select case (collision_mode)
      case (1) ! Continuous ionization
        call Do_Continuous_Ionization_tts(step, nrCollisions, nrIonizations)
      case (2) ! Continuous ionization and discrete recombination
        call Do_Continuous_Ionization_tts(step, nrCollisions, nrIonizations)
        call Do_Discrete_Recombination_tts(step, nrRecombinations)
        nrCollisions = nrCollisions + nrRecombinations
      case (3) ! Discrete ionization
        call Do_Discrete_Ionization_tts(step, nrIonizations)
        nrCollisions = nrIonizations
      case (4) ! Discrete ionization and discrete recombination
        call Do_Discrete_Ionization_tts(step, nrIonizations)
        call Do_Discrete_Recombination_tts(step, nrRecombinations)
        nrCollisions = nrIonizations + nrRecombinations
    end select
  else
    call Update_Collision_Data_All_ots()

    select case (collision_mode)
      case (1) ! Continuous ionization
        call Do_Continuous_Ionization_ots(step, nrCollisions, nrIonizations)
      case (2) ! Continuous ionization and discrete recombination
        call Do_Continuous_Ionization_ots(step, nrCollisions, nrIonizations)
        call Do_Discrete_Recombination_ots(step, nrRecombinations)
        nrCollisions = nrCollisions + nrRecombinations
      case (3) ! Discrete ionization
        call Do_Discrete_Ionization_ots(step, nrIonizations)
        nrCollisions = nrIonizations
      case (4) ! Discrete ionization and discrete recombination
        call Do_Discrete_Ionization_ots(step, nrIonizations)
        call Do_Discrete_Recombination_ots(step, nrRecombinations)
        nrCollisions = nrIonizations + nrRecombinations
    end select
  end if
  end if

  write(ud_coll, '(i6,tr2,i6,tr2,i6,tr2,i6)', iostat=IFAIL) &
          step, nrCollisions, nrIonizations, nrRecombinations
end subroutine Do_Electron_Atom_Collisions

! -----------------------------------------------------------------------------
! --------------------------- RECOMBINATION -----------------------------------
! -----------------------------------------------------------------------------
! We go through all existing ions and electron pairs, fetch the kramers cross section
! radius for the electron's energy and check if the electron will, within the time step,
! enter within one kramers radius of the ion.
! If it does, recombination happens and we remove the electron and ion and add an atom.

subroutine Do_Discrete_Recombination_ots(step,nrRecombinations)
  integer, intent(in)               ::  step
  integer, intent(out)              ::  nrRecombinations
  ! Ion and electron
  double precision, dimension(1:3)  ::  ion_cur_pos, atom_cur_pos, atom_cur_vel
  double precision, dimension(1:3)  ::  elec_cur_pos, elec_cur_vel, elec_cur_acc, elec_next_pos, rel_pos
  double precision                  ::  cur_dist2, elec_cur_speed, dist
  ! Kramers
  double precision, parameter       ::  multiplicator = 1.0d0
  double precision                  ::  recom_rad2, recom_rad
  ! Polynomial solver
  double precision                  ::  a, b, cc, dd, e, t, t1_r, t1_i, t2_r, t2_i, t3_r, t3_i, t4_r, t4_i
  complex(8)                        ::  root1, root2, root3, root4
  integer                           ::  code
  ! Logicals
  logical                           ::  coll_happens
  ! Misc
  integer                           ::  i, j

  nrRecombinations = 0
  
  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP& PRIVATE(i,j,elec_cur_pos, elec_cur_vel, elec_cur_acc, elec_next_pos, rel_pos, cur_dist2) &
  !$OMP& PRIVATE(ion_cur_pos, atom_cur_pos, atom_cur_vel) &
  !$OMP& PRIVATE(elec_cur_speed, dist, recom_rad2, recom_rad) &
  !$OMP& PRIVATE(a, b, cc, dd, e, t, t1_r, t1_i, t2_r, t2_i, t3_r, t3_i, t4_r, t4_i, code, root1, root2, root3, root4) &
  !$OMP& PRIVATE(coll_happens) &
  !$OMP& SHARED(nrPart,nrIon,nrElec,particles_ion_pointer,particles_elec_pointer, particles_species, particles_life) &
  !$OMP& SHARED(particles_cur_pos, particles_mask, particles_step, particles_cur_vel, particles_cur_accel) &
  !$OMP& SHARED(step, time_step, particles_emitter, n_d) &
  !$OMP& SHARED(particles_cur_energy, particles_recom_cross_rad, collision_mode) &
  !$OMP& REDUCTION(+:nrRecombinations)

  ! Cycle through all ions
  ion: do i = 1, nrPart
    if ((particles_species(i) /= species_ion) .or. (particles_mask(i) .eqv. .false.)) cycle ion
    ! Remove ion if its life time is over
    if (step >= particles_life(i)) then ! End of life
      call Mark_Particles_Remove(i, remove_top)
      cycle
    end if 
    
    ion_cur_pos = particles_cur_pos(:,i)

    ! Cycle through all electrons
    elec: do j=1, nrPart
      if ((particles_species(j) /= species_elec) .or. (particles_mask(j) .eqv. .false.)) cycle elec

      ! Fetch electron data
      elec_cur_pos = particles_cur_pos(:,j)
      elec_cur_vel = particles_cur_vel(:,j) 
      elec_cur_acc = particles_cur_accel(:,j)
      elec_cur_speed = norm2(elec_cur_vel)

      rel_pos = elec_cur_pos - ion_cur_pos
      cur_dist2 = dot_product(rel_pos,rel_pos)

      recom_rad = particles_recom_cross_rad(j)*multiplicator
      recom_rad2 = recom_rad**2

      ! print*, 'Recombination radius = ', recom_rad
    
      coll_happens = .false.

      if (cur_dist2 <= recom_rad2) then
        coll_happens = .true.
        t = 0
        dist = norm2(elec_cur_pos - ion_cur_pos)
      else
        a = 0.25d0*(dot_product(elec_cur_acc,elec_cur_acc))
        b = dot_product(elec_cur_vel,elec_cur_acc)
        cc = dot_product(elec_cur_vel,elec_cur_vel) + dot_product(rel_pos,elec_cur_acc)
        dd = 2.0d0*dot_product(rel_pos,elec_cur_vel)
        e = cur_dist2 - recom_rad2

        call SolvePolynomial(a,b,cc,dd,e,code,root1,root2,root3,root4)

        if (code /= 44 .and. code /= 23) then
          t1_r = root1%re
          t1_i = root1%im
          t2_r = root2%re
          t2_i = root2%im
          t3_r = root3%re
          t3_i = root3%im
          t4_r = root4%re
          t4_i = root4%im

          if (code == 31) then
            if (t1_i == 0 .and. t1_r > 0 .and. t1_r <= time_step) then
              coll_happens = .true.
              t = t1_r
            else if (t2_i == 0 .and. t2_r > 0 .and. t2_r <= time_step) then
              coll_happens = .true.
              t = t2_r
            else if (t3_i == 0 .and. t3_r > 0 .and. t3_r <= time_step) then
              coll_happens = .true.
              t = t3_r
            else if (t4_i == 0 .and. t4_r > 0 .and. t4_r <= time_step) then
              coll_happens = .true.
              t = t4_r
            end if            
          else if (code == 42) then 
            if (t1_i == 0 .and. t1_r > 0 .and. t1_r <= time_step) then
              coll_happens = .true.
              t = t1_r
            else if (t2_i == 0 .and. t2_r > 0 .and. t2_r <= time_step) then
              coll_happens = .true.
              t = t2_r
            end if
          end if
        end if
      end if

      ! Recombination
      if (coll_happens .eqv. .true.) then
        elec_next_pos = elec_cur_pos + elec_cur_vel*t+0.5d0*elec_cur_acc*t**2
        dist = norm2(elec_next_pos - ion_cur_pos)
        nrRecombinations = nrRecombinations + 1

        ! Remove ion always
        call Mark_Particles_Remove(i, remove_recom)

        ! Remove electron always
        call Mark_Particles_Remove(j, remove_recom)

        ! Add atom only if ionization is discrete 
        if (collision_mode == 4) then
          atom_cur_pos = ion_cur_pos
          atom_cur_vel = (/0.0d0,0.0d0,0.0d0/)
          !$OMP CRITICAL
          call Add_Particle(atom_cur_pos,atom_cur_vel,species_atom,step,recom_emitter,-1)
          !$OMP END CRITICAL
        end if

        ! Write recombination data
        call Write_Recombination_Data(step, ion_cur_pos, elec_cur_speed, dist, recom_rad, j, i, particles_emitter(j))
        cycle ion       
      end if
    end do elec
  end do ion
  !$OMP END PARALLEL DO
end subroutine Do_Discrete_Recombination_ots

subroutine Do_Discrete_Recombination_tts(step,nrRecombinations)
  integer, intent(in)               ::  step
  integer, intent(out)              ::  nrRecombinations
  ! Ion and electron
  double precision, dimension(1:3)  ::  ion_cur_pos, atom_cur_pos, atom_cur_vel
  double precision, dimension(1:3)  ::  elec_cur_pos, elec_cur_vel, elec_cur_acc, elec_next_pos, rel_pos
  double precision                  ::  cur_dist2, elec_cur_speed, dist
  ! Kramers
  double precision, parameter       ::  multiplicator = 1.0d0
  double precision                  ::  recom_rad2, recom_rad
  ! Polynomial solver
  double precision                  ::  a, b, cc, dd, e, t, t1_r, t1_i, t2_r, t2_i, t3_r, t3_i, t4_r, t4_i
  complex(8)                        ::  root1, root2, root3, root4
  integer                           ::  code
  ! Logicals
  logical                           ::  coll_happens
  ! Misc
  integer                           ::  i, j, k, l

  nrRecombinations = 0
  
  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP& PRIVATE(i,j,k,l,elec_cur_pos, elec_cur_vel, elec_cur_acc, elec_next_pos, rel_pos, cur_dist2) &
  !$OMP& PRIVATE(ion_cur_pos, atom_cur_pos, atom_cur_vel) &
  !$OMP& PRIVATE(elec_cur_speed, dist, recom_rad2, recom_rad) &
  !$OMP& PRIVATE(a, b, cc, dd, e, t, t1_r, t1_i, t2_r, t2_i, t3_r, t3_i, t4_r, t4_i, code, root1, root2, root3, root4) &
  !$OMP& PRIVATE(coll_happens) &
  !$OMP& SHARED(nrPart,nrIon,nrElec,particles_ion_pointer,particles_elec_pointer, particles_species, particles_life) &
  !$OMP& SHARED(particles_cur_pos, particles_mask, particles_step, particles_cur_vel, particles_cur_accel) &
  !$OMP& SHARED(step, time_step, particles_emitter, n_d) &
  !$OMP& SHARED(particles_cur_energy, particles_recom_cross_rad, collision_mode) &
  !$OMP& REDUCTION(+:nrRecombinations)

  ! Cycle through all ions
  ion: do k = 1, nrIon
    i = particles_ion_pointer(k)
    if (particles_mask(i) .eqv. .false.) cycle ion
    ! Remove ion if its life time is over
    if (step >= particles_life(i)) then ! End of life
      call Mark_Particles_Remove(i, remove_top)
      cycle ion
    end if 
    
    ion_cur_pos = particles_cur_pos(:,i)

    ! Cycle through all electrons
    elec: do l=1, nrElec
      j = particles_elec_pointer(l)
      if ((particles_mask(j) .eqv. .false.)) cycle elec

      ! Fetch electron data
      elec_cur_pos = particles_cur_pos(:,j)
      elec_cur_vel = particles_cur_vel(:,j) 
      elec_cur_acc = particles_cur_accel(:,j)
      elec_cur_speed = norm2(elec_cur_vel)

      rel_pos = elec_cur_pos - ion_cur_pos
      cur_dist2 = dot_product(rel_pos,rel_pos)

      recom_rad = particles_recom_cross_rad(j)*multiplicator
      recom_rad2 = recom_rad**2

      ! print*, 'Recombination radius = ', recom_rad
    
      coll_happens = .false.

      if (cur_dist2 <= recom_rad2) then
        coll_happens = .true.
        t = 0
        dist = norm2(elec_cur_pos - ion_cur_pos)
      else
        a = 0.25d0*(dot_product(elec_cur_acc,elec_cur_acc))
        b = dot_product(elec_cur_vel,elec_cur_acc)
        cc = dot_product(elec_cur_vel,elec_cur_vel) + dot_product(rel_pos,elec_cur_acc)
        dd = 2.0d0*dot_product(rel_pos,elec_cur_vel)
        e = cur_dist2 - recom_rad2

        call SolvePolynomial(a,b,cc,dd,e,code,root1,root2,root3,root4)

        if (code /= 44 .and. code /= 23) then
          t1_r = root1%re
          t1_i = root1%im
          t2_r = root2%re
          t2_i = root2%im
          t3_r = root3%re
          t3_i = root3%im
          t4_r = root4%re
          t4_i = root4%im

          if (code == 31) then
            if (t1_i == 0 .and. t1_r > 0 .and. t1_r <= time_step) then
              coll_happens = .true.
              t = t1_r
            else if (t2_i == 0 .and. t2_r > 0 .and. t2_r <= time_step) then
              coll_happens = .true.
              t = t2_r
            else if (t3_i == 0 .and. t3_r > 0 .and. t3_r <= time_step) then
              coll_happens = .true.
              t = t3_r
            else if (t4_i == 0 .and. t4_r > 0 .and. t4_r <= time_step) then
              coll_happens = .true.
              t = t4_r
            end if            
          else if (code == 42) then 
            if (t1_i == 0 .and. t1_r > 0 .and. t1_r <= time_step) then
              coll_happens = .true.
              t = t1_r
            else if (t2_i == 0 .and. t2_r > 0 .and. t2_r <= time_step) then
              coll_happens = .true.
              t = t2_r
            end if
          end if
        end if
      end if

      ! coll_happens = .false. ! REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Recombination
      if (coll_happens .eqv. .true.) then
        elec_next_pos = elec_cur_pos + elec_cur_vel*t+0.5d0*elec_cur_acc*t**2
        dist = norm2(elec_next_pos - ion_cur_pos)
        nrRecombinations = nrRecombinations + 1

        ! Remove ion always
        call Mark_Particles_Remove(i, remove_recom)

        ! Remove electron always
        call Mark_Particles_Remove(j, remove_recom)

        ! Add atom only if ionization is discrete 
        if (collision_mode == 4) then
          atom_cur_pos = ion_cur_pos
          atom_cur_vel = (/0.0d0,0.0d0,0.0d0/)
          !$OMP CRITICAL
          call Add_Particle(atom_cur_pos,atom_cur_vel,species_atom,step,recom_emitter,-1)
          !$OMP END CRITICAL
        end if

        ! Write recombination data
        call Write_Recombination_Data(step, ion_cur_pos, elec_cur_speed, dist, recom_rad, k, i, particles_emitter(j))
        cycle ion       
      end if
    end do elec
  end do ion
  !$OMP END PARALLEL DO
end subroutine Do_Discrete_Recombination_tts

! -----------------------------------------------------------------------------
! --------------------------- IONIZATION -------------------------------------
! -----------------------------------------------------------------------------

subroutine Do_Continuous_Ionization_ots_old(step, nrIonizations, nrCollisions)
  integer, intent(in)               ::  step
    integer, intent(out)              ::  nrIonizations, nrCollisions
    ! Parameters
    double precision, parameter       ::  elec_max_speed2 = (2.0d0*q_0*5000.0d0/m_0)
    ! Electrons and ion
    double precision, dimension(1:3)  ::  elec_cur_pos, elec_prev_pos, elec_cur_vel, direct_vec
    double precision, dimension(1:3)  ::  ejec_elec_pos, ejec_elec_vel
    double precision, dimension(1:3)  ::  ion_pos, ion_vel
    double precision                  ::  elec_energy, elec_cur_path, elec_cur_speed
    integer                           ::  newID, ionID
    ! Ionization
    double precision                  ::  E1, E2, collE, ejecE, inSpeed, outSpeed, newSpeed
    integer                           ::  count_n
    double precision                  ::  cross_tot, cross_ion, mean_path, mean_path_avg
    ! Misc
    integer                           ::  i
    double precision                  ::  rnd, alpha


    nrIonizations = 0
    nrCollisions = 0
    count_n = 0

    if (nrElec == 0) return
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP& PRIVATE(elec_cur_pos,elec_prev_pos,elec_cur_vel,direct_vec,ejec_elec_pos,ejec_elec_vel,ion_pos,ion_vel) &
    !$OMP& PRIVATE(elec_energy,elec_cur_path,elec_cur_speed,newID,ionID,E1,E2,collE,ejecE,inSpeed,outSpeed,newSpeed) &
    !$OMP& PRIVATE(cross_tot,cross_ion,mean_path,mean_path_avg,i,rnd,alpha) &
    !$OMP& SHARED(particles_species,particles_mask,particles_cur_pos,particles_prev_pos,particles_cur_vel) &
    !$OMP& SHARED(particles_cur_energy,particles_ion_cross_sec,particles_tot_cross_sec,nrID,step) &
    !$OMP& SHARED(n_d,nrPart,ion_life_time,particles_emitter,emitters_dim) &
    !$OMP& REDUCTION(+:nrIonizations,nrCollisions,count_n)

    do i=1,nrPart
      if ((particles_species(i) == species_elec) .and. (particles_mask(i) .eqv. .true.)) then
        elec_cur_pos(:) = particles_cur_pos(:, i)
        elec_prev_pos(:) = particles_prev_pos(:, i)
        elec_cur_vel = particles_cur_vel(:, i)

        if (norm2(elec_cur_pos(1:2)) <= emitters_dim(1,1)) then

          elec_energy = particles_cur_energy(i)
          elec_cur_path = norm2(elec_cur_pos - elec_prev_pos)
          elec_cur_speed = norm2(elec_cur_vel)

          if (elec_energy > N_bind) then

            cross_tot = particles_tot_cross_sec(i)

            mean_path = 1.0d0/(n_d*cross_tot)
            mean_path_avg = mean_path_avg + mean_path
            count_n = count_n + 1

            ! Calculate the ratio between the distance travled and the mean free path
            alpha = elec_cur_path/mean_path

            ! Check if we do a collision or not
            call random_number(rnd)
            if (rnd < alpha) then

              !---------------------------------------------
              ! Check if the collision ionizes the N2 or not  
              cross_ion = particles_ion_cross_sec(i)

              alpha = cross_ion/cross_tot

              call random_number(rnd)
              if (rnd < alpha) then ! Check if we ionize or not

                ! --------------------- Energy conservation ---------------------
                E1 = elec_energy ! Starting energy
                E2 = E1 - N_bind ! Energy after ionization
                call random_number(rnd)
                collE = E2*rnd ! Energy of the colliding electron
                ejecE = E2 - collE ! Energy of the ejected electron

                ! --------------------- Colliding electron ----------------------
                ! New direction
                direct_vec = Get_Injected_Vec(elec_energy, elec_cur_vel)
                direct_vec = direct_vec / norm2(direct_vec)
                ! New velocity vector
                particles_cur_vel(:, i) = direct_vec*sqrt(2.0d0*q_0*collE/m_0)
                !$OMP CRITICAL
                call Update_Collision_Data(i)
                !$OMP END CRITICAL
                ! Some variables for the file
                inSpeed = elec_cur_speed
                outSpeed = norm2(particles_cur_vel(:,i))

                ! --------------------- Ejected electron ---------------------   
                ! New position
                call random_number(ejec_elec_pos)
                ejec_elec_pos = ejec_elec_pos - 0.5d0
                ejec_elec_pos = elec_cur_pos + ejec_elec_pos*length_scale
                ! New direction
                direct_vec = Get_Ejected_Vec(elec_energy, elec_energy, elec_cur_vel)
                direct_vec = direct_vec / norm2(direct_vec)
                ! New velocity vector
                ejec_elec_vel = direct_vec*sqrt(2.0d0*q_0*ejecE/m_0)
                ! Some data for the file
                newSpeed = norm2(ejec_elec_vel)

                ! ----------------------- Created ion -----------------------
                ion_pos = elec_cur_pos
                ion_vel = 0.0d0

                ! --------------------- Ionization ---------------------------
                ! Add the new electron to the system
                !$OMP CRITICAL
                newID = nrID
                call Add_Particle(ejec_elec_pos, ejec_elec_vel, species_elec, step, ion_emitter, -1) ! Electron
                call Update_Collision_Data(newID)
                !$OMP END CRITICAL

                ! Add the new positively charged ion to the system
                !$OMP CRITICAL
                ionID = nrID
                call Add_Particle(ion_pos, ion_vel, species_ion, step, ion_emitter, step+ion_life_time) ! Ion
                call Write_Ionization_Data(step,elec_cur_pos,inSpeed,outSpeed,newSpeed,0.0d0, &
                & 0.0d0,i,newID,ionID,particles_emitter(i))
                !$OMP END CRITICAL

                nrIonizations = nrIonizations + 1
              end if

              ! Update the number of collisions
              nrCollisions = nrCollisions + 1
            end if
          end if
        end if
      end if
    end do
  !$OMP END PARALLEL DO
end subroutine Do_Continuous_Ionization_ots_old

subroutine Do_Continuous_Ionization_ots(step,nrCollisions,nrIonizations)
  integer, intent(in)               ::  step
  integer, intent(out)              ::  nrIonizations, nrCollisions
  ! Parameters
  double precision, parameter       ::  elec_max_speed2 = (2.0d0*q_0*5000.0d0/m_0)
  ! Electrons and ion
  double precision, dimension(1:3)  ::  elec_cur_pos, elec_prev_pos, elec_cur_vel, direct_vec
  double precision, dimension(1:3)  ::  ejec_elec_pos, ejec_elec_vel
  double precision, dimension(1:3)  ::  ion_pos, ion_vel
  double precision                  ::  elec_energy, elec_cur_path, elec_cur_speed
  integer                           ::  newID, ionID
  ! Ionization
  double precision                  ::  E1, E2, collE, ejecE, inSpeed, outSpeed, newSpeed
  integer                           ::  count_n
  double precision                  ::  cross_tot, cross_ion, mean_path, mean_path_avg
  ! Misc
  integer                           ::  i
  double precision                  ::  rnd, alpha

  nrCollisions = 0
  nrIonizations = 0
  count_n = 0

  if (nrElec == 0) return
  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP& PRIVATE(i,elec_cur_pos,elec_prev_pos,elec_cur_vel,direct_vec,ejec_elec_pos,ejec_elec_vel,ion_pos,ion_vel) &
  !$OMP& PRIVATE(elec_energy,elec_cur_path,elec_cur_speed,newID,ionID,E1,E2,collE,ejecE,inSpeed,outSpeed,newSpeed) &
  !$OMP& PRIVATE(cross_tot,cross_ion,mean_path,mean_path_avg,rnd,alpha) &
  !$OMP& SHARED(nrPart,particles_elec_pointer,particles_species,particles_mask) &	
  !$OMP& SHARED(particles_cur_pos,particles_prev_pos,particles_cur_vel) &
  !$OMP& SHARED(particles_cur_energy,particles_ion_cross_sec,particles_tot_cross_sec,nrID,step) &
  !$OMP& SHARED(n_d,ion_life_time,particles_emitter,emitters_dim) &
  !$OMP& REDUCTION(+:nrIonizations,nrCollisions,count_n)

  elec: do i=1,nrPart

    if ((particles_species(i) /= species_elec) .or. (particles_mask(i) .eqv. .false.)) cycle elec

    elec_cur_pos(:) = particles_cur_pos(:, i)
    elec_prev_pos(:) = particles_prev_pos(:, i)
    elec_cur_vel = particles_cur_vel(:, i)

    if (norm2(elec_cur_pos(1:2)) <= emitters_dim(1,1)) then

      elec_energy = particles_cur_energy(i)
      elec_cur_path = norm2(elec_cur_pos - elec_prev_pos)
      elec_cur_speed = norm2(elec_cur_vel)

      if (elec_energy > N_bind) then

        cross_tot = particles_tot_cross_sec(i)

        mean_path = 1.0d0/(n_d*cross_tot)
        mean_path_avg = mean_path_avg + mean_path
        count_n = count_n + 1

        ! Calculate the ratio between the distance travled and the mean free path
        alpha = elec_cur_path/mean_path

        ! Check if we do a collision or not
        call random_number(rnd)
        if (rnd < alpha) then

          !---------------------------------------------
          ! Check if the collision ionizes the N2 or not  
          cross_ion = particles_ion_cross_sec(i)
          
          alpha = cross_ion/cross_tot
          
          call random_number(rnd)
          if (rnd < alpha) then ! Check if we ionize or not

            ! --------------------- Energy conservation ---------------------
            E1 = elec_energy ! Starting energy
            E2 = E1 - N_bind ! Energy after ionization
            call random_number(rnd)
            collE = E2*rnd ! Energy of the colliding electron
            ejecE = E2 - collE ! Energy of the ejected electron

            ! --------------------- Colliding electron ----------------------
            ! New direction
            direct_vec = Get_Injected_Vec(elec_energy, elec_cur_vel)
            direct_vec = direct_vec / norm2(direct_vec)
            ! New velocity vector
            particles_cur_vel(:, i) = direct_vec*sqrt(2.0d0*q_0*collE/m_0)
            !$OMP CRITICAL
            call Update_Collision_Data(i)
            !$OMP END CRITICAL
            ! Some variables for the file
            inSpeed = elec_cur_speed
            outSpeed = norm2(particles_cur_vel(:,i))
            
            ! --------------------- Ejected electron ---------------------   
            ! New position
            call random_number(ejec_elec_pos)
            ejec_elec_pos = 2.0d0*(ejec_elec_pos - 0.5d0) ! allow negative_values
            ejec_elec_pos = elec_cur_pos + ejec_elec_pos*length_scale
            ! New direction
            direct_vec = Get_Ejected_Vec(elec_energy, elec_energy, elec_cur_vel)
            direct_vec = direct_vec / norm2(direct_vec)
            ! New velocity vector
            ejec_elec_vel = direct_vec*sqrt(2.0d0*q_0*ejecE/m_0)
            ! Some data for the file
            newSpeed = norm2(ejec_elec_vel)

            ! ----------------------- Created ion -----------------------
            call random_number(ion_pos)
            ion_pos = 2.0d0*(ion_pos - 0.5d0) ! allow negative values
            ion_pos = elec_cur_pos + ion_pos*length_scale
            ion_vel = 0.0d0

            ! --------------------- Ionization ---------------------------
            ! Add the new electron to the system
            !$OMP CRITICAL
            newID = nrID
            call Add_Particle(ejec_elec_pos, ejec_elec_vel, species_elec, step, ion_emitter, -1) ! Electron
            !$OMP END CRITICAL
            !$OMP CRITICAL
            call Update_Collision_Data(newID)
            !$OMP END CRITICAL

            ! Add the new positively charged ion to the system
            !$OMP CRITICAL
            ionID = nrID
            call Add_Particle(ion_pos, ion_vel, species_ion, step, ion_emitter, step+ion_life_time) ! Ion
            !$OMP END CRITICAL
            !$OMP CRITICAL
            call Write_Ionization_Data(step,elec_cur_pos,inSpeed,outSpeed,newSpeed,0.0d0, &
            & 0.0d0,i,newID,ionID,particles_emitter(i))
            !$OMP END CRITICAL

            nrIonizations = nrIonizations + 1

            particles_emitter(i) = ion_emitter
          end if

          ! Update the number of collisions
          nrCollisions = nrCollisions + 1
        end if
      end if
    end if
  end do elec
  !$OMP END PARALLEL DO
end subroutine Do_Continuous_Ionization_ots

subroutine Do_Continuous_Ionization_tts(step,nrCollisions,nrIonizations)
  integer, intent(in)               ::  step
  integer, intent(out)              ::  nrIonizations, nrCollisions
  ! Parameters
  double precision, parameter       ::  elec_max_speed2 = (2.0d0*q_0*5000.0d0/m_0)
  ! Electrons and ion
  double precision, dimension(1:3)  ::  elec_cur_pos, elec_prev_pos, elec_cur_vel, direct_vec
  double precision, dimension(1:3)  ::  ejec_elec_pos, ejec_elec_vel
  double precision, dimension(1:3)  ::  ion_pos, ion_vel
  double precision                  ::  elec_energy, elec_cur_path, elec_cur_speed
  integer                           ::  newID, ionID
  ! Ionization
  double precision                  ::  E1, E2, collE, ejecE, inSpeed, outSpeed, newSpeed
  integer                           ::  count_n
  double precision                  ::  cross_tot, cross_ion, mean_path, mean_path_avg
  ! Misc
  integer                           ::  i, k
  double precision                  ::  rnd, alpha

  nrCollisions = 0
  nrIonizations = 0
  count_n = 0

  if (nrElec == 0) return
  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP& PRIVATE(k,i,elec_cur_pos,elec_prev_pos,elec_cur_vel,direct_vec,ejec_elec_pos,ejec_elec_vel,ion_pos,ion_vel) &
  !$OMP& PRIVATE(elec_energy,elec_cur_path,elec_cur_speed,newID,ionID,E1,E2,collE,ejecE,inSpeed,outSpeed,newSpeed) &
  !$OMP& PRIVATE(cross_tot,cross_ion,mean_path,mean_path_avg,rnd,alpha) &
  !$OMP& SHARED(nrPart,nrElec,particles_elec_pointer,particles_species,particles_mask) &	
  !$OMP& SHARED(particles_cur_pos,particles_prev_pos,particles_cur_vel) &
  !$OMP& SHARED(particles_cur_energy,particles_ion_cross_sec,particles_tot_cross_sec,nrID,step) &
  !$OMP& SHARED(n_d,ion_life_time,particles_emitter,emitters_dim) &
  !$OMP& REDUCTION(+:nrIonizations,nrCollisions,count_n)

  elec: do k=1,nrElec
    i = particles_elec_pointer(k)

    if (particles_mask(i) .eqv. .false.) cycle elec
    elec_cur_pos(:) = particles_cur_pos(:, i)
    elec_prev_pos(:) = particles_prev_pos(:, i)
    elec_cur_vel = particles_cur_vel(:, i)

    if (norm2(elec_cur_pos(1:2)) <= emitters_dim(1,1)) then

      elec_energy = particles_cur_energy(i)
      elec_cur_path = norm2(elec_cur_pos - elec_prev_pos)
      elec_cur_speed = norm2(elec_cur_vel)

      if (elec_energy > N_bind) then

        ! print *, 'passed'

        cross_tot = particles_tot_cross_sec(i)

        mean_path = 1.0d0/(n_d*cross_tot)
        mean_path_avg = mean_path_avg + mean_path
        count_n = count_n + 1

        ! Calculate the ratio between the distance travled and the mean free path
        alpha = elec_cur_path/mean_path

        ! Check if we do a collision or not
        call random_number(rnd)
        if (rnd < alpha) then

          !---------------------------------------------
          ! Check if the collision ionizes the N2 or not  
          cross_ion = particles_ion_cross_sec(i)
          
          alpha = cross_ion/cross_tot
          
          call random_number(rnd)
          if (rnd < alpha) then ! Check if we ionize or not

            ! print *, 'mean_path = ', mean_path, 'elec_cur_path = ', elec_cur_path

            ! --------------------- Energy conservation ---------------------
            E1 = elec_energy ! Starting energy
            E2 = E1 - N_bind ! Energy after ionization
            call random_number(rnd)
            collE = E2*rnd ! Energy of the colliding electron
            ejecE = E2 - collE ! Energy of the ejected electron

            ! --------------------- Colliding electron ----------------------
            ! New direction
            direct_vec = Get_Injected_Vec(elec_energy, elec_cur_vel)
            direct_vec = direct_vec / norm2(direct_vec)
            ! New velocity vector
            particles_cur_vel(:, i) = direct_vec*sqrt(2.0d0*q_0*collE/m_0)
            call Update_Collision_Data(i)
            ! Some variables for the file
            inSpeed = elec_cur_speed
            outSpeed = norm2(particles_cur_vel(:,i))
            
            ! --------------------- Ejected electron ---------------------   
            ! New position
            call random_number(ejec_elec_pos)
            ejec_elec_pos = 2.0d0*(ejec_elec_pos - 0.5d0) ! allow negative_values
            ejec_elec_pos = elec_cur_pos + ejec_elec_pos*length_scale
            ! New direction
            direct_vec = Get_Ejected_Vec(elec_energy, elec_energy, elec_cur_vel)
            direct_vec = direct_vec / norm2(direct_vec)
            ! New velocity vector
            ejec_elec_vel = direct_vec*sqrt(2.0d0*q_0*ejecE/m_0)
            ! Some data for the file
            newSpeed = norm2(ejec_elec_vel)

            ! ----------------------- Created ion -----------------------
            call random_number(ion_pos)
            ion_pos = 2.0d0*(ion_pos - 0.5d0) ! allow negative values
            ion_pos = elec_cur_pos + ion_pos*length_scale
            ion_vel = 0.0d0

            ! --------------------- Ionization ---------------------------
            ! Add the new electron to the system
            !$OMP CRITICAL
            newID = nrID
            call Add_Particle(ejec_elec_pos, ejec_elec_vel, species_elec, step, ion_emitter, -1) ! Electron
            !$OMP END CRITICAL
            !$OMP CRITICAL
            call Update_Collision_Data(newID)
            !$OMP END CRITICAL

            ! Add the new positively charged ion to the system
            !$OMP CRITICAL
            ionID = nrID
            call Add_Particle(ion_pos, ion_vel, species_ion, step, ion_emitter, step+ion_life_time) ! Ion
            !$OMP END CRITICAL
            !$OMP CRITICAL
            call Write_Ionization_Data(step,elec_cur_pos,inSpeed,outSpeed,newSpeed,0.0d0, &
            & 0.0d0,i,newID,ionID,particles_emitter(i))
            !$OMP END CRITICAL

            nrIonizations = nrIonizations + 1

            particles_emitter(i) = ion_emitter
          end if

          ! Update the number of collisions
          nrCollisions = nrCollisions + 1
        end if
      end if
    end if
  end do elec
  !$OMP END PARALLEL DO
end subroutine Do_Continuous_Ionization_tts

subroutine Do_Discrete_Ionization_ots(step,nrIonizations)
  integer, intent(in)               ::  step
  integer, intent(out)              ::  nrIonizations
  ! Parameters
  double precision, parameter       ::  elec_max_speed2 = (2.0d0*q_0*5000.0d0/m_0)
  ! Atom, ion, and electron
  double precision, dimension(1:3)  ::  atom_cur_pos,elec_cur_pos,elec_next_pos,elec_cur_vel,elec_cur_acc
  double precision, dimension(1:3)  ::  ejec_elec_pos, ejec_elec_vel, ion_pos, ion_vel, relative_pos,direct_vec
  double precision                  ::  elec_cur_speed, elec_energy, cur_dist2, ionization_dist
  ! Ionization
  logical                           ::  ionization
  double precision                  ::  ion_rad2, ion_rad, E1, E2, collE, ejecE, inSpeed, outSpeed, newSpeed
  integer                           ::  newID, ionID
  ! Polynomial solver
  double precision                  ::  a,b,cc,dd,e,t,t1_r,t1_i,t2_r,t2_i,t3_r,t3_i,t4_r,t4_i
  complex(8)                        ::  root1,root2,root3,root4
  integer                           ::  code
  ! Misc
  integer                           ::  i,j
  double precision                  ::  rnd

  nrIonizations = 0

  if (nrElec == 0) return
  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP& PRIVATE(i,j,atom_cur_pos,elec_cur_pos,elec_next_pos,elec_cur_vel,elec_cur_acc) &
  !$OMP& PRIVATE(ejec_elec_pos,ejec_elec_vel,ion_pos,ion_vel,relative_pos,direct_vec) &
  !$OMP& PRIVATE(elec_cur_speed,elec_energy,cur_dist2,ionization_dist,ionization,ion_rad2,ion_rad) &
  !$OMP& PRIVATE(E1,E2,collE,ejecE,inSpeed,outSpeed,newSpeed,newID,ionID) &
  !$OMP& PRIVATE (a,b,cc,dd,e,t,t1_r,t1_i,t2_r,t2_i,t3_r,t3_i,t4_r,t4_i,root1,root2,root3,root4,code,rnd) &
  !$OMP& SHARED(nrPart,step,time_step,nrID,ion_life_time,particles_emitter,particles_step) &
  !$OMP& SHARED(particles_species,particles_mask,particles_cur_pos,particles_cur_vel,particles_cur_accel) &
  !$OMP& SHARED(particles_cur_energy, particles_ion_cross_rad) &
  !$OMP& REDUCTION(+:nrIonizations)

  ! Cycle through all atoms
  atom: do i = 1,nrPart
    if ((particles_species(i) /= species_atom) .or. (particles_mask(i) .eqv. .false.)) cycle atom

    atom_cur_pos = particles_cur_pos(:,i)

    ! Cycle through all electrons
    elec: do j = 1,nrPart

      if ((particles_species(j) /= species_elec) .or. (particles_mask(j) .eqv. .false.)) cycle elec
        
      ! Fetch electron data
      elec_cur_pos = particles_cur_pos(:,j)
      elec_cur_vel = particles_cur_vel(:,j)
      elec_cur_acc = particles_cur_accel(:,j)
      elec_cur_speed = norm2(elec_cur_vel)
      elec_energy = particles_cur_energy(j)
      
      relative_pos = elec_cur_pos - atom_cur_pos
      cur_dist2 = dot_product(relative_pos,relative_pos)

      ! Check ionization only if electron has enough energy to ionize
      if (elec_energy >= N_bind) then
        ! print *, 'Ionization energy reached'
        ion_rad = particles_ion_cross_rad(j)
        ion_rad2 = ion_rad**2
        ! print *, 'Ionization radius = ', ion_rad 

        ionization = .false.

        ! Check if ionization happens
        if (cur_dist2 <= ion_rad2) then
          ionization = .true.
          t = 0
          ionization_dist = norm2(elec_cur_pos - atom_cur_pos)
        else
          a = 0.25d0*(dot_product(elec_cur_acc,elec_cur_acc))
          b = dot_product(elec_cur_vel,elec_cur_acc)
          cc = dot_product(elec_cur_vel,elec_cur_vel) + dot_product(relative_pos,elec_cur_acc)
          dd = 2.0d0*dot_product(relative_pos,elec_cur_vel)
          e = cur_dist2 - ion_rad2

          call SolvePolynomial(a,b,cc,dd,e,code,root1,root2,root3,root4)

          if (code /= 44 .and. code /= 23) then
            t1_r = root1%re
            t1_i = root1%im
            t2_r = root2%re
            t2_i = root2%im
            t3_r = root3%re
            t3_i = root3%im
            t4_r = root4%re
            t4_i = root4%im

            if ((abs(t1_i) <= 2.0d-16) .and. (t1_r > 0) .and. (t1_r <= time_step)) then
              ionization = .true.
              t = t1_r
            else if ((abs(t2_i) <= 2.0d-16) .and. (t2_r > 0) .and. (t2_r <= time_step)) then
              ionization = .true.
              t = t2_r
            else if ((abs(t3_i) <= 2.0d-16) .and. (t3_r > 0) .and. (t3_r <= time_step)) then
              ionization = .true.
              t = t3_r
            else if ((abs(t4_i) <= 2.0d-16) .and. (t4_r > 0) .and. (t4_r <= time_step)) then
              ionization = .true.
              t = t4_r
            end if

          end if
        end if

        ! Do ionization
        if (ionization .eqv. .true.) then
          ! print *, 'Atom species = ', particles_species(i), '=', species_atom
          ! print *, 'Ionization distance reached'
          ! Some data for the file
          elec_next_pos = elec_cur_pos + elec_cur_vel*t+0.5d0*elec_cur_acc*t**2
          ionization_dist = norm2(elec_next_pos - atom_cur_pos)
          
          ! --------------------- Energy conservation ---------------------
          E1 = elec_energy ! Starting energy
          E2 = E1 - N_bind ! Energy after ionization
          call random_number(rnd)
          collE = E2*rnd ! Energy of the colliding electron
          ejecE = E2 - collE ! Energy of the ejected electron

          ! --------------------- Colliding electron ----------------------
          ! New direction
          direct_vec = Get_Injected_Vec(elec_energy, elec_cur_vel)
          direct_vec = direct_vec / norm2(direct_vec)
          ! New velocity vector
          particles_cur_vel(:, j) = direct_vec*sqrt(2.0d0*q_0*collE/m_0)
          !$OMP CRITICAL
          call Update_Collision_Data(j)
          !$OMP END CRITICAL
          ! Some variables for the file
          inSpeed = elec_cur_speed
          outSpeed = norm2(particles_cur_vel(:,j))

          ! --------------------- Ejected electron ---------------------
          
          ! New position
          call random_number(ejec_elec_pos)
          ejec_elec_pos = ejec_elec_pos - 0.5d0
          ejec_elec_pos = atom_cur_pos + ejec_elec_pos*length_scale
          ! New direction
          direct_vec = Get_Ejected_Vec(elec_energy, elec_energy, elec_cur_vel)
          direct_vec = direct_vec / norm2(direct_vec)
          ! New velocity vector
          ejec_elec_vel = direct_vec*sqrt(2.0d0*q_0*ejecE/m_0)
          ! Some data for the file
          newSpeed = norm2(ejec_elec_vel)

          ! ----------------------- Created ion -----------------------
          ion_pos = atom_cur_pos
          ion_vel = 0.0d0

          ! --------------------- Ionization ---------------------------
          !$OMP CRITICAL
          ! Add the new electron to the system
          newID = nrID
          call Add_Particle(ejec_elec_pos, ejec_elec_vel, species_elec, step, ion_emitter, -1) ! Electron
          !$OMP END CRITICAL
          !$OMP CRITICAL
          call Update_Collision_Data(newID)
          !$OMP END CRITICAL

          !$OMP CRITICAL
          ! Remove the atom from the system
          call Mark_Particles_Remove(i,remove_ion)
          !$OMP END CRITICAL

          !$OMP CRITICAL
          ! Add the new positively charged ion to the system
          ionID = nrID
          call Add_Particle(ion_pos, ion_vel, species_ion, step, ion_emitter, step+ion_life_time) ! Ion
          !$OMP END CRITICAL
          !$OMP CRITICAL
          call Write_Ionization_Data(step,atom_cur_pos,inSpeed,outSpeed,newSpeed,ionization_dist, &
          & ion_rad,i,newID,ionID,particles_emitter(i))
          !$OMP END CRITICAL

          nrIonizations = nrIonizations + 1

          particles_emitter(j) = ion_emitter

          cycle atom
        end if
      end if
    end do elec
  ! end if
  end do atom
  !$OMP END PARALLEL DO
end subroutine Do_Discrete_Ionization_ots

subroutine Do_Discrete_Ionization_tts(step,nrIonizations)
  integer, intent(in)               ::  step
  integer, intent(out)              ::  nrIonizations
  ! Parameters
  double precision, parameter       ::  elec_max_speed2 = (2.0d0*q_0*5000.0d0/m_0)
  ! Atom, ion, and electron
  double precision, dimension(1:3)  ::  atom_cur_pos,elec_cur_pos,elec_next_pos,elec_cur_vel,elec_cur_acc
  double precision, dimension(1:3)  ::  ejec_elec_pos, ejec_elec_vel, ion_pos, ion_vel, relative_pos,direct_vec
  double precision                  ::  elec_cur_speed, elec_energy, cur_dist2, ionization_dist
  ! Ionization
  logical                           ::  ionization
  double precision                  ::  ion_rad2, ion_rad, E1, E2, collE, ejecE, inSpeed, outSpeed, newSpeed
  integer                           ::  newID, ionID
  ! Polynomial solver
  double precision                  ::  a,b,cc,dd,e,t,t1_r,t1_i,t2_r,t2_i,t3_r,t3_i,t4_r,t4_i
  complex(8)                        ::  root1,root2,root3,root4
  integer                           ::  code
  ! Misc
  integer                           ::  i,j,k,l
  double precision                  ::  rnd

  nrIonizations = 0

  if (nrElec == 0) return
  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP& PRIVATE(i,j,k,l,atom_cur_pos,elec_cur_pos,elec_next_pos,elec_cur_vel,elec_cur_acc) &
  !$OMP& PRIVATE(ejec_elec_pos,ejec_elec_vel,ion_pos,ion_vel,relative_pos,direct_vec) &
  !$OMP& PRIVATE(elec_cur_speed,elec_energy,cur_dist2,ionization_dist,ionization,ion_rad2,ion_rad) &
  !$OMP& PRIVATE(E1,E2,collE,ejecE,inSpeed,outSpeed,newSpeed,newID,ionID) &
  !$OMP& PRIVATE (a,b,cc,dd,e,t,t1_r,t1_i,t2_r,t2_i,t3_r,t3_i,t4_r,t4_i,root1,root2,root3,root4,code,rnd) &
  !$OMP& SHARED(nrAtom,nrElec,particles_atom_pointer,particles_elec_pointer) &
  !$OMP& SHARED(step,time_step,nrID,ion_life_time,particles_emitter,particles_step) &
  !$OMP& SHARED(particles_species,particles_mask,particles_cur_pos,particles_cur_vel,particles_cur_accel) &
  !$OMP& SHARED(particles_cur_energy, particles_ion_cross_rad) &
  !$OMP& REDUCTION(+:nrIonizations)

  ! Cycle through all atoms
  atom: do k = 1,nrAtom
    i = particles_atom_pointer(k)

    ! Move to next if atom has already been ionized
    if (particles_mask(i) .eqv. .false.) cycle atom

    atom_cur_pos = particles_cur_pos(:,i)

    ! Cycle through all electrons
    elec: do l = 1,nrElec
      j = particles_elec_pointer(l)

      if (particles_mask(j) .eqv. .false.) cycle elec
        
      ! Fetch electron data
      elec_cur_pos = particles_cur_pos(:,j)
      elec_cur_vel = particles_cur_vel(:,j)
      elec_cur_acc = particles_cur_accel(:,j)
      elec_cur_speed = norm2(elec_cur_vel)
      elec_energy = particles_cur_energy(j)
      
      relative_pos = elec_cur_pos - atom_cur_pos
      cur_dist2 = dot_product(relative_pos,relative_pos)

      ! Check ionization only if electron has enough energy to ionize
      if (elec_energy >= N_bind) then
        ! print *, 'Ionization energy reached'
        ion_rad2 = particles_ion_cross_rad(j)
        ! print *, 'Ionization radius = ', ion_rad2 

        ionization = .false.

        ! Check if ionization happens
        call random_number(rnd)
        if (cur_dist2 <= ion_rad2) then
        ! if  (rnd <= ion_rad2/cur_dist2) then
          ionization = .true.
          t = 0
          ionization_dist = norm2(elec_cur_pos - atom_cur_pos)
        else
          a = 0.25d0*(dot_product(elec_cur_acc,elec_cur_acc))
          b = dot_product(elec_cur_vel,elec_cur_acc)
          cc = (dot_product(elec_cur_vel,elec_cur_vel) + dot_product(relative_pos,elec_cur_acc))
          dd = 2.0d0*dot_product(relative_pos,elec_cur_vel)
          e = (cur_dist2 - ion_rad2)

          call SolvePolynomial(a,b,cc,dd,e,code,root1,root2,root3,root4)

          if (code /= 44 .and. code /= 23) then
            t1_r = root1%re
            t1_i = root1%im
            t2_r = root2%re
            t2_i = root2%im
            t3_r = root3%re
            t3_i = root3%im
            t4_r = root4%re
            t4_i = root4%im
            
            if ((abs(t1_i) <= 2.22d-16) .and. (t1_r > 0) .and. (t1_r <= time_step)) then
              ionization = .true.
              t = t1_r
            else if ((abs(t2_i) <= 2.22d-16) .and. (t2_r > 0) .and. (t2_r <= time_step)) then
              ionization = .true.
              t = t2_r
            else if ((abs(t3_i) <= 2.22d-16) .and. (t3_r > 0) .and. (t3_r <= time_step)) then
              ionization = .true.
              t = t3_r
            else if ((abs(t4_i) <= 2.22d-16) .and. (t4_r > 0) .and. (t4_r <= time_step)) then
              ionization = .true.
              t = t4_r
            end if

          end if
        end if

      
        ! Do ionization
        if (ionization .eqv. .true.) then
          ! print *, 'Ionization distance reached'
          ! Some data for the file
          elec_next_pos = elec_cur_pos + elec_cur_vel*t+0.5d0*elec_cur_acc*t**2
          ionization_dist = norm2(elec_next_pos - atom_cur_pos)
          
          ! --------------------- Energy conservation ---------------------
          E1 = elec_energy ! Starting energy
          E2 = E1 - N_bind ! Energy after ionization
          call random_number(rnd)
          collE = E2*rnd ! Energy of the colliding electron
          ejecE = E2 - collE ! Energy of the ejected electron

          ! --------------------- Colliding electron ----------------------
          ! New direction
          direct_vec = Get_Injected_Vec(elec_energy, elec_cur_vel)
          direct_vec = direct_vec / norm2(direct_vec)
          ! New velocity vector
          particles_cur_vel(:, j) = direct_vec*sqrt(2.0d0*q_0*collE/m_0)
          !$OMP CRITICAL
          call Update_Collision_Data(j)
          !$OMP END CRITICAL
          ! Some variables for the file
          inSpeed = elec_cur_speed
          outSpeed = norm2(particles_cur_vel(:,j))

          ! --------------------- Ejected electron ---------------------
          
          ! New position
          call random_number(ejec_elec_pos)
          ejec_elec_pos = ejec_elec_pos - 0.5d0
          ejec_elec_pos = atom_cur_pos + ejec_elec_pos*length_scale
          ! New direction
          direct_vec = Get_Ejected_Vec(elec_energy, elec_energy, elec_cur_vel)
          direct_vec = direct_vec / norm2(direct_vec)
          ! New velocity vector
          ejec_elec_vel = direct_vec*sqrt(2.0d0*q_0*ejecE/m_0)
          ! Some data for the file
          newSpeed = norm2(ejec_elec_vel)

          ! ----------------------- Created ion -----------------------
          ion_pos = atom_cur_pos
          ion_vel = 0.0d0

          ! --------------------- Ionization ---------------------------
          !$OMP CRITICAL
          ! Add the new electron to the system
          newID = nrID
          call Add_Particle(ejec_elec_pos, ejec_elec_vel, species_elec, step, ion_emitter, -1) ! Electron
          !$OMP END CRITICAL
          !$OMP CRITICAL
          call Update_Collision_Data(newID)
          !$OMP END CRITICAL

          ! Remove the atom from the system
          call Mark_Particles_Remove(i,remove_ion)

          !$OMP CRITICAL
          ! Add the new positively charged ion to the system
          ionID = nrID
          call Add_Particle(ion_pos, ion_vel, species_ion, step, ion_emitter, step+ion_life_time) ! Ion
          !$OMP END CRITICAL
          !$OMP CRITICAL
          call Write_Ionization_Data(step,atom_cur_pos,inSpeed,outSpeed,newSpeed,ionization_dist, &
          & ion_rad,i,newID,ionID,particles_emitter(i))
          !$OMP END CRITICAL

          nrIonizations = nrIonizations + 1

          cycle atom
        end if
      end if
    end do elec
  end do atom
  !$OMP END PARALLEL DO
end subroutine Do_Discrete_Ionization_tts

subroutine Do_Discrete_Ionization_old(step,nrIonizations)
  integer, intent(in)               ::  step
  integer, intent(out)              ::  nrIonizations
  ! Parameters
  double precision, parameter       ::  elec_max_speed2 = (2.0d0*q_0*5000.0d0/m_0)
  ! Atom, ion, and electron
  double precision, dimension(1:3)  ::  atom_cur_pos,elec_cur_pos,elec_next_pos,elec_cur_vel,elec_cur_acc
  double precision, dimension(1:3)  ::  ejec_elec_pos, ejec_elec_vel, ion_pos, ion_vel, relative_pos,direct_vec
  double precision                  ::  elec_cur_speed, elec_energy, cur_dist2, ionization_dist
  ! Ionization
  logical                           ::  ionization
  double precision                  ::  ion_rad2, ion_rad, E1, E2, collE, ejecE, inSpeed, outSpeed, newSpeed
  integer                           ::  newID, ionID
  ! Polynomial solver
  double precision                  ::  a,b,cc,dd,e,t,t1_r,t1_i,t2_r,t2_i,t3_r,t3_i,t4_r,t4_i
  complex(8)                        ::  root1,root2,root3,root4
  integer                           ::  code
  ! Misc
  integer                           ::  i,k
  double precision                  ::  rnd

  
  nrIonizations = 0
  if (nrElec == 0) return
  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP& PRIVATE(atom_cur_pos,elec_cur_pos,elec_next_pos,elec_cur_vel,elec_cur_acc) &
  !$OMP& PRIVATE(ejec_elec_pos,ejec_elec_vel,ion_pos,ion_vel,relative_pos,direct_vec) &
  !$OMP& PRIVATE(elec_cur_speed,elec_energy,cur_dist2,ionization_dist,ionization,ion_rad2,ion_rad) &
  !$OMP& PRIVATE(E1,E2,collE,ejecE,inSpeed,outSpeed,newSpeed,newID,ionID) &
  !$OMP& PRIVATE (a,b,cc,dd,e,t,t1_r,t1_i,t2_r,t2_i,t3_r,t3_i,t4_r,t4_i,root1,root2,root3,root4,code,i,k,rnd) &
  !$OMP& SHARED(step,nrPart,time_step,nrID,ion_life_time,particles_emitter,particles_step) &
  !$OMP& SHARED(particles_species,particles_mask,particles_cur_pos,particles_cur_vel,particles_cur_accel) &
  !$OMP& SHARED(particles_cur_energy, particles_ion_cross_rad) &
  !$OMP& REDUCTION(+:nrIonizations)

  ! Go through all particles to find atoms
  do i = 1,nrPart
    if ((particles_species(i) == species_atom) .and. (particles_mask(i) .eqv. .true.)) then
      atom_cur_pos = particles_cur_pos(:,i)
      
      ! Go through all particles to find electrons
      do k = 1,nrPart
        if ((particles_species(k) == species_elec) .and. (particles_mask(k) .eqv. .true.)) then
          
          ! Fetch electron data
          elec_cur_pos = particles_cur_pos(:,k)
          elec_cur_vel = particles_cur_vel(:,k)
          elec_cur_acc = particles_cur_accel(:,k)
          elec_cur_speed = norm2(elec_cur_vel)
          elec_energy = particles_cur_energy(k)
          
          relative_pos = elec_cur_pos - atom_cur_pos
          cur_dist2 = dot_product(relative_pos,relative_pos)

          ! Check ionization only if electron has enough energy to ionize
          ! We have a lot to do, lets not waste time
          if (elec_energy >= N_bind) then
            ion_rad = particles_ion_cross_rad(k)
            ion_rad2 = ion_rad**2

            ionization = .false.

            ! Check if ionization happens
            if (cur_dist2 <= ion_rad2) then
              ionization = .true.
              t = 0
              ionization_dist = norm2(elec_cur_pos - atom_cur_pos)
            else
              a = 0.25d0*(dot_product(elec_cur_acc,elec_cur_acc))
              b = dot_product(elec_cur_vel,elec_cur_acc)
              cc = dot_product(elec_cur_vel,elec_cur_vel) + dot_product(relative_pos,elec_cur_acc)
              dd = 2.0d0*dot_product(relative_pos,elec_cur_vel)
              e = cur_dist2 - ion_rad2

              call SolvePolynomial(a,b,cc,dd,e,code,root1,root2,root3,root4)

              if (code /= 44 .and. code /= 23) then
                t1_r = root1%re
                t1_i = root1%im
                t2_r = root2%re
                t2_i = root2%im
                t3_r = root3%re
                t3_i = root3%im
                t4_r = root4%re
                t4_i = root4%im

                if (code == 31) then
                  if (t1_i == 0 .and. t1_r > 0 .and. t1_r <= time_step) then
                    ionization = .true.
                    t = t1_r
                  else if (t2_i == 0 .and. t2_r > 0 .and. t2_r <= time_step) then
                    ionization = .true.
                    t = t2_r
                  else if (t3_i == 0 .and. t3_r > 0 .and. t3_r <= time_step) then
                    ionization = .true.
                    t = t3_r
                  else if (t4_i == 0 .and. t4_r > 0 .and. t3_r <= time_step) then
                    ionization = .true.
                    t = t4_r
                  end if            
                else if (code == 42) then 
                  if (t1_i == 0 .and. t1_r > 0 .and. t1_r <= time_step) then
                    ionization = .true.
                    t = t1_r
                  else if (t2_i == 0 .and. t2_r > 0 .and. t2_r <= time_step) then
                    ionization = .true.
                    t = t2_r
                  end if
                end if
              end if
            end if

            ! Do ionization
            if (ionization .eqv. .true.) then
              ! Some data for the file
              elec_next_pos = elec_cur_pos + elec_cur_vel*t+0.5d0*elec_cur_acc*t**2
              ionization_dist = norm2(elec_next_pos - atom_cur_pos)
              
              ! --------------------- Energy conservation ---------------------
              E1 = elec_energy ! Starting energy
              E2 = E1 - N_bind ! Energy after ionization
              call random_number(rnd)
              collE = E2*rnd ! Energy of the colliding electron
              ejecE = E2 - collE ! Energy of the ejected electron

              ! --------------------- Colliding electron ----------------------
              ! New direction
              direct_vec = Get_Injected_Vec(elec_energy, elec_cur_vel)
              direct_vec = direct_vec / norm2(direct_vec)
              ! New velocity vector
              particles_cur_vel(:, k) = direct_vec*sqrt(2.0d0*q_0*collE/m_0)
              call Update_Collision_Data(k)
              ! Some variables for the file
              inSpeed = elec_cur_speed
              outSpeed = norm2(particles_cur_vel(:,k))

              ! --------------------- Ejected electron ---------------------
              
              ! New position
              call random_number(ejec_elec_pos)
              ejec_elec_pos = ejec_elec_pos - 0.5d0
              ejec_elec_pos = atom_cur_pos + ejec_elec_pos*length_scale
              ! New direction
              direct_vec = Get_Ejected_Vec(elec_energy, elec_energy, elec_cur_vel)
              direct_vec = direct_vec / norm2(direct_vec)
              ! New velocity vector
              ejec_elec_vel = direct_vec*sqrt(2.0d0*q_0*ejecE/m_0)
              ! Some data for the file
              newSpeed = norm2(ejec_elec_vel)

              ! ----------------------- Created ion -----------------------
              ion_pos = atom_cur_pos
              ion_vel = 0.0d0

              ! --------------------- Ionization ---------------------------
              !$OMP CRITICAL
              ! Add the new electron to the system
              newID = nrID
              call Add_Particle(ejec_elec_pos, ejec_elec_vel, species_elec, step, ion_emitter, -1) ! Electron
              call Update_Collision_Data(newID)
              !$OMP END CRITICAL

              ! Remove the atom from the system
              call Mark_Particles_Remove(i,remove_ion)

              !$OMP CRITICAL
              ! Add the new positively charged ion to the system
              ionID = nrID
              call Add_Particle(ion_pos, ion_vel, species_ion, step, ion_emitter, step+ion_life_time) ! Ion
              call Write_Ionization_Data(step,atom_cur_pos,inSpeed,outSpeed,newSpeed,ionization_dist, &
              & ion_rad,i,newID,ionID,particles_emitter(i))
              !$OMP END CRITICAL

              nrIonizations = nrIonizations + 1

              exit
            end if
          end if
        end if
      end do
    end if
  end do
  !$OMP END PARALLEL DO
end subroutine Do_Discrete_Ionization_old

! -------------------------------------------------------------------------------------
! --------------------------- RECOM CROSS SECTION -------------------------------------
! -------------------------------------------------------------------------------------

function Calculate_Kramers_Cross_Section(energy)
  double precision              ::  Calculate_Kramers_Cross_Section
  double precision, intent(in)  ::  energy

  Calculate_Kramers_Cross_Section = 2.105d-26 * (Ryd**2) * (Z_eff**4) / ( N_n*energy * ((N_n**2)*energy + Ryd*Z_eff**2) )
end function Calculate_Kramers_Cross_Section

! -------------------------------------------------------------------------------------

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

! -------------------------------------------------------------------------------------

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

! subroutine Do_Ion_Collisions(step)
!   integer, intent(in)              :: step
!   double precision, dimension(1:3) :: cur_pos, prev_pos, par_vec, old_vel
!   double precision, parameter      :: v2_min      = (2.0d0*q_0*0.1d0/m_0) ! Minimum velocity squared
!   double precision, parameter      :: v2_max      = (2.0d0*q_0*5000.0d0/m_0) ! Maximum velocity squared
!   double precision                 :: mean_path, mean_path_avg, mean_actual_avg ! Mean free path
!   integer                          :: count_n
!   double precision                 :: cross_tot, cross_ion
!   double precision                 :: d ! The distance traveled
!   double precision                 :: rnd, alpha
!   double precision                 :: vel2             ! Squared velocity of the current particle
!   double precision, parameter      :: e_max = 0.10d0    ! Max value of the coefficient of restitution
!   integer                          :: i, nrColl, IFAIL, nrIon
!   double precision                 :: KE ! Kinetic energy
!   double precision                 :: cur_time
!   integer                          :: life_time

!   nrColl = 0
!   nrIon = 0
!   count_n = 0
!   mean_path_avg = 0.0d0
!   mean_actual_avg = 0.0d0

!   if (nrPart > 0) then

!   !$OMP PARALLEL DO DEFAULT(NONE) &
!   !$OMP& PRIVATE(i, cur_pos, prev_pos, d, alpha, rnd, par_vec, vel2, KE, mean_path, cross_tot, cross_ion) &
!   !$OMP& PRIVATE(old_vel, life_time) &
!   !$OMP& SHARED(nrPart, particles_species, particles_life, step, particles_cur_pos) &
!   !$OMP& SHARED(particles_prev_pos, particles_cur_vel, time_step, n_d, ion_life_time) &
!   !$OMP& REDUCTION(+:nrColl, count_n, mean_path_avg, nrIon)
!   do i = 1, nrPart
!     if (particles_species(i) /= species_elec) then
!       if (step >= particles_life(i)) then
!         call Mark_Particles_Remove(i, remove_top) ! Mark the ion to be removed, it has reached the end of its life.
!       end if
!     else

!       cur_pos(:) = particles_cur_pos(:, i)
!       prev_pos(:) = particles_prev_pos(:, i)
!       !prev_pos(:) = particles_last_col_pos(:, i)
!       old_vel = particles_cur_vel(:, i)

!       d = sqrt( (cur_pos(1) - prev_pos(1))**2 + (cur_pos(2) - prev_pos(2))**2 + (cur_pos(3) - prev_pos(3))**2 )
!       vel2 = old_vel(1)**2 + old_vel(2)**2 + old_vel(3)**2
    
!       ! The velocity should be above the minimum velocity we set.
!       !if ((vel2 > v2_min) .and. (vel2 < v2_max)) then
!       if (vel2 > v2_min) then

!         if (vel2 > v2_max) then
!           KE = 0.5d0*m_0*v2_max/q_0
!         else
!           KE = 0.5d0*m_0*vel2/q_0 ! Energy in eV
!         end if

!         cross_tot = Find_Cross_tot_data(KE)

!         mean_path = 1.0d0/(n_d*cross_tot)
!         mean_path_avg = mean_path_avg + mean_path
!         count_n = count_n + 1

!         ! Calculate the ratio between the distance travled and the mean free path
!         alpha = d/mean_path
!         ! if (alpha > 1.0d0) then
!         !   print *, 'WARNING: alpha > 1 in mean path'
!         !   print *, mean_path
!         !   print *, mean_path/1.0E-9
!         !   print *, n_d
!         !   print *, cross_tot
!         !   print *, KE
!         !   print *, sqrt(vel2)
!         !   print *, d
!         !   !pause
!         ! end if

!         ! Check if we do a collision or not
!         call random_number(rnd)
!         if (rnd < alpha) then
!           ! Pick a new random direction for the particle
!           par_vec = Get_Injected_Vec(KE, old_vel)
!           !call random_number(par_vec)
!           !par_vec(1:2) = par_vec(1:2) - 0.5d0
!           !par_vec(3) = par_vec(3) - 0.25d0
!           !par_vec = par_vec - 0.5d0
!           par_vec = par_vec / sqrt(par_vec(1)**2 + par_vec(2)**2 + par_vec(3)**2)

!           ! Set the new velocity
!           call random_number(rnd)

!           particles_cur_vel(:, i) = par_vec*sqrt(vel2)*rnd*1.0d0

!           call random_number(rnd)
!           !par_vec = par_vec*sqrt(vel2)*(1.0d0 - rnd)

!           rnd = rnd*0.80d0 + (1.0d0-0.80d0)
!           par_vec = par_vec*sqrt(2.0d0*q_0*rnd*40.0d0/m_0) ! Give the electron an initial energy between 8 and 40 eV

!           !prev_pos(:) = particles_last_col_pos(:, i)
!           !d = sqrt( (cur_pos(1) - prev_pos(1))**2 + (cur_pos(2) - prev_pos(2))**2 + (cur_pos(3) - prev_pos(3))**2 )
!           !particles_last_col_pos(:, i) = cur_pos(:)
!           !mean_actual_avg = mean_actual_avg + d

!           !---------------------------------------------
!           ! Check if the collision ionizes the N2 or not  
!           cross_ion = Find_Cross_ion_data(KE)
          
!           alpha = cross_ion/cross_tot
!           !alpha = 0.20d0
!           ! if (alpha > 1.0d0) then
!           !   print *, 'WARNING: alpha > 1 in cross section'
!           ! end if
          
!           call random_number(rnd)
!           if (rnd < alpha) then ! Check if we ionize or not

!             ! We ionize
!             ! Pick a position for the new electron to appear at some where close to where the collision occurred
!             call random_number(prev_pos)
!             prev_pos = prev_pos - 0.5d0
!             cur_pos = cur_pos + prev_pos*length_scale

!             ! Pick a direction for the new electron to go in
!             par_vec = Get_Ejected_Vec(KE, KE, old_vel)

!             !$OMP CRITICAL
!             ! Add the new electron to the system
!             call Add_Particle(cur_pos, par_vec, species_elec, step, 1, -1) ! Electron
!             !$OMP END CRITICAL

!             par_vec = 0.0d0

!             ! Pick a position for the ion to appear at some where close to where the collision occurred
!             call random_number(prev_pos)
!             prev_pos = prev_pos - 0.5d0
!             cur_pos = cur_pos + prev_pos*length_scale

!             ! Calculate the life time of the Ion
!             ! This number is drawn from an exponential distribution
!             ! The half life is 0.13843690559395497 ps
!             !alpha = -0.13843690559395497E-12 / time_step
!             !call random_number(rnd)
!             !life_time = NINT(alpha*log(1.0d0 - rnd))
!             life_time = ion_life_time

!             !$OMP CRITICAL
!             ! Add the new positively charged ion to the system
!             call Add_Particle(cur_pos, par_vec, species_ion, step, 1, step+life_time) ! Ion
!             !$OMP END CRITICAL

!             nrIon = nrIon + 1
!           end if

!           ! Update the number of collisions
!           nrColl = nrColl + 1
!         end if
!       else
!         KE = 0.5d0*m_0*vel2/q_0 ! Energy in eV
!         ! if (KE > 0.1d0) then
!         !   print *, 'Outside range'
!         !   print *, KE
!         !   print *, ''
!         ! end if
!       end if
!     end if
!   end do
!   !$OMP END PARALLEL DO

!   end if

!   if (count_n > 1) then
!     mean_path_avg = mean_path_avg / count_n
!   end if
!   mean_actual_avg = 0.0d0

!   cur_time = time_step * step / time_scale ! Scaled in units of time_scale

!   ! Write data
!   write(ud_coll, '(i6, tr2, ES12.4, tr2, i6, tr2, i6, tr2, i6, tr2, ES12.4, tr2, ES12.4)', iostat=IFAIL) &
!            step, cur_time, nrColl, nrIon, count_n, &
!            (mean_path_avg/length_scale), (mean_actual_avg/length_scale)
! end subroutine Do_Ion_Collisions

! function Get_Injected_Vec(T, par_vel)
!   double precision, dimension(1:3)             :: Get_Injected_Vec
!   double precision, dimension(1:3), intent(in) :: par_vel
!   double precision, intent(in)                 :: T ! Energy in eV
!   double precision, parameter                  :: mu = 5.0d0, sigma = 25.0d0
!   double precision                             :: angle
!   double precision                             :: m_factor, rnd, alpha
!   double precision                             :: dot_p, len_vec, len_vel
!   double precision, dimension(1:3)             :: par_vec
!   integer                                      :: n_tries

!   !m_factor = 1.0d0/(sqrt(2.0d0*pi)*sigma)
!   m_factor = folded_normal_dist(mu, sigma, mu)
!   n_tries = 0

!   do
!     ! Get random numbers between 0 and 1
!     call random_number(par_vec)
!     ! Convert this to random numbers between -0.5 and 0.5
!     par_vec = par_vec - 0.5d0
!     ! Dot product
!     dot_p = par_vel(1)*par_vec(1) + par_vel(2)*par_vec(2) + par_vel(3)*par_vec(3)

!     ! Length of the vectors
!     len_vec = sqrt(par_vec(1)**2 + par_vec(2)**2 + par_vec(3)**2)
!     len_vel = sqrt(par_vel(1)**2 + par_vel(2)**2 + par_vel(3)**2)

!     ! Calculate the angle between the two vector using the dot product.
!     ! We want the angle in degrees.
!     angle = acos(dot_p/(len_vec*len_vel)) * 180.0d0/pi

!     ! Get the acceptance/rejection criteria. We use a folded normal distribution
!     ! because the angle between two vector is always positive. [0, pi]/[0, 180]. 
!     alpha = folded_normal_dist(mu, sigma, angle) / m_factor

!     call random_number(rnd)
!     if (rnd < alpha) then
!       exit ! Exit the loop
!     else
!       n_tries = n_tries + 1
!       if (n_tries >= 1000000) then ! This should never take this long
!         print *, 'n_tries > 100000 (Injected)'
!         print *, 'angle = ', angle
!         print *, 'alpha = ', alpha
!         print *, 'n_tries = ', n_tries
!         print *, ''
!         exit
!       end if
!     end if
!   end do

!   ! Return the normalized vector
!   Get_Injected_Vec = (par_vec / len_vec)
! end function Get_Injected_Vec

! ! This function returns a normalized direction vector for the ejected electron
! ! The angle between this vector and the velocity vector of the incident electron
! ! is approximatly distrubted according the experimental values.
! ! See:
! ! Dobly Differential Cross Section for Electron Scattered by Nitrogen
! ! J. C. Nogueira, M. A. Eschiapati Ferreira and Ronaldo S. Barbieri
! function Get_Ejected_Vec(W, T, par_vel)
!   double precision, dimension(1:3)             :: Get_Ejected_Vec
!   double precision, intent(in)                 :: T, W
!   double precision, dimension(1:3), intent(in) :: par_vel
!   double precision, parameter                  :: a = -430.5d0, b = -0.5445d0, c = 89.32d0
!   double precision, parameter                  :: sigma = 48.0d0
!   double precision                             :: angle_max, angle
!   double precision                             :: m_factor, rnd, alpha
!   double precision                             :: dot_p, len_vec, len_vel
!   double precision, dimension(1:3)             :: par_vec
!   integer                                      :: n_tries 

!   if (T < 100.d0) then
!     angle_max = a*100.0d0**b + c
!   else
!     angle_max = a*T**b + c
!   end if

!   !m_factor = 1.0d0/(sqrt(2.0d0*pi)*sigma)
!   m_factor = folded_normal_dist(angle_max, sigma, angle_max)
!   n_tries = 0

!   do
!     call random_number(par_vec)
!     par_vec = par_vec - 0.5d0
!     dot_p = par_vel(1)*par_vec(1) + par_vel(2)*par_vec(2) + par_vel(3)*par_vec(3)
!     len_vec = sqrt(par_vec(1)**2 + par_vec(2)**2 + par_vec(3)**2)
!     len_vel = sqrt(par_vel(1)**2 + par_vel(2)**2 + par_vel(3)**2)

!     angle = acos(dot_p/(len_vec*len_vel)) * 180.0d0/pi

!     alpha = folded_normal_dist(angle_max, sigma, angle) / m_factor

!     call random_number(rnd)
!     if (rnd < alpha) then
!       exit ! Exit the loop
!     else
!       n_tries = n_tries + 1
!       if (n_tries >= 1000000) then ! This should never take this long
!         print *, 'n_tries > 100000 (Ejected)'
!         exit
!       end if
!     end if
!   end do

!   Get_Ejected_Vec = (par_vec / len_vec)
! end function Get_Ejected_Vec

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
    print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file ', filename_tot
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
    print '(a)', 'RUMDEED: ERROR UNABLE TO OPEN file ', filename_ion
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

! --------------------------------------------------------------------------
! Distribute atoms uniformly inside the cylinder for discrete ionization
! --------------------------------------------------------------------------

subroutine Place_Atoms_Cylinder()
  double precision                  :: volume, rnd, cylindR, cylindR2, R2
  double precision, dimension(1:3)  :: atom_pos, atom_vel
  integer                           :: nrStartAtoms, nrStartIons, i

  cylindR = emitters_dim(1, 1) ! Radius of the gas cylinder
  cylindR2 = cylindR**2 ! Radius of cylinder squared
  
  ! Calculate number of atoms
  volume = pi*cylindR2*box_dim(3)
  nrStartAtoms = int(volume*n_d)
  nrStartIons = int(nrStartAtoms*ion_atom_ratio)

  print *, 'Max particles = ', MAX_PARTICLES
  print *, 'Nr atoms = ', nrStartAtoms
  print *, 'Nr ions = ', nrStartIons
  
  if ((collision_mode == 3) .or. (collision_mode == 4)) then
  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP& PRIVATE(i, rnd, atom_pos, atom_vel, R2) &
  !$OMP& SHARED(cylindR, cylindR2, nrStartAtoms, nrStartIons, box_dim, emitters_pos, emitters_dim, ion_life_time, collision_mode)

  ! Place atoms
  do i=1,nrStartAtoms
    ! Find random position in the cylinder
    R2 = cylindR2 + 1.0d0 ! Must be larger than emitR2 for the do while loop to run
    do while (R2 > cylindR2)
      call random_number(atom_pos(1:2)) ! Gives a random number [0,1]
      atom_pos(1:2) = 2.0d0*(atom_pos(1:2) - 0.5d0)*cylindR ! Range is -cylindR to +cylindR
      R2 = atom_pos(1)**2 + atom_pos(2)**2 ! Radius squared of our random point
    end do
    atom_pos(1:2) = emitters_pos(1:2,1) + atom_pos(1:2) + emitters_dim(1:2,1)
    call random_number(rnd)
    atom_pos(3) = box_dim(3)*rnd

    ! Set velocity to zero
    atom_vel(:) = 0.0d0

    if (nrStartIons > 0) then
      !$OMP CRITICAL
      call Add_Particle(atom_pos, atom_vel, species_ion, 0, ion_emitter, ion_life_time)
      !$OMP END CRITICAL
      nrStartIons = nrStartIons - 1
    else
      !$OMP CRITICAL
      call Add_Particle(atom_pos, atom_vel, species_atom, 0, ion_emitter, -1)
      !$OMP END CRITICAL
    end if
  end do
  !$OMP END PARALLEL DO
  end if
end subroutine Place_Atoms_Cylinder

! --------------------------------------------------------------------------
! Update electron energy, ionization and recombination cross section radii
! --------------------------------------------------------------------------

subroutine Update_Collision_Data(i)
  integer, intent(in)         ::  i
  double precision, parameter ::  elec_max_speed2 = (2.0d0*q_0*5000.0d0/m_0)
  double precision            ::  elec_cur_speed2, elec_energy
  double precision            ::  ion_cross_sec, ion_cross_rad, recom_cross_sec, recom_cross_rad, tot_cross_sec


  elec_cur_speed2 = norm2(particles_cur_vel(:,i))**2
  
  if (elec_cur_speed2 > elec_max_speed2) then
    elec_energy = 0.5d0*m_0*elec_max_speed2/q_0
  else
    elec_energy = 0.5d0*m_0*elec_cur_speed2/q_0
  end if
  ion_cross_sec = Find_Cross_ion_data(elec_energy)
  ion_cross_rad = ion_cross_sec/pi
  tot_cross_sec = Find_Cross_tot_data(elec_energy)

  elec_energy = 0.5d0*m_0*elec_cur_speed2/q_0
  recom_cross_sec = Calculate_Kramers_Cross_Section(elec_energy)
  recom_cross_rad = sqrt(Calculate_Kramers_Cross_Section(elec_energy)/pi)

  particles_cur_energy(i) = elec_energy
  particles_ion_cross_sec(i) = ion_cross_sec
  particles_ion_cross_rad(i) = ion_cross_rad
  particles_recom_cross_rad(i) = recom_cross_rad
  particles_tot_cross_sec(i) = tot_cross_sec

end subroutine Update_Collision_Data

! --------------------------------------------------------------------------
! Update collision data for all electrons
! --------------------------------------------------------------------------

subroutine Update_Collision_Data_All_tts()
  integer ::  i,k

  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP& PRIVATE(i,k) &
  !$OMP& SHARED(particles_species, particles_elec_pointer, nrElec)
  do k=1,nrElec
    ! print *, 'Update collision data for electron ', k
    i = particles_elec_pointer(k)
    call Update_Collision_Data(i)
  end do
  !$OMP END PARALLEL DO

end subroutine Update_Collision_Data_All_tts

subroutine Update_Collision_Data_All_ots()
  integer ::  i

  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP& PRIVATE(i) &
  !$OMP& SHARED(particles_species, particles_elec_pointer, nrPart)
  do i=1,nrPart
    ! print *, 'Update collision data for electron ', i
    if (particles_species(i) == species_elec) then
      call Update_Collision_Data(i)
    end if
  end do
  !$OMP END PARALLEL DO

end subroutine Update_Collision_Data_All_ots

end module mod_collisions
