!-------------------------------------------!
! Module for image charge                   !
! Kristinn Torfason                         !
! 16.05.15                                  !
!-------------------------------------------!

Module mod_ic
  use mod_global
  use mod_hyperboloid_tip
  use mod_verlet
  implicit none
contains
  function Accel_on_Tip(pos)
    double precision, dimension(1:3)             :: Accel_on_Tip
    double precision, dimension(1:3), intent(in) :: pos
    double precision, dimension(1:3)             :: pos_a, pos_c, pos_b, n, n_a
    double precision                             :: a, b, phi_a, phi_c
    double precision                             :: xi, R_sphere, dis_a, dis_b, tmp_fac
    double precision, dimension(1:3)             :: E, E_1, E_2
    integer                                      :: i
    double precision, parameter :: pre_fac_E = (-1.0d0*q_0) / m_0 ! Electric field prefactor
    double precision, parameter :: pre_fac_a = m_0 / (-1.0d0*q_0) ! Acceleration prefactor
    

    Accel_on_Tip = 0.0d0

    if (image_charge .eqv. .false.) then

      !call Calculate_Acceleration_particle(par_elec)
      !Accel_on_Tip = par_elec%accel !+ field_E_Hyperboloid(pos)
      Accel_on_Tip = pre_fac_E*Calc_Field_at(pos)

      return
    end if

    if (nrElec > 0) then
      xi = xi_coor(pos(1), pos(2), pos(3))

      tmp_fac = eta_1 / ( sqrt(1.0d0 - eta_1**2) * sqrt(pos(1)**2 + pos(2)**2 + a_foci**2*(1.0d0-eta_1**2)) )
      n(1) = pos(1) * tmp_fac
      n(2) = pos(2) * tmp_fac
      n(3) = -1.0d0

      !n(1) = eta_1 / sqrt(1.0d0 - eta_1**2)*pos(1)/sqrt(pos(1)**2+pos(2)**2+a_foci**2*(1.0d0-eta_1**2))
      !n(2) = eta_1 / sqrt(1.0d0 - eta_1**2)*pos(2)/sqrt(pos(1)**2+pos(2)**2+a_foci**2*(1.0d0-eta_1**2))
      !n(3) = -1.0d0
      n = n / norm2(n)

      R_sphere = abs(a_foci/eta_1 * ( (xi**2 - eta_1**2)**(3.0d0/2.0d0) ) / (sqrt(1.0d0 - eta_1**2)))

      pos_c = n*R_sphere + pos

      do i = 1, nrElec
        pos_a = particles_cur_pos(:, i)

        a = norm2(pos_a - pos_c)
        !if (a < 10.0d0*length) then
          phi_a = phi_coor(pos_a(1), pos_a(2), pos_a(3))
          phi_c = phi_coor(pos(1), pos(2), pos(3))

          !if (abs(phi_a - phi_c) < pi) then
            b = R_sphere**2/a

            n_a = (pos_a - pos_c) / a
            pos_b = pos_c + n_a*b

            dis_a = norm2(pos - pos_a)
            dis_b = norm2(pos - pos_b)

            E_1 = (-1.0d0*q_0)/(4.0d0*pi*epsilon_0) * ( pos - pos_a ) / dis_a**3
            E_2 = (-1.0d0*q_0)/(4.0d0*pi*epsilon_0) * ( R_sphere*( pos - pos_b ) ) / (a*dis_b**3)
            E = E_1 - E_2

            Accel_on_Tip = Accel_on_Tip + E
          !end if
        !end if
      end do

    end if

    Accel_on_Tip = pre_fac_E*( Accel_on_Tip + field_E_Hyperboloid(pos) )
    !print *, 'A = ', Accel_on_Tip

  end function Accel_on_Tip
end module mod_ic
