!-------------------------------------------!
! Module for hyperboloid tip                !
! Kristinn Torfason                         !
! 05.04.13                                  !
!-------------------------------------------!

Module mod_hyperboloid_tip
  use mod_global
  implicit none

  !Hyperbolid tip
  double precision            :: a_foci ! a_foci = d + 2.5 Foci
  double precision            :: theta_tip ! theta = arccos(d/a_foci) Radius of curvature of the tip
  double precision            :: r_tip ! r_tip = a_foci*sin(theta)*tan(theta)
  double precision            :: eta_1 ! eta_1 = -cos(theta) defines the surface of the tip
  double precision, parameter :: eta_2 = 0.0d0 ! Defines the anode plane (should be 0)
  double precision            :: max_xi ! Max value for xi variable for the surface (tip height)
  double precision            :: shift_z ! shift_z = a_foci*eta_1*max_xi Shift the z-coordinate upwards. Base of tip at z = 0
  double precision            :: h_tip ! Height of the tip from base
  double precision            :: R_base ! Base radius
  double precision            :: pre_fac_E_tip, pre_fac_E_tip_unit_voltage
  !double precision            :: xi_0 = 1.0124d0 ! 99.9% Emission area

contains
  pure function xyz_corr(xi, eta, phi)
    double precision, intent(in)     :: xi, eta, phi
    double precision, dimension(1:3) :: xyz_corr
    double precision                 :: xy

    xy = a_foci * sqrt( (xi**2 - 1.0d0)*(1.0d0 - eta**2) )
    xyz_corr(1) = xy * cos(phi)
    xyz_corr(2) = xy * sin(phi)
    xyz_corr(3) = a_foci * xi * eta + shift_z
  end function xyz_corr

  double precision elemental function x_coor(xi, eta, phi)
    double precision, intent(in) :: xi, eta, phi

    x_coor = a_foci * sqrt( (xi**2 - 1.0d0)*(1.0d0 - eta**2) ) * cos(phi)
  end function x_coor

  double precision elemental function y_coor(xi, eta, phi)
    double precision, intent(in) :: xi, eta, phi

    y_coor = a_foci * sqrt( (xi**2 - 1.0d0)*(1.0d0 - eta**2) ) * sin(phi)
  end function y_coor

  double precision elemental function z_coor(xi, eta, phi)
    double precision, intent(in) :: xi, eta, phi

    z_coor = a_foci * xi * eta + shift_z
  end function z_coor

  double precision elemental function xi_coor(x, y, z)
    double precision, intent(in) :: x, y, z

    xi_coor = 1.0d0/(2.0d0*a_foci) * ( sqrt(x**2 + y**2 + (z + a_foci - shift_z)**2) &
          & + sqrt(x**2 + y**2 + (z - a_foci - shift_z)**2) )
  end function xi_coor

  double precision elemental function eta_coor(x, y, z)
    double precision, intent(in) :: x, y, z

    eta_coor = 1.0d0/(2.0d0*a_foci) * ( sqrt(x**2 + y**2 + (z + a_foci - shift_z)**2) &
          & - sqrt(x**2 + y**2 + (z - a_foci - shift_z)**2) )
  end function eta_coor

  double precision elemental function phi_coor(x, y, z)
    double precision, intent(in) :: x, y, z

    if ((abs(x) < 1.0d-18) .and. (abs(y) < 1.0d-18)) then
      phi_coor = 0.0d0
    else
      phi_coor = atan2(y, x)
    end if
  end function phi_coor

  pure function surface_normal(par_pos)
    double precision, dimension(1:3)             :: surface_normal
    double precision, dimension(1:3), intent(in) :: par_pos
    double precision                             :: my_xi, eta_fac, div_fac

    !my_xi = xi_coor(par_pos(1), par_pos(2), par_pos(3))

    !surface_normal(1) = -1.0d0 * par_pos(1) * (par_pos(3) - shift_z) / my_xi
    !surface_normal(2) = par_pos(2) * (par_pos(3) - shift_z) / (eta_1)
    !surface_normal(3) = a_foci**2 * my_xi * (1 - eta_1**2)

    !surface_normal = surface_normal / norm2(surface_normal)

    eta_fac = eta_1 / sqrt(1 - eta_1**2)
    div_fac = -1.0d0/sqrt(par_pos(1)**2 + par_pos(2)**2 + a_foci**2*(1-eta_1**2))

    surface_normal(1) = eta_fac * par_pos(1) * div_fac
    surface_normal(2) = eta_fac * par_pos(2) * div_fac
    surface_normal(3) = 1.0d0

    surface_normal = surface_normal / norm2(surface_normal)
  end function surface_normal


  ! Gives the component of the field normal to the surface
  ! In the prolate spheroidal coordinates this would be the \eta component

  pure function Field_normal(par_pos, field_E)
    double precision, dimension(1:3), intent(in) :: par_pos, field_E
    double precision                             :: Field_normal
    double precision, dimension(1:3)             :: unit_vec

    unit_vec = surface_normal(par_pos)
    Field_normal = dot_product(unit_vec, field_E)
  end function Field_normal

  !pure function field_E_Hyperboloid(pos_xyz)
  pure function field_E_Hyperboloid(pos_xyz, org_pos, is_surface)
    double precision, dimension(1:3), intent(in)           :: pos_xyz
    double precision, dimension(1:3), intent(in), optional :: org_pos
    logical, intent(in), optional                          :: is_surface
    double precision, dimension(1:3)                       :: pos_pro !xi, eta, phi
    double precision, dimension(1:3)                       :: field_E_Hyperboloid
    double precision                                       :: pre_fac_E_tip_xyz, fac_E_xy

    pos_pro(1) = xi_coor(pos_xyz(1), pos_xyz(2), pos_xyz(3))
    pos_pro(2) = eta_coor(pos_xyz(1), pos_xyz(2), pos_xyz(3))
    pos_pro(3) = phi_coor(pos_xyz(1), pos_xyz(2), pos_xyz(3))

!     if ((isnan(pos_pro(1)) == .true.) .or. (isnan(pos_pro(2)) == .true.) .or. (isnan(pos_pro(3)) == .true.)) then
!       print *, 'pos_pro is NaN'
!       print *, pos_pro
!       stop
!     end if

    ! Electric field for the hyperboloid geometry
    pre_fac_E_tip_xyz = pre_fac_E_tip * 1.0d0/(pos_pro(1)**2 - pos_pro(2)**2)

    ! if we are near the tip then the x and y components should be zero
    if (abs(pos_pro(1) - 1.0d0) < 1.0d-6) then
      fac_E_xy = 0.0d0
    else
      fac_E_xy = pos_pro(2) * sqrt( (pos_pro(1)**2 - 1.0d0)/(1.0d0 - pos_pro(2)**2) )
    end if

    field_E_Hyperboloid(1) = -1.0d0*pre_fac_E_tip_xyz * fac_E_xy * cos(pos_pro(3))
    field_E_Hyperboloid(2) = -1.0d0*pre_fac_E_tip_xyz * fac_E_xy * sin(pos_pro(3))
    field_E_Hyperboloid(3) = pre_fac_E_tip_xyz * pos_pro(1)

    !print *, field_E_Hyperboloid

!     if ((isnan(field_E_Hyperboloid(1)) == .true.) .or. (isnan(field_E_Hyperboloid(2)) == .true.) .or. (isnan(field_E_Hyperboloid(3)) == .true.)) then
!       print *, 'field_E_Hyperboloid is NaN'
!       print *, field_E_Hyperboloid
!       stop
!     end if
  end function field_E_Hyperboloid

  double precision pure function Tip_Area(xi_1, xi_2, phi_1, phi_2)
    double precision, intent(in) :: xi_1, xi_2, phi_1, phi_2
    double precision             :: fac_1, fac_2

    fac_1 = xi_1*sqrt(xi_1**2 - eta_1**2) - eta_1**2*log(xi_1 + sqrt(xi_1**2 - eta_1**2))
    fac_2 = xi_2*sqrt(xi_2**2 - eta_1**2) - eta_1**2*log(xi_2 + sqrt(xi_2**2 - eta_1**2))
    Tip_Area = 0.5d0*a_foci**2 * sqrt(1.0d0-eta_1**2) * (phi_2 - phi_1) * (fac_2 - fac_1)
  end function Tip_Area
  !-----------------------------------------------------------------------------
  ! Image charge effects of a sphere placed at the top of the tip.
  ! (Moved here from mod_emission_tip so that mod_verlet can reference it for
  ! the OpenACC geometry dispatch.)
  function Sphere_IC_field(pos_1, pos_2)
    double precision, dimension(1:3)             :: Sphere_IC_field
    double precision, dimension(1:3), intent(in) :: pos_1, pos_2
    double precision                 :: x_a, y_a, z_a
    double precision                 :: x_b, y_b, z_b
    double precision                 :: x, y, z
    double precision                 :: z_0, Sphere_R
    double precision                 :: dis_a, tmp_dis_a, tmp_dis_b

    if (image_charge .eqv. .true.) then

      z_0 = h_tip - r_tip
      Sphere_R = r_tip

      !print *, 'Sphere_R = ', Sphere_R
      !print *, 'z_0 = ', z_0

      x_a = pos_1(1)
      y_a = pos_1(2)
      z_a = pos_1(3)

      x = pos_2(1)
      y = pos_2(2)
      z = pos_2(3)

      dis_a = sqrt(x_a**2 + y_a**2 + (z_a - z_0)**2)

      z_b = z_0 + Sphere_R**2 / ( sqrt(1 + x_a**2/(z_a - z_0)**2 + y_a**2/(z_a - z_0)**2) * dis_a )
      x_b = (z_b - z_0) * x_a / (z_a - z_0)
      y_b = (z_b - z_0) * y_a / (z_a - z_0)

      tmp_dis_a = ( (x - x_a)**2 + (y - y_a)**2 + (z - z_a)**2 )**(3.0d0/2.0d0)
      tmp_dis_b = ( (x - x_b)**2 + (y - y_b)**2 + (z - z_b)**2 )**(3.0d0/2.0d0)

      Sphere_IC_field(1) = 1.0d0*q_0/(4.0d0*pi*epsilon_0) * ( (x_a - x)/tmp_dis_a - (Sphere_R*(x_b - x))/(dis_a*tmp_dis_b) )
      Sphere_IC_field(2) = 1.0d0*q_0/(4.0d0*pi*epsilon_0) * ( (y_a - y)/tmp_dis_a - (Sphere_R*(y_b - y))/(dis_a*tmp_dis_b) )
      Sphere_IC_field(3) = 1.0d0*q_0/(4.0d0*pi*epsilon_0) * ( (z_a - z)/tmp_dis_a - (Sphere_R*(z_b - z))/(dis_a*tmp_dis_b) )

    else
      Sphere_IC_field = 0.0d0
    end if

  end function Sphere_IC_field
end Module mod_hyperboloid_tip
