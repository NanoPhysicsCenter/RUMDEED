! ----------------------------------------------------------------------------
! Vacuum electric field of the hyperboloid tip at (x_1, y_1, z_1): the body
! of field_E_Hyperboloid (mod_hyperboloid_tip, see there for the physics),
! written in the prolate spheroidal coordinates xi / eta / phi. Inlined into
! the OpenACC kernels and shared between Calculate_Acceleration_Particles_ACC
! and Calc_Field_at_Batch via #include.
!
! Inputs:  x_1, y_1, z_1      evaluation position
!          af, sz, pre_E_tip  a_foci, shift_z and pre_fac_E_tip of the tip
! Outputs: fE_x, fE_y, fE_z   the vacuum field
! Clobbers: r_p, r_m, xi, eta, phi, pre_fac_xyz, fac_xy

r_p = sqrt(x_1**2 + y_1**2 + (z_1 + af - sz)**2)
r_m = sqrt(x_1**2 + y_1**2 + (z_1 - af - sz)**2)
xi  = (r_p + r_m) / (2.0d0*af)
eta = (r_p - r_m) / (2.0d0*af)
if ((abs(x_1) < 1.0d-18) .and. (abs(y_1) < 1.0d-18)) then
  phi = 0.0d0
else
  phi = atan2(y_1, x_1)
end if

pre_fac_xyz = pre_E_tip * 1.0d0/(xi**2 - eta**2)

! If we are near the tip then the x and y components should be zero
if (abs(xi - 1.0d0) < 1.0d-6) then
  fac_xy = 0.0d0
else
  fac_xy = eta * sqrt( (xi**2 - 1.0d0)/(1.0d0 - eta**2) )
end if

fE_x = -1.0d0*pre_fac_xyz * fac_xy * cos(phi)
fE_y = -1.0d0*pre_fac_xyz * fac_xy * sin(phi)
fE_z = pre_fac_xyz * xi
