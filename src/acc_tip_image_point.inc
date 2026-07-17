! ----------------------------------------------------------------------------
! Location of the image charge of a particle at (x_a, y_a, z_a) in the
! sphere at the top of the hyperboloid tip: the first half of Sphere_IC_field
! (mod_hyperboloid_tip, see there for the physics). The sphere has center
! (0, 0, z_0_tip) and radius Rsp_tip. Inlined into the OpenACC kernels and
! shared between Calculate_Acceleration_Particles_ACC and Calc_Field_at_Batch
! via #include; acc_tip_ic_force.inc consumes the outputs.
!
! Inputs:  x_a, y_a, z_a      position of the particle being imaged
!          z_0_tip, Rsp_tip   sphere center height and radius
! Outputs: dis_a              distance of the particle from the sphere center
!          x_im, y_im, z_im   position of the image charge

dis_a = sqrt(x_a**2 + y_a**2 + (z_a - z_0_tip)**2)

z_im = z_0_tip + Rsp_tip**2 / ( sqrt(1 + x_a**2/(z_a - z_0_tip)**2 + y_a**2/(z_a - z_0_tip)**2) * dis_a )
x_im = (z_im - z_0_tip) * x_a / (z_a - z_0_tip)
y_im = (z_im - z_0_tip) * y_a / (z_a - z_0_tip)
