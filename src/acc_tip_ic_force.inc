! ----------------------------------------------------------------------------
! Field of the particle at (x_a, y_a, z_a) and of its image charge in the
! sphere at the top of the hyperboloid tip, evaluated at (x_b, y_b, z_b):
! the second half of Sphere_IC_field (mod_hyperboloid_tip, see there for the
! physics). Expects dis_a and (x_im, y_im, z_im) from acc_tip_image_point.inc.
! Note that unlike the planar series this includes the q_0/(4*pi*epsilon_0)
! factor, exactly like the original function does.
!
! Inputs:  x_a, y_a, z_a          position of the imaged particle
!          x_b, y_b, z_b          evaluation position
!          dis_a, x_im, y_im, z_im  from acc_tip_image_point.inc
!          Rsp_tip                sphere radius
! Outputs: ic_x, ic_y, ic_z       the contribution, without the charge
!                                 prefactor (the caller applies pre_fac_c)
! Clobbers: tmp_dis_a, tmp_dis_b

tmp_dis_a = ( (x_b - x_a)**2 + (y_b - y_a)**2 + (z_b - z_a)**2 )**(3.0d0/2.0d0)
tmp_dis_b = ( (x_b - x_im)**2 + (y_b - y_im)**2 + (z_b - z_im)**2 )**(3.0d0/2.0d0)

ic_x = 1.0d0*q_0/(4.0d0*pi*epsilon_0) * ( (x_a - x_b)/tmp_dis_a - (Rsp_tip*(x_im - x_b))/(dis_a*tmp_dis_b) )
ic_y = 1.0d0*q_0/(4.0d0*pi*epsilon_0) * ( (y_a - y_b)/tmp_dis_a - (Rsp_tip*(y_im - y_b))/(dis_a*tmp_dis_b) )
ic_z = 1.0d0*q_0/(4.0d0*pi*epsilon_0) * ( (z_a - z_b)/tmp_dis_a - (Rsp_tip*(z_im - z_b))/(dis_a*tmp_dis_b) )
