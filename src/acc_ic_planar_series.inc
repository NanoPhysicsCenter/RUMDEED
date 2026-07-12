! ----------------------------------------------------------------------------
! Planar image charge partner series (the body of Force_Image_charges_v2 in
! mod_verlet, see that function for the physics). Inlined into the OpenACC
! kernels because procedure calls do not perform inside them; shared between
! Calculate_Acceleration_Particles_ACC and Calc_Field_at_Batch via #include
! so the formula only exists in one place.
!
! Inputs:  z_a            evaluation height
!          z_b            height of the source particle (its partners are at
!                         -z_b and 2*n*d_loc -+ z_b)
!          diff_x, diff_y lateral offsets (evaluation - source), the partners
!                         share the x and y coordinates of their source
!          dxy2           diff_x**2 + diff_y**2
!          Nic, d_loc     N_ic_max and the gap spacing d
! Outputs: ic_x, ic_y, ic_z  the summed contribution, without the charge
!                            prefactor (the caller applies pre_fac_c)
! Clobbers: z_ic, diff_z, r, inv_r3, n

! n = 0: Opposite charge partner below the bottom plane
z_ic = -1.0d0*z_b
diff_z = z_a - z_ic
r = sqrt( dxy2 + diff_z**2 ) + length_scale**2
inv_r3 = 1.0d0 / (r*r*r)
ic_x = -diff_x*inv_r3
ic_y = -diff_y*inv_r3
ic_z = -diff_z*inv_r3

!$acc loop seq
do n = 1, Nic
  ! Opposite charge partners, plus and minus n
  z_ic = 2.0d0*n*d_loc - z_b
  diff_z = z_a - z_ic
  r = sqrt( dxy2 + diff_z**2 ) + length_scale**2
  inv_r3 = 1.0d0 / (r*r*r)
  ic_x = ic_x - diff_x*inv_r3
  ic_y = ic_y - diff_y*inv_r3
  ic_z = ic_z - diff_z*inv_r3

  z_ic = -2.0d0*n*d_loc - z_b
  diff_z = z_a - z_ic
  r = sqrt( dxy2 + diff_z**2 ) + length_scale**2
  inv_r3 = 1.0d0 / (r*r*r)
  ic_x = ic_x - diff_x*inv_r3
  ic_y = ic_y - diff_y*inv_r3
  ic_z = ic_z - diff_z*inv_r3

  ! Same charge partners, plus and minus n
  z_ic = 2.0d0*n*d_loc + z_b
  diff_z = z_a - z_ic
  r = sqrt( dxy2 + diff_z**2 ) + length_scale**2
  inv_r3 = 1.0d0 / (r*r*r)
  ic_x = ic_x + diff_x*inv_r3
  ic_y = ic_y + diff_y*inv_r3
  ic_z = ic_z + diff_z*inv_r3

  z_ic = -2.0d0*n*d_loc + z_b
  diff_z = z_a - z_ic
  r = sqrt( dxy2 + diff_z**2 ) + length_scale**2
  inv_r3 = 1.0d0 / (r*r*r)
  ic_x = ic_x + diff_x*inv_r3
  ic_y = ic_y + diff_y*inv_r3
  ic_z = ic_z + diff_z*inv_r3
end do
