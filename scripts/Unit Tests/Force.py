import numpy as np

length_scale = 1.0E-9 # 1 nm

m_0 = 9.10938356E-31
q_0 = 1.6021766208E-19
pi = 3.14159265358979324 # Pi
c = 299792458.0 # Speed of light (m/s)
mu_0 = 4.0*pi * 1.0E-7 # Vacuum permeability (H/m)
epsilon_0 = 1.0/(mu_0 * c**2) # Vacuum permittivity (F/m)
pre_fac_c = (-q_0)*(-q_0) / (4 * pi * epsilon_0)

pos_1 = np.array([1.0, -2.34, 1.0]) * length_scale
pos_2 = np.array([1.0, -2.34, -1.0]) * length_scale
pos_3 = np.array([1.0, -2.34, 1999.0]) * length_scale

pos_4 = np.array([-4.4, -1.44, 2.0]) * length_scale
pos_5 = np.array([-4.4, -1.44, -2.0]) * length_scale
pos_6 = np.array([-4.4, -1.44, 1998.0]) * length_scale

# Acceleration on particle 1 from others
# Particle 2
diff = pos_1 - pos_2
r = np.sqrt( np.sum(diff**2) ) + length_scale**3 # + length_scale**3 to avoid division with zero
force_2 = pre_fac_c * diff / r**3

# Particle 3
diff = pos_1 - pos_3
r = np.sqrt( np.sum(diff**2) ) + length_scale**3 # + length_scale**3 to avoid division with zero
force_3 = pre_fac_c * diff / r**3

# Particle 4
diff = pos_1 - pos_4
r = np.sqrt( np.sum(diff**2) ) + length_scale**3 # + length_scale**3 to avoid division with zero
force_4 = pre_fac_c * diff / r**3

# Particle 5
diff = pos_1 - pos_5
r = np.sqrt( np.sum(diff**2) ) + length_scale**3 # + length_scale**3 to avoid division with zero
force_5 = pre_fac_c * diff / r**3

# Particle 6
diff = pos_1 - pos_2
r = np.sqrt( np.sum(diff**2) ) + length_scale**3 # + length_scale**3 to avoid division with zero
force_6 = pre_fac_c * diff / r**3

# Total force

force_TOT = force_2 + force_3 + force_4 + force_5 + force_6
accel_TOT = force_TOT/m_0
print(accel_TOT)