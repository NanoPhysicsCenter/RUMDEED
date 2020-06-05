import numpy as np

length_scale = 1.0E-9 # 1 nm

m_0 = 9.10938356E-31
q_0 = 1.6021766208E-19
pi = 3.14159265358979324 # Pi
c = 299792458.0 # Speed of light (m/s)
mu_0 = 4.0*pi * 1.0E-7 # Vacuum permeability (H/m)
epsilon_0 = 1.0/(mu_0 * c**2) # Vacuum permittivity (F/m)
pre_fac_c = q_0*q_0 / (4 * pi * epsilon_0)

V_0 = 1000.0
d = 1000.0E-9

pos_1 = np.array([1.0, -2.34, 1.0]) * length_scale # Negative
pos_2 = np.array([-4.4, -1.44, 2.0]) * length_scale # Negative

pos_3 = np.array([1.0, -2.34, -1.0]) * length_scale # Positive
pos_4 = np.array([1.0, -2.34, 1999.0]) * length_scale # Positive

pos_5 = np.array([-4.4, -1.44, -2.0]) * length_scale # Positive
pos_6 = np.array([-4.4, -1.44, 1998.0]) * length_scale # Positive

# Acceleration on particle 1 from others
# Particle 2 (Negative)
diff = pos_1 - pos_2
r = np.sqrt( np.sum(diff**2) ) #+ length_scale**3 # + length_scale**3 to avoid division with zero
force_2 = (-1.0)*(-1.0)*pre_fac_c * diff / r**3


# Particle 3 (Positive)
diff = pos_1 - pos_3
r = np.sqrt( np.sum(diff**2) ) #+ length_scale**3 # + length_scale**3 to avoid division with zero
force_3 = (-1.0)*(+1.0)*pre_fac_c * diff / r**3
force_3 = 0.0

# Particle 4 (Positive)
diff = pos_1 - pos_4
r = np.sqrt( np.sum(diff**2) ) #+ length_scale**3 # + length_scale**3 to avoid division with zero
force_4 = (-1.0)*(+1.0)*pre_fac_c * diff / r**3
force_4 = 0.0


# Particle 5 (Positive)
diff = pos_1 - pos_5
r = np.sqrt( np.sum(diff**2) ) #+ length_scale**3 # + length_scale**3 to avoid division with zero
force_5 = (-1.0)*(+1.0)*pre_fac_c * diff / r**3

# Particle 6 (Positive)
diff = pos_1 - pos_6
r = np.sqrt( np.sum(diff**2) ) #+ length_scale**3 # + length_scale**3 to avoid division with zero
force_6 = (-1.0)*(+1.0)*pre_fac_c * diff / r**3

# Total force

force_TOT = force_2 + force_3 + force_4 + force_5 + force_6
force_TOT[2] = force_TOT[2] + (-q_0)*(-V_0/d) # Vacuum field
accel_TOT = force_TOT/m_0
#print((-q_0)*(-V_0/d)/m_0)
print(accel_TOT)