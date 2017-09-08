#!/usr/bin/python3
# Unti test for the acceleration

import numpy as np

# Physical constants
pi = 3.14159265358979324 # Pi
h = 6.626070040E-34 # Planck's constant (Js)
k_b = 1.38064852E-23 # Boltzmann constant (J/K)
c = 299792458.0 # Speed of light (m/s)
mu_0 = 4.0*pi * 1E-7 # Vacuum permeability (H/m)
epsilon_0 = 1.0/(mu_0 * c**2) # Vacuum permittivity (F/m)
m_0 = 9.10938356E-31 # Free electron mass (kg)
q_0 = 1.6021766208E-19 # Standard charge (C)

# Defined constants
pre_fac_c = q_0**2/(4*pi*m_0*epsilon_0)
pre_fac_E = q_0/m_0

# System parameters
d = 100.0 # Gap spacing (nm)
V = 2.0 # Voltage over gap (V)
delta_t = 0.25 # Time step (ps)

# Scale factors
time_scale = 1.0E-12
length_scale = 1.0E-9

a_scale_c = time_scale**2/length_scale**3
a_scale_E = time_scale**2/length_scale**3

# Electric field in system
E = np.array([0.0, 0.0, -V/d])

# Particle locations
# Electrons
R_1 = np.array([3.0, -10.0, 101.0])
R_2 = np.array([-9.0, 26.0, 80.0])

# Holes
R_3 = np.array([6.0, -24.0, 118.0])


# Acceleration particles 1 ----------------------------------------------------
a_12 = +1.0*pre_fac_c * a_scale_c * (R_1 - R_2) / np.linalg.norm(R_1 - R_2)**3
a_13 = -1.0*pre_fac_c * a_scale_c * (R_1 - R_3) / np.linalg.norm(R_1 - R_3)**3
a_1E = -1.0*pre_fac_E * E * a_scale_E

a_1 = a_12 + a_13 + a_1E

print('Particle 1')
print(a_1)
print('')


# Acceleration particles 2 -----------------------------------------------------
a_21 = +1.0*pre_fac_c * a_scale_c * (R_2 - R_1) / np.linalg.norm(R_2 - R_1)**3
a_23 = -1.0*pre_fac_c * a_scale_c * (R_2 - R_3) / np.linalg.norm(R_2 - R_3)**3
a_2E = -1.0*pre_fac_E * E * a_scale_E

a_2 = a_21 + a_23 + a_2E

print('Particle 2')
print(a_2)
print('')


# Acceleration particles 2 -----------------------------------------------------
a_31 = -1.0*pre_fac_c * a_scale_c * (R_3 - R_1) / np.linalg.norm(R_3 - R_1)**3
a_32 = -1.0*pre_fac_c * a_scale_c * (R_3 - R_2) / np.linalg.norm(R_3 - R_2)**3
a_3E = +1.0*pre_fac_E * E * a_scale_E

a_3 = a_31 + a_32 + a_3E

print('Particle 3')
print(a_3)
print('')
