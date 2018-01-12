#!/usr/bin/env python3
# Unit test for the acceleration

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
pre_fac_F = q_0/(4*pi*epsilon_0)
pre_fac_E = q_0/m_0


# Scale factors
time_scale = 1.0E-12
length_scale = 1.0E-9

# System parameters
d = 100.0*length_scale # Gap spacing (nm)
V = 2.0 # Voltage over gap (V)
delta_t = 0.25E-3*time_scale # Time step (ps)

#a_scale_c = time_scale**2/length_scale**3
#a_scale_E = time_scale**2/length_scale**3

# Electric field in system
E = np.array([0.0, 0.0, -V/d])

# Particle locations
# Electrons
R_1 = np.array([3.0, -10.0, 101.0])*length_scale
R_2 = np.array([-9.0, 26.0, 80.0])*length_scale

# Holes
R_3 = np.array([6.0, -24.0, 118.0])*length_scale


# Acceleration particles 1 ----------------------------------------------------
a_12 = +1.0*pre_fac_c * (R_1 - R_2) / np.linalg.norm(R_1 - R_2)**3
a_13 = -1.0*pre_fac_c * (R_1 - R_3) / np.linalg.norm(R_1 - R_3)**3
a_1E = -1.0*pre_fac_E * E

a_1 = a_12 + a_13 + a_1E

print('Particle 1')
print(a_1)
print('')


# Acceleration particles 2 -----------------------------------------------------
a_21 = +1.0*pre_fac_c * (R_2 - R_1) / np.linalg.norm(R_2 - R_1)**3
a_23 = -1.0*pre_fac_c * (R_2 - R_3) / np.linalg.norm(R_2 - R_3)**3
a_2E = -1.0*pre_fac_E * E

a_2 = a_21 + a_23 + a_2E

print('Particle 2')
print(a_2)
print('')


# Acceleration particles 2 -----------------------------------------------------
a_31 = -1.0*pre_fac_c * (R_3 - R_1) / np.linalg.norm(R_3 - R_1)**3
a_32 = -1.0*pre_fac_c * (R_3 - R_2) / np.linalg.norm(R_3 - R_2)**3
a_3E = +1.0*pre_fac_E * E

a_3 = a_31 + a_32 + a_3E

print('Particle 3')
print(a_3)
print('')


# Electric field at position ---------------------------------------------------
R_4 = np.array([-4.55, -2.34, 96.44])*length_scale

E_1 = -1.0*pre_fac_F * (R_4 - R_1) / np.linalg.norm(R_4 - R_1)**3
E_2 = -1.0*pre_fac_F * (R_4 - R_2) / np.linalg.norm(R_4 - R_2)**3
E_3 = +1.0*pre_fac_F * (R_4 - R_3) / np.linalg.norm(R_4 - R_3)**3

print('pre_fac_F = ', pre_fac_F)
print('diff = ', (R_4 - R_1))
print('r = ', np.linalg.norm(R_4 - R_1))

E_tot = E_1 + E_2 + E_3 + E
print('Electric field')
print(E_tot)
print('')
