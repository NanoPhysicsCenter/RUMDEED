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

# Scale factors
time_scale = 1.0E-12
length_scale = 1.0E-9

# System parameters
d = 100.0*length_scale # Gap spacing (nm)
V = 2.0 # Voltage over gap (V)
delta_t = 0.25E-3*time_scale # Time step (ps)

R_1 = np.array([0.0, 0.0, 1.0])*length_scale

time_exp = np.sqrt(2.0*d*(d - R_1[2])*m_0/(q_0*V))
time_exp_int = time_exp / delta_t

print(time_exp)
print(time_exp_int)
print(np.ceil(time_exp_int))
