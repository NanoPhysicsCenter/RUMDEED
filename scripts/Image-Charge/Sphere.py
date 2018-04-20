# Calculate the Electric Field on a sphere using the image charge method
# The field is calculated on the sphere as a function of distance of the
# particle.

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import pi, epsilon_0, e

#-------------------------------------------------------------------------------
# Constants the define the tip geometry

#Base radius
R = 250.0E-9 # [nm]

# Height
h = 500E-9 # [nm]

# Gap spacing, measured from the top of the tip.
d = 1000.0E-9 # [nm]
d_plane = d + h # Distance from bottom of tip [nm]

# Voltage
V_0 = 500.0 # [V]

# Tip parameters
a = np.sqrt( d**2*R**2/(h**2+2*d*h) + d**2 )
xi_max = h/d + 1.0
eta_1 = -d/a


max_xi = h/d + 1.0
max_x = a*np.sqrt(max_xi**2-1.0)*np.sqrt(1.0-eta_1**2)
x_tip = np.linspace(-max_x, max_x, 1000)

y_tip = 0.0

xi = np.sqrt(x_tip**2/(a**2*(1.0-eta_1**2)) + 1.0)
z_tip = a * xi * eta_1

#-------------------------------------------------------------------------------
# Radius of sphere
R_a = np.abs(a**2/d-d)

# Center of the sphere. We place it at the top of the tip.
r_c = np.array([0.0, 0.0, h-R_a])

# Size of charge
q = +e # in C

# Position of the charged particle outside the sphere (tip).
r_p = np.array([0.0, 0.0, h+10.0E-9])

# Vector from center of sphere r_c = (x_c, y_c, z_c) to particle r_p = (x_p, y_p, z_p)
R_cp = r_p - r_c

# Distance of the particle from the center of the sphere.
d_p = np.sqrt(np.inner(R_cp, R_cp))

# Scale it to unit size.
R_cp = R_cp / d

# Distance of the image charge partner from the center of the sphere.
d_ic = R_a**2/d

# The position of the image charge partner.
r_ic = r_c + R_cp*d_ic

# Now calculate the size of the charge of the image charge partner.
q_ic = -1.0*R_a/d*q

#-------------------------------------------------------------------------------
E_tot = np.array([])

for i, (x, z) in enumerate(zip(x_tip, z_tip)):
    r = np.array([x, 0.0, z])
    r_kp  = r - r_p
    r_kic = r - r_ic

    d_kp = np.sqrt( np.inner(r_kp, r_kp) )
    d_kic = np.sqrt( np.inner(r_kic, r_kic) )

    E_p = q/(4.0*pi*epsilon_0)*r_kp / d_kp**3
    E_ic = q_ic/(4.0*pi*epsilon_0)*r_kic / d_kic**3

    E = E_p + E_ic
    E_tot = np.append(E_tot, np.sqrt( np.inner(E, E) ))

plt.plot(x_tip/1E-9, E_tot)
plt.show()
