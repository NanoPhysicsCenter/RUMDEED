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
h = 500.0E-9 # [nm]

# Gap spacing, measured from the top of the tip.
d = 1000.0E-9 # [nm]
d_plane = d + h # Distance from bottom of tip [nm]

# Voltage
V_0 = 500.0 # [V]

# Tip parameters
a = np.sqrt( d**2*R**2/(h**2+2*d*h) + d**2 )
xi_max = h/d + 1.0
eta_1 = -d/a

#-------------------------------------------------------------------------------
# Radius of sphere
theta = np.arccos(d/a)
r_tip = a*np.sin(theta)*np.tan(theta)
print(r_tip/1E-9)
#a = 4.999999999999993E-08 # in m

# Size of charge
q = +e # in C

# distance from sphere
#d = np.linspace(r_tip+0.0001E-9, r_tip+10.0E-9, 10000)
d = 50.18E-9

# (a-d)/abs(a-d)**3
term_1 = (r_tip - d)/np.abs(r_tip-d)**3

# (1/(a*d) * (1-a/d)/abs(1-a/d)**3)
d_p = r_tip**2/d
#term_2 = -1.0/(r_tip*d) * (1.0-r_tip/d)/np.abs(1.0-r_tip/d)**3
term_2 = -r_tip/d*(r_tip - d_p)/np.abs(d_p - r_tip)**3

E = q/(4.0*pi*epsilon_0) * ( term_1 + term_2 )

#plt.plot(d/1E-9, np.abs(E))
#plt.show()
print(E)
print(term_1+term_2)
print(term_1)
print(term_2)
