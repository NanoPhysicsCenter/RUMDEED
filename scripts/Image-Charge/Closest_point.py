# Find the closest point on the prolate spheroidal Tip

import numpy as np
import matplotlib.pyplot as plt

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

x_p = 4.621511002908029e-08
y_p = 0.0
z_p = -0.518536659185618e-06

def Find_Point(x_p, y_p, z_p):
    x = x_p # Initial guess for x

    for i in range(1000):
        fac = 1.0 + eta_1**2/(1.0-eta_1**2) - eta_1/np.sqrt(1.0-eta_1**2)*z_p/np.sqrt(x**2*(1.0+y_p**2/x_p**2) + a*(1.0-eta_1**2))
        x = x_p / fac

    y = y_p/x_p*x
    z = eta_1/np.sqrt(1.0-eta_1**2)*np.sqrt(x**2 + y**2 + a*(1.0-eta_1**2))

    return x, y, z

x_a, y_a, z_a = Find_Point(x_p, y_p, z_p)
print(x_a)
print(y_a)
print(z_a)

x_tip = np.linspace(-R, R, 1000)
y_tip = 0.0
xi = np.sqrt(x_tip**2/(a**2*(1-eta_1**2)) + 1)
z_tip = a * xi * eta_1

x_plane = np.linspace(-R, R, 10)
z_plane = x_plane*0.0

plt.plot(x_tip/1E-9, z_tip/1E-9)
plt.plot(x_plane/1E-9, z_plane/1E-9)
#plt.scatter(x_a/1E-9, z_a/1E-9)
plt.scatter(x_p/1E-9, z_p/1E-9)
plt.xlabel('x [nm]')
plt.ylabel('z [nm]')
plt.show()
