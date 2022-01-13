import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import pi, hbar, e, m_e, epsilon_0
from scipy.integrate import nquad, quad

# Fluid model code
# Kristinn Torfason (29.11.2017) (Converted from Current_density_scaled.m)
# Revised (5.11.2018)
#
# This code calculates the fluid model for the field emission
# The continuity equation, \rho(z) p/m_e = J
# Conservation of energy, p^2 / (2m_e) = eV_0 z/d
# Combined they give, \rho(z) = J/\sqrt(m_e d / (2eV_0 z) )
#

# System Parameters
V_0 = 1250.0 # Voltage
d = 1000.0E-9 # Gap spaceing
N_L = 2
L = np.linspace(100.0, 1000.0, N_L) * 1.0E-9 # Length of emitter

# Mixing weight
x = 0.15

# Constants
a_FN = e**2/(16*pi**2*hbar)
b_FN = -4/(3*hbar) * np.sqrt(2*m_e*e)
l_const = e / (4*pi*epsilon_0)
w_theta = 2.00

E_vac = V_0 / d
F = 1.0 # F is scaled in E_vac

# Maximum number of iterations
N_max = 500

J_keep = np.zeros(N_max)
F_keep = np.zeros(N_max)

J_L = np.zeros(N_L)
F_L = np.zeros(N_L)
LD_C = np.zeros(N_L)

# Scale the current density with Child-Langmuir
J_CL = 4/9*epsilon_0*np.sqrt(2*e/m_e)*V_0**(3/2)/d**2

# Approximations to the eleptical integrals for the image-charge effect
# in the Fowler-Nordheim equation
def t_y(l: float) -> float:
    val = 1 + l*(1/9 - 1/18*np.log(l))
    if l > 1:
        print('Warning: l > 1') # This should not happen
    return val

def v_y(l: float) -> float:
    val = 1 - l + 1/6 * l * np.log(l)
    return val

# We want to integrate the charge density over x, y and z
# to get the electric field in center of the emitter.
# The function to integrate is
# J/(9\pi) * z / (\sqrt(z) * (x^2 + y^2 + z^2)^(3/2)  )
# First we do the x integration by hand and then do the
# y and z integrations numerically.
# When we do the x integration we obtain, ( a = L/(2*d) )
# J/(9\pi) * 2*a \sqrt(z) / ( (y^2 + z^2) \sqrt(y^2 + z^2 + a^2) ).
# This is the function we integrate numerically below. z = [0, 1], y = [0, a]
def int_fun(y: float, z: float, a: float) -> float:
    val = np.sqrt(z) / ((y**2 + z**2)*np.sqrt(a**2 + y**2 + z**2))
    return val

# We can also do the y integration. Matlab gives the answer below.
def int_fun_z(z: float, x_c: float, y_c: float, x_1: float, x_2: float, y_1: float, y_2: float) -> float:
    def fun_z(z: float, x_c: float, y_c: float, x_p: float, y_p: float) -> float:
        val = np.arctan( (x_c - x_p)*(y_c - y_p)/(z*np.sqrt(z**2 + (x_c - x_p)**2 + (y_c - y_p)**2 )) )
        return val

    val = 1.0/np.sqrt(z)*( fun_z(z, x_c, y_c, x_2, y_2) - fun_z(z, x_c, y_c, x_1, y_2) - fun_z(z, x_c, y_c, x_2, y_1) + fun_z(z, x_c, y_c, x_1, y_1) )
    return val

print('Starting calculations')
print('k upto ' + str(N_L))
# Loop over L values
for k in range(N_L):

    # Set initial values
    LD_C[k] = L[k]/(2*d)

    J_old = 0.0
    J_keep[:] = 0.0
    F_keep[:] = 0.0

    # Set a = L / (2*d)
    a = LD_C[k]

    x_c = 0.0/d
    y_c = 0.0/d
    x_1 = -L[k]/(2*d)
    x_2 = L[k]/(2*d)
    y_1 = -L[k]/(2*d)
    y_2 = L[k]/(2*d)

    options={'limit': 100} # Increase the number of subintervals in the integration
    # nquad does an N-dimensional integration
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.nquad.html
    # integrate int_fun, y=0..a, z=0..1, with extra arguments (a,) passed to the function
    # We change the integration on y from [-a, a] to [0, a] and multiply with 2 instead (4a = 2*2a).
    #E_const, abserr = nquad(int_fun, [[0.0, a], [0.0, 1.0]], args=(a,), opts=[options, options])
    E_const, abserr = quad(int_fun_z, a=0.0, b=1.0, args=(x_c, y_c, x_1, x_2, y_1, y_2))

    # Loop for the iteration to converge the current density
    for i in range(N_max):
        F = E_vac * F
        l = l_const * F / w_theta**2

        # Calculate the current density from the Fowler-Nordheim
        elec_supply = a_FN * F**2 / (w_theta * t_y(l)**2)
        esc_prob = np.exp(b_FN * w_theta**(3/2) * v_y(l) / F)

        J_new = elec_supply * esc_prob / J_CL
        J = x*J_new + (1-x)*J_old # Mixing for better convergance

        #E_z = J/(9*pi) * 4*a * E_const
        E_z = J/(9*pi) * E_const

        # Set F to new value. The factor 2 is due to image-charge effects
        # F = F_vac - 2*E_z,
        # but we scale the field with F_vac so
        # F = 1 - 2*E_z
        F = 1 - 2*E_z

        # Store the values that we want to keep
        F_keep[i] = F
        J_keep[i] = J
        J_old = J

        if (i >= 2):
            if np.abs((J_keep[i-1] - J_keep[i-2])) < 1.0E-12:
                print('k = ' + str(k) + ' done ' + str(i))
                break

    # Check the convergance
    if np.abs((J_keep[i-1] - J_keep[i-2])) > 1.0E-12:
        print('Warning error > 1E-12')

    # Set values that we want to keep
    J_L[k] = J
    F_L[k] = F

# Plot the results
print('Ploting results')
plt.plot(L/1.0E-9, J_L*J_CL)
plt.xlabel('L [nm]')
plt.ylabel('J [A/m^2]')
plt.show()

print(L/1.0E-9)
print(J_L*J_CL)
