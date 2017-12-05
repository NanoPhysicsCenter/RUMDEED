import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import *
from scipy.integrate import nquad

# Fluid model code
# Kristinn Torfason (29.11.2017) (Converted from Current_density_scaled.m)
#
# This code calculates the fluid model for the field emission
# The continuity equation, \rho(z) p/m_e = J
# Conservation of energy, p^2 / (2m_e) = eV_0 z/d
# Combined they give, \rho(z) = J/\sqrt(m_e d / (2eV_0 z) )
#

# System Parameters
V_0 = 20E3 # Voltage
d = 2500E-9 # Gap spaceing
N_L = 10
L = np.linspace(10, 100, N_L) * 1E-9 # Length of emitter

# Mixing weight
x = 0.15

# Constants
a_FN = e**2/(16*pi**2*hbar)
b_FN = -4/(3*hbar) * np.sqrt(2*m_e*e)
l_const = e / (4*pi*epsilon_0)
w_theta = 4.7

E_vac = V_0 / d
F = 1.0 # F is scaled in E_vac

# Number of iterations
N = 50

J_keep = np.zeros(N)
F_keep = np.zeros(N)

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

# The function to integrate
# J/(9\pi) * z / (\sqrt(z) * (x^2 + y^2 + z^2)^(3/2)  )
# Fyrst we integrate over x to obtain, ( a = L/(2*d) )
# J/(9\pi) * 2*a \sqrt(z) / ( (x^2 + z^2) \sqrt(y^2 + z^2 + a^2) )
# This is the function we integrate numerically below. z = [0, 1], y = [0, a]
# We change the integration on y from [-a, a] to [0, a] and multiply with 2 instead.
def int_fun(y: float, z: float, a: float) -> float:
    val = np.sqrt(z) / ((y**2 + z**2)*np.sqrt(a**2 + y**2 + z**2))
    return val


# Loop over L values
for k in range(N_L):

    # Set initial values
    LD_C[k] = L[k]/(2*d)

    J_old = 0.0
    J_keep[:] = 0.0
    F_keep[:] = 0.0

    # Loop for the iteration to converge the current density
    for i in range(N):
        F = E_vac * F
        l = l_const * F / w_theta**2

        # Calculate the current density from the Fowler-Nordheim
        elec_supply = a_FN * F**2 / (w_theta * t_y(l)**2)
        esc_prob = np.exp(b_FN * w_theta**(3/2) * v_y(l) / F)

        J_new = elec_supply * esc_prob / J_CL
        J = x*J_new + (1-x)*J_old # Mixing for better convergance

        # Set a = L / (2*d)
        a = LD_C[k]

        options={'limit': 100} # Increase the number of subintervals in the integration
        # nquad does an N-dimensional integration
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.nquad.html
        E_z, abserr = nquad(int_fun, [[0.0, a], [0.0, 1.0]], args=(a,), opts=[options, options])
        E_z = J/(9*pi) * 4*a * E_z

        # Set F to out new value. The factor 2 is for image-charge effects
        F = 1 - 2*E_z

        # Set values that we want to keep
        F_keep[i] = F
        J_keep[i] = J
        J_old = J

    # Check the convergance
    if np.abs((J_keep[N-1] - J_keep[N-2])) > 1E-3:
        print('Warning error > 1E-3');

    # Set values that we want to keep
    J_L[k] = J
    F_L[k] = F

# Plot the results
plt.plot(L/1E-9, J_L*J_CL)
plt.xlabel('L [nm]')
plt.ylabel('J [A/m^2]')
plt.show()
