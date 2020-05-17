import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import pi, hbar, e, m_e, epsilon_0
from scipy.integrate import nquad, quad

# Fluid model code
# Kristinn Torfason (29.11.2017) (Converted from Current_density_scaled.m)
# Revised (17.05.2020). Here multiple emitters are added.
#
# This code calculates the fluid model for the field emission
# The continuity equation, \rho(z) p/m_e = J
# Conservation of energy, p^2 / (2m_e) = eV_0 z/d
# Combined they give, \rho(z) = J/\sqrt(m_e d / (2eV_0 z) )
#

#-------------------------------------------------------------------------------
# Define emitters

def Checkerboard(N, L):
    # Side lengths of emitters.
    emitter_Lx = np.ones(N*N) * L/6.0
    emitter_Ly = np.ones(N*N) * L/6.0

    # Work function of each emitter
    emitter_w = np.zeros(N*N)

    # Center coordinates of all emitters
    emitter_x = np.zeros(N*N)
    emitter_y = np.zeros(N*N)
    for j in range(N):

        for i in range(N):
            y = L/(2*N) * (j+1)
            x = L/(2*N) * (i+1)

            emitter_x[i+j] = x
            emitter_y[i+j] = y

            if (((i+j) % 2) == 0):
                emitter_w[i+j] = 2.5
            else:
                emitter_w[i+j] = 2.0

    return emitter_x, emitter_y, emitter_Lx, emitter_Ly, emitter_w

# Checkerboard with N = 6 and L = 1000.0 nm
N = 6
L = 1000.0
emitter_x, emitter_y, emitter_Lx, emitter_Ly, emitter_w = Checkerboard(N, L)


plt.scatter(emitter_x, emitter_y)
for i in range(N*N):
    plt.annotate('{:.2f}'.format(emitter_w[i]), (emitter_x[i], emitter_y[i]))

plt.show()
exit()

#-------------------------------------------------------------------------------
# Fluid model
# This code calculates the fluid model for the field emission.
# The continuity equation, $\rho(z) \frac{p}{m_e} = J$
# Conservation of energy, $\frac{p^2}{2m_e} = eV_0 \frac{z}{d}$
# Combined they give, $\rho(z) = \frac{J}{\sqrt{z}} \sqrt{\frac{md}{2eV}}$

# Constants
a_FN = e**2/(16*pi**2*hbar) # First Fowler Nodheim constant
b_FN = -4/(3*hbar) * np.sqrt(2*m_e*e) # Second Fowler Nordheim constant
l_const = e / (4*pi*epsilon_0) # Constant for image charge calculations
x = 0.15 # Mixing weight in iteration
MAX_ITER = 500 # Maximum number of iterations before giving upp

# System parameters
V = 1.0 # Voltage [V]
d = 1.0 # Gap spacing [nm]

E_vac = V / (d*1.0E-9) # [V/m]

# Scale the current density with the Child-Langmuir law in 1D
J_CL1D = 4/9*epsilon_0*np.sqrt(2*e/m_e)*V**(3/2)/(d*1.0E-9)**2

def Set_System(V_: float, d_: float):
    global V, d, E_vac, J_CL1D
    V = V_
    d = d_
    E_vac = V/(d*1.0E-9)
    J_CL1D = 4/9*epsilon_0*np.sqrt(2*e/m_e)*V**(3/2)/(d*1.0E-9)**2

# Functions

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

# Fowler-Nordheim current density
def F_J(F_scaled: float, w_theta: float) -> float:
    F = E_vac * F_scaled
    l = l_const * F / w_theta**2
    elec_supply = a_FN * F**2 / (w_theta * t_y(l)**2)
    esc_prob = np.exp(b_FN * w_theta**(3/2) * v_y(l) / F)
    J = elec_supply * esc_prob / J_CL1D
    return J

# We want to integrate the charge density over x, y and z
# to get the electric field in center of the emitter.
# The function to integrate over x', y' and z' is
# z' / ( \sqrt(z') * ((x_c - x')^2 + (y_c - y')^2 + z'^2)^(3/2) )
# We can do the x and y integration, Matlab gives that the answer is
# 1/\sqrt(z') * atan( (x_c - x')(y_c - y') / (z \sqrt(z^2 + (x_c - x')^2 + (y_c - y')^2)) )
# Here x_c and y_c are the center coordinates of the emitter we are calculating the field over.
# x' and y' are the integration variables, i.e. the emitter we are integrating over.
# all the spatial coordinates are scaled with the gap spacing, d. 
# We then do the z integration numerically, from 0 to 1.
def int_fun_z(z: float, x_c: float, y_c: float, x_1: float, x_2: float, y_1: float, y_2: float) -> float:
    # The indefinite integral over x and y that Matlab/Wolfram give
    def fun_z(z: float, x_c: float, y_c: float, x_p: float, y_p: float) -> float:
        val = 1.0/np.sqrt(z)*np.arctan( (x_c - x_p)*(y_c - y_p)/(z*np.sqrt(z**2 + (x_c - x_p)**2 + (y_c - y_p)**2 )) )
        return val
    
    # Double integral over x and y
    val = fun_z(z, x_c, y_c, x_2, y_2) - fun_z(z, x_c, y_c, x_1, y_2) - fun_z(z, x_c, y_c, x_2, y_1) + fun_z(z, x_c, y_c, x_1, y_1)   
    return val

@lru_cache(maxsize=None)
def do_int(x_c: float, y_c: float, x_cp: float, y_cp: float, L_x: float, L_y: float) -> float:
    # Calculate the coordinates of the emitter
    #           |-----------| (x_2, y_2)
    #           |           |
    #           |(x_c, y_c) L_y
    #           |           |
    # (x_1, y_1)|----L_x----|
    x_1 = x_cp - L_x/2.0
    x_2 = x_cp + L_x/2.0
    
    y_1 = y_cp - L_y/2.0
    y_2 = y_cp + L_y/2.0
    
    # Do the z integration from 0 to 1
    E_val, abserr = quad(int_fun_z, a=0.0, b=1.0, args=(x_c, y_c, x_1, x_2, y_1, y_2), epsabs=0.5E-12)
    
    # Return the results of the integration
    val = E_val
    return val


# Use this function to set system parameters
# Set_System(Voltage, Gap spacing)
Set_System(2250, 1000.0)


#-------------------------------------------------------------------------------
# Iteration: The function that does the iterations.

def Do_Fluid_Model(emitter_x, emitter_y, emitter_Lx, emitter_Ly, emitter_w):
    
    # Get the number of emitters
    nrEmit = len(emitter_x)
    
    # Define array's that are needed for calculations
    emitter_E    = np.zeros(nrEmit)
    emitter_J    = np.zeros(nrEmit)
    emitter_Jnew = np.zeros(nrEmit)

    # Set the initial field and calculate the inital current density from each emitter
    for i in range(nrEmit):
        # Start with the vacuum field
        emitter_E[i]    = 1.0 # Scaled in E_vac
        
        # Start at full blast and work our way down, could also start at 0 and work uppwards but that takes longer
        emitter_J[i]    = F_J(emitter_E[i], emitter_w[i])
        
    #-----------------------------------------------------------
    abserr = np.finfo(np.float64).max # Absolute error
    relerr = 0.0 # Relative error

    for k in range(MAX_ITER):
        for i in range(nrEmit):

            # We are calculating the field and current density for this emitter
            x_c = emitter_x[i]
            y_c = emitter_y[i]

            E_z = 0.0
            for j in range(nrEmit):

                # We are calculating the influence of this emitter on the emitter at x_c and y_c 
                x_cp = emitter_x[j]
                y_cp = emitter_y[j]

                L_x = emitter_Lx[j]
                L_y = emitter_Ly[j]

                w_theta = emitter_w[j]
                J       = emitter_J[j]

                E_int = do_int(x_c, y_c, x_cp, y_cp, L_x, L_y)
                E_z = E_z + J / (9*pi) * E_int

            F = 1.0 - 2.0*E_z
            emitter_E[i] = F

            J_old = emitter_J[i]
            J_new = F_J(F, w_theta)
            J = x*J_new + (1-x)*J_old # Mixing for better convergance
            emitter_Jnew[i] = J

            # Calculate the absolute error and relative error
            abserr = np.minimum(abserr, np.abs(J-J_old))
            relerr = abserr/J

        emitter_J = emitter_Jnew
        if (abserr < 0.5E-12):
            break
        if (k == MAX_ITER-1):
            print('Warning reached MAX_ITER')
    
    # Return the current density and the surface field
    return emitter_J, emitter_E


#-------------------------------------------------------------------------------
# Results

# Scale the variables with d
emitter_x  = emitter_x/d
emitter_y  = emitter_y/d
emitter_Lx = emitter_Lx/d
emitter_Ly = emitter_Ly/d
nrEmit     = len(emitter_x)

J = np.zeros(nrEmit)
E = np.zeros(nrEmit)
    
# Do the fluid model calculations
J, E = Do_Fluid_Model(emitter_x, emitter_y, emitter_Lx, emitter_Ly, emitter_w)