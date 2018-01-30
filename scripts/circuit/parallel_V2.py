#!/usr/bin/env python3
#
#-------------------------------------------------------------------------------
# Test of node analysis for a resistor and capacitor in parallel with the diode
# Kristinn Torfason
# 25.01.2018
# Nodal Analysis version with finite difference
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, m_e, epsilon_0, pi

#-------------------------------------------------------------------------------
# Child Langmuire limited System
# Set I_D = Current_Diode_Child

# Scales
time_scale = 1.0E-12    # Time scale [s]
length_scale = 1.0E-9  # Length scale [m]

# System variables
d         = 500.0*length_scale # Gap spacing [nm]
L         = 250.0*length_scale  # Side length of emitter [nm]
time_step = 0.25E-3*time_scale     # Size of the time step [s]
steps     = 1000000             # Number of time steps [time_step]

# Curcuit elements
V_0 = 2.0       # Voltage source [V]
R_D = 1.0E6     # Diode resistor [Ohm]
R_C = 1.0E6     # Capacitor resistor [Ohm]
R   = 1.0E6     # Voltage source resistor [Ohm]
#C   = 10.0E-18  # Capacitor [Farad]

# Parallel plate capacitor C = ε_0 A / d
C  = epsilon_0*L**2/d # Capacitor [Farad]

#-------------------------------------------------------------------------------
# Test System
# Set I_D = Current_Diode

# # Scales
# time_scale = 1.0E-6    # Time scale [s]
# length_scale = 1.0E-9  # Length scale [m]
#
# # System variables
# d         = 200.0*length_scale  # Gap spacing [nm]
# L         = 100.0*length_scale  # Side length of emitter [nm]
# time_step = 0.25E0*time_scale     # Size of the time step [s]
# steps     = 1000000             # Number of time steps [time_step]
#
# # Curcuit elements
# V_0 = 2.0     # Volts (Voltage source)
# R_D = 1.0E3   # Ohm   (Diode resistor)
# R_C = 1.0E3   # Ohm   (Capacitor resistor)
# R   = 1.0E3   # Ohm   (Voltage source resistor)
# C   = 10.0E-6 # Farad (Capacitor)
#-------------------------------------------------------------------------------

# To store the previous value of the current density from Child's law
J_cl_prev = 0.0


""" Returns the current from the diode
Assumes that the diode is limited by Child's law and takes the tranist time
to reach full value.

Keyword arguments:
V_diode -- The voltage of the diode
step    -- The current time step
"""
def Current_Diode_Child(V_diode: float, step: int) -> float:
    global J_cl_prev # Store the previous value

    if V_diode > 0.0:

        # Calculate the propagation time and use that as the time it takes to turn on the current, i.e. reach max value
        t_on = d*np.sqrt(2.0*m_e/(e*V_0)) # Turn on time to reach limit [s]
        t_on_step = np.rint(t_on / time_step) # Turn on time to reach limit [time_step]

        J_cl = 4.0*epsilon_0/9.0 * np.sqrt(2.0*e/m_e)*V_diode**(3.0/2.0)/d**2
        if step <= t_on_step:
            J_cl = J_cl/t_on_step * step
    else:
        # If the voltage is less then zero then no current
        J_cl = 0.0

    # The area of the emitter
    Area = L**2

    # Mix the new and old value. To large jumps cause instabillity.
    x = 0.7 # Use 70% of the old value
    J_cl = J_cl_prev*x + (1.0-x)*J_cl

    # Calculate the current
    I_cl = J_cl*Area

    # Keep the previous value
    J_cl_prev = J_cl

    if step > 500000:
        #I_cl = 0.0
        I_cl = I_cl*( 1.0 + 0.1*np.sin(2.0*pi/(10*time_scale)*step*time_step) )

    return I_cl


""" Returns the Current from the diode. This functions assumes a constant
current source.

Keyword arguments:
V_diode -- Voltage of the diode
step    -- The current time step
"""
def Current_Diode(V_diode: float, step: int) -> float:
    return 10.0E-3

# ------------------------------------------------------------------------------
# Start of the simulation
# Initialize variables

print('Starting Simulation')
print('V_0 = {:.2G} V'.format(V_0))
print('R = {:.2E} Ω'.format(R))
print('R_C = {:.2E} Ω'.format(R_C))
print('R_D = {:.2E} Ω'.format(R_D))
print('C = {:.2E} F'.format(C))
print('time_step = {:.4G} ps'.format(time_step/time_scale))
print('d = {:.4G} nm'.format(d/length_scale))
print('L = {:.4G} nm'.format(L/length_scale))
print('')

V_cur  = np.zeros(4)
V_prev = np.zeros(4)

V_C  = np.zeros(steps)
V_D  = np.zeros(steps)
V_0C = np.zeros(steps)

I_T  = np.zeros(steps)
I_DC = np.zeros(steps)
I_C  = np.zeros(steps)
I_D  = 0.0

time = np.zeros(steps)

A = np.zeros((4, 4))
b = np.zeros(4)

#-------------------------------------------------------------------------------
# Main loop

print('Starting main loop with {:d} steps.'.format(steps))
for step in range(steps):
    A[0, 0] = 1.0/R_D + 1.0/R_C
    A[0, 1] = -1.0/R_D
    A[0, 2] = -1.0/R_C
    A[0, 3] = 1.0/R
    b[0]    = 0.0

    A[1, 0] = -1.0/R_D
    A[1, 1] = 1.0/R_D
    b[1]    = -1.0*I_D

    A[2, 0] = -1.0/R_C
    A[2, 2] = 1.0/R_C + C/time_step
    b[2]    = C/time_step*V_prev[2]

    A[3, 0] = 1.0
    A[3, 3] = -1.0
    b[3]    = V_0

    # Store previous values of the voltages
    V_prev[:] = V_cur

    # Solve the system of equations
    V_cur = np.linalg.solve(A, b)

    # Calculate voltages
    V_D[step] = V_cur[1]             # Voltage over the diode
    V_C[step] = V_cur[2]             # Voltage of the capacitor
    V_0C[step] = V_cur[0] - V_cur[3] # Source voltage (Should be equal to V_0)

    # Calculate currents from voltages over resistors
    I_T[step]  = ( -1.0*V_cur[3] / R ) / 1.0E-6            # Total current
    I_C[step]  = ( (V_cur[0] - V_cur[2]) / R_C ) / 1.0E-6  # Capacitor current
    I_DC[step] = ( (V_cur[0] - V_cur[1]) / R_D ) / 1.0E-6  # Diode current

    I_D = Current_Diode_Child(V_D[step], step)
    #I_D = Current_Diode(V_D[step], step)

    # Store the current time
    time[step] = step*time_step / time_scale

print('')
print('Simulation done')

# ------------------------------------------------------------------------------
# Plot results
print('Ploting results')
#print('V_cur = ', V_cur)

plt.figure()
plt.plot(time, V_C, time, V_D)
plt.title('Voltage')
plt.legend(['Capacitor', 'Diode'])
plt.xlabel('Time [ps]')
plt.ylabel('V [V]')
#plt.show()

plt.figure()
plt.plot(time, I_C, '-', time, I_DC, '-', time, I_T, ':')
plt.legend(['Capacitor', 'Diode', 'Total'])
plt.title('Current')
plt.xlabel('Time [ps]')
plt.ylabel('I [uA]')
#plt.show()

plt.show()
