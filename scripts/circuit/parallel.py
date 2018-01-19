#!/usr/bin/env python3
#
#-------------------------------------------------------------------------------
# Test of equations for a resistor and capacitor in parallel with diode
# Kristinn Torfason
# 18.01.2018
# See section Circuit elemens in notes
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, m_e, epsilon_0

# Scales
time_scale = 1.0E-12 # Time scale [S]
length_scale = 1.0E-9 # Length scale [m]

# System variables
V_S = 2.0 # Voltage of source [V]
d = 200.0*length_scale # Gap spacing
L = 100.0*length_scale # Side length of emitter
R = 0.5E6 # Resistance of the resistor in the circuit [Ohm]
R_C = 0.5E6 # Resistance of the resistor in series with the capcitor
C = 1.0E-18 # Capacitance of the capacitor in series [F]
time_step = 0.25E-3*time_scale # Size of the time step [S]
steps = 1000000 # Number of time steps

# Variables for functions
ramo_integral = 0.0
I_Cprev = 0.0

def Voltage_Parallel_Capacitor_Resistor(I):
    global ramo_integral
    global I_Cprev
    C   = 10.0E-18 # Farad
    R_C = 0.5E6 # Ohm
    R   = 0.5E6 # Ohm
    ib   = 1.0/(C*(R+R_C))

    # Caclulate the voltage over the diode
    ramo_integral = ramo_integral + time_step * ( I_Cprev*np.exp(-1.0*time_step*ib) + I ) * 0.5
    I_Cprev = I

    return (R**2)/(C*(R+R_C)**2) * ramo_integral - R*R_C/(R+R_C)*I

# Returns the current from the diode
# Assumes that the diode is limited by Child's law and takes the tranist time
# to reach full value.
def Current_Diode_Child(V_d, step):
    if (V_d <= 0.0):
        t_on = d*np.sqrt(2.0*m_e/(e*V_S)) # Turn on time to reach limit [s]
    else:
        t_on = d*np.sqrt(2.0*m_e/(e*V_S)) # Turn on time to reach limit [s]
    t_on_step = np.rint(t_on / time_step) # Turn on time to reach limit [steps]

    if (V_d > 0.0):
        J_cl = 4.0*epsilon_0/9.0 * np.sqrt(2.0*e/m_e)*V_d**(3.0/2.0)/d**2
        if (step <= t_on_step):
            J_cl = J_cl/t_on_step * step
    else:
        J_cl = 0.0

    A = L**2
    I_cl = J_cl*A

    return I_cl

I = 0.0
V = V_S
current = np.zeros(steps)
voltage = np.zeros(steps)
time    = np.zeros(steps)

for step in range(steps):
    I = Current_Diode_Child(V, step) # Current from the diode
    V_C = Voltage_Parallel_Capacitor_Resistor(I)

    V = V_S - V_C # Voltage over diode

    voltage[step] = V
    current[step] = I / 1.0E-6
    time[step] = step*time_step / time_scale

plt.figure()
plt.plot(time, current)
plt.title('Current')
plt.xlabel('t [ps]')
plt.ylabel('I [uA]')

plt.figure()
plt.plot(time, voltage)
plt.title('Voltage')
plt.xlabel('t [ps]')
plt.ylabel('V [V]')

plt.show()
