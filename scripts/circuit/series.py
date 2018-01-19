#!/usr/bin/env python3
#
#-------------------------------------------------------------------------------
# Test of equations for a resistor and capacitor in series with diode
# Kristinn Torfason
# 17.01.2018
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
R = 0.5E6 # Resistance of the resistor in series [Ohm]
C = 1.0E-18 # Capacitance of the capacitor in series [F]
time_step = 0.25E-3*time_scale # Size of the time step [S]
steps = 100000 # Number of time steps

# Variables for functions
ramo_integral = 0.0
I_Cprev = 0.0

# Gives the voltage over the resistor
def Voltage_Resistor(I_R):
    return R*I_R

# Gives the voltage of the capacitor
def Voltage_Capacitor(I_C):
    global ramo_integral, I_Cprev

    ramo_integral = ramo_integral + time_step * (I_Cprev + I_C) * 0.5
    I_Cprev = I_C
    return 1.0/C * ramo_integral

# Returns the current from the diode
# Assumes that the diode is a resistor with a turn on time t_on to reach
# full value.
def Current_Diode_Resistor(V_d, step):
    R_D = 0.5E6
    t_on = 100 # Turn on time to reach limit [steps]

    if (step > t_1):
        return V_d/R_D
    else:
        return V_d/(R_D*t_on)*step

# Returns the current from the diode
# Assumes that the diode is limited by Child's law and takes the tranist time
# to reach full value.
def Current_Diode_Child(V_d, step):
    t_on = d*np.sqrt(2.0*m_e/(e*V_d)) # Turn on time to reach limit [s]
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
    V_R = Voltage_Resistor(I)
    V_C = Voltage_Capacitor(I)

    V = V_S - V_R - V_C # Voltage over diode
    I = Current_Diode_Child(V, step) # Current from the diode

    #print(V_d)
    #print(I)
    #print('')
    voltage[step] = V
    current[step] = I / 1.0E-6
    time[step] = step*time_step / time_scale
    #if (step == 0):
    #print('V = ', V)
    #print('V_S = ', V_S)
    #print('V_R = ', V_R)
    #print('V_C = ', V_C)
    #print('I = ', I)
    #print('')

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
