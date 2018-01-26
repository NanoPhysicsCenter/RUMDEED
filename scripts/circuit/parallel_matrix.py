#!/usr/bin/env python3
#
#-------------------------------------------------------------------------------
# Test of equations for a resistor and capacitor in parallel with diode
# Kristinn Torfason
# 25.01.2018
# See section Circuit elemens in notes
# Matrix Version
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
time_step = 0.25E-3*time_scale # Size of the time step [S]
steps = 1000000 # Number of time steps

V_0 = 2.0     # Volts (Voltage source)
R_D = 1.0E6   # Ohm   (Diode resistor)
R_C = 1.0E6   # Ohm   (Capacitor resistor)
R   = 1.0E6   # Ohm   (Voltage source resistor)
C   = 10.0E-18 # Farad (Capacitor)

#V_0 = 2.0     # Volts (Voltage source)
#R_D = 1.0E3   # Ohm   (Diode resistor)
#R_C = 1.0E3   # Ohm   (Capacitor resistor)
#R   = 1.0E3   # Ohm   (Voltage source resistor)
#C   = 10.0E-6 # Farad (Capacitor)

J_cl_prev = 0.0 # To store the previous value of the current density from Child's law

# Returns the current from the diode
# Assumes that the diode is limited by Child's law and takes the tranist time
# to reach full value.
def Current_Diode_Child(V_diode: float, step: int) -> float:
    global J_cl_prev # Store the previous value

    if (V_diode > 0.0):

        # Calculate the propagation time and use that as the time it takes to turn on the current, i.e. reach max value
        t_on = d*np.sqrt(2.0*m_e/(e*V_0)) # Turn on time to reach limit [s]
        t_on_step = np.rint(t_on / time_step) # Turn on time to reach limit [steps]

        J_cl = 4.0*epsilon_0/9.0 * np.sqrt(2.0*e/m_e)*V_diode**(3.0/2.0)/d**2
        if (step <= t_on_step):
            J_cl = J_cl/t_on_step * step
    else:
        # If the voltage
        J_cl = 0.0

    # The area of the emitter
    A = L**2

    # Mix the new and old value. To large jumps cause instabillity.
    x = 0.7 # Use 70% of the old value
    J_cl = J_cl_prev*x + (1.0-x)*J_cl

    # Calculate the current
    I_cl = J_cl*A

    # Keep the previous value
    J_cl_prev = J_cl

    if (step > 500000):
        I_cl = 0.0

    return I_cl

def Current_Diode(V_diode: float, step: int) -> float:
    return 10.0E-3

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

    V_prev[:] = V_cur
    V_cur = np.linalg.solve(A, b)

    V_D[step] = V_cur[1]
    V_C[step] = V_cur[2]

    I_T[step]  = ( -1.0*V_cur[3] / R ) / 1.0E-6
    I_C[step]  = ( (V_cur[0] - V_cur[2]) / R_C ) / 1.0E-6
    I_DC[step] = ( (V_cur[0] - V_cur[1]) / R_D ) / 1.0E-6
    V_0C[step] = V_cur[0] - V_cur[3]

    I_D = Current_Diode_Child(V_D[step], step)
    #I_D = Current_Diode(V_D[step], step)

    time[step] = step*time_step / time_scale

print('')
print('V_cur = ', V_cur)

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
