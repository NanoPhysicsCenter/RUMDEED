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
R = 0.5E6 # Resistance of the resistor in the circuit [Ohm]
R_C = 0.5E6 # Resistance of the resistor in series with the capcitor
C = 1.0E-18 # Capacitance of the capacitor in series [F]
time_step = 0.25E-3*time_scale # Size of the time step [S]
steps = 1000000 # Number of time steps

V_0 = 2.0     # Volts (Voltage source)
R_D = 1.0E6   # Ohm   (Diode resistor)
R_C = 1.0E6   # Ohm   (Capacitor resistor)
R   = 1.0E6   # Ohm   (Voltage source resistor)
C   = 0.0E-18 # Farad (Capacitor)

V_cur = np.zeros(4)
V_prev = np.zeros(4)

V_D = V_cur[1]
V_C = V_cur[2]
I   = V_cur[3] / R

I_D = 10.0E-3

A = np.zeros((4, 4))
b = np.zeros(4)

A[0, 0] = 1.0/R_D + 1.0/R_C
A[0, 1] = -1.0/R_D
A[0, 2] = -1.0/R_C
A[0, 3] = 1.0/R
b[0]    = 0.0

A[1, 0] = -1.0/R_D
A[1, 1] = 1.0/R_D
b[1]    = -1.0*I_D

A[2, 0] = -1.0/R_C
A[2, 2] = 1.0/R_C - C/time_step
b[2]    = C/time_step*V_prev[2]

A[3, 0] = 1.0/R
A[3, 3] = -1.0/R
b[3]    = V_0/R

V_prev = V_cur
print('A = ')
print(A)
print('')
print('b = ')
print(b)
V_cur = np.linalg.solve(A, b)

V_D = V_cur[1]
V_C = V_cur[2]
I   = -V_cur[3] / R
I_C = (V_cur[0] - V_cur[2]) / R_C
I_D = (V_cur[0] - V_cur[1]) / R_D
V_0C = V_cur[0] - V_cur[3]

print('')
print('V_cur = ', V_cur)
print('I = ', I)
print('I_D = ', I_D)
print('I_C = ', I_C)
print('V_D = ', V_D)
print('V_C = ', V_C)
print('V_0 = ', V_0C)
