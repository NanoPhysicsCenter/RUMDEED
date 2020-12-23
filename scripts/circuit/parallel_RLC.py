# Parallel RLC circuit
# Kristinn Torfason
# 20.12.2020

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import m_e, epsilon_0, e

#N = 1000000
N = 10000

time_scale = 1.0E-12
time_step = 0.25*time_scale # Sek
#time_step = 5.0E-6 # sek
time_step2 = time_step**2

A = (100.0E-9)**2
V_s = 1.0
d = 1000.0E-9

L = 1.04E-9 # Henry
C = 5.52E-17 # Farad
#C = epsilon_0*A/d
R = 1000.0
#R = 0.5*np.sqrt(L/C)*1.5 # Ohm
#R = 2000.0

#R = 1.0E3 # Ohm
#L = 1.0 # H
#C = 10.0E-6 # Farad

# ##### Luginsland
# A = 4.7E-6
# Q = 150.0
# R = 13.5
# w_0 = 2*np.pi*3.1E9
# d = 3.6E-2
# V_s = 250.0E3

# C = Q/(R*w_0)
# L = 1.0/(C*w_0**2)


############

J_CL = 4*epsilon_0/9*np.sqrt(2*e/m_e)*V_s**(3/2)/d**2
I_CL = A*J_CL

w_0 = 1/np.sqrt(L*C)
Q = R*np.sqrt(C/L)
alpha = 1.0/(2*R*C)

if (alpha**2 < w_0**2):
    print('Underdamped')

if (abs(alpha**2 - w_0**2) < 1.0E-6):
    print('Critically damped')

if (alpha**2 > w_0**2):
    print('Overdamped')

f = w_0/(2*np.pi)
T = 2*np.pi/w_0

print('R = {} Ω'.format(R))
print('L = {} H'.format(L))
print('C = {} F'.format(C))

print('')
print('ω₀ = {} rad/s'.format(w_0))
print('α = {} rad/s'.format(alpha))
print('f = {} Hz'.format(f))
print('T = {} s'.format(T))
print('Q = {}'.format(Q))

print('')
print('I_CL = {} A'.format(I_CL))

V_arr = np.zeros(N)
I_arr = np.zeros(N)
t_arr = np.zeros(N)
N_arr = np.zeros(N)

def V_next(I_cur, I_prev, V_cur, V_prev):
    V_n = time_step/C * I_cur + V_cur * (2.0 + time_step/(R*C) - time_step2/(L*C)) - V_prev - time_step/C * I_prev
    V_n = V_n / (1.0 + time_step/(R*C))

    return V_n

def Child_Lang_Current(V_):
    return A*4*epsilon_0/9*np.sqrt(2*e/m_e)*V_**(3/2)/d**2

def I_next(t, V_cur):
    #I_n = 10.0E-3 # Amper
    if (V_cur <= 0):
        I_n = 0.0
    else:
        I_n = Child_Lang_Current(V_cur)
    #I_n = 10.0E-3*np.sin(2*np.pi/(T)*t)
    #I_n = 0.0
    return I_n


V_cur = V_s
V_prev = 0.0
I_cur = 0.0
I_prev = 0.0
for k in range(N):
    t = k*time_step

    I = I_next(t, V_cur+V_s)
    I_prev = I_cur
    I_cur = I

    V = V_next(I_cur, I_prev, V_cur, V_prev)
    #print(V)

    V_prev = V_cur
    V_cur = V

    V_arr[k] = V
    I_arr[k] = I
    N_arr[k] = k
    t_arr[k] = t


print()

plt.plot(t_arr/time_scale, V_arr+V_s)
plt.show()

print(V_arr[0:10])