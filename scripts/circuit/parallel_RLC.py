# Parallel RLC circuit
# Kristinn Torfason
# 20.12.2020

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import m_e, epsilon_0, e

#N = 1000000
N = 10000

time_step = 0.25E-15 # Sek
#time_step = 5.0E-6 # sek
time_step2 = time_step**2

A = (2500.0E-9)**2
V = 6.0
d = 1000.0E-9

J_CL = 4*epsilon_0/9*np.sqrt(2*e/m_e)*V**(3/2)/d**2
I_CL = A*J_CL

print('I_CL = {} A'.format(I_CL))

R = 1300.5 # Ohm
L = 500*1.04E-7 # Henry
C = 1/5*2.54E-14 # Farad

#R = 1.0E3 # Ohm
#L = 1.0 # H
#C = 10.0E-6 # Farad

w_0 = 1/np.sqrt(L*C)
Q = 1/R*np.sqrt(L/C)

f = w_0/(2*np.pi)
T = 2*np.pi/w_0

print('R = {} Ω'.format(R))
print('L = {} H'.format(L))
print('C = {} F'.format(C))

print('ω₀ = {} rad/s'.format(w_0))
print('f = {} Hz'.format(f))
print('T = {} s'.format(T))
print('Q = {}'.format(Q))

V_arr = np.zeros(N)
I_arr = np.zeros(N)
t_arr = np.zeros(N)
N_arr = np.zeros(N)

def V_next(I_cur, I_prev, V_cur, V_prev):
    V_n = time_step/C * I_cur + V_cur * (2.0 + time_step/(R*C) - time_step2/(L*C)) - V_prev - time_step/C * I_prev
    V_n = V_n / (1.0 + time_step/(R*C))

    return V_n

def I_next(t):
    #I_n = 10.0E-3 # Amper
    I_n = I_CL
    #I_n = 10.0E-3*np.sin(2*np.pi/(T)*t)
    return I_n


V_cur = 0.0
V_prev = 0.0
I_cur = 0.0
I_prev = 0.0
for k in range(N):
    t = k*time_step

    I = I_next(t)
    I_prev = I_cur
    I_cur = I

    V = V_next(I_cur, I_prev, V_cur, V_prev)

    V_prev = V_cur
    V_cur = V

    V_arr[k] = V
    I_arr[k] = I
    N_arr[k] = k
    t_arr[k] = t


print()

plt.plot(N_arr, V_arr)
plt.show()

print(V_arr[0:10])