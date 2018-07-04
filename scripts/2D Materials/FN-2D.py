#------------------------------------------------------------------------------#
# 2D Materials Fowler-Nordheim law                                             #
# See Ang, Y. S., Zubair, M., Ooi, K. J. A., & Ang, L. K. (n.d.).              #
# Generalized Fowler-Nordheim Model of Field-Induced Vertical Electron         #
# Emission from Two-Dimensional Dirac Material.                                #
#                                                                              #
# Kristinn Torfason                                                            #
# 03.07.2018                                                                   #
#                                                                              #
#------------------------------------------------------------------------------#

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, m_e, pi, hbar

# V_⊥, L_⊥ and e_⊥
e_o = 38.3*e #Bound state energy [eV -> J]
v_o = np.sqrt(2.0*e_o/m_e) # Cross plane velocity [m/s]
L_o = 0.335E-9 # 2D material thickness [m] [0.335nm, 0.670nm]

g_sv = 4.0 # Spin-valley degeneracy []
m_eff = 0.03*m_e # Electron effective mass [kg]

e_f = 0.1*e # Fermi-energy [eV -> J]
#v_f = 10.0E6 # Fermi-velocity [m/s]

phi_B0 = 4.5*e # [eV -> J]
phi_B = phi_B0 - e_f # [eV -> J]
#F = 50.0*1.0/1.0E-9 # Electric field [V/m]
F_nm = np.linspace(0.1, 50.0, 100) # [V/nm]
F = F_nm * 1.0/(1.0E-9) # [V/nm * 1nm/1.0E-9m = V/m]

d_f = e*hbar*F/(2.0*np.sqrt(2*m_e*phi_B))
Df = np.exp(-2.0*phi_B/(3.0*d_f))

J = v_o/L_o * g_sv*e*m_eff/(2.0*pi*hbar**2) * e_f * Df * np.exp(-e_f/d_f)
J_cm = J * ( (1.0E-2)/(1.0) )**2

plt.figure()
plt.semilogy(F_nm, J_cm)
plt.xlabel('F [V/nm]')
plt.ylabel('J [A/cm^2]')
plt.xlim([0.0, 50.0])
plt.ylim([10.0E-3, 10.0E9])
#plt.show()

plt.figure()
plt.plot(1/F_nm, np.log(J/F**2))
plt.xlabel('1/F [V/nm]')
plt.ylabel('ln(J/F^2)')
plt.xlim([0.2, 0.4])
plt.ylim([-40, 0])
#plt.show()

plt.figure()
plt.plot(F_nm, Df*np.exp(-e_f/d_f))


plt.show()
