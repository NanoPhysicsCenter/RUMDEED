# Kristinn Torfason
# 21.03.2020
# Script to calculate emittance and brightness

import numpy as np
import pandas as pd
import Vacuum
import os.path as path
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

filepath = './'
filename = path.join(filepath, 'density_absorb_top.bin')

# Binary file layout
# float64 (double precision numbers)
# int32 (32bit integers)
dt = np.dtype([('x', np.float64), ('y', np.float64), ('vx', np.float64), ('vy', np.float64), ('vz', np.float64), ('emit', np.int32), ('sec', np.int32)])

# Memory map the file
# mode=r (Read only)
# order=F (Fortran style array)
emittance = np.array([])
sigma_w = np.array([])
sigma_wp = np.array([])
theta_ell = np.array([])
w_theta = np.array([0.00, 0.10, 0.20, 0.30, 0.40, 0.50])
data_mem = np.memmap(filename, dtype=dt, mode='r', order='F')

# Read the data into dataframe
df = pd.DataFrame.from_records(data=data_mem, columns=data_mem.dtype.names)

df['v'] = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2)

df["x'"] = (df['vx']/df['vz'])/1.0E-3
df['x'] = df['x']

df["y'"] = (df['vy']/df['vz'])/1.0E-3
df['y'] = df['y']

df["r"] = np.sqrt(df["x"]**2 + df["y"]**2)
df["r'"] = (np.sqrt(df['vx']**2 + df['vy']**2) / df['vz']) / 1.0E-3

e_x, sw, swp, th = Vacuum.Calc_Emittance(df, "x", "x'")
e_y, _, _, _ = Vacuum.Calc_Emittance(df, "y", "y'")

sigma_w = sw
sigma_wp = swp
theta_ell = th

print('Emittance in x-direction {:6.2f} nm-mrad'.format(e_x))
print('Emittance in y-direction {:6.2f} nm-mrad'.format(e_y))

# Read data for current
filename_ramo = path.join(filepath, 'ramo_current.dt') # Ramo current
df_cur = pd.read_csv(filepath_or_buffer=filename_ramo, index_col=1, delim_whitespace=True, \
                    header=None, names=['time', 'step', 'current', 'volt', 'nrPart', 'nrElec', 'nrHole', 'avg_mob', 'avg_speed', 'ramo_1', 'ramo_2'])
        
df_cur['cur_roll'] = df_cur['current'].rolling(1000).mean()
cur = df_cur['cur_roll'].iloc[-1]

print('Current is {:6.2f} mA'.format(cur/1.0E-3))

# Calculate brightness
B = 2*cur/(np.pi**2*e_x*e_y * 1.0E-12)
print('Brightness is {:6.2f} A/(m-rad)^2'.format(B))

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=8)
#plt.rcParams['figure.figsize'] = [10, 7.5]
#plt.rcParams['figure.figsize'] = [2.75, 2.0625]
plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
plt.rc('ytick', labelsize=8)    # fontsize of the tick labels
plt.rc('axes', labelsize=8)    # fontsize of the x and y labels

#plt.rcParams['axes.linewidth'] = 5

ell_0 = Ellipse([0.0, 0.0], 2.0*sigma_w, 2.0*sigma_wp, theta_ell, zorder=1, linewidth=0.25)
ax = plt.subplot(111)

for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(0.25)

ax.grid(linestyle='-', linewidth='0.25')

ax.plot(df['x'], df["x'"], 'o', marker='o', markersize=2.5, zorder=0)
plt.xlabel('$x$ [nm]')
plt.ylabel(r'$x^\prime$ [mrad]')
ax.add_artist(ell_0)
ell_0.set_facecolor('red')
ell_0.set_alpha(0.5)

plt.rcParams['savefig.pad_inches'] = 0.0
#plt.savefig('emittance.png', bbox_inches='tight', transparent=True, dpi=600)

plt.show()
