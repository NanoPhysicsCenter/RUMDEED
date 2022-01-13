# Plot the current in ramo_current.dt
# Kristinn Torfason

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os.path as path
import f90nml
import sys
sys.path.append('/home/kristor/Work/Vacuum-MD/scripts/python_package')
import Vacuum

filepath = './out/'
filename_ramo = path.join(filepath, 'ramo_current.dt') # Ramo current

# Read the data into a pandas dataframe
#cur_time, step, ramo_cur, V_d, nrPart, nrElec, nrHole
df_ramo = pd.read_csv(filepath_or_buffer=filename_ramo, index_col=1, delim_whitespace=True, \
                        header=None, names=['time', 'step', 'current', 'volt', 'nrPart', 'nrElec', 'nrHole', 'avg_mob', 'avg_speed', 'ramo_1', 'ramo_2'])

plt.plot(df_ramo.time, df_ramo.current/1.0E-3) # Plot in mA and ps
plt.xlabel('Time [ps]')
plt.ylabel('Current [mA]')

plt.show()