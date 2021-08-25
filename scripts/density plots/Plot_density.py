# Plot the density_emit.bin and density_absorb_top.bin files
# Kristinn Torfason
# 16.01.2020

import numpy as np
import pandas as pd
import os.path as path
import matplotlib.pyplot as plt
import seaborn as sns

# Path to the folder with the density_emit.bin and density_absorb_top.bin files
filepath = '/home/hakon/Documents/FE Brynjar/input7/out/'

# Full path with filename
filename_emit = path.join(filepath, 'density_emit.bin') # density_emit.bin
filename_abs = path.join(filepath, 'density_absorb_top.bin') # density_absorb_top.bin

# File maps
dt_layout_emit = [('x', np.float64), \
                  ('y', np.float64), \
                  ('emit', np.int32), \
                  ('sec', np.int32)]

dt_layout_abs = [('x', np.float64), \
                 ('y', np.float64), \
                 ('v_x', np.float64),\
                 ('v_y', np.float64),\
                 ('v_z', np.float64),\
                 ('emit', np.int32),\
                 ('sec', np.int32)]

data_mem_emit = np.memmap(filename_emit, dtype=dt_layout_emit, mode='r', order='F')
data_mem_abs = np.memmap(filename_abs, dtype=dt_layout_abs, mode='r', order='F')

# Read the data from the files
df_density_emit = pd.DataFrame.from_records(data=data_mem_emit, columns=data_mem_emit.dtype.names)
df_density_abs = pd.DataFrame.from_records(data=data_mem_abs, columns=data_mem_abs.dtype.names)

# Plot emission data
g = sns.JointGrid(x="x", y="y", data=df_density_emit, space=0.125, height=10)
g = g.plot_joint(sns.scatterplot)
g = g.plot_joint(sns.kdeplot, cmap="Blues_d", alpha=0.75)
g = g.plot_marginals(sns.kdeplot, shade=True)
g.set_axis_labels('x [nm]', "y [nm]")
plt.show()

# Plot absorption data
g = sns.JointGrid(x="x", y="y", data=df_density_abs, space=0.125, height=10)
g = g.plot_joint(sns.scatterplot)
g = g.plot_joint(sns.kdeplot, cmap="Blues_d", alpha=0.75)
g = g.plot_marginals(sns.kdeplot, shade=True)
g.set_axis_labels('x [nm]', "y [nm]")
plt.show()