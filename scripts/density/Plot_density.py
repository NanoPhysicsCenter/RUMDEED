# Plot the density_emit.bin and density_absorb_top.bin files
# Kristinn Torfason
# 16.01.2020

import sys
import os.path as path
import matplotlib.pyplot as plt
import seaborn as sns

# The file layouts live in scripts/python_package/rumdeed_io.py
sys.path.insert(0, path.join(path.dirname(path.abspath(__file__)), '..', 'python_package'))
import rumdeed_io

# Directory with the output files, defaults to out/ in the current directory.
# Usage: python Plot_density.py [path/to/out]
filepath = sys.argv[1] if len(sys.argv) > 1 else 'out'

# Full path with filename
filename_emit = path.join(filepath, 'density_emit.bin') # density_emit.bin
filename_abs = path.join(filepath, 'density_absorb_top.bin') # density_absorb_top.bin

# Read the data from the files
df_density_emit = rumdeed_io.read_density_emit(filename_emit)
df_density_abs = rumdeed_io.read_density_absorb_top(filename_abs)

# Plot emission data
g = sns.JointGrid(x="x", y="y", data=df_density_emit, space=0.125, height=10)
g = g.plot_joint(sns.scatterplot)
g = g.plot_joint(sns.kdeplot, cmap="Blues_d", alpha=0.75)
g = g.plot_marginals(sns.kdeplot, fill=True)
g.set_axis_labels('x [nm]', "y [nm]")
plt.show()

# Plot absorption data
g = sns.JointGrid(x="x", y="y", data=df_density_abs, space=0.125, height=10)
g = g.plot_joint(sns.scatterplot)
g = g.plot_joint(sns.kdeplot, cmap="Blues_d", alpha=0.75)
g = g.plot_marginals(sns.kdeplot, fill=True)
g.set_axis_labels('x [nm]', "y [nm]")
plt.show()
