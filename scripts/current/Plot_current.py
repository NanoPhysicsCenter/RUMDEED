# Plot the current in ramo_current.dt
# Kristinn Torfason

import sys
import os.path as path
import matplotlib.pyplot as plt

# The file layouts live in scripts/python_package/rumdeed_io.py
sys.path.insert(0, path.join(path.dirname(path.abspath(__file__)), '..', 'python_package'))
import rumdeed_io

# Directory with the output files, defaults to out/ in the current directory.
# Usage: python Plot_current.py [path/to/out]
filepath = sys.argv[1] if len(sys.argv) > 1 else 'out'
filename_ramo = path.join(filepath, 'ramo_current.dt') # Ramo current

# Read the data into a pandas dataframe
df_ramo = rumdeed_io.read_ramo_current(filename_ramo)

plt.plot(df_ramo.time, df_ramo.current/1.0E-3) # Plot in mA and ps
plt.xlabel('Time [ps]')
plt.ylabel('Current [mA]')

plt.show()
