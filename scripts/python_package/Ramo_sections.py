import numpy as np
import pandas as pd
import os.path as path
import f90nml
import matplotlib.pyplot as plt
import rumdeed_io

filepath = './'

# Directory of the data
filename_ramo = path.join(filepath, 'out/ramo_current.dt') # Ramo current

# Read the data into a pandas dataframe
df = rumdeed_io.read_ramo_current(filename_ramo)

# Max and min time values
x_min = df['time'].min()
x_max = df['time'].max()

# Number of time steps
steps = df['time'].count()

# Read the input file for parameters used
filename_input = path.join(filepath, 'input') # input file
input_nml = f90nml.read(filename_input)

# Read the init.bin file for system parameters
filename_initbin = path.join(filepath, 'out/init.bin')
input_nml['system'] = rumdeed_io.read_init(filename_initbin)

# Set max parameters
MAX_SECTIONS = int(input_nml['system']['MAX_SECTIONS'])
MAX_EMITTERS = int(input_nml['system']['MAX_EMITTERS'])

# Read the ramo current binary file. That is the ramo current broken down into
# sections and emitters. Layout of the array is [Section, Emitter, time step].
# Only written when write_ramo_sec = .true. in the input namelist.
filename_ramo_sec = path.join(filepath, 'out/ramo_current.bin') # Ramo current by emitters and sections
data_sec = rumdeed_io.read_ramo_sections(filename_ramo_sec, MAX_SECTIONS, MAX_EMITTERS, steps)

# Example from a run with 16 sections (MAX_SECTIONS in mod_global.F90 sets the actual number):
# a center section with work function 4.60 eV and and outer section around it with 4.70 eV.
# Sections 5, 6, 9, 10 are in the center
#|---|---|---|---|
#| 0 | 1 | 2 | 3 |
#|---|---|---|---|
#| 4 | 5 | 6 | 7 |
#|---|---|---|---|
#| 8 | 9 |10 |11 |
#|---|---|---|---|
#|12 |13 |14 |15 |
#|---|---|---|---|
# In the Fortran program the first section is labeled 1 and the last one 16.

data_sec_2D = data_sec[:, 0, :]

# One row per time step, one tab-separated column per section
np.savetxt(fname='ramo_sections.dt', X=data_sec_2D.transpose(), fmt='%f', delimiter='\t')

#plt.plot(data_sec[0, 0, :]/1.0E-3)
#plt.show()