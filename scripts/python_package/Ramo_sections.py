import numpy as np
import pandas as pd
import os.path as path
import f90nml
import matplotlib.pyplot as plt

filepath = './'

# Directory of the data
filename_ramo = path.join(filepath, 'out/ramo_current.dt') # Ramo current

# Read the data into a pandas dataframe
#cur_time, step, ramo_cur, V_d, nrPart, nrElec, nrHole
df = pd.read_table(filepath_or_buffer=filename_ramo, index_col=1, delim_whitespace=True, \
                   header=None, names=['time', 'step', 'current', 'volt', 'nrPart', 'nrElec', 'nrHole', 'avg_mob', 'avg_speed', 'ramo_1', 'ramo_2'])

# Max and min time values
x_min = df['time'].min()
x_max = df['time'].max()

# Number of time steps
steps = df['time'].count()

# Read the input file for parameters used
filename_input = path.join(filepath, 'input') # input file
input_nml = f90nml.read(filename_input)

# Read the init.bin file for system parameters
#epsilon_r, m_eeff, m_heff, length_scale, time_scale, vel_scale, cur_scale, MAX_PARTICLES, MAX_EMITTERS, MAX_SECTIONS, MAX_LIFE_TIME
dt = np.dtype([('epsilon_r', np.float64), \
               ('m_eeff', np.float64), \
               ('m_heff', np.float64), \
               ('length_scale', np.float64), \
               ('time_scale', np.float64), \
               ('vel_scale', np.float64), \
               ('cur_scale', np.float64), \
               ('MAX_PARTICLES', np.int32), \
               ('MAX_EMITTERS', np.int32), \
               ('MAX_SECTIONS', np.int32), \
               ('MAX_LIFE_TIME', np.int32) ])
filename_initbin = path.join(filepath, 'out/init.bin')
data_sys = np.memmap(filename_initbin, dtype=dt, mode='r', order='F')
input_nml['system'] = dict(zip(data_sys.dtype.names, data_sys[0]))

# Set max parameters
MAX_SECTIONS = input_nml['system']['MAX_SECTIONS']
MAX_EMITTERS = input_nml['system']['MAX_EMITTERS']

# Read the ramo current binary file. That is the ramo current broken down into sections and emitters
filename_ramo_sec = path.join(filepath, 'out/ramo_current.bin') # Ramo current by emitters and sections

# Layout of the array is [Section, Emitter, time step]
data_sec = np.memmap(filename_ramo_sec, dtype=np.float64, mode='r', order='F', shape=(MAX_SECTIONS, MAX_EMITTERS, steps))

# In this case we have a center section with work function 4.60 eV and and outer section around it with 4.70 eV.
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

f = open("ramo_sections.dt", "w")

for step in range(steps):
    for sec in range(MAX_SECTIONS):
        if (sec == (MAX_SECTIONS-1)):
            f.write('{:f}'.format(data_sec[sec, 0, step]))
        else:
            f.write('{:f}\t'.format(data_sec[sec, 0, step]))
    f.write('\n')

f.close()

#np.savetxt(fname='ramo_sections.dt', X=data_sec_2D.transpose())

#plt.plot(data_sec[0, 0, :]/1.0E-3)
#plt.show()