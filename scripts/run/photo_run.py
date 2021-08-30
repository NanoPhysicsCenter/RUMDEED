# Hákon Örn Árnason
# 28.08.21

import numpy as np
import os
import f90nml
import shutil
import subprocess

# Write work function file
# Format:
# Type of work function method to use
# Size of the matrix / Size of the circle / Gaussian parameters / Vornoi pattern
# Work function in eV, as matrix, circle etc.
def write_w_theta(workfunction):
    filename = 'w_theta'
    np.savetxt(filename, workfunction, fmt='%2.3f', header='1\n1 1', comments='') # header defines method and matrix size
    return None

# Write laser file
# Format:
# Gauss pulse = 1, Continous = 2, ; Initial velocity, zero = 1, norm dist work function dependant = 2
# Laser energy in eV ; std of laser
# Gauss pulse parameters: mu (center); sigma (width); Amplitude
def write_laser(laser):
    filename = 'laser'
    np.savetxt(filename, laser, fmt='%2.3f', header='1 2\n4.7 0.02\n10000 1000 2', comments='')
    #https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
    return None

# Write Fortran input file
def Write_Input():
    nml = f90nml.Namelist({'input': {'v_s': 2500.0}})
    nml['input']['box_dim'] = [0.0, 0.0, 0.0]
    nml['input']['time_step'] = 0.25E-3
    nml['input']['steps'] = 60000
    nml['input']['emission_mode'] = 1
    nml['input']['NrEmit'] = 1
    nml['input']['image_charge'] = True
    nml['input']['emitters_dim'] = [[500.0, 500.0, 0.0]]
    nml['input']['emitters_pos'] = [[-250.0, -250.0, 0.0]]
    nml['input']['emitters_type'] = [2]
    nml['input']['emitters_delay'] = [0]

    f90nml.write(nml, 'input', force=True)
    
    return None


# ----------------------------------------------------------------------------------------
# #N = 12*12 + 1
# N = 2 # Number of runs
# folders = []
# for i in range(N):
#     folders.append(str(i))

# Define work function matrix
W = np.ones((1, 1))*4.70
# for i in range(N):
#     print('')
#     print('Creating directory')
#     print(folders[i])
#     try:
#         os.mkdir(folders[i])
#     except:
#         pass
#     os.chdir(folders[i])

#     print('Writing input file')
#     Write_Input()

#     print('Writing laser file')
#     write_laser(laser_pulse)

#     print('Writing work function file')
#     write_w_theta(flat_function)

#     print('Running Vacuum-MD')
#     shutil.copy2('../Vacuum-MD.out', '.')
#     subprocess.run('./Vacuum-MD.out')
#     os.chdir('..')
