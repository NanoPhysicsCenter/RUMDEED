# Hákon Örn Árnason
# 19.07.22
# Run a series of photo emission simulations, writing the input, work and
# laser files for every combination of the parameters below.
#
# Usage: python Serial_run_script.py <run directory>
# The run directory must contain the RUMDEED.out executable. One
# sub-directory is created per run.

import sys
import numpy as np
import os
import f90nml
import shutil
import subprocess

if len(sys.argv) != 2:
    raise SystemExit('Usage: python Serial_run_script.py <run directory with RUMDEED.out>')
os.chdir(sys.argv[1])

# Write work function file
# Format:
# Type of work function method to use
# Size of the matrix / Size of the circle / Gaussian parameters / Vornoi pattern
# Work function in eV, as matrix, circle etc.
def write_w_theta(workfunction):
    filename = 'work'
    np.savetxt(filename, workfunction, fmt='%2.3f', header='1\n1 1', comments='') # header defines method and matrix size
    return None

# Write laser file
# Format:
# Gauss pulse; on = 1, off = 2
# Laser energy; Fixed = 1, Possion = 2
# Initial velocity; zero = 1, Work function dependant = 2
# Laser energy in eV ; std of laser
# Gauss pulse parameters: mu (center); sigma (width); Amplitude
def write_laser(laser):
    filename = 'laser'
    np.savetxt(filename, laser , fmt='%2.3f', header='2 1 2', comments='')
    return None

# Write Fortran input file
def Write_Input(in_volt):
    nml = f90nml.Namelist({'input': {'v_s': in_volt}})
    nml['input']['box_dim']             = [0.0, 0.0, 10000.0]    # Box dimentions in nano meters
    nml['input']['time_step']           = 0.25E-4                # Time of each step in pico seconds
    nml['input']['steps']               = 20000                  # Number of steps
    nml['input']['emission_mode']       = 1                      # 1 = Photo, 10 = Field emission (see Examples/README.md)
    nml['input']['NrEmit']              = 1                      # If more then 1, remember to add them below
    nml['input']['image_charge']        = True                   # Turn image charge calculations on/off
    nml['input']['N_IC_MAX']            = 1                      # Number of image charge partners to use
    nml['input']['collision_mode']      = 0                      # 0 = Off, 1 = Continuous ionization, ... (see Examples)
    nml['input']['T_temp']              = 293.15                 # Temperature for ion gas or cathode if thermal-field emission
    nml['input']['P_abs']               = 1.0                    # Pressure for ion gas
    nml['input']['write_position_file'] = True                   # Do you want the position file (it is HUGE)
    nml['input']['emitters_dim']        = [[250.0, 250.0, 0.0]]  # Emitter dimentions in nano meters | see docs
    nml['input']['emitters_pos']        = [[-250.0, -250.0, 0.0]] # Emitter position | rectangle: left bottom corner | circle: center
    nml['input']['emitters_type']       = [1]                    # 1 -> Circle; 2 -> Rectangle; 3 -> Rectangle with spots
    nml['input']['emitters_delay']      = [0]                    # Emission delay in time steps

    f90nml.write(nml, 'input', force=True)

    return None

# ----------------------------------------------------------------------------------------
# Input and laser file parameters | BEWARE: number of runs is a multiple of all input variables
diode_volt    = [100]    # Volt  | Voltage across the gap
mu            = 20000        # steps | Center position of pulse
sigma = [1000.000] #list(range(5, 100 + 1, 5))   # steps | Pulse width,
amplitude = [10.000] #10*[5] + 10*[10]  # a.u.  | Amplitude of gauss pulse
photon_energy = 4.800         # eV    | Energy of input laser
photon_std    = 0.020        # eV    | Standard deviation of input laser energy
wf_matrix_size = [1] # list(range(2, 10+1, 2))
run_num = 0 # Make unique number for this series of runs
w_base = 4.70
w_add = 2.0


# Define work function matrix
W = np.ones((1, 1)) * 4.70


# Run script
run_num = 0 # (int) Make unique start number for this series of runs
for k in range(len(diode_volt)):
    for i in range(len(sigma)):
        for j in range(len(amplitude)):
            print('')
            print('Creating directory')
            folder = str(diode_volt[k]) + 'V_' + str(sigma[i]) + 'Sigma_' + str(amplitude[j])[0:4]  + 'Amp_Run' + str(run_num)
            print(folder)
            os.makedirs(folder, exist_ok=True)
            os.chdir(folder)

            print('Writing input file')
            Write_Input(diode_volt[k])

            laser_input = np.array([photon_energy, photon_std, mu, sigma[i], amplitude[j]])

            print('Writing laser file')
            write_laser(laser_input)

            print('Writing work function file')
            write_w_theta(W)

            print('Running RUMDEED')
            shutil.copy2('../RUMDEED.out', '.')
            subprocess.run('./RUMDEED.out')
            os.chdir('..')
            run_num += 1

# Other useful definitions for work function patterns

def Half_fill(w_base, w_add, N):
    A = np.ones((N, N))*w_base
    for i in range(N):

        if (i % 2 == 0):
            j_start = 1
        else:
            j_start = 0

        for j in range(j_start, N, 2):
            A[i, j] = w_base + w_add

    return A
