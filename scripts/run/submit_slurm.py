#!/users/home/kristinnt/anaconda3/bin/python3

# JOB name
#SBATCH -J Vacuum-MD-Entropy_28

# Number of nodes
###SBATCH -N 1
#SBATCH --array=0-50

# Number of cores
#SBATCH --tasks-per-node=8
###SBATCH --exclusive

# Number of GPUs (must use GPU partition)
###SBATCH --gres=gpu:2

# Wall time 
#SBATCH --time=15-00:00:00

# Q
#SBATCH -p normal,omnip,himem

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=kristinnt@ru.is

import os
import f90nml
import numpy as np
import subprocess
import shutil

def Write_W(w):
    filename = 'w_theta'
    sh = w.shape
    np.savetxt(filename, w, fmt='%2.3f', header='1\n'+ str(sh[0]) + ' ' + str(sh[1]), comments='')
    return None

def Fill_One_Spot(W):
    # Check if we have filled all spots
    if (np.abs(W-2.0) < 1.0E-6).all():
        print('All spots are filled!!!')
        return W
    
    Found = False
    while (Found == False):
        # Find a spot
        i = int(np.rint((11- 0) * np.random.random_sample() + 0))
        j = int(np.rint((11 - 0) * np.random.random_sample() + 0))
        
        # Check if we have already change this spot
        if (np.abs(W[i, j] - 2.50) < 1.0E-6):
            Found = True # We found an unchanged spot
        else:
            Found = False # Look for another spot
    
    A = W.copy() # Make a copy of the input matrix
    A[i, j] = 2.0
    return A

def Fill_N_Spots(N, W):
    for i in range(N):
        W = Fill_One_Spot(W)

    return W

def Write_Input():
    nml = f90nml.Namelist({'input': {'v_s': 2000.0}})
    nml['input']['box_dim'] = [0.0, 0.0, 1000.0]
    nml['input']['time_step'] = 1.00E-4
    nml['input']['steps'] = 10000
    nml['input']['emission_mode'] = 99
    nml['input']['NrEmit'] = 1
    nml['input']['image_charge'] = True
    nml['input']['N_IC_MAX'] = 1
    nml['input']['emitters_dim'] = [[1000.0, 1000.0, 0.0]]
    nml['input']['emitters_pos'] = [[-500.0, -500.0, 0.0]]
    nml['input']['emitters_type'] = [2]
    nml['input']['emitters_delay'] = [0]

    f90nml.write(nml, 'input', force=True)
    
    return None


print('This job runs on')
print(os.environ['SLURM_JOB_NODELIST'])
print('Number of cores')
print(os.environ['SLURM_CPUS_ON_NODE'])
print('')


# Get variables
USER = os.environ['USER']
JOB_ID = os.getenv('SLURM_JOB_ID', '1')
TASK_ID = os.getenv('SLURM_ARRAY_TASK_ID', '')

# Create a folder for this task, copy the executable file and change into it
if (TASK_ID != ''):
  os.makedirs(TASK_ID)
  shutil.copy2('Vacuum-MD.out', TASK_ID)
  os.chdir(TASK_ID)
  ORG_PATH = os.getcwd()

# Make the scratch directory
SCRATCH_PATH = os.path.join('/scratch', USER, JOB_ID, TASK_ID)

try:
  print('Creating scratch directory ' + SCRATCH_PATH)
  os.makedirs(SCRATCH_PATH)
except FileExistsError:
  pass

# Copy the program into the scratch directory
shutil.copy2('Vacuum-MD.out', SCRATCH_PATH)

# Create the input files for the program in the scratch directory
W = np.ones((12, 12))*2.5

os.chdir(SCRATCH_PATH)
Write_Input()

#W = np.ones((5, 5))*2.50
#for i in range(int(TASK_ID)):
#  W = Fill_One_Spot(W)

A = Fill_N_Spots(28, W.copy())
Write_W(A)

# Start the program
my_env = os.environ.copy()
my_env['OMP_NUM_THREADS'] = '8' 
subprocess.run('./Vacuum-MD.out', env=my_env)

# Move data from the scratch directory
shutil.move(SCRATCH_PATH, ORG_PATH)

# All done
print('Done')
