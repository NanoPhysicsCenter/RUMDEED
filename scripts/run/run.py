# Kristinn Torfason
# 11.12.18

import numpy as np
import os
import f90nml
import shutil
import subprocess

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

def Fill_One_Spot(W):
    # Check if we have filled all spots
    if (np.abs(W-2.0) < 1.0E-6).all():
        print('All spots are filled!!!')
        return W
    
    Found = False
    while (Found == False):
        # Find a spot
        i = int(np.rint((10 - 0) * np.random.random_sample() + 0))
        j = int(np.rint((10 - 0) * np.random.random_sample() + 0))
        
        # Check if we have already change this spot
        if (np.abs(W[i, j] - 2.50) < 1.0E-6):
            Found = True # We found an unchanged spot
        else:
            Found = False # Look for another spot
    
    A = W.copy() # Make a copy of the input matrix
    A[i, j] = 2.0
    return A

def Write_W(W):
    filename = 'w_theta'
    np.savetxt(filename, W, fmt='%2.3f', header='1\n11 11', comments='')
    return None

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

# ----------------------------------------------------------------------------------------
N = 11*11 + 1
folders = []
for i in range(N):
    folders.append(str(i))

W = np.ones((11, 11))*2.50
A = W.copy()
for i in range(N):
    print('')
    print('Creating directory')
    print(folders[i])
    try:
        os.mkdir(folders[i])
    except:
        pass
    os.chdir(folders[i])

    print('Writing input file')
    Write_Input()

    print('Writing work function file')
    Write_W(A)

    print('Running Vacuum-MD')
    shutil.copy2('../Vacuum-MD.out', '.')
    subprocess.run('./Vacuum-MD.out')
    os.chdir('..')
    A = Fill_One_Spot(A)
