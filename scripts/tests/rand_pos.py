#!/usr/bin/env python3
# Unit test many particles
# This will generate a random position file for the inital test

import numpy as np

nr = 1000 # Number of particles

length_scale = 1.0E-9
length_x = 100.0*length_scale
length_y = 100.0*length_scale
length_z = 100.0*length_scale

# Get random numbers for nr number of particles
# for each space dimension x, y, z
print('Generating random positions')
pos = np.random.rand(nr, 3)

pos[:, 0] = pos[:, 0]*length_x
pos[:, 1] = pos[:, 1]*length_y
pos[:, 2] = pos[:, 2]*length_z

print('Saving to file')
np.savetxt(fname='rand_pos_init.dt', X=pos, fmt='%.18e', delimiter='    ')

print('Done')
