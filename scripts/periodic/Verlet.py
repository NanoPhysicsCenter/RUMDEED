# Script taken from Unit Tests to test periodic model

import numpy as np
#import pandas as pd

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Constants
length_scale = 1.0E-9 # Length scale (1 nanometer)
time_scale = 1.0E-12 # Time scale (1 ps)
m_0 = 9.10938356E-31
q_0 = 1.6021766208E-19
pi = 3.14159265358979324 # Pi
h = 6.626070040E-34 # Planck's constant, h (Js)
k_b = 1.38064852E-23 # Boltzmann constant (J/K)
c = 299792458.0 # Speed of light (m/s)
mu_0 = 4.0*pi * 1.0E-7 # Vacuum permeability (H/m)
epsilon_0 = 1.0/(mu_0 * c**2) # Vacuum permittivity (F/m)

# Numerical parameters
steps = 250
maxElec = 1000
nrElec = 0
nrRemove = 0
time_step = 0.25E-3*time_scale
time_step2 = time_step**2

# System parameters
V_0 = 1000.0 # Voltage
d = 1000.0*length_scale # Gap spacing
E = -V_0/d # Vacuum field
L = 100.0*length_scale # Periodic
num_per = 1


# Work arrays
particles_cur_pos     = np.zeros((3, maxElec*3))
particles_cur_vel     = np.zeros((3, maxElec*3))
particles_cur_accel   = np.zeros((3, maxElec*3))
particles_prev_accel  = np.zeros((3, maxElec*3))
particles_prev2_accel = np.zeros((3, maxElec*3))
particles_mask        = np.zeros(maxElec*3, dtype=bool)

#filename = 'Unit_Test_rand.bin'
#dt = np.dtype([('x', np.float64), ('y', np.float64), ('z', np.float64), ('step', np.int32), ('species', np.int32)])
#data = np.memmap(filename, dtype=dt, mode='r', order='F')

plot_data = False

if (plot_data == True):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

# -------------------------------------------------------------------------------------------------------------------------------------------------
# We read from a file the particles to emit in each time step
def Do_Emission(step):
    global nrElec

    x_0 = 5.5*length_scale
    y_0 = -70.4*length_scale
    z_0 = 36.2*length_scale

    for i in range(-num_per, num_per+1):
        for j in range(-num_per, num_per+1):
            particles_cur_pos[0, nrElec] = x_0 + i*num_per
            particles_cur_pos[1, nrElec] = y_0 + j*num_per
            particles_cur_pos[2, nrElec] = z_0


    # for (x, y, z, emit_step, species) in data:
    #     #print('emit_step = {:d}, step = {:d}'.format(emit_step, step))
    #     if (step == emit_step):
    #         #print('Particle {:d} at x = {:.2f}, y = {:.2f}, z = {:.2f}, {:d} {:d} added'.format(nrElec, x/length_scale, y/length_scale, z/length_scale, emit_step, species))

    #         # Position
    #         particles_cur_pos[0, nrElec] = x.copy()
    #         particles_cur_pos[1, nrElec] = y.copy()
    #         particles_cur_pos[2, nrElec] = z.copy()

    #         # Velocity
    #         particles_cur_vel[:, nrElec] = 0.0

    #         # Acceleration
    #         particles_cur_accel[:, nrElec] = 0.0
    #         particles_prev_accel[:, nrElec] = 0.0
    #         particles_prev2_accel[:, nrElec] = 0.0

    #         # Increase number of electrons
    #         nrElec = nrElec+1

    return None

# -------------------------------------------------------------------------------------------------------------------------------------------------
# Update the position
def Update_Position(step):
    global nrElec, nrRemove

    particles_mask[:] = True # Set all to be removed (Including image charge partners).
    nrRemove = 0

    for i in range(nrElec):
        particles_cur_pos[:, i] = particles_cur_pos[:, i] + particles_cur_vel[:, i]*time_step + 0.5*particles_cur_accel[:, i]*time_step2

        # Store acceleration
        particles_prev2_accel[:, i] = particles_prev_accel[:, i].copy()
        particles_prev_accel[:, i] = particles_cur_accel[:, i].copy()
        particles_cur_accel[:, i] = 0.0

        if ((particles_cur_pos[2, i] < 0.0) or (particles_cur_pos[2, i] > d)):
            particles_mask[i] = True
            nrRemove = nrRemove + 1
            #print('Particle {:d} at z = {:.2f} marked for removal'.format(i, particles_cur_pos[2, i]/length_scale))
        else:
            particles_mask[i] = False

        if (step == 102):
            if (i == 100):
                print('Particle 100 at step 102')
                print(particles_cur_pos[:, i]/length_scale)
                input()

        #x = particles_cur_pos[0, i]/length_scale
        #y = particles_cur_pos[1, i]/length_scale
        #z = particles_cur_pos[2, i]/length_scale

        #print('New position for {:d} x = {:.2f}, y = {:.2f}, z = {:.2f}'.format(i, x, y, z))
        #print(particles_cur_vel[:, i]*time_step)
        #print(0.5*particles_cur_accel[:, i]*time_step2)

    return None

def Calculate_Acceleration(step):
    global nrElec

    pre_fac_c = (-q_0)*(-q_0) / (4 * pi * epsilon_0)

    #print('Calculate Acceleration nrElec = {:d}'.format(nrElec))
    #for i in range(nrElec):
        #print('i = {:d}'.format(i))
        #print(particles_cur_pos[:, i]/length_scale)
    
    #print('')
    # Loop over all particles
    for i in range(nrElec):
        if (particles_mask[i] == True):
            continue
        pos_1 = particles_cur_pos[:, i]

        particles_cur_accel[0:2, i] = 0.0
        particles_cur_accel[2, i] = -q_0*E/m_0
        #print('Vacuum Accel')
        #print(particles_cur_accel[:, i])

        # Loop over all particles and also all image charge particles
        for j in range(nrElec*3):
            if (i == j):
                continue # Skip this
            if (j == (i+nrElec)): # Skip image charge on self
                continue
            if (j == (i+2*nrElec)): # Skip image charge on self
                continue
            if (particles_mask[j] == True) and (j < nrElec):
                continue
            if ((particles_mask[i] == True) and ( (j == i+nrElec) or (j == i+2*nrElec))):
                continue

            pos_2 = particles_cur_pos[:, j]

            diff = pos_1 - pos_2
            r = np.sqrt( np.sum(diff**2) ) + length_scale**3 # + length_scale**3 to avoid division with zero
            force_c = pre_fac_c * diff / r**3


            if (j >= nrElec): # Check if this is an image charge, they have +q
                force_c = (-1.0)*force_c # Flip the sign if image charge

            particles_cur_accel[:, i] = particles_cur_accel[:, i] + force_c/m_0
            if (np.any(np.isinf(particles_cur_accel[:,i])) == True):
                print('INF')
                print('i = {:d}, j = {:d}'.format(i, j))
                print('pos_1 [nm]')
                print(pos_1/length_scale)
                print('pos_2 [nm]')
                print(pos_2/length_scale)
                print('Accel')
                print(force_c/m_0)
                print('')

            if (np.linalg.norm(particles_cur_accel) < 1.0E-3):
                print('NORM IS SMALL')
                print('i = {:d}, j = {:d}'.format(i, j))
                print('pos_1 [nm]')
                print(pos_1/length_scale)
                print('pos_2 [nm]')
                print(pos_2/length_scale)
                print('Accel')
                print(force_c/m_0)
                print('')
                input()

    return None

def Update_Velocity(step):
    global nrElec

    for i in range(nrElec):
        # Beemann
        particles_cur_vel[:, i] = particles_cur_vel[:, i] + 1.0/6.0*( 2.0*particles_cur_accel[:, i] + 5.0*particles_prev_accel[:, i] - particles_prev2_accel[:, i] )*time_step
    return None

def Update_Imagecharge_Positions(step):
    global nrElec

    for i in range(nrElec):
        # Mirror about z = 0
        particles_cur_pos[0:2, i+nrElec] = particles_cur_pos[0:2, i].copy()
        particles_cur_pos[2, i+nrElec] = -1.0*particles_cur_pos[2, i].copy()

        # Mirror about z = d
        particles_cur_pos[0:2, i+2*nrElec] = particles_cur_pos[0:2, i].copy()
        particles_cur_pos[2, i+2*nrElec] = 2*d - particles_cur_pos[2, i].copy()

    return None

def Remove_Particles(step):
    global nrElec, nrRemove

    #for i in range(nrElec):
    #    if (particles_cur_pos[2, i] < 0.0 ):
    #        Remove_Particle_nr(i)

    #    if (particles_cur_pos[2, i] > d):
    #        Remove_Particle_nr(i)

    if (nrRemove > 0):
        #print('Removing particles {:d}'.format(nrRemove))
        nrElec = nrElec - nrRemove

        # np.extract works in the same way as the pack command in Fortran
        particles_cur_pos[0, 0:nrElec] = np.extract(~particles_mask, particles_cur_pos[0, :])
        particles_cur_pos[1, 0:nrElec] = np.extract(~particles_mask, particles_cur_pos[1, :])
        particles_cur_pos[2, 0:nrElec] = np.extract(~particles_mask, particles_cur_pos[2, :])

        particles_cur_vel[0, 0:nrElec] = np.extract(~particles_mask, particles_cur_vel[0, :])
        particles_cur_vel[1, 0:nrElec] = np.extract(~particles_mask, particles_cur_vel[1, :])
        particles_cur_vel[2, 0:nrElec] = np.extract(~particles_mask, particles_cur_vel[2, :])

        particles_cur_accel[0, 0:nrElec] = np.extract(~particles_mask, particles_cur_accel[0, :])
        particles_cur_accel[1, 0:nrElec] = np.extract(~particles_mask, particles_cur_accel[1, :])
        particles_cur_accel[2, 0:nrElec] = np.extract(~particles_mask, particles_cur_accel[2, :])

        particles_prev_accel[0, 0:nrElec] = np.extract(~particles_mask, particles_prev_accel[0, :])
        particles_prev_accel[1, 0:nrElec] = np.extract(~particles_mask, particles_prev_accel[1, :])
        particles_prev_accel[2, 0:nrElec] = np.extract(~particles_mask, particles_prev_accel[2, :])

        particles_prev2_accel[0, 0:nrElec] = np.extract(~particles_mask, particles_prev2_accel[0, :])
        particles_prev2_accel[1, 0:nrElec] = np.extract(~particles_mask, particles_prev2_accel[1, :])
        particles_prev2_accel[2, 0:nrElec] = np.extract(~particles_mask, particles_prev2_accel[2, :])

    return None

def Plot_Particles(step):
    global nrElec

    if (plot_data == True):
        ax.cla()
        ax.set_xlim3d(-50.0, 50.0)
        ax.set_ylim3d(-50.0, 50.0)
        ax.set_zlim3d(0.0, 1000.0)
        ax.view_init(elev=0.0, azim=0.0)
        ax.scatter(particles_cur_pos[0, 0:nrElec]/length_scale, particles_cur_pos[1, 0:nrElec]/length_scale, particles_cur_pos[2, 0:nrElec]/length_scale, marker='o', c='blue')
        #plt.show()
        fig.canvas.draw()
        plt.pause(0.0001)
        #fig.canvas.flush_events()
    return None

def Compair_With_Fortran(step):
    filename_accel = 'accel/accel-{:d}.bin'.format(step) # The filename to read
    filename_pos = 'pos/pos-{:d}.bin'.format(step)

    dt_accel = np.dtype([('x', np.float64), ('y', np.float64), ('z', np.float64)])
    data_accel = np.memmap(filename_accel, dtype=dt_accel, mode='r', order='F')

    dt_pos = np.dtype([('x', np.float64), ('y', np.float64), ('z', np.float64)])
    data_pos = np.memmap(filename_pos, dtype=dt_pos, mode='r', order='F')

    #print('Accel diff')
    i = 0
    for (fortran_data) in zip(data_accel, data_pos):
        axyz_F = np.array(fortran_data[0].tolist())
        pxyz_F = np.array(fortran_data[1].tolist())

        if (i >= nrElec):
            print('Error number of particles is inconsitant')
            raise RuntimeError
        axyz_P = particles_cur_accel[:, i]
        pxyz_P = particles_cur_pos[:, i]
        i = i + 1
        #print(axyz_F)
        #print(axyz_P)
        err_a = np.max(np.divide(np.abs(axyz_F - axyz_P), axyz_P))
        err_p = np.max(np.divide(np.abs(pxyz_F - pxyz_P), pxyz_P))
        if (np.isinf(err_a) == True):
            print('err_a is inf')
            print('nrElec = {:d}, step = {:d}, i = {:d}'.format(nrElec, step, i))
            print(particles_cur_pos[:, i]/length_scale)
            print(particles_cur_vel[:, i])
            print(particles_cur_accel[:, i])
            input()
        #print(err)
        if (err_a > 1.0E-1):
            print('Acceleration Diffrance is more than 0.1')
            print(err_a)
            print(axyz_F)
            print(axyz_P)
            #input()
            print('')
        if (err_p > 1.0E-3):
            print('Position Diffrance is more than 0.001')
            print(err_p)
            print(pxyz_F)
            print(pxyz_P)
            input()
            print('')
        #print('')

    #print('')
    return None

# ----------------------------------------------------------------------------------------------------------
print('Starting')
for step in range(1, steps+1):
    #print(step)
    Do_Emission(step)
    Update_Position(step)
    Update_Imagecharge_Positions(step)
    Calculate_Acceleration(step)
    Update_Velocity(step)
    #Compair_With_Fortran(step)
    Remove_Particles(step)
    Plot_Particles(step)
    #print('')

print('Finished')
if (plot_data == True):
    plt.show(block=True)
