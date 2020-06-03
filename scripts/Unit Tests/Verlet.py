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


# Work arrays
particles_cur_pos     = np.zeros((3, maxElec*3))
particles_cur_vel     = np.zeros((3, maxElec*3))
particles_cur_accel   = np.zeros((3, maxElec*3))
particles_prev_accel  = np.zeros((3, maxElec*3))
particles_prev2_accel = np.zeros((3, maxElec*3))
particles_mask        = np.zeros(maxElec*3, dtype=bool)

filename = 'Unit_Test_rand.bin'
dt = np.dtype([('x', np.float64), ('y', np.float64), ('z', np.float64), ('step', np.int32), ('species', np.int32)])
data = np.memmap(filename, dtype=dt, mode='r', order='F')

plot_data = True

if (plot_data == True):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

# -------------------------------------------------------------------------------------------------------------------------------------------------
# We read from a file the particles to emit in each time step
def Do_Emission(step):
    global nrElec

    for (x, y, z, emit_step, species) in data:
        if (step == emit_step):
            print('Particle {:d} at x = {:.2f}, y = {:.2f}, z = {:.2f}, {:d} {:d} added'.format(nrElec, x/length_scale, y/length_scale, z/length_scale, emit_step, species))

            # Position
            particles_cur_pos[0, nrElec] = x
            particles_cur_pos[1, nrElec] = y
            particles_cur_pos[2, nrElec] = z

            # Velocity
            particles_cur_vel[:, nrElec] = 0.0

            # Acceleration
            particles_cur_accel[:, nrElec] = 0.0
            particles_prev_accel[:, nrElec] = 0.0
            particles_prev2_accel[:, nrElec] = 0.0

            # Increase number of electrons
            nrElec = nrElec+1

    return None

# -------------------------------------------------------------------------------------------------------------------------------------------------
# Update the position
def Update_Position(step):
    global nrElec

    particles_mask[:] = False
    nrRemove = 0

    for i in range(nrElec):
        particles_cur_pos[:, i] = particles_cur_pos[:, i] + particles_cur_vel[:, i]*time_step + 0.5*particles_cur_accel[:, i]*time_step2

        # Store acceleration
        particles_prev2_accel[:, i] = particles_prev_accel[:, i]
        particles_prev_accel[:, i] = particles_cur_accel[:, i]
        particles_cur_accel[:, i] = 0.0

        if ((particles_cur_pos[2, i] < 0.0) or (particles_cur_pos[2, i] > d)):
            particles_mask[i] = True
            nrRemove = nrRemove + 1
            print('Particle {:d} at z = {:.2f} marked for removal'.format(i, particles_cur_pos[2, i]/length_scale))
        else:
            particles_mask[i] = False

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

    # Loop over all particles
    for i in range(nrElec):
        pos_1 = particles_cur_pos[:, i]

        particles_cur_accel[2, i] = -q_0*E/m_0

        # Loop over all particles and also all image charge particles
        for j in range(nrElec*3):
            pos_2 = particles_cur_pos[:, j]

            diff = pos_1 - pos_2
            r = np.sqrt( np.sum(diff**2) ) + length_scale**3 # + length_scale**3 to avoid division with zero
            force_c = pre_fac_c * diff / r**3

            particles_cur_accel[:, i] = particles_cur_accel[:, i] + force_c/m_0


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
        particles_cur_pos[0:1, i+maxElec] = particles_cur_pos[0:1, i]
        particles_cur_pos[2, i+maxElec] = -1.0*particles_cur_pos[2, i]

        # Mirror about z = d
        particles_cur_pos[0:1, i+2*maxElec] = particles_cur_pos[0:1, i]
        particles_cur_pos[2, i+2*maxElec] = 2*d - particles_cur_pos[2, i]

    return None

def Remove_Particles(step):
    global nrElec

    #for i in range(nrElec):
    #    if (particles_cur_pos[2, i] < 0.0 ):
    #        Remove_Particle_nr(i)

    #    if (particles_cur_pos[2, i] > d):
    #        Remove_Particle_nr(i)

    particles_cur_pos[0, :] = np.resize(particles_cur_pos[0, ~particles_mask], (maxElec*3))
    particles_cur_pos[1, :] = np.resize(particles_cur_pos[1, ~particles_mask], (maxElec*3))
    particles_cur_pos[2, :] = np.resize(particles_cur_pos[2, ~particles_mask], (maxElec*3))

    particles_cur_vel[0, :] = np.resize(particles_cur_vel[0, ~particles_mask], (maxElec*3))
    particles_cur_vel[1, :] = np.resize(particles_cur_vel[1, ~particles_mask], (maxElec*3))
    particles_cur_vel[2, :] = np.resize(particles_cur_vel[2, ~particles_mask], (maxElec*3))

    particles_cur_accel[0, :] = np.resize(particles_cur_accel[0, ~particles_mask], (maxElec*3))
    particles_cur_accel[1, :] = np.resize(particles_cur_accel[1, ~particles_mask], (maxElec*3))
    particles_cur_accel[2, :] = np.resize(particles_cur_accel[2, ~particles_mask], (maxElec*3))

    particles_prev_accel[0, :] = np.resize(particles_prev_accel[0, ~particles_mask], (maxElec*3))
    particles_prev_accel[1, :] = np.resize(particles_prev_accel[1, ~particles_mask], (maxElec*3))
    particles_prev_accel[2, :] = np.resize(particles_prev_accel[2, ~particles_mask], (maxElec*3))

    particles_prev2_accel[0, :] = np.resize(particles_prev2_accel[0, ~particles_mask], (maxElec*3))
    particles_prev2_accel[1, :] = np.resize(particles_prev2_accel[1, ~particles_mask], (maxElec*3))
    particles_prev2_accel[2, :] = np.resize(particles_prev2_accel[2, ~particles_mask], (maxElec*3))

    nrElec = nrElec - nrRemove
    return None

def Remove_Particle_nr(i):
    global nrElec

    print('Particle {:d} at z = {:.2f} removed'.format(i, particles_cur_pos[2, i]/length_scale))
    if (nrElec == 1):
        nrElec = 0
    else:
        particles_cur_pos[:, i] = particles_cur_pos[:, nrElec-1]
        particles_cur_vel[:, i] = particles_cur_vel[:, nrElec-1]
        particles_cur_accel[:, i] = particles_cur_accel[:, nrElec-1]
        particles_prev_accel[:, i] = particles_prev_accel[:, nrElec-1]
        particles_prev2_accel[:, i] = particles_prev2_accel[:, nrElec-1]

        nrElec = nrElec - 1

    return None

def Plot_Particles(step):
    global nrElec

    if (plot_data == True):
        ax.cla()
        ax.set_xlim3d(-50.0, 50.0)
        ax.set_ylim3d(-50.0, 50.0)
        ax.set_zlim3d(0.0, 1000.0)
        ax.scatter(particles_cur_pos[0, 0:nrElec]/length_scale, particles_cur_pos[1, 0:nrElec]/length_scale, particles_cur_pos[2, 0:nrElec]/length_scale, marker='o', c='blue')
        #plt.show()
        fig.canvas.draw()
        plt.pause(0.0001)
        #fig.canvas.flush_events()
    return None

# ----------------------------------------------------------------------------------------------------------
print('Starting')
for step in range(1, steps+1):
    print(step)
    Do_Emission(step)
    Update_Position(step)
    Remove_Particles(step)
    Plot_Particles(step)
    Update_Imagecharge_Positions(step)
    Calculate_Acceleration(step)
    Update_Velocity(step)
    print('')

print('Finished')
if (plot_data == True):
    plt.show(block=True)
