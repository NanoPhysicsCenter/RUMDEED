#!/bin/python3
#-------------------------------------#
#  Simple Molecular Dynamics          #
#  Using Verlet Algorithm             #
#  Kristinn Torfason                  #
#  08.09.23                           #
#-------------------------------------#

import numpy as np
import pandas as pd

# Parameters
time_scale = 1.0E-12  # Time scale 1 ps
length_scale = 1.0E-9 # Length scale 1 nm

dt = 0.5*time_scale   # Time step

N = 1000              # Number of particles
steps = 1000          # Number of steps

d = 100.0*length_scale # Gap spacing in z direction
L = 500.0*length_scale # Box length in x and y direction

# Constants
m = 9.1093837015E-31  # Mass of particles
q = 1.602176634E-19   # Charge of particles
epsilon_0 = 8.8541878128E-12 # Permittivity of free space
pre_factor = q**2/(4 * np.pi * epsilon_0) # Pre-factor for acceleration


# Position of particles
x = np.array([])
y = np.array([])
z = np.array([])

# Velocity of particles
vx = np.array([])
vy = np.array([])
vz = np.array([])

# Acceleration of particles
ax = np.array([])
ay = np.array([])
az = np.array([])
ax_new = np.array([])
ay_new = np.array([])
az_new = np.array([])


# Update positions
def Update_Positions():
    x = x + vx*dt + 0.5*ax*dt**2
    y = y + vy*dt + 0.5*ay*dt**2
    z = z + vz*dt + 0.5*az*dt**2

    return None

# Update velocities
def Update_Velocities():
    vx = vx + 0.5*(ax + ax_new)*dt
    vy = vy + 0.5*(ay + ay_new)*dt
    vz = vz + 0.5*(az + az_new)*dt

    return None

# Update accelerations
def Update_Accelerations():
    ax = ax_new
    ay = ay_new
    az = az_new

    return None

# Calculate accelerations using Coulomb's law
def Calculate_Accelerations():
    # Calculate acceleration for each particle
    for i in range(N):
        # Reset accelerations
        ax_new[i] = 0.0
        ay_new[i] = 0.0
        az_new[i] = 0.0

        # Loop over all other particles
        for j in range(N):
            if i != j: # Don't calculate acceleration of particle on itself
                r = np.sqrt((x[i]-x[j])**2 + (y[i]-y[j])**2 + (z[i]-z[j])**2) # Distance between particles
                
                # Calculate acceleration using Coulomb's law
                ax_new[i] += pre_factor * (x[i]-x[j])/(m*r**3)
                ay_new[i] += pre_factor * (y[i]-y[j])/(m*r**3)
                az_new[i] += pre_factor * (z[i]-z[j])/(m*r**3)
    
    return None

# Initialize positions
def Initialize():
    # Initialize positions
    x = L*(2*np.random.rand(N) - 1) # Random positions between -L and L
    y = L*(2*np.random.rand(N) - 1) # Random positions between -L and L
    d = np.random.rand(N)*d # Random positions between 0 and d

    return None

#-------------------------------------#
#  Main program                       #
#-------------------------------------#

# Initialize positions
Initialize()

for i in range(steps):
    # Calculate accelerations
    Calculate_Accelerations()

    # Update positions
    Update_Positions()

    # Calculate accelerations
    Calculate_Accelerations()

    # Update velocities
    Update_Velocities()

    # Update accelerations
    Update_Accelerations()

    # Print progress
    print('Step: ', i, ' of ', steps)

print('Hello World!')