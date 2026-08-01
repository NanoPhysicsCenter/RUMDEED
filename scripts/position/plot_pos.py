# Plot the position file
# Kristinn Torfason
# 05.06.2018
# Mod. Hákon Örn 05.06.2021

import sys
import os.path as path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# The file layouts live in scripts/python_package/rumdeed_io.py
sys.path.insert(0, path.join(path.dirname(path.abspath(__file__)), '..', 'python_package'))
import rumdeed_io

# Path to the position file, defaults to out/position.bin in the current directory.
# Usage: python plot_pos.py [path/to/position.bin]
filename = sys.argv[1] if len(sys.argv) > 1 else "out/position.bin"

with open(filename, 'rb') as f:
    # Header is two int32's: total number of time steps and the output interval N_steps
    max_steps, N_steps = rumdeed_io.read_position_header(f)
    print('max_steps')
    print(max_steps)
    print('N_steps')
    print(N_steps)

    #------------------------------------------------------------
    # set up figure and animation
    plt.close("all")
    fig = plt.figure()

    # Make x and y subplot
    pad = 500
    #ax1 = fig.add_subplot(121, aspect='equal', autoscale_on=True, xlim=(-500, 500+pad), ylim=(-500, 500+pad))
    # Emission tip plot
    ax1 = fig.add_subplot(121, aspect='equal', autoscale_on=True, xlim=(-1000, 1000), ylim=(-1000, 1000))
    plt.xlabel('x [nm]')
    plt.ylabel('y [nm]')

    # Make r and z subplot
    #ax2 = fig.add_subplot(122, autoscale_on=True, xlim=(-500, 500+pad), ylim=(-1, 2501))
    # Plot for emission tip
    ax2 = fig.add_subplot(122, autoscale_on=True, xlim=(-1000, 1000), ylim=(-1, 2500))
    
    plt.xlabel('x [nm]')
    plt.ylabel('z [nm]')

    # Set the layout to tight
    fig.tight_layout()

    # particles1 and 2 hold the positions of the particles in (x, y) and (r, z) respectively
    particles1, = ax1.plot([], [], 'ko', ms=1) # x and y
    particles2, = ax2.plot([], [], 'ko', ms=1) # r and z

    # Draw the cylinders
    #rect    = plt.Rectangle((0, 1), 2, 3, fc='none')
    #circle1 = plt.Circle((0, 0), 1, color='r', fill=False)
    #circle2 = plt.Circle((0, 0), 2, color='b', fill=False)

    #ax1.add_patch(circle1)
    #ax1.add_patch(circle2)
    #ax2.add_patch(rect)

    # Init function, set data to empty
    def init():
        particles1.set_data([], [])
        particles2.set_data([], [])
        #return particles2 #, circle1, circle2, rect
        return particles1, particles2 #, circle1, circle2, rect

    # Read in data for time step i
    def animate(i):
        step, particles = rumdeed_io.read_position_step(f)

        # Scale to nm
        x = particles['x'] / 1.0E-9
        y = particles['y'] / 1.0E-9
        z = particles['z'] / 1.0E-9

        #r = np.sqrt(np.square(x) + np.square(y))
        particles1.set_data(x, y)
        particles2.set_data(x, z)
        #return particles2 #, circle1, circle2, rect
        return particles1, particles2 #, circle1, circle2, rect

    # Create the animation
    anim = animation.FuncAnimation(fig, animate, frames=max_steps//N_steps, interval=1, blit=True, repeat=False, init_func=init)

    # Show the animation
    plt.show()

    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=100, metadata=dict(artist='Me')) #bitrate=1800)

    anim.save('vid.mp4', writer=writer, dpi=100)
