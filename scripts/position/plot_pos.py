# Plot the position file
# Kristinn Torfason
# 05.06.2018
# Mod. Hákon Örn 05.06.2021

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

with open("/home/hakon/Documents/PE Simulations/Aug10/10sigma/19 (1)/out/position.bin", 'rb') as f:
    max_steps = np.fromfile(file=f, count=1, dtype=np.int32)
    max_steps = max_steps[0]
    print('max_steps')
    print(max_steps)


    # for i in range(steps):
    #     step, nrPart = np.fromfile(file=f, count=2, dtype=np.int32)
    #     print(step)
    #     print(nrPart)
    #     for j in range(nrPart):
    #         x, y, z = np.fromfile(file=f, count=3, dtype=np.float64)
    #         emit = np.fromfile(file=f, count=1, dtype=np.int32)
    #         print(x/1.0E-9)
    #         print(y/1.0E-9)
    #         print(z/1.0E-9)
    #         print(emit)
    #         print('')
    #         plt.plot(x/1.0E-9, z/1.0E-9, 'r.')
    #         plt.xlim([-200.0, 200.0])
    #         plt.ylim([0.0, 1500.0])
    #
    #     plt.show()



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
        step, nrPart = np.fromfile(file=f, count=2, dtype=np.int32)

        x = np.zeros(nrPart)
        y = np.zeros(nrPart)
        z = np.zeros(nrPart)
        #r = np.zeros(nrPart)

        for j in range(nrPart):
            x[j], y[j], z[j] = np.fromfile(file=f, count=3, dtype=np.float64)
            emit = np.fromfile(file=f, count=1, dtype=np.int32)
            sec = np.fromfile(file=f, count=1, dtype=np.int32)

        # Scale to nm
        x = x / 1.0E-9
        y = y / 1.0E-9
        z = z / 1.0E-9

        #r = np.sqrt(np.square(x) + np.square(y))
        particles1.set_data(x, y)
        particles2.set_data(x, z)
        #return particles2 #, circle1, circle2, rect
        return particles1, particles2 #, circle1, circle2, rect

    # Create the animation
    anim = animation.FuncAnimation(fig, animate, frames=max_steps, interval=1, blit=True, repeat=False, init_func=init)

    # Show the animation
    plt.show()

    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=100, metadata=dict(artist='Me')) #bitrate=1800)

    anim.save('vid.mp4', writer=writer, dpi=100)
