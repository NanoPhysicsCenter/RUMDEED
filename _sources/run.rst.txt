.. _run:

Input files and running the code
================================

Once you have :ref:`downloaded the source code and build the executable <build>` you can run program using the instructions below.

Input files
-----------
The RUMDEED program is executed from the command line and retrieves input parameters from various input files.
These files should be located in the same directory as the RUMDEED executable file. The primary input file, named **input**, is always necessary.
Additional input files are only required when utilizing specific functionalities of the program.

Input file
++++++++++
The program will read the file called **input** when it starts. This file defines the simulation parameters to be used and is a Fortran namelist file.
The parameters are described below, with an example given at the end.

V_S
    A number in units of Volts [V] that specifies the voltage from the source.
BOX_DIM(1:3)
    Dimensions of the simulations box given as three numbers in nano-meters [nm]. Should specify x, y and z. For now only the z value is used.
TIME_STEP
    The size of the timestep Δt in the simulation given in pico-seconds [ps].
STEPS
    Number of timesteps to do in the simulations. Should be an integer number large than zero.
EMISSION_MODE
    The emission mode controls the type of emission from the cathode. It should given as one of the following integer numbers:
    
    1: :ref:`Photo emission <photo>`.

    3: :ref:`Field emission from a hyperboloid tip <field-tip>`.

    9: :ref:`Thermal-field emission <thermal-field>`.

    10: :ref:`Planar field emission <field>`.
NREMIT
    Number of emitters in the system. Should be an integer larger than zero. Note that not all emission modes support multiple emitters.
IMAGE_CHARGE
    Boolean (.TRUE. / .FALSE.) if system should include image charge effects.
EMITTERS_DIM(1:3, <EMITTER NUMBER>)
    Dimensions of the emitter given as three numbers in nano-meters [nm].
    If the emitter is a circle the first number is interpreted as the radius of the emitter, 
    second is used for workfunction definitions and should be identical to the first number.
    For a rectangular emitter the first two numbers represent the length in the x and y directions respectively.
    Is the emitter is a hyperboloid tip the first number distance from the peak of the tip to the anode. The second number is the base radius of the tip
    and the last number is the height of the tip. See :ref:`Field emission from a hyperboloid tip <field-tip>` in the code description for more details.
    Here **<EMITTER NUMBER>** should be replaced by the number of the emitter.
    This line should be given for all emitters in the system.
EMITTERS_POS(1:3, <EMITTER NUMBER>)
    Position of the emitter given as three numbers in nano-meters [nm]. Here **<EMITTER NUMBER>** should be replaced by the number of the emitter.
    This line should be given for all emitters in the system.
    The position for circle and rectangle is measured from lower left corner of the rectangle (enclosing the circle).
EMITTERS_TYPE(<EMITTER NUMBER>)
    Geometry of the emitter is given as an integer. It should be given as one of the following integer numbers:
    
    1: Circular emitter. Only photo emission supports circular emitters.

    2: Rectangle.
    
    3: Rectangle spots.
    
    Here **<EMITTER NUMBER>** should be replaced by the number of the emitter. This line should be given for all emitters in the system.
EMITTERS_DELAY(<EMITTER NUMBER>)
    The timestep the emitter should become active and start emitting. Here **<EMITTER NUMBER>** should be replaced by the number of the emitter.
    This line should be given for all emitters in the system.

An example :download:`input file <files/input>` with one emitter doing planar field emission can be seen below:

.. code-block:: text

  &INPUT
    V_S = 1.0d3,
    BOX_DIM = 0.0d0, 0.0d0, 1000.0d0,
    TIME_STEP = 0.25d-3,
    STEPS = 20000,
    EMISSION_MODE = 10,
    NREMIT = 1,
    IMAGE_CHARGE = .TRUE.

    EMITTERS_DIM(1:3, 1) 1000.0d0, 1000.0d0, 0.0d0,
    EMITTERS_POS(1:3, 1) = -500.0d0, -500.0d0, 0.0d0,
    EMITTERS_TYPE(1) = 2,
    EMITTERS_DELAY(1) = 0,
  /

Work function
+++++++++++++
To specify the work function on the emitter surface for field emission and thermal-field emission, an input file called **work** is utilized.
For the hyperboloid tip surface, the work function needs to be specified within the file **mod_emission_top.f90** using the work variable.

The first line must contain an integer number. Currently, only the number 1 is supported, indicating a checkerboard work function surface.
The second line should consist of two integer numbers, representing the number of rows and columns in the checkerboard, respectively.
The remaining portion of the file should consist of a matrix of numbers that correspond to the work function values within the checkerboard.
It is important to note that the checkerboard functionality is only supported for a single emitter.

An example checkerboard :download:`work input function file <files/work>` can be seen below:

.. code-block:: text

  1
  2 2
  2.5d0 3.0d0
  3.0d0 2.5d0  

..
  Collisions
  ++++++++++
  N\ :sub:`2` files To do later when collisions are fully implemented

Laser
+++++

Laser file

.. code-block:: console

    1 2 2
    4.7 0.02
    10000 1000 5

.. Photoemission Input Warning
.. +++++++++++++++++++++++++++

The header/first line sets parameters;
The first number enables Gaussian electron emission pulse, 1 = on, 2 = off.
Second number selects type of laser input, 1 for fixed photon energy, 2 for Poisson distributed photon energy.
Third number picks velocity profile for electrons, 1 being zero initial velocity, 2 for work function dependant inital velocity.
    
Second line is laser (photon) energy and variation, first being the laser 'mean' energy level in electronVolts (eV) and second being standard deviation of the laser (in eV's as well). 
This is normal distribution with Box-Muller method.
For work function dependant initial velocity the energy is compared to the work function with the excess making way for Newtonian velocity given to the electrons.

Third line is gauss pulse parameters, center (mu), width (sigma) and A(mplitude) of the pulse. 
The gaussian pulse is simulated with output restriction of electrons according to normal distribution.
This should in theory simulate the Quantum Efficiency and Intensity via amplitude modulation.

Examples
--------

A few examples of input files can be found in the Examples folder.

Space charge limited photo emission
+++++++++++++++++++++++++++++++++++

Photoemission

Planar field emission
+++++++++++++++++++++

The examples for planar field emission can be found in Examples/Planar-FE/ folder. Examples for homogeneous work function of 2.0 eV is given along with an alternating 2.0 and 2.5 eV
checkerboard work function pattern.

Input file

.. literalinclude:: files/Examples/Planar-FE/input

Work function file

.. literalinclude:: files/Examples/Planar-FE/work

Field emission from a prolate spheroidal tip
++++++++++++++++++++++++++++++++++++++++++++

Tip

Planar thermal-field emission
+++++++++++++++++++++++++++++

Thermal-field


Running the code
-----------------

To run the code place the necessary input files into the same folder as the executable file :ref:`"RUMDEED.out" <build>`. The execute the code by running the following command
inside the folder

.. code-block:: console
   
   ./RUMDEED.out

Output files will be placed in a folder called out and are described in :ref:`output`.

.. index:: Collisions, N₂, input, work, work function, time_step, time step, box_dim, steps, emission_mode, nremit