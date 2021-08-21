Input files and running the code
================================

To do run the code

Input files
-----------
input, w_theta, collisions

Input file
++++++++++
The program will read a file called **input** when it starts. This file defines the simulation parameters to be used and is a Fortran namelist file.
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
EMITTERS_DIM(1:3, <EMITTER NUMBER>)
    Dimensions of the emitter given as three numbers in nano-meters [nm].
    If the emitter is a circle only the first number is relevant and is interpreted as the radius of the emitter.
    For a rectangular emitter the first two numbers represent the length in the x and y directions respectively.
    Is the emitter is a hyperboloid tip the first number distance from the peak of the tip to the anode. The second number is the base radius of the tip
    and the last number is the height of the tip. See :ref:`Field emission from a hyperboloid tip <field-tip>` in the code description for more details.
    Here **<EMITTER NUMBER>** should be replaced by the number of the emitter.
    This line should be given for all emitters in the system.
EMITTERS_POS(1:3, <EMITTER NUMBER>)
    Position of the emitter given as three numbers in nano-meters [nm]. Here **<EMITTER NUMBER>** should be replaced by the number of the emitter.
    This line should be given for all emitters in the system.
EMITTERS_TYPE(<EMITTER NUMBER>)
    Geometry of the emitter given as an integer. It should be given as one of the following integer numbers:
    
    1: Circular emitter. The dimensions **to do**

    2: Rectangle. The dimensions ... **to do**
    
    Here **<EMITTER NUMBER>** should be replaced by the number of the emitter. This line should be given for all emitters in the system.
EMITTERS_DELAY(<EMITTER NUMBER>)
    The timestep the emitter should become active and start emitting. Here **<EMITTER NUMBER>** should be replaced by the number of the emitter.
    This line should be given for all emitters in the system.

Example input file with one emitter doing planar field emission:

.. code-block:: text

  &INPUT
    V_S = 1.0d3,
    BOX_DIM = 0.0d0, 0.0d0, 1000.0d0,
    TIME_STEP = 0.25d-3,
    STEPS = 20000,
    EMISSION_MODE = 10,
    NREMIT = 1,

    EMITTERS_DIM(1:3, 1) 1000.0d0, 1000.0d0, 1000.0d0,
    EMITTERS_POS(1:3, 1) = -500.0d0, -500.0d0, 0.0d0,
    EMITTERS_TYPE(1) = 1,
    EMITTERS_DELAY(1) = 0,
  /

Work function
+++++++++++++
work function

Example work function file::

  1 1
  1
  2.5d0

Collisions
++++++++++
N\ :sub:`2` files


Running the code
-----------------

To do

.. index:: Collisions, N₂, input, w_theta, work function, time_step, time step, box_dim, steps, emission_mode, nremit