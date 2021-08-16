Input files and running the code
================================

To do run the code

Input files
-----------
input, w_theta, collisions

Input file
++++++++++
The program will read a file called **input** when it starts. This file defines the simulation parameters to be used. The parameters are described
in the table below.

===================  =====  ===========
Variable             Units  Description
===================  =====  ===========
V_S                  V      Voltage from the source.
BOX_DIM(1:3)         nm     Dimensions of the simulations box. Should specify x, y and z. For now only the z value is used.
TIME_STEP            ps     The size of the timestep Δt in the simulation.
STEPS                #      Number of timesteps to do in the simulations.
EMISSION_MODE        #      The emission mode. Can take on the values of 1 for photo emission, 3 for field emission from a hyperboloid tip, 10 for planar field emission and 10 for thermal-field emission. 
NREMIT               #      Number of emitters in the system.
EMITTERS_DIM(1:3,:)  nm     Dimensions of the emitter.
EMITTERS_POS(1:3,:)  nm     Position of the emitter.
EMITTERS_TYPE(:)     #      Geometry of the emitter.
EMITTERS_DELAY(:)    #      The timestep the emitter should become active and start emitting.
===================  =====  ===========

.. bibliographic fields (which also require a transform):

:V_s: Voltage from the source.
:BOX_DIM(1:3): Dimensions of the simulations box. Should specify x, y and z. For now only the z value is used.

Example input file:

.. code-block:: text

  &INPUT
    V_S = 1.0d3,
    BOX_DIM = 0.0d0, 0.0d0, 1000.0d0,
    TIME_STEP = 0.25d-3,
    STEPS = 20000,
    EMISSION_MODE = 1,
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

Definition Lists
----------------

Term
    Definition
Term : classifier
    Definition paragraph 1.

    Definition paragraph 2.
Term
    Definition

I have no clue why the definition list below is classified as a different style
of definition list than the one above.

Is it the spaces in the term?
    Maybe it was the multiple line paragraph
    in the line below that caused this?

Is it the paragraph above the list maybe?
    I guess a lot of these lists don't have leading paragraphs?

Is it everything all at once?
    Who knows?!

.. index:: Collisions, N₂, input, w_theta, work function, time_step, time step, box_dim, steps, emission_mode, nremit