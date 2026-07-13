.. _run:

Input files and running the code
================================

Once you have :ref:`downloaded the source code and built the executable <build>` you can run the program using the instructions below.

Input files
-----------
The RUMDEED program is executed from the command line and retrieves input parameters from various input files.
These files should be located in the same directory as the RUMDEED executable file. The primary input file, named **input**, is always necessary.
Additional input files are only required when utilizing specific functionalities of the program.

Input file
++++++++++
The program will read the file called **input** when it starts. This file defines the simulation parameters to be used and is a Fortran namelist file.
The parameters are described below, with an example given at the end.

.. rubric:: General setup

V_S
    A number in units of Volts [V] that specifies the voltage from the source.
BOX_DIM(1:3)
    Dimensions of the simulation box given as three numbers in nanometers [nm]. Should specify x, y and z. For now only the z value is used.
TIME_STEP
    The size of the time step Δt in the simulation given in picoseconds [ps].
STEPS
    Number of time steps to do in the simulation. Should be an integer number larger than zero.
EMISSION_MODE
    The emission mode controls the type of emission from the cathode. It should be given as one of the following integer numbers:

    0: Unit and integration tests.

    1: :ref:`Photo emission <photo>`.

    2: Planar field emission (older model, superseded by mode 10).

    3: :ref:`Field emission from a hyperboloid tip <field-tip>`.

    4: Thermionic emission (currently disabled in the code).

    5, 6, 7, 8: Field emission from a 2D material (variants of the 2DEG and Dirac electron-supply models, see **mod_field_emission_2D.f90**).

    9: :ref:`Thermal-field emission <thermal-field>`.

    10: :ref:`Planar field emission <field>`.

    11: Field emission from a cylindrical tip.

    12: Field emission from a torus.

    999: Manual placement of electrons for testing and debugging (see **mod_manual_emission.f90**).
NREMIT
    Number of emitters in the system. Should be an integer larger than zero. Note that not all emission modes support multiple emitters.
    The code is compiled with room for a single emitter (**MAX_EMITTERS** in mod_global.F90); to use more than one emitter that
    constant must be increased and the code recompiled.

.. rubric:: Emitters

EMITTERS_DIM(1:3, <EMITTER NUMBER>)
    Dimensions of the emitter given as three numbers in nanometers [nm].
    If the emitter is a circle, the first number is interpreted as the radius of the emitter,
    the second is used for work function definitions and should be identical to the first number.
    For a rectangular emitter the first two numbers represent the length in the x and y directions respectively.
    If the emitter is a hyperboloid tip, the first number is the distance from the peak of the tip to the anode. The second number is the base radius of the tip
    and the last number is the height of the tip. See :ref:`Field emission from a hyperboloid tip <field-tip>` in the code description for more details.
    Here **<EMITTER NUMBER>** should be replaced by the number of the emitter.
    This line should be given for all emitters in the system.
EMITTERS_POS(1:3, <EMITTER NUMBER>)
    Position of the emitter given as three numbers in nanometers [nm]. Here **<EMITTER NUMBER>** should be replaced by the number of the emitter.
    This line should be given for all emitters in the system.
    The position for circle and rectangle is measured from the lower left corner of the rectangle (enclosing the circle).
EMITTERS_TYPE(<EMITTER NUMBER>)
    Geometry of the emitter is given as an integer. It should be given as one of the following integer numbers:

    1: Circular emitter. Only photo emission supports circular emitters.

    2: Rectangle.

    3: Rectangle spots.

    4: Ring.

    Here **<EMITTER NUMBER>** should be replaced by the number of the emitter. This line should be given for all emitters in the system.
EMITTERS_DELAY(<EMITTER NUMBER>)
    The time step at which the emitter becomes active and starts emitting. Here **<EMITTER NUMBER>** should be replaced by the number of the emitter.
    This line should be given for all emitters in the system.

.. rubric:: Image charge

IMAGE_CHARGE
    Boolean (``.TRUE.`` / ``.FALSE.``) if the system should include image charge effects.
N_IC_MAX
    Number of image charge partners to use when **IMAGE_CHARGE** is ``.TRUE.``. The default 0 uses a single image charge partner,
    while 1 uses 5 partners. See the function Force_Image_charges_v2 in mod_verlet.F90 for details.

.. rubric:: Emission sampling and integration

MH_BATCH
    Boolean, defaults to ``.FALSE.``. Only used in planar field emission (mode 10). When set to
    ``.TRUE.`` all Metropolis-Hastings chains of a time step advance in lockstep and the surface
    field is evaluated for all of them in batches, which is much faster. This is also required for the
    emission to benefit from OpenACC GPU offload, so set it to ``.TRUE.`` in GPU runs. All chains
    then sample the particle configuration from the start of the time step, whereas with the
    default ``.FALSE.`` each chain runs to completion in turn and sees the electrons emitted
    before it within the same step. The two settings agree statistically, but not run for
    run; the default ``.FALSE.`` keeps the serial sampling behaviour (and its results) of
    versions that predate this option.
CUBA_METHOD
    Integration method the Cuba library uses for the emission-current integrals: 1 for Suave, 2 for Divonne (the default) or 3 for Cuhre.
CUBA_EPSABS
    Requested absolute error of the integration. The integrals are in units of electrons per time step, so the default 0.5 means
    half an electron.
CUBA_EPSREL
    Requested relative error of the integration. The default is 1.0d-3. The integration stops when either the absolute or the
    relative error is reached.
CUBA_MINEVAL
    Minimum number of integrand evaluations, default 1000.
CUBA_MAXEVAL
    Maximum number of integrand evaluations, default 5000000.

.. rubric:: Collisions and background gas

COLLISION_MODE
    Electron collisions with neutral N\ :sub:`2` gas. It should be given as one of the following integer numbers (see **mod_collisions.F90**):

    0: No collisions (default).

    1: Continuous ionization.

    2: Continuous ionization and discrete recombination.

    3: Discrete ionization.

    4: Discrete ionization and discrete recombination.
COLLISION_DELAY
    The time step at which collisions become active. The default is 0.
ION_ATOM_RATIO
    Ratio of ions to neutral atoms when the background gas is initialized. The default is 0.
ION_LIFE_TIME
    Lifetime of ions in time steps before they are removed from the system. The default is 100000000.
T_TEMP
    Temperature of the background gas in Kelvin. The default is 293.15 K.
P_ABS
    Pressure of the background gas as a fraction of normal pressure. The default is 1.
ATOM_TIME_INTERVAL
    If given a value larger than zero, neutral atoms are advanced on a coarser time step of ATOM_TIME_INTERVAL × TIME_STEP
    (the two-time-step scheme used with collisions). Set it to 0 to advance all particles every time step. This parameter should
    always be set explicitly when running with collisions, as it has no default value.

.. rubric:: External circuit

R_S, R_P, L_P, C_P
    Series resistance [Ω] and parallel resistance [Ω], inductance [H] and capacitance [F] of an external circuit model.
    These parameters are read from the input file but the circuit model in mod_verlet.F90 is currently disabled, so they have no effect.

.. rubric:: Additional output files

PLANES_N
    Number of imaginary recording planes in the system, an integer from 0 to 10. Particles passing through a plane are recorded in
    the :ref:`planes-?.bin <output>` files. The default is 10 planes at z = 5, 10, 25, 50, 75, 100, 125, 250, 500 and 750 nm.
PLANES_Z(1:PLANES_N)
    The z positions of the recording planes in nanometers [nm]. Only planes with a position larger than zero produce an output file.
WRITE_RAMO_SEC
    Boolean, defaults to ``.FALSE.``. Write the Ramo current broken down into emitter sections to the file ramo_current.bin every time step.
WRITE_POSITION_FILE
    Boolean, defaults to ``.FALSE.``. Write the position of every particle in the system to the file position.bin every time step.
WRITE_PARTICLE_DATA_FILE
    Boolean, defaults to ``.FALSE.``. Write a text file out/particles-<step>.dt every time step with the ID, position, emitter and
    creation step of every particle.
WRITE_ELECTRON_DATA_FILE, WRITE_ION_DATA_FILE
    Booleans, default to ``.FALSE.``. Like WRITE_PARTICLE_DATA_FILE but write only electrons (out/electrons-<step>.dt) or only
    ions (out/ions-<step>.dt). They are only used when WRITE_PARTICLE_DATA_FILE is ``.FALSE.``.
WRITE_RECOMBINATION_FILE
    Boolean, defaults to ``.FALSE.``. Read from the input file but currently unused.
SAMPLE_ATOM_FILE, SAMPLE_ATOM_RATE
    If SAMPLE_ATOM_FILE is ``.TRUE.`` the positions of the neutral atoms are written to a binary file out/atom-<step>.bin every
    SAMPLE_ATOM_RATE time steps (default 500).
SAMPLE_ELEC_FILE, SAMPLE_ELEC_RATE
    If SAMPLE_ELEC_FILE is ``.TRUE.`` the positions of the electrons are written to a binary file out/elec-<step>.bin every
    SAMPLE_ELEC_RATE time steps (default 500).

An example :download:`input file <files/input>` with one emitter doing planar field emission can be seen below:

.. code-block:: text

  &INPUT
    V_S = 1.0d3,
    BOX_DIM = 0.0d0, 0.0d0, 1000.0d0,
    TIME_STEP = 0.25d-3,
    STEPS = 20000,
    EMISSION_MODE = 10,
    NREMIT = 1,
    IMAGE_CHARGE = .TRUE.,

    EMITTERS_DIM(1:3, 1) = 1000.0d0, 1000.0d0, 0.0d0,
    EMITTERS_POS(1:3, 1) = -500.0d0, -500.0d0, 0.0d0,
    EMITTERS_TYPE(1) = 2,
    EMITTERS_DELAY(1) = 0,
  /

Work function
+++++++++++++
To specify the work function on the emitter surface for field emission and thermal-field emission, an input file called **work** is utilized.
For the hyperboloid tip surface, the work function needs to be specified within the file **mod_emission_tip.f90** using the ``w_theta`` variable.

The first line must contain an integer number that selects the work function model (see **mod_work_function.F90**):

1: Checkerboard.

2: Gaussian spots.

4: Voronoi regions.

The number 3 is reserved for a circle model that is not implemented.

For the **checkerboard** model, the second line should consist of two integer numbers, representing the number of rows and columns
in the checkerboard, respectively. The remaining portion of the file should consist of a matrix of numbers that correspond to the
work function values in eV within the checkerboard.

For the **Gaussian spots** model, the second line gives the base work function in eV and the third line the number of Gaussian spots.
Each following line then describes one spot with five numbers: the amplitude added to the base value in eV, the center coordinates
x and y in nm, and the standard deviations in the x and y directions in nm.

For the **Voronoi** model, the second line gives the number of Voronoi sites. Each following line then describes one site with four
numbers: the x and y coordinates of the site on a scale from 0 to 1 over the emitter area, the work function value in eV, and the
section number of the site.

It is important to note that these work function models are only supported for a single emitter.

An example checkerboard :download:`work function input file <files/work>` can be seen below:

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

The photo emission mode reads the laser pulse parameters from an input file called **laser**. An example is shown below:

.. code-block:: console

    1 2 2
    4.7 0.02
    10000 1000 5

The first line sets the following parameters:

- The first number enables the Gaussian electron emission pulse: 1 = on, 2 = off.
- The second number selects the type of laser input: 1 for fixed photon energy, 2 for Poisson distributed photon energy.
- The third number picks the velocity profile for electrons: 1 for zero initial velocity, 2 for work function-dependent initial velocity.

The second line specifies the laser (photon) energy and its spread. The first value is the laser mean photon energy in electronvolts (eV) and the second is the standard deviation (also in eV).
This uses a normal distribution with the Box-Muller method.
For work function-dependent initial velocity, the photon energy is compared to the work function, with the excess converted to Newtonian velocity given to the electrons.

The third line contains the Gaussian pulse parameters: center (mu), width (sigma), and amplitude (A) of the pulse.
The Gaussian pulse is simulated by restricting electron output according to a normal distribution.
This controls the quantum efficiency and intensity via amplitude modulation.

Examples
--------

A few examples of input files can be found in the Examples folder.

Space charge limited photo emission
+++++++++++++++++++++++++++++++++++

Photoemission

Planar field emission
+++++++++++++++++++++

The examples for planar field emission can be found in the Examples/Planar-FE folder. An example for a homogeneous work function of 2.0 eV is given along with an alternating 2.0 and 2.5 eV
checkerboard work function pattern.

Input file

.. literalinclude:: files/Examples/Planar-FE/input

Work function file

.. literalinclude:: files/Examples/Planar-FE/work

Field emission from a hyperboloid tip
+++++++++++++++++++++++++++++++++++++

Tip

Planar thermal-field emission
+++++++++++++++++++++++++++++

Thermal-field


Running the code
-----------------

To run the code, place the necessary input files into the same folder as the executable file :ref:`RUMDEED.out <build>`. Then execute the code by running the following command
inside the folder:

.. code-block:: console

   ./RUMDEED.out

Output files will be placed in a folder called out and are described in :ref:`output`.

.. index:: Collisions, N₂, input, work, work function, time_step, time step, box_dim, steps, emission_mode, nremit