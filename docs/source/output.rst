.. _output:

=====================
Output and data files
=====================

Output files after running the code can be found in the out folder. Files either have an ending in .dt (**D**\ ata **T**\ ext) or .bin (**Bin**\ ary).
Text files ending in .dt have data organized as columns. Files ending in .bin are binary files. Both file types are described below.

----------

| **init.dt**
| The file init.dt contains the parameters used in the program. It is a Fortran namelist file where important parameters are shown.

----------

| **init.bin**
| The file init.bin contains system parameters and scales in binary format, for use by analysis scripts. It consists of a single record with the layout below. See the data analysis examples section for how to read it in Python.

.. list-table:: init.bin
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * - Data
    - type
    - Units
  * - epsilon_r, m_eeff, m_heff
    - 3 × float64
    - #
  * - length_scale
    - float64
    - m
  * - time_scale
    - float64
    - s
  * - vel_scale
    - float64
    - m/s
  * - cur_scale
    - float64
    - A
  * - MAX_PARTICLES, MAX_EMITTERS, MAX_SECTIONS, MAX_LIFE_TIME
    - 4 × int32
    - #

----------

| **emitted.dt**
| The file emitted.dt contains information about the number of particles emitted from the cathode.
| The data in the columns is

.. list-table:: emitted.dt
   :widths: auto
   :header-rows: 1
   :stub-columns: 1

   * - 
     - Time
     - Time step
     - Number of emitted particles
     - Total number of electrons in the system
     - Number of emitted particles from emitter number 1
     - Number of emitted particles from emitter number 2
     - ...
     - Number of emitted particles from emitter number N
   * - Units
     - ps
     - #
     - #
     - #
     - #
     - #
     - ...
     - #
   * - Type
     - float
     - integer
     - integer
     - integer
     - integer
     - integer
     - ...
     - integer

----------

| **absorbed.dt**, **absorbed_bot.dt**, **absorbed_top.dt** and **absorbed_recom.dt**
| The file absorbed.dt contains information about the number of particles absorbed by the cathode or anode, while absorbed_bot.dt contains data only for the cathode and absorbed_top.dt only for the anode. The file absorbed_recom.dt contains the number of particles removed through recombination.
| The data in the columns is

.. list-table:: absorbed.dt, absorbed_bot.dt, absorbed_top.dt, absorbed_recom.dt
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * -
    - Time
    - Time step
    - Number of absorbed particles
    - Number of electrons absorbed
    - Number of ions absorbed
  * - Units
    - ps
    - #
    - #
    - #
    - #
  * - Type
    - float
    - integer
    - integer
    - integer
    - integer

| The file absorbed_top.dt contains additional columns after the ones above: the number of electrons absorbed from emitter number 1 through N (integers).

----------

| **field.dt**
| The file field.dt contains information about the electric field.
| The data in the columns is

.. list-table:: field.dt
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * -
    - Time step
    - Average of the x component of the Electric field on the cathode
    - Average of the y component of the Electric field on the cathode
    - Average of the z component of the Electric field on the cathode
    - Unrounded results of the integration of how many electrons to emit from the cathode
    - Debug information for the MH algorithm
    - Debug information for the MH algorithm
    - Debug information for the MH algorithm
  * - Units
    - #
    - V/m
    - V/m
    - V/m
    - #
    - #
    - #
    - #
  * - Type
    - integer
    - float
    - float
    - float
    - float
    - float
    - float
    - float

----------

| **integration.dt**
| The file integration.dt contains information about the numerical integration performed by the Cuba library :cite:p:`HAHN200578`. The default integration method is Divonne; it can be changed with the **cuba_method** parameter in the input file. Please refer to the Cuba library manual for more information. It can be found in the Cuba directory.
| The data in the columns is

.. list-table:: integration.dt
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * -
    - Number of the emitter
    - Number of regions used in the integration
    - Number of evaluations
    - Fail flag
    - Integral results
    - Error
    - Probability
  * - Units
    - #
    - #
    - #
    - #
    - #
    - #
    - #
  * - Type
    - integer
    - integer
    - integer
    - integer
    - float
    - float
    - float

| Note that the cylindrical tip and torus emission modes omit the first column (the number of the emitter).

----------

| **ramo_current.dt**
| The file ramo_current.dt contains information about the Ramo current.
| The data in the columns is

.. list-table:: ramo_current.dt
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * -
    - Time
    - Time step
    - Ramo current
    - Voltage over the gap
    - Number of particles in the system
    - Number of electrons in the system
    - Number of ions in the system
    - Average mobility of the particles in the system
    - Average speed of the particles in the system
    - Average speed of the electrons in the system
    - Average speed of the ions in the system
    - Ramo current for electrons
    - Ramo current for ions
    - Ramo current for atoms
  * - Units
    - ps
    - #
    - A
    - V
    - #
    - #
    - #
    - m²/Vs
    - m/s
    - m/s
    - m/s
    - A
    - A
    - A
  * - Type
    - float
    - integer
    - float
    - float
    - integer
    - integer
    - integer
    - float
    - float
    - float
    - float
    - float
    - float
    - float

----------

| **volt.dt**
| The file volt.dt contains information about the voltage over the gap. The voltage is currently constant during the simulation and the time dependent column is always zero.
| The data in the columns is

.. list-table:: volt.dt
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * -
    - Time
    - Time step
    - Voltage over the gap
    - Time dependent voltage (Currently set to zero)
  * - Units
    - ps
    - #
    - V
    - V
  * - Type
    - float
    - integer
    - float
    - float

----------

| **gauss.dt**
| The file gauss.dt contains information about the Gaussian distribution used in the photoemission model. The Gaussian emission is turned on and off with the first number in the :ref:`laser file <run>`, where its parameters are also set.
| The data in the columns is

.. list-table:: gauss.dt
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * -
    - Time step
    - Number of electrons allowed by the Gaussian distribution
  * - Units
    - #
    - #
  * - Type
    - integer
    - integer

----------

| **collisions.dt**
| The file collisions.dt is written every time step when collisions are enabled (**collision_mode** > 0 in the input file).
| The data in the columns is

.. list-table:: collisions.dt
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * -
    - Time step
    - Number of collisions
    - Number of ionizations
    - Number of recombinations
  * - Units
    - #
    - #
    - #
    - #
  * - Type
    - integer
    - integer
    - integer
    - integer

----------

| **density_emit_elec.bin**, **density_emit_ion.bin** and **density_emit_atom.bin**
| These files contain the positions where particles were added to the system, split by species: electrons in density_emit_elec.bin, ions in density_emit_ion.bin and neutral atoms in density_emit_atom.bin. The files are written in binary format and the record below is repeated for each particle.
| Note that the files density_emit.bin and density_ion.bin are still created but are no longer written to; the data now goes to the per-species files described here.
| An example of how to read this kind of data is shown in the Python Jupyter notebooks found in the examples folder.

.. list-table:: density_emit_elec.bin, density_emit_ion.bin, density_emit_atom.bin
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * - Data
    - type
    - Units
  * - x position
    - float64
    - nm
  * - y position
    - float64
    - nm
  * - z position
    - float64
    - nm
  * - Emitter ID
    - int32
    - #
  * - Particle ID
    - int32
    - #

----------

| **density_absorb_top.bin** and **density_absorb_bot.bin**
| The files density_absorb_top.bin and density_absorb_bot.bin contains information about the density of the absorbed particles. The files are written in binary format and contain the x and y coordinates of the absorbed particles, along with absorbed id, section and the ID number of the particle. This data is repeated for each absorbed particle. The top file also contains the velocity of the absorbed particles.
| An example of how to read the data is shown in the Python Jupyter notebooks found in the examples folder.

.. list-table:: density_absorb_top.bin
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * - Data
    - type
    - Units
  * - x position
    - float64
    - nm
  * - y position
    - float64
    - nm
  * - x velocity
    - float64
    - m/s
  * - y velocity
    - float64
    - m/s
  * - z velocity
    - float64
    - m/s
  * - Emitter ID
    - int32
    - #
  * - Emitter section
    - int32
    - #
  * - Particle ID
    - int32
    - #

.. list-table:: density_absorb_bot.bin
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * - Data
    - type
    - Units
  * - x position
    - float64
    - nm
  * - y position
    - float64
    - nm
  * - Emitter ID
    - int32
    - #
  * - Emitter section
    - int32
    - #
  * - Particle ID
    - int32
    - #

----------

| **ionization_data.bin**
| The file ionization_data.bin contains one binary record for each ionization event.

.. list-table:: ionization_data.bin
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * - Data
    - type
  * - Time step
    - int32
  * - Ion position (x, y, z)
    - 3 × float64
  * - Speed of the incoming electron
    - float64
  * - Speed of the incoming electron after the collision
    - float64
  * - Speed of the new electron
    - float64
  * - Distance between electron and atom
    - float64
  * - Ionization cross section radius
    - float64
  * - ID of the incoming electron
    - int32
  * - ID of the new electron
    - int32
  * - ID of the ion
    - int32
  * - Emitter of the incoming electron
    - int32

----------

| **recombination_data.bin**
| The file recombination_data.bin contains one binary record for each recombination event.

.. list-table:: recombination_data.bin
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * - Data
    - type
  * - Time step
    - int32
  * - Ion position (x, y, z)
    - 3 × float64
  * - Speed of the electron
    - float64
  * - Distance between electron and ion
    - float64
  * - Kramers cross section radius
    - float64
  * - ID of the electron
    - int32
  * - ID of the ion
    - int32
  * - Emitter of the electron
    - int32
  * - Lifetime of the ion in time steps
    - int32

----------

| **position.bin**
| The file position.bin is only written when **write_position_file** is set to ``.TRUE.`` in the input file. It starts with a header of two int32 numbers: the total number of time steps and the sampling interval. Then for every time step a block follows consisting of the time step number and the number of particles (2 × int32), followed by one record per particle with the position x, y, z (3 × float64, in meters), the emitter number, the emitter section and the particle ID (3 × int32).

----------

| **ramo_current.bin**
| The file ramo_current.bin is only written when **write_ramo_sec** is set to ``.TRUE.`` in the input file. Every time step the full array of Ramo currents per emitter section is written as MAX_SECTIONS × MAX_EMITTERS float64 numbers in amperes (column-major order, sections varying fastest).

----------

| **particles-<step>.dt**, **electrons-<step>.dt** and **ions-<step>.dt**
| These text files are written every time step when **write_particle_data_file**, **write_electron_data_file** or **write_ion_data_file** is set to ``.TRUE.`` in the input file. The first line contains the number of particles. Each following line describes one particle: ID, position x, y, z, emitter number and the time step the particle was created. In particles-<step>.dt the positions are in nm, in the other two files they are in meters.

----------

| **atom-<step>.bin** and **elec-<step>.bin**
| Snapshot files of atom and electron positions written every **sample_atom_rate** / **sample_elec_rate** time steps when **sample_atom_file** / **sample_elec_file** is set to ``.TRUE.`` in the input file.

----------

| **density_emit.bin**, **density_ion.bin**, **mh.bin**, **laplace_average_field.dt** and **laplace_grid.dt**
| These files are created in the out folder but the current version of the code does not write any data to them.

----------

| **density_absorb_recom.bin**
| The file density_absorb_recom.bin contains information about particles removed through recombination. The file is written in binary format and the record below is repeated for each removed particle. Note that unlike the other density files, the positions in this file are given in meters.

.. list-table:: density_absorb_recom.bin
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * - Data
    - type
    - Units
  * - x position
    - float64
    - m
  * - y position
    - float64
    - m
  * - z position
    - float64
    - m
  * - Emitter ID
    - int32
    - #
  * - Emitter section
    - int32
    - #
  * - Particle ID
    - int32
    - #
  * - Species (1 = electron, 2 = ion)
    - int32
    - #
  * - Time step
    - float64
    - #

----------

| **planes-?.bin**
| The files planes-?.bin contain information about particles when they pass through imaginary planes in the system. The number of planes-?.bin files is set using the **planes_N** variable in the input file. The positions of the imaginary planes are specified using the **planes_z** variable. Up to 10 planes can be specified. The files are numbered from 1 to planes_N, and only planes with a position larger than zero produce a file. If the input file does not set these variables, the default is 10 planes at z = 5, 10, 25, 50, 75, 100, 125, 250, 500 and 750 nm. The planes are specified in the input file as

.. code-block:: text
  
    planes_N = 2
    planes_z = 50.0d0, 100.0d0

| The files are written in binary format and contain the x and y coordinates of the particles, the velocity of the particles, along with the emitter id, emitter section and particle id. This data is repeated for each particle that passes through the plane.

.. list-table:: planes-?.bin
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * - Data
    - type
    - Units
  * - x position
    - float64
    - nm
  * - y position
    - float64
    - nm
  * - x velocity
    - float64
    - m/s
  * - y velocity
    - float64
    - m/s
  * - z velocity
    - float64
    - m/s
  * - Emitter ID
    - int32
    - #
  * - Emitter section
    - int32
    - #
  * - Particle ID
    - int32
    - #
