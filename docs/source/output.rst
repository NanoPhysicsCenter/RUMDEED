.. _output:

=====================
Output and data files
=====================

Output files after running the code can be found in the out folder. Files either have an ending in .dt (**D**\ ata **T**\ ext) or .bin (**Bin**\ ary).
Text file ending in .dt have to data orginzed as columns. Files ending in .bin are binary files. Both files types are decsribed below.

----------

| **init.dt**
| The file init.dt contains the parameters used in the program. It is a Fortran namelist file were importand parameters are shown.

----------

| **emitted.dt**
| The file emitted.dt contains information about the number of particles emitted from the cathode.
| The date in the columns is

.. list-table:: emitted.dt
   :widths: auto
   :header-rows: 1
   :stub-columns: 1

   * - 
     - Time
     - Time step
     - Number of emitted particles
     - Total number of electron in the system
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

| **absorbed.dt**, **absorb_bot.dt** and **absorb_top.dt**
| The file absorbed.dt contains information about the number of particles absorbed by the cathode or anode. While absorb_bot.dt contains data only for the cathode and absorb_top.dt for the anode.
| The date in the columns is

.. list-table:: absorbed.dt, absorb_bot.dt, absorb_top.dt
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * -
    - Time
    - Time step
    - Number of absorbed particles
    - Number of electrons absorbed
    - Number of holes absorbed
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

----------

| **field.dt**
| The file field.dt contains information about the electric field.
| The date in the columns is

.. list-table:: field.dt
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * -
    - Time step
    - Average of the x compontent of the Electric field on the cathode
    - Average of the y compontent of the Electric field on the cathode
    - Average of the z compontent of the Electric field on the cathode
    - Unrounded results of the MC integration of how many electrons to emitt from the cathode
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

---------

| **integration.dt** 
| The file integration.dt contains information about the Monte Carlo integration performed by the Cuba library :cite:p:`HAHN200578`. The default integration method used is call Divone. Please refer to the Cuba library manual for more information. It can be found in the Cuba directory.
| The date in the columns is

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

---------

| **ramo_current.dt**
| The file ramo_current.dt contains information about the Ramo current.
| The date in the columns is

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
    - Number of holes in the system
    - Average mobility of the particles in the system
    - Average speed of the particles in the system
    - Ramo current for species 1
    - ...
    - Ramo current for species N
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
    - A
    - ...
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
    - ...
    - float

---------

| **volt.dt**
| The file volt.dt contains information about the voltage over the gap. In future releases it will be possible to have voltage that depends on time or other factors.
| The date in the columns is

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

---------

| **gauss.dt**
| The file gauss.dt contains information about the Gaussian distribution used in the photoemission model. The Gaussian emission can be turned on and off in the mod_photo_emission.f90 file by setting the variable **EmitGauss** to true or false. Parameters for the Gaussian emission can be set in the :ref:`laser file <run>`.
| The date in the columns is

.. list-table:: gauss.dt
  :widths: auto
  :header-rows: 1
  :stub-columns: 1

  * -
    - Time step
    - Number of electrons allowed by the gaussian distribution
  * - Units
    - #
    - #
  * - Type
    - integer
    - integer

---------

| **density_emit.bin**
| The file density_emit.bin contains information about the density of the emitted particles. The file is written in binary format and contains the x and y coordinates of the emitted particles, along with emitted id, section and the ID number of the particle. This data is repated for each emitted particle.
| An example of how to read the data is shown in the Python Jupiter notebooks found in the exmaples folder.

.. list-table:: density_emit.bin
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

---------

| **density_absorb_top.bin** and **density_absorb_bot.bin**
| The files density_absorb_top.bin and density_absorb_bot.bin contains information about the density of the absorbed particles. The files are written in binary format and contain the x and y coordinates of the absorbed particles, along with absorbed id, section and the ID number of the particle. This data is repated for each absorbed particle. The top file also contains the velocity of velocity of the absorbed particles.
| An example of how to read the data is shown in the Python Jupiter notebooks found in the exmaples folder.

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

---------

| **planes-?.bin**
| The files planes-?.bin contain information about particles when they pass through imaginary planes in the system. The number of planes-?.dt files is set using the **planes_N** variable in the input file. The positions of the imagninary planes are the specified using the **planes_z** variable. Up to 10 planes can be specified. The planes are numbered from 0 to 9. The planes are specified in the input file as

.. code-block:: text
  
    planes_N = 2
    planes_z = 0.0d0, 100.0d0

| The files are written in binary format and contain the x and y coordinates of the particles, the velocity of the particles, along with the emitter id, emitter section and particle id. This data is repated for each particle that passes through the plane.

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
