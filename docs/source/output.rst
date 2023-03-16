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
| The file emitted.dt contains information about the number of particles emitted from the cathode. The date in is column is

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

Scripts
-------

Describe some of the scripts in the scripts/ folder. How to use them and what they do.