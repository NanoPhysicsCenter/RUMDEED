Building the code
=================

If you do not already have to code you can clone it from the NanoPhysics Centers Github repository.
You will need to have **git** installed and then you can clone the repository using the command:

.. code-block:: console

   git clone -b master https://github.com/NanoPhysicsCenter/Vacuum-MD.git Vacuum-MD

Where the last "Vacuum-MD" is the name of the directory you want to create and store the code in.
Once you have to code you will need the GNU Fortran compiler (**gfortran**) and **make** to compile the it.
Go into the directory "Vacuum-MD" and enter the command:

.. code-block:: console

   make

This will compile the code and an executable named "Vacuum-MD.out" will be created inside the build directory.

.. index:: build, make, git, compiler, gfortran, github
