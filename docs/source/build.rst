Building the code
=================

.. _installation:

Installation
------------

If you do not already have to code you can clone it from the NanoPhysics Centers Github repository.
You will need to have have **git** installed and then to clone the repository use the command:

.. code-block:: console

   git clone -b master https://github.com/NanoPhysicsCenter/Vacuum-MD.git Vacuum-MD

Where last the "Vacuum-MD" is the name of the directory you want to create and store the code in.
Once you have to code you will need **gfortran** and **make** to compile the code. Go into the directory "Vacuum-MD" and enter the command:

.. code-block:: console

   make

This will then compile the code and an executable named "Vacuum-MD" will be created inside the build directory.