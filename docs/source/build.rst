.. _build:

Building the code
=================

The code was built to run on Linux systems. For Windows users it's recommended to use the
`Windows Subsystem for Linux (WSL) <https://docs.microsoft.com/en-us/windows/wsl/install>`_.
The code can be obtained from the NanoPhysics Centers Github repository.
You will need to have **git** installed to clone the repository using the command:

.. code-block:: console

   git clone -b main https://github.com/NanoPhysicsCenter/RUMDEED.git RUMDEED

Where the last "RUMDEED" is the name of the directory you want to create and store the code in.
Once you have to code you will need the GNU Compiler Collection (**gfortran** and **gcc**) and **make** to
compile the code into an executable file. The Intel compilers can be used instead of the GNU compilers (see below),
but the code has mainly been tested with the GNU compilers. Once the compilers and make are installed you can
go into the directory "RUMDEED" and enter the command:

.. code-block:: console
   
   cd RUMDEED
   make

This will compile the code and create an executable named "RUMDEED.out" inside the build directory.
When you have the executable built you can :ref:`run the code <run>`.

To use a different compiler the variables **CC** and **FC** can be passed to the make command.
For example to compile with the Intel compilers:

.. code-block:: console

   make CC=icc FC=ifort

This documentation can also be generate using the make file with:

.. code-block:: console

   make docs

This will produce html and pdf files in the docs/build directory. You must have **Python** and **venv** installed to generate the HTML documentation.
To make the pdf file **xelatex** and **xindy** must be installed.

.. index:: build, make, git, compiler, gfortran, github, Intel, ifort
