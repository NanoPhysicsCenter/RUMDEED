.. _build:

Building the code
=================

The code was built to run on Linux systems. For Windows users it's recommended to use the
`Windows Subsystem for Linux (WSL) <https://docs.microsoft.com/en-us/windows/wsl/install>`_.

The RUMDEED code can be obtained from the NanoPhysics Centers Github repository.
You will need to have **git** installed to clone the repository using the command:

.. code-block:: console

   git clone -b main https://github.com/NanoPhysicsCenter/RUMDEED.git RUMDEED

Where the last "RUMDEED" is the name of the directory you want to create and store the code in.
Once you have to code you will need the GNU Compiler Collection (**gfortran** and **gcc**) and **make** to
compile the code into an executable file. The Intel compilers can be used instead of the GNU compilers (see below),
but the code has mainly been tested with the GNU compilers. Once the necessary build utilities are installed you can
go into the directory "RUMDEED" and enter the command:

.. code-block:: console
   
   make

This will compile the code and create an executable named "RUMDEED.out" inside the build directory.
When you have the executable built you can :ref:`run the code <run>`.

To use a different compiler the variables **CC** and **FC** can be passed to the make command.
For example to compile with the Intel compilers:

.. code-block:: console

   make CC=icc FC=ifort

Cuba library
------------
The code uses the `Cuba library <https://feynarts.de/cuba/>`_ for numerical integration :cite:p:`HAHN200578`. The make file will automatically download and build the code
for the Cuba library. To download the code **wget** needs to be installed. If the download fails the file 'Cuba-4.2.2.tar.gz' can be download manually from the Cuba website and
placed into the cuba folder. The script 'install-cuba.sh' can then be run to build the library.

Documentation
-------------
This documentation can also be generate using the make file with the command:

.. code-block:: console

   make docs

This will produce html and pdf files in the docs/build directory. You must have **Python** and **venv** installed to generate the HTML documentation.
To make the pdf file **xelatex** and **xindy** must be installed.

.. index:: build, make, git, compiler, gfortran, github, Intel, ifort, Cuba library