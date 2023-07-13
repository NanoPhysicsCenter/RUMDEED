.. _build:

Building the code
=================

The code has been developed to run specifically on Linux systems. If you are a Windows user,
it is recommended to utilize the Windows Subsystem for Linux (WSL)
for running the code. You can find installation instructions for WSL at the following link:
`Windows Subsystem for Linux (WSL) <https://docs.microsoft.com/en-us/windows/wsl/install>`_.

To obtain the RUMDEED code, you can access it from the NanoPhysics Centers Github repository.
Prior to cloning the repository, ensure that you have the **git** tool installed.
Use the following command to clone the repository:

.. code-block:: console

   git clone -b main https://github.com/NanoPhysicsCenter/RUMDEED.git RUMDEED

In the command mentioned above, "RUMDEED" represents the name of the directory in which you wish to create and store the code.
You can replace this with any other name of your choice.

Once you have acquired the code, you will need to have the GNU Compiler Collection (**gfortran** and **gcc**) and **make**
installed on your system in order to compile the code and generate an executable file.
Although it is possible to use the Intel compilers instead of the GNU compilers,
it is important to note that the code has primarily been tested with the GNU compilers.

After installing the necessary build utilities, navigate to the "RUMDEED" directory and execute the following command:

.. code-block:: console
   
   make

This command will compile the code and generate an executable file named *RUMDEED.out* within the build directory.
Once the compilation is complete, you can proceed to run the code as explained in the documentation.

If you prefer to use a different compiler, you can specify the **CC** and **FC** variables when executing the make command.
For instance, if you want to compile using the Intel compilers, utilize the following command:

.. code-block:: console

   make CC=icx FC=ifx

Please refer to the :ref:`"Input files and running the code" <run>` in the documentation for further instructions on running the code.

Cuba library
------------
The code uses the `Cuba library <https://feynarts.de/cuba/>`_ for numerical integration :cite:p:`HAHN200578`.
When executing the make file, you will be prompted to download and build the Cuba library.
It is necessary to have **wget** installed on your system to have the make file automatically download and build the library.
In case the download encounters any issues or if you prefer to manually download the library,
you can obtain the file 'Cuba-4.2.2.tar.gz' directly from the Cuba website.
After acquiring the file, please place it into the cuba folder. To build the library, simply execute the 'install-cuba.sh' script.

Documentation
-------------
This documentation can also be generate using the make file with the command:

.. code-block:: console

   make docs

This will produce html and pdf files in the docs/build directory. You must have **Python** and **venv** installed to generate the HTML documentation.
To make the pdf file **xelatex** and **xindy** must be installed.

.. index:: build, make, git, compiler, gfortran, github, Intel, ifort, ifx, Cuba library