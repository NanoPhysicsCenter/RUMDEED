.. _build:

Building the code
=================

The code has been developed to run specifically on Linux systems. If you are a Windows user,
it is recommended to utilize the Windows Subsystem for Linux (WSL)
for running the code. You can find installation instructions for WSL at the following link:
`Windows Subsystem for Linux (WSL) <https://docs.microsoft.com/en-us/windows/wsl/install>`_.

To obtain the RUMDEED code, you can access it from the NanoPhysics Center's GitHub repository.
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

This command will compile the code and generate an executable file named **RUMDEED.out** within the build directory.
Once the compilation is complete, you can proceed to run the code as explained in the documentation.
The unit and integration tests can be run with the command ``make test``.

If you prefer to use a different compiler, you can specify the **CC** and **FC** variables when executing the make command.
For instance, if you want to compile using the Intel compilers, utilize the following command:

.. code-block:: console

   make CC=icx FC=ifx

Please refer to the :ref:`Input files and running the code <run>` section in the documentation for further instructions on running the code.

Building for NVIDIA GPUs (OpenACC)
----------------------------------
The particle-particle interaction kernel (the :math:`O(N^2)` Coulomb and image charge
force calculation in *mod_verlet.F90*) can be offloaded to an NVIDIA GPU using OpenACC.
This requires the **nvfortran** compiler from the
`NVIDIA HPC SDK <https://developer.nvidia.com/hpc-sdk>`_. To build with GPU offload, use:

.. code-block:: console

   make CC=gcc FC=nvfortran ACC=gpu

The GPU compute capability can optionally be given with the **GPU_CC** variable, for
example ``make ACC=gpu GPU_CC=cc80`` for an A100.
This can be combined with OpenMP (which is enabled by default): the offloaded kernel runs
on the GPU while the remaining parallel loops (emission, position and velocity updates)
keep running on the CPU cores. The regular CPU-only builds are unaffected: when compiling
with gfortran or without the **ACC** variable the OpenACC directives are ignored and the
OpenMP code path is used as before.

The GPU kernels implement the planar geometry (the vacuum field *field_E_planar* together
with the image charge partners of *Force_Image_charges_v2*), which covers the planar
emission modes (photo, field, thermionic, field-thermo, 2D and manual emission), and the
hyperboloid tip geometry (*field_E_Hyperboloid* with *Sphere_IC_field*) used by emission
mode 3. Runs using the cylindrical tip or torus geometries automatically fall back to the
OpenMP code path at runtime: those geometries interpolate their fields from a mesh with
kd-tree searches, which cannot run inside a GPU kernel.

For planar field emission (mode 10), also set **MH_BATCH = .TRUE.** in the input file so
that the emission-side field evaluations run in batches on the GPU as well (see the
:ref:`Input files and running the code <run>` section); without it the emission remains
on the CPU and dominates the run time.

For testing the OpenACC code path on a machine without a GPU, the kernels can be compiled
to run on the CPU cores instead:

.. code-block:: console

   make CC=gcc FC=nvfortran ACC=multicore OPENMP=no

POLARSO Laplace solver (optional)
---------------------------------
The POLARSO Laplace solver solves the Poisson equation on a grid with Intel MKL's
PARDISO. Because of the MKL dependency it is a compile-time option that is disabled by
default; without it a stub module is compiled and MKL is not needed at all. To enable it:

.. code-block:: console

   make POLARSO=yes

The MKL installation is found through the **MKLROOT** variable (default
``/opt/intel/oneapi/mkl/latest``), which can be overridden on the make command line or
from the environment. Builds with **POLARSO=yes** use MKL's ILP64 interface and therefore
8-byte default integers. At runtime the solver is activated by placing a ``polarso``
namelist file with ``use_polarso = .true.`` in the run directory; a binary built without
POLARSO stops with an error message if such a file is present. The ``polarso`` namelist
also accepts the grid parameters ``polarso_dim``, ``polarso_pos``, ``polarso_step`` and
``polarso_padding``, and the flag ``write_field_files`` (see mod_polarso.F90).

Cuba library
------------
The code uses the `Cuba library <https://feynarts.de/cuba/>`_ for numerical integration :cite:p:`HAHN200578`.
When executing the makefile, you will be prompted to download and build the Cuba library.
It is necessary to have **wget** installed on your system to have the makefile automatically download and build the library.
In case the download encounters any issues or if you prefer to manually download the library,
you can obtain the file **Cuba-4.2.2.tar.gz** directly from the Cuba website.
After acquiring the file, please place it into the cuba folder. To build the library, simply execute the **install-cuba.sh** script.

Documentation
-------------
This documentation can also be generated using the makefile with the command:

.. code-block:: console

   make docs

This will produce HTML and PDF files in the docs/build directory. You must have **Python** and **venv** installed to generate the HTML documentation.
To make the PDF file, **xelatex** and **xindy** must be installed.

.. index:: build, make, git, compiler, gfortran, github, Intel, ifort, ifx, Cuba library, Windows Subsystem for Linux (WSL), Python, venv, xelatex, xindy, GPU, OpenACC, nvfortran, NVIDIA HPC SDK