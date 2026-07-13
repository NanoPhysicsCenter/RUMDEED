.PHONY: all clean install

# makefile for RUMDEED
# Kristinn Torfason
# 12.10.09

# Compiler to use
# Override default values if nothing else is set
ifeq ($(origin FC),default)
FC = gfortran
endif
ifeq ($(origin CC),default)
CC = gcc
endif

GIT_VERSION := $(shell git describe --dirty --always --tags || echo `date +'%d/%m/%Y'`)

CUBA_LIB=../cuba/lib/libcuba.a

ifndef OPENMP
OPENMP = yes # Set to something else to disable OpenMP
$(info Using OpenMP)
else
$(info OpenMP disabled)
endif

TESTING_MODE ?= 2 # Set to 1 to enable testing mode

# Compile flags
# UNKNOWN COMPILERS
ifndef FFLAGS
FFLAGS=-O2
endif

# INTEL FORTRAN CLASSIC
ifeq ($(FC),ifort)
$(info Using Intel Fortran Classic)
FFLAGS=-O2 -xHost -ip -ipo #-fp-model precise -g
#-fltconsistency
#-prec-div -prec-sqrt
#-fp-model precise
#FFLAGS=-check all -g -warn interfaces -warn all -traceback -ftrapuv
ifneq (,$(findstring yes,$(OPENMP)))
FFLAGS+=-qopenmp
endif
endif

# NEW INTEL FORTRAN
ifeq ($(FC),ifx)
$(info Using New Intel Fortran Compiler)
FFLAGS = -O2 -xHost -ipo
ifneq (,$(findstring yes,$(OPENMP)))
FFLAGS += -fopenmp
endif
endif

# GNU FORTRAN
# ifeq ($(FC),gfortran)
ifneq (,$(findstring gfortran,$(FC)))
$(info Using GNU Fortran)
FFLAGS = -O3 -march=native -flto=auto -ffree-line-length-none
#FFLAGS=-g -fcheck=all -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow
ifneq (,$(findstring yes,$(OPENMP)))
FFLAGS+=-fopenmp
else
FFLAGS+=-fno-openmp
endif
$(info $$OPENMP is [${OPENMP}])
endif

# PGI FORTRAN
ifeq ($(FC),pgfortran)
$(info Using PGI Fortran)
FFLAGS=-fast
#FFLAGS=-fast -Minfo -Mllvm
#FFLAGS=-g -gopt -Mbounds -Mchkfpstk -Mchkfpstk -Mchkfpstk -Melf -traceback
ifneq (,$(findstring yes,$(OPENMP)))
FFLAGS+=-mp
endif
endif

# NV FORTRAN
ifeq ($(FC),nvfortran)
$(info Using NVidia Fortran)
#FFLAGS=-O3 -fastsse -Mipa=fast,inline -Munroll
FFLAGS=-O3 -fastsse
ifneq (,$(findstring yes,$(OPENMP)))
FFLAGS+=-mp
endif
# OpenACC offload of the N-body kernels (see mod_verlet.F90).
# ACC=gpu       Offload to an NVIDIA GPU (can be combined with OPENMP=yes,
#               the OpenMP loops stay on the CPU).
# ACC=multicore Run the OpenACC kernels on the CPU cores instead, useful for
#               testing on machines without a GPU. Use with OPENMP=no.
ifeq ($(ACC),gpu)
$(info Using OpenACC on GPU)
FFLAGS+=-acc=gpu -Minfo=accel
# Optional GPU compute capability, e.g. make ACC=gpu GPU_CC=cc80 for an A100.
# Needed when compiling on a machine without a GPU (e.g. a cluster login
# node), where nvfortran cannot auto-detect the target architecture.
ifdef GPU_CC
FFLAGS+=-gpu=$(GPU_CC)
endif
endif
ifeq ($(ACC),multicore)
$(info Using OpenACC on CPU cores)
FFLAGS+=-acc=multicore -Minfo=accel
endif
endif

ifeq ($(FC),flang)
$(info Using AMD Fortran)
FFLAGS=-O3 -march=native -D_AOCC_=1
ifneq (,$(findstring yes,$(OPENMP)))
FFLAGS+=-fopenmp
endif
endif

# POLARSO Laplace solver (mod_polarso.F90). Optional because it needs the
# Intel MKL library (PARDISO); disabled by default. Enable with:
#   make POLARSO=yes [MKLROOT=/path/to/mkl]
# The solver uses MKL's ILP64 interface, so builds with POLARSO=yes use
# 8-byte default integers. Without it a stub module is compiled instead and
# MKL is not needed at all.
POLARSO ?= no
ifeq ($(POLARSO),yes)
$(info Using POLARSO (requires Intel MKL))
MKLROOT ?= /opt/intel/oneapi/mkl/latest
FFLAGS += -D_POLARSO_=1
MKL_DEPS = mkl_pardiso.mod mkl_pardiso.o \
           mkl_spblas.mod mkl_spblas.o \
           mkl_rci.mod mkl_rci.o
ifneq (,$(findstring gfortran,$(FC)))
FFLAGS += -fdefault-integer-8
MKL_LIBS = -m64 -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp
else ifeq ($(FC),nvfortran)
FFLAGS += -i8
MKL_LIBS = -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5
else
# Intel ifx / ifort
FFLAGS += -i8
MKL_LIBS = -qmkl-ilp64=parallel ${MKLROOT}/lib/libmkl_blas95_ilp64.a ${MKLROOT}/lib/libmkl_lapack95_ilp64.a
endif
else
MKL_DEPS =
MKL_LIBS =
endif

# Chunk size for OpenMP
#FFLAGS+=-DCHUNK_SIZE=250

# Name of the output file
PROJECT_NAME=RUMDEED

# makefile
all: $(PROJECT_NAME).out

$(PROJECT_NAME).out: ../src/main.F90 ../cuba/lib/libcuba.a \
			mod_global.mod mod_global.o \
			mod_cuba_integration.mod mod_cuba_integration.o \
			mod_verlet.mod mod_verlet.o \
			mod_pair.mod mod_pair.o	\
			mod_tests.mod mod_tests.o \
			mod_photo_emission.mod mod_photo_emission.o \
			mod_field_emission.mod mod_field_emission.o \
			mod_emission_tip.mod mod_emission_tip.o \
			mod_field_emission_v2.mod mod_field_emission_v2.o \
			mod_field_thermo_emission.mod mod_field_thermo_emission.o \
			mod_field_emission_2D.mod mod_field_emission_2D.o \
			mod_hyperboloid_tip.mod mod_hyperboloid_tip.o \
			mod_ic.mod mod_ic.o \
			mod_cylindrical_tip.mod mod_cylindrical_tip.o \
			mod_torus.mod mod_torus.o \
			mod_work_function.mod mod_work_function.o \
			mod_polynomialroots.mod mod_polynomialroots.o \
			mod_collisions.mod mod_collisions.o \
			mod_kevin_rjgtf_v2.mod mod_kevin_rjgtf_v2.o \
			mod_velocity.mod mod_velocity.o \
			mod_manual_emission.mod mod_manual_emission.o \
			kdtree2.o kdtree2_module.mod \
			mod_polarso.mod mod_polarso.o \
			$(MKL_DEPS)
		$(FC) $(FFLAGS) -D_GIT_VERSION_=\"${GIT_VERSION}\" -D_TESTING_MODE_=${TESTING_MODE} \
		    ../src/main.F90 mod_global.o mod_cuba_integration.o mod_verlet.o mod_pair.o \
			mod_tests.o mod_photo_emission.o mod_field_emission.o mod_field_emission_v2.o \
			mod_field_thermo_emission.o mod_kevin_rjgtf_v2.o \
			mod_emission_tip.o mod_hyperboloid_tip.o mod_ic.o mod_field_emission_2D.o \
			mod_cylindrical_tip.o mod_torus.o \
			mod_polynomialroots.o \
			mod_work_function.o mod_collisions.o \
			mod_velocity.o mod_manual_emission.o \
			mod_polarso.o \
			kdtree2.o \
			$(CUBA_LIB) $(MKL_LIBS) -lpthread -lm -ldl -o $(PROJECT_NAME).out

../cuba/Cuba-4.2.2.tar.gz:
		cd ../cuba/; ./get-cuba.sh; md5sum -c md5hash.md5

../cuba/Cuba-4.2.2/cuba.h: ../cuba/Cuba-4.2.2.tar.gz
		cd ../cuba/; tar -xvzf Cuba-4.2.2.tar.gz

../cuba/lib/libcuba.a: ../cuba/Cuba-4.2.2/cuba.h
		cd ../cuba/; chmod +x ./install-cuba.sh; ./install-cuba.sh ${CC} ${FC}

mod_global.o: ../src/mod_global.F90
		$(FC) -c $(FFLAGS) ../src/mod_global.F90

mod_global.mod: mod_global.o
	@true

mod_cuba_integration.o: ../src/mod_cuba_integration.F90 mod_global.mod
		$(FC) -c $(FFLAGS) ../src/mod_cuba_integration.F90

mod_cuba_integration.mod: mod_cuba_integration.o
	@true

mod_ic.o: ../src/mod_ic.f90 mod_global.o mod_hyperboloid_tip.o mod_verlet.o
		$(FC) -c $(FFLAGS) ../src/mod_ic.f90

mod_ic.mod: mod_ic.o
	@true

mod_verlet.o: ../src/mod_verlet.F90 mod_global.mod mod_pair.mod mod_collisions.mod mod_hyperboloid_tip.mod \
		../src/acc_ic_planar_series.inc ../src/acc_tip_image_point.inc ../src/acc_tip_ic_force.inc ../src/acc_tip_field_E.inc
		$(FC) -c $(FFLAGS) -I../src ../src/mod_verlet.F90

mod_verlet.mod: mod_verlet.o
	@true

mod_pair.o: ../src/mod_pair.F90 mod_global.mod 
		$(FC) -c $(FFLAGS) ../src/mod_pair.F90

mod_pair.mod: mod_pair.o
	@true

mod_tests.o: ../src/mod_tests.F90 mod_global.mod mod_verlet.mod mod_pair.mod mod_field_emission_v2.mod mod_emission_tip.mod \
             mod_work_function.mod mod_collisions.mod mod_cuba_integration.mod
		$(FC) -c $(FFLAGS) ../src/mod_tests.F90

mod_tests.mod: mod_tests.o
	@true

mod_photo_emission.o: ../src/mod_photo_emission.f90 mod_global.mod mod_verlet.mod mod_pair.mod mod_velocity.mod mod_polarso.mod
		$(FC) -c $(FFLAGS) $(MKL_LIBS) ../src/mod_photo_emission.f90

mod_photo_emission.mod: mod_photo_emission.o
	@true

mod_field_emission.o: ../src/mod_field_emission.f90 mod_global.mod mod_verlet.mod mod_pair.mod mod_velocity.mod
		$(FC) -c $(FFLAGS) ../src/mod_field_emission.f90

mod_field_emission.mod: mod_field_emission.o
	@true

mod_field_emission_v2.o: ../src/mod_field_emission_v2.F90 mod_global.mod mod_cuba_integration.mod mod_verlet.mod mod_pair.mod mod_work_function.mod mod_velocity.mod mod_polarso.mod
		$(FC) -c $(FFLAGS) $(MKL_LIBS) ../src/mod_field_emission_v2.F90

mod_field_emission_v2.mod: mod_field_emission_v2.o
	@true

mod_field_thermo_emission.o: ../src/mod_field_thermo_emission.F90 mod_global.mod mod_cuba_integration.mod mod_verlet.mod mod_pair.mod mod_work_function.mod mod_kevin_rjgtf_v2.mod mod_velocity.mod
		$(FC) -c $(FFLAGS) ../src/mod_field_thermo_emission.F90

mod_field_thermo_emission.mod: mod_field_thermo_emission.o
	@true

mod_kevin_rjgtf_v2.o: ../src/mod_kevin_rjgtf_v2.f90
		$(FC) -c $(FFLAGS) ../src/mod_kevin_rjgtf_v2.f90

mod_kevin_rjgtf_v2.mod: mod_kevin_rjgtf_v2.o
	@true

mod_field_emission_2D.o: ../src/mod_field_emission_2D.f90 mod_global.mod mod_verlet.mod mod_pair.mod mod_velocity.mod
		$(FC) -c $(FFLAGS) ../src/mod_field_emission_2D.f90

mod_field_emission_2D.mod: mod_field_emission_2D.o
	@true

mod_emission_tip.o: ../src/mod_emission_tip.f90 mod_global.mod mod_cuba_integration.mod mod_verlet.mod mod_pair.mod mod_hyperboloid_tip.mod mod_ic.mod mod_velocity.mod mod_kevin_rjgtf_v2.mod
		$(FC) -c $(FFLAGS) ../src/mod_emission_tip.f90

mod_emission_tip.mod: mod_emission_tip.o
	@true

mod_hyperboloid_tip.o: ../src/mod_hyperboloid_tip.f90 mod_global.mod
		$(FC) -c $(FFLAGS) ../src/mod_hyperboloid_tip.f90

mod_hyperboloid_tip.mod: mod_hyperboloid_tip.o
	@true

mod_cylindrical_tip.o: ../src/mod_cylindrical_tip.F90 mod_global.mod mod_cuba_integration.mod mod_verlet.mod mod_pair.mod mod_velocity.mod kdtree2_module.mod
		$(FC) -c $(FFLAGS) ../src/mod_cylindrical_tip.F90

mod_cylindrical_tip.mod: mod_cylindrical_tip.o
	@true

mod_torus.o: ../src/mod_torus.F90 mod_global.mod mod_cuba_integration.mod mod_verlet.mod mod_pair.mod mod_velocity.mod
		$(FC) -c $(FFLAGS) ../src/mod_torus.F90

mod_torus.mod: mod_torus.o
	@true

mod_manual_emission.o: ../src/mod_manual_emission.f90 mod_global.mod mod_verlet.mod mod_pair.mod mod_field_emission_v2.mod
		$(FC) -c $(FFLAGS) ../src/mod_manual_emission.f90

mod_manual_emission.mod: mod_manual_emission.o
	@true

mod_work_function.o: ../src/mod_work_function.F90 mod_global.mod
		$(FC) -c $(FFLAGS) ../src/mod_work_function.F90

mod_work_function.mod: mod_work_function.o
	@true

mod_polynomialroots.o: ../src/mod_polynomialroots.F90
		$(FC) -c $(FFLAGS) ../src/mod_polynomialroots.F90

mod_polynomialroots.mod: mod_polynomialroots.o
	@true

mod_polarso.o: ../src/mod_polarso.F90 mod_global.mod mod_verlet.mod mod_pair.mod mod_work_function.mod $(MKL_DEPS)
		$(FC) -c $(FFLAGS) $(MKL_LIBS) ../src/mod_polarso.F90

mod_polarso.mod: mod_polarso.o
	@true

mod_velocity.o: ../src/mod_velocity.f90 mod_global.mod
		$(FC) -c $(FFLAGS) ../src/mod_velocity.f90

mod_velocity.mod: mod_velocity.o
	@true

mod_collisions.o: ../src/mod_collisions.F90 mod_global.mod mod_pair.mod mod_velocity.mod mod_polynomialroots.mod
		$(FC) -c $(FFLAGS) ../src/mod_collisions.F90

mod_collisions.mod: mod_collisions.o
	@true

kdtree2.o: ../src/kdtree2/kdtree2.f90
		$(FC) -c $(FFLAGS) ../src/kdtree2/kdtree2.f90

kdtree2_module.mod: kdtree2.o
	@true

# mkl
mkl_pardiso.o: $(MKLROOT)/include/mkl_pardiso.f90
		$(FC) -c $(FFLAGS) -I$(MKLROOT)/include $< -o mkl_pardiso.o
		cp mkl_pardiso.mod ../src/ 

mkl_pardiso.mod: mkl_pardiso.o
	@true

mkl_spblas.o: $(MKLROOT)/include/mkl_spblas.f90
		$(FC) -c $(FFLAGS) -I$(MKLROOT)/include $< -o mkl_spblas.o
		cp mkl_spblas.mod ../src/

mkl_spblas.mod: mkl_spblas.o
	@true

mkl_rci.o: $(MKLROOT)/include/mkl_rci.f90
		$(FC) -c $(FFLAGS) -I$(MKLROOT)/include $< -o mkl_rci.o
		cp mkl_rci.mod ../src/

mkl_rci.mod: mkl_rci.o
	@true

clean:
		rm -f *.o
		rm -f *.mod
		rm -f *.smod
		rm -f $(PROJECT_NAME).out
		rm -f $(CUBA_LIB)
		rm -f ../cuba/Cuba-4.2.2/libcuba.a
		cd ../cuba/Cuba-4.2.2/ && make clean

install:
		cp $(PROJECT_NAME).out ../data/
