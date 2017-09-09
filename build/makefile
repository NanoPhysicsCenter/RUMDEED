# makefile for Vacuum
# Kristinn Torfason
# 12.10.09

# Compiler to use
ifndef FCOMPILER
#FCOMPILER=ifort
FCOMPILER=gfortran
#FCOMPILER=pgfortran
endif

UNIT_TEST = 1

# Compile flags
# UNKNOWN COMPILERS
ifndef FFLAGS
FFLAGS=-O2
endif

# INTEL FORTRAN
ifeq ($(FCOMPILER),ifort)
$(info Using Intel Fortran)
#FFLAGS=-O2 -xHost -qopenmp -ip -ipo
FFLAGS=-check all -g -warn interfaces -warn all -qopenmp
endif

# GNU FORTRAN
ifeq ($(FCOMPILER),gfortran)
$(info Using GNU Fortran)
#FFLAGS=-O3 -fopenmp -ffast-math -m64 -funroll-loops -flto -march=native -fno-range-check
FFLAGS=-g -fcheck=all -Wall -fopenmp
endif

# PGI FORTRAN
ifeq ($(FCOMPILER),pgfortran)
$(info Using PGI Fortran)
#FFLAGS=-O2 -fast -mp -pgf90libs
FFLAGS=-g -gopt -Mbounds -mp
endif

ifdef UNIT_TEST
FFLAGS+=-D_UNIT_TEST_
endif

# Name of the output file
PROJECT_NAME=Vacuum-MD

# makefile
all: $(PROJECT_NAME).out

$(PROJECT_NAME).out: ../src/main.F90 mod_global.mod mod_verlet.mod mod_pair.mod	mod_unit_tests.mod
		$(FCOMPILER) $(FFLAGS) ../src/main.F90 mod_global.o mod_verlet.o mod_pair.o mod_unit_tests.o \
			-lpthread -lm -ldl -o $(PROJECT_NAME).out

mod_global.mod: ../src/mod_global.F90
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_global.F90

mod_verlet.mod: ../src/mod_verlet.F90 mod_global.mod mod_pair.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_verlet.F90

mod_pair.mod: ../src/mod_pair.F90 mod_global.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_pair.F90

mod_unit_tests.mod: ../src/mod_unit_tests.f90 mod_global.mod mod_verlet.mod mod_pair.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_unit_tests.f90

clean:
		rm -f *.o
		rm -f *.mod
		rm -f $(PROJECT_NAME).out

install:
		cp $(PROJECT_NAME).out ../data/
