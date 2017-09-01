# makefile for Vacuum
# Kristinn Torfason
# 12.10.09

# Compiler to use
ifndef FCOMPILER
FCOMPILER=ifort
#FCOMPILER=gfortran
#FCOMPILER=pgfortran
endif


# Compile flags
# UNKNOWN COMPILERS
FFLAGS=-O2

# INTEL FORTRAN
ifeq ($(FCOMPILER),ifort)
$(info Using Intel Fortran)
FFLAGS=-O2 -xHost -qopenmp -ip -ipo
#FFLAGS=-check all -g -warn all -qopenmp
endif

# GNU FORTRAN
ifeq ($(FCOMPILER),gfortran)
$(info Using GNU Fortran)
FFLAGS=-O3 -fopenmp -ffast-math -m64 -funroll-loops -flto -march=native
endif

# PGI FORTRAN
ifeq ($(FCOMPILER),pgfortran)
$(info Using PGI Fortran)
FFLAGS=-O2 -mp
endif

# Name of the output file
PROJECT_NAME=Vacuum-MD


# makefile
all: $(PROJECT_NAME).out

$(PROJECT_NAME).out: ../src/main.f90 mod_global.mod
		$(FCOMPILER) $(FFLAGS) ../src/main.f90 mod_global.o\
			-lpthread -lm -ldl -o $(PROJECT_NAME).out

mod_global.mod: ../src/mod_global.f90
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_global.f90

clean:
		rm -f *.o
		rm -f *.mod
		rm -f $(PROJECT_NAME).out

install:
		cp $(PROJECT_NAME).out ../data/
