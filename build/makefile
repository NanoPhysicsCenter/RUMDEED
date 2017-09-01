# makefile for Vacuum
# Kristinn Torfason
# 12.10.09

# Compiler to use
FCOMPILER=ifort

# Compile flags
#FFLAGS=-O2 -parallel -xHost -openmp -mcmodel=large -shared-intel
FFLAGS=-O2 -xHost -qopenmp -ip -ipo
#FFLAGS=-check all -g -warn all -openmp

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
