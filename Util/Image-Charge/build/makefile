# makefile for IC
# Kristinn Torfason
# 02.07.14

# Compiler to use
FCOMPILER=ifort
#FCOMPILER=gfortran

# Compile flags
#FFLAGS=-O2 -parallel -xHost -openmp -mcmodel=large -shared-intel
FFLAGS=-O3 -xHost -qopenmp -ip -ipo -par-affinity=granularity=core,scatter
#FFLAGS=-check all -g -traceback -warn all -openmp
#FFLAGS=-O3 -m64 -fopenmp -flto -march=native

# Name of the output file
PROJECT_NAME=IC

# makefile
all: $(PROJECT_NAME).out

$(PROJECT_NAME).out: main.o conic.o cons.o auxil.o cgamma.o
		$(FCOMPILER) $(FFLAGS) main.o conic.o cons.o auxil.o cgamma.o \
		-o $(PROJECT_NAME).out

main.o: ../main.f90 auxil.o cons.o cgamma.o conic.o
		$(FCOMPILER) -c $(FFLAGS) ../main.f90
conic.o: ../conic.f90 cons.o auxil.o
		$(FCOMPILER) -c $(FFLAGS) ../conic.f90
cons.o: ../cons.f90
		$(FCOMPILER) -c $(FFLAGS) ../cons.f90
auxil.o: ../auxil.f90
		$(FCOMPILER) -c $(FFLAGS) ../auxil.f90
cgamma.o: ../cgamma.f
		$(FCOMPILER) -c $(FFLAGS) ../cgamma.f

install:
		mkdir -p ../run/
		cp $(PROJECT_NAME).out ../run/

clean:
		rm -f *.o
		rm -f *.mod
		rm -f *__genmod.f90
		rm -f $(PROJECT_NAME).out
