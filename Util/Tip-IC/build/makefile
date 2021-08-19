# makefile for IC
# Kristinn Torfason
# 02.07.14

# Compiler to use
FCOMPILER=ifort
#FCOMPILER=gfortran

# Compile flags
#FFLAGS=-O2 -parallel -xHost -openmp -mcmodel=large -shared-intel
#FFLAGS=-O3 -xHost -qopenmp -ip -ipo -par-affinity=granularity=core,scatter
FFLAGS=-check all -g -traceback -warn all
#FFLAGS=-O3 -m64 -fopenmp -flto -march=native

# Name of the output file
PROJECT_NAME=IC

# makefile
all: $(PROJECT_NAME).out

$(PROJECT_NAME).out: main.o conicpr.o cons.o auxil.o cgamma.o
		$(FCOMPILER) $(FFLAGS) main.o conicpr.o cons.o auxil.o bessel.o gammaCHI.o air.o cgamma.o \
		-o $(PROJECT_NAME).out

main.o: ../src/main.f90 auxil.o cons.o cgamma.o conicpr.o
		$(FCOMPILER) -c $(FFLAGS) ../src/main.f90
conicpr.o: ../src/conicpr.f90 cons.o auxil.o gammaCHI.o bessel.o
		$(FCOMPILER) -c $(FFLAGS) ../src/conicpr.f90
cons.o: ../src/cons.f90 ../src/constants
		$(FCOMPILER) -c $(FFLAGS) ../src/cons.f90
auxil.o: ../src/auxil.f90
		$(FCOMPILER) -c $(FFLAGS) ../src/auxil.f90
bessel.o: ../src/bessel.f90 air.o
		$(FCOMPILER) -c $(FFLAGS) ../src/bessel.f90
gammaCHI.o: ../src/gammaCHI.f90
		$(FCOMPILER) -c $(FFLAGS) ../src/gammaCHI.f90
air.o: ../src/air.f90
		$(FCOMPILER) -c $(FFLAGS) ../src/air.f90
cgamma.o: ../src/cgamma.f
		$(FCOMPILER) -c $(FFLAGS) ../src/cgamma.f

install:
		mkdir -p ../run/
		cp $(PROJECT_NAME).out ../run/

clean:
		rm -f *.o
		rm -f *.mod
		rm -f *__genmod.f90
		rm -f $(PROJECT_NAME).out
