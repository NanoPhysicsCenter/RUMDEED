# makefile for Vacuum
# Kristinn Torfason
# 12.10.09

# Compiler to use
ifndef FCOMPILER
FCOMPILER=ifort
#FCOMPILER=gfortran
#FCOMPILER=pgfortran
endif

# Include and Library paths
MKLPATH=$(MKLROOT)/lib/intel64
MKLINCLUDE=$(MKLROOT)/include


# Compile flags
# UNKNOWN COMPILERS
FFLAGS=-O2

# INTEL FORTRAN
ifeq ($(FCOMPILER),ifort)
$(info Using Intel Fortran)
FFLAGS=-O2 -xHost -qopenmp -ip -ipo
MKL_LIBS=-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread
#FFLAGS=-check all -g -warn all -qopenmp
endif

# GNU FORTRAN
ifeq ($(FCOMPILER),gfortran)
$(info Using GNU Fortran)
FFLAGS=-O3 -fopenmp -ffast-math -m64 -funroll-loops -flto -march=native -fno-range-check
MKL_LIBS=-Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread

endif

# PGI FORTRAN
ifeq ($(FCOMPILER),pgfortran)
$(info Using PGI Fortran)
FFLAGS=-O2 -fast -mp -pgf90libs -I${MKLROOT}/include
MKL_LIBS=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_pgi_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group
endif

# Name of the output file
PROJECT_NAME=Vacuum-MD

# makefile
all: $(PROJECT_NAME).out

$(PROJECT_NAME).out: ../src/main.F90 mod_global.mod mod_verlet.mod mod_pair.mod
		$(FCOMPILER) $(FFLAGS) ../src/main.F90 mod_global.o mod_verlet.o mod_pair.o \
			$(MKL_LIBS) -lpthread -lm -ldl -o $(PROJECT_NAME).out

mod_global.mod: ../src/mod_global.F90
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_global.F90 -I$(MKLINCLUDE)

mod_verlet.mod: ../src/mod_verlet.F90 mod_global.mod mod_pair.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_verlet.F90 -I$(MKLINCLUDE)

mod_pair.mod: ../src/mod_pair.F90 mod_global.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_pair.F90 -I$(MKLINCLUDE)

clean:
		rm -f *.o
		rm -f *.mod
		rm -f $(PROJECT_NAME).out

install:
		cp $(PROJECT_NAME).out ../data/
