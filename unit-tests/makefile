# makefile for Vacuum
# Kristinn Torfason
# 06.11.2018

# Compiler to use
ifndef FCOMPILER
#FCOMPILER=ifort
FCOMPILER=gfortran
#FCOMPILER=pgfortran
endif

GIT_VERSION := $(shell git describe --dirty --always --tags)

CUBA_LIB=../cuba/lib/libcuba.a

#UNIT_TEST = 1
OPENMP = 1

# INTEL FORTRAN
ifeq ($(FCOMPILER),ifort)
$(info Using Intel Fortran)
FFLAGS=-O2 -xHost -ip -ipo -prec-div -prec-sqrt -fp-model precise -g
#-fltconsistency
#-prec-div -prec-sqrt
#-fp-model precise
#FFLAGS=-check all -g -warn interfaces -warn all
ifdef OPENMP
FFLAGS+=-qopenmp
endif
endif

# GNU FORTRAN
# ifeq ($(FCOMPILER),gfortran)
ifneq (,$(findstring gfortran,$(FCOMPILER)))
$(info Using GNU Fortran)
#FFLAGS=-Ofast -m64 -flto -march=native
FFLAGS=-g -fcheck=all -Wall -fbacktrace -ffpe-trap=zero,overflow,underflow,denormal
ifdef OPENMP
FFLAGS+=-fopenmp
endif
endif

# PGI FORTRAN
ifeq ($(FCOMPILER),pgfortran)
$(info Using PGI Fortran)
FFLAGS=-O2 -fast
#FFLAGS=-g -gopt -Mbounds -Mchkfpstk -Mchkfpstk -Mchkfpstk -Melf -traceback
ifdef OPENMP
FFLAGS+=-mp
endif
endif

all: verlet-tests.out

verlet-tests.out: verlet-tests.f90 mod_global.o mod_verlet.o mod_pair.o ../cuba/lib/libcuba.a
				$(FCOMPILER) $(FFLAGS) -D_GIT_VERSION_=\"${GIT_VERSION}\" verlet-tests.f90 \
				mod_global.o mod_verlet.o mod_pair.o $(CUBA_LIB) -lpthread -lm -ldl \
				-o verlet-tests.out

mod_global.mod mod_global.o: ../src/mod_global.F90
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_global.F90

mod_verlet.mod mod_verlet.o: ../src/mod_verlet.F90 mod_global.mod mod_pair.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_verlet.F90

mod_pair.mod mod_pair.o: ../src/mod_pair.F90 mod_global.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_pair.F90
