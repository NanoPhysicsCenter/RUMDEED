# makefile for RUMDEED
# Kristinn Torfason
# 12.10.09

# Compiler to use
ifeq ($(origin FC),default)
FC = gfortran
endif
ifeq ($(origin CC),default)
CC = gcc
endif

GIT_VERSION := $(shell git describe --dirty --always --tags || echo `date +'%d/%m/%Y'`)

CUBA_LIB=../cuba/lib/libcuba.a

#UNIT_TEST = 1
OPENMP = 1 # Comment out to disable OpenMP | changing value does nothing

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
ifdef OPENMP
FFLAGS+=-qopenmp
endif
endif

# NEW INTEL FORTRAN
ifeq ($(FC),ifx)
$(info Using New Intel Fortran Compiler)
FFLAGS=-O2 -xHost -ipo
ifdef OPENMP
FFLAGS+=-qopenmp
endif
endif

# GNU FORTRAN
# ifeq ($(FC),gfortran)
ifneq (,$(findstring gfortran,$(FC)))
$(info Using GNU Fortran)
FFLAGS=-O3 -march=native -flto=auto
#FFLAGS=-g -fcheck=all -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow
ifdef OPENMP
FFLAGS+=-fopenmp
endif
endif

# PGI FORTRAN
ifeq ($(FC),pgfortran)
$(info Using PGI Fortran)
FFLAGS=-fast
#FFLAGS=-fast -Minfo -Mllvm
#FFLAGS=-g -gopt -Mbounds -Mchkfpstk -Mchkfpstk -Mchkfpstk -Melf -traceback
ifdef OPENMP
FFLAGS+=-mp
endif
endif

# NV FORTRAN
ifeq ($(FC),nvfortran)
$(info Using NVidia Fortran)
#FFLAGS=-O3 -fastsse -Mipa=fast,inline -Munroll
FFLAGS=-O3 -fastsse
ifdef OPENMP
FFLAGS+=-mp
endif
endif

ifeq ($(FC),flang)
$(info Using AMD Fortran)
FFLAGS=-O3 -march=native -D_AOCC_=1
ifdef OPENMP
FFLAGS+=-fopenmp
endif
endif

ifdef UNIT_TEST
FFLAGS+=-D_UNIT_TEST_
endif

# Chunk size for OpenMP
#FFLAGS+=-DCHUNK_SIZE=250

# Name of the output file
PROJECT_NAME=RUMDEED

# makefile
all: $(PROJECT_NAME).out unit_test.out

unit_test.out: $(PROJECT_NAME).out ../src/unit_tests.F90
			$(FC) $(FFLAGS) -D_GIT_VERSION_=\"${GIT_VERSION}\" \
			../src/unit_tests.F90 mod_global.o mod_verlet.o mod_pair.o \
			mod_unit_tests.o mod_photo_emission.o mod_field_emission.o mod_field_emission_v2.o \
			mod_field_thermo_emission.o mod_kevin_rjgtf.o \
			mod_emission_tip.o mod_hyperboloid_tip.o mod_field_emission_2D.o \
			mod_work_function.o mod_collisions.o \
			mod_velocity.o \
			$(CUBA_LIB) -lpthread -lm -ldl -o unit_test.out

$(PROJECT_NAME).out: ../src/main.F90 ../cuba/lib/libcuba.a \
			mod_global.mod mod_global.o \
			mod_verlet.mod mod_verlet.o \
			mod_pair.mod mod_pair.o	\
			mod_unit_tests.mod mod_unit_tests.o \
			mod_photo_emission.mod mod_photo_emission.o \
			mod_field_emission.mod mod_field_emission.o \
			mod_emission_tip.mod mod_emission_tip.o \
			mod_field_emission_v2.mod mod_field_emission_v2.o \
			mod_field_thermo_emission.mod mod_field_thermo_emission.o \
			mod_field_emission_2D.mod mod_field_emission_2D.o \
			mod_work_function.mod mod_work_function.o \
			mod_collisions.mod mod_collisions.o \
			mod_kevin_rjgtf.mod mod_kevin_rjgtf.o \
			mod_velocity.mod mod_velocity.o \
			mod_manual_emission.mod mod_manual_emission.o
		$(FC) $(FFLAGS) -D_GIT_VERSION_=\"${GIT_VERSION}\" \
		    ../src/main.F90 mod_global.o mod_verlet.o mod_pair.o \
			mod_unit_tests.o mod_photo_emission.o mod_field_emission.o mod_field_emission_v2.o \
			mod_field_thermo_emission.o mod_kevin_rjgtf.o \
			mod_emission_tip.o mod_hyperboloid_tip.o mod_field_emission_2D.o \
			mod_work_function.o mod_collisions.o \
			mod_velocity.o mod_manual_emission.o \
			$(CUBA_LIB) -lpthread -lm -ldl -o $(PROJECT_NAME).out

../cuba/Cuba-4.2.2.tar.gz:
		cd ../cuba/; ./get-cuba.sh; md5sum -c md5hash.md5

../cuba/Cuba-4.2.2/cuba.h: ../cuba/Cuba-4.2.2.tar.gz
		cd ../cuba/; tar -xvzf Cuba-4.2.2.tar.gz

../cuba/lib/libcuba.a: ../cuba/Cuba-4.2.2/cuba.h
		cd ../cuba/; chmod +x ./install-cuba.sh; ./install-cuba.sh ${CC} ${FC}

mod_global.mod mod_global.o: ../src/mod_global.F90
		$(FC) -c $(FFLAGS) ../src/mod_global.F90

mod_ic.mod mod_ic.o: ../src/mod_ic.f90 mod_global.o mod_hyperboloid_tip.o mod_verlet.o
		$(FC) -c $(FFLAGS) ../src/mod_ic.f90

mod_verlet.mod mod_verlet.o: ../src/mod_verlet.F90 mod_global.mod mod_pair.mod mod_collisions.mod
		$(FC) -c $(FFLAGS) ../src/mod_verlet.F90

mod_pair.mod mod_pair.o: ../src/mod_pair.F90 mod_global.mod
		$(FC) -c $(FFLAGS) ../src/mod_pair.F90

mod_unit_tests.mod mod_unit_tests.o: ../src/mod_unit_tests.f90 mod_global.mod mod_verlet.mod mod_pair.mod mod_field_emission_v2.mod
		$(FC) -c $(FFLAGS) ../src/mod_unit_tests.f90

mod_photo_emission.mod mod_photo_emission.o: ../src/mod_photo_emission.f90 mod_global.mod mod_verlet.mod mod_pair.mod mod_velocity.mod
		$(FC) -c $(FFLAGS) ../src/mod_photo_emission.f90

mod_field_emission.mod mod_field_emission.o: ../src/mod_field_emission.f90 mod_global.mod mod_verlet.mod mod_pair.mod mod_velocity.mod
		$(FC) -c $(FFLAGS) ../src/mod_field_emission.f90

mod_field_emission_v2.mod mod_field_emission_v2.o: ../src/mod_field_emission_v2.F90 mod_global.mod mod_verlet.mod mod_pair.mod mod_work_function.mod mod_velocity.mod
		$(FC) -c $(FFLAGS) ../src/mod_field_emission_v2.F90

mod_field_thermo_emission.mod mod_field_thermo_emission.o: ../src/mod_field_thermo_emission.F90 mod_global.mod mod_verlet.mod mod_pair.mod mod_work_function.mod mod_kevin_rjgtf.mod mod_velocity.mod
		$(FC) -c $(FFLAGS) ../src/mod_field_thermo_emission.F90

mod_kevin_rjgtf.mod mod_kevin_rjgtf.o: ../src/mod_kevin_rjgtf.F90
		$(FC) -c $(FFLAGS) ../src/mod_kevin_rjgtf.F90

mod_field_emission_2D.mod mod_field_emission_2D.o: ../src/mod_field_emission_2D.f90 mod_global.mod mod_verlet.mod mod_pair.mod mod_velocity.mod
		$(FC) -c $(FFLAGS) ../src/mod_field_emission_2D.f90

mod_emission_tip.mod mod_emission_tip.o: ../src/mod_emission_tip.f90 mod_global.mod mod_verlet.mod mod_pair.mod mod_hyperboloid_tip.mod mod_ic.mod mod_velocity.mod
		$(FC) -c $(FFLAGS) ../src/mod_emission_tip.f90

mod_hyperboloid_tip.mod mod_hyperboloid_tip.o: ../src/mod_hyperboloid_tip.f90 mod_global.mod
		$(FC) -c $(FFLAGS) ../src/mod_hyperboloid_tip.f90

mod_manual_emission.mod mod_manual_emission.o: ../src/mod_manual_emission.f90 mod_global.mod mod_verlet.mod mod_pair.mod
		$(FC) -c $(FFLAGS) ../src/mod_manual_emission.f90

mod_work_function.mod mod_work_function.o: ../src/mod_work_function.F90
		$(FC) -c $(FFLAGS) ../src/mod_work_function.F90

mod_velocity.mod mod_velocity.o: ../src/mod_velocity.f90
		$(FC) -c $(FFLAGS) ../src/mod_velocity.f90

mod_collisions.mod mod_collisions.o: ../src/mod_collisions.F90 mod_global.mod mod_pair.mod mod_velocity.mod
		$(FC) -c $(FFLAGS) ../src/mod_collisions.F90

clean:
		rm -f *.o
		rm -f *.mod
		rm -f *.smod
		rm -f $(PROJECT_NAME).out
		rm -f unit_test.out
		rm -f $(CUBA_LIB)
		rm -f ../cuba/Cuba-4.2.2/libcuba.a
		cd ../cuba/Cuba-4.2.2/ && make clean

install:
		cp $(PROJECT_NAME).out ../data/
