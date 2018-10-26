# makefile for Vacuum
# Kristinn Torfason
# 12.10.09

# Compiler to use
ifndef FCOMPILER
#FCOMPILER=ifort
FCOMPILER=gfortran-8
#FCOMPILER=pgfortran
endif

CUBA_LIB=/home/kristor/software/cuba-4.2/lib/libcuba.a

#UNIT_TEST = 1
OPENMP = 1

# Compile flags
# UNKNOWN COMPILERS
ifndef FFLAGS
FFLAGS=-O2
endif

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
FFLAGS=-Ofast -m64 -flto -march=native -g
#FFLAGS=-g -fcheck=all -Wall -fbacktrace -ffpe-trap=zero,overflow,underflow,denormal
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

ifdef UNIT_TEST
FFLAGS+=-D_UNIT_TEST_
endif

# Name of the output file
PROJECT_NAME=Vacuum-MD

# makefile
all: $(PROJECT_NAME).out

$(PROJECT_NAME).out: ../src/main.F90 mod_global.mod mod_global.o \
			mod_verlet.mod mod_verlet.o \
			mod_pair.mod mod_pair.o	\
			mod_unit_tests.mod mod_unit_tests.o \
			mod_photo_emission.mod mod_photo_emission.o \
			mod_field_emission.mod mod_field_emission.o \
			mod_field_emission_tip.mod mod_field_emission_tip.o \
			mod_field_emission_v2.mod mod_field_emission_v2.o \
			mod_field_emission_2D.mod mod_field_emission_2D.o \
			mod_field_emission_v2@smod_work_function.smod smod_work_function.o \
			mod_field_emission_v2@smod_mc_integration.smod smod_mc_integration.o \
			mod_field_emission_v2@smod_metro_algo.smod smod_metro_algo.o mod_ziggurat.o
		$(FCOMPILER) $(FFLAGS) ../src/main.F90 mod_global.o mod_verlet.o mod_pair.o \
			mod_unit_tests.o mod_photo_emission.o mod_field_emission.o mod_field_emission_v2.o \
			mod_field_emission_tip.o mod_hyperboloid_tip.o mod_field_emission_2D.o \
			smod_work_function.o smod_mc_integration.o smod_metro_algo.o mod_ziggurat.o \
			$(CUBA_LIB) -lpthread -lm -ldl -o $(PROJECT_NAME).out

mod_global.mod mod_global.o: ../src/mod_global.F90
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_global.F90

mod_verlet.mod mod_verlet.o: ../src/mod_verlet.F90 mod_global.mod mod_pair.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_verlet.F90

mod_pair.mod mod_pair.o: ../src/mod_pair.F90 mod_global.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_pair.F90

mod_unit_tests.mod mod_unit_tests.o: ../src/mod_unit_tests.f90 mod_global.mod mod_verlet.mod mod_pair.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_unit_tests.f90

mod_photo_emission.mod mod_photo_emission.o: ../src/mod_photo_emission.f90 mod_global.mod mod_verlet.mod mod_pair.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_photo_emission.f90

mod_field_emission.mod mod_field_emission.o: ../src/mod_field_emission.f90 mod_global.mod mod_verlet.mod mod_pair.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_field_emission.f90

mod_field_emission_v2.mod mod_field_emission_v2.o: ../src/mod_field_emission_v2.f90 mod_global.mod mod_verlet.mod mod_pair.mod ziggurat.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_field_emission_v2.f90

mod_field_emission_2D.mod mod_field_emission_2D.o: ../src/mod_field_emission_2D.f90 mod_global.mod mod_verlet.mod mod_pair.mod
			$(FCOMPILER) -c $(FFLAGS) ../src/mod_field_emission_2D.f90

mod_field_emission_tip.mod mod_field_emission_tip.o: ../src/mod_field_emission_tip.f90 mod_global.mod mod_verlet.mod mod_pair.mod mod_hyperboloid_tip.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_field_emission_tip.f90

mod_hyperboloid_tip.mod mod_hyperboloid_tip.o: ../src/mod_hyperboloid_tip.f90 mod_global.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_hyperboloid_tip.f90

mod_therminoic_emission.mod mod_therminoic_emission.o: ../src/mod_thermionic_emission.f90 mod_global.mod mod_verlet.mod mod_pair.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_thermionic_emission.f90

mod_field_emission_v2@smod_work_function.smod smod_work_function.o: ../src/smod_work_function.F90 mod_field_emission_v2.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/smod_work_function.F90

mod_field_emission_v2@smod_mc_integration.smod smod_mc_integration.o: ../src/smod_mc_integration.F90 mod_field_emission_v2.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/smod_mc_integration.F90

mod_field_emission_v2@smod_metro_algo.smod smod_metro_algo.o: ../src/smod_metro_algo.F90 mod_field_emission_v2.mod
		$(FCOMPILER) -c $(FFLAGS) ../src/smod_metro_algo.F90

ziggurat.o ziggurat.mod: ../src/mod_ziggurat.f90
		$(FCOMPILER) -c $(FFLAGS) ../src/mod_ziggurat.f90

clean:
		rm -f *.o
		rm -f *.mod
		rm -f *.smod
		rm -f $(PROJECT_NAME).out

install:
		cp $(PROJECT_NAME).out ../data/
