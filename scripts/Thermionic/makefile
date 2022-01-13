FCOMPILER=gfortran-8
#FFLAGS=-Ofast -flto -march=native
FFLAGS=-O3 -fdefault-double-8 -march=native

all: GenThermField.out
GenThermField.out: GenThermField.f90
		$(FCOMPILER) $(FFLAGS) GenThermField.f90 -lpthread -lm -ldl -o GenThermField.out

clean:
		rm -f *.o
		rm -f *.mod
		rm -f *.smod
		rm -f GenThermField.out
