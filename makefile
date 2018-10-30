# makefile for Vaccum
# Kristinn Torfason
# 02.09.09

.PHONY: all clean test install doc cuba

all:
	cd build; make all
clean:
	cd build; make clean; cd ../doc; make clean
test:
	cd tests; ./compilers.py
install:
	cd build; make install
doc:
	cd doc; make all
cuba:
	cd cuba; ./install-cuba.sh
