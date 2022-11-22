# makefile for RUMDEED
# Kristinn Torfason
# 02.09.09

.PHONY: all clean test install doc docs cuba

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
docs:
	cd docs; make html; make latexpdf
