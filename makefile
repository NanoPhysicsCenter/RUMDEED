# makefile for Vaccum
# Kristinn Torfason
# 02.09.09

all:
	cd build; make all
clean:
	cd build; make clean
test:
	cd tests; ./compilers.py
install:
	cd build; make install
doc:
	cd doc; make all
