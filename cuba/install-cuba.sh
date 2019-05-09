#!/bin/bash

installdir=$PWD
cd 4.2
./configure --prefix="${installdir// /\\ }" CC=gcc-9.1 FC=gfortran-9.1
make clean
make
make install
make clean
#rm config.h
