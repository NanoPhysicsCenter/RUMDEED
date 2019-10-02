#!/bin/bash

installdir=$PWD
cd 4.2
./configure --prefix="${installdir// /\\ }" CC=gcc-6 FC=gfortran-6
make clean
make
make install
make clean
#rm config.h
