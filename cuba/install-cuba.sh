#!/bin/bash

installdir=$PWD
cd 4.2
./configure --prefix="${installdir// /\\ }" CC=gcc-8 FC=gfortran-8
make clean
make
make install
make clean
#rm config.h
