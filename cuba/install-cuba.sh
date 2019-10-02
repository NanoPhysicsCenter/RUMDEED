#!/bin/bash

installdir=$PWD
cd 4.2
./configure --prefix="${installdir// /\\ }" CC=gcc FC=gfortran
make clean
make
make install
make clean
#rm config.h
