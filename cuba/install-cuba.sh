#!/bin/bash

installdir=$PWD
cd 4.2.1
./configure --prefix="${installdir// /\\ }" CC=gcc-10 FC=gfortran-10
make clean
make lib -j8
make install
make clean
#rm config.h
