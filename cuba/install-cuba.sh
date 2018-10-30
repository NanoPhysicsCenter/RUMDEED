#!/bin/bash

installdir=$PWD
cd 4.2
./configure --prefix="$installdir"
make
make install
make clean
rm config.h
