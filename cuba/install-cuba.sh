#!/bin/bash

installdir=~/software/cuba-4.2
eval installdir=$installdir
cd 4.2
./configure --prefix=$installdir
make
make install
