#!/bin/bash

installdir=$PWD
cd 4.2.1
./configure --prefix="${installdir// /\\ }" CC=$1 FC=$2 && make clean && make lib -j 8 && make install && make clean
