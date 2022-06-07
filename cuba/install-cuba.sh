#!/bin/bash

installdir=$PWD
cd Cuba-4.2.2
./configure --prefix="${installdir// /\\ }" CC=$1 FC=$2 && make clean && make lib -j 8 && make install && make clean
