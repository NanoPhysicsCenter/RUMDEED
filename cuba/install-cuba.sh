#!/bin/bash

installdir=$PWD
cd Cuba-4.2.2
# -std=gnu11 is needed because Cuba defines its own bool type, which clashes
# with the bool keyword of C23 (the default of GCC 15 and later).
./configure --prefix="${installdir// /\\ }" CC=$1 FC=$2 CFLAGS="-O3 -fomit-frame-pointer -std=gnu11" && make clean && make lib && make install && make clean
