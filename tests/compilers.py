#!/usr/bin/python3

import os
import subprocess

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

os.chdir('../build/')
make_process = subprocess.Popen("make clean all FCOMPILER=ifort", shell=True, stderr=subprocess.STDOUT)
if make_process.wait() != 0:
    print(bcolors.FAIL + 'Intel FAILED' + bcolors.ENDC)
else:
    print(bcolors.OKGREEN + 'Intel OK' + bcolors.ENDC)

print('')

make_process = subprocess.Popen("make clean all FCOMPILER=gfortran", shell=True, stderr=subprocess.STDOUT)
if make_process.wait() != 0:
    print(bcolors.FAIL + 'GNU FAILED' + bcolors.ENDC)
else:
    print(bcolors.OKGREEN + 'GNU OK' + bcolors.ENDC)

print('')

make_process = subprocess.Popen("make clean all FCOMPILER=pgfortran", shell=True, stderr=subprocess.STDOUT)
if make_process.wait() != 0:
    print(bcolors.FAIL + 'PGI FAILED' + bcolors.ENDC)
else:
    print(bcolors.OKGREEN + 'PGI OK' + bcolors.ENDC)

print('')
