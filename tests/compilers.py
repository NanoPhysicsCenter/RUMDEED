#!/usr/bin/python3

import os
import subprocess
import time

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
make_process = subprocess.Popen("make clean all install FCOMPILER=ifort UNIT_TEST=1", shell=True, stderr=subprocess.STDOUT)
if make_process.wait() != 0:
    print(bcolors.FAIL + 'Intel FAILED' + bcolors.ENDC)
else:
    print(bcolors.OKGREEN + 'Intel OK' + bcolors.ENDC)
    os.chdir('../data/')
    Vacuum_Process = subprocess.Popen("./Vacuum-MD.out", shell=True, stderr=subprocess.STDOUT)
    Vacuum_Process.wait()
    time.sleep(5)

#print('')


os.chdir('../build/')
make_process = subprocess.Popen("make clean all install FCOMPILER=gfortran-8 UNIT_TEST=1", shell=True, stderr=subprocess.STDOUT)
if make_process.wait() != 0:
    print(bcolors.FAIL + 'GNU FAILED' + bcolors.ENDC)
else:
    print(bcolors.OKGREEN + 'GNU OK' + bcolors.ENDC)
    os.chdir('../data/')
    Vacuum_Process = subprocess.Popen("./Vacuum-MD.out", shell=True, stderr=subprocess.STDOUT)
    Vacuum_Process.wait()
    time.sleep(5)

print('')

os.chdir('../build/')
make_process = subprocess.Popen("make clean all install FCOMPILER=pgfortran UNIT_TEST=1", shell=True, stderr=subprocess.STDOUT)
if make_process.wait() != 0:
    print(bcolors.FAIL + 'PGI FAILED' + bcolors.ENDC)
else:
    print(bcolors.OKGREEN + 'PGI OK' + bcolors.ENDC)
    os.chdir('../data/')
    Vacuum_Process = subprocess.Popen("./Vacuum-MD.out", shell=True, stderr=subprocess.STDOUT)
    Vacuum_Process.wait()
    #time.sleep(10)

print('')
