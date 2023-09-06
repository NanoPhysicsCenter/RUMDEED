#!/usr/bin/python3

import os
import subprocess
import time
import shutil

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# Intel New
os.chdir('../build/')
print('Compiling with Intel OneAPI')
make_process = subprocess.Popen(["make", "clean", "all", "CC=icx", "FC=ifx", "UNIT_TEST=1"], shell=False, stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL)
if make_process.wait() != 0:
    print(bcolors.FAIL + 'Intel Classic FAILED' + bcolors.ENDC)
else:
    print(bcolors.OKGREEN + 'Intel Classic OK' + bcolors.ENDC)
    os.chdir('../tests/')
    shutil.rmtree('run_intel', ignore_errors=True)
    os.mkdir('run_intel')
    os.chdir('run_intel')
    os.mkdir('pos')
    shutil.copy('../../build/RUMDEED.out', './RUMDEED.out')
    shutil.copy('../input', './input')
    Vacuum_Process = subprocess.Popen("./RUMDEED.out", shell=True, stderr=subprocess.STDOUT)
    Vacuum_Process.wait()
    os.chdir('../')
    time.sleep(5)

print('')
##exit()

os.chdir('../build/')
print('Compiling with GNU')
make_process = subprocess.Popen(["make", "clean", "all", "CC=gcc", "FC=gfortran", "UNIT_TEST=1"], shell=False, stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL)
if make_process.wait() != 0:
    print(bcolors.FAIL + 'GNU FAILED' + bcolors.ENDC)
else:
    print(bcolors.OKGREEN + 'GNU OK' + bcolors.ENDC)
    os.chdir('../tests/')
    shutil.rmtree('run_gnu', ignore_errors=True)
    os.mkdir('run_gnu')
    os.chdir('run_gnu')
    os.mkdir('pos')
    os.mkdir('accel')
    shutil.copy('../../build/RUMDEED.out', './RUMDEED.out')
    shutil.copy('../input', './input')
    Vacuum_Process = subprocess.Popen("./RUMDEED.out", shell=True, stderr=subprocess.STDOUT)
    Vacuum_Process.wait()
    os.chdir('../')
    time.sleep(5)

print('')

""" os.chdir('../build/')
make_process = subprocess.Popen("make clean all install CC=gcc FC=pgfortran UNIT_TEST=1", shell=True, stderr=subprocess.STDOUT)
if make_process.wait() != 0:
    print(bcolors.FAIL + 'PGI FAILED' + bcolors.ENDC)
else:
    print(bcolors.OKGREEN + 'PGI OK' + bcolors.ENDC)
    os.chdir('../data/')
    Vacuum_Process = subprocess.Popen("./RUMDEED.out", shell=True, stderr=subprocess.STDOUT)
    Vacuum_Process.wait()
    #time.sleep(10)

print('') """
