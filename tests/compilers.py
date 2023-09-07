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

# -------------------------------------------------------------------
# Intel OneAPI SERIAL
""" os.chdir('../build/')
print('Compiling using Intel OneAPI Serial')
make_process = subprocess.Popen(["make", "clean", "all", "CC=icx", "FC=ifx", "TESTING_MODE=1", "OPENMP=no"], shell=False, stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL)
if make_process.wait() != 0:
    print(bcolors.FAIL + 'Intel OneAPI Serial FAILED' + bcolors.ENDC)
else:
    print(bcolors.OKGREEN + 'Intel OneAPI Serial compiled' + bcolors.ENDC)
    os.chdir('../tests/')
    shutil.rmtree('run_intel', ignore_errors=True)
    os.makedirs('run_intel/serial')
    os.chdir('run_intel/serial')
    os.mkdir('pos')
    shutil.copy('../../../build/RUMDEED.out', './RUMDEED.out')
    shutil.copy('../../input', './input')
    Vacuum_Process = subprocess.Popen("./RUMDEED.out", shell=True, stderr=subprocess.STDOUT)
    Vacuum_Process.wait()
    os.chdir('../../')

print('') """
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Intel OneAPI OpenMP
""" os.chdir('../build/')
print('Compiling using Intel OneAPI with OpenMP')
make_process = subprocess.Popen(["make", "clean", "all", "CC=icx", "FC=ifx", "TESTING_MODE=1", "OPENMP=yes"], shell=False, stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL)
if make_process.wait() != 0:
    print(bcolors.FAIL + 'Intel OneAPI OpenMP FAILED' + bcolors.ENDC)
else:
    print(bcolors.OKGREEN + 'Intel OneAPI OpenMP compiled' + bcolors.ENDC)
    os.chdir('../tests/')
    #shutil.rmtree('run_intel', ignore_errors=True)
    os.makedirs('run_intel/openmp')
    os.chdir('run_intel/openmp')
    os.mkdir('pos')
    shutil.copy('../../../build/RUMDEED.out', './RUMDEED.out')
    shutil.copy('../../input', './input')
    Vacuum_Process = subprocess.Popen("./RUMDEED.out", shell=True, stderr=subprocess.STDOUT)
    Vacuum_Process.wait()
    os.chdir('../../')

print('') """
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# GNU SERIAL
os.chdir('../build/')
print('Compiling using GNU Serial')
make_process = subprocess.Popen(["make", "clean", "all", "CC=gcc-12", "FC=gfortran-12", "TESTING_MODE=1", "OPENMP=no"], shell=False, stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL)
if make_process.wait() != 0:
    print(bcolors.FAIL + 'GNU SERIAL FAILED' + bcolors.ENDC)
else:
    print(bcolors.OKGREEN + 'GNU SERIAL compiled' + bcolors.ENDC)
    os.chdir('../tests/')
    shutil.rmtree('run_gnu', ignore_errors=True)
    os.makedirs('run_gnu/serial')
    os.chdir('run_gnu/serial')
    os.mkdir('pos')
    os.mkdir('accel')
    shutil.copy('../../../build/RUMDEED.out', './RUMDEED.out')
    shutil.copy('../../input', './input')
    Vacuum_Process = subprocess.Popen("./RUMDEED.out", shell=True, stderr=subprocess.STDOUT)
    Vacuum_Process.wait()
    os.chdir('../../')

print('')
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# GNU SERIAL
os.chdir('../build/')
print('Compiling using GNU OpenMP')
make_process = subprocess.Popen(["make", "clean", "all", "CC=gcc-12", "FC=gfortran-12", "TESTING_MODE=1", "OPENMP=yes"], shell=False, stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL)
if make_process.wait() != 0:
    print(bcolors.FAIL + 'GNU OpenMP FAILED' + bcolors.ENDC)
else:
    print(bcolors.OKGREEN + 'GNU OpenMP compiled' + bcolors.ENDC)
    os.chdir('../tests/')
    #shutil.rmtree('run_gnu', ignore_errors=True)
    os.makedirs('run_gnu/openmp')
    os.chdir('run_gnu/openmp')
    os.mkdir('pos')
    os.mkdir('accel')
    shutil.copy('../../../build/RUMDEED.out', './RUMDEED.out')
    shutil.copy('../../input', './input')
    Vacuum_Process = subprocess.Popen("./RUMDEED.out", shell=True, stderr=subprocess.STDOUT)
    Vacuum_Process.wait()
    os.chdir('../../')

print('')
# -------------------------------------------------------------------