# Hákon Örn Árnason
# 30.10.2022
# Compile statistical data (emittance, brightness, pulse shape) from a large
# number of photo emission runs into one CSV file.
#
# Usage: python PE_Main_script.py <folder with run directories> <output.csv>
# Each run directory must contain the input, work and laser files along with
# the out/ directory from the simulation.

import sys
import numpy as np
import pandas as pd
import glob
import os
import os.path as path
import f90nml
from scipy.constants import pi, e, hbar, m_e, epsilon_0

# The file layouts live in scripts/python_package
sys.path.insert(0, path.join(path.dirname(path.abspath(__file__)), '..', 'python_package'))
import Vacuum
import rumdeed_io

if len(sys.argv) != 3:
    raise SystemExit('Usage: python PE_Main_script.py <folder with run directories> <output.csv>')
folderpath = sys.argv[1]
output_csv = sys.argv[2]

# Define external calculation parameters
pulse_width_multiplier = 16          # Input pulse length multiplier (int)

# Create Dataset for input variables and calculated outputs
Main_dataset = pd.DataFrame()

for filepath in sorted(glob.glob(path.join(folderpath, '*'))):
    if not path.isdir(filepath):
        continue
    print(filepath)
    # Load in RUMDEED input file
    input_file          = f90nml.read(path.join(filepath, 'input'))
    # Aquire Input voltage from input file
    input_voltage       = input_file['input']['v_s']                # Voltage (float)
    box_dimensions      = input_file['input']['box_dim']            # Dimensions of diode (list[float, float, float])
    time_step           = input_file['input']['time_step']          # Time step of simulation (float)
    emission_mode       = input_file['input']['emission_mode']      # Type of emission (int)
    emitter_dim         = input_file['input']['emitters_dim']       # Emitter dimentions (list[float, float, float])
    emitter_type        = input_file['input']['emitters_type']      # Emitter type (int)

    # Open the workfunction input file (the file the code reads is called 'work')
    with open(path.join(filepath, 'work')) as wf_file:
        array               = [[float(x) for x in line.split()] for line in wf_file]
        #wf_decision_operator  = int(array[0][0]) # Program load operator?
        wf_matrix_column    = int(array[1][0])    # Matrix number of columns
        wf_matrix_row       = int(array[1][1])    # Matrix number of rows
        wf_emitter          = array[2][0]         # Emitter workfunction
        wf_block            = array[3][0]         # Blocked emitter workfunction

    # Open the laser input file (Photoemission specific)
    with open(path.join(filepath, 'laser')) as laser_file:
        array               = [[float(x) for x in line.split()] for line in laser_file]
        #decision_operators  = int(array[0][0]), int(array[0][1]), int(array[0][2]) # Program load operators ([1, 1, 2])
        laser_energy        = array[1][0]       # Energy of the laser in eV (float)
        laser_variation     = array[2][0]       # Variation of laser energy in eV (float)
        pulse_center        = int(array[3][0])  # Center of pulse in steps (int)
        pulse_width         = int(array[4][0])  # Width of the pulse in steps (Gaussian) (int)
        pulse_amplitude     = array[5][0]       # Amplitude of the pulse (float)

    # Load absorbed electrons file
    filename_absorbed_top = path.join(filepath, 'out/absorbed_top.dt')
    data                  = np.loadtxt(filename_absorbed_top) # Read absorbed at anode file
    absor_elec            = data[:, 3] # Data column for absorbed electrons at the anode
    tot_elec_top          = 0          # Temp var for total electron count in system (float64)
    max_elec_top          = 0          # Temp var. for e count to locate max (float64)
    out_pulse_start       = np.nan     # Stays NaN if no electrons were absorbed
    out_pulse_end         = np.nan

    # Count the absorbed electrons at anode
    for k in range(len(absor_elec)):
        tot_elec_top = tot_elec_top + absor_elec[k]
    tot_charge_top = tot_elec_top * e
    # Check when anode pulse starts (first electron exits)
    for g in range(len(data[:, 2])):
        if (data[g, 2] != 0):
            out_pulse_start = data[g, 0]
            break
    # Check when anode pulse ends (last electron exits)
    for k in range(len(absor_elec)):
        max_elec_top = max_elec_top + absor_elec[k]
        if (tot_elec_top == max_elec_top):
            out_pulse_end = data[k, 0]
            break

    # Calculate the length of the outgoing (anode) pulse with regards to incoming (laser) pulse. (float64)
    out_pulse_ratio = (out_pulse_end - out_pulse_start) / (pulse_width * pulse_width_multiplier)

    # Load emitted electron file
    filename_emitted = path.join(filepath, 'out/emitted.dt')
    data                  = np.loadtxt(filename_emitted) # Read emitted electrons file
    emitted_elec          = data[:, 3] # Data column for emitted electrons at the cathode
    tot_elec_emit         = 0          # Temp var for total electron count in system (float64)
    max_elec_emit         = 0          # Temp var. for e count to locate max (float64)
    in_pulse_start        = np.nan
    in_pulse_end          = np.nan

    # Count the emitted electrons at cathode
    for k in range(len(emitted_elec)):
        tot_elec_emit = tot_elec_emit + emitted_elec[k]
    tot_charge_emit = tot_elec_emit * e
    # Check when cathode pulse starts (first electron emitted)
    for g in range(len(data[:, 2])):
        if (data[g, 2] != 0):
            in_pulse_start = data[g, 0]
            break
    # Check when cathode pulse ends (last electron emitted)
    for k in range(len(emitted_elec)):
        max_elec_emit = max_elec_emit + emitted_elec[k]
        if (tot_elec_emit == max_elec_emit):
            in_pulse_end = data[k, 0]
            break


    # Density file
    filename = path.join(filepath, 'out/density_absorb_top.bin')

    try:
        # Current format (particle id added in 2022)
        df = rumdeed_io.read_density_absorb_top(filename)
    except ValueError:
        # Pre-2022 version without the id field
        dt_old  = np.dtype([('x' , np.float64), ('y', np.float64) , ('vx', np.float64), ('vy', np.float64),
                            ('vz', np.float64), ('emit', np.int32), ('sec', np.int32)                     ])
        data_mem = np.memmap(filename, dtype=dt_old, mode='r', order='F')
        df       = pd.DataFrame.from_records(data=data_mem, columns=data_mem.dtype.names)

    df['v']     = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2)

    df["x'"]    = (df['vx']/df['vz'])/1.0E-3
    df["y'"]    = (df['vy']/df['vz'])/1.0E-3

    df["r"]     = np.sqrt(df["x"]**2 + df["y"]**2)
    df["r'"]    = (np.sqrt(df['vx']**2 + df['vy']**2) / df['vz']) / 1.0E-3

    e_x, sigma_w, sigma_wp, theta_ell = Vacuum.Calc_Emittance(df, "x", "x'")
    e_y,       _,        _,         _ = Vacuum.Calc_Emittance(df, "y", "y'")

    # Read data for current
    filename_ramo   = path.join(filepath, 'out/ramo_current.dt') # Ramo current
    df_cur          = rumdeed_io.read_ramo_current(filename_ramo)

    cur             = df_cur['current'].max() # Find the max Ramo current

    # Calculate brightness
    brightness      = 2 * cur / (np.pi**2 * e_x * e_y * 1.0E-12)

    # Write this simulation run into dataframe
    temp_set        = pd.DataFrame({'Voltage'            : [input_voltage],
                                    'Box_dim'            : [box_dimensions],
                                    'Time_step'          : [time_step],
                                    'Emission_mode'      : [emission_mode],
                                    'Emitter_dim'        : [emitter_dim],
                                    'Emitter_type'       : [emitter_type],
                                    'Pulse_width'        : [pulse_width],
                                    'Amplitude'          : [pulse_amplitude],
                                    'Laser_energy'       : [laser_energy],
                                    'Laser_var'          : [laser_variation],
                                    'Pulse_center'       : [pulse_center],
                                    'Brightness'         : [brightness],
                                    'Emittance_x'        : [e_x],
                                    'Emittance_y'        : [e_y],
                                    'Pulse_ratio'        : [out_pulse_ratio],
                                    'Max_ramo'           : [cur],
                                    'Total_charge'       : [tot_charge_top],
                                    'Total_charge_emit'  : [tot_charge_emit],
                                    'Inpulse_start'      : [in_pulse_start],
                                    'Inpulse_end'        : [in_pulse_end],
                                    'Outpulse_start'     : [out_pulse_start],
                                    'Outpulse_end'       : [out_pulse_end],
                                    'WF_Column'          : [wf_matrix_column],
                                    'WF_Row'             : [wf_matrix_row],
                                    'WF_Emitter'         : [wf_emitter],
                                    'WF_Blocker'         : [wf_block],
                                  })

    # Write the simulation run into the Main dataset
    Main_dataset    = pd.concat([temp_set, Main_dataset])

Main_dataset.to_csv(output_csv)
print('Wrote', output_csv)
