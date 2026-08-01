# Shared reader module for RUMDEED output files.
# Kristinn Torfason
#
# This module is the single place that encodes the layout of the files
# written by the Fortran code. If an output format changes in src/, update
# the corresponding dtype or column list here and all scripts follow.
#
# Format sources in the Fortran code:
#   position.bin            -> Write_Position in src/mod_pair.F90
#   ramo_current.dt         -> Write_Ramo_Current in src/mod_pair.F90
#   ramo_current.bin        -> Write_Ramo_Current in src/mod_pair.F90
#                              (only written when write_ramo_sec = .true.)
#   density_absorb_top.bin  -> Mark_Particles_Remove in src/mod_pair.F90
#   density_absorb_bot.bin  -> Mark_Particles_Remove in src/mod_pair.F90
#   density_emit.bin        -> Add_Particle in src/mod_pair.F90
#                              (all emitted particles, includes section and species)
#   density_emit_*.bin      -> Add_Particle in src/mod_pair.F90
#                              (one file per species, no section)
#   planes-N.bin            -> Check_Planes in src/mod_verlet.F90
#   init.bin                -> Write_Initial_Variables in src/main.F90

import numpy as np
import pandas as pd
import os.path as path

# ----------------------------------------------------------------
# init.bin: 7 float64 followed by 4 int32 (72 bytes)
INIT_DTYPE = np.dtype([('epsilon_r', np.float64),
                       ('m_eeff', np.float64),
                       ('m_ieff', np.float64),
                       ('length_scale', np.float64),
                       ('time_scale', np.float64),
                       ('vel_scale', np.float64),
                       ('cur_scale', np.float64),
                       ('MAX_PARTICLES', np.int32),
                       ('MAX_EMITTERS', np.int32),
                       ('MAX_SECTIONS', np.int32),
                       ('MAX_LIFE_TIME', np.int32)])

# density_absorb_top.bin and planes-N.bin (52 bytes per record)
# x, y are in units of length_scale, velocities are in m/s
DENSITY_ABSORB_TOP_DTYPE = np.dtype([('x', np.float64),
                                     ('y', np.float64),
                                     ('vx', np.float64),
                                     ('vy', np.float64),
                                     ('vz', np.float64),
                                     ('emit', np.int32),
                                     ('sec', np.int32),
                                     ('id', np.int32)])
PLANES_DTYPE = DENSITY_ABSORB_TOP_DTYPE

# density_absorb_bot.bin (28 bytes per record), x, y in units of length_scale
DENSITY_ABSORB_BOT_DTYPE = np.dtype([('x', np.float64),
                                     ('y', np.float64),
                                     ('emit', np.int32),
                                     ('sec', np.int32),
                                     ('id', np.int32)])

# density_emit.bin, every emitted particle (40 bytes per record)
# x, y, z in units of length_scale; species: 1 = electron, 2 = ion, 3 = atom
DENSITY_EMIT_DTYPE = np.dtype([('x', np.float64),
                               ('y', np.float64),
                               ('z', np.float64),
                               ('emit', np.int32),
                               ('sec', np.int32),
                               ('id', np.int32),
                               ('species', np.int32)])

# density_emit_elec.bin / density_emit_ion.bin / density_emit_atom.bin
# (32 bytes per record), x, y, z in units of length_scale
DENSITY_EMIT_SPECIES_DTYPE = np.dtype([('x', np.float64),
                                       ('y', np.float64),
                                       ('z', np.float64),
                                       ('emit', np.int32),
                                       ('id', np.int32)])

# ramo_current.dt columns, 11 fixed plus one ramo current per species
# (nrSpecies = 3: electrons, ions, atoms)
RAMO_COLUMNS = ['time', 'step', 'current', 'volt', 'nrPart', 'nrElec', 'nrIon',
                'avg_mob', 'avg_part_speed', 'avg_elec_speed', 'avg_ion_speed',
                'ramo_elec', 'ramo_ion', 'ramo_atom']

# position.bin per particle record: x, y, z in meters (not scaled)
POSITION_PARTICLE_DTYPE = np.dtype([('x', np.float64),
                                    ('y', np.float64),
                                    ('z', np.float64),
                                    ('emit', np.int32),
                                    ('sec', np.int32),
                                    ('id', np.int32)])


# ----------------------------------------------------------------
# Read init.bin and return the system parameters as a dictionary
def read_init(filename='out/init.bin'):
    data = np.fromfile(filename, dtype=INIT_DTYPE, count=1)
    return dict(zip(INIT_DTYPE.names, data[0]))


# ----------------------------------------------------------------
# Read ramo_current.dt into a pandas dataframe indexed by time step
def read_ramo_current(filename='out/ramo_current.dt'):
    return pd.read_csv(filepath_or_buffer=filename, index_col=1, sep=r'\s+',
                       header=None, names=RAMO_COLUMNS)


# ----------------------------------------------------------------
# Memory map a binary file of fixed size records and return a dataframe.
# An empty file (e.g. no particles absorbed yet) gives an empty dataframe.
def _read_records(filename, dtype):
    if path.getsize(filename) == 0:
        return pd.DataFrame(np.empty(0, dtype=dtype))
    data_mem = np.memmap(filename, dtype=dtype, mode='r', order='F')
    return pd.DataFrame.from_records(data=data_mem, columns=dtype.names)


# Read density_absorb_top.bin (or a file with the same layout)
def read_density_absorb_top(filename='out/density_absorb_top.bin'):
    return _read_records(filename, DENSITY_ABSORB_TOP_DTYPE)


# Read density_absorb_bot.bin
def read_density_absorb_bot(filename='out/density_absorb_bot.bin'):
    return _read_records(filename, DENSITY_ABSORB_BOT_DTYPE)


# Read density_emit.bin (all emitted particles with section and species)
def read_density_emit(filename='out/density_emit.bin'):
    return _read_records(filename, DENSITY_EMIT_DTYPE)


# Read one of the density_emit_elec/ion/atom.bin files
def read_density_emit_species(filename='out/density_emit_elec.bin'):
    return _read_records(filename, DENSITY_EMIT_SPECIES_DTYPE)


# Read a planes-N.bin file (same record layout as density_absorb_top.bin)
def read_plane(filename):
    return _read_records(filename, PLANES_DTYPE)


# ----------------------------------------------------------------
# Read ramo_current.bin, the Ramo current broken down into sections and
# emitters. Layout of the returned array is [Section, Emitter, time step].
# max_sections and max_emitters come from init.bin, steps is the number of
# lines in ramo_current.dt (see read_ramo_sections_auto below).
def read_ramo_sections(filename, max_sections, max_emitters, steps):
    # The file is only written when write_ramo_sec = .true. in the input
    # namelist. Check the size so we fail with a clear message instead of
    # an obscure memmap error.
    expected_size = max_sections * max_emitters * steps * 8 # float64
    actual_size = path.getsize(filename)
    if actual_size != expected_size:
        raise SystemExit('{:s} is {:d} bytes but expected {:d} '
                         '(MAX_SECTIONS={:d} x MAX_EMITTERS={:d} x steps={:d} x 8). '
                         'Was the simulation run with write_ramo_sec = .true. in the input file?'
                         .format(filename, actual_size, expected_size,
                                 max_sections, max_emitters, steps))

    return np.memmap(filename, dtype=np.float64, mode='r', order='F',
                     shape=(max_sections, max_emitters, steps))


# Convenience version that pulls the dimensions from init.bin and
# ramo_current.dt in the same output directory.
def read_ramo_sections_auto(outdir='out'):
    sys_params = read_init(path.join(outdir, 'init.bin'))
    steps = read_ramo_current(path.join(outdir, 'ramo_current.dt'))['time'].count()
    return read_ramo_sections(path.join(outdir, 'ramo_current.bin'),
                              int(sys_params['MAX_SECTIONS']),
                              int(sys_params['MAX_EMITTERS']), steps)


# ----------------------------------------------------------------
# position.bin readers. The file starts with two int32's, the total number
# of time steps and the output interval N_steps. Then for every recorded
# step: step and nrPart as int32's followed by nrPart particle records.

# Read the header from an open file object
def read_position_header(f):
    max_steps, N_steps = np.fromfile(f, count=2, dtype=np.int32)
    return max_steps, N_steps


# Read one time step from an open file object.
# Returns (step, particles) where particles is a structured array with
# fields x, y, z (meters), emit, sec, id. Returns (None, None) at EOF.
def read_position_step(f):
    header = np.fromfile(f, count=2, dtype=np.int32)
    if header.size < 2:
        return None, None
    step, nrPart = header
    particles = np.fromfile(f, count=nrPart, dtype=POSITION_PARTICLE_DTYPE)
    return step, particles


# Iterate over all time steps in a position.bin file.
# Yields (step, particles) tuples, see read_position_step.
def iter_position(filename='out/position.bin'):
    with open(filename, 'rb') as f:
        read_position_header(f)
        while True:
            step, particles = read_position_step(f)
            if step is None:
                return
            yield step, particles
