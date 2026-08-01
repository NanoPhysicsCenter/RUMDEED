# RUMDEED examples

Each folder holds a ready-to-run `input` file (Fortran namelist, see
`namelist /input/` in `src/mod_global.F90` for all available variables)
plus whatever auxiliary files that emission mode reads. Run an example by
executing `RUMDEED.out` inside its folder; output is written to `out/`
(created automatically).

```sh
cd Examples/Planar-FE/2.0eV
../../../build/RUMDEED.out
```

| Example | `EMISSION_MODE` | Description | Extra files read |
|---|---|---|---|
| `Planar-FE/2.0eV` | 10 | Planar field emission, 2.0 eV work function | `work` |
| `Planar-FE/2.5eV` | 10 | Planar field emission, 2.5 eV work function | `work` |
| `Planar-FE/4.7eV` | 10 | Planar field emission, 4.7 eV work function | `work` |
| `Ion` | 10 | Planar field emission with electron-N2 collisions (`COLLISION_MODE = 2`) | `work`, `N2-tot-cross.txt`, `N2-ion-cross.txt` |
| `Photo` | 1 | Planar photo emission | `work`, `laser` |
| `Checkerboard-TFE` | 9 | Thermal-field (GTF) emission from a 4x4 checkerboard of 2.0/2.5 eV patches | `work` |
| `Tip-FE` | 3 | Field emission from a hyperboloid tip (`EMITTERS_TYPE = 1`) | — |
| `Tip-GTF` | 3 | Thermal-field (GTF) emission from a hyperboloid tip (`EMITTERS_TYPE = 4`) | — |
| `Torus` | 12 | Field emission from a torus / looped CNT, field interpolated from mesh data | `Torus_mesh_data.txt`, `Torus_mesh_data_1V.txt`, `image_charge_data_torus.txt` (symlinks into `data/Looped-CNT/`) |
| `GPU-Planar-FE` | 10 | Planar field emission on an NVIDIA GPU (OpenACC build, `MH_BATCH = .True.`) | `work` |

## The `work` file

Planar emission modes (1, 9, 10) read the position dependent work function
from a file named `work` (`mod_work_function.F90`). First line is the model
type: 1 = checkerboard, 2 = Gaussians, 3 = circles, 4 = Voronoi. For the
checkerboard the next line is the grid size (`ny nx`) followed by the work
function values in eV; a 1x1 grid gives a uniform work function. The tip
(mode 3) and torus (mode 12) modes use hard coded work functions instead
(4.7 eV and 4.65 eV).

## Collisions

`COLLISION_MODE` 1-4 requires the N2 cross section data files
`N2-tot-cross.txt` and `N2-ion-cross.txt` in the run directory (see the
`Ion` example; the data files are also in `ion-data/`). Modes 3 and 4 place
neutral atoms in a cylinder of radius `EMITTERS_DIM(1,1)` at startup — at
NTP density this exceeds `MAX_PARTICLES` for micron-sized systems, so lower
`P_abs` accordingly.

## Running on GPUs

The N-body and emission-field kernels can be offloaded to NVIDIA GPUs with
OpenACC: build with `make FC=nvfortran ACC=gpu` (add `GPU_CC=cc80` etc. when
compiling on a machine without a visible GPU) and set `MH_BATCH = .True.`
in mode-10 decks. The planar modes (1, 9, 10) and the hyperboloid tip
(mode 3) run on the GPU; the cylindrical tip (11) and torus (12) fall back
to the OpenMP CPU path automatically. See `GPU-Planar-FE/README.md` for
details, including a Slurm job script example.

## Notes

- `NREMIT` is capped by `MAX_EMITTERS` in `src/mod_global.F90` (currently 1);
  increase it there and rebuild for multi-emitter runs
  (`Planar-FE/4.7eV multi`).
- Old input files that use `COLLISIONS = .True./.False.` no longer parse;
  the variable was replaced by the integer `COLLISION_MODE`.
