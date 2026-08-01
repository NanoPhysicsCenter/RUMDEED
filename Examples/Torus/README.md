# Torus (looped CNT) field emission

Field emission from a torus shaped emitter (a looped carbon nanotube standing
on the cathode plane), `EMISSION_MODE = 12` (`mod_torus.F90`).

## Required data files

The electric field is not computed analytically. It is interpolated (inverse
distance weighting over the 4 nearest neighbours of a kd-tree) from field maps
solved externally (COMSOL) and read from the run directory:

- `Torus_mesh_data.txt` — mesh points and field at the operating voltage.
  Symlink to `data/Looped-CNT/Torus_mesh_data-160kv.txt` here. Other voltages
  (110, 180, 200 and 500 kV) are available in `data/Looped-CNT/`; if you switch
  the file, set `V_S` in `input` to match.
- `Torus_mesh_data_1V.txt` — field map at 1 V, used for the Ramo current.
- `image_charge_data_torus.txt` — tabulated image charge positions/magnitudes.

File format: 6 columns per line, `x y z Ex Ey Ez`, coordinates in mm (converted
to m on read). The maximum z coordinate of the mesh must equal `BOX_DIM(3)`
(here 0.4 mm = 4.0e5 nm) or the run aborts at startup.

## Notes

- The work function is hard coded to 4.65 eV (CNT) in `mod_torus.F90`
  (`w_theta_pos_tip`); the `work` file is not read in this mode.
- The `EMITTERS_*` variables are ignored; the emission area is a patch at the
  top of the torus (`rad_patch` in `mod_torus.F90`).
- At startup the module writes the surface field along/around the top of the
  torus (`Calc_E_Along_Top` / `Calc_E_Around_Top`) which is useful to verify
  the mesh data was read correctly.

## Run

```sh
../../build/RUMDEED.out
```
