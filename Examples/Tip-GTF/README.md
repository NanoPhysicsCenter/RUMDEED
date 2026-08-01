# Tip GTF (general thermal-field) emission

Emission from a hyperboloid tip (`EMISSION_MODE = 3`, `mod_emission_tip.f90`)
using the General Thermal-Field (GTF) current density of Jensen's oGTF model
(`mod_kevin_rjgtf_v2.f90`), selected with `EMITTERS_TYPE(1) = 4`.

The GTF model covers the full range from pure field emission to pure
thermionic emission, so the run can be moved between regimes with `T_temp`
and `V_S` alone.

## Geometry

`EMITTERS_DIM(1:3, 1) = d_tip, R_base, h_tip` in nano-meters:

- `d_tip` — distance from the tip apex to the anode,
- `R_base` — radius of the tip at its base,
- `h_tip` — height of the tip.

The gap spacing is `d_tip + h_tip`; `BOX_DIM` is not used for the geometry.

## Notes

- The work function for the tip is hard coded to 4.7 eV (`w_theta` in
  `mod_emission_tip.f90`); the `work` file is not read in this mode.
- The cathode temperature is `T_temp`.
- The number of emitted electrons per step comes from integrating the GTF
  current density over the tip surface with the CUBA library; the
  `cuba_method`, `cuba_epsabs`, `cuba_epsrel`, `cuba_mineval` and
  `cuba_maxeval` namelist variables tune this integration.

## Run

```sh
../../build/RUMDEED.out
```
