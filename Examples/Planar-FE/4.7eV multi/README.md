# Planar FE, 4.7 eV, multiple emitters

This example is a placeholder: the current build caps the number of emitters
at `MAX_EMITTERS = 1` in `src/mod_global.F90`. To run a multi-emitter
simulation, increase `MAX_EMITTERS`, rebuild, and use an input file like the
one in `../4.7eV/` with `NREMIT > 1` and one `EMITTERS_DIM` / `EMITTERS_POS`
/ `EMITTERS_TYPE` / `EMITTERS_DELAY` block per emitter.
