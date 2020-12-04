# Molecular Dynamics code
    
    This branch is focused on photoemission.
    + Across a large vacuum gap (500 micro meters)
    + Differing workfunction pattern across the emitter
      + checkerboard patter
      + random generated pattern 
    
    See master branch for more info


## Input file
```
&
INPUT
  V_S = 1.0d1,
  BOX_DIM = 200.0d0, 200.0d0, 500000.0d0,
  TIME_STEP = 0.25d-3,
  STEPS = 20000,
  EMISSION_MODE = 1,
  NREMIT = 1,

  EMITTERS_DIM(1:3, 1) = 1000.0d0, 1000.0d0, 1000.0d0,
  EMITTERS_POS(1:3, 1) = 100.0d0, 100.0d0, 100.0d0,
  EMITTERS_TYPE(1) = 1,
  EMITTERS_DELAY(1) = 0,

/
```

## Notes
  Needs fixing?
  See [notes.pdf](doc/notes.pdf)?

## Photo-Emission
    Photo
      Add detail about function

## Field-Emission
    Field
      Where is this located?
## Acknowledgments

* Rannís grant nr. ?
* Rannís grant nr. 174127-052
