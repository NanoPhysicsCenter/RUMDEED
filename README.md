# Molecular Dynamics code
    
    Photoemission Branch
    + Workfunction pattern
      + checkerboard
      + random generated pattern - WIP
      + Emission tip
      + Circle - WIP
    + Velocity distibution is Gaussian
      + Velocity has Newtonian calculations
    + Output pulse can be Gaussian
      + Quantum efficiency is controlled by amplitude modulation
      + Pulse repetition is possible - WIP
    + Input file for laser and pulse - WIP
    
    See master branch for more info


## Input file
```
&
INPUT
  V_S = 1.0d1,
  BOX_DIM = 0.0d0, 0.0d0, 2500.0d0,
  TIME_STEP = 0.25d-3,
  STEPS = 60000,
  EMISSION_MODE = 1,
  NREMIT = 1,
  IMAGE_CHARGE = .TRUE.,

  EMITTERS_DIM(1:3, 1) = 500.0d0, 500.0d0, 0.0d0,
  EMITTERS_POS(1:3, 1) = -250.0d0, -250.0d0, 0.0d0,
  EMITTERS_TYPE(1) = 2,
  EMITTERS_DELAY(1) = 0,

/
```

## laser file
    1 2
    4.7 0.02
    10000 1000 5

### Input Warning
  The Gauss_Emission needs to be enabled if such input/output is desired, "1" in the laser file, "2" is instant start and continous.
  The second number picks velocity profile, "1" being zero initial velocity while "2" will give inital velocity dependant on workfunction.
  Second line is laser (photon) energy and variation, with mean controlling the energy level in electronVolts (eV) and std being standard deviation of the laser (in eV's as well). 
  This is normal distribution with Box-Muller method.
  This energy is compared to the work function with the excess making way for Newtonian inital velocity given to the electrons.
  Third line is gauss pulse parameters, center (mu), width (sigma) and A(mplitude) of the pulse. 
  
  The gaussian pulse is simulated with output restriction of electrons according to normal distribution. This should in theory simulate the Quantum Efficiency and Intensity via amplitude modulation.

  
      

## Notes
  Needs fixing?
  See [notes.pdf](doc/notes.pdf)?

## Acknowledgments

* Rannís grant nr. 174512-051
* Rannís grant nr. 218029-051
* Rannís grant nr. 174127-052
