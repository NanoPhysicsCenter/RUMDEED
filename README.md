# Molecular Dyanmics code
  Hello

## Input file
```
&INPUT                                                                                                                                                                                                                                                                         
  V_S = 2.0d3,                                                                                                                                                                                                                                                                 
  BOX_DIM = 0.0d0, 0.0d0, 1000.0d0,                                                                                                                                                                                                                                            
  TIME_STEP = 0.25d-3,                                                                                                                                                                                                                                                         
  STEPS = 1000,                                                                                                                                                                                                                                                                
  NREMIT = 1,                                                                                                                                                                                                                                                                  
  EMITTERS_DIM(1:3, 1) = 2500.0d0, 2500.0d0, 0.0d0,                                                                                                                                                                                                                            
  EMITTERS_POS(1:3, 1) = 0.0d0, 0.0d0, 0.0d0,                                                                                                                                                                                                                                  
  EMITTERS_TYPE(1) = 2,                                                                                                                                                                                                                                                        
  EMITTERS_DELAY(1) = 0,                                                                                                                                                                                                                                                       
  EMISSION_MODE = 2,                                                                                                                                                                                                                                                           
/ 
```

## Notes
  See [notes.pdf](Vacuum-MD/doc/notes.pdf)

## Photo-Emission
    Photo

## Field-Emission
    Field

## Acknowledgments

* Rannís grant nr. 174127-052
