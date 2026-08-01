# Hákon Örn Árnason
# Build a statistical summary (means and standard deviations per parameter
# combination) from the dataset compiled by PE_Main_script.py.
#
# Usage: python statistical_set.py <main dataset.csv> <output.csv>

import sys
import numpy as np
import pandas as pd
from scipy.constants import pi, e, hbar, m_e, epsilon_0

if len(sys.argv) != 3:
    raise SystemExit('Usage: python statistical_set.py <main dataset.csv> <output.csv>')

# Read main data set from file
Main_dataset = pd.read_csv(sys.argv[1])
output_csv = sys.argv[2]

# Set pulse width multiplier if used
pulse_width_multiplier = 16

# Create new dataset for statistical work from Main_dataset
statistical_set = pd.DataFrame({
    'Voltage'       : [],
    'Pulse_width'   : [],
    'B_mean'        : [],
    'B_std'         : [],
    'Q_mean'        : [],
    'Q_std'         : [],
    'Amplitude'     : [],
    'norm_tp'       : [],
    'norm_Q'        : [],
    'norm_I'        : [],
    'Ex_mean'       : [],
    'Ex_std'        : [],
    'Ey_mean'       : [],
    'Ey_std'        : [],
    'R_mean'        : [],
    'R_std'         : [],
    'P_mean'        : [],
    'P_std'         : [],
    'Emitter'       : [],
    })
    # <<< Use the below for emitter dimention size comparison >>>
    # emitter_values       = np.unique(Main_dataset['Emitter_dim'].values).tolist()
    # for emit in emitter_values:
    #     emit_set         = voltage_set[ voltage_set['Emitter_dim'] == emit ]
# >>> Use for Checkerboard <<<
#    emitter_values       = np.unique(Main_dataset['WF_Column'].values)
#    for emit in emitter_values:
#        emit_set         = voltage_set[ voltage_set['WF_Column'] == emit ]


VOLT_values              = np.unique(Main_dataset['Voltage'].values)
for volt in VOLT_values:
    voltage_set          = Main_dataset[Main_dataset['Voltage']    == volt ]
    emitter_values       = np.unique(Main_dataset['Emitter_dim'].values).tolist()
    for emit in emitter_values:
        emit_set         = voltage_set[ voltage_set['Emitter_dim'] == emit ]
        AMP_values       = np.unique(emit_set['Amplitude'].values)
        for amp in AMP_values:
            amp_set      = emit_set[emit_set['Amplitude']          == amp  ]
            PW_values    = np.unique(amp_set['Pulse_width'].values)
            for pw in PW_values:
                pw_set   = amp_set[amp_set['Pulse_width']          == pw   ]
                norm_tp  = pw * np.sqrt( (volt * e ) / (2 * m_e) ) * pulse_width_multiplier * 1E-10 # <-- 2.5E-4 / 2500E-9 * 1E-12
                # sqrt(qV/2m) * sigma * FWHM (pwm) * delta_t * 1 / D
                mean_B   = pw_set['Brightness'  ].mean()
                std_B    = pw_set['Brightness'  ].std()
                mean_Q   = pw_set['Total_charge'].mean()
                std_Q    = pw_set['Total_charge'].std()
                mean_R   = pw_set['Max_ramo'    ].mean()
                std_R    = pw_set['Max_ramo'    ].std()
                mean_Ex  = pw_set['Emittance_x' ].mean()
                std_Ex   = pw_set['Emittance_x' ].std()
                mean_Ey  = pw_set['Emittance_y' ].mean()
                std_Ey   = pw_set['Emittance_y' ].std()
                mean_P   = pw_set['Pulse_ratio' ].mean()
                std_P    = pw_set['Pulse_ratio' ].std()
                # Normalized charge for a square !!! NEED TO ADD CIRCLE !!!
                norm_Q   = (mean_Q * 2500E-9 )/ (500E-9 * 500E-9 * epsilon_0 * volt)
                # Normalized current for a square
                #norm_I   = np.sqrt(( 2 * m_e ) / e) * (( mean_R * (2500E-9 * 2500E-9 )) / ( epsilon_0 * 500E-9 * 500E-9 * volt**(3/2) ))
                # For a circle - normalized current - NEED TO FIX
                norm_I   = np.sqrt(( 2 * m_e ) / e) * (( mean_R * (2500E-9 * 2500E-9 )) / ( epsilon_0 * ((250E-9)**2) * volt**(3/2) ))
                # wrkfnct  = np.unique(emit_set['WF_Column'].values) # Workfunction values if based on checkerboard
                emitter  = np.unique(emit_set['Emitter_dim'].values) # Emitter_dim
                #P_start  = pw_set['Outpulse_start']
                add_series_set = pd.DataFrame({
                    'Voltage'       : [volt     ],
                    'Pulse_width'   : [pw       ],
                    'B_mean'        : [mean_B   ],
                    'B_std'         : [std_B    ],
                    'Q_mean'        : [mean_Q   ],
                    'Q_std'         : [std_Q    ],
                    'Amplitude'     : [amp      ],
                    'norm_tp'       : [norm_tp  ],
                    'norm_Q'        : [norm_Q   ],
                    'norm_I'        : [norm_I   ],
                    'Ex_mean'       : [mean_Ex  ],
                    'Ex_std'        : [std_Ex   ],
                    'Ey_mean'       : [mean_Ey  ],
                    'Ey_std'        : [std_Ey   ],
                    'R_mean'        : [mean_R   ],
                    'R_std'         : [std_R    ],
                    'P_mean'        : [mean_P   ],
                    'P_std'         : [std_P    ],
                    'Emitter'       : [emitter  ],
                    #'WF_col'        : [wrkfnct  ],
                    #'P_start'       : [P_start  ],
                    })
                # Add calculated set to statistical set
                statistical_set = pd.concat([statistical_set, add_series_set])

# To save statistical set as a separate file
statistical_set.to_csv(output_csv)
print('Wrote', output_csv)
