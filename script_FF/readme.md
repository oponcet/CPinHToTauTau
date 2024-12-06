fake_factor_derivation/
│
├── inputs/
│   ├── inputs_pckl_region*.json # containts the path to the *.pck files
│
├── outputs/                    # Store the output JSON and ROOT files
│   ├── region # ex pi_pi_2jets
│       ├── cat_plots_region.root/A/TH1D(pt1).... #TH1D, TH2D of the category ABCD before any correction, will serve as inputs for next steps
│       ├── FF_pt1_region.json # json file containing Fake Factor 
│       ├── CR_QCD_region.json # json file containing the QCD closure correction
│       ├── ER_dR.json # json file containing the dR extrapolation correction
│       ├── extrapolation_plots/ 
│       ├── final_control_plots/
│
├── src/
│   ├── __init__.py
│   ├── input_processing.py       # Script to convert all needed pck file in ROOT files by region
│   ├── derive_FF.py              # Script to caclulate Fake factors and error (FF and FF0). Call fit_function.py and save results in json.
│   ├── fit_functions.py          # Fitting functions for FF, CR, ER, etc.
│   ├── apply_FF_control.py       # Function to apply fake factor and create control plots, especially QCD.     
│   ├── derive_corr_CR.py         # Function to derive CR correction from QCD control plot. Call fit_function.py and save results in json.
│   ├── apply_FF_CR_control.py    # Function to apply fake factor and QCD CR correction and create control plots.
│   ├── apply_CR_ABCD0.py         # Function to apply QCD CR correction and create plots. 
│   ├── derive_corr_ER.py         # Function to derive ER correction from dR control plot. Call fit_function.py and save results in json.
│   ├── apply_FF_CR_ER_control.py # Function to apply fake factor and QCD CR correction and create control plots.
│   ├── utils.py                  # Utility functions (e.g., for JSON saving, data manipulation)
│
│
├── README.md
└── requirements.txt



+---------------------------+
|         Main              |
|---------------------------|
| + main()                 |
|---------------------------|
| Calls various processing |
| functions in sequence    |
+---------------------------+
          |
          v
+---------------------------+
|     input_processing.py   |
|---------------------------|
| + convert_pck_to_root()  |
|---------------------------|
| Converts pck to ROOT,    |
| generates input histos   |
+---------------------------+
          |
          v
+---------------------------+
|       derive_FF.py        |
|---------------------------|
| + derive_fake_factor()   |
| + save_ff_to_json()      |
|---------------------------|
| Fits FF from A/B         |
| and saves to JSON.       |
+---------------------------+
          |
          v
+---------------------------+
|     apply_FF_control.py   |
|---------------------------|
| + apply_ff_to_histos()   |
| + generate_control_plots()|
|---------------------------|
| Applies FF and generates |
|     control plots.       |
+---------------------------+
          |
          v
+---------------------------+
|     derive_corr_CR.py     |
|---------------------------|
| + derive_cr_correction() |
| + save_cr_to_json()      |
|---------------------------|
| Fits CR corrections and  |
| saves to JSON.           |
+---------------------------+
          |
          v
+---------------------------+
| apply_FF_CR_control.py    |
|---------------------------|
| + apply_ff_cr_to_histos()|
| + generate_cr_plots()    |
|---------------------------|
| Applies FF + CR to ABCD  |
| and     ABCD0   and      |
| generates control plots. |
+---------------------------+
          |
          v
+---------------------------+
|       derive_FF.py        |
|---------------------------|
| + derive_fake_factor()   |
| + save_ff_to_json()      |
|---------------------------|
| Fits FF0 from A0/B0 (CR)  |
| and saves to JSON.       |
+---------------------------+
          |
          v
+---------------------------+
|     apply_FF_control.py   |
|---------------------------|
| + apply_ff_to_histos()   |
| + generate_control_plots()|
|---------------------------|
| Applies FF0 and generates|
|     control plots.       |
+---------------------------+
          |
          v
+---------------------------+
|     derive_corr_ER.py     |
|---------------------------|
| + derive_er_correction() |
| + save_er_to_json()      |
|---------------------------|
| Fits ER corrections and  |
| saves to JSON.           |
+---------------------------+
          |
          v
+---------------------------+
| apply_FF_CR_ER_control.py |
|---------------------------|
| + apply_all_corrections()|
| + generate_final_plots() |
|---------------------------|
| Applies FF, CR, ER, and  |
| generates final control  |
| plots to ABCD and ABCD0. |
+---------------------------+

+---------------------------+
|       fit_functions.py     |
|---------------------------|
| + fit_pol_landau()        |
| + fit_cr_function()       |
| + fit_er_function()       |
|---------------------------|
| Provides reusable fitting |
| methods for FF, CR, ER.   |
+---------------------------+

+---------------------------+
|         utils.py           |
|---------------------------|
| + save_to_json()          |
| + load_from_json()        |
| + manipulate_histos()     |
|---------------------------|
| Common utilities for JSON |
| handling, data operations |
+---------------------------+
