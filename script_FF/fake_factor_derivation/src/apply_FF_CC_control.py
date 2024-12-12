'''
Author : @oponcet
Date : 12-12-2024
Script to apply pt dependent Fake Factors and met_var_qcd_h1 closure correction on pt 
distribution of the cat C to get D category.

Input:
    - FF root file (TH1D)
    - CC root file (TH1D)
    - 2D histogram of C met_var_qcd_h1, pt (TH2D)
    - 1D histogram of D pt (THStack) to build the final D histogram
                         
'''

import os
import ROOT

def apply_FF_CC_C(input_file, catC, catD, dm, njet)

    output_root_file_path = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF_CC/{dm}/{dm}_{njet}_D.root'
    output_root_file_png = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF_CC/{dm}/{dm}_{njet}_D.png'
    # Ensure directories exist if not create them
    os.makedirs(os.path.dirname(output_root_file_path), exist_ok=True)
    os.makedirs(os.path.dirname(output_root_file_png), exist_ok=True)

    # Load input ROOT files and histograms 
    





