
'''
Author : @oponcet
Date : 10-12-2024
Main script calling the function of Fake Factor derivation, closure and extrapolation 
correction, control plots ... 

CPPYY_BACKEND_LIBRARY=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/lib/libcppyy_backend3_9.so
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/bin/thisroot.sh 


Input: 
    - inputs/inputs_json/region.json (json file with the region configuration)

Output: see called function
'''

import os
import json
import ROOT
import uproot

from derive_FF import calculate_fake_factor
from apply_FF import *
from apply_FF_CCorr import *
from plotting_hep import plot_stack_from_root_file




def main(config_path, dm):

    # ROOT configuration
    ROOT.gROOT.SetBatch(True) # Do not display canvas

    # Load the JSON configuration
    with open(config_path, 'r') as f:
        config = json.load(f)

    categories = config['categories']
    variable = "hcand_1_pt"

    catA, catB = None, None  # Initialize to None

    if dm == 'pi_1':
        njet = 'allnjet'
        catA = 'tautau__real_1__hadA__allnjet__pi_1'
        catB = 'tautau__real_1__hadB__allnjet__pi_1'

        catA0 = 'tautau__real_1__hadA0__allnjet__pi_1'
        catB0 = 'tautau__real_1__hadB0__allnjet__pi_1'

        input_file_region = f'script_FF/fake_factor_derivation/inputs/inputs_rootfile/{dm}/{dm}_allnjet.root'   # script_FF/fake_factor_derivation/inputs/inputs_rootfile/pi_1/pi_1_has_0j.root

        calculate_fake_factor(input_file_region, catA, catB, dm, njet)
        # calculate_fake_factor(input_file_region, catA0, catB0, dm, njet)


        # apply_fake_factor_CD(input_file_region, catA, catB, dm, njet)
    else:
        for njet in config['categories'][dm]:
            # Loop over ABCD categories
            for cat in config['categories'][dm][njet]['ABCD']:
                if '_hadA__' in cat:
                    catA = cat
                if '_hadB__' in cat:
                    catB = cat
                if '_hadC__' in cat:
                    catC = cat
                if '_hadD__' in cat:
                    catD = cat
            for cat in config['categories'][dm][njet]['A0B0C0D0']:
                if '_hadA0__' in cat:
                    catA0 = cat
                if '_hadB0__' in cat:
                    catB0 = cat
                if '_hadC0__' in cat:
                    catC0 = cat
                if '_hadD0__' in cat:
                    catD0 = cat
                
                
            if not catA or not catB:
                raise ValueError(f"Categories _hadA__ or _hadB__ not found for {dm} in njet {njet}")   
            
            # Define input file
            
            input_file_region = f'script_FF/fake_factor_derivation/inputs/inputs_rootfile/{dm}/{dm}_{njet}.root'   # script_FF/fake_factor_derivation/inputs/inputs_rootfile/pi_1/pi_1_has_0j.root

        
        #########################################
        ###        PLOT THE HCAND_1_PT        ###
        #########################################

        # Replace the following paths with actual file paths and names
        # variable = "hcand_1_pt"
        # canvas_name = f"{catD}/{variable}/canvas"
        # mc_stack = canvas_name["mc_stack"]  # THStack
        # data_hist = canvas_name["data_hist"]  # TH1D


        # output_dir = f"script_FF/fake_factor_derivation/outputs/outputs_applyFF_CCorr/dm/pi_1"  # output/plots

        # # # Ensure directories exist if not create them
        # # os.makedirs(output_dir, exist_ok=True)

        # root_file_to_plot = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF_CCorr/dm/pi_1/pi_1_has_0j_hcand_1_pt_D.root'

        # plot_stack_from_root_file(root_file_to_plot, output_dir)
       
        # #########################################
        # ###         DERIVE FAKE FACTORS       ###
        # #########################################
        # # # Calculate fake factor

            # # 0 ABCD categories
            # calculate_fake_factor(input_file_region, catA0, catB0, dm, njet)

            ## ABCD categories
            calculate_fake_factor(input_file_region, catA, catB, dm, njet)

        

        # #########################################
        # ###          APPLY FAKE FACTORS       ###
        # ###             FF x C -> D           ###
        #########################################
        ## Apply the fake factor to the region C and D to get the final D region. Correction are not yet apply.
        ## It's a first set of control plot on pt distribution of the region C and D.


            ## 0 ABCD categories
            # apply_fake_factor_CD(input_file_region, catC0, catD0, dm, njet)

            ## ABCD categories
            # apply_fake_factor_CD(input_file_region, catC, catD, dm, njet)


        # # #########################################
        # # ###          APPLY FAKE FACTORS       ###
        # # ###             FF x B -> A           ###
        # # #########################################
        # # # Apply the fake factor to the region B and A to get the final A region. Correction are not yet apply. It's to get the 
        # # # closure correction of met_var_qcd_h1 in the region A.

        # variables = config['variables']
        # vars2D = [[var['var1'], var['var2']] for var in variables]
        # #print(f"Variables 2D: {vars2D}")

        # for var2D in vars2D:
        #     apply_fake_factor_AB(input_file_region, catA, catB, dm, njet, var2D)

        # # #########################################
        # # ###  APPLY FAKE FACTORS AND CCorr     ###
        # # ###             FF x C -> D           ###
        # # #########################################
        # # ## Apply the fake factor and the colsure correction to the region C and D to get the final D region.
        # # ## It's a second set of control plot on pt distribution of the region C and D.

        # apply_fake_factor_CCorr_CD(input_file_region, catC, catD, dm, njet)





if __name__ == "__main__":
    dms = ['pi_1', 'rho_1', 'a1dm2_1', 'a1dm10_1', 'a1dm11_1']
    # dms = ['pi_1']
    for dm in dms:
        config_path = f'script_FF/fake_factor_derivation/inputs/inputs_json/fake_factors_{dm}.json'
        main(config_path, dm)
