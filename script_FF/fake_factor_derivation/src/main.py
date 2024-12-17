
'''
Author : @oponcet
Date : 10-12-2024
Main script calling the function of Fake Factor derivation, closure and extrapolation 
correction, control plots ... 

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




def main(config_path, dm):

    # ROOT configuration
    ROOT.gROOT.SetBatch(True) # Do not display canvas

    # Load the JSON configuration
    with open(config_path, 'r') as f:
        config = json.load(f)

    categories = config['categories']
    variable = "hcand_1_pt"

    catA, catB = None, None  # Initialize to None

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
            
        if not catA or not catB:
            raise ValueError(f"Categories _hadA__ or _hadB__ not found for {dm} in njet {njet}")   
        

        # Define input file
        input_file_region = f'script_FF/fake_factor_derivation/inputs/inputs_rootfile/{dm}/{dm}_{njet}.root'   # script_FF/fake_factor_derivation/inputs/inputs_rootfile/pi_1/pi_1_has_0j.root

        
        #########################################
        ###        PLOT THE HCAND_1_PT        ###
        #########################################
        ''' WIP
        
        # Replace the following paths with actual file paths and names
        variable = "hcand_1_pt"
        canvas_name = f"{catD}/{variable}/canvas"
        # mc_stack = canvas_name["mc_stack"]  # THStack
        # data_hist = canvas_name["data_hist"]  # TH1D


        output_dir = f"script_FF/fake_factor_derivation/inputs/inputs_rootfile/{dm}/plots"  # output/plots

        # Ensure directories exist if not create them
        os.makedirs(output_dir, exist_ok=True)

        # Extract histograms from the ROOT file
        with uproot.open(input_file_region) as root_file:
            # Extract MC histograms
            mc_stack = root_file[f"{catD}/{variable}/canvas/mc_stack"]
            print(f"mc_stack: {mc_stack}")
            mc_histograms = {
                hist.name: (hist.values(), hist.axis().edges()) for hist in mc_stack.values()
            }

            # Extract Data histogram (if available)
            data_histogram = None
            if data_histogram_name in root_file:
                data_hist = root_file[f"{catD}/{variable}/canvas/data_hist"]
                print(f"data_hist: {data_hist}")
                data_histogram = (data_hist.values(), data_hist.axis().edges())


        # Create the plot
        create_stack_plot_hep(mc_stack, data_hist, variable, output_dir)
        '''
        #########################################
        ###         DERIVE FAKE FACTORS       ###
        #########################################
        # # Calculate fake factor

        calculate_fake_factor(input_file_region, catA, catB, dm, njet)
        

        #########################################
        ###          APPLY FAKE FACTORS       ###
        ###             FF x C -> D           ###
        #########################################
        ## Apply the fake factor to the region C and D to get the final D region. Correction are not yet apply.
        ## It's a first set of control plot on pt distribution of the region C and D.

        apply_fake_factor_CD(input_file_region, catC, catD, dm, njet)


        #########################################
        ###          APPLY FAKE FACTORS       ###
        ###             FF x B -> A           ###
        #########################################
        # Apply the fake factor to the region B and A to get the final A region. Correction are not yet apply. It's to get the 
        # closure correction of met_var_qcd_h1 in the region A.

        variables = config['variables']
        vars2D = [[var['var1'], var['var2']] for var in variables]
        #print(f"Variables 2D: {vars2D}")

        for var2D in vars2D:
            apply_fake_factor_AB(input_file_region, catA, catB, dm, njet, var2D)

        #########################################
        ###  APPLY FAKE FACTORS AND CCorr     ###
        ###             FF x C -> D           ###
        #########################################
        ## Apply the fake factor and the colsure correction to the region C and D to get the final D region.
        ## It's a second set of control plot on pt distribution of the region C and D.

        apply_fake_factor_CCorr_CD(input_file_region, catC, catD, dm, njet)





if __name__ == "__main__":
    # dms = ['pi_1', 'rho_1', 'a1dm2_1', 'a1dm10_1', 'a1dm11_1']
    dms = ['pi_1']
    for dm in dms:
        config_path = f'script_FF/fake_factor_derivation/inputs/inputs_json/fake_factors_{dm}.json'
        main(config_path, dm)
