
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

from derive_FF import calculate_fake_factor
from apply_FF import *



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
        
        #########################################
        ###         DERIVE FAKE FACTORS       ###
        #########################################
        
        # Define input file
        input_file_region = f'script_FF/fake_factor_derivation/inputs/inputs_rootfile/{dm}/{dm}_{njet}.root'   # script_FF/fake_factor_derivation/inputs/inputs_rootfile/pi_1/pi_1_has_0j.root

        """
        # Calculate fake factor
        calculate_fake_factor(input_file_FF, catA, catB, dm, njet)
        """

        #########################################
        ###          APPLY FAKE FACTORS       ###
        ###             FF x C -> D           ###
        #########################################

        apply_fake_factor_CD(input_file_region, catC, catD, dm, njet)


        #########################################
        ###          APPLY FAKE FACTORS       ###
        ###             FF x B -> A           ###
        #########################################

        # variables = config['variables']
        # vars2D = [[var['var1'], var['var2']] for var in variables]
        # #print(f"Variables 2D: {vars2D}")

        # for var2D in vars2D:
        #     apply_fake_factor_AB(input_file_region, catA, catB, dm, njet, var2D)









if __name__ == "__main__":
    # dms = ['pi_1', 'rho_1', 'a1dm2_1', 'a1dm10_1', 'a1dm11_1']
    dms = ['pi_1']
    for dm in dms:
        config_path = f'script_FF/fake_factor_derivation/inputs/inputs_json/fake_factors_{dm}.json'
        main(config_path, dm)
