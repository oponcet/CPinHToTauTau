'''
Author : @oponcet
Date : 06-12-2024
Script to apply pt-depedent Fake Factors to distribution of any other variable 
(met_var_qcd_h1 for instance). 
Apply FF to B data_min_mc and use it in A region as Fake jets -> tau_h 

Input: 
    - inputs/inputs_json/region.json (json file with the region configuration)
    - outputs/FakeFactors/dm/region.json (json file with the fake factors) het the function to apply 
    - inputs/inputs_rootfile/dm/region.root/A/met_var_qcd_h1_hcand_1_pt/data_minus_mc (TH2D)
    - inputs/inputs_rootfile/dm/region.root/B/met_var_qcd_h1_hcand_1_pt/data_minus_mc (TH2D)

Output:
    - outputs/outputs_applyFF/dm/region.root/A/met_var_qcd_h1/fakes_jets (TH1D)
                                                             /canvas_FF (THstack)  
    - outputs/outputs_applyFF/dm/region_FF_apply_catA.png (png/pdf)                                    

'''
import os
import json
import ROOT

def apply_fake_factor(input_file, catA, catB, dm, njet,var2D):
    """
    Parameters:
    - input_file: Path to input ROOT file for region A. 
    - catA: category A (numerator) -> TH1D is in catA/variable/data_minus_mc
    - catB: category B (denominator) -> TH1D is in catB/variable/data_minus_mc
    - dm: decay mode
    - njet: number of jets
    - var2D: variable to apply the fake factor

    """
    # Define output paths
    output_root_file = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/dm/{dm}/{dm}_{njet}.root'
    output_png_file = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/dm/{dm}/{dm}_{njet}.png'

     # Ensure directories exist if not create them
    os.makedirs(os.path.dirname(output_root_file), exist_ok=True)
    

    # Load input ROOT files and histograms
    root_file = ROOT.TFile.Open(input_file, "READ")

    
    if not root_file or root_file.IsZombie():
        raise FileNotFoundError("Input ROOT files could not be opened.")

    var1 = var2D[0]
    var2 = var2D[1]

    # Get the 2D histograms
    hist_b = root_file.Get(f'{catB}/{var1}-{var2}/data_minus_mc').Clone()

    hist_b.SetDirectory(0) # Detach histogram from ROOT file to avoid deletion

    if not hist_b:
        raise KeyError(f"Histogram '{catA}/{var1}-{var2}/data_minus_mc not found in one of the input files.")

    # close the input ROOT file
    root_file.Close()


    # Get the fake factor
    fake_factor_file = f'script_FF/fake_factor_derivation/outputs/FakeFactors/{dm}/{dm}_{njet}.root'


    # Load input ROOT files and histograms
    root_fake_factor_file = ROOT.TFile.Open(fake_factor_file, "READ")
    print(f"Fake factor file: {fake_factor_file}")

    # Get the TF1 fake factor from the fit 
    FF_function = root_fake_factor_file.Get(f'FitResult')

    print(f"Fake factor function: {FF_function}")

    # ## Apply FF_function to hist_b
    hist_b_fake = hist_b.Clone()
    
    for i in range(hist_b_fake.GetNbinsX()):
        for j in range(hist_b_fake.GetNbinsY()):
            pt = hist_b_fake.GetYaxis().GetBinCenter(i)
            eta = hist_b_fake.GetXaxis().GetBinCenter(j)
            fake_factor = FF_function.Eval(pt)
            hist_b_fake.SetBinContent(i,j,hist_b_fake.GetBinContent(i,j)*fake_factor)




    ## Project the 2D histogram ton the hcand_1_pt axis
    hist_b_fake_proj = hist_b_fake.ProjectionY() 

    # Save the output ROOT file
    output_root_file = ROOT.TFile(output_root_file, "RECREATE")
    output_root_file.cd()
    hist_b_fake_proj.Write("fakes_jets")
    output_root_file.Close()





