'''
Author : @oponcet
Date : 06-12-2024
Script calculate the fakfe factor for the different regions (DM/DM_Njets). Call fit fit_functions.py to fit the fake factor. 
Retrun the fake factor in a root file and in json file.

Input: 
    - inputs/inputs_json/region.json (json file with the region configuration)
    - inputs/inputs_rootfile/dm/region.root/A/pt1/data_minus_mc (TH1D)
    - inputs/inputs_rootfile/dm/region.root/B/pt1/data_minus_mc (TH1D)

Output:
    - outputs/FakeFactors/dm/region.root (fake factor in TH1D) 
    - outputs/FakeFactors/dm/region.json (fake factor in json)
'''

import os
import json
import ROOT
from correctionlib import schemav2 as cs
from array import array

from fit_functions import fit_fake_factor


def calculate_fake_factor(input_file, catA, catB, dm, njet):
    """
    Calculates the fake factor from two input histograms and saves the results.

    Parameters:
    - input_file: Path to input ROOT file for region A. 
    - catA: category A (numerator) -> TH1D is in catA/variable/data_minus_mc
    - catB: category B (denominator) -> TH1D is in catB/variable/data_minus_mc
    - dm: decay mode
    - njet: number of jets

    """

    # Reset the style of the canvas to use default style
    ROOT.gROOT.SetStyle("Plain")
    

    # Define output paths
    output_root_file = f'script_FF/fake_factor_derivation/outputs/FakeFactors/{dm}/{dm}_{njet}.root'
    output_json_file = f'script_FF/fake_factor_derivation/outputs/FakeFactors/{dm}/{dm}_{njet}.json'

    # Ensure directories exist if not create them
    os.makedirs(os.path.dirname(output_root_file), exist_ok=True)
    os.makedirs(os.path.dirname(output_json_file), exist_ok=True)
    

    # Load input ROOT files and histograms
    root_file = ROOT.TFile.Open(input_file, "READ")
   

    if not root_file or root_file.IsZombie():
        raise FileNotFoundError("Input ROOT files could not be opened.")

    hist_a = root_file.Get(f'{catA}/hcand_1_pt/data_minus_mc')
    hist_b = root_file.Get(f'{catB}/hcand_1_pt/data_minus_mc')

    if not hist_a or not hist_b:
        raise KeyError(f"Histogram '{catA}/hcand_1_pt/data_minus_mc not found in one of the input files.")

    # Calculate fake factor (histogram division)
    fake_factor_hist = hist_a.Clone("fake_factor")
    fake_factor_hist.Divide(hist_b)


    # Rebin the fake factor histogram like this [40,45,55,60,65,70,80,90,100,120,140,200]
    custom_bins = [40,45,50,55,60,65,70,80,90,100,120,140,200]
    n_bins = len(custom_bins) - 1

    # # print(array('d', custom_bins))

    fake_factor_hist_rebinned = fake_factor_hist.Rebin(n_bins, "fake_factor_rebinned", array('d', custom_bins))

    # normalize to the bin width
    for i in range(1, fake_factor_hist_rebinned.GetNbinsX() + 1):
        bin_width = fake_factor_hist_rebinned.GetBinWidth(i)/5. # divide by 5 because the bin width is 5 GeV
        print(f"bin_width: {bin_width}")
        bin_content = fake_factor_hist_rebinned.GetBinContent(i)
        bin_error = fake_factor_hist_rebinned.GetBinError(i)
        fake_factor_hist_rebinned.SetBinContent(i, bin_content / bin_width)
        fake_factor_hist_rebinned.SetBinError(i, bin_error / bin_width)

    fake_factor_hist = fake_factor_hist_rebinned


    # Fit the fake factor histogram using `fit_functions.py`
    # fit_result = fit_fake_factor(fake_factor_hist)
    #fit_result, h_uncert, fake_factor_hist, fit_details = fit_fake_factor(fake_factor_hist, 40, 200, usePol1=True) # polOnly = None

    # if last bin = 0 
    if fake_factor_hist.GetBinContent(fake_factor_hist.GetNbinsX()) == 0:
        print("Last bin is 0")
        xmax = 140
        ff_last_bin = fake_factor_hist.GetBinContent(fake_factor_hist.GetNbinsX() - 1)
    else:
        ff_last_bin = fake_factor_hist.GetBinContent(fake_factor_hist.GetNbinsX())
        xmax = 200

    print(f"dms: {dm}, njet: {njet}")
    if dm == 'a1dm2_1' and njet == "has_2j":
        fit_method = 3
    else:
        fit_method = 2

    fit_result, h_uncert, fake_factor_hist, fit_details = fit_fake_factor(fake_factor_hist, 40, xmax, usePol1=False, polOnly=fit_method) # polOnly = None


    # Save the fake factor and fit results to a ROOT file
    os.makedirs(os.path.dirname(output_root_file), exist_ok=True)
    output_file = ROOT.TFile.Open(output_root_file, "RECREATE")
    fake_factor_hist.Write("FakeFactor")
    fit_result.Write("FitResult")
    h_uncert.Write("Uncertainties")

    #### CANVAS FULL PLOT #### 

    # Plot the uncertainties
    uncert_canvas = ROOT.TCanvas("Fake_Factor", "Fake Factor", 800, 600)
    
   
    # Draw the fake factor histogram
    fake_factor_hist.Draw("EP")
    fake_factor_hist.SetLineColor(ROOT.kBlack)
    fake_factor_hist.SetTitle("Fit Fake Factor;p_{T} (GeV) ;Fake Factor")
    fake_factor_hist.GetYaxis().SetRangeUser(0, 2.5)

    # Draw the uncertainties as a filled area
    h_uncert.Draw("E3 SAME")
    h_uncert.SetFillColorAlpha(ROOT.kAzure + 7 , 0.3)
    h_uncert.SetLineColor(ROOT.kAzure + 7)
    h_uncert.SetLineWidth(2)

    # Draw the fit result
    fit_result.Draw("SAME")
    fit_result.SetLineColor(ROOT.kAzure + 7)

    # Add a legend for clarity
    legend = ROOT.TLegend(0.15, 0.75, 0.35, 0.9) #  
    legend.AddEntry(fake_factor_hist, "Fake Factor", "EP")
    legend.AddEntry(fit_result, "Fit Result", "L")
    legend.AddEntry(h_uncert, "95% CL (Uncertainties)", "F")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.Draw()

    # remove stat box
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)

    # Save the canvas to an image file
    output_image_path = output_root_file.replace(".root", "_FullPlot.png")
    uncert_canvas.SaveAs(output_image_path)

    #### PNG FILE: FIT DETAIL #### uncert_canvas
    # Create a text box to show the fit details
    fit_details_text = ROOT.TPaveText(0.35, 0.69, 0.6, 0.89, "NDC") # x1, y1, x2, y2, option
    fit_details_text.SetBorderSize(0)
    fit_details_text.SetFillColor(0)
    fit_details_text.SetTextAlign(12)
    fit_details_text.SetTextSize(0.03)

    # Add fit statistics and parameters to the text box
    fit_details_text.AddText(f"Chi2 = {fit_details['Chi2']:.4f}")
    fit_details_text.AddText(f"NDf = {fit_details['NDf']}")

    # Add parameter values with their errors
    for i, param in enumerate(fit_details['Parameters']):
        fit_details_text.AddText(f"p{i} = {param['p']:.4f} \pm {param['error']:.4f}")

    # Draw the text box
    fit_details_text.Draw()

    # Save the fit details canvas as a separate PNG file
    output_details_image_path = output_root_file.replace(".root", "_FitDetails.png")
    uncert_canvas.SaveAs(output_details_image_path)

    # Save the canvas to the ROOT file for later use
    uncert_canvas.Write("Fullplot")

    output_file.Close()


    #### JSON FILE ####

    # Save the fake factor to a JSON file
    os.makedirs(os.path.dirname(output_json_file), exist_ok=True)

    fit_formula = str(fit_result.GetExpFormula("P"))  # Explicitly cast to a Python string

    save_to_correctionlib_with_fit(fake_factor_hist, fit_result, output_json_file, dm, njet, fit_formula, [fit_result.GetParameter(i) for i in range(fit_result.GetNpar())], correction_name="fake_factor", variable_name="pt", pt_max=xmax, pt_default=ff_last_bin, range_speficied=False)

    # fake_factor_json = {
    #     "bin_edges": list(fake_factor_hist.GetXaxis().GetXbins()),
    #     "bin_values": [fake_factor_hist.GetBinContent(i) for i in range(1, fake_factor_hist.GetNbinsX() + 1)],
    #     "bin_errors": [fake_factor_hist.GetBinError(i) for i in range(1, fake_factor_hist.GetNbinsX() + 1)]
    # }

    # with open(output_json_file, "w") as json_file:
    #     json.dump(fake_factor_json, json_file, indent=4)

    # # Close input files
    # root_file.Close()

    # print(f"Fake factor saved to {output_root_file} and {output_json_file}")

   
def save_to_correctionlib_with_fit(fake_factor_hist, fit_result, output_json_file, dm, njet, fit_formula, fit_params, correction_name="fake_factor", variable_name="pt", pt_max=200, pt_default=1.0, range_speficied=False):
    """
    Converts ROOT histogram with fake factors and a fit function into CorrectionLib JSON format.

    Parameters:
    - fake_factor_hist: ROOT.TH1D object containing fake factors.
    - fit_result: ROOT.TF1 fit result.
    - output_json_file: Path to save the CorrectionLib-compatible JSON file.
    - dm: Decay mode.
    - njet: Number of jets.
    - fit_formula: Formula string of the fit (TF1.GetExpFormula()).
    - fit_params: List of fit parameters (TF1.GetParameter()).
    """
    # Extract bin edges and values from the histogram
    bin_edges = []
    bin_values = []
    uncertainties = []
    
    for i in range(1, fake_factor_hist.GetNbinsX() + 1):
        bin_edges.append(fake_factor_hist.GetBinLowEdge(i))
        bin_values.append(fake_factor_hist.GetBinContent(i))
        uncertainties.append(fake_factor_hist.GetBinError(i))
    
    # Add the upper edge of the last bin
    bin_edges.append(fake_factor_hist.GetBinLowEdge(fake_factor_hist.GetNbinsX() + 1))
    
    # # Create CorrectionLib format with both binning and formula
    # correction_binned = {
    #     "name": f"fake_factor_{dm}_{njet}_binned",
    #     "description": f"Binned fake factor for decay mode {dm}, {njet} jets",
    #     "version": 1,
    #     "inputs": [{"name": "pt", "type": "real"}],
    #     "output": {"name": "fake_factor", "type": "real"},
    #     "data": {
    #         "nodetype": "binning",
    #         "inputs": ["pt"],
    #         "edges": bin_edges,
    #         "content": [{"value": val, "uncertainty": unc} for val, unc in zip(bin_values, uncertainties)],
    #         "flow": "clamp"
    #     }
    # }

    if range_speficied:
        fit_expression = f"""
            if(pt > {pt_max}) {{ fake_factor = {pt_default}; }}
            else {{ fake_factor = {fit_formula}; }}"""
    else:
        fit_expression = fit_formula
    
    # Create CorrectionLib format with the formula
    correction_fit = {
        "name": f"{correction_name}_{dm}_{njet}_fit",
        "description": f"Fit fake factor for decay mode {dm}, {njet} jets",
        "version": 1,
        "inputs": [{"name": "pt", "type": "real", "description": "Transverse momentum of the hadronic tau"}],
        "output": {"name": "fake_factor", "type": "real", "description": "Fake factor for the hadronic tau"},
        "data": {
            "nodetype": "formula",
            "expression": fit_expression,
            "parameters": [{"name": f"p{i}", "value": fit_params[i]} for i in range(len(fit_params))],
            "variables": [{"name": variable_name, "type": "real"}]
        }
    }

    # Combine both corrections
    correctionlib_json = {
        "schema_version": 2,
        "description": f"Fake factor for decay mode {dm}, {njet} jets for 2022_preEE",
        "corrections": [correction_fit]
    }

    # Save to JSON file
    with open(output_json_file, "w") as f:
        json.dump(correctionlib_json, f, indent=4)
    
    print(f"CorrectionLib JSON with fit formula saved to {output_json_file}")




# if __name__ == "__main__":
#     # Define input and output paths
#     input_file_a = "inputs/inputs_rootfile/dm/region.root/A/pt1/data_minus_mc"
#     input_file_b = "inputs/inputs_rootfile/dm/region.root/B/pt1/data_minus_mc"
#     output_root_file = "outputs/FakeFactors/dm/region.root"
#     output_json_file = "outputs/FakeFactors/dm/region.json"

#     # Calculate fake factor
#     calculate_fake_factor(input_file_a, input_file_b, output_root_file, output_json_file)

# def main(config_path, dm):

#     # ROOT configuration
#     ROOT.gROOT.SetBatch(True) # Do not display canvas

#     # Load the JSON configuration
#     with open(config_path, 'r') as f:
#         config = json.load(f)

#     categories = config['categories']
#     variable = "hcand_1_pt"

#     catA, catB = None, None  # Initialize to None

#     for njet in config['categories'][dm]:
#         for cat in config['categories'][dm][njet]['ABCD']:

#             if '_hadA__' in cat:
#                 catA = cat
#             if '_hadB__' in cat:
#                 catB = cat
            
#         if not catA or not catB:
#             raise ValueError(f"Categories _hadA__ or _hadB__ not found for {dm} in njet {njet}")   
        
#         # Define input file
#         input_file = f'script_FF/fake_factor_derivation/inputs/inputs_rootfile/{dm}/{dm}_{njet}.root'   # script_FF/fake_factor_derivation/inputs/inputs_rootfile/pi_1/pi_1_has_0j.root

#         # Calculate fake factor
#         calculate_fake_factor(input_file, catA, catB, dm, njet)

# if __name__ == "__main__":
#     # dms = ['pi_1', 'rho_1', 'a1dm2_1', 'a1dm10_1', 'a1dm11_1']
#     dms = ['pi_1']
#     for dm in dms:
#         config_path = f'script_FF/fake_factor_derivation/inputs/inputs_json/fake_factors_{dm}.json'
#         main(config_path, dm)


