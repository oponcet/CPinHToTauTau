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
import re
import cmsstyle as CMS


from fit_functions import fit_fake_factor



def rebin_to_custom_bins(hist, custom_bins):
    """
    Rebins a histogram to custom bin ranges and normalizes content to bin widths.

    Args:
        hist (ROOT.TH1): The original histogram to rebin.
        custom_bins (list): List of custom bin edges (e.g., [40, 44, 48, ...]).
    
    Returns:
        ROOT.TH1: A new histogram with custom bins and normalized bin content.
    """

    print("###############################################")
    # Ensure custom_bins is a valid list
    if len(custom_bins) < 2:
        raise ValueError("custom_bins must contain at least two values to define bin ranges.")
    
    # Convert to array format for ROOT
    bin_edges_array = array('d', custom_bins)
    n_bins = len(custom_bins) - 1

    # Create a new histogram with custom binning
    rebinned_hist = ROOT.TH1F(
        f"{hist.GetName()}_rebinned",
        hist.GetTitle(),
        n_bins,
        bin_edges_array
    )

    for i in range(0,n_bins):
        new_bin_low = bin_edges_array[i] # low edge of the bin
        new_bin_high = bin_edges_array[i + 1] # high edge of the bin
        print(f"new_bin_low: {new_bin_low}, new_bin_high: {new_bin_high}")
        # Find the corresponding bin in the original histogram
        new_bin_content = 0
        n_bin_merged = 0
        new_bin_error_sq = 0  # Sum of squared errors for the new bin

        # Loop over the original histogram bins
        for j in range(1, hist.GetNbinsX() + 1):
            bin_center = hist.GetBinCenter(j)
            bin_content = hist.GetBinContent(j)
            bin_error = hist.GetBinError(j)

            # Check if the bin center is within the new bin range
            if new_bin_low <= bin_center < new_bin_high and bin_content != 0:
                # Add bin content to the new histogram
                n_bin_merged += 1
                new_bin_content += bin_content
                new_bin_error_sq += bin_error ** 2
                print(f"bin_center: {bin_center}, bin_content: {bin_content}")

        # Skip empty bins
        if n_bin_merged == 0:
            continue
        bin_width = new_bin_high - new_bin_low    
        new_bin_content /= n_bin_merged # normalize to the bin width
        new_bin_error = (new_bin_error_sq**0.5) / n_bin_merged  # Propagate uncertainty
        print(f"new_bin_content: {new_bin_content} for bin center : {rebinned_hist.GetBinCenter(i + 1)} ")
        rebinned_hist.SetBinContent(i + 1, new_bin_content)
        rebinned_hist.SetBinError(i + 1, new_bin_error)

    return rebinned_hist


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
    # ROOT.gROOT.SetStyle("Plain")
    

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
    custom_bins = [40,45,50,55,60,65,70,80,120,200]
    # custom_bins = [40,45,50,55,60,65,70,80,100,120,200]

    # custom_bins = [40,44,50,54,58,64,70,80,120,200]
    # custom_bins = [40,44,48,52,56,60,70,80,120,200]
    
    fake_factor_hist = rebin_to_custom_bins(fake_factor_hist, custom_bins)


    # Fit the fake factor histogram using `fit_functions.py`
    # fit_result = fit_fake_factor(fake_factor_hist)
    #fit_result, h_uncert, fake_factor_hist, fit_details = fit_fake_factor(fake_factor_hist, 40, 200, usePol1=True) # polOnly = None

    # if last bin = 0 
    if fake_factor_hist.GetBinContent(fake_factor_hist.GetNbinsX()) == 0:
        print("Last bin is 0")
        xmax = 120
        ff_last_bin = fake_factor_hist.GetBinContent(fake_factor_hist.GetNbinsX() - 1)
    else:
        ff_last_bin = fake_factor_hist.GetBinContent(fake_factor_hist.GetNbinsX())
        xmax = 160

    print(f"dms: {dm}, njet: {njet}")

    if dm == "pi_111" :
        fit_result, h_uncert, fake_factor_hist, fit_details = fit_fake_factor(fake_factor_hist, 40, xmax, usePol1=True) # polOnly = None

    else:
        if dm == "a1dm10_1" and njet == "has_0j":
            fit_method = 3
        elif (dm == "rho_1" and njet == "has_1j"):
            fit_method = 3
        elif dm == "pi_1":
            fit_method = 3
        else:
            fit_method = 2

        fit_result, h_uncert, fake_factor_hist, fit_details = fit_fake_factor(fake_factor_hist, 40, xmax, usePol1=False, polOnly=fit_method) # polOnly = None
    # fit_result, h_uncert, fake_factor_hist, fit_details = fit_fake_factor(fake_factor_hist, 40, xmax, usePol1=True) # polOnly = None


    # Save the fake factor and fit results to a ROOT file
    os.makedirs(os.path.dirname(output_root_file), exist_ok=True)
    output_file = ROOT.TFile.Open(output_root_file, "RECREATE")
    fake_factor_hist.Write("FakeFactor")
    fit_result.Write("FitResult")
    h_uncert.Write("Uncertainties")

    #### CANVAS FULL PLOT #### 

    # Plot the uncertainties
    uncert_canvas = ROOT.TCanvas("Fake_Factor", "Fake Factor", 800, 600)

    # Set style plain 
     
   
    # Draw the fake factor histogram
    fake_factor_hist.Draw("EP")
    fake_factor_hist.SetLineColor(ROOT.kBlack)
    fake_factor_hist.SetMarkerStyle(20)
    title = f"Fake Factor for {dm} {njet} jets"
    fake_factor_hist.SetTitle(title)    
    fake_factor_hist.GetYaxis().SetTitle("Fake Factor")
    fake_factor_hist.GetXaxis().SetTitle("p_{T} (GeV)")
    # fake_factor_hist.SetTitle("{};p_{T} (GeV) ;Fake Factor")
    # fake_factor_hist.GetYaxis().SetRangeUser(0, 1.5)
    fake_factor_hist.GetYaxis().SetRangeUser(0, 0.5)


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
    legend.AddEntry(h_uncert, "68% CL (Uncertainties)", "F")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)
    legend.Draw()

    # remove stat box
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
  

    # Save the canvas to an image file
    output_image_path = output_root_file.replace(".root", "_FullPlot.png")
    uncert_canvas.SaveAs(output_image_path)
    output_image_path = output_root_file.replace(".root", "_FullPlot.pdf")
    uncert_canvas.SaveAs(output_image_path)

    # #### PNG FILE: FIT DETAIL #### uncert_canvas
    # Display fit info like NDF, Chi2, etc.
    ROOT.gStyle.SetOptFit(1)
    # change stat box position and size
    ROOT.gStyle.SetStatY(0.9)
    ROOT.gStyle.SetStatX(0.78)
    ROOT.gStyle.SetStatW(0.15)
    ROOT.gStyle.SetStatH(0.15)
    # # Create a text box to show the fit details
    # fit_details_text = ROOT.TPaveText(0.35, 0.69, 0.6, 0.89, "NDC") # x1, y1, x2, y2, option
    # fit_details_text.SetBorderSize(0)
    # fit_details_text.SetFillColor(0)
    # fit_details_text.SetTextAlign(12)
    # fit_details_text.SetTextSize(0.03)

    # # Add fit statistics and parameters to the text box
    # fit_details_text.AddText(f"Chi2 = {fit_details['Chi2']:.4f}")
    # fit_details_text.AddText(f"NDf = {fit_details['NDf']}")

    # # Add parameter values with their errors
    # for i, param in enumerate(fit_details['Parameters']):
    #     fit_details_text.AddText(f"p{i} = {param['p']:.4f} \pm {param['error']:.4f}")

    # # Draw the text box
    # fit_details_text.Draw()

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

    save_to_correctionlib_with_fit(fake_factor_hist, fit_result, output_json_file, dm, njet, fit_formula, [fit_result.GetParameter(i) for i in range(fit_result.GetNpar())], correction_name="fake_factor", variable_name="pt",pt_min=40, pt_max=xmax)

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

   
def save_to_correctionlib_with_fit(fake_factor_hist, fit_result, output_json_file, dm, njet, fit_formula, fit_params, correction_name="fake_factor", variable_name="pt", pt_min=40, pt_max=200):
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

    fit_formula_converted = convert_fit_formula_to_correctionlib(fit_formula, min_value=pt_min, max_value=pt_max)

    if correction_name == "fake_factor":
        main_description = "Fake factors for the httcp analysis"
        name = "fake_factor"
        output_description = "Fake factor to apply to data-MC"
    else:
        main_description = "Closure correction for the httcp analysis"
        name = "closure_correction"
        output_description = "Closure correction to apply to data-MC"

    if njet == "has_0j":
        njet_int = 0
    elif njet == "has_1j":
        njet_int = 1
    elif njet == "has_2j":
        njet_int = 2
    else:
        njet_int = -1
    
    if dm == "pi_1":
        dm_int = 0
    elif dm == "rho_1":
        dm_int = 1
    elif dm == "a1dm2_1":
        dm_int = 2
    elif dm == "a1dm10_1":
        dm_int = 10
    elif dm == "a1dm11_1":
        dm_int = 11
    else:
        dm_int = -1


    correctionlib_json = {
        "schema_version": 2,
        "description": main_description,
        "corrections": [
            {
                "name": name,
                "description": main_description,
                "version": 1,
                "inputs": [
                    {
                        "name": variable_name,
                        "type": "real",
                        "description": "Transverse momentum of the tau"
                    },
                    {
                        "name": "dm",
                        "type": "int",
                        "description": "Reconstructed tau decay mode of leading tau: 0,1,2,10,11"
                    },
                    {
                        "name": "njets",
                        "type": "int",
                        "description": "Number of jets in the event 0,1,2"
                    }
                ],
                "output": {
                    "name": name,
                    "type": "real",
                    "description": output_description
                },
                "data": {
                    "nodetype": "category",
                    "input": "dm",
                    "content": [
                    {
                        "key": dm_int,
                        "value": {
                            "nodetype": "category",
                            "input": "njets",
                            "content": [
                                {
                                    "key": njet_int,
                                    "value": {
                                        "nodetype": "formula",
                                        "expression": fit_formula_converted,
                                        "parser": "TFormula",
                                        "variables": [variable_name]
                                    }
                                }
                            ]
                        }
                    }]
                }
            }
        ]
    }








    # Save to JSON file
    with open(output_json_file, "w") as f:
        json.dump(correctionlib_json, f, indent=4)
    
    print(f"CorrectionLib JSON with fit formula saved to {output_json_file}")


def convert_fit_formula_to_correctionlib(fit_formula, min_value=None, max_value=None):
    """
    Convert a ROOT.TF1 formula to a CorrectionLib-compatible formula.

    Parameters:
        fit_formula (str): Formula string extracted from ROOT.TF1 (e.g., TF1.GetExpFormula()).
        min_value (float): Minimum value for the variable (applied with max()).
        max_value (float): Maximum value for the variable (applied with min()).
    
    Returns:
        str: Converted formula compatible with CorrectionLib.
    """
    # Replace TMath functions with standard operations
    replacements = {
        "TMath::Sq(x)": "x*x",
        "TMath::Sqrt(x)": "sqrt(x)",
        "TMath::Power(x,": "pow(x,",
        "TMath::Exp(x)": "exp(x)",
        "TMath::Log(x)": "log(x)"
    }
    print("min_value: ", min_value)
    print("max_value: ", max_value)
    
    print(f"Original fit formula: {fit_formula}")
    for key, value in replacements.items():
        fit_formula = fit_formula.replace(key, value)

    # Add bounds using min and max if specified
    if min_value is not None or max_value is not None:
        var = "x"  # Default variable name
        bounded_variable = var
        if min_value is not None:
            bounded_variable = f"max({bounded_variable},{min_value})"
        if max_value is not None:
            bounded_variable = f"min({bounded_variable},{max_value})"
        # Replace the variable with its bounded version
        fit_formula = fit_formula.replace(var, bounded_variable)

    # Clean up extra spaces (optional, for neatness)
    fit_formula = re.sub(r'\s+', ' ', fit_formula).strip()

    print(f"Converted fit formula: {fit_formula}")
    
    return fit_formula