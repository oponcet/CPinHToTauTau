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
from save_correctionlib import *



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
        xmax = 120

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


    fit_up, fit_down = get_variated_function(fit_result)

    fit_up_formula = str(fit_up.GetExpFormula("P"))
    fit_down_formula = str(fit_down.GetExpFormula("P"))

    print(f"Fit formula: {fit_formula}")
    print(f" fit up formula: {fit_up_formula}")
    print(f" fit down formula: {fit_down_formula}")


    save_to_correctionlib_with_fit(fake_factor_hist, output_json_file, dm, njet, fit_formula, fit_up_formula, fit_down_formula, "fake_factor", "pt", 40, xmax)
