import ROOT
import pickle
import numpy as np
import os
import json
import sys
from array import array
from argparse import ArgumentParser
from itertools import product  # To combine dm and n_jets lists
from fit_functions import fit_fake_factor
from save_correctionlib import *
import cmsstyle as CMS
 


# Set ROOT to batch mode (no GUI popping up)
ROOT.gROOT.SetBatch(True)

# -------------------------------
# Utility Functions
# -------------------------------
def create_th1d_histogram(name, bin_edges, bin_contents, bin_uncertainties):
    """
    Create a ROOT TH1D histogram from bin contents and uncertainties.
    """
    n_bins = len(bin_contents)
    th1d = ROOT.TH1D(name, name, n_bins, bin_edges[0], bin_edges[-1])
    
    for i, (content, uncertainty) in enumerate(zip(bin_contents, bin_uncertainties), start=1):
        if content > 1.0e-5 and uncertainty < 1:
            th1d.SetBinContent(i, content)
            th1d.SetBinError(i, uncertainty)
    
    return th1d


def load_pkl_file(file_path):
    """Load a pickle file containing the ratio data."""
    with open(file_path, "rb") as file:
        return pickle.load(file)


def ensure_directory(directory):
    """Ensure a directory exists."""
    os.makedirs(directory, exist_ok=True)


def configure_directories(config):
    """Configure and return directories and necessary paths for each DM and Njet combination."""
    ERA = config["era"]
    CORRECTION_TYPE = config["correction_type"] 
    VARIABLE = config["variable"]
    HIST_NAME = config["hist_name"]
    DM_LIST = config["dm"]
    N_JETS_LIST = config["n_jets"]
    
    if VARIABLE == "met_var_qcd_h1":
        PT_RANGE = (-1.5, 1.5)
    else:
        PT_RANGE = (35, 100)

    OUTPUT_DIR_BASE = f"/afs/cern.ch/user/o/oponcet/private/analysis/CPinHToTauTau/script_FF/fake_factor_derivation/outputs/{ERA}/{CORRECTION_TYPE}"
    ensure_directory(OUTPUT_DIR_BASE)

    # Handle all combinations of DM and n_jets
    categories = list(product(DM_LIST, N_JETS_LIST))
    
    input_files = []
    output_dirs = []
    categories_str = []

    for dm, n_jets in categories:
        CATEGORY, INPUT_FILE = get_input_file(dm, n_jets, config, VARIABLE,CORRECTION_TYPE)
        output_dir = os.path.join(OUTPUT_DIR_BASE, f"{CATEGORY}")
        ensure_directory(output_dir)

        input_files.append(INPUT_FILE)
        output_dirs.append(output_dir)
        categories_str.append(CATEGORY)

    return input_files, output_dirs, HIST_NAME, categories_str, PT_RANGE


def get_input_file(dm, n_jets, config, VARIABLE,CORRECTION_TYPE):
    """Return the input file path based on dm, n_jets, and the given config."""
    if dm == -1 and n_jets == -1:
        CATEGORY = "inclusive"
        INPUT_FILE = f"{config['input_base_path']}/RATIO_{VARIABLE}_{CATEGORY}.pkl"
    elif dm == -1 and n_jets != -1:
        CATEGORY = n_jets
        INPUT_FILE = f"{config['input_base_path']}/RATIO_{VARIABLE}_{CATEGORY}.pkl"
    elif dm != -1 and n_jets == -1:
        CATEGORY = dm
        INPUT_FILE = f"{config['input_base_path']}/RATIO_{VARIABLE}_{CATEGORY}.pkl"
    elif dm != -1 and n_jets != -1:
        CATEGORY = f"{dm}_{n_jets}"
        INPUT_FILE = f"{config['input_base_path']}/{CORRECTION_TYPE}_{VARIABLE}_dm_{dm}_njet_{n_jets}.pkl" 
    else:
        raise ValueError("Invalid DM and N_JETS values")

    return CATEGORY, INPUT_FILE


def save_root_file(OUTPUT_DIR, th1d, HIST_NAME, fit, h_uncert, ratio_hist, config, CATEGORY, combine=False):
    """Save the histogram, fit results, and uncertainties to a ROOT file."""
    if combine == True:
        era = config["era"].split("_")[0]  #only keep the year
    else :
        era = config["era"]
    output_root_file = os.path.join(OUTPUT_DIR, f"{config['correction_type']}_{era}_{CATEGORY}.root")
    output_file = ROOT.TFile(output_root_file, "RECREATE")
    th1d.Write(HIST_NAME)
    fit.Write("fit")
    h_uncert.Write("h_uncert")
    ratio_hist.Write("ratio_hist")
    output_file.Close()
    return output_root_file


def plot_results(fit, h_uncert, ratio_hist, OUTPUT_DIR, CATEGORY, output_root_file, lumi=1):
    """Create and save the plot results."""
    canvas = ROOT.TCanvas("Extrapolation Correction", "Extrapolation Correction", 800, 600)
    ratio_hist.Draw("EP")
    ratio_hist.SetLineColor(ROOT.kBlack)
    ratio_hist.SetMarkerStyle(20)
    
    ratio_hist.SetTitle("Extrapolation Correction for inclusive in DM and Njets")
    ratio_hist.GetYaxis().SetTitle("Extrapolation Correction")
    ratio_hist.GetXaxis().SetTitle("p_{T} (GeV)")
    ratio_hist.GetYaxis().SetRangeUser(0, 1.0)

    h_uncert.Draw("E3 SAME")
    h_uncert.SetFillColorAlpha(ROOT.kAzure + 7, 0.3)
    h_uncert.SetLineColor(ROOT.kAzure + 7)
    h_uncert.SetLineWidth(2)

    fit.Draw("SAME")
    fit.SetLineColor(ROOT.kAzure + 7)

    legend = ROOT.TLegend(0.15, 0.75, 0.35, 0.9)
    legend.AddEntry(ratio_hist, "Extrapolation Correction", "EP")
    legend.AddEntry(fit, "Fit Result", "L")
    legend.AddEntry(h_uncert, "68% CL (Uncertainties)", "F")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)
    legend.Draw()

    ROOT.gStyle.SetOptFit(0) 

    # Add luminosity label with cms style 
    luminosity_label = ROOT.TLatex()
    luminosity_label.SetNDC()
    luminosity_label.SetTextFont(42)
    luminosity_label.SetTextSize(0.03)
    # display only 3 digits
    lumi = "{:.2f}".format(lumi)
    luminosity_label.DrawLatex(0.72, 0.91, f"{lumi} fb^{{-1}} (13.6 TeV)")

    output_image_base = output_root_file.replace(".root", "")
    canvas.SaveAs(f"{output_image_base}_FullPlot.root")
    canvas.SaveAs(f"{output_image_base}_FullPlot.png")
    canvas.SaveAs(f"{output_image_base}_FullPlot.pdf")

    # -------------------------------
    # Adjust Fit Box and Save Details
    # -------------------------------
    ROOT.gStyle.SetOptStat(0)  # Disable stat box
    ROOT.gStyle.SetStatY(0.8)
    ROOT.gStyle.SetStatX(0.7)
    ROOT.gStyle.SetStatW(0.15)
    ROOT.gStyle.SetStatH(0.15)
    ROOT.gStyle.SetOptFit(1)   # Show fit info

    canvas.SaveAs(f"{output_image_base}_FitDetails.root")
    canvas.SaveAs(f"{output_image_base}_FitDetails.png")
    canvas.SaveAs(f"{output_image_base}_FitDetails.pdf")

    return canvas

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
    rebinned_hist = ROOT.TH1D(
        f"{hist.GetName()}_rebinned",
        hist.GetTitle(),
        n_bins,
        bin_edges_array
    )

    for i in range(0,n_bins):
        new_bin_low = bin_edges_array[i] # low edge of the bin
        new_bin_high = bin_edges_array[i + 1] # high edge of the bin
        # print(f"new_bin_low: {new_bin_low}, new_bin_high: {new_bin_high}")
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
                # print(f"bin_center: {bin_center}, bin_content: {bin_content}")

        # Skip empty bins
        if n_bin_merged == 0:
            continue
        bin_width = new_bin_high - new_bin_low    
        new_bin_content /= n_bin_merged # normalize to the bin width
        new_bin_error = (new_bin_error_sq**0.5) / n_bin_merged  # Propagate uncertainty
        # print(f"new_bin_content: {new_bin_content} for bin center : {rebinned_hist.GetBinCenter(i + 1)} ")
        rebinned_hist.SetBinContent(i + 1, new_bin_content)
        rebinned_hist.SetBinError(i + 1, new_bin_error)

    return rebinned_hist

def save_json_correction(fit, fit_up, fit_down, ratio_hist, output_root_file, config, PT_RANGE, dm, njet):
    """Save the fit results to a JSON file using the correctionlib format."""
    output_json_file = output_root_file.replace(".root", ".json")
    fit_formula = str(fit.GetExpFormula("P"))  # Explicitly cast to a Python string

    fit_up_formula = str(fit_up.GetExpFormula("P"))
    fit_down_formula = str(fit_down.GetExpFormula("P"))

    print(f"Fit formula: {fit_formula}")
    print(f" fit up formula: {fit_up_formula}")
    print(f" fit down formula: {fit_down_formula}")

    save_to_correctionlib_with_fit(ratio_hist, output_json_file, dm, njet, fit_formula, fit_up_formula, fit_down_formula, 
                                   config['correction_type'], config['variable'], PT_RANGE[0], PT_RANGE[1])

def merge_histograms_lumi_weighted(hist_list, lumi_list):
    """Merge histograms using luminosity-weighted averaging."""
    total_lumi = sum(lumi_list)

    # verify that number of hsit un nb files
    if len(hist_list) != len(lumi_list):
        raise ValueError("Number of histograms and luminosities must be equal.")

    merged_hist = hist_list[0].Clone()  # Clone first histogram to maintain structure
    
    # Reset values to zero
    merged_hist.Reset()

    hist1 = hist_list[0]
    hist2 = hist_list[1]

    lumi1 = lumi_list[0]
    lumi2 = lumi_list[1]

    # Test if they have the same binning
    if hist1.GetNbinsX() != hist2.GetNbinsX():
        raise ValueError("Histograms must have the same number of bins.")

    # should have the same binning
    for bin in range(1, hist1.GetNbinsX() + 1):
        bin_content_hist1 = hist1.GetBinContent(bin)
        bin_content_hist2 = hist2.GetBinContent(bin)

        bin_uncert_hist1 = hist1.GetBinError(bin)
        bin_uncert_hist2 = hist2.GetBinError(bin)

        # Calculate the weighted average
        merged_bin_content = (lumi1 * bin_content_hist1 + lumi2 * bin_content_hist2) / total_lumi
        merged_bin_uncert = np.sqrt((lumi1 * bin_uncert_hist1**2 + lumi2 * bin_uncert_hist2**2) / total_lumi)

        merged_hist.SetBinContent(bin, merged_bin_content)
        merged_hist.SetBinError(bin, merged_bin_uncert)

    return merged_hist


# -------------------------------
# Main Function
# -------------------------------
def main(args):

    config_files = args.config

    nb_files = len(config_files)

    combined_ratios = {}

    for config_file in config_files:

        if not os.path.exists(config_file):
            raise ValueError(f"Configuration file {config_file} does not exist.")
        
        print("Using configuration file: %s"%(config_file))

        # Load configuration from JSON file
        with open(config_file, "r") as config_file:
            config = json.load(config_file)
        
        lumi = config.get("luminosity", 1)

        # Assign configuration values and directories
        input_files, output_dirs, HIST_NAME, categories, PT_RANGE = configure_directories(config)
        
        for INPUT_FILE, OUTPUT_DIR, CATEGORY in zip(input_files, output_dirs, categories):
            print(f"Processing category: {CATEGORY}")

            dm = CATEGORY.split("_")[0]
            njet = CATEGORY.split("_")[1]


            # -------------------------------
            # Load Data
            # -------------------------------
            ratio = load_pkl_file(INPUT_FILE)
            h = ratio[0, :, :]

            bin_edges = h.axes[1].edges
            bin_contents = h.values().flatten()
            bin_uncertainties = np.sqrt(h.variances()).flatten()

            print(f"Loaded data from {INPUT_FILE}")

            th1d = create_th1d_histogram(HIST_NAME, bin_edges, bin_contents, bin_uncertainties)
            th1d.SetMinimum(0)
            th1d.SetMaximum(1.6)
   

            custom_bins = [35,40,45,50,55,60,65,70,80,120,200] # Custom binning for hcand_1_pt
            th1d = rebin_to_custom_bins(th1d, custom_bins)


            if nb_files > 1:
                # Store results for merging
                key = (dm, njet)
                if key not in combined_ratios:
                    combined_ratios[key] = {"hists": [], "lumis": []}
                 
                combined_ratios[key]["hists"].append(th1d)
                combined_ratios[key]["lumis"].append(lumi)


            # -------------------------------
            # Fit Fake Factor & Save Fit Outputs
            # -------------------------------


            fit_range = (35, 200)
            fit, h_uncert, ratio_hist, fit_up, fit_down = fit_fake_factor(th1d, *fit_range, usePol1=False, polOnly=3)

            # Save ROOT file
            output_root_file = save_root_file(OUTPUT_DIR, th1d, HIST_NAME, fit, h_uncert, ratio_hist, config, CATEGORY)
            
            # Plot Results
            canvas = plot_results(fit, h_uncert, ratio_hist, OUTPUT_DIR, CATEGORY, output_root_file, lumi=lumi)

            # Save to JSON
            save_json_correction(fit, fit_up, fit_down, ratio_hist, output_root_file, config, PT_RANGE, dm, njet)

    if nb_files > 1:
        # Merge histograms for the same (dm, njet) category
        for (dm, njet), data in combined_ratios.items():
            print(f"Combining histograms for dm={dm}, njet={njet} using luminosity weighting")

            print("data: ", data)

            merged_hist = merge_histograms_lumi_weighted(data["hists"], data["lumis"])

            combined_era = config["era"].split("_")[0]  #only keep the year
            CORRECTION_TYPE = config["correction_type"] 

            OUTPUT_DIR_BASE = f"/afs/cern.ch/user/o/oponcet/private/analysis/CPinHToTauTau/script_FF/fake_factor_derivation/outputs/{combined_era}/{CORRECTION_TYPE}"
            ensure_directory(OUTPUT_DIR_BASE)

            OUTPUT_DIR = os.path.join(OUTPUT_DIR_BASE, f"{dm}_{njet}")
            ensure_directory(OUTPUT_DIR)


            # Fit Fake Factor & Save Fit Outputs
            fit_range = (35, 200)
            fit, h_uncert, ratio_hist, fit_up, fit_down = fit_fake_factor(merged_hist, *fit_range, usePol1=False, polOnly=3)

            # Save ROOT file
            output_root_file = save_root_file(OUTPUT_DIR, merged_hist, HIST_NAME, fit, h_uncert, ratio_hist, config, f"{dm}_{njet}", combine=True)

            # Plot Results
            canvas = plot_results(fit, h_uncert, ratio_hist, OUTPUT_DIR, f"{dm}_{njet}", output_root_file, lumi=sum(data["lumis"]))

            # Save to JSON
            save_json_correction(fit, fit_up, fit_down, ratio_hist, output_root_file, config, PT_RANGE, dm, njet)
# -------------------------------
# Call the main function if script is executed
# -------------------------------
if __name__ == '__main__':

    argv = sys.argv
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', dest='config', type=str, nargs='+', default=['script_FF/fake_factor_derivation/src/config/config_default.json'], action='store', help="set config file (multiple files supported)")
    args = parser.parse_args()

    main(args)

    print(">>>\n>>> done\n")
