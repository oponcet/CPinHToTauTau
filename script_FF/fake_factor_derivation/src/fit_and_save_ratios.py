import ROOT
import pickle
import numpy as np
import os
import json
import sys
from argparse import ArgumentParser
from itertools import product  # To combine dm and n_jets lists
from fit_functions import fit_fake_factor
from save_correctionlib import *


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
        PT_RANGE = (35, 200)

    OUTPUT_DIR_BASE = f"script_FF/fake_factor_derivation/outputs/{ERA}/{CORRECTION_TYPE}"
    ensure_directory(OUTPUT_DIR_BASE)

    # Handle all combinations of DM and n_jets
    categories = list(product(DM_LIST, N_JETS_LIST))
    
    input_files = []
    output_dirs = []
    categories_str = []

    for dm, n_jets in categories:
        CATEGORY, INPUT_FILE = get_input_file(dm, n_jets, config, VARIABLE)
        output_dir = os.path.join(OUTPUT_DIR_BASE, f"{CATEGORY}")
        ensure_directory(output_dir)

        input_files.append(INPUT_FILE)
        output_dirs.append(output_dir)
        categories_str.append(CATEGORY)

    return input_files, output_dirs, HIST_NAME, categories_str, PT_RANGE


def get_input_file(dm, n_jets, config, VARIABLE):
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
        INPUT_FILE = f"{config['input_base_path']}/ratio_{VARIABLE}_dm_{dm}_njet_{n_jets}.pkl" 
    else:
        raise ValueError("Invalid DM and N_JETS values")

    return CATEGORY, INPUT_FILE


def save_root_file(OUTPUT_DIR, th1d, HIST_NAME, fit, h_uncert, ratio_hist, config, CATEGORY):
    """Save the histogram, fit results, and uncertainties to a ROOT file."""
    output_root_file = os.path.join(OUTPUT_DIR, f"{config['correction_type']}_{config['era']}_{CATEGORY}.root")
    output_file = ROOT.TFile(output_root_file, "RECREATE")
    th1d.Write(HIST_NAME)
    fit.Write("fit")
    h_uncert.Write("h_uncert")
    ratio_hist.Write("ratio_hist")
    output_file.Close()
    return output_root_file


def plot_results(fit, h_uncert, ratio_hist, OUTPUT_DIR, CATEGORY, output_root_file):
    """Create and save the plot results."""
    canvas = ROOT.TCanvas("Extrapolation Correction", "Extrapolation Correction", 800, 600)
    ratio_hist.Draw("EP")
    ratio_hist.SetLineColor(ROOT.kBlack)
    ratio_hist.SetMarkerStyle(20)
    
    ratio_hist.SetTitle("Extrapolation Correction for inclusive in DM and Njets")
    ratio_hist.GetYaxis().SetTitle("Extrapolation Correction")
    ratio_hist.GetXaxis().SetTitle("p_{T} (GeV)")
    ratio_hist.GetYaxis().SetRangeUser(0, 1.6)

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

    output_image_base = output_root_file.replace(".root", "")
    canvas.SaveAs(f"{output_image_base}_FullPlot.png")
    canvas.SaveAs(f"{output_image_base}_FullPlot.pdf")

    # -------------------------------
    # Adjust Fit Box and Save Details
    # -------------------------------
    ROOT.gStyle.SetOptStat(0)  # Disable stat box
    ROOT.gStyle.SetStatY(0.4)
    ROOT.gStyle.SetStatX(0.7)
    ROOT.gStyle.SetStatW(0.15)
    ROOT.gStyle.SetStatH(0.15)
    ROOT.gStyle.SetOptFit(1)   # Show fit info

    canvas.SaveAs(f"{output_image_base}_FitDetails.png")
    canvas.SaveAs(f"{output_image_base}_FitDetails.pdf")

    return canvas


def save_json_correction(fit, fit_up, fit_down, ratio_hist, output_root_file, config, PT_RANGE):
    """Save the fit results to a JSON file using the correctionlib format."""
    output_json_file = output_root_file.replace(".root", ".json")
    fit_formula = str(fit.GetExpFormula("P"))  # Explicitly cast to a Python string

    fit_up_formula = str(fit_up.GetExpFormula("P"))
    fit_down_formula = str(fit_down.GetExpFormula("P"))

    print(f"Fit formula: {fit_formula}")
    print(f" fit up formula: {fit_up_formula}")
    print(f" fit down formula: {fit_down_formula}")

    save_to_correctionlib_with_fit(ratio_hist, output_json_file, -1, -1, fit_formula, fit_up_formula, fit_down_formula, 
                                   config['correction_type'], config['variable'], PT_RANGE[0], PT_RANGE[1])

# -------------------------------
# Main Function
# -------------------------------
def main(args):

    config = args.config
    if not os.path.exists(args.config):
        raise ValueError(f"Configuration file {args.config} does not exist.")
    
    print("Using configuration file: %s"%(args.config))

    # Load configuration from JSON file
    with open(args.config, "r") as config_file:
        config = json.load(config_file)

    # Assign configuration values and directories
    input_files, output_dirs, HIST_NAME, categories, PT_RANGE = configure_directories(config)
    
    for INPUT_FILE, OUTPUT_DIR, CATEGORY in zip(input_files, output_dirs, categories):
        print(f"Processing category: {CATEGORY}")

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

        # -------------------------------
        # Fit Fake Factor & Save Fit Outputs
        # -------------------------------
        fit, h_uncert, ratio_hist, fit_up, fit_down = fit_fake_factor(th1d, *PT_RANGE, usePol1=False, polOnly=3)

        # Save ROOT file
        output_root_file = save_root_file(OUTPUT_DIR, th1d, HIST_NAME, fit, h_uncert, ratio_hist, config, CATEGORY)
        
        # Plot Results
        canvas = plot_results(fit, h_uncert, ratio_hist, OUTPUT_DIR, CATEGORY, output_root_file)

        # Save to JSON
        save_json_correction(fit, fit_up, fit_down, ratio_hist, output_root_file, config, PT_RANGE)


# -------------------------------
# Call the main function if script is executed
# -------------------------------
if __name__ == '__main__':

    argv = sys.argv
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', dest='config', type=str, default='script_FF/fake_factor_derivation/src/config/config_default.json', action='store', help="set config file")
    args = parser.parse_args()

    main(args)

    print(">>>\n>>> done\n")
