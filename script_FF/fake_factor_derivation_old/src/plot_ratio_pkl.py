'''
Author : @oponcet
Date : 29-01-2025

Script to plot the ratio from the pickle file and fit the fake factor with a Landau function and/or polynomial functions and save in a root file.

'''

import ROOT
import pickle
import numpy as np
import os
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


# -------------------------------
# Configuration
# -------------------------------
CORRECTION_TYPE = "extrapolation_correction" # "extrapolation_correction" or "closure_correction" or "FakeFactor"
CATEGORY = "inclusive"  # "inclusive" or "pi_1_has_Ojet" ...
VARIABLE = "met_var_qcd_h1"
INPUT_FILE = f"/afs/cern.ch/user/o/oponcet/private/analysis/CPinHToTauTau/script_FF/fake_factor_derivation/outputs/extrapolation_correction/RATIO_{VARIABLE}_{CATEGORY}.pkl"
OUTPUT_DIR = f"script_FF/fake_factor_derivation/outputs/{CORRECTION_TYPE}/"
HIST_NAME = "ratio"
if VARIABLE == "met_var_qcd_h1":
    PT_RANGE = (-1.5, 1.5)  # Fit range
else:
    PT_RANGE = (35, 200)  # Fit range

ensure_directory(OUTPUT_DIR)

# -------------------------------
# Load Data
# -------------------------------
ratio = load_pkl_file(INPUT_FILE)
h = ratio[0, :, :]

bin_edges = h.axes[1].edges
bin_contents = h.values().flatten()
bin_uncertainties = np.sqrt(h.variances()).flatten()

th1d = create_th1d_histogram(HIST_NAME, bin_edges, bin_contents, bin_uncertainties)
th1d.SetMinimum(0)
th1d.SetMaximum(1.6)

# -------------------------------
# Save Histogram to ROOT File
# -------------------------------
output_root_file = os.path.join(OUTPUT_DIR, f"{CORRECTION_TYPE}_{CATEGORY}.root")
output_file = ROOT.TFile(output_root_file, "RECREATE")
th1d.Write(HIST_NAME)

# -------------------------------
# Fit Fake Factor & Save Fit Outputs
# -------------------------------
fit_result, h_uncert, ratio_hist, _ = fit_fake_factor(th1d, *PT_RANGE, usePol1=False, polOnly=3)
fit_result.Write("fit_result")
h_uncert.Write("h_uncert")
ratio_hist.Write("ratio_hist")

# -------------------------------
# Plot Fit Results
# -------------------------------
def create_canvas():
    """Create and configure a ROOT TCanvas for plotting."""
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
    
    fit_result.Draw("SAME")
    fit_result.SetLineColor(ROOT.kAzure + 7)
    
    legend = ROOT.TLegend(0.15, 0.75, 0.35, 0.9)
    legend.AddEntry(ratio_hist, "Extrapolation Correction", "EP")
    legend.AddEntry(fit_result, "Fit Result", "L")
    legend.AddEntry(h_uncert, "68% CL (Uncertainties)", "F")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)
    legend.Draw()
    
    return canvas

# Create and save the full plot
uncert_canvas = create_canvas()

output_image_base = output_root_file.replace(".root", "")
uncert_canvas.SaveAs(f"{output_image_base}_FullPlot.png")
uncert_canvas.SaveAs(f"{output_image_base}_FullPlot.pdf")

# -------------------------------
# Adjust Fit Box and Save Details
# -------------------------------
ROOT.gStyle.SetOptStat(0)  # Disable stat box
ROOT.gStyle.SetStatY(0.4)
ROOT.gStyle.SetStatX(0.7)
ROOT.gStyle.SetStatW(0.15)
ROOT.gStyle.SetStatH(0.15)
ROOT.gStyle.SetOptFit(1)   # Show fit info


uncert_canvas.SaveAs(f"{output_image_base}_FitDetails.png")
uncert_canvas.SaveAs(f"{output_image_base}_FitDetails.pdf")
uncert_canvas.Write("Fullplot")

output_file.Close()
print(f"Output saved to {output_root_file}")

# -------------------------------
# Save Correction to CorrectionLib
# -------------------------------

# Save the correction to the CorrectionLib
#### JSON FILE ####

output_json_file = output_root_file.replace(".root", ".json")

# Save the fake factor to a JSON file

fit_formula = str(fit_result.GetExpFormula("P"))  # Explicitly cast to a Python string


fit_up, fit_down = get_variated_function(fit_result)

fit_up_formula = str(fit_up.GetExpFormula("P"))
fit_down_formula = str(fit_down.GetExpFormula("P"))

print(f"Fit formula: {fit_formula}")
print(f" fit up formula: {fit_up_formula}")
print(f" fit down formula: {fit_down_formula}")


save_to_correctionlib_with_fit(ratio_hist, output_json_file, -1, -1, fit_formula, fit_up_formula, fit_down_formula, CORRECTION_TYPE, VARIABLE, PT_RANGE[0], PT_RANGE[1])

