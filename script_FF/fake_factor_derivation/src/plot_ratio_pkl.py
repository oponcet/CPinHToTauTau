import ROOT
import pickle
import numpy as np
import matplotlib.pyplot as plt
import scinum as sn
import hist
import cmsstyle as CMS
import os
from fit_functions import fit_fake_factor


# dont open the canvas
ROOT.gROOT.SetBatch(True)

def create_th1d_histogram(var, bin_edges, bin_contents, bin_uncertainties):
    """
    Create a ROOT TH1D histogram from bin contents and uncertainties.
    """
    dataset = "dataset"
    n_bins = len(bin_contents)
    th1d = ROOT.TH1D(dataset, var, n_bins, bin_edges[0], bin_edges[-1])
    
    for i, (content, uncertainty) in enumerate(zip(bin_contents, bin_uncertainties), start=1):
        th1d.SetBinContent(i, content)
        th1d.SetBinError(i, uncertainty)

    return th1d


# File path to the .pkl file
file_path = "/eos/user/g/gsaha/CPinHToTauTauOutput/cf_store/analysis_httcp/cf.PlotVariables1D/QCD0/Ratio/RATIO_met_var_qcd_h1_inclusive.pkl"

# Load the scinum.Number object
with open(file_path, "rb") as file:
    ratio = pickle.load(file)

# Inspect the structure
print("Loaded object:", ratio)
print("Type:", type(ratio))


category_axis = ratio.axes['category']
# category_index = category_axis.index(410009000)
category_index = 0

h = ratio[category_index, :, :]


# convert to ROOT histogram TH1D 

bin_edges = h.axes[1].edges
print("bin_edges:", bin_edges)
bin_contents = h.values().flatten()
print("bin_contents:", bin_contents)
bin_uncertainties = (h.variances()).flatten()
print("bin_uncertainties:", bin_uncertainties)
var = "met_var_qcd_h1"
th1d = create_th1d_histogram(var, bin_edges, bin_contents, bin_uncertainties)

# set thi1 y axis range 0 to 1.5
th1d.SetMinimum(0)
th1d.SetMaximum(1.6)


# SAve the th1d in a root file 
output_floder = "script_FF/fake_factor_derivation/outputs/extrapolation_correction/"
# Ensure the output folder exists
os.makedirs(output_floder, exist_ok=True)

output_file_path = output_floder + "extrapolation_correction_inclusive.root"
output_file = ROOT.TFile(output_file_path, "RECREATE")
th1d.Write("ratio")

# fit the fake factor

fit_result, h_uncert, fake_factor_hist, fit_details = fit_fake_factor(th1d, -1.5, 1.5, usePol1=False, polOnly=4) # polOnly = None

# Save the fit result in a root file
fit_result.Write("fit_result")
h_uncert.Write("h_uncert")
fake_factor_hist.Write("fake_factor_hist")

#### CANVAS FULL PLOT #### 

# Plot the uncertainties
uncert_canvas = ROOT.TCanvas("Extrapolation Correction", "Extrapolation Correction", 800, 600)
    

# Draw the fake factor histogram
fake_factor_hist.Draw("EP")
fake_factor_hist.SetLineColor(ROOT.kBlack)
fake_factor_hist.SetMarkerStyle(20)
title = f"Extrapolation Correction for inclusive in DM and Njets"
fake_factor_hist.SetTitle(title)    
fake_factor_hist.GetYaxis().SetTitle("Extrapolation Correction")
fake_factor_hist.GetXaxis().SetTitle("p_{T} (GeV)")
fake_factor_hist.GetYaxis().SetRangeUser(0, 1.6)


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
legend.AddEntry(fake_factor_hist, "Extrapolation Correction", "EP")
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
output_image_path = output_file_path.replace(".root", "_FullPlot.png")
uncert_canvas.SaveAs(output_image_path)
output_image_path = output_file_path.replace(".root", "_FullPlot.pdf")
uncert_canvas.SaveAs(output_image_path)

# #### PNG FILE: FIT DETAIL #### uncert_canvas
# change fit box position and size
ROOT.gStyle.SetStatY(0.4)
ROOT.gStyle.SetStatX(0.8)
ROOT.gStyle.SetStatW(0.15)
ROOT.gStyle.SetStatH(0.15)
# Display fit info like NDF, Chi2, etc.
ROOT.gStyle.SetOptFit(1)





# Save the fit details canvas as a separate PNG file
output_details_image_path = output_file_path.replace(".root", "_FitDetails.png")
uncert_canvas.SaveAs(output_details_image_path)

# Save the fit details canvas as a separate PDF file
output_details_image_path = output_file_path.replace(".root", "_FitDetails.pdf")
uncert_canvas.SaveAs(output_details_image_path)


# Save the canvas to the ROOT file for later use
uncert_canvas.Write("Fullplot")

output_file.Close()