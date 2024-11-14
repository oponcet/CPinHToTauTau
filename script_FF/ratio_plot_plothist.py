import ROOT
import numpy as np
from boost_histogram import Histogram, axis
from plothist import plot_two_hist_comparison

ROOT.gROOT.SetBatch(True)

def ratio_plot(file1_path, file2_path, hist_name):
    # Open the ROOT files
    file1 = ROOT.TFile.Open(file1_path)
    file2 = ROOT.TFile.Open(file2_path)
    
    # Get the histograms
    hist1 = file1.Get(hist_name)
    hist2 = file2.Get(hist_name)
    
    # Check if histograms are successfully loaded
    if not hist1 or not hist2:
        print(f"Failed to load histograms {hist_name}")
        return
    
    # Convert ROOT histograms to boost_histogram objects
    def convert_root_to_boost_hist(root_hist):
        nbins = root_hist.GetNbinsX()
        edges = [root_hist.GetBinLowEdge(i + 1) for i in range(nbins + 1)]
        values = [root_hist.GetBinContent(i + 1) for i in range(nbins)]
        
        hist = Histogram(axis.Variable(edges))
        for i, value in enumerate(values):
            hist[i] = value
        
        return hist
    
    hplot1 = convert_root_to_boost_hist(hist1)
    hplot2 = convert_root_to_boost_hist(hist2)
    
    # Ensure the bin edges match
    if not np.array_equal(hplot1.axes[0].edges, hplot2.axes[0].edges):
        print("Histogram bin edges do not match!")
        return
    
    # Use plothist to create comparison plots
    fig, ax_main, ax_comparison = plot_two_hist_comparison(
        hplot1,
        hplot2,
        xlabel="leading $ tau p_{T}$ [GeV]",
        ylabel="Entries/5.0 GEV",
        h1_label="SS Iso",
        h2_label="SS anti Iso",
    )
    
    # Create the ratio plot
    ratio = np.divide(hplot1.view(), hplot2.view(), out=np.zeros_like(hplot1.view()), where=hplot2.view()!=0)
    
    # Plot the ratio on the comparison axis
    edges = hplot1.axes[0].edges
    # ax_comparison.step(edges[:-1], ratio, where='mid', color='black')
    ax_comparison.set_ylim(0, 0.2)  # Set the y-axis limits to be between 0 and 0.2
    ax_comparison.legend()
    
    # Save the figure
    fig.savefig("1d_comparison_ratio.pdf", bbox_inches="tight")
    
    # Close the ROOT files
    file1.Close()
    file2.Close()

# Example usage
file1_path = "/afs/cern.ch/user/o/oponcet/private/analysis/CPinHToTauTau/histogram_cat509_tau_1_pt.root"
file2_path = "/afs/cern.ch/user/o/oponcet/private/analysis/CPinHToTauTau/histogram_cat510_tauantiiso_1_pt.root"
hist_name = "hist_name"

ratio_plot(file1_path, file2_path, hist_name)
