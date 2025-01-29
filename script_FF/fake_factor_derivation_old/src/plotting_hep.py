
import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

def plot_stack_from_root_file(root_file_path, output_dir):
    # Open the ROOT file
    root_file = ROOT.TFile(root_file_path, "READ")
    
    # Retrieve the canvas
    canvas = root_file.Get("canvas")  # Ensure you know the correct canvas name
    if not canvas:
        print("Canvas not found in the ROOT file.")
        return

    # Retrieve the stack (assuming it's in the canvas)
    mc_stack = canvas.GetListOfPrimitives()  # This contains the stack (and possibly other objects)
    
    histograms = []
    # Loop through primitives to get the TH1D histograms from the stack
    for obj in mc_stack:
        if isinstance(obj, ROOT.TH1D):
            histograms.append(obj)

    if not histograms:
        print("No TH1D histograms found in the canvas!")
        return

    # Set up the CMS style using mplhep
    hep.style.use("CMS")

    # Create a figure using matplotlib for plotting
    fig, ax = plt.subplots(figsize=(8, 6))

    # Prepare data for the stacked plot
    bin_edges = [histograms[0].GetBinLowEdge(i) for i in range(1, histograms[0].GetNbinsX() + 2)]

    # Convert ROOT histograms to numpy arrays for matplotlib/ mplhep
    stack_data = []
    for hist in histograms:
        hist_data = np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)])
        stack_data.append(hist_data)

    # Stack and plot histograms
    for i, hist_data in enumerate(stack_data):
        if i == 0:
            ax.bar(bin_edges[:len(hist_data)], hist_data, label=histograms[i].GetName(), color='C0', zorder=3)
        else:
            ax.bar(bin_edges[:len(hist_data)], hist_data, label=histograms[i].GetName(),
                   bottom=np.sum(stack_data[:i], axis=0), color=f"C{i}", zorder=3)

    # Add labels, title, and legend
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Events")
    ax.set_title("Stacked Histogram")
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))

    # Add CMS style
    # hep.cms.lumitext(ax, lumi=35, year=2018)

    # Optional: Add a grid
    ax.grid(True)

    # Show the plot
    plt.tight_layout()
    plt.savefig(f"{output_dir}/stacked_plot.png")

    # Close the ROOT file after plotting
    root_file.Close()
