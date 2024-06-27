import pickle
import hist
import uproot
import numpy as np
import matplotlib.pyplot as plt

import ROOT
ROOT.gROOT.SetBatch(True)
# Path to your pickle file
pickle_file_path = '/eos/user/o/oponcet/code/analysis/CP/cf_store/analysis_httcp/cf.CreateHistograms/run2_UL2017_nano_tau_v10_limited/data_tau_b/nominal/calib__main/sel__main_FF/prod__main_FF/DRiso_DRantiIso_testWorking1/histograms__vars_tau_1_pt__0.pickle'

# Load the histogram data from pickle
with open(pickle_file_path, 'rb') as f:
    hist_data = pickle.load(f)

# Extract the desired histogram (assuming 'tau_1_pt' is the key)
histogram = hist_data['tau_1_pt']


# Iterate over the axes to get categories
category_axis = histogram.axes[0]
print(f"category_axis : {category_axis}") # IntCategory([1, 100, 101, 301, 509, 510, 201], growth=True, name='category')
for category in category_axis:
    print(f"category : {category}")
    # Extract histogram for the current category
    hist_category = histogram[hist.loc(category),:,:,:]

    # Convert hist_category to numpy arrays
    h = hist_category.to_numpy()
    hist_values = h[0][0].flatten()
    bin_edges = h[3]

    # Create numpy histogram
    # numpy_hist = np.histogram(bin_edges[:-1], bins=bin_edges, weights=hist_values)
    numpy_hist = (hist_values, bin_edges)
    print(f"numpy_hist : {numpy_hist}")



    # Extract the number of bins, bin edges, and histogram values
    num_bins = len(hist_values)
    bin_edges = np.array(bin_edges)
    hist_values = np.array(hist_values)

    # Create a ROOT histogram
    # ROOT.TH1F(name, title, number of bins, array of bin edges)
    root_hist = ROOT.TH1F("hist_name", "hist_title", num_bins, bin_edges)

    # Fill the ROOT histogram with the values from the numpy histogram
    for i in range(num_bins):
        root_hist.SetBinContent(i + 1, hist_values[i])

    # Optionally, you can draw the histogram to visualize it
    c1 = ROOT.TCanvas("c1", "Canvas", 800, 600)
    # root_hist.Draw()
    root_hist.SaveAs(f"histogram_cat{category}_tau_1_pt.root")  # Save the histogram as a ROOT file
    c1.SaveAs("histogram.png")  # Save the histogram as a PNG file
    

    # # Convert numpy histogram to aghast histogram
    # aghast_hist = aghast.from_numpy(numpy_hist)
    # aghast_hist.dump()
    # print(f"aghast_hist : {aghast_hist}")


    # Plotting (example)
    plt.figure(figsize=(10, 6))
    plt.bar(bin_edges[:-1], hist_values, width=np.diff(bin_edges), align='edge', edgecolor='black')
    plt.xlabel('Bins')
    plt.ylabel('Counts')
    plt.title('Histogram Counts')
    plt.grid(True)
    plt.tight_layout()
    # plt.show()

    # # Save to ROOT file
    # root_hist = aghast.to_root(aghast_hist, "root_hist")
    
    
    # canvas = ROOT.TCanvas()
    # root_hist.Draw()
    # canvas.Draw()

    # Save to ROOT file

        



            # root_file[f"{category_axis.name}_{category.label}"] = aghast.to_root(aghast_hist, f"{category_axis.name}_{category.label}")

# print(f"Histograms saved to {root_file_path}")
