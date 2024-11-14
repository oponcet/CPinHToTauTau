# import pickle
# import uproot
# import numpy as np
# from hist import Hist
# from hist.axis import IntCategory, StrCategory, Variable
# import ROOT

# # Step 1: Load data from pickle file
# pickle_file_path = '/eos/user/o/oponcet/code/analysis/CP/cf_store/analysis_httcp/cf.CreateHistograms/run2_UL2017_nano_tau_v10_limited/data_mu_b/nominal/calib__main/sel__main_FF/prod__main_FF/DR1_test4/histograms__vars_tau_1_pt__0.pickle'
# with open(pickle_file_path, 'rb') as f:
#     hist_data = pickle.load(f)

# print(f"data : {hist_data}")

# # Ensure that hist_data is of the correct type
# if not isinstance(hist_data, Hist):
#     raise TypeError("Loaded data is not of type 'Hist'")

# # Get the number of bins
# num_bins = hist_data.axes[0].size
# print(f"Number of bins: {num_bins}")

# # Convert Hist to a NumPy array
# bins = np.array([axis.centers for axis in hist_data.axes])
# weights = np.ravel(hist_data.values(flow=False))
# print(f"weight of bin : {weights}")

# # Create a ROOT histogram
# root_hist = ROOT.TH1D("histogram", "Histogram", num_bins, 0, num_bins)

# # Fill the ROOT histogram
# for bin_idx in range(num_bins):
#     root_hist.SetBinContent(bin_idx + 1, weights[bin_idx])

# # Write the ROOT histogram to a ROOT file
# output_file_path = "output_file.root"
# with ROOT.TFile(output_file_path, "RECREATE") as f:
#     root_hist.Write()

# print(f"Histogram has been successfully converted and saved to {output_file_path}")
