# import pickle
# import uproot
# import ROOT
# import numpy as np
# import hist
# from hist import Hist
# import matplotlib.pyplot as plt
# import aghast



# # Step 1: Load data from pickle file
# pickle_file_path = '/eos/user/o/oponcet/code/analysis/CP/cf_store/analysis_httcp/cf.CreateHistograms/run2_UL2017_nano_tau_v10_limited/data_tau_b/nominal/calib__main/sel__main_FF/prod__main_FF/DR1_DR2_test2/histograms__vars_tau_1_pt__0.pickle'
# with open(pickle_file_path, 'rb') as f:
#     hist_data = pickle.load(f)

# hist_data = hist_data['tau_1_pt']

# print(f"data : {hist_data}")

# hist = hist_data.profile(axis=0)
# print(f"hist : {hist}")

# the_hist = hist[0,0,:] # Select only data hist withoud shift

# # print(f"the_hist : {the_hist}")

# # Extracting data from the_hist
# # print(f"the_hist.axes : {the_hist.axes}")
# # print(f"the_hist.axes[0].edges[:-1] : {the_hist.axes[0].edges[:-1]}")
# # bins = the_hist.axes[0].edges[:-1]  # bin edges
# # values = the_hist.values()          # bin values
# # print(f"values : {values}")

# print(">>>>>>>>>>>>>>>>>>")
# np_array = the_hist.view()
# print(f"np_array : {np_array}")


# # Extract values for plotting
# values = np_array['value']
# counts = np_array['count']




# # Plotting
# plt.figure(figsize=(10, 6))
# plt.bar(np.arange(len(np_array)), counts, align='center', alpha=0.5)
# plt.xlabel('Bins')
# plt.ylabel('Counts')
# plt.title('Histogram Counts')
# plt.xticks(np.arange(len(np_array)))
# plt.grid(True)
# plt.tight_layout()

# Display the plot
# plt.show()



# # 
# ghastly_hist = aghast.from_numpy(numpy_hist)

# root_hist = aghast.to_root(ghastly_hist, "root_hist")


# # # Save the_hist using pickle (optional)
# # output_file_path = 'the_hist.pickle'
# # with open(output_file_path, 'wb') as f:
# #     pickle.dump(the_hist, f)

# # print(f"Saved the_hist to {output_file_path}")




import pickle
import hist
import aghast
import uproot

# Path to your pickle file
pickle_file_path = '/eos/user/o/oponcet/code/analysis/CP/cf_store/analysis_httcp/cf.CreateHistograms/run2_UL2017_nano_tau_v10_limited/data_tau_b/nominal/calib__main/sel__main_FF/prod__main_FF/DR1_DR2_test2/histograms__vars_tau_1_pt__0.pickle'

# Load the histogram data from pickle
with open(pickle_file_path, 'rb') as f:
    hist_data = pickle.load(f)

# Extract the desired histogram (assuming 'tau_1_pt' is the key)
histogram = hist_data['tau_1_pt']

print(f"data : {hist_data}")


# Output ROOT file path
root_file_path = '/eos/user/o/oponcet/code/analysis/CP/cf_store/analysis_httcp/cf.CreateHistograms/run2_UL2017_nano_tau_v10_limited/data_tau_b/nominal/calib__main/sel__main_FF/prod__main_FF/DR1_DR2_test2/rootfile.root'

# Create a ROOT file and save histograms using aghast and uproot
with uproot.recreate(root_file_path) as root_file:
    # Iterate over the axes to get categories
    for category_axis in histogram.axes:
        print(f"category_axis : {category_axis}") # IntCategory([1, 100, 101, 301, 509, 510, 201], growth=True, name='category')
        for category in category_axis:
            print(f"category : {category}")
            print("num_categories = len(axis)=", len(category_axis))
            # Extract histogram for the current category
            hist_category = histogram[hist.loc(category),:,:,:]

            print(f"hist_category : {hist_category}")


            h = hist_category.to_numpy()
            print(f"h : {h}")

            hist_values = h[0][0].flatten()
            print(f"hist_values : {hist_values}")
            bin_edges = h[3]
            print(f"bin_edges : {bin_edges}")
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            numpy_hist = np.histogram(bin_edges[:-1], bins=bin_edges, weights=hist_values)
            print(f"numpy_hist : {numpy_hist}")



            # Plotting
            plt.figure(figsize=(10, 6))
            plt.bar(bin_edges[:-1], hist_values, width=np.diff(bin_edges), align='edge', edgecolor='black')
            plt.xlabel('Bins')
            plt.ylabel('Counts')
            plt.title('Histogram Counts')
            plt.grid(True)
            plt.tight_layout()
            # plt.show()


            # print type of hist_category
            print(f"type(h) : {type(h)}")
            print(f"type(hist_category) : {type(hist_category)}")
            print(f"type(numpy_hist) : {type(numpy_hist)}")


            # Convert hist_category to aghast histogram
            # ghastly_hist = aghast.from_numpy(h)
            aghast_hist = aghast.from_numpy(numpy_hist)
            

            print(f"ghastly_hist : {ghastly_hist}")

            # file = ROOT.TFile("demo_root_file.root", "RECREATE")
            # root_hist = aghast.to_root(ghastly_hist, "root_hist")
            # file.Write()

            # # Save to ROOT file
            # root_file[f"{category_axis.name}_{category.label}"] = hist_aghast.to_rootio()


print(f"Histograms saved to {root_file_path}")
