import pickle
import hist
import numpy as np
import matplotlib.pyplot as plt
import mplhep


def load_histogram_data(pickle_file_path):
    """
    Load histogram data from a specified pickle file.

    Parameters:
    pickle_file_path (str): The path to the pickle file.

    Returns:
    dict or object: The loaded histogram data from the pickle file.

    Raises:
    FileNotFoundError: If the specified file is not found.
    Exception: For other errors during the loading process.
    """
    try:
        with open(pickle_file_path, 'rb') as f:
            hist_data = pickle.load(f)
            print("Histogram data loaded successfully!")
            return hist_data
    except FileNotFoundError:
        print(f"Error: The file '{pickle_file_path}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred while loading the file: {e}")
        raise


def plot_histogram(hist_data, variable_name):
    """
    Plot the histogram data using matplotlib.

    Parameters:
    hist_data (Hist): The histogram data object.
    variable_name (str): The name of the variable for the plot title.
    """
    # Extract the Variable object from hist_data

    # # Extract the histogram bins and weights
    # bins = hist_data.axes[0].edges  # Bin edges
    # weights = hist_data.values  # Histogram weights

    print(f"hist_data.axes[1] :  {hist_data.axes[1]}")


    # Convert to numpy arrays for plotting
    values = hist_data.values()  # returns an array where each element corresponds to the weight or count in a specific bin of the histogram. This is the data that represents the frequency or the weighted count of observations in each bin.
    # errors = np.sqrt(hist_data.values().sum(axis=0)) # Calculates the errors (standard deviations) for each bin.
    edges = hist_data.axes[1].edges  # Bin edges
    variances = hist_data.variances() 

    values = hist_data.values()
    print(f"values : {values}")


    # # If you want to get the sum of weights
    # weight_sum = values.Weight()
    # print("Sum of weights:", weight_sum)

    

    print(f"values  : {values}")
    # print(f"errors : {errors}")
    print(f"edges  : {edges}")
    print(f"variances  : {variances}")

    # Calculate the bin centers
    bin_centers = (edges[:-1] + edges[1:]) / 2
    print(f"bin_centers : {bin_centers}")


    # Calculate Poisson intervals
    yerr = poisson_interval(values, variances)
    yerr[np.isnan(yerr)] = 0
    yerr[0] = values - yerr[0]
    yerr[1] -= values
    yerr[yerr < 0] = 0

    
    # # Plot the histogram
    # plt.figure(figsize=(10, 6))
    # plt.errorbar(edges[:-1], counts, yerr=errors, fmt='o', label='Histogram Data', capsize=5)
    # plt.title(f'Histogram of {variable_name}')
    # plt.xlabel(variable_name)
    # plt.ylabel('Counts')
    # plt.grid(True)
    # plt.show()





mplhep.style.use("CMS")

# Load histogram data from the pickle file
cf_store_path = "/eos/user/o/oponcet2/analysis/FakeFactor/analysis_httcp"

variable = "tau_1_pt"
category = "FFDRantiIso_tautau"
print(f"Category Name: {category}")
print(f"Variable Name: {variable}")

pickle_file_path = cf_store_path + '/cf.DataDrivenEstimation/run3_2022_preEE_nano_cp_tau_v12/nominal/calib__main/sel__main_FF/prod__main_FF/datasets_19_130ad08fc8/htcondor2/qcd_histogram__'+category+'_'+variable+'.pickle'
print(f"Loading histogram data from: {pickle_file_path}")

hist_data = load_histogram_data(pickle_file_path)
print(f"hist_data : {hist_data}")

#print type of hist_data
print(f"Type of hist_data : {type(hist_data)}")

# # Plot the histogram data
# hist_data.plot1d(hist_data)

# Plot the histogram data
fig, ax = plt.subplots()
# hist.plot1d(ax=ax, label='Histogram', color='blue')

hist_data.plot1d(ax=ax, color="black", lw=3, ls="-")

# Add labels and title
plt.xlabel('Tau $p_{T}$ / GeV')
plt.ylabel('Events')
plt.title(f'Histogram of {variable} for {category}')
# plt.show()

plt.savefig(f"script_FF/plots/plot_{category}_{variable}.png", format='png')

    # # Plot the histogram for the variable
    # if hist_data:
    #     plot_histogram(hist_data, variable)

# Extract the desired histogram (assuming 'tau_1_pt' is the key)
# histogram = hist_data[variable]


# # Iterate over the axes to get categories
# category_axis = histogram.axes[0]
# print(f"category_axis : {category_axis}") # IntCategory([1, 100, 101, 301, 509, 510, 201], growth=True, name='category')
# for cat in category_axis:
#     print(f"category : {cat}")
#     # Extract histogram for the current category
#     hist_category = histogram[hist.loc(category),:,:,:]

#     # Convert hist_category to numpy arrays
#     h = hist_category.to_numpy()
#     hist_values = h[0][0].flatten()
#     bin_edges = h[3]

#     # Create numpy histogram
#     # numpy_hist = np.histogram(bin_edges[:-1], bins=bin_edges, weights=hist_values)
#     numpy_hist = (hist_values, bin_edges)
#     print(f"numpy_hist : {numpy_hist}")



#     # Extract the number of bins, bin edges, and histogram values
#     num_bins = len(hist_values)
#     bin_edges = np.array(bin_edges)
#     hist_values = np.array(hist_values)

#     # Plotting (example)
#     plt.figure(figsize=(10, 6))
#     plt.bar(bin_edges[:-1], hist_values, width=np.diff(bin_edges), align='edge', edgecolor='black')
#     plt.xlabel('Bins')
#     plt.ylabel('Counts')
#     plt.title('Histogram Counts')
#     plt.grid(True)
#     plt.tight_layout()
#     # plt.show()
