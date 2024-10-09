'''
Script to calculate the fake factor and plot the data-MC subtraction and the ratio of two histograms
Date : 2021-09-30
Author : @oponcet
'''

import pickle
import hist
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import scinum as sn
import hist
from hist import Hist
import matplotlib.cm as cm  # Importing color maps



def hist_to_num(h: hist.hist, unc_name=str(sn.DEFAULT)) -> sn.Number:
        '''
        Return an `sn.Number` object, where:
        - The nominal values (i.e., bin contents) are `h.values()`.
        - The uncertainties are stored in a dictionary where the key is `unc_name` and the value 
          is the square root of the variances (h.variances()**0.5).
        '''
        return sn.Number(h.values(), {unc_name: h.variances()**0.5})
        return sn.Number(h.values(), {unc_name: h.variances()**0.5})



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


def get_all_hist():
    ''' 
    Load all histograms from the pickle files for data and MC samples
    Returns: dictionnary of histograms for data and MC samples
    '''

    ## Define the paths to the pickle files
    eos_path = "/eos/user/o/oponcet2/analysis/FakeFactor/analysis_httcp/"
    task = "cf.MergeHistograms/run3_2022_preEE_nano_cp_tau_v12/"
    dataset_data = ["data_tau_C", "data_tau_D"]
    dataset_mc = ["dy_lep_m50", "h_ggf_tautau_prod_cp_even_sm", "st_tchannel_t", "st_tchannel_tbar", "st_tw_tb_dl", "st_tw_tb_fh", "st_tw_tb_sl", "st_tw_t_dl", "st_tw_t_fh", "st_tw_t_sl", "tt_dl", "tt_fh", "tt_sl", "wj_incl", "ww", "wz", "zz"]
    hist_path = "nominal/calib__main/sel__main_FF/prod__main_FF/weight__all_weights/htcondor2/hist__tau_1_pt.pickle"

    ## Load data histogram from the pickle file
    hists_data = {}
    for dataset in dataset_data:
        pickle_file_path = eos_path + task + dataset + "/" + hist_path
        # print(f"Loading histogram data from: {pickle_file_path}")
        hists_data[dataset] = load_histogram_data(pickle_file_path)
        # print(f"hist_data : {hists_data[dataset]}")

    ## Load MC histogram from the pickle file
    hists_mc = {}
    for dataset in dataset_mc:
        pickle_file_path = eos_path + task + dataset + "/" + hist_path
        # print(f"Loading histogram data from: {pickle_file_path}")
        hists_mc[dataset] = load_histogram_data(pickle_file_path)
        # print(f"hist_data : {hists_mc[dataset]}")

    # print(f">>>>>>>>>>>>>>>>>>>>> Data histograms: \n{hists_data}")

    # print(f">>>>>>>>>>>>>>>>>>>>> MC histograms: \n{hists_mc}")

    return hists_data, hists_mc

def get_hist_by_category_var(hists_data, hists_mc, category):
    '''
    Get histograms for a specific category and variable
    Returns: dictionnary of histograms for data and MC samples by category
    '''

    # Initialize the dictionaries to store the histograms for the category
    hists_data_cat = {}
    hists_mc_cat = {}

    # Loop over the datasets to get the histograms for the category
    for dataset in hists_data.keys():
        category_axis = hists_data[dataset].axes['category']
        if category not in category_axis:
            print(f"Error: Category {category} not found in the histogram data for dataset {dataset}")
            continue
        category_index = category_axis.index(category)
        hists_data_cat[dataset] = hists_data[dataset][category_index, :, :, :].project('tau_1_pt')
    
    # Loop over the datasets to get the histograms for the category
    for dataset in hists_mc.keys():
        category_axis = hists_mc[dataset].axes['category']
        if category not in category_axis:
            print(f"Error: Category {category} not found in the histogram data for dataset {dataset}")
            continue
        category_index = category_axis.index(category)
        hists_mc_cat[dataset]  = hists_mc[dataset][category_index, :, :, :].project('tau_1_pt')

    return hists_data_cat, hists_mc_cat


def convert_all_hist_tonumpy(hists_data, hists_mc, variable,category):
    '''
    Convert all histograms to numpy arrays
    Returns: dictionnary of histograms for data and MC samples in numpy arrays
    '''

    hists_data_np = {}
    hists_data_variance_np = {}
    hists_mc_np = {}
    hists_mc_variance_np = {}


    # Loop over the datasets to convert the histograms to numpy arrays
    for dataset in hists_data.keys():
        category_axis = hists_data[dataset].axes['category']
        if category not in category_axis:
            print(f"Error: Category {category} not found in the histogram data for dataset {dataset}")
            continue
        category_index = category_axis.index(category)
        hists_data_np[dataset] = hists_data[dataset][category_index, :, :, :].project('tau_1_pt').values()
        hists_data_variance_np[dataset] = hists_data[dataset][category_index, :, :, :].project('tau_1_pt').variances()**0.5
    
    # Loop over the datasets to convert the histograms to numpy arrays
    for dataset in hists_mc.keys():
        category_axis = hists_mc[dataset].axes['category']
        if category not in category_axis:
            print(f"Error: Category {category} not found in the histogram data for dataset {dataset}")
            continue
        category_index = category_axis.index(category)
        hists_mc_np[dataset] = hists_mc[dataset][category_index, :, :, :].project('tau_1_pt').values()
        hists_mc_variance_np[dataset] = hists_mc[dataset][category_index, :, :, :].project('tau_1_pt').variances()**0.5

    # get varaible edges : 
    variable_bin_edges = hists_data['data_tau_C'].axes[variable].edges
    variable_centers = (variable_bin_edges[:-1] + variable_bin_edges[1:]) / 2  # Get the bin centers for plotting


    # print(f">>>>>>>>>>>>>>>>>>>>> Data histograms in numpy: \n{hists_data_np}")
    # print(f">>>>>>>>>>>>>>>>>>>>> Data histograms variance in numpy: \n{hists_data_variance_np}")
    # print(f">>>>>>>>>>>>>>>>>>>>> MC histograms in numpy: \n{hists_mc_np}")
    # print(f">>>>>>>>>>>>>>>>>>>>> MC histograms variance in numpy: \n{hists_mc_variance_np}")

    # print(f">>>>>>>>>>>>>>>>>>>>> Variable bin edges: \n{variable_bin_edges}")
    # print(f">>>>>>>>>>>>>>>>>>>>> Variable centers: \n{variable_centers}")

    return hists_data_np, hists_mc_np, hists_data_variance_np, hists_mc_variance_np, variable_bin_edges, variable_centers


def plot_hist_cat(hists_data_cat, hists_mc_cat, cat):
    '''
    Plot histograms for a specific category
    '''

    hep.style.use(hep.style.CMS)

    # Initialize the sum variable as None
    hists_data_cat_sum = None

    # Loop over the datasets to sum histograms
    for dataset in hists_data_cat:
        print(dataset)
        hist = hists_data_cat[dataset]

        # Initialize the sum histogram if it's None
        if hists_data_cat_sum is None:
            hists_data_cat_sum = hist  # Start with the first histogram
        else:
            hists_data_cat_sum += hist  # Accumulate the histograms

    # Prepare the figure and axis for plotting
    fig, ax = plt.subplots()

    # Use a color map to dynamically generate colors based on the number of MC categories
    num_mc_categories = len(hists_mc_cat)
    cmap = cm.get_cmap('tab10', num_mc_categories)  # Use 'tab10', 'viridis', or any other colormap

    # Plot stacked MC histograms
    mc_total = None  # This will store the total MC histogram
    mc_labels = list(hists_mc_cat.keys())

    for idx, dataset in enumerate(mc_labels):
        mc_hist = hists_mc_cat[dataset]

        # Ensure no negative values and only work with non-negative values
        mc_values = np.maximum(mc_hist.values(), 0)  # Ensure MC values are non-negative

        if mc_total is None:
            mc_total = mc_values.copy()  # Start with the first histogram
        else:
            mc_total += mc_values  # Accumulate the histograms
        
        # Assign color from the color map
        color = cmap(idx / num_mc_categories)

        # Plot the stacked MC histogram for each dataset
        ax.bar(mc_hist.axes[0].centers, mc_values, width=np.diff(mc_hist.axes[0].edges), 
               label=dataset, align='center', color=color, alpha=0.7, edgecolor='black')

    # Plot the data histogram with error bars (as points with errors)
    data_hist = hists_data_cat_sum
    data_values = data_hist.values()
    data_variances = data_hist.variances()

    # Ensure no negative data values
    data_values = np.maximum(data_values, 0)  # Ensure data values are non-negative

    # Plot data as points with error bars
    ax.errorbar(data_hist.axes[0].centers, data_values, yerr=np.sqrt(data_variances), fmt='o', 
                color='black', label='Data', capsize=3, markersize=5)

    # Add labels, title, and legend
    ax.set_xlabel("Tau $p_T$ / GeV")
    ax.set_ylabel("Events")
    ax.set_title(f"Stacked MC vs Data - Category {cat}")
    ax.legend()

    # Show the plot and save it
    plt.tight_layout()
    plt.savefig(f"script_FF/plots/stacked_plot_cat_{cat}.png", dpi=140)
    # plt.show()


def calculate_data_MC_substraction(hists_data_cat, hists_mc_cat, cat):
    '''
    Calculate the data-MC subtraction for a specific category and plot the resulting histogram
    Returns: histogram of the data-MC subtraction
    '''

    hists_data_cat_sum = None
    hists_mc_cat_sum = None

    print(hists_data_cat)

    # from IPython import embed; embed()

    for dataset in hists_data_cat:
        hist = hists_data_cat[dataset]
        if hists_data_cat_sum is None:
            hists_data_cat_sum = hist
        else:
            hists_data_cat_sum += hist

    for dataset in hists_mc_cat:
        hist = hists_mc_cat[dataset]
        if hists_mc_cat_sum is None:
            hists_mc_cat_sum = hist
        else:
            hists_mc_cat_sum += hist
            
    print(f">>>>>>>>>>>>>>>>>>>>> Data histograms sum: \n{hists_data_cat_sum}")
    print(f">>>>>>>>>>>>>>>>>>>>> MC histograms sum: \n{hists_mc_cat_sum}")

    # Guarantee non-negative values
    hists_data_cat_sum.view().value = np.maximum(hists_data_cat_sum.values(), 0)
    hists_mc_cat_sum.view().value = np.maximum(hists_mc_cat_sum.values(), 0)

    print(f">>>>>>>>>>>>>>>>>>>>> Data histograms sum after non-0: \n{hists_data_cat_sum}")
    print(f">>>>>>>>>>>>>>>>>>>>> MC histograms sum after non-0: \n{hists_mc_cat_sum}")


    hists_data_cat_sum_num = hist_to_num(hists_data_cat_sum)
    hists_mc_cat_sum_num = hist_to_num(hists_mc_cat_sum)
    

    # Subtract MC from data
    hists_data_mc_sub_num = hists_data_cat_sum_num - hists_mc_cat_sum_num


    # combine uncertainties and store values in bare arrays
    hists_data_mc_sub_values = hists_data_mc_sub_num()
    hists_data_mc_sub_variances = hists_data_mc_sub_num(sn.UP, sn.ALL, unc=True)**2

    print(f">>>>>>>>>>>>>>>>>>>>> Data - MC subtraction values: \n{hists_data_mc_sub_values}")
    print(f">>>>>>>>>>>>>>>>>>>>> Data - MC subtraction variances: \n{hists_data_mc_sub_variances}")    

    # # Create hist from the subtraction
    # hists_data_mc_sub_hist = hist.hist(hists_data_mc_sub_num.axes[0], hists_data_mc_sub_values, hists_data_mc_sub_variances)



    # Create histo copy of the subtraction
    hists_data_mc_sub = hists_mc_cat_sum.copy().reset()

    print(f">>>>>>>>>>>>>>>>>>>>> Data - MC subtraction hist: \n{hists_data_mc_sub}")

    hists_data_mc_sub.view().value[...] = hists_data_mc_sub_values
    hists_data_mc_sub.view().variance[...] = hists_data_mc_sub_variances

    # Ensure non-negative values
    hists_data_mc_sub.view().value = np.maximum(hists_data_mc_sub.values(), 0)


    print(f">>>>>>>>>>>>>>>>>>>>> Data - MC subtraction hist: \n{hists_data_mc_sub}")
    

    # Plot the data-MC subtraction
    fig, ax = plt.subplots()
    ax.bar(hists_data_mc_sub.axes[0].centers, hists_data_mc_sub.values(), width=np.diff(hists_data_mc_sub.axes[0].edges))
    ax.set_ylabel("Events")
    ax.set_title(f"Data - MC Subtraction - Category {cat}")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"script_FF/plots/data_mc_subtraction_cat_{cat}.png", dpi=140)
    # plt.show()
    return hists_data_mc_sub


def calculate_ratio_hist(hist_iso, hist_antiiso):
    '''
    Calculate the ratio of two histograms and plot the resulting histogram
    Returns: ratio values and variances
    '''

    # Calculate the ratio of two histograms

    hist_iso_num = hist_to_num(hist_iso)
    hist_antiiso_num = hist_to_num(hist_antiiso)


    ratio = hist_iso_num / hist_antiiso_num

    # Create hist ratio with same x range as the hist iso
    ratio_hist = hist_iso.copy().reset()

    # Get the ratio values and variances
    ratio_values = ratio()
    ratio_variances = ratio(sn.UP, sn.ALL, unc=True)**2

    # Fill the ratio hist
    ratio_hist.view().value[...] = ratio_values
    ratio_hist.view().variance[...] = ratio_variances


    # Plot the ratio
    fig, ax = plt.subplots()
    ax.errorbar(ratio_hist.axes[0].centers, ratio_hist.values(), yerr=np.sqrt(ratio_hist.variances()), fmt='o', label='Ratio', color='black')
    ax.set_ylabel("Ratio")
    ax.set_title("Ratio of Iso / AntiIso")
    ax.legend()
    plt.savefig("script_FF/plots/ratio_hist.png", dpi=140)


    return ratio_values, ratio_variances

if __name__ == "__main__":
    '''
    Main function to run the script
    '''

    hists_data, hists_mc =  get_all_hist()



    ## plot all histograms by category
    category = [600,500] # FFDRantiIso_tautau = 600, FFDRIso_tautau = 500, category that i want 

    hists_mc_data_sub = {}
    for cat in category:
        hists_data_cat, hists_mc_cat = get_hist_by_category_var(hists_data, hists_mc, cat)
        plot_hist_cat(hists_data_cat, hists_mc_cat,cat)
        hists_mc_data_sub[cat] = calculate_data_MC_substraction(hists_data_cat, hists_mc_cat, cat)

    calculate_ratio_hist(hists_mc_data_sub[500], hists_mc_data_sub[600])



