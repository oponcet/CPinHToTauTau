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
import copy
import json



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
            # print("Histogram data loaded successfully!")
            return hist_data
    except FileNotFoundError:
        print(f"Error: The file '{pickle_file_path}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred while loading the file: {e}")
        raise


def get_all_hist(eos_path, task, dataset_data, dataset_mc, hist_path):
    ''' 
    Load all histograms from the pickle files for data and MC samples
    Returns: dictionnary of histograms for data and MC samples
    '''

    
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

def get_hist_by_category_var(hists_data, hists_mc, category, var = 'hcand_1_pt'):
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
            # print(f"Warning: Category {category} not found in the histogram data for dataset {dataset}")
            continue
        category_index = category_axis.index(category)
        hists_data_cat[dataset] = hists_data[dataset][category_index, :, :, :].project(var)[40j:200j:1j]

    
    # Loop over the datasets to get the histograms for the category
    for dataset in hists_mc.keys():
        category_axis = hists_mc[dataset].axes['category']
        if category not in category_axis:
            # print(f"Warning: Category {category} not found in the histogram data for dataset {dataset}")
            continue
        category_index = category_axis.index(category)
        hists_mc_cat[dataset]  = hists_mc[dataset][category_index, :, :, :].project(var)[40j:200j:1j]

    return hists_data_cat, hists_mc_cat

def get_2D_hist_by_category_var(hists_data, hists_mc, category, var_x='hcand_1_pt', var_y='hcand_2_pt'):
    '''
    Get 2D histograms for a specific category and variables.
    
    Parameters:
    hists_data : dict
        Dictionary containing histograms for data samples.
    hists_mc : dict
        Dictionary containing histograms for MC samples.
    category : str
        The category to extract histograms for.
    var_x : str
        The variable for the x-axis in the 2D histogram.
    var_y : str
        The variable for the y-axis in the 2D histogram.
        
    Returns:
    hists_data_cat : dict
        Dictionary of 2D histograms for data samples by category.
    hists_mc_cat : dict
        Dictionary of 2D histograms for MC samples by category.
    '''

    # Initialize dictionaries to store histograms for the specified category
    hists_data_cat = {}
    hists_mc_cat = {}

    # Loop over the data samples to retrieve and project histograms for the specified category
    for dataset in hists_data.keys():
        category_axis = hists_data[dataset].axes['category']
        if category not in category_axis:
            print(f"Warning: Category {category} not found in the histogram data for dataset {dataset}")
            continue
        category_index = category_axis.index(category)
        hists_data_cat[dataset] = hists_data[dataset][category_index, :, :, :, :].project(var_x, var_y)[:,40j:200j:1j]
    
    # Loop over the MC samples to retrieve and project histograms for the specified category
    for dataset in hists_mc.keys():
        category_axis = hists_mc[dataset].axes['category']
        if category not in category_axis:
            print(f"Warning: Category {category} not found in the histogram data for dataset {dataset}")
            continue
        category_index = category_axis.index(category)
        hists_mc_cat[dataset] = hists_mc[dataset][category_index, :, :, :, :].project(var_x, var_y)[:,40j:200j:1j]

    return hists_data_cat, hists_mc_cat


def plot_hist_cat(hists_data_cat_, hists_mc_cat_, cat, hists_mc_data_sub_):
    '''
    Plot histograms for a specific category, including MC uncertainties
    '''

    hep.style.use(hep.style.CMS)

    hists_data_cat = copy.deepcopy(hists_data_cat_)
    hists_mc_cat = copy.deepcopy(hists_mc_cat_)

    hists_mc_data_sub = hists_mc_data_sub_.copy() if hists_mc_data_sub_ is not None else None

    # Initialize the sum variable for the data histogram
    hists_data_cat_sum = None

    # Loop over the datasets to sum data histograms
    for dataset in hists_data_cat:
        hist = hists_data_cat[dataset]

        if hists_data_cat_sum is None:
            hists_data_cat_sum = hist  # Start with the first histogram
        else:
            hists_data_cat_sum += hist  # Accumulate the histograms

    # Prepare the figure and axis for plotting
    fig, ax = plt.subplots()

    # Define colors for specific datasets
    colors = {
        'tt': "#9e9ac8",  # Violet
        'dy': "#feb24c",  # Orange
        'diboson': "#a96b59",   # Brown for Diboson (WW)
        'wj': "#d73027", # Red for W+jets
        'higgs': "#253494", # dark blue
        'fake': "#a1d99b" # green
    }

    # Use a color map to dynamically generate colors based on the number of MC categories
    num_mc_categories = len(hists_mc_cat)
    cmap = cm.get_cmap('tab10', num_mc_categories)

    # Plot stacked MC histograms
    mc_total = None  # This will store the total MC histogram
    mc_uncertainty_squared = None  # To accumulate uncertainties squared
    mc_cumulative = None  # This will keep track of the cumulative height of the stack
    legend_entries = set()  # Track legend entries

    mc_labels = list(hists_mc_cat.keys())
    # print(f"MC labels: {mc_labels}")

    for idx, dataset in enumerate(mc_labels):
        # print(f"Processing dataset: {dataset}")
        mc_hist = hists_mc_cat[dataset]

        mc_values = np.maximum(mc_hist.values(), 0)  # Ensure MC values are non-negative
        mc_variances = mc_hist.variances()

        # Initialize cumulative and uncertainty sum
        if mc_cumulative is None:
            mc_cumulative = mc_values.copy()
            mc_uncertainty_squared = mc_variances.copy()
        else:
            mc_cumulative += mc_values
            mc_uncertainty_squared += mc_variances  # Sum uncertainties in quadrature

        # Assign color from the color map
        color = cmap(idx / num_mc_categories)
        if dataset.startswith("st") or dataset.startswith("tt"):
            label = "t/tbar"
            color = colors['tt']
        elif dataset in ["ww", "wz", "zz"]:
            label = "Diboson"
            color = colors["diboson"]
        elif dataset.startswith('dy'):
            label = "Drell-Yan"
            color = colors['dy']
        elif dataset.startswith('wj'):
            label = "W+jets"
            color = colors['wj']
        else:
            label = dataset  # Fallback for any other dataset

        # Plot the stacked MC histogram
        ax.bar(mc_hist.axes[0].centers, mc_values, width=np.diff(mc_hist.axes[0].edges),
               bottom=mc_cumulative - mc_values, label=label if label not in legend_entries else "",
               align='center', color=color, alpha=0.9)

        legend_entries.add(label)

        # If data-MC subtraction histograms are available, plot them
        if hists_mc_data_sub is not None and dataset == list(hists_mc_cat.keys())[-1]:
            if category == 4307 or category == 4305: # D or D0
                color = colors['fake']
                label = "Fake jets"
            else:
                color = '#f1b6da'
                label = 'Data - MC'
            data_mc_values = hists_mc_data_sub.values()
            mc_cumulative += data_mc_values
            ax.bar(hists_mc_data_sub.axes[0].centers, data_mc_values,
                   width=np.diff(hists_mc_data_sub.axes[0].edges),
                   bottom=mc_cumulative - data_mc_values, color=color, alpha=0.9, label=label, align='center')

    # Plot MC uncertainty as a shaded band
    mc_uncertainty = np.sqrt(mc_uncertainty_squared)  # Take square root to get standard deviation
    ax.fill_between(mc_hist.axes[0].centers,
                    mc_cumulative - mc_uncertainty,
                    mc_cumulative + mc_uncertainty,
                    step='mid', hatch='///', label='MC Uncertainty',facecolor='none',linewidth=0)

    # Plot the data histogram with error bars
    data_hist = hists_data_cat_sum
    print("data_hist ", data_hist)
    data_values = np.maximum(data_hist.values(), 0)  # Ensure data values are non-negative
    data_variances = data_hist.variances()

    ax.errorbar(data_hist.axes[0].centers, data_values, yerr=np.sqrt(data_variances),
                fmt='o', color='black', label='Data', capsize=3, markersize=5)

    # Add labels, title, and legend
    ax.set_xlabel("Tau $p_T$ / GeV")
    ax.set_ylabel("Events")
    ax.set_title(f"Stacked MC vs Data - Category {cat}")
    ax.legend()

    # Show the plot and save it
    plt.tight_layout()
    if hists_mc_data_sub is not None:
        plt.savefig(f"script_FF/plots/stacked_plot_cat_with_sub_{cat}.png", dpi=140)
    else:
        plt.savefig(f"script_FF/plots/stacked_plot_cat_{cat}.png", dpi=140)


def calculate_data_MC_substraction(hists_data_cat_, hists_mc_cat_, cat,cat_title=""):
    '''
    Calculate the data-MC subtraction for a specific category and plot the resulting histogram
    Returns: histogram of the data-MC subtraction
    '''
    # Copy
    hists_data_cat = copy.deepcopy(hists_data_cat_)
    hists_mc_cat = copy.deepcopy(hists_mc_cat_)

    hists_data_cat_sum = None
    hists_mc_cat_sum = None

    # print(hists_data_cat)

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
            
    # print(f">>>>>>>>>>>>>>>>>>>>> Data histograms sum: \n{hists_data_cat_sum}")
    # print(f">>>>>>>>>>>>>>>>>>>>> MC histograms sum: \n{hists_mc_cat_sum}")

    # Guarantee non-negative values
    hists_data_cat_sum.view().value = np.maximum(hists_data_cat_sum.values(), 0)
    hists_mc_cat_sum.view().value = np.maximum(hists_mc_cat_sum.values(), 0)

    # print(f">>>>>>>>>>>>>>>>>>>>> Data histograms sum after non-0: \n{hists_data_cat_sum}")
    # print(f">>>>>>>>>>>>>>>>>>>>> MC histograms sum after non-0: \n{hists_mc_cat_sum}")


    hists_data_cat_sum_num = hist_to_num(hists_data_cat_sum)
    hists_mc_cat_sum_num = hist_to_num(hists_mc_cat_sum)
    

    # Subtract MC from data
    hists_data_mc_sub_num = hists_data_cat_sum_num - hists_mc_cat_sum_num


    # combine uncertainties and store values in bare arrays
    hists_data_mc_sub_values = hists_data_mc_sub_num()
    hists_data_mc_sub_variances = hists_data_mc_sub_num(sn.UP, sn.ALL, unc=True)**2

    # print(f">>>>>>>>>>>>>>>>>>>>> Data - MC subtraction values: \n{hists_data_mc_sub_values}")
    # print(f">>>>>>>>>>>>>>>>>>>>> Data - MC subtraction variances: \n{hists_data_mc_sub_variances}")    

    # # Create hist from the subtraction
    # hists_data_mc_sub_hist = hist.hist(hists_data_mc_sub_num.axes[0], hists_data_mc_sub_values, hists_data_mc_sub_variances)



    # Create histo copy of the subtraction
    hists_data_mc_sub = hists_mc_cat_sum.copy().reset()

    # print(f">>>>>>>>>>>>>>>>>>>>> Data - MC subtraction hist: \n{hists_data_mc_sub}")

    hists_data_mc_sub.view().value[...] = hists_data_mc_sub_values
    hists_data_mc_sub.view().variance[...] = hists_data_mc_sub_variances

    # Ensure non-negative values
    hists_data_mc_sub.view().value = np.maximum(hists_data_mc_sub.values(), 0)


    # print(f">>>>>>>>>>>>>>>>>>>>> Data - MC subtraction hist: \n{hists_data_mc_sub}")
    

    fig, ax = plt.subplots()
    ax.bar(hists_data_mc_sub.axes[0].centers, hists_data_mc_sub.values(), width=np.diff(hists_data_mc_sub.axes[0].edges), color='#f1b6da')
    ax.set_ylabel("Events")
    ax.set_title(f"Data - MC Subtraction - Category {cat_title}")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"script_FF/plots/data_mc_subtraction_cat_{cat_title}.png", dpi=140)
    # plt.show()
    return hists_data_mc_sub


def calculate_data_MC_subtraction_2D(hists_data_cat, hists_mc_cat, cat, var1, var2):
    '''
    Calculate the data-MC subtraction for a specific 2D category and plot the resulting histogram.
    Returns: 2D histogram of the data-MC subtraction.
    '''

    # Initialize sums for data and MC
    hists_data_cat_sum = None
    hists_mc_cat_sum = None

    # Sum data histograms
    for dataset in hists_data_cat:
        hist = hists_data_cat[dataset]
        hists_data_cat_sum = hist if hists_data_cat_sum is None else hists_data_cat_sum + hist

    # Sum MC histograms
    for dataset in hists_mc_cat:
        hist = hists_mc_cat[dataset]
        hists_mc_cat_sum = hist if hists_mc_cat_sum is None else hists_mc_cat_sum + hist

    # print(f">>>>>>>>>>>>>>>>>>>>> Data histograms sum: \n{hists_data_cat_sum}")
    # print(f">>>>>>>>>>>>>>>>>>>>> MC histograms sum: \n{hists_mc_cat_sum}")

    # Ensure non-negative values in histograms
    hists_data_cat_sum.view().value = np.maximum(hists_data_cat_sum.values(), 0)
    hists_mc_cat_sum.view().value = np.maximum(hists_mc_cat_sum.values(), 0)

    # Convert histograms to sn.Number objects for arithmetic operations
    hists_data_cat_sum_num = hist_to_num(hists_data_cat_sum)
    hists_mc_cat_sum_num = hist_to_num(hists_mc_cat_sum)

    # Perform Data - MC subtraction
    hists_data_mc_sub_num = hists_data_cat_sum_num - hists_mc_cat_sum_num


    # combine uncertainties and store values in bare arrays
    hists_data_mc_sub_values = hists_data_mc_sub_num()
    hists_data_mc_sub_variances = hists_data_mc_sub_num(sn.UP, sn.ALL, unc=True)**2


    # Create histo copy of the subtraction
    hists_data_mc_sub = hists_mc_cat_sum.copy().reset()

    # print(f">>>>>>>>>>>>>>>>>>>>> Data - MC subtraction hist: \n{hists_data_mc_sub}")

    hists_data_mc_sub.view().value[...] = hists_data_mc_sub_values
    hists_data_mc_sub.view().variance[...] = hists_data_mc_sub_variances

    # Ensure non-negative values
    hists_data_mc_sub.view().value = np.maximum(hists_data_mc_sub.values(), 0)

    # Prepare category title mapping
    category_titles = {
        410101: "A_pi_1",
        410103: "B_pi_1",
    }
    cat_title = category_titles.get(cat, str(cat))

    # Plot using hist.plot2d_full for 2D plotting
    fig, ax = plt.subplots(figsize=(8, 6))
    hists_data_mc_sub.plot2d(ax=ax)
    # ax.set_xlabel(f"{var1}")
    # ax.set_ylabel(f"{var2}")
    ax.set_title(f"Data - MC Subtraction - Category {cat_title}")

    plt.tight_layout()
    plt.savefig(f"script_FF/plots/data_mc_subtraction_2D_cat_{cat_title}.png", dpi=140)


    return hists_data_mc_sub


def calculate_ratio_hist(hist_iso, hist_antiiso,title=""):
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

    # Calculate bin widths
    bin_edges = ratio_hist.axes[0].edges  # Get bin edges
    bin_widths = bin_edges[1:] - bin_edges[:-1]  # Calculate widths of bins

    # Plot the ratio
    fig, ax = plt.subplots()
    ax.errorbar(ratio_hist.axes[0].centers, ratio_hist.values(), yerr=np.sqrt(ratio_hist.variances()),xerr=bin_widths / 2, fmt='o', label='Ratio', color='black')
    ax.set_ylabel("Ratio")
    ax.set_xlabel("Tau $p_T$ / GeV")
    # set ylim
    ax.set_ylim(0, 0.5)
    ax.set_title(f"Ratio of Iso / AntiIso {title}")
    ax.legend()
    plt.savefig(f"script_FF/plots/FakeFactor_{title}.png", dpi=140)

    # Extract data from the histogram
    bin_centers = ratio_hist.axes[0].centers
    values = ratio_hist.values()
    variances = ratio_hist.variances()
    yerr = np.sqrt(variances)
    xerr = bin_widths / 2

    # Replace any Infinity values in `values` and `yerr` with 0
    values = np.where(np.isinf(values), 0, values)
    yerr = np.where(np.isinf(yerr), 0, yerr)

    # Extract data for JSON
    # Convert data to a dictionary format
    data_dict = {
        "title": title,
        "bin_centers": bin_centers.tolist(),
        "values": values.tolist(),
        "yerr": yerr.tolist(),
        "xerr": xerr.tolist()
    }

    # Save the data to a JSON file
    with open(f"script_FF/data/FakeFactor_{title}.json", 'w') as json_file:
        json.dump(data_dict, json_file, indent=4)

    return ratio_values, ratio_variances, ratio_hist

def plot_hist_cat_with_ratio(hists_data_cat_, hists_mc_cat_, cat, hists_mc_data_sub_,fake=False, var = 'hcand_1_pt', cat_title=""):
    '''
    Plot histograms for a specific category, including MC uncertainties and a ratio plot
    '''

    hep.style.use(hep.style.CMS)

    hists_data_cat = copy.deepcopy(hists_data_cat_)
    hists_mc_cat = copy.deepcopy(hists_mc_cat_)

    hists_mc_data_sub = hists_mc_data_sub_.copy() if hists_mc_data_sub_ is not None else None

    # Initialize the sum variable for the data histogram
    hists_data_cat_sum = None

    # Loop over the datasets to sum data histograms
    for dataset in hists_data_cat:
        hist = hists_data_cat[dataset]

        if hists_data_cat_sum is None:
            hists_data_cat_sum = hist  # Start with the first histogram
        else:
            hists_data_cat_sum += hist  # Accumulate the histograms

    # Prepare the figure with ratio subplot
    fig, (ax_main, ax_ratio) = plt.subplots(nrows=2, ncols=1, gridspec_kw={'height_ratios': [3, 1]}, sharex=True)


    # Define colors for specific datasets
    colors = {
        'tt': "#9e9ac8",  # Violet
        'dy': "#feb24c",  # Orange
        'diboson': "#a96b59",   # Brown for Diboson (WW)
        'wj': "#d73027", # Red for W+jets
        'higgs': "#253494", # dark blue 
        'fake': "#addd8e" # green
    }

    # Use a color map to dynamically generate colors based on the number of MC categories
    num_mc_categories = len(hists_mc_cat)
    cmap = cm.get_cmap('tab10', num_mc_categories)

    # Plot stacked MC histograms
    mc_total = None  # This will store the total MC histogram
    mc_uncertainty_squared = None  # To accumulate uncertainties squared
    mc_cumulative = None  # This will keep track of the cumulative height of the stack
    legend_entries = set()  # Track legend entries

    mc_labels = list(hists_mc_cat.keys())
    # print(f"MC labels: {mc_labels}")

    for idx, dataset in enumerate(mc_labels):
        # print(f"Processing dataset: {dataset}")
        mc_hist = hists_mc_cat[dataset]

        mc_values = np.maximum(mc_hist.values(), 0)  # Ensure MC values are non-negative
        mc_variances = mc_hist.variances()

        # Initialize cumulative and uncertainty sum
        if mc_cumulative is None:
            mc_cumulative = mc_values.copy()
            mc_uncertainty_squared = mc_variances.copy()
        else:
            mc_cumulative += mc_values
            mc_uncertainty_squared += mc_variances  # Sum uncertainties in quadrature

        # Assign color from the color map
        color = cmap(idx / num_mc_categories)
        if dataset.startswith("st") or dataset.startswith("tt"):
            label = "t/tbar"
            color = colors['tt']
        elif dataset in ["ww", "wz", "zz"]:
            label = "Diboson"
            color = colors["diboson"]
        elif dataset.startswith('dy'):
            label = "Drell-Yan"
            color = colors['dy']
        elif dataset.startswith('wj'):
            label = "W+jets"
            color = colors['wj']
        elif dataset.startswith('h_gg'):
            label = "Higgs"
            color = colors['higgs']
        else:
            label = dataset  # Fallback for any other dataset

        # Plot the stacked MC histogram
        ax_main.bar(mc_hist.axes[0].centers, mc_values, width=np.diff(mc_hist.axes[0].edges),
                    bottom=mc_cumulative - mc_values, label=label if label not in legend_entries else "",
                    align='center', color=color, alpha=0.9)

        legend_entries.add(label)

      # If data-MC subtraction histograms are available, plot them
    if hists_mc_data_sub is not None and dataset == list(hists_mc_cat.keys())[-1]:
        if fake: # D or D0
            color = colors['fake']
            label = "Fake jets"
        else:
            color = '#f1b6da'
            label = 'Data - MC'
        data_mc_values = hists_mc_data_sub.values()
        mc_cumulative += data_mc_values
        ax_main.bar(hists_mc_data_sub.axes[0].centers, data_mc_values,
                width=np.diff(hists_mc_data_sub.axes[0].edges),
                bottom=mc_cumulative - data_mc_values, color=color, alpha=0.9, label=label, align='center')



    # Plot the data histogram with error bars
    data_hist = hists_data_cat_sum
    data_values = np.maximum(data_hist.values(), 0)  # Ensure data values are non-negative
    data_variances = data_hist.variances()

    ax_main.errorbar(data_hist.axes[0].centers, data_values, yerr=np.sqrt(data_variances),
                     fmt='o', color='black', label='Data', capsize=3, markersize=5)

    # Plot the MC uncertainty as a shaded band
    mc_uncertainty = np.sqrt(mc_uncertainty_squared)  # Take square root to get standard deviation
    ax_main.fill_between(mc_hist.axes[0].centers,
                         mc_cumulative - mc_uncertainty,
                         mc_cumulative + mc_uncertainty,
                         step='mid', hatch='///', label='MC Uncertainty', facecolor='none', linewidth=0)

    # Add labels, title, and legend to the main plot
    ax_main.set_ylabel("Events")
    ax_main.set_title(f"Stacked MC vs Data - Category {cat_title}")
    ax_main.legend()

    # Plot the ratio (data/MC) on the ratio subplot
    ratio_values = data_values / mc_cumulative
    ratio_errors = np.sqrt(data_variances) / mc_cumulative

    ax_ratio.errorbar(data_hist.axes[0].centers, ratio_values, yerr=ratio_errors,
                      fmt='o', color='black', label='Data/MC', capsize=3, markersize=5)
    ax_ratio.axhline(1, color='red', linestyle='--')  # Reference line at ratio = 1
    ax_ratio.set_ylabel("Data/MC")
    ax_ratio.set_ylim(0, 2)
    ax_ratio.set_xlabel(var)

    # Show the plot and save it
    plt.tight_layout()

    # Save in Json 
    # Extract data centers, values, and errors
    bin_centers = data_hist.axes[0].centers  # Bin centers for the ratio plot
    # Replace NaN or infinite values with zero to avoid issues
    ratio_values = np.nan_to_num(ratio_values, nan=0.0, posinf=0.0, neginf=0.0)
    ratio_errors = np.nan_to_num(ratio_errors, nan=0.0, posinf=0.0, neginf=0.0)

    # Create dictionary for JSON output
    data_dict = {
    "title": cat_title,  # Category title
    "bin_centers": bin_centers.tolist(),    # Bin centers for each histogram bin
    "ratio_values": ratio_values.tolist(),   # Data/MC ratio values
    "ratio_errors": ratio_errors.tolist()    # Errors in the Data/MC ratio
    }


    # Save data_dict to JSON
    output_path = f"script_FF/data/control_plots/{var}_DataMC_Ratio_{cat_title}.json"
    with open(output_path, "w") as json_file:
        json.dump(data_dict, json_file, indent=4)


    if hists_mc_data_sub is not None:
        if fake:
            plt.savefig(f"script_FF/plots/{var}_stacked_plot_cat_with_sub_{cat_title}_ratio_fake.png", dpi=140)
        else:
            plt.savefig(f"script_FF/plots/{var}_stacked_plot_cat_with_sub_{cat_title}_ratio.png", dpi=140)
    else:
        plt.savefig(f"script_FF/plots/{var}_stacked_plot_cat_{cat_title}_ratio.png", dpi=140)

def process_histograms(base_path, task, datasets_data, datasets_mc, hist_path, category, cats_title):
    hists_data, hists_mc = get_all_hist(base_path, task, datasets_data, datasets_mc, hist_path)
    hists_mc_data_sub = {}

    # Define the new binning    
    
    for cat, cat_title in zip(category, cats_title):
        hists_data_cat, hists_mc_cat = get_hist_by_category_var(hists_data, hists_mc, cat,var = 'hcand_1_pt')

        plot_hist_cat_with_ratio(hists_data_cat, hists_mc_cat, cat, None,fake=False, cat_title=cat_title)
        hists_mc_data_sub[cat] = calculate_data_MC_substraction(hists_data_cat, hists_mc_cat, cat, cat_title=cat_title)
        plot_hist_cat_with_ratio(hists_data_cat, hists_mc_cat, cat, hists_mc_data_sub[cat],fake=False, cat_title=cat_title)

    
    return hists_mc_data_sub


def apply_fake_factor(hists_mc_data_sub_SR, fakefator_hist):

    hists_mc_data_sub_SR_num = hist_to_num(hists_mc_data_sub_SR)
    fakefator_hist_num = hist_to_num(fakefator_hist)
    
    hists_mc_data_sub_SR_num = hists_mc_data_sub_SR_num * fakefator_hist_num  # SR = AR * FF

    hists_mc_data_sub_SR_num_values = hists_mc_data_sub_SR_num()
    hists_mc_data_sub_SR_num_variances = hists_mc_data_sub_SR_num(sn.UP, sn.ALL, unc=True)**2
    
    return hists_mc_data_sub_SR_num_values, hists_mc_data_sub_SR_num_variances



def apply_fake_factor_var(hist_mc_data_sub_var_2D, fakefactor_hist, var_axis):
    '''
    Apply a pt-dependent fake factor to the data-MC subtraction histogram for a specific variable axis 
    (e.g., MET_pt, Njets) based on the corresponding pt distribution.
    
    Parameters:
    hist_mc_data_sub_pt (2D Histogram): The pt-var-based data-MC subtraction histogram (used for the fake factor calculation).

    fakefactor_hist (Histogram): The fake factor histogram, typically derived from region A/B or C/D, based on pt.
    var_axis (str): The name of the variable axis (e.g., 'MET_pt') to which the fake factor is being applied.
    
    Returns:
    hist_fake (1D Histogram): The data-MC subtraction histogram of the var with the pt-dependent fake factor applied.

    bin_0 (var) = bin_0_0 (var,pt) * bin_0 (FF) + ... + bin_0_n (var,pt) * bin_n (FF)
    '''
    
      # Initialize a 1D histogram along the `var_axis` with the same binning as the original
    hist_fake = hist_mc_data_sub_var_2D.project(var_axis).copy()
    hist_fake.reset()  # Clear values for accumulation

    # Calculate the pt integral for normalization
    pt_hist = hist_mc_data_sub_var_2D.project('hcand_1_pt').copy()
    pt_integral = np.sum(pt_hist.values())
    # pt_integral = np.sum(pt_hist.view().value, axis=1)
    # print(f"pt_integral: {pt_integral}")

    # Initialize arrays to hold values and variances for final assignment
    var_values = np.zeros(hist_mc_data_sub_var_2D.axes[0].size)
    var_variances = np.zeros(hist_mc_data_sub_var_2D.axes[0].size)

    # Iterate over each bin in the `var_axis` and `pt` axis to apply the fake factor
    for var_bin in range(hist_mc_data_sub_var_2D.axes[0].size):
        total_fake_content = 0
        total_fake_variance = 0

        for pt_bin in range(hist_mc_data_sub_var_2D.axes[1].size):
            # Get the bin content and variance from the 2D data-MC histogram
            var_pt_content = hist_mc_data_sub_var_2D.view().value[var_bin, pt_bin]
            var_pt_variance = hist_mc_data_sub_var_2D.view().variance[var_bin, pt_bin]

            # avoid nan values

            # Get the fake factor and its uncertainty for the current pt bin
            fake_factor = fakefactor_hist.view().value[pt_bin]
            fake_factor_variance = fakefactor_hist.view().variance[pt_bin]

            # Set NaN values to 0
            var_pt_content = np.nan_to_num(var_pt_content, nan=0.0)
            var_pt_variance = np.nan_to_num(var_pt_variance, nan=0.0)
            fake_factor = np.nan_to_num(fake_factor, nan=0.0)
            fake_factor_variance = np.nan_to_num(fake_factor_variance, nan=0.0)


            # Calculate the contribution to the current `var` bin by applying the fake factor
            fake_content = var_pt_content * fake_factor
            fake_variance = (var_pt_content ** 2 * fake_factor_variance) + (var_pt_variance * fake_factor ** 2)

            # Accumulate total content and variance for the current var bin
            total_fake_content += fake_content
            total_fake_variance += fake_variance

        # Store calculated content and variance for the `var` bin
        var_values[var_bin] = total_fake_content
        var_variances[var_bin] = total_fake_variance

    # Assign all calculated values and variances at once
    hist_fake.view().value[...] = var_values
    hist_fake.view().variance[...] = var_variances
    # print(f"hist_fake: {hist_fake}")


    # plot the hist_fake

    fig, ax = plt.subplots()
    ax.errorbar(hist_fake.axes[0].centers, hist_fake.values(), yerr=np.sqrt(hist_fake.variances()), fmt='o', label='Fake Factor Applied', color='black')
    ax.set_ylabel("Events")
    ax.set_xlabel(f"{var_axis}")
    ax.set_title(f"Fake Factor Applied to {var_axis}")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"script_FF/plots/FakeFactor_{var_axis}_applied.png", dpi=140)


    return hist_fake




def plot_fake_factors(fakefator_hist, fakefator_values, fakefator_variance, fakefator_hist_0, fakefator_values_0, fakefator_variance_0):
    fig, ax = plt.subplots()
    ax.errorbar(fakefator_hist.axes[0].centers, fakefator_values, yerr=np.sqrt(fakefator_variance), fmt='o', label='Fake Factor A/B', color='#d73027')
    ax.errorbar(fakefator_hist_0.axes[0].centers, fakefator_values_0, yerr=np.sqrt(fakefator_variance_0), fmt='o', label='Fake Factor A0/BO', color='#c994c7')
    ax.set_ylabel("Fake Factor")
    ax.set_xlabel("Tau $p_T$ / GeV")
    ax.set_title("Fake Factor")
    ax.legend()
    plt.savefig("script_FF/plots/FakeFactor_comparaison.png", dpi=140)




def main(config_path):
    # Load the JSON configuration file
    with open(config_path, 'r') as f:
        config = json.load(f)

    # Extract configuration parameters
    eos_path = config['paths']['eos_path']
    task = config['paths']['task']
    hist_path_base = config['paths']['hist_path']
    hist_path_FF = hist_path_base + "hist__hcand_1_pt.pickle"
    dataset_data = config['datasets']['data']
    dataset_mc = config['datasets']['mc']
    categories = config['categories']
    variables = config['variables']

    # Load histograms for data and MC from the specified paths
    hists_data, hists_mc = get_all_hist(eos_path, task, dataset_data, dataset_mc, hist_path_FF)



    # Loop through each DM category
    for dm_name, dm_categories in categories.items():
        print(f"Processing DM category {dm_name}")
        # ABCD or ABCD0
        for grouping_name, category_dict in dm_categories.items(): 
            print(f"Processing grouping {grouping_name}")

            if grouping_name == "ABCD":
                catA, catB, catC, catD = category_dict['A'], category_dict['B'], category_dict['C'], category_dict['D']
                cats_titleABC = ["A_"+dm_name, "B_"+dm_name, "C_"+dm_name]
                cats_titleD = "D_"+dm_name
            elif grouping_name == "A0B0C0D0":
                catA, catB, catC, catD = category_dict['A0'], category_dict['B0'], category_dict['C0'], category_dict['D0']
                cats_titleABC = ["A0_"+dm_name, "B0_"+dm_name, "C0_"+dm_name]
                cats_titleD = "D0_"+dm_name
            else:
                raise ValueError(f"Invalid grouping name: {grouping_name}")


            category = [catA, catB, catC, catD]
            categoryABC = [catA, catB, catC]
    
            hists_mc_data_sub_ABC = process_histograms(eos_path, task, dataset_data, dataset_mc, hist_path_FF, categoryABC, cats_titleABC)

            # Derive Fake factors
            title = grouping_name + "_"+ dm_name
            print(f"Processing Fake Factor for {title}")
            fakefactor_values_ABCD, fakefactor_variance_ABCD, fakefactor_hist_ABCD = calculate_ratio_hist(hists_mc_data_sub_ABC[catA], hists_mc_data_sub_ABC[catB], title)   
        
            hists_mc_data_sub_C = hists_mc_data_sub_ABC[catC]

            ## SR : region D

            ## Calculate the fake factor
            ## Create hist from the subtraction which will be fake process in the SR
            hists_mc_data_sub_D = hists_mc_data_sub_C.copy().reset()

            # Get histograms for the category and plot them
            hists_data_cat_D, hists_mc_cat_D = get_hist_by_category_var(hists_data, hists_mc, catD ,var = 'hcand_1_pt')
            plot_hist_cat_with_ratio(hists_data_cat_D, hists_mc_cat_D, catD,None,fake=False, cat_title=cats_titleD)

            # Apply fake factor to the D region from the C region
            hists_mc_data_sub_D_num_values, hists_mc_data_sub_D_num_variances = apply_fake_factor(hists_mc_data_sub_C, fakefactor_hist_ABCD)

            # Create hist from the subtraction
            hists_mc_data_sub_D.view().value[...] = hists_mc_data_sub_D_num_values
            hists_mc_data_sub_D.view().variance[...] = hists_mc_data_sub_D_num_variances

            # Plot
            plot_hist_cat_with_ratio(hists_data_cat_D, hists_mc_cat_D ,catD , hists_mc_data_sub_D,fake=True,cat_title=cats_titleD)


            # Loop through each variable pair
            for variable_pair in variables:
                var1, var2 = variable_pair['var1'], variable_pair['var2']
                print(f"Processing variable pair: {var1}, {var2}")

                # Check if it's closure correction (ABCD) or extrapolation correction (A0B0C0D0)
                # cat1 = FF X cat2
                if grouping_name == "ABCD":
                    cat1, cat2 = category_dict['A'], category_dict['B']
                    cat1_title = "A_"+dm_name
                elif grouping_name == "A0B0C0D0":
                    cat1, cat2 = category_dict['D0'], category_dict['C0']
                    cat1_title = "D0_"+dm_name
                else:
                    raise ValueError(f"Invalid grouping name: {grouping_name}")


                # Define histogram path for 2D variables
                hist_path_2D = f"{hist_path_base}hist__{var1}-{var2}.pickle"
                hists_data_2D, hists_mc_2D = get_all_hist(eos_path, task, dataset_data, dataset_mc, hist_path_2D)


                # Get 2d histogram for the variable and category
                hists_data_cat_1, hists_mc_cat_1 = get_2D_hist_by_category_var(hists_data_2D, hists_mc_2D, cat1,var1, var2)
                hists_data_cat_2, hists_mc_cat_2 = get_2D_hist_by_category_var(hists_data_2D, hists_mc_2D, cat2,var1, var2)


                # ## Get hist sub of each category and
                hists_mc_data_sub_1 = calculate_data_MC_subtraction_2D(hists_data_cat_1, hists_mc_cat_1, cat1, var1, var2)
                hists_mc_data_sub_2 = calculate_data_MC_subtraction_2D(hists_data_cat_2, hists_mc_cat_2, cat2, var1, var2)

                # Get var hist with fake factor applied
                hist_fake = apply_fake_factor_var(hists_mc_data_sub_2, fakefactor_hist_ABCD, var1)

                # Project the 2D histogram to 1D for plotting for var1 using hist_fake
                # Initialize a new dictionary to hold the projected 1D histograms
                hists_data_cat_1_1D = {}
                hists_mc_cat_1_1D = {}

                # Project each histogram in hists_data_cat_D to 1D for var1
                for key, hist in hists_data_cat_1.items():
                    if hist:  # Check if the histogram is not None
                        # Project the 2D histogram to 1D
                        hists_data_cat_1_1D[key] = hist.project(var1)

                # Project each histogram in hists_mc_cat_D to 1D for var1
                for key, hist in hists_mc_cat_1.items():
                    if hist:  # Check if the histogram is not None
                        # Project the 2D histogram to 1D
                        hists_mc_cat_1_1D[key] = hist.project(var1)


                # Plot the 1D histograms for the variable
                plot_hist_cat_with_ratio(hists_data_cat_1_1D, hists_mc_cat_1_1D, cat1, hist_fake,fake=True, var=var1, cat_title=cat1_title)


if __name__ == "__main__":
    main("script_FF/fake_factors.json")
