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
            print(f"Warning: Category {category} not found in the histogram data for dataset {dataset}")
            continue
        category_index = category_axis.index(category)
        hists_data_cat[dataset] = hists_data[dataset][category_index, :, :, :].project('tau_1_pt')
    
    # Loop over the datasets to get the histograms for the category
    for dataset in hists_mc.keys():
        category_axis = hists_mc[dataset].axes['category']
        if category not in category_axis:
            print(f"Warning: Category {category} not found in the histogram data for dataset {dataset}")
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
            print(f"Warning: Category {category} not found in the histogram data for dataset {dataset}")
            continue
        category_index = category_axis.index(category)
        hists_data_np[dataset] = hists_data[dataset][category_index, :, :, :].project('tau_1_pt').values()
        hists_data_variance_np[dataset] = hists_data[dataset][category_index, :, :, :].project('tau_1_pt').variances()**0.5
    
    # Loop over the datasets to convert the histograms to numpy arrays
    for dataset in hists_mc.keys():
        category_axis = hists_mc[dataset].axes['category']
        if category not in category_axis:
            print(f"Warning: Category {category} not found in the histogram data for dataset {dataset}")
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


# def plot_hist_cat(hists_data_cat_, hists_mc_cat_, cat, hists_mc_data_sub_):
#     '''
#     Plot histograms for a specific category
#     '''

#     hep.style.use(hep.style.CMS)

#     hists_data_cat = copy.deepcopy(hists_data_cat_)
#     hists_mc_cat = copy.deepcopy(hists_mc_cat_)

#     hists_mc_data_sub = hists_mc_data_sub_.copy() if hists_mc_data_sub_ is not None else None
    
#     # print(f">>>>>>>>>>>>>>>>>>>>> Data histograms: \n{hists_data_cat}")
#     # print(f">>>>>>>>>>>>>>>>>>>>> MC histograms: \n{hists_mc_cat}")
#     # print(f">>>>>>>>>>>>>>>>>>>>> Data - MC subtraction histograms: \n{hists_mc_data_sub}")

#     # Initialize the sum variable as None
#     hists_data_cat_sum = None

#     # Loop over the datasets to sum histograms
#     for dataset in hists_data_cat:
#         # print(dataset)
#         hist = hists_data_cat[dataset]

#         # Initialize the sum histogram if it's None
#         if hists_data_cat_sum is None:
#             hists_data_cat_sum = hist  # Start with the first histogram
#         else:
#             hists_data_cat_sum += hist  # Accumulate the histograms

#     # Prepare the figure and axis for plotting
#     fig, ax = plt.subplots()

#      # Define colors for specific datasets
#     colors = {
#         'tt': "#9e9ac8",  # Violet
#         'dy': "#feb24c",  # Orange
#         'diboson': "#a96b59",   # Brown for Diboson (WW)
#         'wj': "#d73027", # Red for W+jets
#     }

#     # Use a color map to dynamically generate colors based on the number of MC categories
#     num_mc_categories = len(hists_mc_cat)
#     cmap = cm.get_cmap('tab10', num_mc_categories)  # Use 'tab10', 'viridis', or any other colormap

#     # Plot stacked MC histograms
#     mc_total = None  # This will store the total MC histogram
#     mc_labels = list(hists_mc_cat.keys())
#     print(f"MC labels: {mc_labels}")

#     # Prepare a variable to keep track of the cumulative height of the stack
#     mc_cumulative = None  
#     # Keep track of legend entries
#     legend_entries = set()

#     for idx, dataset in enumerate(mc_labels):
#         print(f"Processing dataset: {dataset}")
#         mc_hist = hists_mc_cat[dataset]

#         # Ensure no negative values and only work with non-negative values
#         mc_values = np.maximum(mc_hist.values(), 0)  # Ensure MC values are non-negative

#         # print(f"MC values: {mc_values}")

#         # If mc_cumulative is None, start with the first histogram
#         if mc_cumulative is None:
#             mc_cumulative = mc_values.copy()  # Copy the first histogram values
#         else:
#             mc_cumulative += mc_values  # Add to the cumulative histogram

    
#         # Create a dictionary for legend labels
#         legend_labels = {}

#         # Assign color from the color map
#         color = cmap(idx / num_mc_categories)

#         # Determine the legend label based on dataset name
#         if dataset.startswith("st") or dataset.startswith("tt"):
#             label = "t/tbar"
#             color = colors['tt']
#         elif dataset in ["ww", "wz", "zz"]:
#             color = colors["diboson"]
#             label = "Diboson"
#         elif dataset.startswith('dy'):
#             label = "Drell-Yan"
#             color = colors['dy']  # Use the color for Drell-Yan
#         elif dataset.startswith('wj'):
#             label = "W+jets"
#             color = colors['wj']
#         else:
#             label = dataset  # Fallback for any other dataset
        
#         # Add to legend labels dictionary
#         legend_labels[label] = color

#         # Plot the stacked MC histogram for each dataset
#         ax.bar(mc_hist.axes[0].centers, mc_values, width=np.diff(mc_hist.axes[0].edges),
#                bottom=mc_cumulative - mc_values,  # Adjust the bottom for stacking
#                label=label if label not in legend_entries else "", align='center', color=color, alpha=0.9)

#         # Update legend entries
#         if label not in legend_entries:
#             legend_entries.add(label)  # Mark this label as added

        
#         if hists_mc_data_sub is not None and dataset == list(hists_mc_cat.keys())[-1]:
#             # Plot the data-MC subtraction histogram as part of the stack
#             color = '#f1b6da'  # Use a custom color for the data-MC subtraction (pink)
#             data_mc_values = hists_mc_data_sub.values()
#             mc_cumulative += data_mc_values
#             # Plot it on top of the cumulative stack
#             ax.bar(hists_mc_data_sub.axes[0].centers, data_mc_values, width=np.diff(hists_mc_data_sub.axes[0].edges),bottom=mc_cumulative - data_mc_values, color=color, alpha=0.9, label='Data - MC', align='center') 



#     # Plot the data histogram with error bars (as points with errors)
#     data_hist = hists_data_cat_sum
#     data_values = data_hist.values()
#     data_variances = data_hist.variances()

#     # Ensure no negative data values
#     data_values = np.maximum(data_values, 0)  # Ensure data values are non-negative

#     # Plot data as points with error bars
#     ax.errorbar(data_hist.axes[0].centers, data_values, yerr=np.sqrt(data_variances), fmt='o', 
#                 color='black', label='Data', capsize=3, markersize=5)

#     # Add labels, title, and legend
#     ax.set_xlabel("Tau $p_T$ / GeV")
#     ax.set_ylabel("Events")
#     ax.set_title(f"Stacked MC vs Data - Category {cat}")
#     ax.legend()

#     # Show the plot and save it
#     plt.tight_layout()
#     if hists_mc_data_sub is not None:
#         plt.savefig(f"script_FF/plots/stacked_plot_cat_with_sub_{cat}.png", dpi=140)
#     else:
#         plt.savefig(f"script_FF/plots/stacked_plot_cat_{cat}.png", dpi=140)
#     # plt.show()

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
    print(f"MC labels: {mc_labels}")

    for idx, dataset in enumerate(mc_labels):
        print(f"Processing dataset: {dataset}")
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
            color = '#f1b6da'  # Custom color for data-MC subtraction (pink)
            data_mc_values = hists_mc_data_sub.values()
            mc_cumulative += data_mc_values
            ax.bar(hists_mc_data_sub.axes[0].centers, data_mc_values,
                   width=np.diff(hists_mc_data_sub.axes[0].edges),
                   bottom=mc_cumulative - data_mc_values, color=color, alpha=0.9, label='Data - MC', align='center')

    # Plot MC uncertainty as a shaded band
    mc_uncertainty = np.sqrt(mc_uncertainty_squared)  # Take square root to get standard deviation
    ax.fill_between(mc_hist.axes[0].centers,
                    mc_cumulative - mc_uncertainty,
                    mc_cumulative + mc_uncertainty,
                    step='mid', hatch='///', label='MC Uncertainty',facecolor='none',linewidth=0)

    # Plot the data histogram with error bars
    data_hist = hists_data_cat_sum
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


def calculate_data_MC_substraction(hists_data_cat_, hists_mc_cat_, cat):
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
    

    # Plot the data-MC subtraction
    fig, ax = plt.subplots()
    ax.bar(hists_data_mc_sub.axes[0].centers, hists_data_mc_sub.values(), width=np.diff(hists_data_mc_sub.axes[0].edges), color='#f1b6da')
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


    return ratio_values, ratio_variances, ratio_hist


if __name__ == "__main__":
    '''
    Main function to run the script
    '''


    # Derivative region

    ## Define the paths to the pickle files
    eos_path_DR = "/eos/user/o/oponcet2/analysis/FakeFactor/analysis_httcp/"
    task_DR = "cf.MergeHistograms/run3_2022_preEE_nano_cp_tau_v12/"
    dataset_data_DR = ["data_tau_C", "data_tau_D"]
    dataset_mc_DR = ["dy_lep_m50", "h_ggf_tautau_prod_cp_even_sm", "st_tchannel_t", "st_tchannel_tbar", "st_tw_tb_dl", "st_tw_tb_fh", "st_tw_tb_sl", "st_tw_t_dl", "st_tw_t_fh", "st_tw_t_sl", "tt_dl", "tt_fh", "tt_sl", "wj_incl", "ww", "wz", "zz"]
    hist_path_DR = "nominal/calib__main/sel__main_FF/prod__main_FF/weight__all_weights/htcondor2/hist__tau_1_pt.pickle"

    hists_data_DR, hists_mc_DR =  get_all_hist(eos_path_DR, task_DR, dataset_data_DR, dataset_mc_DR, hist_path_DR)


    ## plot all histograms by category
    category = [600,500] # FFDRantiIso_tautau = 600, FFDRIso_tautau = 500, category that i want 

    hists_mc_data_sub_DR = {}
    for cat in category:
        hists_data_cat_DR, hists_mc_cat_DR = get_hist_by_category_var(hists_data_DR, hists_mc_DR, cat)
        plot_hist_cat(hists_data_cat_DR, hists_mc_cat_DR,cat,None)
        hists_mc_data_sub_DR[cat] = calculate_data_MC_substraction(hists_data_cat_DR, hists_mc_cat_DR, cat)
        # plot_hist_cat(hists_data_cat, hists_mc_cat,cat,None)
        plot_hist_cat(hists_data_cat_DR, hists_mc_cat_DR, cat, hists_mc_data_sub_DR[cat])

        
    fakefator_values, fakefator_variance, fakefator_hist = calculate_ratio_hist(hists_mc_data_sub_DR[500], hists_mc_data_sub_DR[600])


    ## AR region
    ## Define the paths to the pickle files
    eos_path_AR = "/eos/user/o/oponcet2/analysis/FakeFactor/analysis_httcp/"
    task_AR = "cf.MergeHistograms/run3_2022_preEE_nano_cp_tau_v12/"
    dataset_data_AR = ["data_tau_C", "data_tau_D"]
    dataset_mc_AR = ["dy_lep_m50", "h_ggf_tautau_prod_cp_even_sm", "st_tchannel_t", "st_tchannel_tbar", "st_tw_tb_dl", "st_tw_tb_fh", "st_tw_tb_sl", "st_tw_t_dl", "st_tw_t_fh", "st_tw_t_sl", "tt_dl", "tt_fh", "tt_sl", "wj_incl", "ww", "wz", "zz"]
    hist_path_AR = "nominal/calib__main/sel__main_FF/prod__main_FF/weight__all_weights/htcondor4/hist__tau_1_pt.pickle"

    hists_data_AR, hists_mc_AR =  get_all_hist(eos_path_AR, task_AR, dataset_data_AR, dataset_mc_AR, hist_path_AR)
    category = 700 # tautau_antiIso = 700, category that i want 

    hists_data_cat_AR, hists_mc_cat_AR = get_hist_by_category_var(hists_data_AR, hists_mc_AR, category)
    plot_hist_cat(hists_data_cat_AR, hists_mc_cat_AR,category,None)
    hists_mc_data_sub_AR = calculate_data_MC_substraction(hists_data_cat_AR, hists_mc_cat_AR, category)
    plot_hist_cat(hists_data_cat_AR, hists_mc_cat_AR, category, hists_mc_data_sub_AR)

    ## SR 
    ## Define the paths to the pickle files
    hists_data_SR, hists_mc_SR =  hists_data_DR, hists_mc_DR
    category = 301 # tautau
    hists_data_cat_SR, hists_mc_cat_SR = get_hist_by_category_var(hists_data_SR, hists_mc_SR, category)
    plot_hist_cat(hists_data_cat_SR, hists_mc_cat_SR,category,None)

    from IPython import embed; embed()

    ## Calculate the fake factor
    hists_mc_data_sub_SR = hists_mc_data_sub_AR.copy().reset()

    hists_mc_data_sub_SR_num = hist_to_num(hists_mc_data_sub_SR)
    hists_mc_data_sub_AR_num = hist_to_num(hists_mc_data_sub_AR)
    fakefator_hist_num = hist_to_num(fakefator_hist)

    # ABCD method for the fake factor
    # shape: (SHIFT, VAR) 
    hists_mc_data_sub_SR_num = hists_mc_data_sub_AR_num *  fakefator_hist_num # SR = AR * FF = AR * (ss_iso/ss_noniso)

    # combine uncertainties and store values in bare arrays
    hists_mc_data_sub_SR_num_values = hists_mc_data_sub_SR_num()
    hists_mc_data_sub_SR_num_variances = hists_mc_data_sub_SR_num(sn.UP, sn.ALL, unc=True)**2

    print(f">>>>>>>>>>>>>>>>>>>>> Data - MC subtraction values: \n{hists_mc_data_sub_SR_num_values}")
    print(f">>>>>>>>>>>>>>>>>>>>> Data - MC subtraction variances: \n{hists_mc_data_sub_SR_num_variances}") 

    # Create hist from the subtraction
    hists_mc_data_sub_SR.view().value[...] = hists_mc_data_sub_SR_num_values
    hists_mc_data_sub_SR.view().variance[...] = hists_mc_data_sub_SR_num_variances

    # Plot
    plot_hist_cat(hists_data_cat_SR, hists_mc_cat_SR,category, hists_mc_data_sub_SR)
