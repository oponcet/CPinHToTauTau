'''
Author : @oponcet
Date : 2021-06-01
Script to get the histograms from the pickle files and convert them to ROOT histograms and save them in a ROOT file per region

Output: Example 
├── outputs/                    # Store the output JSON and ROOT files
│   ├── pi_1/
│   ├── pi_1_Oj.root 
│       ├── A/
│           ├── hcand_1_pt/
│               ├── THStack hcand_1_pt                # THStack of the histograms for the variable pt1 with all datasets    
│               ├── TH1D hcand_1_pt_data_tau_C        # TH1D of the histogram for the variable pt1 with the dataset tau_C
│               ... 
│               ├── TH1D hcand_1_pt_datasets          # TH1D of the histogram for the variable pt1 with all datasets     
│           ├── hcand_1_pt__met_var_qcd_h1/
│               ├── TH2D hcand_1_pt__met_var_qcd_h1_data_tau_C        # TH1D of the histogram for the variable pt1 with the dataset tau_C
                ... 
│               ├── TH2D hcand_1_pt__met_var_qcd_h1_datasets          # TH1D of the histogram for the variable pt1 with all datasets         
│       ├── B # same for B
│       ├── C
│       ├── D
│       ├── A0
│       ├── B0
│       ├── C0
│       ├── D0

'''

import pickle
import hist
import numpy as np
import scinum as sn
from hist import Hist
import json
import ROOT
import os



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


def get_all_hist(eos_path, task, dataset_data, dataset_mc, hist_path, var, cat_id):
    ''' 
    Load all histograms from the pickle files for data and MC samples
    Returns: dictionnary of histograms for data and MC samples
    '''
    hist_path_var = hist_path + "hist__" + var + ".pickle"
    ## Load data histogram from the pickle file
    hists_data = {}
    for dataset in dataset_data:
        pickle_file_path = eos_path + task + dataset + "/" + hist_path_var
        # print(f"Loading histogram data from: {pickle_file_path}")
        hists_data[dataset] = load_histogram_data(pickle_file_path)
        # print(f"hist_data : {hists_data[dataset]}")

        category_axis = hists_data[dataset].axes['category']
        if cat_id not in category_axis:
            # print(f"Warning: The category {cat_id} is not found in the histogram {dataset}.")
            hists_data[dataset] = None
            continue
        else: 
            cat_index = hists_data[dataset].axes['category'].index(cat_id)
            hists_data[dataset] = hists_data[dataset][cat_index, :, :, :].project(var)
        # print(f"hist_data : {hists_data[dataset]}")

    ## Load MC histogram from the pickle file
    hists_mc = {}
    for dataset in dataset_mc:

        pickle_file_path = eos_path + task + dataset + "/" + hist_path_var
        # print(f"Loading histogram data from: {pickle_file_path}")
        hists_mc[dataset] = load_histogram_data(pickle_file_path)
        # print(f"hist_mc : {hists_mc[dataset]}")

        category_axis = hists_mc[dataset].axes['category']
        if cat_id not in category_axis:
            # print(f"Warning: The category {cat_id} is not found in the histogram {dataset}.")
            hists_mc[dataset] = None
            continue
        else:
            cat_index = hists_mc[dataset].axes['category'].index(cat_id)
            hists_mc[dataset] = hists_mc[dataset][cat_index, :, :, :].project(var)
            #hists_mc[dataset] = hists_mc[dataset][cat_index, :, :, :].project(var)[40j:200j:1j]

        # print(f"hist_mc : {hists_mc[dataset]}

    # print(f">>>>>>>>>>>>>>>>>>>>> Data histograms: \n{hists_data}")

    # print(f">>>>>>>>>>>>>>>>>>>>> MC histograms: \n{hists_mc}")

    return hists_data, hists_mc


def main(config_path, dm):
    # Load the JSON configuration file containing the paths, categories, and variables
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


    # For 1D histograms
    vars1D = ["hcand_1_pt"]
    for var in variables:
        var1 = var['var1']
        vars1D.append(var1)
    
    print(f"Variables 1D : {vars1D}")
    
    # For 2D histograms
    vars2D = []
    for var in variables:
        var1 = var['var1']
        var2 = var['var2']
        vars2D.append([var1, var2])
    
    print(f"Variables 2D: {vars2D}")

    print("DM = ", dm)

    # Create the output directory if it does not exist
    if not os.path.exists('outputs'):
        os.makedirs('outputs')
    
    output_path = f'script_FF/fake_factor_derivation/outputs/{dm}'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Loop over the njet
    for njet in categories[dm].keys():
        print(f"njet: {njet}")

        # Create the ROOT file in the output directory  
        output_file = ROOT.TFile(f"{output_path}/{dm}_{njet}.root", "RECREATE")
        print(f"Output file: {output_file}")

        # Create TDirectory for each region in the ROOT file
        for cat_group in categories[dm][njet].keys(): # ABCD or A0B0C0D0
            for cat in categories[dm][njet][cat_group].keys(): # A, B, C, D, A0, B0, C0, D0
                # print(f"Processing cat: {cat }")
                cat_dir = output_file.mkdir(cat)
                cat_dir.cd()

                cat_id = categories[dm][njet][cat_group][cat]
                print(f"cat_id: {cat_id}")

                # Iterate over the variables 
                for var1D in vars1D:
                    # Create a directory for the var1D within the ROOT file
                    var1D_dir = cat_dir.mkdir(var1D)
                    var1D_dir.cd()
                
                    # Save all datasets in TH1D histograms

                    # Get the hist from the pickle file
                    hists_data, hists_mc = get_all_hist(eos_path, task, dataset_data, dataset_mc, hist_path_base, var1D, cat_id)

                    # Convert each dataset hist to TH1D and save it in the ROOT file
                    for dataset in dataset_data:
                        hist_data = hists_data[dataset]

                        if hist_data is None:
                            # print(f"Warning: No histogram data found for dataset {dataset}")
                            continue

                        else: 
                            # Extract bin edges, contents, and uncertainties
                            bin_edges = hist_data.axes[0].edges
                            bin_contents = hist_data.values()
                            bin_uncertainties = hist_data.variances() ** 0.5
                            
                            # Create ROOT TH1D
                            n_bins = len(bin_contents)
                            th1d = ROOT.TH1D(dataset, var1D, n_bins, bin_edges[0], bin_edges[-1])
                            
                            # Fill TH1D with contents and uncertainties
                            for i, (content, uncertainty) in enumerate(zip(bin_contents, bin_uncertainties), start=1):
                                th1d.SetBinContent(i, content)
                                th1d.SetBinError(i, uncertainty)

                            th1d.Write()
                            # print(f"TH1D histogram saved for dataset {dataset}")

                    for dataset in dataset_mc:
                        hist_mc = hists_mc[dataset]

                        if hist_mc is None:
                            # print(f"Warning: No histogram data found for dataset {dataset}")
                            continue

                        else: 
                            # Extract bin edges, contents, and uncertainties
                            bin_edges = hist_mc.axes[0].edges
                            bin_contents = hist_mc.values()
                            bin_uncertainties = hist_mc.variances() ** 0.5
                            
                            # Create ROOT TH1D
                            n_bins = len(bin_contents)
                            th1d = ROOT.TH1D(dataset, var1D, n_bins, bin_edges[0], bin_edges[-1])
                            
                            # Fill TH1D with contents and uncertainties
                            for i, (content, uncertainty) in enumerate(zip(bin_contents, bin_uncertainties), start=1):
                                th1d.SetBinContent(i, content)
                                th1d.SetBinError(i, uncertainty)

                            th1d.Write()
                            # print(f"TH1D histogram saved for dataset {dataset}")
    
        output_file.Close()

if __name__ == "__main__":
    # Define the path to the configuration file

    # dms = ['pi_1', 'rho_1', 'a1dm2_1', 'a1dm10_1', 'a1dm11_1']
    dms = ['pi_1']

    for dm in dms:
        config_path = f'script_FF/fake_factor_derivation/inputs/fake_factors_{dm}.json'

        # Run the main function
        main(config_path, dm)
