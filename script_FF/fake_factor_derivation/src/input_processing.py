'''
Author : @oponcet
Date : 06-12-24
Script to get the histograms from the pickle files and convert them to ROOT histograms and save them in a ROOT file per region

Input: json file containing the paths to the pickle files, the list the datasets, categories and variables

Output: Example 
├── inputs/inputs_rootfile                    # Store the output JSON and ROOT files
│   ├── pi_1/
│   ├── pi_1_Oj.root 
│       ├── A/
│           ├── hcand_1_pt/
│               ├── THStack hcand_1_pt                                # THStack of the histograms for the variable pt1 with all datasets    
│               ├── TH1D data_tau_C                                   # TH1D of the histogram for the variable pt1 with the dataset tau_C
│               ... 
│               ├── TH1D datasets                                     # TH1D of the histogram for the variable pt1 with all datasets     
│           ├── met_var_qcd_h1_met_var_qcd_h1/
│               ├── TH2D data_tau_C                                   # TH1D of the histogram for the variable pt1 with the dataset tau_C
│               ... 
│               ├── TH2D datasets                                     # TH1D of the histogram for the variable pt1 with all datasets         
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
import json
import ROOT
from ROOT import TCanvas, THStack, TLegend, TColor
import os
from plotting import *
import cmsstyle as CMS

def hist_to_num(h: hist.hist, unc_name=str(sn.DEFAULT)) -> sn.Number:
    """
    Convert hist object to sn.Number containing values and uncertainties.
    """
    return sn.Number(h.values(), {unc_name: h.variances()**0.5})


def load_histogram_data(pickle_file_path):
    """
    Load histogram data from a pickle file.
    """
    try:
        with open(pickle_file_path, 'rb') as f:
            return pickle.load(f)
    except (FileNotFoundError, Exception) as e:
        print(f"Error loading {pickle_file_path}: {e}")
        return None


def get_histograms(eos_path, task, dataset_data, dataset_mc, hist_path, var, cat_id, twoD=False):
    """
    Load histograms for data and MC samples and filter by category.
    """
    hist_path_var = hist_path + "hist__" + var + ".pickle"

    # Load data histograms
    hists_data = load_datasets_histograms(eos_path, task, dataset_data, hist_path_var, cat_id, var, twoD)
    # Load MC histograms
    hists_mc = load_datasets_histograms(eos_path, task, dataset_mc, hist_path_var, cat_id, var, twoD)
    
    
    return hists_data, hists_mc


def load_datasets_histograms(eos_path, task, datasets, hist_path_var, cat_id, var, twoD=False):
    """
    Helper function to load histograms for a list of datasets.
    """
    hists = {}

    # print("var ", var)
    for dataset in datasets:
        pickle_file_path = eos_path + task + dataset + "/" + hist_path_var
        hist_data = load_histogram_data(pickle_file_path)
        
        if hist_data:
            category_axis = hist_data.axes['category']
            if cat_id not in category_axis:
                hists[dataset] = None
            else:
                cat_index = hist_data.axes['category'].index(cat_id)
                if twoD == True:
                    var1 = var.split('-')[0]
                    var2 = var.split('-')[1]
                    hists[dataset] = hist_data[cat_index, :, :, :, :].project(var1, var2)
                else:
                    hists[dataset] = hist_data[cat_index, :, :, :].project(var)
        else:
            hists[dataset] = None
    return hists


def create_th1d_histogram(dataset, var, bin_edges, bin_contents, bin_uncertainties):
    """
    Create a ROOT TH1D histogram from bin contents and uncertainties.
    """
    n_bins = len(bin_contents)
    th1d = ROOT.TH1D(dataset, var, n_bins, bin_edges[0], bin_edges[-1])
    
    for i, (content, uncertainty) in enumerate(zip(bin_contents, bin_uncertainties), start=1):
        th1d.SetBinContent(i, content)
        th1d.SetBinError(i, uncertainty)

    return th1d


def save_histograms_to_root(hists_data, hists_mc, var, cat_dir):
    """
    Convert histograms to ROOT TH1D format, save them in the ROOT file,
    and create a stack plot for MC and a summed histogram for data.
    """
    # Create a list to store TH1D histograms for the stack plot
    mc_histograms = []
    data_histogram = None

    mc_sum_histogram = None
    
    # Save data histograms and sum them
    for dataset, hist_data in hists_data.items():
        if hist_data:
            bin_edges = hist_data.axes[0].edges
            bin_contents = hist_data.values()
            bin_uncertainties = hist_data.variances() ** 0.5
            th1d = create_th1d_histogram(dataset, var, bin_edges, bin_contents, bin_uncertainties)
            th1d.Write()
            
            # Sum the data histograms
            if data_histogram is None:
                data_histogram = th1d.Clone(f"{dataset}_data")
            else:
                data_histogram.Add(th1d)
            # rename the data histogram
            data_histogram.SetName("data_hist")

    # Save MC histograms and add them to the stack
    for dataset, hist_mc in hists_mc.items():
        if hist_mc:
            bin_edges = hist_mc.axes[0].edges
            bin_contents = hist_mc.values()
            bin_uncertainties = hist_mc.variances() ** 0.5
            th1d = create_th1d_histogram(dataset, var, bin_edges, bin_contents, bin_uncertainties)
            th1d.Write()

            # Add to the stack of MC histograms
            mc_histograms.append(th1d)

            # Sum the MC histograms
            if mc_sum_histogram is None:
                mc_sum_histogram = th1d.Clone(f"{dataset}_mc")
            else:
                mc_sum_histogram.Add(th1d)

    # Create data minus MC histogram
    if data_histogram and mc_sum_histogram:
        data_minus_mc = data_histogram.Clone(f"data_minus_mc")
        data_minus_mc.Add(mc_sum_histogram, -1)
        # set to zero the negative values
        for i in range(data_minus_mc.GetNbinsX()):
            if data_minus_mc.GetBinContent(i + 1) < 0:
                data_minus_mc.SetBinContent(i + 1, 0)
        data_minus_mc.Write()


    # Create a stack plot and a summary plot
    create_stack_plot_and_summary(mc_histograms, data_histogram, var, cat_dir)


def create_stack_plot_and_summary(mc_histograms, data_histogram, var, cat_dir, colors=None):
    """
    Create a stack plot for the MC histograms and a summary plot for Data vs MC.
    Also allows customization of colors for different MC processes.

    Parameters:
    - mc_histograms: List of MC histograms.
    - data_histogram: Data histogram.
    - var: Variable name to label the plot.
    - cat_dir: Directory where the plot and ROOT file will be saved.
    - colors: Dictionary containing colors for different MC processes (optional).
    """
    # Create a canvas for the Data vs MC plot
    # canvas = ROOT.TCanvas("canvas", "Data vs MC", 800, 600)
    CMS.SetExtraText("Private work")
    CMS.SetCmsTextFont(52)
    CMS.SetCmsTextSize(0.75*0.76)
    
    # canvas = CMS.cmsCanvas('Data/MC',40,200,0,1,"Energy [GeV]","Events/5 GeV",square=CMS.kSquare,extraSpace=0.05,iPos=0)
    canvas = CMS.cmsCanvas('canvas', 0, 1, 0, 1, '', '', square = CMS.kSquare, extraSpace=0.01, iPos=0)

    # Create stack for MC histograms
    stack = ROOT.THStack("mc_stack", "MC Stack")

    # Add each MC histogram to the stack with custom colors
    for i, hist in enumerate(mc_histograms):
        # Set custom fill color based on the given dictionary
        color, label = assign_color_and_label(hist.GetName()) 
        hist.SetFillColor(ROOT.TColor.GetColor(color))
        hist.SetLineColor(ROOT.kBlack)  # Set outline color to black
        hist.SetLineWidth(1)
        # Add histogram to stack
        stack.Add(hist)

    # Draw data histogram
    if data_histogram:
        data_histogram.SetLineColor(ROOT.kBlack)
        data_histogram.SetMarkerStyle(20)
        data_histogram.SetMarkerColor(ROOT.kBlack)
        data_histogram.Draw("E1")

    # Draw the stack on top of the data histogram
    stack.Draw("HIST SAME")
    
    # Set axis titles
    data_histogram.GetXaxis().SetTitle(var)
    data_histogram.GetYaxis().SetTitle("Events/5 GeV")

   # Add a legend
    legend = CMS.cmsLeg(0.7, 0.7, 0.9, 0.9, textSize=0.02)
    # legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # Position can be adjusted
    for hist in mc_histograms:
        # Get color and label for each histogram from the assign_color_and_label function
        color, label = assign_color_and_label(hist.GetName())
        # add entry if label not already in the legend
        if label not in [entry.GetLabel() for entry in legend.GetListOfPrimitives()]:
            legend.AddEntry(hist, label, "f")
    
    if data_histogram:
        legend.AddEntry(data_histogram, "Data", "lep")
    
    legend.Draw()

   

    cat = cat_dir.GetName()
    # print("cat ", cat)  

    if var == "hcand_1_pt":
        # remove stat box
        # ROOT.gStyle.SetOptStat(0)

        data_histogram.GetXaxis().SetRangeUser(40, 200)
        data_histogram.GetXaxis().SetTitle("p_{T}^{#tau} [GeV]")

        # Update the canvas to ensure everything is rendered properly
        stack.Modified() # Update the stack 
        canvas.Modified()
        canvas.Update()
        

        
        # Save the canvas to a file
        outdir_plot = f"script_FF/fake_factor_derivation/inputs/inputs_rootfile_ffapply/{dm}/plots"
        if not os.path.exists(outdir_plot):
            os.makedirs(outdir_plot)
        canvas.SaveAs(f"{outdir_plot}/input_stack_hcand_1_pt_{cat}.pdf")

    # Optionally, write the canvas to the ROOT file
    canvas.Write()


def save_histograms_to_root2D(hists_data, hists_mc, var1, var2, cat_dir):
    """
    Convert histograms to ROOT TH2D format and save them in the ROOT file.
    """
    # Sum of data histograms and mc 
    data_sum_histogram = None
    mc_sum_histogram = None

    # Save data histograms
    for dataset, hist_data in hists_data.items():
        if hist_data:
            # Extract bin edges, contents, and uncertainties for the 2D histogram
            bin_edges_x = hist_data.axes[0].edges
            bin_edges_y = hist_data.axes[1].edges
            bin_contents = hist_data.values()
            bin_uncertainties = hist_data.variances() ** 0.5

            th2d = create_th2d_histogram(dataset, var1, var2, bin_edges_x, bin_edges_y, bin_contents, bin_uncertainties)
            th2d.Write()
            # Add to sum 
            if data_sum_histogram is None:
                data_sum_histogram = th2d.Clone(f"{dataset}_data")
            else:
                data_sum_histogram.Add(th2d)

    # Save MC histograms
    for dataset, hist_mc in hists_mc.items():
        if hist_mc:
            # Extract bin edges, contents, and uncertainties for the 2D histogram
            bin_edges_x = hist_mc.axes[0].edges
            bin_edges_y = hist_mc.axes[1].edges
            bin_contents = hist_mc.values()
            bin_uncertainties = hist_mc.variances() ** 0.5

            th2d = create_th2d_histogram(dataset, var1, var2, bin_edges_x, bin_edges_y, bin_contents, bin_uncertainties)
            th2d.Write()
            # Add to sum
            if mc_sum_histogram is None:
                mc_sum_histogram = th2d.Clone(f"{dataset}_mc")
            else:
                mc_sum_histogram.Add(th2d)

    # Create data minus MC histogram
    if data_sum_histogram and mc_sum_histogram:
        data_minus_mc = data_sum_histogram.Clone(f"data_minus_mc")
        data_minus_mc.Add(mc_sum_histogram, -1)
        # set to zero the negative values
        for i in range(data_minus_mc.GetNbinsX()):
            for j in range(data_minus_mc.GetNbinsY()):
                if data_minus_mc.GetBinContent(i + 1, j + 1) < 0:
                    data_minus_mc.SetBinContent(i + 1, j + 1, 0)
       
        data_minus_mc.Write()

    

def create_th2d_histogram(dataset, var1, var2, bin_edges_x, bin_edges_y, bin_contents, bin_uncertainties):
    """
    Create a ROOT TH2D histogram from bin contents and uncertainties.
    """
    n_bins_x = len(bin_edges_x) - 1
    n_bins_y = len(bin_edges_y) - 1
    th2d = ROOT.TH2D(dataset, f"{var1}-{var2}", n_bins_x, bin_edges_x, n_bins_y, bin_edges_y)

    for i in range(n_bins_x):
            for j in range(n_bins_y):
                th2d.SetBinContent(i + 1, j + 1, bin_contents[i, j])
                th2d.SetBinError(i + 1, j + 1, bin_uncertainties[i, j])

    return th2d
                

def prepare_output_directory(dm):
    """
    Prepare the output directory structure.
    """
    output_path = f'script_FF/fake_factor_derivation/inputs/inputs_rootfile/{dm}'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    return output_path


def process_categories(categories, dm, njet, output_file, vars1D, vars2D, eos_path, task, dataset_data, dataset_mc, hist_path_base):
    """
    Process each category and save the corresponding histograms.
    """
    print(categories[dm][njet])
    for group, cat_dict in categories[dm][njet].items():
        print(f"Group: {group}")
        print(f"Categories: {cat_dict}")
        # for cat_group, cat_dict in njet_categories.items():
        for cat, cat_id in cat_dict.items():
            cat_dir = output_file.mkdir(cat)
            cat_dir.cd()
            print(f"Processing category: {cat} with id {cat_id}")
            
            for var1D in vars1D:
                # if var1D already processed, skip
                if cat_dir.GetListOfKeys().Contains(var1D):
                    continue
                
                var1D_dir = cat_dir.mkdir(var1D)
                var1D_dir.cd()

                hists_data, hists_mc = get_histograms(eos_path, task, dataset_data, dataset_mc, hist_path_base, var1D, cat_id, twoD=False)
                save_histograms_to_root(hists_data, hists_mc, var1D, cat_dir)
                
            cat_dir.cd()
            for var2D in vars2D:
                print(f"Processing 2D variable: {var2D}")
                var2D_name = f"{var2D[0]}-{var2D[1]}"
                var2D_dir = cat_dir.mkdir(var2D_name)
                var2D_dir.cd()

                hists_data, hists_mc = get_histograms(eos_path, task, dataset_data, dataset_mc, hist_path_base, var2D_name, cat_id, twoD=True)
                save_histograms_to_root2D(hists_data, hists_mc, var2D[0], var2D[1], cat_dir)

def sum_data_minus_mc_across_njets(dm):
    """
    Sum the data minus MC histograms across the njet categories.
    """

    # Define the input and output paths
    output_path = f'script_FF/fake_factor_derivation/inputs/inputs_rootfile/{dm}'
    output_file = ROOT.TFile(f"{output_path}/{dm}_allnjet.root", "RECREATE")
    print(f"Output file: {output_file.GetName()}")

    for cat in ["A", "B", "C", "D", "A0", "B0", "C0", "D0"]:
        # Reset data_minus_mc_sum for each category
        data_minus_mc_sum = None
        print(f"Processing category: {cat}")
        
        for njet in ["has_0j", "has_1j", "has_2j"]:
            input_file_name = f"{output_path}/{dm}_{njet}.root"
            input_file = ROOT.TFile(input_file_name, "READ")
            if not input_file or input_file.IsZombie():
                print(f"Error: Could not open file {input_file_name}")
                continue

            # Navigate to the directory
            dir_name = f"tautau__real_1__had{cat}__{njet}__{dm}"
            cat_dir_in = input_file.Get(dir_name)
            if not cat_dir_in:
                print(f"Warning: Directory '{dir_name}' not found in {input_file_name}")
                input_file.Close()
                continue

            var_dir = cat_dir_in.Get("hcand_1_pt")
            if not var_dir:
                print(f"Warning: Directory 'hcand_1_pt' not found in {cat_dir_in.GetName()}")
                input_file.Close()
                continue

            # Retrieve the histogram
            data_minus_mc = var_dir.Get("data_minus_mc")
            if not data_minus_mc:
                print(f"Warning: Histogram 'data_minus_mc' not found in {var_dir.GetName()}")
                input_file.Close()
                continue

            print(f"Adding histogram dm {dm}, cat {cat}, njet {njet}")

            # Sum the histograms
            if data_minus_mc_sum is None:
                data_minus_mc_sum = data_minus_mc.Clone("data_minus_mc_sum")
                data_minus_mc_sum.SetDirectory(0)  # Detach from the file
            else:
                data_minus_mc_sum.Add(data_minus_mc)
                data_minus_mc_sum.SetDirectory(0)

            input_file.Close()

        # Write the summed histogram to the output file
        if data_minus_mc_sum:
            dir_name = f"tautau__real_1__had{cat}__allnjet__{dm}"
            if not output_file.GetDirectory(dir_name):
                output_file.mkdir(dir_name)

            output_file.cd(dir_name)
            # create directory for the variable
            var_dir =  output_file.GetDirectory(dir_name).mkdir("hcand_1_pt")
            output_file.cd(f"{dir_name}/hcand_1_pt")

            data_minus_mc_sum.Write("data_minus_mc")
            print(f"Summed histogram written for category: {cat}")
        else:
            print(f"No histograms found for category: {cat}")

    # Close the output file
    output_file.Close()
    print("All operations completed successfully.")



def create_canvas_stack_all_jets(dm, output_path):
    """
    Create a canvas for stacking histograms from all jet multiplicities and save as a PDF.
    This combines the histograms from multiple jet multiplicities (`0j`, `1j`, `2j`) into a single plot.

    Args:
        dm (str): The decay mode (e.g., 'pi_1').
        output_path (str): The path to save the canvas.
    """
    
    """
    Sum the data minus MC histograms across the njet categories.
    """

    # Define the input and output paths
    output_path = f'script_FF/fake_factor_derivation/inputs/inputs_rootfile/{dm}'
    output_file = ROOT.TFile(f"{output_path}/{dm}_allnjet.root", "UPDATE")
    print(f"Output file: {output_file.GetName()}")

    for cat in ["A", "B", "C", "D", "A0", "B0", "C0", "D0"]:
        # Reset data_minus_mc_sum for each category
        mc_stack_sum = None
        data_hist_sum = None
        print(f"Processing category: {cat}")
        
        for njet in ["has_0j", "has_1j", "has_2j"]:
            input_file_name = f"{output_path}/{dm}_{njet}.root"
            input_file = ROOT.TFile(input_file_name, "READ")
            if not input_file or input_file.IsZombie():
                print(f"Error: Could not open file {input_file_name}")
                continue

            # Navigate to the directory
            dir_name = f"tautau__real_1__had{cat}__{njet}__{dm}"
            cat_dir_in = input_file.Get(dir_name)
            if not cat_dir_in:
                print(f"Warning: Directory '{dir_name}' not found in {input_file_name}")
                input_file.Close()
                continue

            var_dir = cat_dir_in.Get("hcand_1_pt")
            if not var_dir:
                print(f"Warning: Directory 'hcand_1_pt' not found in {cat_dir_in.GetName()}")
                input_file.Close()
                continue

            # Retrieve the canvas
            canvas = var_dir.Get("canvas")
            if not canvas:
                print(f"Warning: canvas  not found in {var_dir.GetName()}")
                input_file.Close()
                continue

        
            # Retrive the THStack and data histogram in the canvas
            mc_stack = canvas.GetPrimitive("mc_stack")
            data_hist = canvas.GetPrimitive("data_hist")

            # Sum the mc histograms of the THStack
            if mc_stack_sum is None:
                mc_stack_sum = mc_stack.Clone("mc_stack_sum")
            else:
                for i in range(mc_stack.GetNhists()):
                    # Get the i-th histogram from mc_stack
                    hist = mc_stack.GetStack().At(i)
                    
                    # Add the histogram to the sum histogram (mc_stack_sum)
                    mc_stack_sum.GetStack().Last().Add(hist)

            # Do the same for the data histogram
            if data_hist_sum is None:
                data_hist_sum = data_hist.Clone("data_minus_mc_sum")
            else:
                data_hist_sum.Add(data_hist)

            data_hist_sum.SetDirectory(0)

            input_file.Close()

        # Write the in the output file
        if data_hist_sum:
            dir_name = f"tautau__real_1__had{cat}__allnjet__{dm}"
            if not output_file.GetDirectory(dir_name):
                output_file.mkdir(dir_name)

            output_file.cd(dir_name)
            # create directory for the variable
            var_dir =  output_file.GetDirectory(dir_name).mkdir("hcand_1_pt")
            output_file.cd(f"{dir_name}/hcand_1_pt")

            #### Create the canvas ### 

            CMS.SetExtraText("Private work")
            CMS.SetCmsTextFont(52)
            CMS.SetCmsTextSize(0.75*0.76)
            
            # canvas = CMS.cmsCanvas('Data/MC',40,200,0,1,"Energy [GeV]","Events/5 GeV",square=CMS.kSquare,extraSpace=0.05,iPos=0)
            canvas = CMS.cmsCanvas('canvas', 0, 1, 0, 1, '', '', square = CMS.kSquare, extraSpace=0.01, iPos=0)

            # Create stack for MC histograms
            stack = ROOT.THStack("mc_stack", "MC Stack")

            # Add each MC histogram to the stack with custom colors
            for i, hist in enumerate(mc_stack_sum.GetHists()):
                # Set custom fill color based on the given dictionary
                color, label = assign_color_and_label(hist.GetName()) 
                hist.SetFillColor(ROOT.TColor.GetColor(color))
                hist.SetLineColor(ROOT.kBlack)  # Set outline color to black
                hist.SetLineWidth(1)
                # Add histogram to stack
                stack.Add(hist)

            # Draw data histogram
            if data_hist_sum:
                data_hist_sum.SetLineColor(ROOT.kBlack)
                data_hist_sum.SetMarkerStyle(20)
                data_hist_sum.SetMarkerColor(ROOT.kBlack)
                data_hist_sum.Draw("E1")

            # rename to data_hist
            data_hist_sum.SetName("data_hist")
            # Draw the stack on top of the data histogram
            stack.Draw("HIST SAME")
            
            # Set axis titles
            data_hist_sum.GetXaxis().SetTitle("hcand_1_pt")
            data_hist_sum.GetYaxis().SetTitle("Events/5 GeV")

        # Add a legend
            legend = CMS.cmsLeg(0.7, 0.7, 0.9, 0.9, textSize=0.02)
            # legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # Position can be adjusted
            for hist in mc_stack_sum.GetHists():
                # Get color and label for each histogram from the assign_color_and_label function
                color, label = assign_color_and_label(hist.GetName())
                # add entry if label not already in the legend
                if label not in [entry.GetLabel() for entry in legend.GetListOfPrimitives()]:
                    legend.AddEntry(hist, label, "f")
            
            if data_hist_sum:
                legend.AddEntry(data_hist_sum, "Data", "lep")
            
            legend.Draw()

     

            data_hist_sum.GetXaxis().SetRangeUser(40, 200)
            data_hist_sum.GetXaxis().SetTitle("p_{T}^{#tau} [GeV]")

            # Update the canvas to ensure everything is rendered properly
            stack.Modified() # Update the stack 
            canvas.Modified()
            canvas.Update()
            canvas.Write()
            print(f"Summed histogram written for category: {cat}")
        else:
            print(f"No histograms found for category: {cat}")

    # Close the output file
    output_file.Close()
    print("All operations completed successfully.")









def main(config_path, dm):

    # ROOT configuration
    ROOT.gROOT.SetBatch(True) # Do not display canvas

    # Load the JSON configuration
    with open(config_path, 'r') as f:
        config = json.load(f)

    eos_path = config['paths']['eos_path']
    task = config['paths']['task']
    hist_path_base = config['paths']['hist_path']
    dataset_data = config['datasets']['data']
    dataset_mc = config['datasets']['mc']
    categories = config['categories']
    variables = config['variables']

    # For 1D histograms
    vars1D = ["hcand_1_pt"] + [var['var1'] for var in variables]
    # print(f"Variables 1D: {vars1D}")

    # vars2D = [[var['var1'], var['var2']] for var in variables]
    vars2D = []
    # print(f"Variables 2D: {vars2D}")

    # Prepare the output directory
    output_path = prepare_output_directory(dm)

    # Loop over the categories and process the histograms
    for njet in categories[dm].keys():
        output_file = ROOT.TFile(f"{output_path}/{dm}_{njet}.root", "RECREATE")
        # print(f"Output file: {output_file}")
        process_categories(categories, dm,njet, output_file, vars1D, vars2D, eos_path, task, dataset_data, dataset_mc, hist_path_base)
        output_file.Close()
    if dm == "pi_1":
            sum_data_minus_mc_across_njets(dm)
            create_canvas_stack_all_jets(dm, output_path)


if __name__ == "__main__":
    dms = ['pi_1', 'rho_1', 'a1dm2_1', 'a1dm10_1', 'a1dm11_1']
    # dms = ['pi_1']
    for dm in dms:
        config_path = f'script_FF/fake_factor_derivation/inputs/inputs_json/fake_factors_{dm}.json'
        main(config_path, dm)
