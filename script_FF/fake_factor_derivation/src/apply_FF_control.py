'''
Author : @oponcet
Date : 06-12-2024
Script to apply pt-depedent Fake Factors to distribution of any other variable 
(met_var_qcd_h1 for instance). 
Apply FF to B data_min_mc and use it in A region as Fake jets -> tau_h 

Input: 
    - inputs/inputs_json/region.json (json file with the region configuration)
    - outputs/FakeFactors/dm/region.json (json file with the fake factors) het the function to apply 
    - inputs/inputs_rootfile/dm/region.root/A/met_var_qcd_h1_hcand_1_pt/data_minus_mc (TH2D)
    - inputs/inputs_rootfile/dm/region.root/B/met_var_qcd_h1_/data_minus_mc (TH1D)

Output:
    - outputs/outputs_applyFF/dm/region.root/A/met_var_qcd_h1/fakes_jets (TH1D)
                                                             /canvas_FF (THstack)  
    - outputs/outputs_applyFF/dm/region_FF_apply_catA.png (png/pdf)                                    

'''
import os
import json
import ROOT

def apply_fake_factor(input_file, catA, catB, dm, njet,var2D):
    """
    Parameters:
    - input_file: Path to input ROOT file for region A. 
    - catA: category A (numerator) -> TH1D is in catA/variable/data_minus_mc
    - catB: category B (denominator) -> TH1D is in catB/variable/data_minus_mc
    - dm: decay mode
    - njet: number of jets
    - var2D: variable to apply the fake factor

    """
    # Define output paths
    output_root_file_path = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/dm/{dm}/{dm}_{njet}.root'
    output_png_file = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/dm/{dm}/{dm}_{njet}.png'

     # Ensure directories exist if not create them
    os.makedirs(os.path.dirname(output_root_file_path), exist_ok=True)
    

    # Load input ROOT files and histograms
    root_file = ROOT.TFile.Open(input_file, "READ")

    
    if not root_file or root_file.IsZombie():
        raise FileNotFoundError("Input ROOT files could not be opened.")

    var1 = var2D[0]
    var2 = var2D[1]

    # Get the 2D histograms
    hist_b = root_file.Get(f'{catB}/{var1}-{var2}/data_minus_mc').Clone()

    hist_b.SetDirectory(0) # Detach histogram from ROOT file to avoid deletion

    if not hist_b:
        raise KeyError(f"Histogram '{catA}/{var1}-{var2}/data_minus_mc not found in one of the input files.")

    # close the input ROOT file
    root_file.Close()


    # Get the fake factor
    fake_factor_file = f'script_FF/fake_factor_derivation/outputs/FakeFactors/{dm}/{dm}_{njet}.root'


    # Load input ROOT files and histograms
    root_fake_factor_file = ROOT.TFile.Open(fake_factor_file, "READ")
    print(f"Fake factor file: {fake_factor_file}")

    # Get the TF1 fake factor from the fit 
    FF_function = root_fake_factor_file.Get(f'FitResult')

    print(f"Fake factor function: {FF_function}")

    # ## Apply FF_function to hist_b
    hist_a_fake = hist_b.Clone()
    
    for i in range(hist_a_fake.GetNbinsX()):
        for j in range(hist_a_fake.GetNbinsY()):
            pt = hist_a_fake.GetYaxis().GetBinCenter(i)
            eta = hist_a_fake.GetXaxis().GetBinCenter(j)
            fake_factor = FF_function.Eval(pt)
            hist_a_fake.SetBinContent(i,j,hist_a_fake.GetBinContent(i,j)*fake_factor)


    ## Project the 2D histogram ton the hcand_1_pt axis
    hist_a_fake_proj = hist_a_fake.ProjectionX() 

    # set name of the histogram
    hist_a_fake_proj.SetName("fakes_jets")

    # Save the output ROOT file
    output_root_file = ROOT.TFile(output_root_file_path, "RECREATE")
    output_root_file.cd()
    hist_a_fake_proj.Write("fakes_jets")
    output_root_file.Close()

    ##########################################################
    ### Build hist stack and data_minus_mc_wfake histogram ###
    ##########################################################

    # Load input ROOT files and histograms
    root_file = ROOT.TFile.Open(input_file, "READ")

    # Get the stack canvas
    canvas_a_stack = root_file.Get(f'{catA}/{var1}/canvas').Clone()
    root_file.Close()

    # Open the output ROOT file
    output_root_file = ROOT.TFile(output_root_file_path, "UPDATE")
    output_root_file.cd()

    # Get THStack mc_stacks in canvas and data
    data_hist = canvas_a_stack.GetPrimitive("data_hist")  # Ensure this matches your histogram name
    mc_stack = canvas_a_stack.GetPrimitive('mc_stack')

    # Add the fake jets to the stack
    mc_stack.Add(hist_a_fake_proj)

    # Save the canvas to the output file
    canvas_a_stack.Write() 

    ##########################################################
    ### Build hist stack and data_minus_mc_wfake histogram ###
    ##########################################################
    # Compute the total MC histogram
    mc_total = mc_stack.GetStack().Last().Clone()  # Get the total MC histogram from the stack

    # # Compute the "Data - MC with fake jets" histogram
    # data_minus_mc_wfake = data_hist.Clone("data_minus_mc_wfake")
    # data_minus_mc_wfake.Add(mc_total, -1)


    # Create a ratio plot
    ratio = data_hist.Clone("ratio")
    ratio.Divide(mc_total)

    print("ratio: ", ratio)

    # Customize the ratio plot
    # Don't show the ratio plot title
    ratio.SetTitle("")
    # ratio.SetTitle("Data / MC with Fake Jets")
    ratio.GetYaxis().SetTitle("Ratio")
    ratio.GetYaxis().SetNdivisions(505)
    # ratio.GetYaxis().SetRangeUser(0.5, 1.5)  # Adjust range as needed
    ratio.GetXaxis().SetTitle("Variable")    # Replace with appropriate axis label

    print("creating canvas")
    # Create a new canvas and draw the plots
    output_canvas_ratio = ROOT.TCanvas("output_canvas_ratio", "Data-MC Ratio", 800, 600)
    output_canvas_ratio.Divide(1, 2)

    # Upper pad for the histograms
    upper_pad = output_canvas_ratio.cd(1)
    upper_pad.SetPad(0, 0.3, 1, 1) # (xlow, ylow, xup, yup)
    upper_pad.SetBottomMargin(0.04)
    data_hist.SetMarkerStyle(20)
    data_hist.Draw("E")
    mc_stack.Draw("HIST SAME")
    data_hist.Draw("E SAME")
    
    # Lower pad for the ratio plot
    lower_pad = output_canvas_ratio.cd(2)
    lower_pad.SetPad(0, 0, 1, 0.3)
    lower_pad.SetTopMargin(0.02)
    lower_pad.SetBottomMargin(0.3)
    ratio.Draw("E")


    # Save the canvas to the output file
    output_canvas_ratio.Write()

    # Save the canvas as a PNG file
    output_canvas_ratio.SaveAs(output_png_file)

    # Also save the ratio in a TH1D
    ratio_hist = ratio.Clone("ratio_hist")
    ratio_hist.Write("ratio_hist")

    output_root_file.Close()








