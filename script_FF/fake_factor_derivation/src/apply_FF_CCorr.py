'''
Author : @oponcet
Date : 12-12-2024
Script to apply pt dependent Fake Factors and met_var_qcd_h1 closure correction on pt 
distribution of the cat C to get D category.

Input:
    - FF root file (TH1D)
    - CCorr root file (TH1D)
    - 2D histogram of C met_var_qcd_h1, pt (TH2D)
    - 1D histogram of D pt (THStack) to build the final D histogram
                         
'''

import os
import ROOT
from plotting import *
import cmsstyle as CMS
from ROOT import TCanvas, THStack, TLegend, TColor
from plotting import assign_color_and_label

def apply_fake_factor_CCorr_CD(input_file, catC, catD, dm, njet):
    """
    Apply FF to C data_min_mc and use it in D region as Fake jets -> tau_h and produce
    the stack plot of the pt_1.

    Parameters:
    - input_file: Path to input ROOT file for region A.
    - catC: category C (numerator) -> TH1D is in catC/variable/data_minus_mc
    - catD: category D (denominator) -> TH1D is in catD/variable/data_minus_mc
    - dm: decay mode
    - njet: number of jets
    """

    output_root_file_path = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF_CCorr/dm/{dm}/{dm}_{njet}_hcand_1_pt_D.root'
    output_png_file = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF_CCorr/dm/{dm}/{dm}_{njet}_hcand_1_pt_D'

    # Ensure directories exist if not create them
    os.makedirs(os.path.dirname(output_root_file_path), exist_ok=True)

    #####################################################################
    ### Load input ROOT files and histograms that contains all histos ###
    #####################################################################

    # Load input ROOT files and histograms that contains all histos
    root_file = ROOT.TFile.Open(input_file, "READ")

    hist_c = root_file.Get(f'{catC}/met_var_qcd_h1-hcand_1_pt/data_minus_mc')

    if not hist_c :
        raise KeyError(f"Histogram '{catC}/met_var_qcd_h1-hcand_1_pt/data_minus_mc not found in one of the input files.")

    ## Detach the histograms from the ROOT file
    hist_c.SetDirectory(0)

    root_file.Close()
    

    #####################################################################
    ###  Load the fake factor and apply it to the histograms          ###
    #####################################################################

    # Get the fake factor
    fake_factor_file = f'script_FF/fake_factor_derivation/outputs/FakeFactors/{dm}/{dm}_{njet}.root'

    # Load input ROOT files and histograms
    root_fake_factor_file = ROOT.TFile.Open(fake_factor_file, "READ")

    # Get the TF1 fake factor from the fit
    FF_function = root_fake_factor_file.Get(f'FitResult')

    ## Detach the histograms from the ROOT file
    # FF_function.SetDirectory(0)

    root_fake_factor_file.Close()

    # Apply the fake factor to the histograms
    hist_2D_fake = hist_c.Clone()



    for i in range(hist_2D_fake.GetNbinsX()):
        met = hist_2D_fake.GetXaxis().GetBinCenter(i) # X axis is met_var_qcd_h1
        for j in range(hist_2D_fake.GetNbinsY()):
            pt = hist_2D_fake.GetYaxis().GetBinCenter(j) # Y axis is hcand_1_pt
            fake_factor = FF_function.Eval(pt) # Evaluate the fake factor at the given pt
            hist_2D_fake.SetBinContent(i,j,hist_2D_fake.GetBinContent(i,j)*fake_factor)


    #####################################################################
    ###  Load the closure correction and apply it to the histograms   ###
    #####################################################################

    # Get the closure correction
    CCorr_file = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/closure/dm/{dm}/{dm}_{njet}_met_var_qcd_h1.root' # script_FF/fake_factor_derivation/outputs/outputs_applyFF/closure/dm/pi_1/pi_1_has_2j_met_var_qcd.root
 

    # Load input ROOT files and histograms
    root_CCorr_file = ROOT.TFile.Open(CCorr_file, "READ")
    
    # Get the TF1 closure correction from the fit
    CCorr_function = root_CCorr_file.Get(f'FitResult')

    root_CCorr_file.Close()

    # Apply the closure correction to the histograms
    for i in range(hist_2D_fake.GetNbinsX()):
        met = hist_2D_fake.GetXaxis().GetBinCenter(i) # X axis is met_var_qcd_h1
        for j in range(hist_2D_fake.GetNbinsY()):
            pt = hist_2D_fake.GetYaxis().GetBinCenter(j)
            CCorr = CCorr_function.Eval(met) # Evaluate the closure correction at the given met
            hist_2D_fake.SetBinContent(i,j,hist_2D_fake.GetBinContent(i,j)*CCorr)


    #####################################################################
    ###  Project the 2D histogram on the hcand_1_pt axis              ###
    #####################################################################

   
    ## Project the 2D histogram ton the hcand_1_pt axis
    hist_2D_fake_proj = hist_2D_fake.ProjectionY()  # Project the 2D histogram on the hcand_1_pt axis

    # set name of the histogram
    hist_2D_fake_proj.SetName("fakes_jets")
    hist_2D_fake.SetName("fakes_jets_2D")
    # Save the output ROOT file
    output_root_file = ROOT.TFile(output_root_file_path, "RECREATE")
    output_root_file.cd()

    # Save th histograms
    hist_2D_fake.Write("fakes_jets_2D")
    hist_2D_fake_proj.Write("fakes_jets")

    # detatch the histograms from the ROOT file
    hist_2D_fake.SetDirectory(0)
    hist_2D_fake_proj.SetDirectory(0)
    
    output_root_file.Close()

    #######################################################################
    ###  Create stack plot of D with hist_c_fake as fake contribution   ###
    #######################################################################

    # Load input ROOT files and histograms
    root_file = ROOT.TFile.Open(input_file, "READ")

    # Get the stack canvas
    canvas_d_stack = root_file.Get(f'{catD}/hcand_1_pt/canvas').Clone()

    # Get THStack mc_stacks in canvas and data
    data_hist = canvas_d_stack.GetPrimitive("data_hist")  # Ensure this matches your histogram name

    mc_stack = canvas_d_stack.GetPrimitive('mc_stack')

    # Add the fake jets to the stack
    mc_stack.Add(hist_2D_fake_proj)

    # Detach the histograms from the ROOT file
    data_hist.SetDirectory(0)
    hist_2D_fake_proj.SetDirectory(0)

    # Close canvas and create a new one
    canvas_d_stack.Close()
    canvas_d_stack = CMS.cmsCanvas('Data-MC', 0, 1, 0, 1, '', '', square = CMS.kSquare, extraSpace=0.01, iPos=0)


    # rename the fake jets histogram
    hist_2D_fake_proj.SetName("fakes_jets")

    # Draw data
    data_hist.Draw("E ")
    mc_stack.Draw("HIST SAME")
    data_hist.Draw("E SAME")

    ### STYLE THE PLOT ###

    # legend
    # Add a legend
    legend = CMS.cmsLeg(0.7, 0.7, 0.9, 0.9, textSize=0.02)
    # legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # Position can be adjusted
    for hist in mc_stack:
        # Get color and label for each histogram from the assign_color_and_label function
        color, label = assign_color_and_label(hist.GetName())
        hist.SetFillColor(ROOT.TColor.GetColor(color))
        hist.SetLineColor(ROOT.kBlack)
        hist.SetLineWidth(1)
        hist.SetMarkerSize(0)
        # add entry if label not already in the legend
        if label not in [entry.GetLabel() for entry in legend.GetListOfPrimitives()]:
            legend.AddEntry(hist, label, "f")
    
    if data_hist:
        legend.AddEntry(data_hist, "Data", "lep")

    # Axis
    data_hist.GetYaxis().SetTitle("Events/5 GeV")
    data_hist.GetXaxis().SetTitle("p_{T} (GeV)")
    data_hist.GetXaxis().SetRangeUser(40, 200)
    

    # Save the canvas to the output file
    output_root_file = ROOT.TFile(output_root_file_path, "RECREATE")
    output_root_file.cd()
    canvas_d_stack.Write()


    # Save the canvas as a PNG file
    canvas_d_stack.SaveAs(f"{output_png_file}.pdf")


    #######################################################################
    ###   Create stack plot of D with hist_c_fake as fake with ratio    ###
    #######################################################################
    # Compute the total MC histogram
    mc_total = mc_stack.GetStack().Last().Clone()  # Get the total MC histogram from the stack

    output_png_file_ratio = f'{output_png_file}_ratio.pdf'
    
    
    # Create a ratio plot
    ratio = data_hist.Clone("ratio")
    ratio.Divide(mc_total)

    # Customize the ratio plot
    # Don't show the ratio plot title
    ratio.SetTitle("")
    # ratio.SetTitle("Data / MC with Fake Jets")
    ratio.GetYaxis().SetTitle("Ratio")
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetYaxis().SetRangeUser(0.5, 2.5)  # Adjust range as needed
    ratio.GetXaxis().SetTitle("p_{T} (GeV)")    # Replace with appropriate axis label

    print("creating canvas")
    # Create a new canvas and draw the plots
    output_canvas_ratio = CMS.cmsDiCanvas('Data-MC Ratio', 40, 200, 0, 1, 0, 2.5, 'p_{T} (GeV)', 'Events/5 GeV', 'Data/MC', square = CMS.kSquare, iPos=0, extraSpace=0.01)

    
    # output_canvas_ratio.Divide(1, 2)

    # Upper pad for the histograms
    upper_pad = output_canvas_ratio.cd(1)
    upper_pad.SetPad(0, 0.3, 1, 1) # (xlow, ylow, xup, yup)
    upper_pad.SetBottomMargin(0.02)
    data_hist.SetMarkerStyle(20)

    data_hist.Draw("E ")
    mc_stack.Draw("HIST SAME")
    data_hist.Draw("E SAME")

    # remove x axis title and labels
    data_hist.GetXaxis().SetTitle("")
    data_hist.GetXaxis().SetLabelSize(0)
    data_hist.GetXaxis().SetRangeUser(40, 200)
    data_hist.GetYaxis().SetTitle("Events/5 GeV")
    

    # legend
    # increase text size
    legend.SetTextSize(0.02)
    legend.Draw()
    
    # Lower pad for the ratio plot
    lower_pad = output_canvas_ratio.cd(2)
    lower_pad.SetPad(0, 0, 1, 0.27)
    lower_pad.SetTopMargin(0.036)
    # lower_pad.SetBottomMargin(0.3)
    # change label size
    ratio.GetYaxis().SetLabelSize(0.1)
    ratio.GetYaxis().SetTitleSize(0.1)
    ratio.GetXaxis().SetLabelSize(0.1)
    ratio.GetXaxis().SetTitleSize(0.1)
    ratio.Draw("E")

    # Draw dash line at y=1
    line = ROOT.TLine(40, 1, 200, 1)
    line.SetLineStyle(2)
    line.Draw()


    # Save the canvas to the output file
    output_canvas_ratio.Write()

    # Save the canvas as a PNG file
    output_canvas_ratio.SaveAs(output_png_file_ratio)

    # detatch the histograms from the ROOT file
    ratio.SetDirectory(0)

    output_root_file.Close()



#####################################################################################################
### Commented because I think that we cannot combine 3 two 2D histograms to get the 3D histogram  ###
### taking into account the correlation between the 3 variables                                   ###
#####################################################################################################

# def apply_fake_factor_CCorr_CD_var(input_file, catC, catD, dm, njet, var):
#     """
#     Apply FF to C data_min_mc and use it in D region as Fake jets -> tau_h and produce
#     the stack plot of the variable.

#     Parameters:
#     - input_file: Path to input ROOT file for region A.
#     - catC: category C (numerator) -> TH1D is in catC/variable/data_minus_mc
#     - catD: category D (denominator) -> TH1D is in catD/variable/data_minus_mc
#     - dm: decay mode
#     - njet: number of jets
#     - var: variable to plot
#     """

#     output_root_file_path = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF_CCorr/dm/{dm}/{dm}_{njet}_{var}_D.root'
#     output_png_file = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF_CCorr/dm/{dm}/{dm}_{njet}_{var}_D'

#     # Ensure directories exist if not create them
#     os.makedirs(os.path.dirname(output_root_file_path), exist_ok=True)


#     #####################################################################
#     ### Load input ROOT files and histograms that contains all histos ###
#     #####################################################################

#     # Load input ROOT files and histograms that contains all histos
#     root_file = ROOT.TFile.Open(input_file, "READ")
#     print("opening file ", input_file)

#     hist_c_qcd_pt1 = root_file.Get(f'{catC}/met_var_qcd_h1-hcand_1_pt/data_minus_mc')
#     print("hist_c_qcd_pt1", hist_c_qcd_pt1)

#     hist_c_var_pt1 = root_file.Get(f'{catC}/{var}-hcand_1_pt/data_minus_mc')
#     print("hist_c_var_pt1", hist_c_var_pt1)


#     hist_c_var_qcd = root_file.Get(f'{catC}/{var}-met_var_qcd_h1/data_minus_mc')

#     if not hist_c_qcd_pt1 or not hist_c_var_pt1 or not hist_c_var_qcd:
#         raise KeyError(f"Histogram '{catC}/met_var_qcd_h1-hcand_1_pt/data_minus_mc' or '{catC}/{var}-hcand_1_pt/data_minus_mc' or '{catC}/{var}-met_var_qcd_h1/data_minus_mc' not found in one of the input files.")
   
   
   
#     ## Detach the histograms from the ROOT file
#     hist_c_qcd_pt1.SetDirectory(0)
#     hist_c_var_pt1.SetDirectory(0)
#     hist_c_var_qcd.SetDirectory(0)

#     root_file.Close()

#     #####################################################################
#     ###  Load the fake factor and apply it to the histograms          ###
#     #####################################################################

#     # Get the fake factor
#     fake_factor_file = f'script_FF/fake_factor_derivation/outputs/FakeFactors/{dm}/{dm}_{njet}.root'

#     # Load input ROOT files and histograms
#     root_fake_factor_file = ROOT.TFile.Open(fake_factor_file, "READ")

#     # Get the TF1 fake factor from the fit
#     FF_function = root_fake_factor_file.Get(f'FitResult')

#     ## Detach the histograms from the ROOT file
#     # FF_function.SetDirectory(0)

#     root_fake_factor_file.Close()

#     #####################################################################
#     ###  Load the closure correction and apply it to the histograms   ###
#     #####################################################################

#     # Get the closure correction
#     CCorr_file = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/closure/dm/{dm}/{dm}_{njet}_met_var_qcd_h1.root' # script_FF/fake_factor_derivation/outputs/outputs_applyFF/closure/dm/pi_1/pi_1_has_2j_met_var_qcd.root
 

#     # Load input ROOT files and histograms
#     root_CCorr_file = ROOT.TFile.Open(CCorr_file, "READ")
    
#     # Get the TF1 closure correction from the fit
#     CCorr_function = root_CCorr_file.Get(f'FitResult')

#     root_CCorr_file.Close()

#     #####################################################################
#     ### Apply the fake factor and closure correction to the histogram ###
#     #####################################################################   

#     # ### Projection histograms
#     # hist_c_var_proj_hist_c_var_pt1 = hist_c_var_pt1.ProjectionX("hist_c_var_proj_hist_c_var_pt1")
#     # hist_c_var_proj_hist_c_var_qcd = hist_c_var_qcd.ProjectionX("hist_c_var_proj_hist_c_var_qcd")

#     # hist_c_pt1_proj_hist_c_var_pt1 = hist_c_var_pt1.ProjectionY("hist_c_pt1_proj_hist_c_var_pt1")
#     # hist_c_pt1_proj_hist_c_qcd_pt1 = hist_c_qcd_pt1.ProjectionY("hist_c_pt1_proj_hist_c_qcd_pt1")

#     # hist_c_qcd_proj_hist_c_var_qcd = hist_c_var_qcd.ProjectionY("hist_c_qcd_proj_hist_c_var_qcd")
#     # hist_c_qcd_proj_hist_c_qcd_pt1 = hist_c_qcd_pt1.ProjectionX("hist_c_qcd_proj_hist_c_qcd_pt1")


#     # ##### HERE big issue, because the projection are slightly different !!! 


#     # # save the histogram
#     # output_root_file = ROOT.TFile(output_root_file_path, "RECREATE")
#     # output_root_file.cd()
#     # hist_c_var_pt1.Write()
#     # hist_c_qcd_pt1.Write()
#     # hist_c_var_qcd.Write()
#     # hist_c_var_proj_hist_c_var_pt1.Write()
#     # hist_c_var_proj_hist_c_var_qcd.Write()
#     # hist_c_pt1_proj_hist_c_var_pt1.Write()
#     # hist_c_pt1_proj_hist_c_qcd_pt1.Write()
#     # hist_c_qcd_proj_hist_c_var_qcd.Write()
#     # hist_c_qcd_proj_hist_c_qcd_pt1.Write()


#     # ### Projection histograms 
#     hist_c_var_proj = hist_c_var_pt1.ProjectionX("hist_c_var_proj_hist_c_var_pt1")
#     hist_c_pt1_proj = hist_c_var_pt1.ProjectionY("hist_c_pt1_proj_hist_c_var_pt1")
#     hist_c_qcd_proj = hist_c_var_qcd.ProjectionY("hist_c_qcd_proj_hist_c_var_qcd")

#     # # save the histogram
#     output_root_file = ROOT.TFile(output_root_file_path, "RECREATE")
#     output_root_file.cd()
#     hist_c_var_proj.Write()
#     hist_c_pt1_proj.Write()
#     hist_c_qcd_proj.Write()


#     # Create the 1D histogram for the variable "var" based on its range and binning
#     hist_c_var = ROOT.TH1D("hist_c_var", f"Weighted {var} Histogram",
#                         hist_c_var_qcd.GetNbinsX(), hist_c_var_qcd.GetXaxis().GetXmin(),
#                         hist_c_var_qcd.GetXaxis().GetXmax())

#     # cret 3d histogram
#     hist_c_var_3D = ROOT.TH3D("hist_c_var_3D", f"Weighted {var} Histogram",
#                         hist_c_var_qcd.GetNbinsX(), hist_c_var_qcd.GetXaxis().GetXmin(),
#                         hist_c_var_qcd.GetXaxis().GetXmax(),
#                         hist_c_qcd_pt1.GetNbinsY(), hist_c_qcd_pt1.GetYaxis().GetXmin(),
#                         hist_c_qcd_pt1.GetYaxis().GetXmax(),
#                         hist_c_var_pt1.GetNbinsY(), hist_c_var_pt1.GetYaxis().GetXmin(),
#                         hist_c_var_pt1.GetYaxis().GetXmax())


#     # Loop through all bins in the 'var' (x-axis of hist_c_var_qcd)
#     for x_var in range(1, hist_c_var_qcd.GetNbinsX() + 1):  # 'var' bins
#         var_value = hist_c_var_qcd.GetXaxis().GetBinCenter(x_var)  # Variable value
#         for y_pt in range(1, hist_c_qcd_pt1.GetNbinsY() + 1):  # 'pt' bins (y-axis in pt histograms)
#             pt_value = hist_c_qcd_pt1.GetYaxis().GetBinCenter(y_pt)  # pt value
#             for x_qcd in range(1, hist_c_qcd_pt1.GetNbinsX() + 1):  # 'qcd' bins (x-axis in qcd histograms)
#                 qcd_value = hist_c_qcd_pt1.GetXaxis().GetBinCenter(x_qcd)  # qcd/met_var value

#                 # Retrieve bin content using `FindBin`-determined values
#                 bin_content_qcd_pt1 = hist_c_qcd_pt1.GetBinContent(hist_c_qcd_pt1.FindBin(qcd_value, pt_value))
#                 bin_content_var_pt1 = hist_c_var_pt1.GetBinContent(hist_c_var_pt1.FindBin(var_value, pt_value))
#                 bin_content_var_qcd = hist_c_var_qcd.GetBinContent(hist_c_var_qcd.FindBin(var_value, qcd_value))


#                 # Normalization
#                 norm_qcd = hist_c_var_proj.GetBinContent(hist_c_var_proj.FindBin(var_value))
#                 norm_pt1 = hist_c_pt1_proj.GetBinContent(hist_c_pt1_proj.FindBin(pt_value))
#                 norm_var = hist_c_qcd_proj.GetBinContent(hist_c_qcd_proj.FindBin(qcd_value))
    
#                 # print("norm_qcd", norm_qcd)
#                 # print("norm_pt1", norm_pt1)
#                 # print("norm_var", norm_var)
                

#                 # Skip empty or invalid bins to avoid unnecessary operations
#                 if bin_content_qcd_pt1 <= 0 or bin_content_var_pt1 <= 0 or bin_content_var_qcd <= 0:
#                     continue

#                 # Apply fake factor (FF) and closure correction (CCorr)
#                 fake_factor_weight = FF_function.Eval(pt_value)
#                 closure_corr_weight = CCorr_function.Eval(qcd_value)

#                 # Combined weight
#                 # combined_weight = fake_factor_weight * closure_corr_weight
#                 combined_weight = 1

#                 # Weighted contribution to the histogram for the variable
#                 weighted_value = bin_content_qcd_pt1 * bin_content_var_pt1 * bin_content_var_qcd * combined_weight
#                 sum_value = bin_content_qcd_pt1 + bin_content_var_pt1 + bin_content_var_qcd

#                 # # # # Normalize 
#                 # weighted_value = weighted_value / (norm_qcd * norm_pt1 * norm_var)
#                 weighted_value = weighted_value 

#                 # Find the corresponding bin in the 1D histogram for `var`
#                 var_bin = hist_c_var.FindBin(var_value)
#                 hist_c_var.AddBinContent(var_bin, weighted_value)

#                 # build 3d histogram
#                 hist_c_var_3D.SetBinContent(x_var, y_pt, x_qcd, weighted_value)
                



#     # print("saving hist_c_var, ", hist_c_var)

#     # Save the histogram
#     hist_c_var.Write()
#     hist_c_var_3D.Write()

#     # project the 3D histogram to get the 1D var histogram
#     hist_c_var_1D = hist_c_var_3D.ProjectionX("hist_c_var_proj_x")
#     hist_c_var_1D.Write()














#     # # Prepare output for dR distribution (weighted)
#     # hist_dR = hist_2D_DR_pt.ProjectionX("weighted_dR")  # Project dR axis from hist_2D_DR_pt
#     # hist_dR.Reset()  # Reset to fill with weighted values

#     # # Loop through all bins in dR for the two histograms
#     # for x_bin in range(1, hist_2D_DR_pt.GetNbinsX() + 1):  # Loop over dR bins
#     #     dR = hist_2D_DR_pt.GetXaxis().GetBinCenter(x_bin)  # dR bin center
        
#     #     for y_bin_pt in range(1, hist_2D_DR_pt.GetNbinsY() + 1):  # Loop over pT bins
#     #         pT = hist_2D_DR_pt.GetYaxis().GetBinCenter(y_bin_pt)  # Bin center for pT
            
#     #         # Get the corresponding var_met_qcd bin for this dR
#     #         for y_bin_met in range(1, hist_2D_DR_var_met_qcd.GetNbinsY() + 1):  # Loop over var_met_qcd bins
#     #             var_met_qcd = hist_2D_DR_var_met_qcd.GetYaxis().GetBinCenter(y_bin_met)  # Bin center for var_met_qcd

#     #             # Get the original bin contents for the given (dR, pT) and (dR, var_met_qcd)
#     #             bin_content_pt = hist_2D_DR_pt.GetBinContent(x_bin, y_bin_pt)  # From hist_2D_DR_pt
#     #             bin_content_met = hist_2D_DR_var_met_qcd.GetBinContent(x_bin, y_bin_met)  # From hist_2D_DR_var_met_qcd

#     #             if bin_content_pt == 0 or bin_content_met == 0:
#     #                 continue  # Skip empty bins

#     #             # Apply the Fake Factor (FF) and Closure Correction (C_Corr)
#     #             fake_factor = FF_function.Eval(pT)  # From pT dimension
#     #             CCorr = CCorr_function.Eval(var_met_qcd)  # From var_met_qcd dimension

#     #             # Weight the bin content
#     #             weight = fake_factor * CCorr
#     #             weighted_content = bin_content_pt * bin_content_met * weight

#     #             # Accumulate into the dR histogram
#     #             existing_content = hist_dR.GetBinContent(x_bin)
#     #             hist_dR.SetBinContent(x_bin, existing_content + weighted_content)

#     # # Save the final histogram
#     # output_root_file_path = f"outputs/dm/{dm}/{dm}_{njet}_weighted_dR.root"
#     # output_root_file = ROOT.TFile(output_root_file_path, "RECREATE")
#     # hist_dR.SetName("dR_distribution")
#     # hist_dR.Write()
#     # output_root_file.Close()
