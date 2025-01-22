'''
Author : @oponcet
Date : 06-12-2024
Script to apply pt-depedent Fake Factors to distribution of any other variable 
(met_var_qcd_h1 for instance). 
Apply FF to B data_min_mc and use it in A region as Fake jets -> tau_h 


'''
import os
import json
import ROOT
from plotting import *
import cmsstyle as CMS
from fit_functions import fit_fake_factor
from derive_FF import save_to_correctionlib_with_fit

def remove_TPaveText(canvas):
    # Remove any existing TPave elements
    for primitive in canvas.GetListOfPrimitives():
        if isinstance(primitive, ROOT.TPaveText):  # Check if it's a TPaveText (or similar)
            primitive.Delete()  # Remove the unwanted TPaveText

def apply_fake_factor_CD(input_file, catC, catD, dm, njet):
    """
    Apply FF to C data_min_mc and use it in D region as Fake jets -> tau_h and produce
    the stack plot of the variable.

    Parameters:
    - input_file: Path to input ROOT file for region A.
    - catC: category C (numerator) -> TH1D is in catC/variable/data_minus_mc
    - catD: category D (denominator) -> TH1D is in catD/variable/data_minus_mc
    - dm: decay mode
    - njet: number of jets
    """

    output_root_file_path = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/dm/{dm}/{dm}_{njet}_hcand_1_pt_D.root'
    output_png_file = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/dm/{dm}/{dm}_{njet}_hcand_1_pt_D'

    # Ensure directories exist if not create them
    os.makedirs(os.path.dirname(output_root_file_path), exist_ok=True)

    #####################################################################
    ### Load input ROOT files and histograms that contains all histos ###
    #####################################################################

    # Load input ROOT files and histograms that contains all histos
    root_file = ROOT.TFile.Open(input_file, "READ")

    hist_c = root_file.Get(f'{catC}/hcand_1_pt/data_minus_mc')

    if not hist_c :
        raise KeyError(f"Histogram '{catC}/hcand_1_pt/data_minus_mc not found in one of the input files.")

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
    hist_c_fake = hist_c.Clone()

    # change range of hist_c_fake to 40-200

    for i in range(hist_c_fake.GetNbinsX()):
        pt = hist_c_fake.GetXaxis().GetBinCenter(i)
        fake_factor = FF_function.Eval(pt)
        hist_c_fake.SetBinContent(i, hist_c_fake.GetBinContent(i)*fake_factor)

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

    # Detach the histograms from the ROOT file
    data_hist.SetDirectory(0)
    hist_c_fake.SetDirectory(0)

    # Close canvas and create a new one
    canvas_d_stack.Close()
    canvas_d_stack = CMS.cmsCanvas('Data-MC', 0, 1, 0, 1, '', '', square = CMS.kSquare, extraSpace=0.01, iPos=0)


    # rename the fake jets histogram
    hist_c_fake.SetName("fakes_jets")

    # Add the fake jets to the stack 
    mc_stack.Add(hist_c_fake)

    # Draw data
    data_hist.Draw("E ")
    mc_stack.Draw("HIST SAME")
    data_hist.Draw("E SAME")

    data_hist.SetTitle("")

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

    # # Save the output ROOT file
    # output_root_file.Close()

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
    ratio.GetYaxis().SetRangeUser(0.5, 2.5)  # Adjust range as needed
    # automatci range 
    # ratio.GetYaxis().SetRangeUser(ratio.GetMinimum()-0.1, ratio.GetMaximum()+0.1)
    ratio.GetYaxis().SetNdivisions(3)
    ratio.GetXaxis().SetTitle("p_{T} (GeV)")    # Replace with appropriate axis label
    ratio.SetMarkerStyle(20)

    print("creating canvas")
    # Create a new canvas and draw the plots
    output_canvas_ratio = CMS.cmsDiCanvas('Data-MC Ratio', 40, 200, 0, 1, 0, 2.5, 'p_{T} (GeV)', 'Events/5 GeV', 'Data/MC', square = CMS.kSquare, iPos=0, extraSpace=0.01)



    # Upper pad for the histograms
    upper_pad = output_canvas_ratio.cd(1)
    upper_pad.SetPad(0, 0.3, 1, 1) # (xlow, ylow, xup, yup)
    upper_pad.SetBottomMargin(0.02)
    data_hist.SetMarkerStyle(20)

    data_hist.Draw("E ")
    mc_stack.Draw("HIST SAME")
    data_hist.Draw("E SAME")

    # data_hist.SetTitle("")
    title = f"Fake Factor for {dm} {njet} jets"
    data_hist.SetTitle(title) 
    


    # remove x axis title and labels
    # data_hist.GetXaxis().SetTitle("")
    data_hist.GetXaxis().SetLabelSize(0)
    data_hist.GetXaxis().SetRangeUser(40, 200)
    data_hist.GetYaxis().SetTitle("Events/5 GeV")

    output_canvas_ratio.Update()
    output_canvas_ratio.Modified()


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
    ratio.GetYaxis().SetTitleOffset(0.5)  # Adjust this value to move the title to the right
    ratio.GetYaxis().SetTitle("Data/MC")
    ratio.GetYaxis().SetRangeUser(0.0, 2.5)
    ratio.GetYaxis().SetNdivisions(3)
    ratio.SetMarkerStyle(20)


    ratio.Draw("E")
    ratio.SetTitle("")

    # Draw dash line at y=1
    line = ROOT.TLine(40, 1, 200, 1)
    line.SetLineStyle(2)
    line.Draw()

    output_canvas_ratio.Update()
    output_canvas_ratio.Modified()

    # Save the canvas to the output file
    output_canvas_ratio.Write()

    # Save the canvas as a PNG file
    output_canvas_ratio.SaveAs(output_png_file_ratio)
    

    # detatch the histograms from the ROOT file
    ratio.SetDirectory(0)

    output_root_file.Close()




def apply_fake_factor_AB(input_file, catA, catB, dm, njet,var2D):
    """
    Apply FF to B data_min_mc and use it in A region as Fake jets -> tau_h and produce
    closure control plots.
    Parameters:
    - input_file: Path to input ROOT file for region A. 
    - catA: category A (numerator) -> TH1D is in catA/variable/data_minus_mc
    - catB: category B (denominator) -> TH1D is in catB/variable/data_minus_mc
    - dm: decay mode
    - njet: number of jets
    - var2D: variable to apply the fake factor

    """

    var1 = var2D[0]
    var2 = var2D[1]

    # Define output paths
    output_root_file_path = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/dm/{dm}/{dm}_{njet}_{var1}_A.root'
    output_png_file = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/dm/{dm}/{dm}_{njet}_{var1}_A'

    # Ensure directories exist if not create them
    os.makedirs(os.path.dirname(output_root_file_path), exist_ok=True)
    os.makedirs(os.path.dirname(output_png_file), exist_ok=True)


    hist_a_fake, hist_a_fake_proj = apply_fake_factor_var(input_file, catB, dm, njet,var2D)

    # print("hist_a_fake_proj",hist_a_fake_proj)
    # print("hist_a_fake",hist_a_fake)

    ##########################################################
    ### Build hist stack and data_minus_mc_wfake histogram ###
    ##########################################################

    # Load input ROOT files and histograms
    root_file = ROOT.TFile.Open(input_file, "READ")

    # Get the stack canvas
    canvas_a_stack = root_file.Get(f'{catA}/{var1}/canvas').Clone()
    root_file.Close()

    # Open the output ROOT file
    output_root_file = ROOT.TFile(output_root_file_path, "RECREATE")
    output_root_file.cd()

    # Get THStack mc_stacks in canvas and data
    data_hist = canvas_a_stack.GetPrimitive("data_hist")  # Ensure this matches your histogram name
    mc_stack = canvas_a_stack.GetPrimitive('mc_stack')

    # Detach the histograms from the ROOT file
    data_hist.SetDirectory(0)

    # Close canvas and create a new one
    canvas_a_stack.Close()
    canvas_a_stack = CMS.cmsCanvas('Data-MC', 0, 1, 0, 1, '', '', square = CMS.kSquare, extraSpace=0.01, iPos=0)

    # Add the fake jets to the stack
    mc_stack.Add(hist_a_fake_proj)

    # Draw 
    data_hist.Draw("E1")
    mc_stack.Draw("HIST SAME")
    data_hist.Draw("E1 SAME")


    ### STYLE THE PLOT ###
    # Add a legend
    legend = CMS.cmsLeg(0.7, 0.7, 0.9, 0.9, textSize=0.02)
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
    data_hist.GetYaxis().SetTitle("Events/Bin GeV")
    data_hist.GetXaxis().SetTitle(var1)


    legend.SetTextSize(0.02)
    legend.Draw()

    # Save the canvas to the output file
    canvas_a_stack.Write() 


    #######################################################################
    ###   Create stack plot of A with hist_B_fake as fake with ratio    ###
    #######################################################################
    # Compute the total MC histogram
    mc_total = mc_stack.GetStack().Last().Clone()  # Get the total MC histogram from the stack

    output_png_file_ratio = f'{output_png_file}_ratio.pdf'
    
    
    # Create a ratio plot
    ratio = data_hist.Clone("ratio")
    ratio.Divide(mc_total)


    print("creating canvas")
    # Create a new canvas and draw the plots
    output_canvas_ratio = CMS.cmsDiCanvas('Data-MC Ratio', 0, 1, 0, 1, 0, 2.5, 'p_{T} (GeV)', 'Events/5 GeV', 'Data/MC', square = CMS.kSquare, iPos=0, extraSpace=0.01)


    # Upper pad for the histograms
    upper_pad = output_canvas_ratio.cd(1)
    upper_pad.SetPad(0, 0.3, 1, 1) # (xlow, ylow, xup, yup)
    upper_pad.SetBottomMargin(0.02)
    data_hist.SetMarkerStyle(20)
    
    data_hist.Draw("E")
    mc_stack.Draw("HIST SAME")
    data_hist.Draw("E SAME")

    # remove x axis title and labels
    data_hist.GetXaxis().SetTitle("")
    data_hist.GetXaxis().SetLabelSize(0)
    data_hist.GetYaxis().SetTitle("Events/Bin GeV")
    minx = data_hist.GetXaxis().GetXmin()
    maxx = data_hist.GetXaxis().GetXmax()

    # Set title
    data_hist.SetTitle("")
    

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
    ratio.GetYaxis().SetTitleOffset(0.5)  # Adjust this value to move the title to the right
    ratio.GetYaxis().SetTitle("Data/MC")
    ratio.GetYaxis().SetRangeUser(0.5, 1.5)
    ratio.GetYaxis().SetNdivisions(3)
    ratio.SetMarkerStyle(20)
    ratio.Draw("E")
    ratio.SetTitle("")

    # Draw dash line at y=1
    line = ROOT.TLine(minx, 1, maxx, 1) # (x1, y1, x2, y2)
    line.SetLineStyle(2)
    line.Draw("SAME")

    output_canvas_ratio.Update()


    # Save the canvas to the output file
    output_canvas_ratio.Write()

    # Save the canvas as a PNG file
    output_canvas_ratio.SaveAs(output_png_file_ratio)

    # detatch the histograms from the ROOT file
    ratio.SetDirectory(0)

    output_root_file.Close()
    ### derive the closure correction for met_var_qcd_h1_hcand_1_pt

    print(f"var1: {var1}")

    if var1 == "met_var_qcd_h1":
        print(f"ratio: {ratio}")
        derive_CCorr(ratio, "met_var_qcd_h1", dm, njet)


def derive_CCorr(ratio_hist, variable, dm, njet):

    # ROOT.gROOT.SetStyle("Plain")

    # output file 
    output_root_file_path = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/closure/dm/{dm}/{dm}_{njet}_{variable}.root'

    # ensure directories exist if not create them
    os.makedirs(os.path.dirname(output_root_file_path), exist_ok=True)

    # Open the output ROOT file
    output_file = ROOT.TFile(output_root_file_path, "RECREATE")

    # Also save the ratio in a TH1D
    ratio_hist.Write("ratio_hist")

    # Divide the number on bin by two by merging the bins
    ratio_hist.Rebin(2)
    # Divive each bin by 2 to get the average
    ratio_hist.Scale(0.5)
    
    

    fit_result, h_uncert, ratio_hist, fit_details = fit_fake_factor(ratio_hist, -1.5, 1.5, polOnly=3, usePol1=False)

    # Save the fake factor and fit results to a ROOT file

    ratio_hist.Write("ClosureCorrection")
    fit_result.Write("FitResult")
    h_uncert.Write("Uncertainties")

    #### CANVAS FULL PLOT #### 

    # Plot the uncertainties as a filled area
    CMS.SetCmsText("")
    CMS.SetExtraText("Private Work")
    CMS.SetLumi("")
    CMS.SetEnergy(13.6, unit='TeV')
    uncert_canvas = CMS.cmsCanvas('Closure Correction', -1.5, 1.4, 0, 2.5, 'MET_var_QCD', 'Closure Correction', square = CMS.kSquare, extraSpace=0.01, iPos=0)
   
    # Draw the uncertainties as a filled area
    # Set label size for axes
    h_uncert.GetXaxis().SetLabelSize(0.02)
    h_uncert.GetYaxis().SetLabelSize(0.02)
    CMS.cmsDraw(h_uncert, "E3", fcolor = ROOT.TColor.GetColor("#85D1FBff"), alpha = 0.8, lwidth = 0, lcolor = ROOT.kAzure + 7)

    # Draw the fake factor histogram
    CMS.cmsDraw(ratio_hist, "Same P")
    ratio_hist.SetLineColor(ROOT.kBlack)
    ratio_hist.SetMarkerStyle(20)
    ratio_hist.GetYaxis().SetRangeUser(0, 2.5)

    # Add a legend for clarity
    legend = ROOT.TLegend(0.15, 0.75, 0.35, 0.9) #  
    legend.AddEntry(ratio_hist, "Fake Factor", "EP")
    legend.AddEntry(fit_result, "Fit Result", "L")
    legend.AddEntry(h_uncert, "95% CL (Uncertainties)", "F")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)
    legend.Draw()


    # remove stat box and fit box
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)



    # Save the canvas to an image file
    output_image_path = output_root_file_path.replace(".root", "_FullPlot.pdf")
    uncert_canvas.SaveAs(output_image_path)

    # #### PNG FILE: FIT DETAIL #### uncert_canvas
    # # Create a text box to show the fit details
    fit_details_text = ROOT.TPaveText(0.35, 0.69, 0.6, 0.89, "NDC") # x1, y1, x2, y2, option
    fit_details_text.SetBorderSize(0)
    fit_details_text.SetFillColor(0)
    fit_details_text.SetTextAlign(12)
    fit_details_text.SetTextSize(0.02)

    # # Add fit statistics and parameters to the text box
    fit_details_text.AddText(f"Chi2 = {fit_details['Chi2']:.4f}")
    fit_details_text.AddText(f"NDf = {fit_details['NDf']}")

    # # Add parameter values with their errors
    # for i, param in enumerate(fit_details['Parameters']):
    #     fit_details_text.AddText(f"p{i} = {param['p']:.4f} \pm {param['error']:.4f}")

    # # Draw the text box
    fit_details_text.Draw()

    # Add fit details to the canvas
    ROOT.gStyle.SetOptFit(1)


    # Save the fit details canvas as a separate PNG file
    output_details_image_path = output_root_file_path.replace(".root", "_FitDetail.pdf")
    uncert_canvas.SaveAs(output_details_image_path)

    # Save the canvas to the ROOT file for later use
    uncert_canvas.Write("Fullplot")

    # detacth the histograms from the ROOT file
    ratio_hist.SetDirectory(0)

    output_file.Close()

    #### JSON FILE ####
    output_json_file = f'script_FF/fake_factor_derivation/outputs/ClosureCorrection/{dm}/{dm}_{njet}.json'


    # Save the fake factor to a JSON file
    os.makedirs(os.path.dirname(output_json_file), exist_ok=True)

    fit_formula = str(fit_result.GetExpFormula("P"))  # Explicitly cast to a Python string

    save_to_correctionlib_with_fit(ratio_hist, fit_result, output_json_file, dm, njet, fit_formula, [fit_result.GetParameter(i) for i in range(fit_result.GetNpar())], correction_name="closure_correction", variable_name="met_var_qcd_h1", pt_min= -1.5, pt_max=1.5)



def apply_fake_factor_var(input_file, cat, dm, njet,var2D):
    """
    Apply FF to any variable Y, use the 2D plot of (Y,pt) to apply the fake factor.
    Return the 2D histogram and 1D with the fake factor applied of data_minus_mc.
    This function is used to derive FF0 of fakes. 
    Parameters:
    - input_file: Path to input ROOT file for region.
    - cat: category 
    - dm: decay mode
    - njet: number of jets
    - var2D: variable to apply the fake factor

    """

    var1 = var2D[0]
    var2 = var2D[1]

    if "hadA" in cat:
        output_root_file_path = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/dm/{dm}/{dm}_{njet}_{var1}_A.root'
    if "hadB" in cat:
        output_root_file_path = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/dm/{dm}/{dm}_{njet}_{var1}_B.root'
    if "hadC" in cat:
        output_root_file_path = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/dm/{dm}/{dm}_{njet}_{var1}_C.root'
    if "hadD" in cat:
        output_root_file_path = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/dm/{dm}/{dm}_{njet}_{var1}_D.root'
    else:
        output_root_file_path = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/dm/{dm}/{dm}_{njet}_{var1}_.root'

    # Define output paths
    # output_root_file_path = f'script_FF/fake_factor_derivation/outputs/outputs_applyFF/dm/{dm}/{cat}/{dm}_{njet}_{var1}_.root'

    # Ensure directories exist if not create them
    os.makedirs(os.path.dirname(output_root_file_path), exist_ok=True)
    

    # Load input ROOT files and histograms
    root_file = ROOT.TFile.Open(input_file, "READ")
    
    if not root_file or root_file.IsZombie():
        raise FileNotFoundError("Input ROOT files could not be opened.")

    # Get the 2D histograms
    hist_2D = root_file.Get(f'{cat}/{var1}-{var2}/data_minus_mc').Clone()

    hist_2D.SetDirectory(0) # Detach histogram from ROOT file to avoid deletion

    if not hist_2D:
        raise KeyError(f"Histogram '{catA}/{var1}-{var2}/data_minus_mc not found in one of the input files.")

    # close the input ROOT file
    root_file.Close()

    # Get the fake factor
    fake_factor_file = f'script_FF/fake_factor_derivation/outputs/FakeFactors/{dm}/{dm}_{njet}.root'

    # Load input ROOT files and histograms
    root_fake_factor_file = ROOT.TFile.Open(fake_factor_file, "READ")

    # Get the TF1 fake factor from the fit
    FF_function = root_fake_factor_file.Get(f'FitResult')

    # Apply FF_function to hist_2D
    hist_2D_fake = hist_2D.Clone()

    for i in range(hist_2D_fake.GetNbinsX()):
        var = hist_2D_fake.GetXaxis().GetBinCenter(i) # X axis is var
        for j in range(hist_2D_fake.GetNbinsY()):
            pt = hist_2D_fake.GetYaxis().GetBinCenter(j) # Y axis is pt
            fake_factor = FF_function.Eval(pt)
            hist_2D_fake.SetBinContent(i,j,hist_2D_fake.GetBinContent(i,j)*fake_factor)

    ## Project the 2D histogram on the VAR axis
    hist_2D_fake_proj = hist_2D_fake.ProjectionX() 

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

    return hist_2D_fake, hist_2D_fake_proj

    


    









