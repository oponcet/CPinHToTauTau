import ROOT
import json
import numpy as np

def FitFF(h, minh, maxh, func='erf',):
    h_uncert = ROOT.TH1D(h.GetName()+'_uncert', "", 1000, minh, maxh)
    
    # Define fit function based on input
    if func == 'erf':
        f2 = ROOT.TF1("f2", "[0]*TMath::Erf((x-[1])/[2])", minh, maxh)
        f2.SetParameter(2, 40)
    elif 'pol' in func:
        f2 = ROOT.TF1("f2", func, minh, maxh)
    else:
        f1 = ROOT.TF1("f1", "landau", minh, maxh)
        f2 = ROOT.TF1("f2", "[0]*TMath::Landau(x,[1],[2]) + [3]", minh, maxh)

    # Fit first with Landau to get initial values for parameters
    if func == 'landau':
        h.Fit("f1", 'IR')
        f2.SetParameter(0, f1.GetParameter(0))
        f2.SetParameter(1, f1.GetParameter(1))
        f2.SetParameter(2, f1.GetParameter(2))
        f2.SetParameter(3, 0)

    # Fit with the full function, retrying if necessary
    rep = True
    count = 0
    maxN = 100
    while rep:
        fitresult = h.Fit("f2", 'SIR')
        rep = int(fitresult) != 0
        if not rep or count > maxN:
            ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(h_uncert, 0.68)
            fit = f2
            break
        count += 1

    # fit.SetName(h.GetName() + '_fit')

    print(f'Chi2/NDF = {f2.GetChisquare():.2f}/{f2.GetNDF():.0f}, p-value = {f2.GetProb():.2f}')
    return fit, h_uncert, h

def PlotFF(f, h, uncertainty_hist, name, title='', output_folder='./', var_title='p_{T} (GeV)'):
    # Create a canvas
    c1 = ROOT.TCanvas("c1", "Canvas", 800, 600)

    # Set axis titles
    f.GetXaxis().SetTitle('p_{T} (GeV)')
    f.GetYaxis().SetTitle('Fake Factor Correction')
    f.SetTitle(title)
    f.SetLineColor(ROOT.kAzure + 7)
    f.SetLineWidth(2)

    # Draw the histogram with error bars
    h.SetStats(0)  # Disable statistics box
    h.SetLineColor(ROOT.kBlack)  # Set fill color
    h.Draw("EP")  # Draw histogram with error bars

    # Draw the fit function
    f.Draw("same")  # Draw fit function

    # Draw the uncertainty histogram (if provided)
    if uncertainty_hist:
        #uncertainty_hist.SetFillColor(ROOT.kBlue - 10)  # Set fill color for uncertainty
        uncertainty_hist.SetFillColorAlpha(ROOT.kAzure + 7 , 0.4)
        uncertainty_hist.Draw("E3 same")  # Draw uncertainty as a filled histogram

   

    # Set title axis
    h.GetXaxis().SetTitle(var_title)
    h.GetYaxis().SetTitle('Correction Factors')
    
    # axis ranges
    # h.GetXaxis().SetRangeUser(40, 200)
    # sh.GetYaxis().SetRangeUser(0, 0.4)

    # Redraw axis
    c1.RedrawAxis()
    


    # Save the canvas as PDF
    c1.Print(f"{output_folder}/{name}.png")
    c1.Print(f"{output_folder}/{name}.root")
    c1.Print(f"{output_folder}/{name}.pdf")

def main():


    dms = ["DM0", "DM1", "DM10", "DM11"]  # Loop over all DMs, e.g.
    variables = ["puppi_met_phi", "puppi_met_pt", "hcand_dr", "hcand_invm"]
    variables_title = ["MET #phi", "MET p_{T}", "#Delta R(h, cand)", "m_{h, cand}"]
    for var in variables:
        for dm in dms:
            # Load the JSON file
            with open(f"script_FF/data/control_plots/{var}_DataMC_Ratio_A_{dm}.json", 'r') as json_file:
                data = json.load(json_file)

            # Extract data from JSON
            bin_centers = np.array(data["bin_centers"])
            values = np.array(data["ratio_values"])
            yerr = np.array(data["ratio_errors"])

            print("bin center = ," , bin_centers)

            
            # Create a ROOT histogram
            n_bins = len(bin_centers)
            h = ROOT.TH1D("h", f"Fake Factor: {data['title']}", n_bins, bin_centers.min(), bin_centers.max())
            
            # Fill histogram with data
            for i in range(n_bins):
                h.SetBinContent(i + 1, values[i])
                h.SetBinError(i + 1, yerr[i])

            minh = bin_centers.min()
            maxh = bin_centers.max()
            # Fit the histogram
            fit, h_uncert, h = FitFF(h,minh, maxh,  func='landau')  # Fit using your chosen function

            # # Plot the results
            PlotFF(fit, h, h_uncert, name=f"{var}_{data['title']}", title='Correction Fit', output_folder='script_FF/plots/fit/control_plots', var_title=variables_title[variables.index(var)])


# Call the main function to execute the script
if __name__ == "__main__":
    main()
