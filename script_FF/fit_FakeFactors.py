import ROOT
import json
import numpy as np

def FitFF(h, func='erf'):
    h_uncert = ROOT.TH1D(h.GetName()+'_uncert', "", 1000, 40, 200)
    
    # Define fit function based on input
    if func == 'erf':
        f2 = ROOT.TF1("f2", "[0]*TMath::Erf((x-[1])/[2])", 40., 200.)
        f2.SetParameter(2, 40)
    elif 'pol' in func:
        f2 = ROOT.TF1("f2", func, 40., 200.)
    else:
        f1 = ROOT.TF1("f1", "landau", 40, 200)
        f2 = ROOT.TF1("f2", "[0]*TMath::Landau(x,[1],[2]) + [3]", 40, 200)

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

def PlotFF(f, h, uncertainty_hist, name, title='', output_folder='./'):
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
    h.GetXaxis().SetTitle('p_{T} (GeV)')
    h.GetYaxis().SetTitle('Fake Factors')
    
    # axis ranges
    h.GetXaxis().SetRangeUser(40, 200)
    # sh.GetYaxis().SetRangeUser(0, 0.4)

    # Redraw axis
    c1.RedrawAxis()
    


    # Save the canvas as PDF
    c1.Print(f"{output_folder}/{name}.png")
    c1.Print(f"{output_folder}/{name}.root")
    c1.Print(f"{output_folder}/{name}.pdf")

def main():

    dms = ["DM0", "DM1", "DM10", "DM11"]  # Loop over all DMs, e.g.
    for dm in dms:
        # Load the JSON file
        with open(f"script_FF/data/FakeFactor_ABCD_{dm}.json", 'r') as json_file:
            data = json.load(json_file)

        # Extract data from JSON
        bin_centers = np.array(data["bin_centers"])
        values = np.array(data["values"])
        yerr = np.array(data["yerr"])
        xerr = np.array(data["xerr"])

        print("bin center = ," , bin_centers)

        rebin = False  # Set to True to merge bins within specified ranges

        if rebin == True:
            # Define new bin edges based on requested merging
            new_bin_edges = np.array([40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 120, 140, 200], dtype=np.float64)  # Non-uniform bin edges

            # Aggregate bin contents and errors within the specified ranges
            new_values = []
            new_yerr = []
            new_bin_centers = []

            # Merge bins within the ranges
            i = 0
            while i < len(bin_centers):
                if 80 <= bin_centers[i] < 90:
                    merged_value = (values[i] + values[i + 1] )/2
                    merged_error = np.sqrt((yerr[i]**2 + yerr[i + 1]**2)/ 2)
                    new_values.append(merged_value)
                    new_yerr.append(merged_error)
                    new_bin_centers.append((bin_centers[i] + bin_centers[i + 1]) / 2)
                    i += 2  # Skip the next bin (since it's merged)
                if 90 <= bin_centers[i] < 100:
                    merged_value = (values[i] + values[i + 1] )/2
                    merged_error = np.sqrt((yerr[i]**2 + yerr[i + 1]**2)/ 2)
                    new_values.append(merged_value)
                    new_yerr.append(merged_error)
                    new_bin_centers.append((bin_centers[i] + bin_centers[i + 1]) / 2)
                    i += 2  # Skip the next bin (since it's merged)
                if 100 <= bin_centers[i] < 120:
                    merged_value = (values[i] + values[i + 1] + values[i + 2] + values[i + 3])/4
                    merged_error = np.sqrt((yerr[i]**2 + yerr[i + 1]**2 + yerr[i + 2]**2 + yerr[i + 3]**2)/4)
                    new_values.append(merged_value)
                    new_yerr.append(merged_error)
                    new_bin_centers.append((bin_centers[i] + bin_centers[i + 1] + bin_centers[i + 2] + bin_centers[i + 3]) / 4)
                    i += 4  # Skip the next bin (since it's merged)
                elif 120 <= bin_centers[i] < 140:
                    merged_value = (values[i] + values[i + 1] + values[i + 2] + values[i + 3])/4 
                    merged_error = np.sqrt((yerr[i]**2 + yerr[i + 1]**2 + yerr[i + 2]**2 + yerr[i + 3]**2)/4)
                    new_values.append(merged_value)
                    new_yerr.append(merged_error)
                    new_bin_centers.append((bin_centers[i] + bin_centers[i + 1] + bin_centers[i + 2] + bin_centers[i + 3]) / 4)
                    i += 4
                elif 140 <= bin_centers[i] < 200:
                    # NB OF BINS TO MERGE
                    n_bins = len(bin_centers[i:])
                    print("n_bins = ," , n_bins)
                    # merged_value = np.mean(values[i:]) # Sum all remaining values
                    remaining_values = [v for v in values[i:] if v != 0]
                    merged_value = np.mean(remaining_values) if remaining_values else 0 
                    print("values = ," , values[i:])
                    merged_error = np.sqrt((np.sum(yerr[i:]**2))/n_bins)  # Sum in quadrature
                    new_values.append(merged_value)
                    new_yerr.append(merged_error)
                    new_bin_centers.append(np.mean(bin_centers[i:]))  # Average bin center
                    print("bin center = ," , np.mean(bin_centers[i:]))
                    print("new bin value = ," , merged_value)
                    break  # End of range
                else:
                    # For bins not being merged, just add them as they are
                    new_values.append(values[i])
                    new_yerr.append(yerr[i])
                    new_bin_centers.append(bin_centers[i])
                    i += 1

            n_bins = len(new_bin_centers)
            print("new bin edges = ," , new_bin_edges)
            h = ROOT.TH1D("h", f"Fake Factor: {data['title']}", n_bins, new_bin_edges)


            # Fill the histogram with merged data
            for j in range(len(new_values)):
                h.SetBinContent(j+1 , new_values[j])
                h.SetBinError(j+1, new_yerr[j]) 
                print("bin", j+1, "value = ," , new_values[j])
                print('bin center = ,' , new_bin_centers[j])
                print('new_bin_edges = ,' , new_bin_edges[j+1])

                print("new_values = ," , new_values)
                print("bin center =  ," , new_bin_centers)
        
        else: 
            print("no rebin")
            # Create a ROOT histogram
            n_bins = len(bin_centers)
            h = ROOT.TH1D("h", f"Fake Factor: {data['title']}", n_bins, bin_centers.min() - xerr[0], bin_centers.max() + xerr[0])
            
            # Fill histogram with data
            for i in range(n_bins):
                h.SetBinContent(i + 1, values[i])
                h.SetBinError(i + 1, yerr[i])

        # Fit the histogram
        fit, h_uncert, h = FitFF(h, func='landau')  # Fit using your chosen function

        # Plot the results
        PlotFF(fit, h, h_uncert, name=data['title'], title='Fake Factor Fit', output_folder='script_FF/plots/fit')


# Call the main function to execute the script
if __name__ == "__main__":
    main()
