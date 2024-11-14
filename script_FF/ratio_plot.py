import ROOT

ROOT.gROOT.SetBatch(True)

def ratio_plot(file1_path, file2_path, hist_name):
    # Open the ROOT files
    file1 = ROOT.TFile.Open(file1_path)
    file2 = ROOT.TFile.Open(file2_path)
    
    # Get the histograms
    hist1 = file1.Get(hist_name)
    hist2 = file2.Get(hist_name)
    
    # Check if histograms are successfully loaded
    if not hist1 or not hist2:
        print(f"Failed to load histograms {hist_name}")
        return
    
    # Create a canvas
    c = ROOT.TCanvas("c", "Ratio Plot", 800, 600)
    
    # Create the ratio histogram
    ratio_hist = hist1.Clone("ratio_hist")
    ratio_hist.Divide(hist2)
    
    # Draw the histograms and the ratio plot
    c.Divide(1, 2)
    
    c.cd(1)
    hist1.SetLineColor(ROOT.kRed)
    hist2.SetLineColor(ROOT.kBlue)
    hist2.Draw("HIST")
    hist1.Draw("HIST SAME")
    
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    legend.AddEntry(hist1, "SS Iso", "l")
    legend.AddEntry(hist2, "SS anti Iso", "l")
    legend.Draw()
    
    c.cd(2)
    ratio_hist.SetLineColor(ROOT.kBlack)
    ratio_hist.SetTitle("Ratio Plot of SS Iso / SS anti Iso")
    ratio_hist.GetYaxis().SetTitle("Ratio")
    ratio_hist.Draw("P")
    ratio_hist.SetMarkerStyle(2)

    
    # Save the canvas
    c.SaveAs("ratio_plot.root")
    
    # Close the files
    file1.Close()
    file2.Close()

# Example usage
file1_path = "/afs/cern.ch/user/o/oponcet/private/analysis/CPinHToTauTau/histogram_cat509_tau_1_pt.root"
file2_path = "/afs/cern.ch/user/o/oponcet/private/analysis/CPinHToTauTau/histogram_cat510_tauantiiso_1_pt.root"
hist_name = "hist_name"

ratio_plot(file1_path, file2_path, hist_name)

