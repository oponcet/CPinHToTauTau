'''
Author : @oponcet
Date : 06-12-2024
Script with different useful function:
    - FitFakeFactors: fit the fake factor with landau or pol1 function
'''

import ROOT

def fit_fake_factor(h, usePol1=False, polOnly=None):
    """
    Fit the fake factor with a Landau function and/or polynomial functions.

    Parameters:
    - h: TH1D histogram to fit.
    - usePol1: bool, optional, default=False. If True, include a pol1 term in the fit.
    - polOnly: int, optional, default=None. If specified, fit only with pol0 or pol1 function:
        - 0: Use pol0 function only.
        - 1: Use pol1 function only.
        - None: Use Landau or Landau+pol1 by default.

    Returns:
    - fit: TF1 fit function with the final parameters.
    - h_uncert: TH1D histogram representing the fit uncertainties.
    - h: TH1D histogram, with bins adjusted during fitting.
    """
    if not h or h.GetNbinsX() == 0:
        raise ValueError("Input histogram is invalid or empty.")
    
    # Prepare an uncertainty histogram
    h_uncert = ROOT.TH1D(h.GetName() + '_uncert', "", 1000, h.GetBinLowEdge(1), h.GetBinLowEdge(h.GetNbinsX() + 1))
    
    # Define the fit functions
    f1 = ROOT.TF1("f1", "landau", 40, 200)
    f2 = ROOT.TF1("f2", "[0]*TMath::Landau(x,[1],[2])+[3]", 40, 200)
    if usePol1:
        f2 = ROOT.TF1("f2", "[0]*TMath::Landau(x,[1],[2])+[3]+[4]*x", 40, 200)
    
    if polOnly == 0:
        f1 = ROOT.TF1("f1", "pol0", 20, 200)
        f2 = ROOT.TF1("f2", "pol0", 20, 200)
    elif polOnly == 1:
        f1 = ROOT.TF1("f1", "pol1", 20, 200)
        f2 = ROOT.TF1("f2", "pol1", 20, 200)
    
    # Reset histogram, keeping only bins with content > 0
    h_clone = h.Clone()
    h_clone.Reset()
    for i in range(1, h.GetNbinsX() + 1):
        content = h.GetBinContent(i)
        error = h.GetBinError(i)
        if content > 0:
            h_clone.SetBinContent(i, content)
            h_clone.SetBinError(i, error)
    h = h_clone

    # Fit logic
    fit = None
    if polOnly is None:
        # Initial Landau fit for parameter seeding
        h.Fit("f1", 'IR')
        f2.SetParameter(0, f1.GetParameter(0))
        f2.SetParameter(1, f1.GetParameter(1))
        f2.SetParameter(2, f1.GetParameter(2))
        f2.SetParameter(3, 0)  # Initial offset
        if usePol1:
            f2.SetParameter(4, 0)  # Initial slope

        # Iterative fitting
        rep = True
        count = 0
        while rep:
            fitresult = h.Fit("f2", 'SIR')
            rep = int(fitresult) != 0  # Repeat if fit fails
            if not rep or count > 100:
                # Generate confidence intervals for uncertainty
                ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(h_uncert, 0.68)
                fit = f2
                break
            count += 1

    else:
        # Direct fit with specified polynomial
        fitresult = h.Fit("f2", 'SIR')
        ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(h_uncert, 0.68)
        fit = f2

    if fit is None:
        raise RuntimeError("Fit did not converge after 100 iterations.")

    # Get fit statistics and parameters
    chi2 = fit.GetChisquare()
    ndf = fit.GetNDF()

    # Get parameter values and errors
    parameters = []
    for i in range(fit.GetNpar()):
        parameters.append({
            "p": fit.GetParameter(i),
            "error": fit.GetParError(i)
        })

    # Fill the fit details dictionary
    fit_details = {
        "Chi2": chi2,
        "NDf": ndf,
        "Parameters": parameters
    }

    # Name and return the fit
    fit.SetName(h.GetName() + '_fit')
    return fit, h_uncert, h, fit_details
