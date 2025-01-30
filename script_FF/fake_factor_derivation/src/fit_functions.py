'''
Author : @oponcet
Date : 06-12-2024
Script with different useful function:
    - FitFakeFactors: fit the fake factor with landau or pol1 function
'''

import ROOT
import math

def fit_fake_factor(h, xmin, xmax, usePol1=False, polOnly=None):
    """
    Fit the fake factor with a Landau function and/or polynomial functions.

    Parameters:
    - h: TH1D histogram to fit.
    - usePol1: bool, optional, default=False. If True, include a pol1 term in the fit.
    - polOnly: int, optional, default=None. If specified, fit only with pol0 or pol1 function:
        - 0: Use pol0 function only.
        - 1: Use pol1 function only.
        - 2: Use pol2 function only.
        - 3: Use pol3 function only.
        - None: Use Landau or Landau+pol1 by default.

    Returns:
    - fit: TF1 fit function with the final parameters.
    - h_uncert: TH1D histogram representing the fit uncertainties.
    - h: TH1D histogram, with bins adjusted during fitting.
    """
    if not h or h.GetNbinsX() == 0:
        raise ValueError("Input histogram is invalid or empty.")

    print(f"Fitting histogram {h}")
    
    # Prepare an uncertainty histogram
    h_uncert = ROOT.TH1D(h.GetName() + '_uncert', "", 1000, h.GetBinLowEdge(1), h.GetBinLowEdge(h.GetNbinsX() + 1))
    
    # Define the fit functions
    f1 = ROOT.TF1("f1", "landau", xmin, xmax)
    f2 = ROOT.TF1("f2", "[0]*TMath::Landau(x,[1],[2])+[3]", xmin, xmax)
    if usePol1:
        f2 = ROOT.TF1("f2", "[0]*TMath::Landau(x,[1],[2])+[3]+[4]*x", xmin, xmax)
    
    if polOnly == 0:
        f1 = ROOT.TF1("f1", "pol0", xmin, xmax)
        f2 = ROOT.TF1("f2", "pol0", xmin, xmax)
    elif polOnly == 1:
        f1 = ROOT.TF1("f1", "pol1", xmin, xmax)
        f2 = ROOT.TF1("f2", "pol1", xmin, xmax)
    elif polOnly == 2:
        f1 = ROOT.TF1("f1", "pol2", xmin, xmax)
        f2 = ROOT.TF1("f2", "pol2", xmin, xmax)
    elif polOnly == 3:
        f1 = ROOT.TF1("f1", "pol3", xmin, xmax)
        f2 = ROOT.TF1("f2", "pol3", xmin, xmax)
    elif polOnly == 4:
        f1 = ROOT.TF1("f1", "pol4", xmin, xmax)
        f2 = ROOT.TF1("f2", "pol4", xmin, xmax)
    elif polOnly == -1: # use Gaussian approximation  exp(− ((x−mean)**2) / (2x sigma**2))
        f1 = ROOT.TF1("f1", "[0]*exp(-0.5*((x-[1])/[2])^2)", xmin, xmax)  # Gaussian function without the linear term
        f2 = ROOT.TF1("f2", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x", xmin, xmax)

    
    
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
                print(f"Fit converged after {count} iterations.")
                print(f"Fit result: {fitresult}")
                # Generate confidence intervals for uncertainty
                ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(h_uncert, 0.68)
                fit = f2
                # print(f"Final fit function: {fit.GetExpFormula('P')}")
                # uncerts, fit_up, fit_down = DecomposeUncerts(fitresult, fit)
                fit_up, fit_down, fit_nom= get_variated_fitfunction(fit, fitresult, h_uncert)
                # print(f"Final fit up: {fit_up.GetExpFormula('P')}")
                # print(f"Final fit down: {fit_down.GetExpFormula('P')}")
                # print(f"Final fit nom: {fit_nom.GetExpFormula('P')}")
                break
            count += 1

    else:
        # Direct fit with specified polynomial
        fitresult = h.Fit("f2", 'SIR')
        ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(h_uncert, 0.68)
        fit = f2
        # uncerts, fit_up, fit_down = DecomposeUncerts(fitresult, fit)
        # print(f"Final fit function: {fit.GetExpFormula('P')}")
        fit_up, fit_down, fit_nom = get_variated_fitfunction(fit, fitresult, h_uncert)
        # print(f"Final fit up: {fit_up.GetExpFormula('P')}")
        # print(f"Final fit down: {fit_down.GetExpFormula('P')}")
        # print(f"Final fit nom: {fit_nom.GetExpFormula('P')}")


    # Set range of h_uncert to match fit range
    h_uncert.SetAxisRange(xmin, xmax)
    f2.SetRange(xmin, xmax)

    if fit is None:
        raise RuntimeError("Fit did not converge after 100 iterations.")

    
    # Name and return the fit
    fit.SetName(h.GetName() + '_fit')
    return fit, h_uncert, h, fit_up, fit_down


def fit_correction(h, func='pol1'):
    """
    Fit the correction with a pol1 function.
    """

    # Prepare an uncertainty histogram
    h_uncert = ROOT.TH1D(h.GetName()+'_uncert',"",1000,h.GetBinLowEdge(1),h.GetBinLowEdge(h.GetNbinsX()+1))
    
    # Define the fit functions
    f1 = ROOT.TF1("f1",func)
    #  fit with the full functions
    # repeat fit up to 100 times until the fit converges properly
    rep = True
    count = 0
    while rep:
        fitresult = h.Fit("f1",'SI')
        rep = int(fitresult) != 0
        if not rep or count>100:
            ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(h_uncert, 0.68)
            fit = f1
            break
        count+=1
    fit.SetName(h.GetName()+'_fit')
    return fit, h_uncert

def DecomposeUncerts(fitresult, fit):
    # Decompose the uncertainties of the fit for each fit parameter
    # from https://github.com/danielwinterbottom/TauSF/blob/main/python/fit_tools.py#L4-L43 thanks Danny

    shifted_functions = []
    # Decompose the covariance matrix into eigenvectors
    cov = ROOT.TMatrixD(fitresult.GetCovarianceMatrix())
    eig = ROOT.TMatrixDEigen(cov)
    eigenvectors = eig.GetEigenVectors()

    # Estimate uncertainty variations based on the eigenvectors
    pars = ROOT.TVectorD(fit.GetNpar())
    for i in range(fit.GetNpar()):
        pars[i] = fit.GetParameter(i)
    variances = eig.GetEigenValues()
    transposed_eigenvectors = eigenvectors.Clone().T()

    # Loop over the parameters
    for i in range(fit.GetNpar()):

        temp = ROOT.TVectorD(fit.GetNpar())
        for j in range(fit.GetNpar()):
            temp[j] = transposed_eigenvectors(i, j)


        temp*=variances(i,i)**0.5
        fit_up=fit.Clone()  
        fit_down=fit.Clone() 

        # shift each parameter of the function in turn
        for j in range(fit.GetNpar()): 
            p_uncert = temp[j] 
            nom = fit.GetParameter(j)
            p_up=nom+p_uncert   
            p_down=nom-p_uncert   
            fit_up.SetParameter(j,p_up)
            fit_down.SetParameter(j,p_down)
            fit_up.SetName(fit.GetName()+'_uncert%i_up' %i)
            fit_down.SetName(fit.GetName()+'_uncert%i_down' %i)

            shifted_functions.append(('uncert%i' % i, fit_up, fit_down))
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>> parameter %i nominal is %f" % (i, fit.GetParameter(i)))
        print(f"decompose nominal = {fit.GetExpFormula('P')}")
        print(f"decompose fit up = {fit_up.GetExpFormula('P')}")
        print(f"decompose fit down = {fit_down.GetExpFormula('P')}")


    print(f"Final fit up = {fit_up.GetExpFormula('P')}")
    print(f"Final fit down = {fit_down.GetExpFormula('P')}")

    print(f"Decomposed uncertainties:" , shifted_functions)
    return shifted_functions, fit_up, fit_down

import ROOT

def get_variated_fitfunction(fit, fitresult, h_uncert):
    '''
    Function that returns the fit_up (TF1) and fit_down (TF1) functions by extracting the upper and lower 
    bands from the uncertainty envelope and fitting them with the same fit function.
    '''
    
    # Number of bins in h_uncert should match the number of parameters in the fit result
    fit_up = fit.Clone("fit_up")  # Clone the original fit for the upper bound
    fit_down = fit.Clone("fit_down")  # Clone the original fit for the lower bound
    fit_nom = fit.Clone("fit_nom")  # Clone the original fit for the nominal value
    
    # We need to create histograms for the upper and lower bounds
    # Create new histograms to store the upper and lower bounds
    h_up = h_uncert.Clone("h_up")
    h_down = h_uncert.Clone("h_down")
    
    # Extract the upper and lower bounds from the uncertainty envelope
    for i in range(1, h_uncert.GetNbinsX() + 1):
        # Get the value of the uncertainty band at bin i
        envelop_value = h_uncert.GetBinContent(i)
        envelop_error = h_uncert.GetBinError(i)
        
        # Set the upper and lower bounds of the uncertainty envelope
        h_up.SetBinContent(i, envelop_value + envelop_error)   # Upper bound
        h_down.SetBinContent(i, envelop_value - envelop_error)  # Lower bound

        # Set error of upper and lower bounds to 0
        h_up.SetBinError(i, 0)
        h_down.SetBinError(i, 0)
    
    # Fit the upper and lower bounds with the same fit function (same functional form as the original fit)
    fit_up_result = h_up.Fit(fit_up.GetName(), "S")  # S for 'Save' fit result
    fit_down_result = h_down.Fit(fit_down.GetName(), "S")  # S for 'Save' fit result

    fit_nom_result = h_uncert.Fit(fit.GetName(), "S")  # S for 'Save' fit result

    # write in a root file
    output_file = ROOT.TFile("fit.root", "RECREATE")
    fit_nom.Write()
    fit_up.Write()
    fit_down.Write()

    # plot theme on the same canvas
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)

    h_uncert.Draw("E3") 
    h_uncert.SetFillColorAlpha(ROOT.kAzure + 7, 0.3)

    fit_nom.Draw("SAME")
    fit_up.Draw("SAME")
    fit_down.Draw("SAME")

    h_up.Draw("P SAME")
    h_down.Draw("P SAME")

    canvas.Write()

    canvas.SaveAs("fit.png")

    output_file.Close()
    
    return fit_up, fit_down, fit_nom


    