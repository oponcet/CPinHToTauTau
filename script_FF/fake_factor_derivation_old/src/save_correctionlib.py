'''
Author : @oponcet
Date : 29-01-2025
Script to save the fit results to a CorrectionLib-compatible JSON file with helpers functions.
'''

import os
import json
import ROOT
from correctionlib import schemav2 as cs
import re

def get_variated_function(fit_result):
    """
    Get up and down variation of the fit function using error on the parameters.

    Parameters:
    - fit_result: TF1 fit result.
    """
    fit_up = fit_result.Clone(f"{fit_result.GetName()}_up")
    fit_down = fit_result.Clone(f"{fit_result.GetName()}_down")

    for i in range(fit_result.GetNpar()):
        p = fit_result.GetParameter(i)
        e = fit_result.GetParError(i)
        fit_up.SetParameter(i, p + e)
        fit_down.SetParameter(i, p - e)        

    return fit_up, fit_down


   
def save_to_correctionlib_with_fit(ratio_hist, output_json_file, dm, njet, fit_formula, fit_up_formula, fit_down_formula, correction, variable, x_min=40, x_max=200):
    """
    Converts ROOT histogram with fake factors and a fit function into CorrectionLib JSON format.

    Parameters:
    - ratio_hist: ROOT.TH1D object containing fake factors.
    - output_json_file: Path to save the CorrectionLib-compatible JSON file.
    - dm: Decay mode.
    - njet: Number of jets.
    - fit_formula: Formula string of the fit (TF1.GetExpFormula()).
    - fit_up_formula: Formula string of the fit up variation.
    - fit_down_formula: Formula string of the fit down variation.
    - correction: Type of correction (fake_factor, extrapolation_correction, closure_correction).
    - variable: Variable name (e.g., "met_var_qcd_h1").
    - x_min: Minimum value for the variable (default: 40).
    - x_max: Maximum value for the variable (default: 200).
    - fit_params: List of fit parameters (TF1.GetParameter()).
    """

    fit_formula_converted = convert_fit_formula_to_correctionlib(fit_formula, min_value=x_min, max_value=x_max)
    fit_up_formula_converted = convert_fit_formula_to_correctionlib(fit_up_formula, min_value=x_min, max_value=x_max)
    fit_down_formula_converted = convert_fit_formula_to_correctionlib(fit_down_formula, min_value=x_min, max_value=x_max)

    if correction == "fake_factor":
        main_description = "Fake factors for the httcp analysis"
        name = "fake_factor"
        output_description = "Fake factor to apply to data-MC"
    elif correction == "extrapolation_correction":
        main_description = "Extrapolation correction for the httcp analysis"
        name = "extrapolation_correction"
        output_description = "Extrapolation correction to apply to data-MC"
    else:
        main_description = "Closure correction for the httcp analysis"
        name = "closure_correction"
        output_description = "Closure correction to apply to data-MC"

    njet_int = {"has_0j": 0, "has_1j": 1, "has_2j": 2}.get(njet, -1)
    dm_int = {"pi_1": 0, "rho_1": 1, "a1dm2_1": 2, "a1dm10_1": 10, "a1dm11_1": 11}.get(dm, -1)

    #########################################
    ###      INCLUSIF IN NJETS AND DM     ###
    #########################################
    if njet_int == -1 and dm_int == -1: # inclusif case for both dm and njet
        correctionlib_json = {
        "schema_version": 2,
        "description": main_description,
        "corrections": [
            {
                "name": name,
                "description": main_description,
                "version": 1,
                "inputs": [
                    {
                        "name": variable,
                        "type": "real",
                        "description": "Transverse momentum of the tau"
                    },
                    { 
                        "name": "syst",
                        "type": "string",
                        "description": "Systematic variations: 'nom', 'up', 'down'"
                    }
                ],
                "output": {
                    "name": name,
                    "type": "real",
                    "description": output_description
                },
                "data": {
                        "nodetype": "category",
                        "input": "syst",
                        "content": [
                            {
                                "key": "nom",
                                "value": {
                                    "nodetype": "formula",
                                    "expression": fit_formula_converted,
                                    "parser": "TFormula",
                                    "variables": [variable]
                                }
                            },
                            {
                                "key": "up",
                                "value": {
                                    "nodetype": "formula",
                                    "expression": fit_up_formula_converted,
                                    "parser": "TFormula",
                                    "variables": [variable]
                                }
                            },
                            {
                                "key": "down",
                                "value": {
                                    "nodetype": "formula",
                                    "expression": fit_down_formula_converted,
                                    "parser": "TFormula",
                                    "variables": [variable]
                                }
                            }
                        ]
                    }
            }
        ]
    }
    #########################################
    ###          INCLUSIF IN NJETS        ###
    #########################################
    elif njet_int == -1 and dm_int != -1: # inclusif case for njet
        correctionlib_json = {
            "schema_version": 2,
            "description": main_description,
            "corrections": [
                {
                    "name": name,
                    "description": main_description,
                    "version": 1,
                    "inputs": [
                        {
                            "name": variable,
                            "type": "real",
                            "description": "Transverse momentum of the tau"
                        },
                        {
                            "name": "dm",
                            "type": "int",
                            "description": "Reconstructed tau decay mode of leading tau: 0,1,2,10,11"
                        },
                        { 
                            "name": "syst",
                            "type": "string",
                            "description": "Systematic variations: 'nom', 'up', 'down'"
                        }
                    ],
                    "output": {
                        "name": name,
                        "type": "real",
                        "description": output_description
                    },
                    "data": {
                        "nodetype": "category",
                        "input": "dm",
                        "content": [
                        {
                            "key": dm_int,
                            "value": {
                                    "nodetype": "category",
                                    "input": "syst",
                                    "content": [
                                        {
                                            "key": "nom",
                                            "value": {
                                                "nodetype": "formula",
                                                "expression": fit_formula_converted,
                                                "parser": "TFormula",
                                                "variables": [variable]
                                            }
                                        },
                                        {
                                            "key": "up",
                                            "value": {
                                                "nodetype": "formula",
                                                "expression": fit_up_formula_converted,
                                                "parser": "TFormula",
                                                "variables": [variable]
                                            }
                                        },
                                        {
                                            "key": "down",
                                            "value": {
                                                "nodetype": "formula",
                                                "expression": fit_down_formula_converted,
                                                "parser": "TFormula",
                                                "variables": [variable]
                                            }
                                        }
                                    ]
                                }
                            }
                        ]
                    }
                }
            ]
        }
    #########################################
    ###           INCLUSIF IN DMS         ###
    #########################################
    elif njet_int != -1 and dm_int == -1: # inclusif case for dm
      correctionlib_json = {
        "schema_version": 2,
        "description": main_description,
        "corrections": [
            {
                "name": name,
                "description": main_description,
                "version": 1,
                "inputs": [
                    {
                        "name": variable,
                        "type": "real",
                        "description": "Transverse momentum of the tau"
                    },
                    {
                        "name": "njets",
                        "type": "int",
                        "description": "Number of jets in the event 0,1,2"
                    },
                    { 
                        "name": "syst",
                        "type": "string",
                        "description": "Systematic variations: 'nom', 'up', 'down'"
                    }
                ],
                "output": {
                    "name": name,
                    "type": "real",
                    "description": output_description
                },
                "data": {
                    "nodetype": "category",
                    "input": "njets",
                    "content": [
                        {
                            "key": njet_int,
                            "value": {
                                    "nodetype": "category",
                                    "input": "syst",
                                    "content": [
                                        {
                                            "key": "nom",
                                            "value": {
                                                "nodetype": "formula",
                                                "expression": fit_formula_converted,
                                                "parser": "TFormula",
                                                "variables": [variable]
                                            }
                                        },
                                        {
                                            "key": "up",
                                            "value": {
                                                "nodetype": "formula",
                                                "expression": fit_up_formula_converted,
                                                "parser": "TFormula",
                                                "variables": [variable]
                                            }
                                        },
                                        {
                                            "key": "down",
                                            "value": {
                                                "nodetype": "formula",
                                                "expression": fit_down_formula_converted,
                                                "parser": "TFormula",
                                                "variables": [variable]
                                            }
                                        }
                                    ]
                                }
                        }
                    ]
                }
            }
        ]
    }
    #########################################
    ###    DMS AND NJETS MULTIPLICITY     ###
    #########################################
    else: # inclusif case for both dm and njet
           correctionlib_json = {
        "schema_version": 2,
        "description": main_description,
        "corrections": [
            {
                "name": name,
                "description": main_description,
                "version": 1,
                "inputs": [
                    {
                        "name": variable,
                        "type": "real",
                        "description": "Transverse momentum of the tau"
                    },
                    {
                        "name": "dm",
                        "type": "int",
                        "description": "Reconstructed tau decay mode of leading tau: 0,1,2,10,11"
                    },
                    {
                        "name": "njets",
                        "type": "int",
                        "description": "Number of jets in the event 0,1,2"
                    },
                    { 
                        "name": "syst",
                        "type": "string",
                        "description": "Systematic variations: 'nom', 'up', 'down'"
                    }
                ],
                "output": {
                    "name": name,
                    "type": "real",
                    "description": output_description
                },
                "data": {
                    "nodetype": "category",
                    "input": "dm",
                    "content": [
                    {
                        "key": dm_int,
                        "value": {
                            "nodetype": "category",
                            "input": "njets",
                            "content": [
                                {
                                    "key": njet_int,
                                    "value": {
                                            "nodetype": "category",
                                            "input": "syst",
                                            "content": [
                                                {
                                                    "key": "nom",
                                                    "value": {
                                                        "nodetype": "formula",
                                                        "expression": fit_formula_converted,
                                                        "parser": "TFormula",
                                                        "variables": [variable]
                                                    }
                                                },
                                                {
                                                    "key": "up",
                                                    "value": {
                                                        "nodetype": "formula",
                                                        "expression": fit_up_formula_converted,
                                                        "parser": "TFormula",
                                                        "variables": [variable]
                                                    }
                                                },
                                                {
                                                    "key": "down",
                                                    "value": {
                                                        "nodetype": "formula",
                                                        "expression": fit_down_formula_converted,
                                                        "parser": "TFormula",
                                                        "variables": [variable]
                                                    }
                                                }
                                            ]
                                        }
                                }
                            ]
                        }
                    }]
                }
            }
        ]
    }

    #########################################
    ###           SAVE JSON FILE          ### 
    #########################################
    with open(output_json_file, "w") as f:
        json.dump(correctionlib_json, f, indent=4)
    
    print(f"CorrectionLib JSON with fit formula saved to {output_json_file}")

        

def convert_fit_formula_to_correctionlib(fit_formula, min_value=None, max_value=None):
    """
    Convert a ROOT.TF1 formula to a CorrectionLib-compatible formula.

    Parameters:
        fit_formula (str): Formula string extracted from ROOT.TF1 (e.g., TF1.GetExpFormula()).
        min_value (float): Minimum value for the variable (applied with max()).
        max_value (float): Maximum value for the variable (applied with min()).
    
    Returns:
        str: Converted formula compatible with CorrectionLib.
    """
    # Replace TMath functions with standard operations
    replacements = {
        "TMath::Sq(x)": "x*x",
        "TMath::Sqrt(x)": "sqrt(x)",
        "TMath::Power(x,": "pow(x,",
        "TMath::Exp(x)": "exp(x)",
        "TMath::Log(x)": "log(x)"
    }
    # print("min_value: ", min_value)
    # print("max_value: ", max_value)
    
    print(f"Original fit formula: {fit_formula}")
    for key, value in replacements.items():
        fit_formula = fit_formula.replace(key, value)

    # Add bounds using min and max if specified
    if min_value is not None or max_value is not None:
        var = "x"  # Default variable name
        bounded_variable = var
        if min_value is not None:
            bounded_variable = f"max({bounded_variable},{min_value})"
        if max_value is not None:
            bounded_variable = f"min({bounded_variable},{max_value})"
        # Replace the variable with its bounded version
        fit_formula = fit_formula.replace(var, bounded_variable)

    # Clean up extra spaces (optional, for neatness)
    fit_formula = re.sub(r'\s+', ' ', fit_formula).strip()

    print(f"Converted fit formula: {fit_formula}")
    
    return fit_formula