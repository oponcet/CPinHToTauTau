### This config is used for listing the variables used in the analysis ###

from columnflow.config_util import add_category

import order as od

from columnflow.columnar_util import EMPTY_FLOAT
from columnflow.util import DotDict
from columnflow.columnar_util import ColumnCollection

def keep_columns(cfg: od.Config) -> None:
    # columns to keep after certain steps
    cfg.x.keep_columns = DotDict.wrap({
        "cf.ReduceEvents": {
            # general event info
            "run", "luminosityBlock", "event",
            "PV.npvs","Pileup.nPU","genWeight", "LHEWeight.originalXWGTUP",
            "deterministic_seed", "process_id", "cutflow.*",
        } | {f"MET.{var}" for var in [
            "pt", "phi", "significance",
            "covXX", "covXY", "covYY",
        ]     
         } | {f"Jet.{var}" for var in [
             "pt", "eta", "phi", "mass", 
             "btagDeepFlavB", "hadronFlavour",
         ] 
          } | {f"Tau.{var}" for var in [
              "pt","eta","phi","mass","dxy","dz", "charge", "rawDeepTau2018v2p5VSjet",
              "idDeepTau2018v2p5VSjet", "idDeepTau2018v2p5VSe", "idDeepTau2018v2p5VSmu", 
          ] 
           } | {f"Muon.{var}" for var in [
               "pt","eta","phi","mass","dxy","dz", "charge", 
               "pfRelIso04_all","mT"
           ] 
            } | {
                ColumnCollection.ALL_FROM_SELECTOR
            },
        "cf.MergeSelectionMasks": {
            "normalization_weight", "cutflow.*", "process_id", "category_ids",
        },
        "cf.UniteColumns": {
            "*",
        },
    })


def add_variables(cfg: od.Config) -> None:
    # add variables
    # (the "event", "run" and "lumi" variables are required for some cutflow plotting task,
    # and also correspond to the minimal set of columns that coffea's nano scheme requires)
    cfg.add_variable(
        name="event",
        expression="event",
        binning=(1, 0.0, 1.0e9),
        x_title="Event number",
        discrete_x=True,
    )
    cfg.add_variable(
        name="run",
        expression="run",
        binning=(1, 100000.0, 500000.0),
        x_title="Run number",
        discrete_x=True,
    )
    cfg.add_variable(
        name="lumi",
        expression="luminosityBlock",
        binning=(1, 0.0, 5000.0),
        x_title="Luminosity block",
        discrete_x=True,
    )
    cfg.add_variable(
        name="n_jet",
        expression="n_jet",
        binning=(11, -0.5, 10.5),
        x_title="Number of jets",
        discrete_x=True,
    )
    cfg.add_variable(
        name="jets_pt",
        expression="Jet.pt",
        binning=(40, 0.0, 100.0),
        unit="GeV",
        x_title=r"$p_{T} of all jets$",
    )
    cfg.add_variable(
        name="jets_eta",
        expression="Jet.eta",
        binning=(20, -3., 3.),
        unit="GeV",
        x_title=r"$\eta of all jets$",
    )
    cfg.add_variable(
        name="jet1_pt",
        expression="Jet.pt[:,0]",
        null_value=EMPTY_FLOAT,
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"Jet 1 $p_{T}$",
    )
    cfg.add_variable(
        name="jet1_eta",
        expression="Jet.eta[:,0]",
        null_value=EMPTY_FLOAT,
        binning=(30, -3.0, 3.0),
        x_title=r"Jet 1 $\eta$",
    )
    
    cfg.add_variable(
        name="jet_raw_DeepJetFlavB",
        expression="Jet.btagDeepFlavB",
        null_value=EMPTY_FLOAT,
        binning=(30, 0,1),
        x_title=r"raw DeepJetFlawB",
    )
    
    
    # cfg.add_variable(
    #     name="ht",
    #     expression=lambda events: ak.sum(events.Jet.pt, axis=1),
    #     binning=(40, 0.0, 800.0),
    #     unit="GeV",
    #     x_title="HT",
    # )
    # weights
    cfg.add_variable(
        name="mc_weight",
        expression="mc_weight",
        binning=(200, -100, 100),
        x_title="MC weight",
    )
    # cutflow variables
    cfg.add_variable(
        name="cf_jet1_pt",
        expression="cutflow.jet1_pt",
        binning=(40, 0.0, 400.0),
        unit="GeV",
        x_title=r"Jet 1 $p_{T}$",
    )
    
    cfg.add_variable(
        name="muon_eta",
        expression="Muon.eta",
        null_value=EMPTY_FLOAT,
        binning=(30, -3.0, 3.0),
        x_title=r"muon $\eta$",
    )
    
    cfg.add_variable(
        name="muon_pt",
        expression="Muon.pt",
        null_value=EMPTY_FLOAT,
        binning=(30, 20.0, 80.0),
        unit="GeV",
        x_title=r"muon $p_{T}$",
    )
    cfg.add_variable(
        name="muon_phi",
        expression="Muon.phi",
        null_value=EMPTY_FLOAT,
        binning=(20, -3.14159, 3.14159),
        x_title=r"muon $\varphi$",
    )
    
    cfg.add_variable(
        name="tau_phi",
        expression="Tau.phi",
        null_value=EMPTY_FLOAT,
        binning=(20, -3.14159, 3.14159),
        x_title=r"tau $\varphi$",
    )
    
    cfg.add_variable(
        name="tau_eta",
        expression="Tau.eta",
        null_value=EMPTY_FLOAT,
        binning=(30, -3.0, 3.0),
        x_title=r"tau $\eta$",
    )
    
    cfg.add_variable(
        name="tau_pt",
        expression="Tau.pt",
        null_value=EMPTY_FLOAT,
        binning=(30, 25, 85),
        unit="GeV",
        x_title=r"tau $p_{T}$",
    )
    
    cfg.add_variable(
        name="muon_mT",
        expression="Muon.mT",
        null_value=EMPTY_FLOAT,
        binning=(40, 0.0, 200.0),
        unit="GeV",
        x_title=r"leading muon $m_{T}$",
    )
    
    cfg.add_variable(
        name="mutau_mass",
        expression="mutau_mass",
        null_value=EMPTY_FLOAT,
        binning=(40, 0.0, 200.0),
        unit="GeV",
        x_title=r"$m_{vis}$",
    )
    
    cfg.add_variable(
        name="pu_weight",
        expression="pu_weight",
        null_value=EMPTY_FLOAT,
        binning=(30, 0,3),
        unit="GeV",
        x_title=r"Pileup weight",
    )
