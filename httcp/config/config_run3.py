# coding: utf-8

"""
Configuration of the higgs_cp analysis.
"""

import functools
import itertools

import os
import re
import law
import yaml
from glob import glob
import order as od
from scinum import Number

from columnflow.util import DotDict, maybe_import, dev_sandbox
from columnflow.config_util import (
    get_root_processes_from_campaign, 
    add_category,
    add_shift_aliases, 
    get_shifts_from_sources,
    verify_config_processes,
)

ak = maybe_import("awkward")

#thisdir = os.path.dirname(os.path.abspath(__file__))
thisdir = "/afs/cern.ch/work/g/gsaha/public/IPHC/Work/ColumnFlowAnalyses/CPinHToTauTau/httcp/config"
print(f"thisdir: {thisdir}")
corrdir = os.path.join(os.path.dirname(thisdir), "data")
print(f"corrdir: {corrdir}")

def add_config (ana: od.Analysis,
                campaign: od.Campaign,
                config_name           = None,
                config_id             = None,
                limit_dataset_files   = None
) -> od.Config :

    print(f"\ncampaign: {campaign.name}")
    print(f"config name: {config_name}, id: {config_id}")
    # gather campaign data
    year = campaign.x.year
    postfix = campaign.x.postfix
    
    # some validations
    assert year in [2022, 2023, 2024]

    # get all root processes
    procs = get_root_processes_from_campaign(campaign)
    
    # create a config by passing the campaign, so id and name will be identical
    cfg = ana.add_config(campaign,
                         name  = config_name,
                         id    = config_id)

    # add processes we are interested in
    process_names = [
        ## Data
        "data",
        ## W + jets
        "w_lnu",
        ## Drell-Yan
        "dy_lep_m50",
        ## TTJets
        "tt",
        ## Single top
        "st_tchannel",
        "st_twchannel",
        ## VV [diboson inclusive]
        "vv",
        ## Signal
        "h_ggf_tautau",
    ]

    for process_name in process_names:
        # development switch in case datasets are not _yet_ there
        if process_name not in procs:
            print(f"WARNING: {process_name} not in cmsdb processes")
            continue
        # add the process
        #if process_name == "h_ggf_tautau":
        #    procs.get(process_name).is_signal = True
        proc = cfg.add_process(procs.get(process_name))
        #if proc.name == "h_ggf_tautau":
        #    proc.is_signal = True
    # configuration of colors, labels, etc. can happen here
    from httcp.config.styles import stylize_processes
    stylize_processes(cfg)    

    # add datasets we need to study
    dataset_names = [
        ##W+jets
        "wj_incl",
        #Drell-Yan
        "dy_lep_m50",
        ## ttbar
        "tt_sl",
        "tt_dl",
        "tt_fh",
        ##single top
        "st_tchannel_t",
        "st_tchannel_tbar",
        "st_tw_t_sl",
        "st_tw_tb_sl",
        "st_tw_t_dl",
        "st_tw_tb_dl",
        "st_tw_t_fh",
        "st_tw_tb_fh",
        ##Diboson
        "ww",
        "wz",
        "zz",
        ## single top tW channel
        #"st_t_wminus_to_2l2nu",
        #"st_tbar_wplus_to_lnu2q",
        #"st_tbar_wplus_to_2l2nu",
        "h_ggf_tautau_prod_cp_even_sm",
    ]
    
    datasets_data = []
    if year == 2022:
        if postfix == "preEE":
            datasets_data = ["data_e_C",   "data_e_D",
                             "data_single_mu_C",
                             "data_mu_C",  "data_mu_D", 
                             "data_tau_C", "data_tau_D"]

        elif postfix == "postEE":
            datasets_data = ["data_e_E",   "data_e_G",
                             "data_mu_E",  "data_mu_G",
                             "data_tau_E", "data_tau_F", "data_tau_G"]
        else:
            raise RuntimeError(f"Wrong postfix: {campaign.x.postfix}")

    elif year == 2023:
        if postfix == "preBPix":
            datasets_data = ["data_e0_C",  "data_e1_C",
                             "data_mu0_C", "data_mu1_C",
                             "data_tau_C"]
        elif postfix == "postBPix":
            datasets_data = ["data_e0_D",  "data_e1_D",
                             "data_mu0_D", "data_mu1_D",
                             "data_tau_D"]
        else:
            raise RuntimeError(f"Wrong postfix: {campaign.x.postfix}")

    elif year == 2024:
        raise RuntimeError("2024? Too Early ...")

    else:
        raise RuntimeError(f"Check Year in __init__.py in cmsdb campaign: {year}")


    dataset_names = datasets_data + dataset_names
    for dataset_name in dataset_names:
        # development switch in case datasets are not _yet_ there
        if dataset_name not in campaign.datasets:
            print(f"WARNING: {dataset_name} not in cmsdb campaign")
            continue
        
        # add the dataset
        dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))
        if re.match(r"^(ww|wz|zz)_.*pythia$", dataset.name):
            dataset.add_tag("no_lhe_weights")
        
        # for testing purposes, limit the number of files to 1
        for info in dataset.info.values():
            if limit_dataset_files:
                info.n_files = min(info.n_files, limit_dataset_files)

    # verify that the root process of all datasets is part of any of the registered processes
    verify_config_processes(cfg, warn=True)

    # default objects, such as calibrator, selector, producer, ml model, inference model, etc
    cfg.x.default_calibrator      = "main"
    cfg.x.default_selector        = "main"
    cfg.x.default_producer        = "main"
    cfg.x.default_ml_model        = None
    cfg.x.default_inference_model = "example"
    cfg.x.default_categories      = ("incl",)
    #cfg.x.default_variables = ("n_jet", "jet1_pt")
    cfg.x.default_variables       = ("event","channel_id")

    # process groups for conveniently looping over certain processs
    # (used in wrapper_factory and during plotting)
    cfg.x.process_groups = {
        "backgrounds": (backgrounds := [
            "w_lnu",
            "tt",
            "dy_lep_m50",
            "st",
            "vv",
        ]),
        "bkg_sig"       : (bkg_sig       := [*backgrounds, "h_ggf_tautau"]),
        "data_bkg"      : (data_bkg      := [*backgrounds, "data"]),
        "data_bkg_sig"  : (data_bkg_sig  := [*backgrounds, "data", "h_ggf_tautau"]),
    }
    cfg.x.process_settings_groups = {
        "unstack_processes": {proc: {"unstack": True, "scale": 10.0} for proc in ("h_ggf_tautau")},
    }
    cfg.x.general_settings_groups = {
        "compare_shapes": {"skip_ratio": False, "shape_norm": False, "yscale": "log"} #"cms_label": "simpw"},
    }

    # dataset groups for conveniently looping over certain datasets
    # (used in wrapper_factory and during plotting)
    cfg.x.dataset_groups = {}

    # category groups for conveniently looping over certain categories
    # (used during plotting)
    cfg.x.category_groups = {}

    # variable groups for conveniently looping over certain variables
    # (used during plotting)
    cfg.x.variable_groups = {}

    # shift groups for conveniently looping over certain shifts
    # (used during plotting)
    cfg.x.shift_groups = {}

    # selector step groups for conveniently looping over certain steps
    # (used in cutflow tasks)
    cfg.x.selector_step_groups = {
        "default": ["json", "met_filter", "dl_res_veto", "trigger", "lepton", "jet"],
    }

    # whether to validate the number of obtained LFNs in GetDatasetLFNs
    # (currently set to false because the number of files per dataset is truncated to 2)
    cfg.x.validate_dataset_lfns = False
    
    # lumi values in inverse pb
    # TODO later: preliminary luminosity using norm tag. Must be corrected, when more data is available
    # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis
    if year == 2022:
        if postfix == "preEE":
            cfg.x.luminosity = Number(7980.4, {
                "total": 0.014j,
            })
        elif postfix == "postEE":
            cfg.x.luminosity = Number(26671.7, {
                "total": 0.014j,
            })
        else:
            raise RuntimeError(f"Wrong postfix: {campaign.x.postfix}")

    elif year == 2023:
        if postfix == "preBPix":
            cfg.x.luminosity = Number(17794 , {
                "lumi_13TeV_correlated": 0.0j,
            })
        elif postfix == "postBPix":
            cfg.x.luminosity = Number(9451, {
                "lumi_13TeV_correlated": 0.0j,
            })
    elif year == 2024:
        cfg.x.luminosity = Number(0, {
            "lumi_13TeV_correlated": 0.0j,
        })
    else:
        raise RuntimeError(f"Wrong year: {year}")


    # minimum bias cross section in mb (milli) for creating PU weights, values from
    # https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData?rev=45#Recommended_cross_section
    # TODO later: Error for run three not available yet. Using error from run 2.
    cfg.x.minbias_xs = Number(69.2, 0.046j)

    year_postfix = ""
    if year == 2022:
        year_postfix = "EE" if postfix == "postEE" else ""
    elif year == 2023:
        year_postfix = "BPix" if postfix == "postBPix" else ""
        
    # b-tag working points
    # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
    # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22EE/
    # TODO later: complete WP when data becomes available
    btag_key = f"{year}{year_postfix}"
    cfg.x.btag_working_points = DotDict.wrap({
        "deepjet": {
            "loose"  : {"2022": 0.0583, "2022EE": 0.0614, "2023": 0.0583, "2023BPix": 0.0583, "2024": 0.0}[btag_key],
            "medium" : {"2022": 0.3086, "2022EE": 0.3196, "2023": 0.3086, "2023BPix": 0.3086, "2024": 0.0}[btag_key],
            "tight"  : {"2022": 0.7183, "2022EE": 0.73,   "2023": 0.7183, "2023BPix": 0.7183, "2024": 0.0}[btag_key],
            "vtight" : {"2022": 0.8111, "2022EE": 0.8184, "2023": 0.8111, "2023BPix": 0.8111, "2024": 0.0}[btag_key],
            "vvtight": {"2022": 0.9512, "2022EE": 0.9542, "2023": 0.9512, "2023BPix": 0.9512, "2024": 0.0}[btag_key],
        },
        "robustParticleTransformer": {
            "loose"  : {"2022": 0.0849, "2022EE": 0.0897, "2023": 0.0849, "2023BPix": 0.0849, "2024": 0.0}[btag_key],
            "medium" : {"2022": 0.4319, "2022EE": 0.451,  "2023": 0.4319, "2023BPix": 0.4319, "2024": 0.0}[btag_key],
            "tight"  : {"2022": 0.8482, "2022EE": 0.8604, "2023": 0.8482, "2023BPix": 0.8482, "2024": 0.0}[btag_key],
            "vtight" : {"2022": 0.9151, "2022EE": 0.9234, "2023": 0.9151, "2023BPix": 0.9151, "2024": 0.0}[btag_key],
            "vvtight": {"2022": 0.9874, "2022EE": 0.9893, "2023": 0.9874, "2023BPix": 0.9874, "2024": 0.0}[btag_key],
        },
        "particleNet": {
            "loose"  : {"2022": 0.047,  "2022EE": 0.0499, "2023": 0.0499, "2023BPix": 0.0499, "2024": 0.0}[btag_key],
            "medium" : {"2022": 0.245,  "2022EE": 0.2605, "2023": 0.2605, "2023BPix": 0.2605, "2024": 0.0}[btag_key],
            "tight"  : {"2022": 0.6734, "2022EE": 0.6915, "2023": 0.6915, "2023BPix": 0.6915, "2024": 0.0}[btag_key],
            "vtight" : {"2022": 0.7862, "2022EE": 0.8033, "2023": 0.8033, "2023BPix": 0.8033, "2024": 0.0}[btag_key],
            "vvtight": {"2022": 0.961,  "2022EE": 0.9664, "2023": 0.9664, "2023BPix": 0.9664, "2024": 0.0}[btag_key],
        },
    })
    # 2023 is dummy ... CORRECT IT LATER ===>> TODO

    cfg.x.deep_tau_tagger = "DeepTau2018v2p5"
    cfg.x.deep_tau_info = DotDict.wrap({
        "DeepTau2018v2p5": {
            "wp": {
                "vs_e": {"VVVLoose": 1, "VVLoose": 2, "VLoose": 3, "Loose": 4, "Medium": 5, "Tight": 6, "VTight": 7, "VVTight": 8},
                "vs_m": {"VVVLoose": 1, "VVLoose": 1, "VLoose": 1, "Loose": 2, "Medium": 3, "Tight": 4, "VTight": 4, "VVTight": 4},
                "vs_j": {"VVVLoose": 1, "VVLoose": 2, "VLoose": 3, "Loose": 4, "Medium": 5, "Tight": 6, "VTight": 7, "VVTight": 8},
            },
            "vs_e": "VVLoose",
            "vs_m": "Tight",
            "vs_j": "Medium"
        },
    })


  
    # Adding triggers
    if year == 2022:
        from httcp.config.triggers import add_triggers_run3_2022
        add_triggers_run3_2022(cfg, postfix)
    elif year == 2023:
        from httcp.config.triggers import add_triggers_run3_2023
        add_triggers_run3_2023(cfg, postfix)
    elif year == 2024:
        raise RuntimeError("Babushcha: too early")
    else:
        raise RuntimeError(f"Wrong year: {year}. Check __init__.py in cmsdb campaign")


    from httcp.config.met_filters import add_met_filters
    add_met_filters(cfg)

    year2 = year%100
    #external_path_parent = os.path.join(os.environ.get('HTTCP_BASE'), f"httcp/data/corrections")
    external_path_parent = os.path.join(corrdir, "corrections")
    external_path_tail   = f"{year}{postfix}" if postfix else f"{year}"
    external_path        = os.path.join(external_path_parent, f"{external_path_tail}")
    normtag_path = "/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags"

    print(f"external_path_parent : {external_path_parent}")
    print(f"external_path        : {external_path}")

    #json_mirror = os.path.join(os.environ.get('HTTCP_BASE'), f"httcp/data/jsonpog-integration")
    json_mirror = os.path.join(corrdir, "jsonpog-integration")
    print(f"json_mirror          : {json_mirror}")

    normtagjson = None
    if year == 2022:
        normtagjson = "/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json"
    elif year == 2023:
        normtagjson = "/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json"
    elif year == 2024:
        raise RuntimeWarning("too early")
    else:
        raise RuntimeError(f"Check year : {year}")

    goldenjson  = glob(f"{external_path}/Lumi/*.json")[0]

    print(f"GoldenJSON           : {goldenjson}")
    print(f"NormtagJSON          : {normtagjson}")

    cfg.x.external_files = DotDict.wrap({
        # lumi files
        "lumi": {
            "golden"  : (goldenjson,  "v1"),  # noqa
            "normtag" : (normtagjson, "v1"),
        },

        # pileup weight corrections
        "pu_sf": (f"{json_mirror}/POG/LUM/{year}_Summer{year2}{year_postfix}/puWeights.json.gz", "v1"),

        # jet energy correction
        "jet_jerc": (f"{json_mirror}/POG/JME/{year}_Summer{year2}{year_postfix}/jet_jerc.json.gz", "v1"),

        # Add Muon POG scale factors
        "muon_sf": (f"{json_mirror}/POG/MUO/{year}_Summer{year2}{year_postfix}/muon_Z.json.gz", "v1"),
        #"muon_sf": (f"{external_path_parent}/muon_SFs_2022_preEE.root", "v1"),
        
        # electron scale factors
        "electron_sf": (f"{json_mirror}/POG/EGM/{year}_Summer{year2}{year_postfix}/electron.json.gz", "v1"),
        
        # tau energy correction and scale factors
        #"tau_sf": (f"{external_path_parent}/tau_DeepTau2018v2p5_2022_preEE.json.gz", "v1"),  # noqa
        "tau_correction": (f"{external_path_parent}/tau_DeepTau2018v2p5_2022_preEE.json.gz", "v1"),  # noqa

        # btag scale factor
        #"btag_sf_corr": (f"{json_mirror}/POG/BTV/{year}_Summer{year2}{year_postfix}/btagging.json.gz", "v1"),

        # met phi corrector 
        # unavailable for Run3
        # "met_phi_corr": (f"{json_mirror}/POG/JME/2018_UL/met.json.gz", "v1"),
    })

    # register shifts
    cfg.add_shift(name="nominal", id=0)

    """
    # ----------------------------------------------------------------------- #
    #                       jec & jer configuration                           #
    # ----------------------------------------------------------------------- #
    # https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC?rev=201
    # TODO later: check this corrections summary correction_file (jet_jerc.json.gz) after setting sandbox_dev
    # https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/JME?ref_type=heads
    cfg.x.jec = DotDict.wrap({
        "campaign": f"Summer{year2}{year_postfix}_22Sep2023",                  # campaign name for this JEC corrections
        "version": {2022: "V2", 2023: "V1"}[year],                             # version of the corrections 
        "jet_type": "AK4PFPuppi",                                              # Type of jets that the corrections should be applied on
        "levels": ["L1FastJet", "L2Relative", "L2L3Residual", "L3Absolute"],   # relevant levels in the derivation process of the JEC
        "levels_for_type1_met": ["L1FastJet"],                                 # relevant levels in the derivation process of the Type 1 MET JEC
        "uncertainty_sources": [                                               # names of the uncertainties to be applied
            "Total",
            "CorrelationGroupMPFInSitu",
            "CorrelationGroupIntercalibration",
            "CorrelationGroupbJES",
            "CorrelationGroupFlavor",
            "CorrelationGroupUncorrelated",
        ],
    })

    # JER
    # https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution?rev=107 # TODO later: check this
    cfg.x.jer = DotDict.wrap({
        "campaign": f"Summer{year2}{year_postfix}_22Sep2023",
        "version": "JR" + {2022: "V1", 2023: "V1"}[year],
        "jet_type": "AK4PFPuppi",
    })
    """

    # ----------------------------------------------------------------------- #
    #                                electron id SF                           #
    # ----------------------------------------------------------------------- #
    # names of electron correction sets and working points
    # (used in the electron_sf producer)    
    electron_sf_tag = ""
    if year == 2022:
        electron_sf_tag = "2022Re-recoE+PromptFG" if year_postfix else "2022Re-recoBCD"
    elif year == 2023:
        electron_sf_tag = "2023PromptC" if year_postfix else "2023PromptD"
    
    cfg.x.electron_sf_names = (
        "Electron-ID-SF",
        electron_sf_tag,
        "wp80iso",
    )

    # ----------------------------------------------------------------------- #
    #                                muon id SF                               #
    # ----------------------------------------------------------------------- #
    # names of muon correction sets and working points
    # (used in the muon producer)
    cfg.x.muon_sf_names = (
        "NUM_TightPFIso_DEN_MediumID", 
        f"{year}_{postfix}"
    )

    """
    # register shifts
    cfg.add_shift(name="nominal", id=0)

    cfg.add_shift(name="minbias_xs_up", id=7, type="shape")
    cfg.add_shift(name="minbias_xs_down", id=8, type="shape")
    add_shift_aliases(
        cfg,
        "minbias_xs",
        {
            "pu_weight": "pu_weight_{name}",
            "normalized_pu_weight": "normalized_pu_weight_{name}",
        },
    )

    # load jec sources
    with open(os.path.join(thisdir, "jec_sources.yaml"), "r") as f:
        all_jec_sources = yaml.load(f, yaml.Loader)["names"]

    for jec_source in cfg.x.jec.uncertainty_sources:
        idx = all_jec_sources.index(jec_source)
        cfg.add_shift(
            name=f"jec_{jec_source}_up",
            id=5000 + 2 * idx,
            type="shape",
            tags={"jec"},
            aux={"jec_source": jec_source},
        )
        cfg.add_shift(
            name=f"jec_{jec_source}_down",
            id=5001 + 2 * idx,
            type="shape",
            tags={"jec"},
            aux={"jec_source": jec_source},
        )
        add_shift_aliases(
            cfg,
            f"jec_{jec_source}",
            {
                "Jet.pt": "Jet.pt_{name}",
                "Jet.mass": "Jet.mass_{name}",
                "MET.pt": "MET.pt_{name}",
                "MET.phi": "MET.phi_{name}",
            },
        )

    cfg.add_shift(name="jer_up", id=6000, type="shape", tags={"jer"})
    cfg.add_shift(name="jer_down", id=6001, type="shape", tags={"jer"})
    add_shift_aliases(
        cfg,
        "jer",
        {
            "Jet.pt": "Jet.pt_{name}",
            "Jet.mass": "Jet.mass_{name}",
            "MET.pt": "MET.pt_{name}",
            "MET.phi": "MET.phi_{name}",
        },
    )
    """
    """
    for i, (match, dm) in enumerate(itertools.product(["jet", "e"], [0, 1, 10, 11])):
        cfg.add_shift(name=f"tec_{match}_dm{dm}_up", id=20 + 2 * i, type="shape", tags={"tec"})
        cfg.add_shift(name=f"tec_{match}_dm{dm}_down", id=21 + 2 * i, type="shape", tags={"tec"})
        add_shift_aliases(
            cfg,
            f"tec_{match}_dm{dm}",
            {
                "Tau.pt": "Tau.pt_{name}",
                "Tau.mass": "Tau.mass_{name}",
                "MET.pt": "MET.pt_{name}",
                "MET.phi": "MET.phi_{name}",
            },
        )
    """
    """
    # start at id=50
    cfg.x.tau_unc_names = [
        "jet_dm0", "jet_dm1", "jet_dm10",
        "e_barrel", "e_endcap",
        "mu_0p0To0p4", "mu_0p4To0p8", "mu_0p8To1p2", "mu_1p2To1p7", "mu_1p7To2p3",
    ]
    for i, unc in enumerate(cfg.x.tau_unc_names):
        cfg.add_shift(name=f"tau_{unc}_up", id=50 + 2 * i, type="shape")
        cfg.add_shift(name=f"tau_{unc}_down", id=51 + 2 * i, type="shape")
        add_shift_aliases(cfg, f"tau_{unc}", {"tau_weight": f"tau_weight_{unc}_" + "{direction}"})
    """
    cfg.add_shift(name="e_up", id=90, type="shape")
    cfg.add_shift(name="e_down", id=91, type="shape")
    add_shift_aliases(cfg, "e", {"electron_weight": "electron_weight_{direction}"})

    cfg.add_shift(name="mu_up", id=100, type="shape")
    cfg.add_shift(name="mu_down", id=101, type="shape")
    add_shift_aliases(cfg, "mu", {"muon_weight": "muon_weight_{direction}"})

    cfg.add_shift(name="pdf_up", id=130, type="shape")
    cfg.add_shift(name="pdf_down", id=131, type="shape")
    add_shift_aliases(
        cfg,
        "pdf",
        {
            "pdf_weight": "pdf_weight_{direction}",
            #"normalized_pdf_weight": "normalized_pdf_weight_{direction}",
        },
    )


    # target file size after MergeReducedEvents in MB
    cfg.x.reduced_file_size = 512.0
    
    from httcp.config.variables import keep_columns
    keep_columns(cfg)


    # event weight columns as keys in an OrderedDict, mapped to shift instances they depend on
    # get_shifts = functools.partial(get_shifts_from_sources, cfg)
    # configurations for all possible event weight columns as keys in an OrderedDict,
    # mapped to shift instances they depend on
    # (this info is used by weight producers)
    get_shifts = functools.partial(get_shifts_from_sources, cfg)
    cfg.x.event_weights = DotDict({
        "normalization_weight": [],
        #"pdf_weight": get_shifts("pdf"),
        #"normalized_pu_weight": get_shifts("minbias_xs"),
        #"electron_weight": get_shifts("e"),
        #"muon_weight": get_shifts("mu"),
        #"tau_weight": get_shifts(*(f"tau_{unc}" for unc in cfg.x.tau_unc_names)),
    })

    cfg.x.default_weight_producer = "all_weights"

    # versions per task family, either referring to strings or to callables receving the invoking
    # task instance and parameters to be passed to the task family
    def set_version(cls, inst, params):
        # per default, use the version set on the command line
        version = inst.version 
        return version if version else 'dev1'
            
        
    cfg.x.versions = {
        "cf.CalibrateEvents"    : set_version,
        "cf.SelectEvents"       : set_version,
        "cf.MergeSelectionStats": set_version,
        "cf.MergeSelectionMasks": set_version,
        "cf.ReduceEvents"       : set_version,
        "cf.MergeReductionStats": set_version,
        "cf.MergeReducedEvents" : set_version,
    }
    # channels
    cfg.add_channel(name="etau",   id=1)
    cfg.add_channel(name="mutau",  id=2)
    cfg.add_channel(name="tautau", id=4)
    
    campaign_tag = cfg.campaign.x("custom").get("creator")
    if campaign_tag == "desy" or campaign_tag == "IPHC":
        def get_dataset_lfns(dataset_inst: od.Dataset, shift_inst: od.Shift, dataset_key: str) -> list[str]:
            # destructure dataset_key into parts and create the lfn base directory
            print(f"Creating custom get_dataset_lfns for {config_name}")   
            try:
               basepath = cfg.campaign.x("custom").get("location")
            except:
                print("Did not find any basebath in the campaigns")
                basepath = "" 
            lfn_base = law.wlcg.WLCGDirectoryTarget(
                f"{basepath}{dataset_key}",
                fs="wlcg_fs_eoscms_redirector",
            )
            print(f"lfn basedir:{lfn_base}")
            # loop though files and interpret paths as lfns
            return [
                lfn_base.child(basename, type="f").path
                for basename in lfn_base.listdir(pattern="*.root")
            ]
        # define the lfn retrieval function
        cfg.x.get_dataset_lfns = get_dataset_lfns
        # define a custom sandbox
        cfg.x.get_dataset_lfns_sandbox = dev_sandbox("bash::$CF_BASE/sandboxes/cf.sh")
        # define custom remote fs's to look at
        cfg.x.get_dataset_lfns_remote_fs =  lambda dataset_inst: "wlcg_fs_eoscms_redirector"
        
    # add categories using the "add_category" tool which adds auto-generated ids
    from httcp.config.categories import add_categories
    add_categories(cfg)
        
    from httcp.config.variables import add_variables
    add_variables(cfg)
    
    
