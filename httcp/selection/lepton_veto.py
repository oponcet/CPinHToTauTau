# coding: utf-8

"""
Extra-Lepton-Veto
https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2325&ancode=HIG-20-006&tp=an&line=HIG-20-006
http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_192_v15.pdf
"""

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import set_ak_column
from columnflow.util import maybe_import, DotDict



np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

@selector(
    uses={
        "Muon.pt", "Muon.eta", "Muon.phi", "Muon.mass",
        "Electron.pt", "Electron.eta", "Electron.phi", "Electron.mass",
    },
    exposed=False,
)
def extra_lepton_veto(
        self: Selector,
        events: ak.Array,
        extra_electron_index: ak.Array,
        extra_muon_index: ak.Array,
        hcand_pair: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    """
    The events are rejected if they contain a (very loosely selected) third lepton, separated from the
    two tau candidates by DR > 0.5
    """
    # extra leptons
    extra_lep  = ak.Array(ak.concatenate([events.Muon[extra_muon_index],
                                          events.Electron[extra_electron_index]], axis=-1),
                          behavior=coffea.nanoevents.methods.nanoaod.behavior)
    extra_lep  = ak.with_name(extra_lep, "PtEtaPhiMLorentzVector")

    # has pair or not
    # e.g.
    # hcand_pair = [
    #                [[p4,p4], [      ], [     ] ],
    #                [[     ], [      ], [p4,p4] ],
    #                [[p4,p4], [      ], [p4,p4] ]
    #              ]
    # is_pair = [
    #             [ True, False, False ],
    #             [ False, False, True ],
    #             [ True, False, True  ]
    #           ]
    # has_single_pair = [ True, True, False ]
    is_pair = ak.num(hcand_pair.pt, axis=-1) == 2
    has_single_pair = ak.sum(ak.num(hcand_pair.pt, axis=-1), axis=1) == 2

    hcand_pair_p4 = 1 * hcand_pair
    hcand_pair_p4_dummy = hcand_pair_p4[:,1][:,:0][:,None]
    masked_hcand_pair_p4 = ak.mask(hcand_pair_p4, is_pair)
    hcand_pair_p4 = ak.where(has_single_pair, masked_hcand_pair_p4, hcand_pair_p4_dummy)
    hcand_pair_p4 = ak.drop_none(hcand_pair_p4)[:,0]
    hcand_lep1 = hcand_pair_p4[:,0:1]
    hcand_lep2 = hcand_pair_p4[:,1:2]

    dr_hlep1_extraleps = extra_lep.metric_table(hcand_lep1)
    dr_hlep2_extraleps = extra_lep.metric_table(hcand_lep2)

    dr_mask = (
        (dr_hlep1_extraleps > 0.5) 
        & (dr_hlep2_extraleps > 0.5)
    )
    
    # this is kept as False intentionally,
    # because the end result will be then true for tautau channel
    dummy = (events.event < 0)[:,None]
    # extra leps must be away from both leptons of the higgs candidates
    has_extra_lepton = ak.where(has_single_pair, 
                                ak.all(dr_mask, axis=-1),
                                dummy)


    has_no_extra_lepton = ak.sum(has_extra_lepton, axis=1) == 0

    # For the purpose of debugging
    if self.config_inst.x.verbose.selection.extra_lep_veto:
        for i in range(100):
            if events.channel_id[i] < 1: continue
            print(f"event : {events.event[i]}")
            print(f"hcand_pairs_pt: {hcand_pair.pt[i]}")
            print(f"extra leps pt: {extra_lep.pt[i]}")
            print(f"h1 pt: {hcand_lep1.pt[i]}, h2 pt: {hcand_lep2.pt[i]}")
            print(f"dr_hlep1_extraleps : {dr_hlep1_extraleps[i]}")
            print(f"dr_hlep2_extraleps : {dr_hlep2_extraleps[i]}")
            print(f"dr_mask : {dr_mask[i]}")
            print(f"has_extra_lepton : {has_extra_lepton[i]}")
            print(f"has_no_extra_lepton : {has_no_extra_lepton[i]}")
            print("\n")

    return events, SelectionResult(
        steps={"extra_lepton_veto": has_no_extra_lepton}
    )



@selector(
    uses={
        "Muon.pt", "Muon.eta", "Muon.phi", "Muon.mass", "Muon.charge",
        "Electron.pt", "Electron.eta", "Electron.phi", "Electron.mass", "Electron.charge",
    },
    exposed=False,
)
def double_lepton_veto(
        self: Selector,
        events: ak.Array,
        double_veto_electron_index: ak.Array,
        double_veto_muon_index: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    """
    A veto on events containing dileptons pairs is applied.
    The oppositesign leptons must be separated by DR > 0.15, 
    and a loose selection on the leptons is applied.
    """
    # get the double_veto_muon from events.Muon using double_veto_muon_index
    double_veto_muon     = ak.Array(events.Muon[double_veto_muon_index], 
                                    behavior=coffea.nanoevents.methods.nanoaod.behavior)
    # get the double_veto_muon p4s
    double_veto_muon     = ak.with_name(double_veto_muon, "PtEtaPhiMLorentzVector")

    # get the double_veto_electron from events.Electron using double_veto_electron_index
    double_veto_electron = ak.Array(events.Electron[double_veto_electron_index], 
                                    behavior=coffea.nanoevents.methods.nanoaod.behavior)
    # get the double_veto_electron p4s
    double_veto_electron = ak.with_name(double_veto_electron, "PtEtaPhiMLorentzVector")

    # get all possible combinatorics separately for double_veto_muons and electrons
    # mu : [ [(m11,m12),(m12,m13),(m11,m13)],
    #        [(m21,m22)],
    #        [],
    #        []
    #      ]
    mu_pair = ak.combinations(double_veto_muon, 2, axis=1)
    # el : [ [],
    #        [],
    #        [e31,e32],
    #        []
    #      ]
    el_pair = ak.combinations(double_veto_electron, 2, axis=1)

    # get the components
    # mu1 : [ [m11,m12,m11],
    #         [m21],
    #         [],
    #         [],
    #       ]
    # mu2 : [ [m12,m13,m13],
    #         [m22],
    #         [],
    #         [],
    #       ]
    mu1,mu2 = ak.unzip(mu_pair)
    # el1 : [ [],
    #         [],
    #         [e31],
    #         [],
    #       ]
    # el2 : [ [],
    #         [],
    #         [e32],
    #         [],
    #       ]
    el1,el2 = ak.unzip(el_pair)

    # main algo
    presel_mask = lambda leps1, leps2: ((leps1.charge * leps2.charge < 0) & (leps1.delta_r(leps2) > 0.15))

    # dlveto masks
    # dlveto_mu_mask : [ [True, False, False],
    #                    [False],
    #                    [],
    #                    []
    #                  ]
    dlveto_mu_mask = presel_mask(mu1,mu2)
    # dlveto_el_mask : [ [],
    #                    [],
    #                    [True],
    #                    [],
    #                  ]
    dlveto_el_mask = presel_mask(el1,el2)

    

    # concatenate the masks
    # It should not contain any True
    # dlveto_mask : [ [True, False, False],
    #                 [False],
    #                 [True]
    #               ]
    dlveto_mask = ak.concatenate([dlveto_mu_mask, dlveto_el_mask], axis=1)
    # this final event level mask should not conta
    # dl_veto : [ False, True, False, True ]
    dl_veto = ak.sum(dlveto_mask, axis=1) == 0
    
    if self.config_inst.x.verbose.selection.dilep_veto:
        for i in range(100):
            if events.channel_id[i] < 1: continue
            print(f"Event : {events.event[i]}")
            print(f"dl_veto_mu_pair_pt : {mu1.pt[i]}, {mu2.pt[i]}")
            print(f"is Z like pair?      {dlveto_mu_mask[i]}")
            print(f"dl_veto_el_pair_pt : {el1.pt[i]}, {el2.pt[i]}")
            print(f"is Z like pair?      {dlveto_el_mask[i]}")
            print(f"concat masks       : {dlveto_mask[i]}")
            print(f"Any True means Z like pair exists")
            print(f"has no Z like pair : {dl_veto[i]}\n")


    """
    # convert to event level masks
    ## dl_mu_veto : [ True, False, False ]
    ## dl_mu_veto = ak.any(dlveto_mu_mask, axis=1)
    ## dl_el_veto : [ False, False, True ]
    ## dl_el_veto = ak.any(dlveto_el_mask, axis=1)


    # dl_mu_veto : [ False, True, True, True ]
    # That means 1st event has Z like muon pair
    dl_mu_veto = ak.sum(dlveto_mu_mask, axis=1) == 0
    # dl_el_veto : [ True, True, False, True ]
    # 3rd events has Z like ele pair
    dl_el_veto = ak.sum(dlveto_el_mask, axis=1) == 0


    # combination of mu and el veto masks
    # So, for the example, 1st and 3rd event should be rejected
    # dl_veto : [ False, True, False, True ]
    dl_veto    = dl_mu_veto & dl_el_veto
    """
    
    return events, SelectionResult(
        steps={"dilepton_veto": dl_veto}
    )
    

