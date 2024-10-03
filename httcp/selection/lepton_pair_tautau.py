# coding: utf-8

"""
Prepare h-Candidate from SelectionResult: selected lepton indices & channel_id [trigger matched] 
"""

from typing import Optional
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column


np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")



def get_sorted_pair(
        dtrpairs: ak.Array,
        dtrpairindices: ak.Array,
)->ak.Array:
    # redundant, because taus were sorted by the deeptau before
    sorted_idx = ak.argsort(dtrpairs["0"].rawDeepTau2018v2p5VSjet, ascending=False)
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]

    # if the deep tau val of tau-0 is the same for the first two pair
    where_same_iso_1 = (
        ak.firsts(dtrpairs["0"].rawDeepTau2018v2p5VSjet[:,:1], axis=1)
        ==
        ak.firsts(dtrpairs["0"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1)
    )

    # if so, sort the pairs according to the deep tau of the 2nd tau
    sorted_idx = ak.where(where_same_iso_1,
                          ak.argsort(dtrpairs["1"].rawDeepTau2018v2p5VSjet, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]
    # if the deep tau val of tau-1 is the same for the first two pair 
    where_same_iso_2 = (
        ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,:1], axis=1)
        ==
        ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1)
    )

    # sort them with the pt of the 1st tau
    sorted_idx = ak.where(where_same_iso_2,
                          ak.argsort(dtrpairs["0"].pt, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]
    # check if the first two pairs have the second tau with same rawDeepTau2017v2p1VSjet
    where_same_pt_1 = (
        ak.firsts(dtrpairs["0"].pt[:,:1], axis=1)
        ==
        ak.firsts(dtrpairs["0"].pt[:,1:2], axis=1)
    )

    # if so, sort the taus with their pt
    sorted_idx = ak.where(where_same_pt_1,
                          ak.argsort(dtrpairs["1"].pt, ascending=False),
                          sorted_idx)

    # finally, the pairs are sorted
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]

    lep1 = ak.singletons(ak.firsts(dtrpairs["0"], axis=1))
    lep2 = ak.singletons(ak.firsts(dtrpairs["1"], axis=1))

    lep1idx = ak.singletons(ak.firsts(dtrpairindices["0"], axis=1))
    lep2idx = ak.singletons(ak.firsts(dtrpairindices["1"], axis=1))
    # print(f"lep1 pt: {lep1.pt}")
    # print(f"lep2 pt: {lep2.pt}")
    dtrpair    = ak.concatenate([lep1, lep2], axis=1) 
    dtrpairidx = ak.concatenate([lep1idx, lep2idx], axis=1)     

    return dtrpairidx



@selector(
    uses={
        "Tau.pt", "Tau.eta", "Tau.phi", "Tau.mass",
        "Tau.charge", "Tau.rawDeepTau2018v2p5VSjet",
        "Tau.idDeepTau2018v2p5VSjet", "Tau.idDeepTau2018v2p5VSe", "Tau.idDeepTau2018v2p5VSmu",
    },
    exposed=False,
)
def tautau_selection(
        self: Selector,
        events: ak.Array,
        lep_indices: ak.Array,
        tau_iso_mask: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult, ak.Array]:

    # # prepare vectors for output vectors
    false_mask = (abs(events.event) < 0)
    tau2_isolated = false_mask
    leptons_os = false_mask
    channel_id = np.uint8(1) * false_mask


    # #from IPython import embed; embed()

    # # Extra channel specific selections on tau
    # # -------------------- #
    # taus            = events.Tau[lep_indices]
    # tau_tagger      = self.config_inst.x.deep_tau_tagger
    # tau_tagger_wps  = self.config_inst.x.deep_tau_info[tau_tagger].wp
    # is_good_tau     = (
    #     (taus.idDeepTau2018v2p5VSjet   >= tau_tagger_wps.vs_j.Medium)
    #     & (taus.idDeepTau2018v2p5VSe   >= tau_tagger_wps.vs_e.VVLoose)
    #     & (taus.idDeepTau2018v2p5VSmu  >= tau_tagger_wps.vs_m.VLoose)
    # )
    # lep_indices    = lep_indices[is_good_tau]
    # # -------------------- # 

    
    # Sorting leps [Tau] by deeptau [descending]
    lep_sort_key       = events.Tau[lep_indices].rawDeepTau2018v2p5VSjet
    lep_sorted_indices = ak.argsort(lep_sort_key, axis=-1, ascending=False)
    lep_indices        = lep_indices[lep_sorted_indices]

    leps_pair        = ak.combinations(events.Tau[lep_indices], 2, axis=1)
    lep_indices_pair = ak.combinations(lep_indices, 2, axis=1)

    
    # pair of leptons: probable higgs candidate -> leps_pair
    # and their indices                         -> lep_indices_pair 
    lep1, lep2 = ak.unzip(leps_pair)
    lep1_idx, lep2_idx = ak.unzip(lep_indices_pair)

    # # determine the os/ss charge sign relation
    is_os = (lep1.charge * lep2.charge) < 0

    # Check if the number of isolated taus per event is 2 or more (passing Medium WPvsJet)
    is_iso = ak.sum(tau_iso_mask, axis=1) >= 2

    tau2_isolated = is_iso

    preselection = {
        "tautau_is_pt_40"      : (lep1.pt > 40) & (lep2.pt > 40),
        "tautau_is_eta_2p1"    : (np.abs(lep1.eta) < 2.1) & (np.abs(lep2.eta) < 2.1),
        # "leptons_os"         : (lep1.charge * lep2.charge) < 0,
        "tautau_dr_0p5"        : (1*lep1).delta_r(1*lep2) > 0.5,  #deltaR(lep1, lep2) > 0.5,
        # "tau2_isolated"        : is_iso
    }

    # Initialize a mask to identify good pairs based on the condition that `lep1_idx` is non-negative.
    good_pair_mask = lep1_idx >= 0

    # Create a dictionary to store the selection steps for pairs.
    pair_selection_steps = {}
    for cut in preselection.keys():
        # Update the good_pair_mask by applying each cut in conjunction with the current mask.
        good_pair_mask = good_pair_mask & preselection[cut]
        # Store the current state of the good_pair_mask in the pair_selection_steps dictionary.
        pair_selection_steps[cut] = good_pair_mask


    # Select pairs of leptons based on the good_pair_mask, filtering the original pairs.    
    leps_pair_sel = leps_pair[good_pair_mask]
    lep_indices_pair_sel = lep_indices_pair[good_pair_mask]

    # Extract the first index of each lepton in the selected pairs, treating them as singletons.
    lep1idx = ak.singletons(ak.firsts(lep_indices_pair_sel["0"], axis=1))
    lep2idx = ak.singletons(ak.firsts(lep_indices_pair_sel["1"], axis=1))

    
    # Concatenate the indices of the selected leptons into a single array.
    lep_indices_pair_sel_single = ak.concatenate([lep1idx, lep2idx], axis=1)

    # Identify pairs with more than one lepton.
    where_many   = ak.num(lep_indices_pair_sel, axis=1) > 1
    # For pairs with multiple leptons, sort them; otherwise, use the single concatenated indices.
    pair_indices = ak.where(where_many, 
                            get_sorted_pair(leps_pair_sel,
                                            lep_indices_pair_sel),
                            lep_indices_pair_sel_single)

    
    # Update the isolation status and OS selection
    tau2_isolated = ak.where(good_pair_mask, is_iso, False)  # Set based on good pairs
    leptons_os = ak.where(good_pair_mask, is_os, False)  # Set based on good pairs
    leptons_os = ak.fill_none(leptons_os, False)

    # Save new columns in the events object
    events = set_ak_column(events, "leptons_os", leptons_os)
    events = set_ak_column(events, "tau2_isolated", tau2_isolated)


    # Return a selection result containing the auxiliary steps and the calculated pair indices.
    return events, SelectionResult(
        steps = pair_selection_steps,
        objects = {
            "Tau": {
                "Tau": 
            }
        }
        aux = {
            "tautau_pair" = ak.concatenate([events.Tau[tautau_indices_pair[:,0:1]], 
                                  events.Tau[tautau_indices_pair[:,1:2]]], 
                                 axis=1)
        },
    ), pair_indices
