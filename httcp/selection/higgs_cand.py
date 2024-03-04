# coding: utf-8

"""
Prepare h-Candidate from SelectionResult: selected lepton indices & channel_id [trigger matched] 
"""

from typing import Optional
from columnflow.selection import SelectionResult
from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column

from hcp.util import invariant_mass, deltaR, transverse_mass

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")


def sort_and_get_pair_semilep(
        dtrpairs: ak.Array,
        dtrpairindices: ak.Array,
        var1: str, var2: str, var3: str, var4: str
)->ak.Array:
    sorted_idx = ak.argsort(dtrpairs["0"][var1], ascending=True)
    # Sort the pairs based on pfRelIso03_all of the first object in each pair
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]
    # Check if there are multiple pairs for each event
    where_many = ak.num(dtrpairs["0"], axis=1) > 1
    # Extract the pfRelIso03_all values for the first object in each pair
    lep1_pfRelIso03_all = dtrpairs["0"][var1]
    # Check if the pfRelIso03_all values are the same for the first two objects in each pair
    where_same_iso_1 = (
        where_many &
        ak.fill_none(
            (
                ak.firsts(dtrpairs["0"][var1][:,:1], axis=1) 
                ==
                ak.firsts(dtrpairs["0"][var1][:,1:2], axis=1)
            ), False
        )
    )
    # Sort the pairs based on pt if pfRelIso03_all is the same for the first two objects
    sorted_idx = ak.where(where_same_iso_1,
                          ak.argsort(dtrpairs["0"][var2], ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]
    # Check if the pt values are the same for the first two objects in each pair    
    where_same_pt_1 = (
        where_many &
        ak.fill_none(
            (
                ak.firsts(dtrpairs["0"][var2][:,:1], axis=1)
                ==
                ak.firsts(dtrpairs["0"][var2][:,1:2], axis=1)
            ), False
        )
    )
    # if so, sort the pairs with tau rawDeepTau2017v2p1VSjet
    sorted_idx = ak.where(where_same_pt_1,
                          ak.argsort(dtrpairs["1"][var3], ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]
    # check if the first two pairs have taus with same rawDeepTau2017v2p1VSjet
    where_same_iso_2 = (
        where_many &
        ak.fill_none(
            (
                ak.firsts(dtrpairs["1"][var3][:,:1], axis=1)
                ==
                ak.firsts(dtrpairs["1"][var3][:,1:2], axis=1)
            ), False
        )
    )
    # Sort the pairs based on pt if rawDeepTau2017v2p1VSjet is the same for the first two objects
    sorted_idx = ak.where(where_same_iso_2,
                          ak.argsort(dtrpairs["1"][var4], ascending=False),
                          sorted_idx)
    # finally, the pairs are sorted
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]

    # Extract the first object in each pair (lep1) and the second object (lep2)
    lep1 = ak.singletons(ak.firsts(dtrpairs["0"], axis=1))
    lep2 = ak.singletons(ak.firsts(dtrpairs["1"], axis=1))

    lep1idx = ak.singletons(ak.firsts(dtrpairindices["0"], axis=1))
    lep2idx = ak.singletons(ak.firsts(dtrpairindices["1"], axis=1))

    # Concatenate lep1 and lep2 to create the final dtrpair
    dtrpair    = ak.concatenate([lep1, lep2], axis=1)
    dtrpairidx = ak.concatenate([lep1idx, lep2idx], axis=1)

    return dtrpair, dtrpairidx



def sort_and_get_pair_fullhad(
        dtrpairs: ak.Array,
        dtrpairindices: ak.Array,
        var1: str, var2: str,
)->ak.Array:
    # redundant, because taus were sorted by the deeptau before
    sorted_idx = ak.argsort(dtrpairs["0"][var1], ascending=True)
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]
    where_many = ak.num(dtrpairs["0"], axis=1) > 1

    # if the deep tau val of tau-0 is the same for the first two pair
    where_same_iso_1 = (
        where_many &
        ak.fill_none(
            (
                ak.firsts(dtrpairs["0"][var1][:,:1], axis=1)
                ==
                ak.firsts(dtrpairs["0"][var1][:,1:2], axis=1)
            ), False
        )
    )

    # if so, sort the pairs according to the deep tau of the 2nd tau
    sorted_idx = ak.where(where_same_iso_1,
                          ak.argsort(dtrpairs["1"][var1], ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]
    # if the deep tau val of tau-1 is the same for the first two pair 
    where_same_iso_2 = (
        where_many &
        ak.fill_none(
            (
                ak.firsts(dtrpairs["1"][var1][:,:1], axis=1)
                ==
                ak.firsts(dtrpairs["1"][var1][:,1:2], axis=1)
            ), False
        )
    )

    # sort them with the pt of the 1st tau
    sorted_idx = ak.where(where_same_iso_2,
                          ak.argsort(dtrpairs["0"][var2], ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]
    # check if the first two pairs have the second tau with same rawDeepTau2017v2p1VSjet
    where_same_pt_1 = (
        where_many &
        ak.fill_none(
            (
                ak.firsts(dtrpairs["0"][var2][:,:1], axis=1)
                ==
                ak.firsts(dtrpairs["0"][var2][:,1:2], axis=1)
            ), False
        )
    )

    # if so, sort the taus with their pt
    sorted_idx = ak.where(where_same_pt_1,
                          ak.argsort(dtrpairs["1"][var2], ascending=False),
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

    return dtrpair, dtrpairidx



@selector(
    uses={
        "Muon.pt",
    },
    exposed=False,
)
def select_higgs_cand(
        self: Selector,
        lep_pairs: ak.Array,
        lep_pair_indices: ak.Array,
        channel: str,
        **kwargs
) -> ak.Array:
    sorted_lep_pair = None
    if channel == "etau" :
        sorted_dtrpair, sorted_dtrpairidx = sort_and_get_pair_semilep(lep_pairs,
                                                                      lep_pair_indices
                                                                      "pfRelIso03_all",
                                                                      "pt",
                                                                      "rawDeepTau2018v2p5VSjet",
                                                                      "pt")
    elif channel == "mutau" :
        sorted_dtrpair, sorted_dtrpairidx = sort_and_get_pair_semilep(lep_pairs, 
                                                                      lep_pair_indices,
                                                                      "pfRelIso03_all",
                                                                      "pt",
                                                                      "rawDeepTau2018v2p5VSjet",
                                                                      "pt")
    elif channel == "tautau" :
        sorted_dtrpair, sorted_dtrpairidx = sort_and_get_pair_fullhad(lep_pairs,
                                                                      lep_pair_indices, 
                                                                      "rawDeepTau2018v2p5VSjet",
                                                                      "pt")
        
    #return sorted_dtrpair, sorted_dtrpairidx
    return sorted_dtrpairidx
