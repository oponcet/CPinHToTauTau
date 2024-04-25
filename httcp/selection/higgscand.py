# coding: utf-8

"""
Prepare h-Candidate from SelectionResult: selected lepton indices & channel_id [trigger matched] 
"""

from typing import Optional
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
#maybe_import("coffea.nanoevents.methods.nanoaod")



@selector(
    uses={
        "channel_id",
    },
    produces={
        "hcand.pt", "hcand.eta", "hcand.phi", "hcand.mass", "hcand.charge", "hcand.decayModePNet",
    },
    exposed=False,
)
def higgscand(
        self: Selector,
        events: ak.Array,
        hcand_pair: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:
    #sel_hcand = ak.sum(ak.num(hcand_pair.pt, axis=-1), axis=1) > 0
    sel_hcand = ak.sum(ak.num(hcand_pair.pt, axis=-1), axis=1) == 2
    #hcand_col = ak.firsts(hcand_pair, axis=1)
    #print(ak.to_list(hcand_col.pt))

    # convert hcand to common p4
    #hcand_col = transform_hcand(hcand_col)

    #events = set_ak_column(events, "hcand", hcand_col)
    #events = set_ak_column(events, "hcand", hcand_pair)
    """
    empty_hcand_pair = hcand_pair[:,:0]
    hcand_pair_concat = ak.where(events.channel_id == 1, hcand_pair[:,0], empty_hcand_pair)
    hcand_pair_concat = ak.where(events.channel_id == 2, hcand_pair[:,1], hcand_pair_concat)
    hcand_pair_concat = ak.where(events.channel_id == 4, hcand_pair[:,2], hcand_pair_concat)
    """
    empty_hcand_pair = hcand_pair[:,:0][:,None]
    hcand_pair_concat = ak.where(events.channel_id == 1, hcand_pair[:,0][:,None], empty_hcand_pair)
    hcand_pair_concat = ak.where(events.channel_id == 2, hcand_pair[:,1][:,None], hcand_pair_concat)
    hcand_pair_concat = ak.where(events.channel_id == 4, hcand_pair[:,2][:,None], hcand_pair_concat)

    hcand_pair_concat = ak.where(events.channel_id == 3, 
                                 ak.concatenate([hcand_pair[:,0][:,None], hcand_pair[:,1][:,None]], axis=1),
                                 hcand_pair_concat)
    hcand_pair_concat = ak.where(events.channel_id == 5,
                                 ak.concatenate([hcand_pair[:,0][:,None], hcand_pair[:,2][:,None]], axis=1),
                                 hcand_pair_concat)
    hcand_pair_concat = ak.where(events.channel_id == 6, 
                                 ak.concatenate([hcand_pair[:,1][:,None], hcand_pair[:,2][:,None]], axis=1),
                                 hcand_pair_concat)
    
    #hcand_array = ak.zip({"pt": ak.values_astype(hcand_pair_concat.pt, "float32"),
    #                      "eta": ak.values_astype(hcand_pair_concat.eta, "float32"),
    #                      "phi": ak.values_astype(hcand_pair_concat.phi, "float32"),
    #                      "mass": ak.values_astype(hcand_pair_concat.mass, "float32"),
    #                      "charge": ak.values_astype(hcand_pair_concat.charge, "float32"),
    #                      "decayModePNet": ak.values_astype(hcand_pair_concat.decayModePNet, "float32")})

    #events = set_ak_column(events, "hcand", ak.zip())
    #from IPython import embed; embed()
    hcand_pair_concat = ak.Array(ak.to_list(hcand_pair_concat))
    events = set_ak_column(events, "hcand", hcand_pair_concat)
    #from IPython import embed; embed()

    sel_hcand = ak.fill_none(ak.num(ak.firsts(hcand_pair_concat.pt, axis=1), axis=1) == 2, False)
    #sel_hcand = ak.fill_none(ak.num(hcand_pair_concat.pt, axis=1) == 2, False)
    
    return events, SelectionResult(
        steps={
            #"atleast_one_higgs_cand_per_event": ak.num(ak.firsts(hcand_pair_concat.pt, axis=1), axis=1) == 2,
            "atleast_one_higgs_cand_per_event": sel_hcand,
        },
    )
