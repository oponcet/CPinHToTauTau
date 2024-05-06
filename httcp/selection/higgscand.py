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


def enforce_hcand_type(hcand_pair_concat, field_type_dict):
    temp = {}
    for field, typename in field_type_dict.items():
        temp[field] = ak.enforce_type(ak.values_astype(hcand_pair_concat[field], typename), f"var * var * {typename}")
    hcand_array = ak.zip(temp)
    return hcand_array
    


@selector(
    uses={
        "channel_id",
    },
    #produces={
    #},
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
    #                      "decayMode": ak.values_astype(hcand_pair_concat.decayMode, "float32")})
    """
    hcand_array = ak.zip(
        {
            "pt"           : ak.enforce_type(ak.values_astype(hcand_pair_concat.pt,    "float64"), "var * var * float64"),
            "eta"          : ak.enforce_type(ak.values_astype(hcand_pair_concat.eta,   "float64"), "var * var * float64"),
            "phi"          : ak.enforce_type(ak.values_astype(hcand_pair_concat.phi,   "float64"), "var * var * float64"),
            "mass"         : ak.enforce_type(ak.values_astype(hcand_pair_concat.mass,  "float64"), "var * var * float64"),
            "charge"       : ak.enforce_type(ak.values_astype(hcand_pair_concat.charge,        "int64"), "var * var * int64"),
            "decayMode": ak.enforce_type(ak.values_astype(hcand_pair_concat.decayMode, "int64"), "var * var * int64"),
            "genPartFlav"  : ak.enforce_type(ak.values_astype(hcand_pair_concat.genPartFlav,   "int64"), "var * var * int64"),
            "genPartIdx"   : ak.enforce_type(ak.values_astype(hcand_pair_concat.genPartIdx,    "int64"), "var * var * int64"),
        }
    )
    """
    hcand_array = enforce_hcand_type(hcand_pair_concat, 
                                     {"pt"            : "float64",
                                      "eta"           : "float64",
                                      "phi"           : "float64",
                                      "mass"          : "float64",
                                      "charge"        : "int64",
                                      "decayMode" : "int64",
                                      "genPartFlav"   : "int64",
                                      "genPartIdx"    : "int64",
                                      "rawIdx"        : "int64"}
    )

    sel_hcand = ak.fill_none(ak.num(ak.firsts(hcand_array.pt, axis=1), axis=1) == 2, False)

    #events = set_ak_column(events, "hcand", hcand_array)
    
    return events, hcand_array, SelectionResult(
        steps={
            #"atleast_one_higgs_cand_per_event": ak.num(ak.firsts(hcand_pair_concat.pt, axis=1), axis=1) == 2,
            "One_higgs_cand_per_event": sel_hcand,
        },
    )


def select_tauprods(hcand_idx, tauprods):
    hcand_idx_brdcst, tauprod_tauIdx = ak.broadcast_arrays(ak.firsts(hcand_idx,axis=1), tauprods.tauIdx)
    #hcand_idx_brdcst, tauprod_tauIdx = ak.broadcast_arrays(hcand_idx, tauprods.tauIdx)
    hcandprod_mask                   = tauprod_tauIdx == hcand_idx_brdcst
    hcandprods                       = tauprods[hcandprod_mask]

    return hcandprods

is_pion   = lambda prods : ((np.abs(prods.pdgId) == 211) | (np.abs(prods.pdgId) == 321))
is_photon = lambda prods : prods.pdgId == 22

has_one_pion    = lambda prods : (ak.sum(is_pion(prods),   axis = 1) == 1)[:,None]
has_three_pions = lambda prods : (ak.sum(is_pion(prods),   axis = 1) == 3)[:,None]
has_photons     = lambda prods : (ak.sum(is_photon(prods), axis = 1) >  0)[:,None]
has_no_photons  = lambda prods : (ak.sum(is_photon(prods), axis = 1) == 0)[:,None]

@selector(
    uses={
        "channel_id", "TauProd.*",#"hcand.rawIdx",
    },
    produces={
        "hcand.pt", "hcand.eta", "hcand.phi", "hcand.mass", "hcand.charge", "hcand.rawIdx",
        "hcand.decayMode", "hcand.genPartFlav", "hcand.genPartIdx",
        "hcandprod.pt", "hcandprod.eta", "hcandprod.phi", "hcandprod.mass", 
        "hcandprod.charge", "hcandprod.pdgId", "hcandprod.tauIdx",
    },
    exposed=False,
)
def higgscandprod(
        self: Selector,
        events: ak.Array,
        hcand_array: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:
    tauprods = events.TauProd
    hcand    = hcand_array #events.hcand
    hcand1 = ak.firsts(hcand[:,:,0:1], axis=1)
    hcand2 = ak.firsts(hcand[:,:,1:2], axis=1)
    #hcand1 = ak.firsts(hcand[:,0:1], axis=1)
    #hcand2 = ak.firsts(hcand[:,1:2], axis=1)
    hcand_concat = ak.concatenate([hcand1, hcand2], axis=1)
    
    
    hcand1_idx = hcand1.rawIdx
    hcand2_idx = hcand2.rawIdx

    hcand1prods = ak.where(events.channel_id == 4,
                           select_tauprods(hcand1_idx, tauprods), 
                           tauprods[:,:0])
    hcand2prods = ak.where(((events.channel_id == 1) | (events.channel_id == 2) | (events.channel_id == 4)), 
                           select_tauprods(hcand2_idx, tauprods),
                           tauprods[:,:0])

    #dummy = ak.full_like(hcand1.pt, False)
    dummy = (events.event >= 0)[:,None]
    hcand1_mask = ak.where(hcand1.decayMode == 0,
                           has_one_pion(hcand1prods),
                           ak.where(((hcand1.decayMode == 1) | (hcand1.decayMode == 2)),
                                    (has_one_pion(hcand1prods) & has_photons(hcand1prods)),
                                    ak.where(hcand1.decayMode == 10,
                                             has_three_pions(hcand1prods),
                                             ak.where(hcand1.decayMode == 11,
                                                      (has_three_pions(hcand1prods) & has_photons(hcand1prods)),
                                                      dummy)
                                         )
                                )
                       )
    hcand2_mask = ak.where(hcand2.decayMode == 0,
                           has_one_pion(hcand2prods),
                           ak.where(((hcand2.decayMode == 1) | (hcand2.decayMode == 2)),
                                    (has_one_pion(hcand2prods) & has_photons(hcand2prods)),
                                    ak.where(hcand2.decayMode == 10,
                                             has_three_pions(hcand2prods),
                                             ak.where(hcand2.decayMode == 11,
                                                      (has_three_pions(hcand2prods) & has_photons(hcand2prods)),
                                                      dummy)
                                         )
                                )
                       )
    
    hcand_prod_mask = ak.concatenate([hcand1_mask, hcand2_mask], axis=1)
    
    hcand_prods = ak.concatenate([hcand1prods[:,None], hcand2prods[:,None]], axis=1)

    hcand_prods_array = enforce_hcand_type(ak.from_regular(hcand_prods),
                                           {"pt"            : "float64",
                                            "eta"           : "float64",
                                            "phi"           : "float64",
                                            "mass"          : "float64",
                                            "charge"        : "int64",
                                            "pdgId"         : "int64",
                                            "tauIdx"        : "int64"}
                                       )

    events = set_ak_column(events, "hcand",     hcand_concat)
    events = set_ak_column(events, "hcandprod", hcand_prods_array)

    #for i in range(500): 
    #    print(f"ch : {events.channel_id[i]}\t{hcand1.decayMode[i]}\t{hcand2.decayMode[i]}\t{hcand1_mask[i]}\t{hcand2_mask[i]}\t{hcand1prods.pdgId[i]}\t{hcand2prods.pdgId[i]}")

    #from IPython import embed; embed()
    #1/0

    return events, SelectionResult(
        steps={
            "has proper tau decay products": ak.sum(hcand_prod_mask, axis=1) == 2,
        },
    )
