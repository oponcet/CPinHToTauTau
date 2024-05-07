# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

from typing import Optional
from columnflow.production import Producer, producer
from columnflow.production.categories import category_ids
from columnflow.production.normalization import normalization_weights
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.production.cms.mc_weight import mc_weight
#from columnflow.production.cms.muon import muon_weights
from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column

from httcp.production.PhiCPNeutralPion import PhiCPNPMethod
from httcp.production.ReconstructPi0 import reconstructPi0

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

@producer(
    uses={
        # nano columns
        "Jet.pt",
    },
    produces={
        # new columns
        "ht", "n_jet",
    },
)
def features(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    events = set_ak_column(events, "ht", ak.sum(events.Jet.pt, axis=1))
    events = set_ak_column(events, "n_jet", ak.num(events.Jet.pt, axis=1), value_type=np.int32)

    return events


@producer(
    uses={
        "channel_id", 
        "hcand.pt", "hcand.eta", "hcand.phi", "hcand.mass", "hcand.decayMode",
        "hcandprod.pt", "hcandprod.eta", "hcandprod.phi", "hcandprod.mass","hcandprod.pdgId",
    },
)
def reArrangeDecayProducts(
        self: Producer,
        events: ak.Array,
        **kwargs
) :
    hcand      = events.hcand
    hcandprod  = events.hcandprod

    hcand1     = hcand[:, 0:1]
    hcand2     = hcand[:, 1:2]
    hcand1prod = ak.firsts(hcandprod[:,0:1], axis=1)
    hcand2prod = ak.firsts(hcandprod[:,1:2], axis=1)

    ispion_pos = lambda prod: ((prod.pdgId ==  211) | (prod.pdgId ==  321))
    ispion_neg = lambda prod: ((prod.pdgId == -211) | (prod.pdgId == -321))
    isphoton   = lambda prod: (prod.pdgId == 22)

    hcand1prod_pions   = hcand1prod[(ispion_pos(hcand1prod) | ispion_neg(hcand1prod))]
    hcand1prod_photons = hcand1prod[isphoton(hcand1prod)]
    hcand2prod_pions   = hcand2prod[(ispion_pos(hcand2prod) | ispion_neg(hcand2prod))]
    hcand2prod_photons = hcand2prod[isphoton(hcand2prod)]

    has_three_pions   = lambda prod : ak.sum((ispion_pos(prod) | ispion_neg(prod)), axis=1) == 3
    has_two_pos_one_neg_pions = lambda prod : has_three_pions(prod) & (ak.sum(ispion_pos(prod), axis=1) == 2) & (ak.sum(ispion_neg(prod), axis=1) == 1)
    has_two_neg_one_pos_pions = lambda prod : has_three_pions(prod) & (ak.sum(ispion_pos(prod), axis=1) == 1) & (ak.sum(ispion_neg(prod), axis=1) == 2)
    get_sorted_pion_indices   = lambda pions: ak.where(has_two_pos_one_neg_pions(pions),
                                                       ak.argsort(pions.pdgId, ascending=True),
                                                       ak.where(has_two_neg_one_pos_pions(pions),
                                                                ak.argsort(pions.pdgId, ascending=False),
                                                                ak.local_index(pions.pdgId)))
    
    hcand1pions_sorted_indices = get_sorted_pion_indices(hcand1prod_pions)
    hcand2pions_sorted_indices = get_sorted_pion_indices(hcand2prod_pions)

    #from IPython import embed; embed()
    #1/0
    
    hcand1prod_pions = hcand1prod_pions[hcand1pions_sorted_indices]
    hcand2prod_pions = hcand2prod_pions[hcand2pions_sorted_indices]
    
    # hcand1 and its decay products
    p4_hcand1     = ak.with_name(hcand1, "PtEtaPhiMLorentzVector")
    p4_hcand1_pi  = ak.with_name(hcand1prod_pions, "PtEtaPhiMLorentzVector")
    p4_hcand1_pi0 = reconstructPi0(hcand1prod_photons)

    # hcand2 and its decay products
    p4_hcand2     = ak.with_name(hcand2, "PtEtaPhiMLorentzVector")
    p4_hcand2_pi  = ak.with_name(hcand2prod_pions, "PtEtaPhiMLorentzVector")
    p4_hcand2_pi0 = reconstructPi0(hcand2prod_photons)

    return events, {"p4_hcand1"     : p4_hcand1, 
                    "p4_hcand1_pi"  : p4_hcand1_pi, 
                    "p4_hcand1_pi0" : p4_hcand1_pi0, 
                    "p4_hcand2"     : p4_hcand2, 
                    "p4_hcand2_pi"  : p4_hcand2_pi, 
                    "p4_hcand2_pi0" : p4_hcand2_pi0}

@producer(
    uses={
        "channel_id", 
        "hcand.pt", "hcand.eta", "hcand.phi", "hcand.mass", "hcand.decayMode",
        "hcandprod.pt", "hcandprod.eta", "hcandprod.phi", "hcandprod.mass","hcandprod.pdgId",
        reArrangeDecayProducts,
    },
    #produces={
    #    "phicp_NP",
    #},
)
def ProducePhiCP(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    events, P4_dict = self[reArrangeDecayProducts](events)

    is_rho_rho = ak.sum(events.hcand.decayMode == 1, axis=1) == 2
    is_a1_rho  = ak.sum((((events.hcand.decayMode[:,0:1]    == 10) 
                         & (events.hcand.decayMode[:,1:2]   == 1))
                        | ((events.hcand.decayMode[:,0:1]   == 1)
                           & (events.hcand.decayMode[:,1:2] == 10))), axis=1) == 1
    is_a1_a1   = ak.sum(events.hcand.decayMode == 10, axis=1) == 2

    #from IPython import embed; embed()
    #1/0

    dummy_phicp = P4_dict["p4_hcand1"].pt[:,:0]
    phicp_NP_obj = PhiCPNPMethod(P4_dict["p4_hcand1_pi"], P4_dict["p4_hcand1_pi0"], 
                                 P4_dict["p4_hcand2_pi"], P4_dict["p4_hcand2_pi0"])    
    #from IPython import embed; embed()
    phicp = ak.where(is_rho_rho, 
                     phicp_NP_obj.comp_PhiCPNP("rhorho", is_rho_rho),
                     ak.where(is_a1_rho,
                              phicp_NP_obj.comp_PhiCPNP("a1rho", is_a1_rho),
                              ak.where(is_a1_a1,
                                       phicp_NP_obj.comp_PhiCPNP("a1a1", is_a1_a1),
                                       dummy_phicp
                                   )
                          )
                 )
    #print(f"phicp: {phicp}")
    #from IPython import embed; embed()
    phicp = ak.enforce_type(phicp, "var * float64")
    #events = set_ak_column(events, "phicp_NP", phicp)
    
    return events, phicp


@producer(
    uses={
        # nano columns
        "hcand.*", "hcandprod.*", #reArrangeDecayProducts,
        ProducePhiCP,
    },
    produces={
        # new columns
        "hcand_invm", "hcand_dr", "phicp_NP",
    },
)
def hcand_features(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:

    events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
    #hcand_ = ak.with_name(ak.firsts(events.hcand, axis=1), "PtEtaPhiMLorentzVector")
    hcand_ = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
    hcand1 = hcand_[:,0:1]
    hcand2 = hcand_[:,1:2]
    
    mass = (hcand1 + hcand2).mass
    dr = ak.firsts(hcand1.metric_table(hcand2), axis=1)
    dr = ak.enforce_type(dr, "var * float32")

    events, phicp_NP = self[ProducePhiCP](events)

    #from IPython import embed; embed()
    #1/0
    events = set_ak_column(events, "hcand_invm", mass)
    events = set_ak_column(events, "hcand_dr", dr)
    events = set_ak_column(events, "phicp_NP", phicp_NP)
    
    return events



@producer(
    uses={
        mc_weight, category_ids,
        # nano columns
        "Jet.pt",
    },
    produces={
        mc_weight, category_ids,
        # new columns
        "cutflow.jet1_pt",
    },
)
def cutflow_features(
    self: Producer,
    events: ak.Array,
    object_masks: dict[str, dict[str, ak.Array]],
    **kwargs,
) -> ak.Array:
    if self.dataset_inst.is_mc:
        events = self[mc_weight](events, **kwargs)

    # apply object masks and create new collections
    reduced_events = create_collections_from_masks(events, object_masks)

    # create category ids per event and add categories back to the
    events = self[category_ids](reduced_events, target_events=events, **kwargs)

    # add cutflow columns
    events = set_ak_column(
        events,
        "cutflow.jet1_pt",
        Route("Jet.pt[:,0]").apply(events, EMPTY_FLOAT),
    )

    return events


@producer(
    uses={
        features, category_ids, normalization_weights, deterministic_seeds, 
        hcand_features, #ProducePhiCP, #muon_weights,
    },
    produces={
        features, category_ids, normalization_weights, deterministic_seeds, hcand_features, #ProducePhiCP,
        #muon_weights,
    },
)
def main(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    # features
    events = self[features](events, **kwargs)

    events = self[hcand_features](events, **kwargs)

    #from IPython import embed; embed()
    #events = self[ProducePhiCP](events, **kwargs)

    # category ids
    events = self[category_ids](events, **kwargs)

    # deterministic seeds
    events = self[deterministic_seeds](events, **kwargs)

    # mc-only weights
    if self.dataset_inst.is_mc:
        # normalization weights
        events = self[normalization_weights](events, **kwargs)

        # muon weights
        #events = self[muon_weights](events, **kwargs)
    
    return events
