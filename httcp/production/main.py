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
from httcp.production.ReArrangeHcandProds import reArrangeDecayProducts, reArrangeGenDecayProducts
#from httcp.production.ReconstructPi0 import reconstructPi0

from httcp.production.dilepton_features import hcand_mass, mT, rel_charge #TODO: rename mutau_vars -> dilepton_vars
from httcp.production.weights import pu_weight, muon_weight, tau_weight
from httcp.production.sample_split import split_dy

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)


@producer(
    uses={
        "channel_id", 
        "hcand.pt", "hcand.eta", "hcand.phi", "hcand.mass", "hcand.decayMode",
        "hcandprod.pt", "hcandprod.eta", "hcandprod.phi", "hcandprod.mass","hcandprod.pdgId",
        #"GenTau.*", "GenPart.*",
        reArrangeDecayProducts, reArrangeGenDecayProducts,
    },
)
def ProducePhiCP(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    events, P4_dict = self[reArrangeDecayProducts](events)

    events, P4_gen_dict = self[reArrangeGenDecayProducts](events)
    #from IPython import embed; embed()
    #1/0

    is_rho_rho = ak.sum(events.hcand.decayMode == 1, axis=1) == 2
    is_a1_rho  = ak.sum((((events.hcand.decayMode[:,0:1]    == 10) 
                         & (events.hcand.decayMode[:,1:2]   == 1))
                        | ((events.hcand.decayMode[:,0:1]   == 1)
                           & (events.hcand.decayMode[:,1:2] == 10))), axis=1) == 1
    is_a1_a1   = ak.sum(events.hcand.decayMode == 10, axis=1) == 2

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
    #from IPython import embed; embed()
    phicp = ak.enforce_type(phicp, "var * float64")
    
    return events, # phicp


@producer(
    uses={
        # nano columns
        "hcand.*", #"hcandprod.*", #reArrangeDecayProducts,
        ProducePhiCP, #"GenTau.*", "GenPart.*",
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

    events = set_ak_column(events, "hcand_invm", mass)
    events = set_ak_column(events, "hcand_dr",   dr)
    events = set_ak_column(events, "phicp_NP",   phicp_NP)
    
    return events

@producer(
    uses={
        rel_charge,
        category_ids,
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        tau_weight,
        hcand_mass,
    },
    produces={
        rel_charge,
        category_ids,
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        tau_weight,
        hcand_mass,
    },
)
def main(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    events = self[rel_charge](events, **kwargs)
    
    events = self[category_ids](events, **kwargs) 
    if self.dataset_inst.is_mc:
        events = self[normalization_weights](events, **kwargs)
        processes = self.dataset_inst.processes.names()
        if ak.any(['dy' in proc for proc in processes]):
            print("Splitting Drell-Yan dataset...")
            events = self[split_dy](events, **kwargs)
        events = self[pu_weight](events, **kwargs)
        events = self[muon_weight](events, **kwargs)
        events = self[tau_weight](events, **kwargs) 
    #events = self[hcand_features](events, **kwargs)       
    # features
    events = self[hcand_mass](events, **kwargs)
   # events = self[mT](events, **kwargs)
    return events
