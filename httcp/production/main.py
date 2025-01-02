# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

import law
from typing import Optional
from columnflow.production import Producer, producer
from columnflow.production.categories import category_ids
from columnflow.production.normalization import normalization_weights
from columnflow.production.normalization import stitched_normalization_weights

from columnflow.production.cms.pileup import pu_weight
from columnflow.production.cms.pdf import pdf_weights
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.production.cms.mc_weight import mc_weight

from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import

from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional

#from httcp.production.PhiCPNeutralPion import PhiCPNPMethod
from httcp.production.ReArrangeHcandProds import reArrangeDecayProducts, reArrangeGenDecayProducts
from httcp.production.PhiCP_Producer import ProduceDetPhiCP, ProduceGenPhiCP
#from httcp.production.weights import tauspinner_weight
from httcp.production.extra_weights import zpt_reweight, ff_weight # ff_weight : dummy
from httcp.production.muon_weights import muon_id_weights, muon_iso_weights, muon_trigger_weights, muon_xtrigger_weights
from httcp.production.electron_weights import electron_idiso_weights, electron_trigger_weights, electron_xtrigger_weights
from httcp.production.tau_weights import tau_weights# , tauspinner_weights


from httcp.production.dilepton_features import hcand_mass, mT, rel_charge #TODO: rename mutau_vars -> dilepton_vars
#from httcp.production.weights import pu_weight, muon_weight, tau_weight
#from httcp.production.weights import tau_weight
from httcp.production.sample_split import split_dy

from httcp.production.angular_features import ProduceDetCosPsi, ProduceGenCosPsi

from httcp.util import IF_DATASET_HAS_LHE_WEIGHTS, IF_DATASET_IS_DY, IF_DATASET_IS_W, IF_DATASET_IS_SIGNAL
from httcp.util import IF_RUN2, IF_RUN3, IF_ALLOW_STITCHING

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

logger = law.logger.get_logger(__name__)

@producer(
    uses={
        # nano columns
        "hcand.*", optional("GenTau.*"), optional("GenTauProd.*"),
        "Jet.pt",
        "PuppiMET.pt", "PuppiMET.phi",
        reArrangeDecayProducts,
        reArrangeGenDecayProducts,
        ProduceGenPhiCP,# ProduceGenCosPsi, 
        ProduceDetPhiCP, #ProduceDetCosPsi,
    },
    produces={
        # new columns
        "hcand_invm",
        "hcand_dr",
        "n_jet",
        ProduceGenPhiCP, #ProduceGenCosPsi,
        ProduceDetPhiCP, #ProduceDetCosPsi,
        "dphi_met_h1", "dphi_met_h2", "met_var_qcd_h1",
    },
)
def hcand_features(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:
    events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
    hcand_ = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
    hcand1 = hcand_[:,0:1]
    hcand2 = hcand_[:,1:2]

    met = ak.with_name(events.PuppiMET, "PtEtaPhiMLorentzVector")

    
    mass = (hcand1 + hcand2).mass
    #dr = ak.firsts(hcand1.metric_table(hcand2), axis=1)
    #dr = ak.enforce_type(dr, "var * float32")
    dr = hcand1.delta_r(hcand2)

    # deltaPhi between MET and hcand1 & 2
    dphi_met_h1 = met.delta_phi(hcand1)
    dphi_met_h2 = met.delta_phi(hcand2)
    met_var_qcd_h1 = met.pt * np.cos(dphi_met_h1)/hcand1.pt
    
    events = set_ak_column_f32(events, "hcand_invm",  mass)
    events = set_ak_column_f32(events, "hcand_dr",    dr)
    events = set_ak_column_f32(events, "dphi_met_h1", dphi_met_h1)
    events = set_ak_column_f32(events, "dphi_met_h2", dphi_met_h2)
    events = set_ak_column_f32(events, "met_var_qcd_h1", met_var_qcd_h1)
    
    events = set_ak_column_i32(events, "n_jet", ak.num(events.Jet.pt, axis=1))

    
   
    events, P4_dict     = self[reArrangeDecayProducts](events)
    events              = self[ProduceDetPhiCP](events, P4_dict)
    #from IPython import embed; embed()
    # events              = self[ProduceDetCosPsi](events, P4_dict)
    
    if self.config_inst.x.extra_tags.genmatch:
        if "is_signal" in list(self.dataset_inst.aux.keys()):
            #    if self.dataset_inst.aux["is_signal"]:
            events, P4_gen_dict = self[reArrangeGenDecayProducts](events)
            events = self[ProduceGenPhiCP](events, P4_gen_dict)
            #events = self[ProduceGenCosPsi](events, P4_gen_dict)
    
    return events


@producer(
    uses={
        ##deterministic_seeds,
        normalization_weights,
        IF_ALLOW_STITCHING(stitched_normalization_weights),
        #split_dy,
        pu_weight,
        IF_DATASET_HAS_LHE_WEIGHTS(pdf_weights),
        # -- muon -- #
        muon_id_weights,
        muon_iso_weights,
        muon_trigger_weights,
        muon_xtrigger_weights,
        # -- electron -- #
        electron_idiso_weights,
        electron_trigger_weights,
        electron_xtrigger_weights,
        # -- tau -- #
        tau_weights,
        # IF_DATASET_IS_SIGNAL(tauspinner_weights),
        IF_DATASET_IS_DY(zpt_reweight),
        hcand_features,
        hcand_mass,
        category_ids,
        #ff_weight,
    },
    produces={
        ##deterministic_seeds,
        normalization_weights,
        IF_ALLOW_STITCHING(stitched_normalization_weights),
        #split_dy,
        pu_weight,
        IF_DATASET_HAS_LHE_WEIGHTS(pdf_weights),
        # -- muon -- #
        muon_id_weights,
        muon_iso_weights,
        muon_trigger_weights,
        muon_xtrigger_weights,
        # -- electron -- #
        electron_idiso_weights,
        electron_trigger_weights,
        electron_xtrigger_weights,
        # -- tau -- #
        tau_weights,
        # IF_DATASET_IS_SIGNAL(tauspinner_weights),
        IF_DATASET_IS_DY(zpt_reweight),
        hcand_features,
        hcand_mass,
        "channel_id",
        "trigger_ids",
        category_ids,
        #ff_weight,
    },
)
def main(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    # deterministic seeds
    ##events = self[deterministic_seeds](events, **kwargs)
    ##events = self[category_ids](events, **kwargs)

    #events = self[ff_weight](events, **kwargs)
    if self.dataset_inst.is_mc:
        # allow stitching is applicable only when datasets are DY or wjets, only if the stitching booleans are true in config
        allow_stitching = bool(ak.any([(self.dataset_inst.has_tag("is_dy") and self.config_inst.x.allow_dy_stitching),
                                       (self.dataset_inst.has_tag("is_w") and self.config_inst.x.allow_w_stitching)]))
        if allow_stitching:
            events = self[stitched_normalization_weights](events, **kwargs)
        else:
            events = self[normalization_weights](events, **kwargs)

        events = self[pu_weight](events, **kwargs)
        #if not self.dataset_inst.has_tag("no_lhe_weights"):
        if self.has_dep(pdf_weights):
            events = self[pdf_weights](events, **kwargs)
        # ----------- Muon weights ----------- #
        events = self[muon_id_weights](events, **kwargs)
        events = self[muon_iso_weights](events, **kwargs)
        events = self[muon_trigger_weights](events, **kwargs)
        events = self[muon_xtrigger_weights](events, **kwargs)
        # ----------- Electron weights ----------- #
        events = self[electron_idiso_weights](events, **kwargs)
        events = self[electron_trigger_weights](events, **kwargs)
        events = self[electron_xtrigger_weights](events, **kwargs)
        # ----------- Tau weights ----------- #        
        events = self[tau_weights](events, do_syst=True, **kwargs)

        #from IPython import embed; embed()
        # if self.has_dep(tauspinner_weights):
        #     events = self[tauspinner_weights](events, **kwargs)

        if self.has_dep(zpt_reweight):
            events = self[zpt_reweight](events, **kwargs)

        #processes = self.dataset_inst.processes.names()
        #if ak.any(['dy_' in proc for proc in processes]):
        #    logger.info("splitting (any) Drell-Yan dataset ... ")
        #    events = self[split_dy](events,**kwargs)
            
    events = self[hcand_features](events, **kwargs)       

    # features
    # events = self[hcand_mass](events, **kwargs)
    # events = self[mT](events, **kwargs)

    return events
