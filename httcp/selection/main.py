# coding: utf-8

"""
Exemplary selection methods.
"""

from operator import and_
from functools import reduce
from collections import defaultdict

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.stats import increment_stats
from columnflow.selection.cms.json_filter import json_filter
from columnflow.selection.cms.met_filters import met_filters

from columnflow.production.processes import process_ids
from columnflow.production.cms.mc_weight import mc_weight

from columnflow.util import maybe_import

from httcp.production.example import cutflow_features

from httcp.selection.physics_objects import *
from httcp.selection.trigger import trigger_selection
from httcp.selection.lepton_pair_etau import etau_selection
from httcp.selection.lepton_pair_mutau import mutau_selection
from httcp.selection.lepton_pair_tautau import tautau_selection
from httcp.selection.event_category import get_categories

np = maybe_import("numpy")
ak = maybe_import("awkward")


# exposed selectors
# (those that can be invoked from the command line)
@selector(
    uses={
        # selectors / producers called within _this_ selector
        json_filter, met_filters, mc_weight, cutflow_features, process_ids,
        trigger_selection, muon_selection, electron_selection, tau_selection, jet_selection,
        etau_selection, mutau_selection, tautau_selection,
        increment_stats,
    },
    produces={
        # selectors / producers whose newly created columns should be kept
        mc_weight, trigger_selection, cutflow_features, process_ids,
    },
    exposed=True,
)
def main(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:

    # prepare the selection results that are updated at every step
    results = SelectionResult()

    # filter bad data events according to golden lumi mask
    if self.dataset_inst.is_data:
        events, json_filter_results = self[json_filter](events, **kwargs)
        results += json_filter_results

    # trigger selection
    events, trigger_results = self[trigger_selection](events, **kwargs)
    results += trigger_results

    # met filter selection
    events, met_filter_results = self[met_filters](events, **kwargs)
    results += met_filter_results

    # double lepton veto
    events, extra_double_lepton_veto_results = self[double_lepton_veto](events,
                                                                        dlveto_electron_indices,
                                                                        dlveto_muon_indices)
    results += extra_double_lepton_veto_results

    # jet selection
    events, bjet_veto_result = self[jet_selection](events, 
                                                   call_force=True, 
                                                   **kwargs)
    results += bjet_veto_results

    # muon selection
    # e.g. mu_idx: [ [0,1], [], [1], [0], [] ] 
    events, muon_results, good_muon_indices, veto_muon_indices, dlveto_muon_indices = self[muon_selection](events,
                                                                                                           call_force=True, 
                                                                                                           **kwargs)
    results += muon_results

    # electron selection
    # e.g. ele_idx: [ [], [0,1], [], [], [1,2] ] 
    events, ele_results, good_ele_indices, veto_ele_indices, dlveto_ele_indices = self[electron_selection](events,
                                                                                                           call_force=True, 
                                                                                                           **kwargs)
    results += ele_results

    # tau selection
    # e.g. tau_idx: [ [1], [0,1], [1,2], [], [0,1] ] 
    events, tau_results, good_tau_indices = self[tau_selection](events,
                                                                good_ele_indices,
                                                                good_muon_indices,
                                                                call_force=True, 
                                                                **kwargs)
    results += tau_results

    # e-tau pair i.e. hcand selection
    # e.g. [ [], [e1, tau1], [], [], [e1, tau2] ]
    events, etau_results, etau_pair_lep1_idx, etau_pair_lep2_idx = self[etau_selection](events,
                                                                                        good_ele_indices,
                                                                                        good_tau_indices,
                                                                                        call_force=True,
                                                                                        **kwargs)
    results += etau_results
    etau_indices_pair = ak.concatenate([etau_pair_lep1_idx, etau_pair_lep2_idx], axis=1)
    etau_pair         = ak.concatenate([events.Electron[etau_pair_lep1_idx], 
                                        events.Tau[etau_pair_lep2_idx]], axis=1)

    # mu-tau pair i.e. hcand selection
    # e.g. [ [mu1, tau1], [], [mu1, tau2], [], [] ]
    events, mutau_results, mutau_pair_lep1_idx, mutau_pair_lep2_idx = self[mutau_selection](events,
                                                                                            good_muon_indices,
                                                                                            good_tau_indices,
                                                                                            call_force=True,
                                                                                            **kwargs)
    results += mutau_results
    mutau_indices_pair = ak.concatenate([mutau_pair_lep1_idx, mutau_pair_lep2_idx], axis=1)
    mutau_pair = ak.concatenate([events.Muon[mutau_pair_lep1_idx], 
                                 events.Tau[mutau_pair_lep2_idx]], axis=1)

    # tau-tau pair i.e. hcand selection
    # e.g. [ [], [tau1, tau2], [], [], [] ]
    events, tautau_results, tautau_pair_lep1_idx, tautau_pair_lep2_idx = self[tautau_selection](events,
                                                                                                good_tau_indices,
                                                                                                call_force=True,
                                                                                                **kwargs)
    results += tautau_results
    tautau_indices_pair = ak.concatenate([tautau_pair_lep1_idx, tautau_pair_lep2_idx], axis=1)
    tautau_pair = ak.concatenate([events.Tau[mutau_pair_lep1_idx], 
                                  events.Tau[mutau_pair_lep2_idx]], axis=1)
    
    
    # make sure events have at least one lepton pair
    # hcand pair: [ [[mu1,tau1]], [[e1,tau1],[tau1,tau2]], [[mu1,tau2]], [], [[e1,tau2]] ]
    hcand_pair = ak.concatenate([etau_pair[:,None], mutau_pair[:,None], tautau_pair[:,None]], axis=1)
    count_higgs_cand_pair = ak.concatenate([etau_pair_lep1_idx[:,None], 
                                            mutau_pair_lep1_idx[:,None], 
                                            tautau_pair_lep1_idx[:,None]], axis=1)
    hcand_results = SelectionResult(
        steps={
            "Atleast_one_higgs_cand": ak.sum(ak.num(count_higgs_cand_pair, axis=-1), axis=1) > 0,
        },
    )
    results += hcand_results


    # channel selection
    # channel_id is now in columns
    events, channel_results = self[get_categories](events, trigger_results, 
                                                   etau_pair_indices, mutau_pair_indices, tautau_pair_indices)
    results += channel_results

    # extra lepton veto
    events, extra_lepton_veto_results = self[extra_lepton_veto](events, 
                                                                veto_ele_indices,
                                                                veto_muon_indices,
                                                                hcand_pair)
    results += extra_lepton_veto_results


    # create process ids
    events = self[process_ids](events, **kwargs)

    # add the mc weight
    if self.dataset_inst.is_mc:
        events = self[mc_weight](events, **kwargs)

    # add cutflow features, passing per-object masks
    events = self[cutflow_features](events, results.objects, **kwargs)

    # increment stats
    weight_map = {
        "num_events": Ellipsis,
        "num_events_selected": results.event,
    }
    group_map = {}
    if self.dataset_inst.is_mc:
        weight_map = {
            **weight_map,
            # mc weight for all events
            "sum_mc_weight": (events.mc_weight, Ellipsis),
            "sum_mc_weight_selected": (events.mc_weight, results.event),
        }
        group_map = {
            # per process
            "process": {
                "values": events.process_id,
                "mask_fn": (lambda v: events.process_id == v),
            },
            # per jet multiplicity
            "njet": {
                "values": results.x.n_jets,
                "mask_fn": (lambda v: results.x.n_jets == v),
            },
        }
    events, results = self[increment_stats](
        events,
        results,
        stats,
        weight_map=weight_map,
        group_map=group_map,
        **kwargs,
    )

    return events, results
