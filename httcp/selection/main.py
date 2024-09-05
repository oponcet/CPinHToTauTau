# coding: utf-8

"""
Exemplary selection methods.
"""

from typing import Optional
from operator import and_
from functools import reduce
from collections import defaultdict, OrderedDict

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.stats import increment_stats
from columnflow.selection.cms.json_filter import json_filter
from columnflow.selection.cms.met_filters import met_filters

from columnflow.production.processes import process_ids
from columnflow.production.cms.mc_weight import mc_weight
from columnflow.production.util import attach_coffea_behavior
from columnflow.production.categories import category_ids

from columnflow.util import maybe_import
from columnflow.columnar_util import optional_column as optional
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column

from httcp.selection.physics_objects import *
from httcp.selection.trigger import trigger_selection
from httcp.selection.lepton_pair_etau import etau_selection
from httcp.selection.lepton_pair_mutau import mutau_selection
from httcp.selection.lepton_pair_tautau import tautau_selection
from httcp.selection.event_category import get_categories
#from httcp.selection.match_trigobj import match_trigobj
from httcp.selection.trigobject_matching import match_trigobj
from httcp.selection.lepton_veto import *
from httcp.selection.higgscand import higgscand, higgscandprod

#from httcp.production.main import cutflow_features
from httcp.production.weights import scale_mc_weight
from httcp.production.dilepton_features import hcand_mass, mT, rel_charge #TODO: rename mutau_vars -> dilepton_vars

from httcp.util import filter_by_triggers, get_objs_p4, trigger_object_matching_deep

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

@selector(uses={"process_id", optional("mc_weight")}) #TODO: Move it to utils
def custom_increment_stats(
    self: Selector,
    events: ak.Array,
    results: SelectionResult,
    stats: dict,
    **kwargs,
) -> ak.Array:
    """
    Unexposed selector that does not actually select objects but instead increments selection
    *stats* in-place based on all input *events* and the final selection *mask*.
    """
    # get event masks
    event_mask = results.event

    # get a list of unique process ids present in the chunk
    unique_process_ids = np.unique(events.process_id)
    # increment plain counts
    n_evt_per_file = self.dataset_inst.n_events/self.dataset_inst.n_files
    stats["num_events"] = n_evt_per_file
    stats["num_events_selected"] += float(ak.sum(event_mask, axis=0))
    #print(f"stats : {stats}")
    #from IPython import embed; embed()
    if self.dataset_inst.is_mc:
        stats[f"sum_mc_weight"] = n_evt_per_file
        stats.setdefault(f"sum_mc_weight_per_process", defaultdict(float))
        for p in unique_process_ids:
            stats[f"sum_mc_weight_per_process"][int(p)] = n_evt_per_file
        
    # create a map of entry names to (weight, mask) pairs that will be written to stats
    weight_map = OrderedDict()
    if self.dataset_inst.is_mc:
        # mc weight for selected events
        weight_map["mc_weight_selected"] = (events.mc_weight, event_mask)

    # get and store the sum of weights in the stats dictionary
    for name, (weights, mask) in weight_map.items():
        joinable_mask = True if mask is Ellipsis else mask

        # sum of different weights in weight_map for all processes
        stats[f"sum_{name}"] += float(ak.sum(weights[mask]))
        # sums per process id
        stats.setdefault(f"sum_{name}_per_process", defaultdict(float))
        for p in unique_process_ids:
            stats[f"sum_{name}_per_process"][int(p)] += float(ak.sum(
                weights[(events.process_id == p) & joinable_mask],
            ))

    return events, results


def get_2n_pairs(etau_indices_pair,
                 mutau_indices_pair,
                 tautau_indices_pair):

    has_one_etau_pair   = ak.num(etau_indices_pair, axis=1) == 2
    has_one_mutau_pair  = ak.num(mutau_indices_pair, axis=1) == 2
    has_one_tautau_pair = ak.num(tautau_indices_pair, axis=1) == 2

    has_at_least_one_pair = (has_one_etau_pair | has_one_mutau_pair | has_one_tautau_pair)

    return has_at_least_one_pair    

"""
def get_event_level_eff(events, results):
    from tabulate import tabulate
    steps_ = results.steps.keys()
    indiv_selections_ = []
    comb_selections_ = []
    indiv_headers_ = ["selections", "nevents", "abs eff"]
    comb_headers_ = ["selections", "nevents", "abs eff", "rel eff"]
    init = len(events)
    comb_count_array = ak.Array(np.ones(init, dtype=bool))
    for step_ in steps_:
        count_array      = results.steps[step_]
        comb_count_array = comb_count_array & count_array
        count            = ak.sum(count_array)
        comb_count       = ak.sum(comb_count_array)
        indiv_selections_.append([step_, count, round(count/init,3)])
        comb_selections_.append([step_, comb_count, round(comb_count/init,3)])
    indiv_table_ = tabulate(indiv_selections_,
                            indiv_headers_,
                            tablefmt="pretty")
    print(f"\n---> Efficiencies of individual selections: \n{indiv_table_}")

    comb_selections_ = np.array(comb_selections_)
    comb_selections_counts_ = comb_selections_[:,1]
    comb_den_ = np.array([init] + comb_selections_counts_[:-1].tolist())
    rel_eff_ = np.round(np.asarray(comb_selections_counts_, float)/np.asarray(comb_den_, float), decimals=3)
    comb_selections_ = np.concatenate([comb_selections_, rel_eff_[:,None]], axis=1).tolist()
    comb_table_ = tabulate(comb_selections_,
                           comb_headers_,
                           tablefmt="pretty")
        
    print(f"\n---> Efficiencies of combined selections: \n{comb_table_}")
    

def get_object_eff(results, tag):
    from tabulate import tabulate
    print(f"\n{tag}")
    keys = [key for key in results.aux.keys() if key.startswith(f"{tag}_")]
    rows = []
    rows_evt_level = []

    n0 = 0
    nevt0 = 0

    for i, key in enumerate(keys):
        mask = results.aux[key]
        n = ak.sum(ak.sum(mask, axis=1))
        nevt = ak.sum(ak.any(mask, axis=1))
        if i == 0:
            n0 = n
            nevt0 = nevt
        rows.append([key, n, round(n/n0, 3)])
        rows_evt_level.append([key, nevt, round(nevt/nevt0, 3)])
    table = tabulate(rows, ["selection", f"n_{tag}", "abseff"], tablefmt="pretty")
    evt_table = tabulate(rows_evt_level, ["selection", f"n_{tag}", "abseff"], tablefmt="pretty")
    
    print(f"object level : \n{table}")
    print(f"event level  : \n{evt_table}")    
"""

    
# exposed selectors
# (those that can be invoked from the command line)
@selector(
    uses={
        "event",
        # selectors / producers called within _this_ selector
        attach_coffea_behavior,
        json_filter, 
        met_filters, 
        scale_mc_weight, 
        process_ids,
        trigger_selection, 
        muon_selection, 
        electron_selection, 
        tau_selection, 
        jet_selection,
        etau_selection, 
        mutau_selection, 
        tautau_selection, 
        get_categories,
        extra_lepton_veto, 
        double_lepton_veto, 
        match_trigobj,
        increment_stats, 
        custom_increment_stats,
        higgscand,
        gentau_selection,
        higgscandprod,
        rel_charge,
        category_ids,
    },
    produces={
        # selectors / producers whose newly created columns should be kept
        scale_mc_weight, 
        trigger_selection, 
        muon_selection, 
        electron_selection, 
        tau_selection, 
        jet_selection,
        etau_selection, 
        mutau_selection, 
        tautau_selection, 
        get_categories, 
        process_ids,
        extra_lepton_veto, 
        double_lepton_veto, 
        match_trigobj,
        higgscandprod,
        gentau_selection,
        rel_charge,
        category_ids,
        "trigger_ids",
    },
    exposed=True,
)
def main(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    
    # ensure coffea behaviors are loaded
    events = self[attach_coffea_behavior](events, **kwargs)

    # prepare the selection results that are updated at every step
    results = SelectionResult()

    # --------- Sel : Good Lumi JSON for DATA only -------- #
    # filter bad data events according to golden lumi mask
    if self.dataset_inst.is_data:
        events, json_filter_results = self[json_filter](events, **kwargs)
        results += json_filter_results
    else:
        results += SelectionResult(steps={"json": np.ones(len(events), dtype=bool)})


    # -------- Sel : Met Filters -------- #
    # met filter selection
    events, met_filter_results = self[met_filters](events, **kwargs)
    results += met_filter_results


    # --------- Sel : Trigger Selection -------- #
    # trigger selection
    events, trigger_results = self[trigger_selection](events, **kwargs)
    results += trigger_results

    
    # -------- Sel : b-veto -------- #
    # jet selection
    events, bjet_veto_result = self[jet_selection](events, 
                                                   call_force=True, 
                                                   **kwargs)
    results += bjet_veto_result


    # electron selection
    # e.g. ele_idx: [ [], [0,1], [], [], [1,2] ] 
    events, ele_results, good_ele_indices, veto_ele_indices, dlveto_ele_indices = self[electron_selection](events,
                                                                                                           call_force=True,
                                                                                                           **kwargs)
    results += ele_results

    
    # muon selection
    # e.g. mu_idx: [ [0,1], [], [1], [0], [] ] 
    events, muon_results, good_muon_indices, veto_muon_indices, dlveto_muon_indices = self[muon_selection](events,
                                                                                                           call_force=True,
                                                                                                           **kwargs)
    results += muon_results


    # tau selection
    # e.g. tau_idx: [ [1], [0,1], [1,2], [], [0,1] ] 
    events, tau_results, good_tau_indices = self[tau_selection](events, call_force=True, **kwargs)
    results += tau_results

    # double lepton veto
    events, extra_double_lepton_veto_results = self[double_lepton_veto](events,
                                                                        dlveto_ele_indices,
                                                                        dlveto_muon_indices)
    results += extra_double_lepton_veto_results
    
    #from IPython import embed; embed()

    
    # check if there are at least two leptons with at least one tau
    _lepton_indices = ak.concatenate([good_muon_indices, good_ele_indices, good_tau_indices], axis=1)
    nlep_mask = ((ak.num(_lepton_indices, axis=1) >= 2) & (ak.num(good_tau_indices, axis=1) >= 1))

    match_nlep = SelectionResult(
        steps = {
            "has_2_or_more_leps_with_at_least_1_tau"  : nlep_mask,
        },
    )
    results += match_nlep
    
    
    # e-tau pair i.e. hcand selection
    # e.g. [ [], [e1, tau1], [], [], [e1, tau2] ]
    etau_results, etau_indices_pair = self[etau_selection](events,
                                                           good_ele_indices,
                                                           good_tau_indices,
                                                           call_force=True,
                                                           **kwargs)
    results += etau_results

    # mu-tau pair i.e. hcand selection
    # e.g. [ [mu1, tau1], [], [mu1, tau2], [], [] ]
    mutau_results, mutau_indices_pair = self[mutau_selection](events,
                                                              good_muon_indices,
                                                              good_tau_indices,
                                                              call_force=True,
                                                              **kwargs)
    results += mutau_results

    # tau-tau pair i.e. hcand selection
    # e.g. [ [], [tau1, tau2], [], [], [] ]
    tautau_results, tautau_indices_pair = self[tautau_selection](events,
                                                                 good_tau_indices,
                                                                 call_force=True,
                                                                 **kwargs)
    results += tautau_results

    
    pre_match_at_least_one_pair = get_2n_pairs(etau_indices_pair,
                                               mutau_indices_pair,
                                               tautau_indices_pair)


    pre_match_pair_result = SelectionResult(
        steps = {
            "has_at_least_1_pair_before_trigobj_matching"  : pre_match_at_least_one_pair,
        },
    )
    results += pre_match_pair_result

    #from IPython import embed; embed()
    #for i in range(1000):
    #  print(etau_indices_pair[i], mutau_indices_pair[i], tautau_indices_pair[i], pre_match_at_least_one_pair[i])
    
    # Trigger objects matching
    events, matchedResults = self[match_trigobj](events,
                                                 trigger_results,
                                                 etau_indices_pair,
                                                 mutau_indices_pair,
                                                 tautau_indices_pair,
                                                 True)
    results += matchedResults

    etau_indices_pair    = matchedResults.x.etau["pair_indices"]
    mutau_indices_pair   = matchedResults.x.mutau["pair_indices"]
    tautau_indices_pair  = matchedResults.x.tautau["pair_indices"]

    
    post_match_at_least_one_pair = get_2n_pairs(etau_indices_pair,
                                                mutau_indices_pair,
                                                tautau_indices_pair)
    post_match_pair_result = SelectionResult(
        steps = {
            "has_at_least_1_pair_after_trigobj_matching"  : post_match_at_least_one_pair,
        },
    )
    results += post_match_pair_result


    #from IPython import embed; embed()
    
    etau_pair    = ak.concatenate([events.Electron[etau_indices_pair[:,0:1]],
                                   events.Tau[etau_indices_pair[:,1:2]]],
                                  axis=1)
    
    mutau_pair   = ak.concatenate([events.Muon[mutau_indices_pair[:,0:1]], 
                                   events.Tau[mutau_indices_pair[:,1:2]]],
                                  axis=1)

    tautau_pair  = ak.concatenate([events.Tau[tautau_indices_pair[:,0:1]], 
                                   events.Tau[tautau_indices_pair[:,1:2]]], 
                                  axis=1)

    # channel selection
    # channel_id is now in columns
    events, channel_results = self[get_categories](events,
                                                   trigger_results,
                                                   etau_indices_pair,
                                                   mutau_indices_pair,
                                                   tautau_indices_pair)
    results += channel_results

    #etau_trigger_names = matchedResults.x.etau["trigger_names"]
    etau_trigger_ids   = matchedResults.x.etau["trigger_ids"]

    #mutau_trigger_names = matchedResults.x.mutau["trigger_names"]
    mutau_trigger_ids   = matchedResults.x.mutau["trigger_ids"]
    
    #tautau_trigger_names = matchedResults.x.tautau["trigger_names"]
    tautau_trigger_ids   = matchedResults.x.tautau["trigger_ids"]
    
    #trigger_names = ak.concatenate([etau_trigger_names[:,None], mutau_trigger_names[:,None], tautau_trigger_names[:,None]], axis=1)
    #trigger_ids   = ak.concatenate([etau_trigger_ids[:,None], mutau_trigger_ids[:,None], tautau_trigger_ids[:,None]], axis=1)
    trigger_ids = ak.concatenate([etau_trigger_ids, mutau_trigger_ids, tautau_trigger_ids], axis=1)
    trigger_ids = ak.values_astype(ak.fill_none(ak.firsts(trigger_ids, axis=1), 999), np.uint64)
    #from IPython import embed; embed()
    # save the trigger_ids column
    events = set_ak_column(events, "trigger_ids", trigger_ids)
    # hcand pair: [ [[mu1,tau1]], [[e1,tau1],[tau1,tau2]], [[mu1,tau2]], [], [[e1,tau2]] ]
    hcand_pairs = ak.concatenate([etau_pair[:,None], mutau_pair[:,None], tautau_pair[:,None]], axis=1)

    
    # extra lepton veto
    # it is only applied on the events with one higgs candidate only
    events, extra_lepton_veto_results = self[extra_lepton_veto](events,
                                                                veto_ele_indices,
                                                                veto_muon_indices,
                                                                hcand_pairs)
    results += extra_lepton_veto_results
    

    #from IPython import embed; embed()
    # hcand results
    events, hcand_array, hcand_results = self[higgscand](events, hcand_pairs)
    results += hcand_results
    
    # hcand prod results
    events, hcandprod_results = self[higgscandprod](events, hcand_array)
    results += hcandprod_results

    # gen particles info
    # hcand-gentau match = True/False
    """
    if "is_signal" in list(self.dataset_inst.aux.keys()):
        if self.dataset_inst.aux["is_signal"]:
            #print("hcand-gentau matching")
            events, gentau_results = self[gentau_selection](events, True)
            results += gentau_results
    """
    #from IPython import embed; embed()
    #1/0

    # create process ids
    #events = self[process_ids](events, **kwargs)

    # combined event selection after all steps
    event_sel = reduce(and_, results.steps.values())
    results.event = event_sel
    
    # add the mc weight
    if self.dataset_inst.is_mc:
        events = self[scale_mc_weight](events, **kwargs)

    # rel-charge
    events = self[rel_charge](events, **kwargs)

    # add cutflow features, passing per-object masks
    #events = self[cutflow_features](events, results.objects, **kwargs)

    events = self[process_ids](events, **kwargs)
    events = self[category_ids](events, **kwargs)     

    events, results = self[custom_increment_stats]( 
        events,
        results,
        stats,
    )

    """
    # inspect cuts
    if self.config_inst.x.verbose.selection.main:
        print(f"\n---> Inspecting event selections : ")
        get_event_level_eff(events, results)
        
        print(f"\n---> Inspecting trigger selections: ")

        from tabulate import tabulate        
        # trigger details
        trigger_names = results.aux["trigger_names"]
        trigger_ids   = results.aux["trigger_ids"]

        HLT_names = [trigger.name for trigger in self.config_inst.x.triggers]
        HLT_ids   = [trigger.id for trigger in self.config_inst.x.triggers]

        trig_info = []
        for i in range(len(HLT_ids)):
            HLT = HLT_names[i]
            ID  = HLT_ids[i]
            nPassed = ak.sum(ak.any(trigger_ids == ID, axis=1))
            trig_info.append([HLT, ID, nPassed, round(nPassed/len(events), 3)])
        trig_table = tabulate(trig_info, ["HLT", "ID", "nEvents Passed", "Efficiency"], tablefmt="pretty")
        print(trig_table)

        print(f"\n---> Inspecting object selections: ")
        
        # muon
        get_object_eff(results, "muon")
        get_object_eff(results, "electron")
        get_object_eff(results, "tau")
        
        print(f"\n---> Inspecting pair selections: ")
        
        # pairs
        print(f"\nBefore trigobj matching : ")
        get_object_eff(results, "etau")
        get_object_eff(results, "mutau")
        get_object_eff(results, "tautau")
        
        print(f"\nAfter trigobj matching: ")
        pairs = []
        pairs.append(["etau", ak.sum(ak.num(results.aux["etau"]["pair_indices"], axis=1) == 2)])
        pairs.append(["mutau", ak.sum(ak.num(results.aux["mutau"]["pair_indices"], axis=1) == 2)])
        pairs.append(["tautau", ak.sum(ak.num(results.aux["tautau"]["pair_indices"], axis=1) == 2)])
        pair_table = tabulate(pairs, ["pair", "nEvents with pair"], tablefmt="pretty")
        
        print(pair_table)

        print(f"\nCategorization: ")
        cats = []
        cats.append(["is_etau", ak.sum(results.aux["cat_is_etau"])])
        cats.append(["is_mutau", ak.sum(results.aux["cat_is_mutau"])])
        cats.append(["is_tautau", ak.sum(results.aux["cat_is_tautau"])])
        cats.append(["is_etau_mutau", ak.sum(results.aux["cat_is_etau_mutau"])])
        cats.append(["is_etau_tautau", ak.sum(results.aux["cat_is_etau_tautau"])])
        cats.append(["is_mutau_tautau", ak.sum(results.aux["cat_is_mutau_tautau"])])
        cats_table = tabulate(cats, ["category", "nEvents"], tablefmt="pretty")

        print(cats_table)

        print(f"\nEvents selected per channel: ")
        sel_ev = ak.sum(events.channel_id > 0)
        print(f"nSelectedEvents : {sel_ev}")
        channels = []
        etau_ev = ak.sum(events.channel_id == 1)
        mtau_ev = ak.sum(events.channel_id == 2)
        ttau_ev = ak.sum(events.channel_id == 4)
        mixed   = ak.sum(~((events.channel_id == 0)
                           | (events.channel_id == 1)
                           | (events.channel_id == 2)
                           | (events.channel_id == 4)))
        channels.append(["etau", etau_ev, round(etau_ev/sel_ev, 3)])
        channels.append(["mutau", mtau_ev, round(mtau_ev/sel_ev, 3)])
        channels.append(["tautau", ttau_ev, round(ttau_ev/sel_ev, 3)])
        channels.append(["other", mixed, round(mixed/sel_ev, 3)])
        channel_table = tabulate(channels, ["channel", "nEvents", "eff"], tablefmt="pretty")
        
        print(channel_table)
        print(f"Total selected events in etau, mutau and tautau chennels : {etau_ev+mtau_ev+ttau_ev}\n\n")
    """

    return events, results
