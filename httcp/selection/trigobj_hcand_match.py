# coding: utf-8

"""
Exemplary selection methods.
"""

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import

from httcp.selection.physics_objects import *
from httcp.selection.trigger import trigger_selection
from httcp.selection.lepton_pair_etau import etau_pair_selection
from httcp.selection.lepton_pair_mutau import mutau_pair_selection
from httcp.selection.lepton_pair_tautau import tautau_pair_selection

np = maybe_import("numpy")
ak = maybe_import("awkward")


@selector(
    uses={
        "Channel_id", "trigger_ids",
    },
    exposed=False
)
def trigobj_lepton_match(
        self: Selector,
        events: ak.Array,
        trigger_result: SelectionResult,
        hcand_pair: ak.Array,
        **kwargs
) -> ak.Array:
    trigger_fired_and_all_legs_match = trigger_result.x.trigger_fired_and_all_legs_match
    # hcand pair: e.g. [ [[mu1,tau1]], 
    #                    [[e1,tau1],[tau1,tau2]], 
    #                    [[mu1,tau2]], 
    #                    [[tau1,tau2]], 
    #                    [[e1,tau2]] ]
    # Trigger_ids: e.g. [ [101,201], [201], [502], [103,301], [501,503] ]
    has_single_pair = ak.sum(ak.num(ak.local_index(hcand_pair), axis=-1), axis=1) == 2
    has_multi_pairs = ak.sum(ak.num(ak.local_index(hcand_pair), axis=-1), axis=1) > 2

    # concat
    # e.g. [ [([mu1,tau1],101),([mu1,tau1],201)], [([e1,tau1],201),([tau1,tau2],201)], [], [], [] ] 
    pair_trigger_ids_hcand = ak.cartesian([hcand_pair, events.trigger_ids])
    # unzip
    # e.g. 
    # hcand : [ [[mu1,tau1],[mu1,tau1]], [[e1,tau1],[tau1,tau2]], [], [], [] ]
    # trigid: [ [101,201],               [201,201],               [], [], [] ]
    hcand, trigid = ak.unzip(pair_trigger_ids_hcand)

    # HLTs         : [HLT1, HLT2, ...., HLT6]
    # HLTIDs       : [201,  203, ....., 501 ]
    # HLT_leg_mask : []
    HLTs = ak.Array([trig_tuple[0] for trig_tuple in trigger_result])
    # broadcast HLTs
    HLTs_brdcst, trigid_brdcst = ak.broadcast_arrays(HLTs[:,...][None][None], trigid[:,:,None])
    is_trig = HLTs_brdcst.id == trigid_brdcst
    HLTs_with_fired_id = HLTs_brdcst[is_trig]



    return events, SelectionResults(
        steps = {
            "trigger_fired_all_legs_match": trigger_fired_and_all_legs_match,
        },
    )
