"""
Exemplary selection methods.
"""

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.util import maybe_import

np = maybe_import("numpy")
ak = maybe_import("awkward")


@selector(
    #uses={#
    #
    #   },
    produces={
        "channel_id",
    },
    exposed=False,
)
def get_categories(
        self: Selector,
        events: ak.Array,
        trigger_results: SelectionResult,
        etau_pair_indices: ak.Array,
        mutau_pair_indices: ak.Array,
        tautau_pair_indices: ak.Array,
        #FFDR_tautau_pair_indices: ak.Array,
        FFDRIso_tautau_pair_indices: ak.Array,
        FFDRantiIso_tautau_pair_indices: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:
    # get channels from the config
    ch_etau   = self.config_inst.get_channel("etau")
    ch_mutau  = self.config_inst.get_channel("mutau")
    ch_tautau = self.config_inst.get_channel("tautau")
    #ch_FFDR_tautau = self.config_inst.get_channel("FFDR_tautau")
    ch_FFDRIso_tautau = self.config_inst.get_channel("FFDRIso_tautau")
    ch_FFDRantiIso_tautau = self.config_inst.get_channel("FFDRantiIso_tautau")

    false_mask       = (abs(events.event) < 0)
    single_triggered = false_mask
    cross_triggered  = false_mask
    empty_indices    = ak.zeros_like(1 * events.event, dtype=np.uint16)[..., None][..., :0]

    channel_selections = {
        "cat_is_etau"           : [ch_etau.id, # 1
                                   ((ak.num(etau_pair_indices, axis=1) == 2) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 0) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 0))],
        "cat_is_mutau"          : [ch_mutau.id,  # 2
                                   ((ak.num(etau_pair_indices, axis=1) == 0) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 2) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 0))],
        "cat_is_tautau"         : [ch_tautau.id, # 4
                                   ((ak.num(etau_pair_indices, axis=1) == 0) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 0) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 2))],
        "cat_is_etau_mutau"     : [ch_etau.id + ch_mutau.id,  # 3
                                   ((ak.num(etau_pair_indices, axis=1) == 2) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 2) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 0))],
        "cat_is_etau_tautau"    : [ch_etau.id + ch_tautau.id,  # 5
                                   ((ak.num(etau_pair_indices, axis=1) == 2) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 0) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 2))], 
        "cat_is_mutau_tautau"   : [ch_mutau.id + ch_tautau.id, # 6
                                   ((ak.num(etau_pair_indices, axis=1) == 0) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 2) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 2))], 
        # "cat_is_FFDR_tautau"    : [ch_FFDR_tautau.id, # 8
        #                               ((ak.num(etau_pair_indices, axis=1) == 0) 
        #                                 & (ak.num(mutau_pair_indices, axis=1) == 0) 
        #                                 & (ak.num(FFDR_tautau_pair_indices, axis=1) == 2))],
        "cat_is_FFDRIso_tautau"    : [ch_FFDRIso_tautau.id, # 9
                                      ((ak.num(etau_pair_indices, axis=1) == 0) 
                                        & (ak.num(mutau_pair_indices, axis=1) == 0) 
                                        & (ak.num(FFDRIso_tautau_pair_indices, axis=1) == 2))], 
        "cat_is_FFDRantiIso_tautau" : [ch_FFDRantiIso_tautau.id, # 10
                                      ((ak.num(etau_pair_indices, axis=1) == 0) 
                                        & (ak.num(mutau_pair_indices, axis=1) == 0) 
                                        & (ak.num(FFDRantiIso_tautau_pair_indices, axis=1) == 2))],                               
    }

    selection_steps  = {}
    channel_id       = np.uint8(1) * false_mask
    for key, val in channel_selections.items():
        selection_steps[key] = val[1]
        channel_id = ak.where(val[1], val[0], channel_id)
        
    # some final type conversions
    channel_id = ak.values_astype(channel_id, np.uint8)
    events = set_ak_column(events, "channel_id", channel_id)

    return events, SelectionResult(aux=selection_steps)
