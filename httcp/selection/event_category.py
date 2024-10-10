"""
Exemplary selection methods.
"""

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.util import maybe_import

from httcp.util import transverse_mass
from httcp.util import IF_RUN2, IF_RUN3

np = maybe_import("numpy")
ak = maybe_import("awkward")


@selector(
    #uses={#
    #
    #   },
    produces={
        "channel_id",
        #"is_os", "is_iso1", "is_iso2",
        #"is_real1", "is_real2",
    },
    exposed=False,
)
def get_categories(
        self: Selector,
        events: ak.Array,
        results: SelectionResult,
        etau_pair_indices: ak.Array,
        mutau_pair_indices: ak.Array,
        tautau_pair_indices: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:
    # get channels from the config
    ch_etau   = self.config_inst.get_channel("etau")
    ch_mutau  = self.config_inst.get_channel("mutau")
    ch_tautau = self.config_inst.get_channel("tautau")

    false_mask       = (abs(events.event) < 0)

    channel_selections = {
        "cat_is_etau"           : [ch_etau.id, 
                                   ((ak.num(etau_pair_indices, axis=1) == 2) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 0) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 0))],
        "cat_is_mutau"          : [ch_mutau.id, 
                                   ((ak.num(etau_pair_indices, axis=1) == 0) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 2) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 0))],
        "cat_is_tautau"         : [ch_tautau.id, 
                                   ((ak.num(etau_pair_indices, axis=1) == 0) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 0) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 2))],
        "cat_is_etau_mutau"     : [ch_etau.id + ch_mutau.id, 
                                   ((ak.num(etau_pair_indices, axis=1) == 2) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 2) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 0))],
        "cat_is_etau_tautau"    : [ch_etau.id + ch_tautau.id, 
                                   ((ak.num(etau_pair_indices, axis=1) == 2) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 0) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 2))], 
        "cat_is_mutau_tautau"   : [ch_mutau.id + ch_tautau.id, 
                                   ((ak.num(etau_pair_indices, axis=1) == 0) 
                                    & (ak.num(mutau_pair_indices, axis=1) == 2) 
                                    & (ak.num(tautau_pair_indices, axis=1) == 2))], 
    }

    selection_steps  = {}
    channel_id       = np.uint8(1) * false_mask
    for key, val in channel_selections.items():
        selection_steps[key] = val[1]
        channel_id = ak.where(val[1], val[0], channel_id)
        
    # some final type conversions
    channel_id = ak.values_astype(channel_id, np.uint8)
    events = set_ak_column(events, "channel_id", channel_id)


    """
    # categoris for SR, AR and DR
    cat_tautau_is_OS   = results.x.cat_tautau_is_OS
    cat_tautau_is_ISO1 = results.x.cat_tautau_is_ISO1
    cat_tautau_is_ISO2 = results.x.cat_tautau_is_ISO2

    is_OS   = ak.where(events.channel_id == ch_tautau.id, cat_tautau_is_OS, False)
    is_ISO1 = ak.where(events.channel_id == ch_tautau.id, cat_tautau_is_ISO1, False)
    is_ISO2 = ak.where(events.channel_id == ch_tautau.id, cat_tautau_is_ISO2, False)

    events = set_ak_column(events, "is_os", is_OS)
    events = set_ak_column(events, "is_iso1", is_ISO1)
    events = set_ak_column(events, "is_iso2", is_ISO2)

    if self.dataset_inst.is_mc:
        cat_tautau_is_REAL1 = results.x.cat_tautau_is_REAL1
        cat_tautau_is_REAL2 = results.x.cat_tautau_is_REAL2

        is_REAL1 = ak.where(events.channel_id == ch_tautau.id, cat_tautau_is_REAL1, False)
        is_REAL2 = ak.where(events.channel_id == ch_tautau.id, cat_tautau_is_REAL2, False)

        events = set_ak_column(events, "is_real1", is_REAL1)
        events = set_ak_column(events, "is_real2", is_REAL2)
    
    """

    return events, SelectionResult(
        #steps={
        #    "category_et_mt_or_tt": ((channel_id == ch_etau.id) | (channel_id == ch_mutau.id) | (channel_id == ch_tautau.id)),
        #},
        aux=selection_steps)





@selector(
    uses={
        "channel_id", "hcand.*",
        IF_RUN2("MET.pt", "MET.phi"),
        IF_RUN3("PuppiMET.pt", "PuppiMET.phi"),
    },
    produces={
        "is_os",
        "is_b_veto",
        "is_low_mt",
        "is_iso_1", "is_iso_2",
        #"is_real_1", "is_real_2",
    },
    exposed=False,
)
def build_abcd_masks(
        self: Selector,
        events: ak.Array,
        bjet_veto_mask: ak.Array,
        **kwargs
) -> ak.Array:

    # get channels from the config
    ch_etau   = self.config_inst.get_channel("etau")
    ch_mutau  = self.config_inst.get_channel("mutau")
    ch_tautau = self.config_inst.get_channel("tautau")

    hcand = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
    h1 = hcand[:,0:1]
    h2 = hcand[:,1:2]


    #from IPython import embed; embed()
    
    
    # tau tagger wp
    tau_tagger      = self.config_inst.x.deep_tau_tagger
    tau_tagger_info = self.config_inst.x.deep_tau_info[tau_tagger]
    #wp_tag = self.config_inst.x.deep_tau_info[tau_tagger].vs_j
    #tau_tagger_wps  = self.config_inst.x.deep_tau_info[tau_tagger].wp
    vs_jet_wp       = lambda tau_tagger_info, ch : tau_tagger_info.wp.vs_j[tau_tagger_info.vs_j[ch]]

    # OS --> opposite sign
    is_os = (h1.charge * h2.charge) < 0
    is_os = ak.fill_none(ak.any(is_os, axis=1), False)

    #from IPython import embed; embed()
    

    # ISO1 --> required for tau-tau channel only to categorise events on the basis of
    # leading tau isolation
    is_iso_1_dummy = h1.rawIdx < 0
    is_iso_1 = ak.where(events.channel_id == ch_tautau.id,
                        h1.idVsJet >= vs_jet_wp(tau_tagger_info, ch_tautau.name),
                        is_iso_1_dummy)
    is_iso_1 = ak.fill_none(ak.any(is_iso_1, axis=1), False)

    # ISO2 --> required for alla channels
    id_etau_pass = h2.idVsJet >= vs_jet_wp(tau_tagger_info, "etau")
    id_mutau_pass = h2.idVsJet >= vs_jet_wp(tau_tagger_info, "mutau")
    id_tautau_pass = h2.idVsJet >= vs_jet_wp(tau_tagger_info, "tautau")
    
    is_iso_2 = ak.where(events.channel_id == ch_tautau.id,
                        id_tautau_pass,
                        ak.where(events.channel_id == ch_mutau.id,
                                 id_mutau_pass,
                                 ak.where(events.channel_id == ch_etau.id,
                                          id_etau_pass,
                                          is_iso_1_dummy)
                                 )
                        )
    is_iso_2 = ak.fill_none(ak.any(is_iso_2, axis=1), False)

    # REAL1 --> to get the contribution of real MC taus only with genPartFlav > 0
    # only required for tau-tau channel
    #is_real_1 = ak.where(events.channel_id == ch_tautau.id, h1.genPartFlav > 0, is_iso_1_dummy)
    #is_real_1 = ak.fill_none(ak.any(is_real_1, axis=1), False)

    # REAL2 --> the same for the sub-leading tau
    # valid for all channels
    #is_real_2 = h2.genPartFlav > 0
    #is_real_2 = ak.fill_none(ak.any(is_real_2, axis=1), False)

    # BVETO --> events with / without bjets
    is_b_veto = bjet_veto_mask

    # LOWMT --> Required for leptonic channels
    met = events.MET if self.config_inst.campaign.x.run == 2 else events.PuppiMET
    is_low_mt = transverse_mass(h1, met) < 50
    is_low_mt = ak.fill_none(ak.any(is_low_mt, axis=1), False)


    # set columns
    events = set_ak_column(events, "is_os",     is_os)
    events = set_ak_column(events, "is_iso_1",  is_iso_1)
    events = set_ak_column(events, "is_iso_2",  is_iso_2)
    #events = set_ak_column(events, "is_real_1", is_real_1)
    #events = set_ak_column(events, "is_real_2", is_real_2)
    events = set_ak_column(events, "is_low_mt", is_low_mt)
    events = set_ak_column(events, "is_b_veto", is_b_veto)

    #from IPython import embed; embed()

    return events
