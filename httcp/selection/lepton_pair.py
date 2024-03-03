# coding: utf-8

"""
Lepton-Pairs Selection
https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2325&ancode=HIG-20-006&tp=an&line=HIG-20-006
http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_192_v15.pdf
"""

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import set_ak_column
from columnflow.util import maybe_import, DotDict

from httcp.selection.physics_object import muon_selection, electron_selection, tau_selection

from httcp.config.trigger_util import Trigger
from httcp.util import deltaR, new_invariant_mass

np = maybe_import("numpy")
ak = maybe_import("awkward")


@selector(
    uses={
        "Muon.pt", "Muon.eta", "Muon.phi", "Muon.charge", "Muon.mass",
        "Electron.pt", "Electron.eta", "Electron.phi", "Electron.charge", "Electron.mass",
        "Tau.pt", "Tau.eta", "Tau.phi", "Tau.charge", "Tau.mass",
        "MET.pt", "MET.phi",
    },
    exposed=False,
)
def get_presel_pairs(
        self: Selector,
        events: ak.Array,
        lep1_indices: ak.Array,
        lep2_indices: ak.Array,
        channel: int,
        **kwargs,
) -> tuple[SelectionResult, ak.Array, ak.Array]:    
    #leps1 = None
    #leps2 = None
    leps_pair        = None
    lep_indices_pair = None

    if channel == self.config_inst.get_channel("etau"):
        lep_indices_pair = ak.cartesian([lep1_indices, 
                                         lep2_indices], axis=1)
        leps_pair        = ak.cartesian([events.Electron[lep1_indices], 
                                         events.Tau[lep2_indices]], axis=1)
    elif channel == self.config_inst.get_channel("mutau"):
        lep_indices_pair = ak.cartesian([lep1_indices, 
                                         lep2_indices], axis=1)
        leps_pair        = ak.cartesian([events.Muon[lep1_indices], 
                                         events.Tau[lep2_indices]], axis=1)
    elif channel == self.config_inst.get_channel("tautau"):
        lep_indices_pair = ak.combinations(lep1_indices, 2, axis=-1)
        leps_pair        = ak.combinations(events.Tau[lep1_indices], 2, axis=-1)
    else:
        raise RuntimeError(f"Mention proper channel: WRONG {channel}")

    
    # pair of leptons: probable higgs candidate -> leps_pair
    # and their indices                         -> lep_indices_pair 
    lep1, lep2 = ak.unzip(leps_pair)
    lep1_idx, lep2_idx = ak.unzip(lep_indices_pair)

    preselection = {
        "is_os"         : (lep1.charge * lep2.charge) < 0,
        "dr_0p5"        : deltaR(lep1, lep2) > 0.5,
        "mT_50"         : transverse_mass(lep1, events.MET) < 50
    }

    good_pair_mask = ak.local_index(lep_indices_pair) >= 0
    pair_selection_steps = {}
    for cut in preselection.keys():
        good_pair_mask = good_pair_mask & preselection[cut]
        pair_selection_steps[cut] = ak.sum(preselection[cut], axis=1) > 0
        

    leps_pair_sel = leps_pair[good_pair_mask]
    lep_indices_pair_sel = lep_indices_pair[good_pair_mask]

    return SelectionResult(steps=pair_selection_steps), leps_pair_sel, lep_indices_pair_sel 



@selector(
    uses={
        electron_selection, muon_selection, tau_selection,
        # nano columns
        "event", 
        "Electron.charge", "Muon.charge", "Tau.charge", 
        "Electron.mass", "Muon.mass", "Tau.mass",
    },
    produces={
        electron_selection, muon_selection, tau_selection,
        # new columns
        "channel_id", "single_triggered", "cross_triggered",
    },
)
def lepton_pair_selection(
        self: Selector,
        events: ak.Array,
        trigger_results: SelectionResult,
        **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    """
    Combined lepton selection.
    """
    # get channels from the config
    ch_etau = self.config_inst.get_channel("etau")
    ch_mutau = self.config_inst.get_channel("mutau")
    ch_tautau = self.config_inst.get_channel("tautau")

    print(f"channels: {ch_etau}, {ch_mutau}, {ch_tautau}")
    
    # prepare vectors for output vectors
    false_mask = (abs(events.event) < 0)
    channel_id = np.uint8(1) * false_mask
    single_triggered = ak.copy(false_mask)
    cross_triggered = ak.copy(false_mask)
    empty_indices = ak.zeros_like(1 * events.event, dtype=np.uint16)[..., None][..., :0]

    sel_electron_indices = empty_indices
    sel_muon_indices = empty_indices
    sel_tau_indices = empty_indices

    # Initialize empty indices for higgs-candidate
    dummy_indices_obj = ak.zeros_like(
        ak.singletons(
            ak.firsts(events.Electrons, axis=1)
        )
    )
    # Initialize an empty h-cand array
    higgs_cand = empty_indices_obj[...,:0]

    # perform each lepton selection step separately per trigger, avoid caching
    sel_kwargs = {**kwargs, "call_force": True}
    sel_mask = {}
    for trigger, trigger_fired, leg_masks in trigger_results.x.trigger_data:
        print(f"trigger: {trigger}")
        print(f"trigger_fired: {trigger_fired}")
        is_single = trigger.has_tag("single_trigger")
        is_cross = trigger.has_tag("cross_trigger")

        print(f"Triggered? is_single: {is_single} :: is_cross: {is_cross} ")
        
        # electron selection
        electron_result, good_electron_indices, veto_electron_indices, dlveto_electron_indices = self[electron_selection](
            events,
            **sel_kwargs,
        )
        
        # muon selection
        muon_result, good_muon_indices, veto_muon_indices, dlveto_muon_indices = self[muon_selection](
            events,
            **sel_kwargs,
        )

        # tau selection
        tau_result, good_tau_indices = self[tau_selection](
            events,
            **sel_kwargs,
        )

        if trigger.has_tag({"single_e", "cross_e_tau"}): 
            # expect at least 1 electron,
            # 0 veto electron (loose ele, but not the one selected already as trigger matched, otherwise loose ele only),
            # 0 veto muons (loose muons),
            # and at least one tau
            is_etau = (
                trigger_fired
                & (ak.num(good_electron_indices, axis=1) >= 1)
                & (ak.num(veto_muon_indices, axis=1) == 0)
                & (ak.num(good_tau_indices, axis=1) >= 1)
            )
            
            where_etau = (channel_id == 0) & is_etau
            channel_id = ak.where(where_etau, ch_etau.id, channel_id)
            
            single_triggered = ak.where(where_etau & is_single, True, single_triggered)
            cross_triggered = ak.where(where_etau & is_cross, True, cross_triggered)

            sel_electron_indices = ak.where(where_etau, good_electron_indices, sel_electron_indices)
            sel_tau_indices = ak.where(where_etau, good_tau_indices, sel_tau_indices)   

                        
            etau_presel_results, etau_pairs, etau_index_pairs = get_presel_pairs(events, 
                                                                                 sel_electron_indices,
                                                                                 sel_tau_indices,
                                                                                 ch_etau)
            
            
            higgs_cand = ak.where(where_etau, 
                                  self[select_higgs_cand](etau_pairs,
                                                          etau_index_pairs,
                                                          "etau"), 
                                  higgs_cand)

            # extra single lepton veto


            sel_mask["etau_selection"] = is_etau



            # hcand selection will be called here

        elif trigger.has_tag({"single_mu", "cross_mu_tau"}): 
            # expect at least 1 muon,
            # 0 veto electron (loose ele, but not the one selected already as trigger matched, otherwise loose ele only),
            # 0 veto muons (loose muons),
            # and at least one tau
            is_mutau = (
                trigger_fired
                & (ak.num(good_muon_indices, axis=1) >= 1)
                & (ak.num(veto_electron_indices, axis=1) == 0)
                & (ak.num(good_tau_indices, axis=1) >= 1)
            )
            
            where_mutau = (channel_id == 0) & is_mutau
            channel_id = ak.where(where_mutau, ch_mutau.id, channel_id)
            
            single_triggered = ak.where(where_mutau & is_single, True, single_triggered)
            cross_triggered  = ak.where(where_mutau & is_cross, True, cross_triggered)

            sel_muon_indices = ak.where(where_mutau, good_muon_indices, sel_muon_indices)
            sel_tau_indices  = ak.where(where_mutau, good_tau_indices, sel_tau_indices)   


            mutau_presel_results, mutau_pairs, mutau_index_pairs = get_presel_pairs(events, 
                                                                                    sel_muon_indices,
                                                                                    sel_tau_indices,
                                                                                    ch_mutau)

            higgs_cand = ak.where(where_mutau, 
                                  self[select_higgs_cand](mutau_pairs,
                                                          mutau_index_pairs,
                                                          "mutau"), 
                                  higgs_cand)

            sel_mask["mutau_selection"] = is_mutau

            # hcand selection will be called here

        elif trigger.has_tag({"cross_tau_tau"}):
            # expect at least 2 taus,
            # 0 veto muon (loose muons)
            # 0 veto electrons (loose electrons),
            is_tautau = (
                trigger_fired &
                (ak.num(veto_electron_indices, axis=1) == 0) &
                (ak.num(veto_muon_indices, axis=1) == 0) &
                (ak.num(good_tau_indices, axis=1) >= 2)
            )

            where_tautau = (channel_id == 0) & is_tautau
            channel_id   = ak.where(where_tautau, ch_tautau.id, channel_id)
            
            single_triggered = ak.where(where_tautau & is_single, True, single_triggered)
            cross_triggered  = ak.where(where_tautau & is_cross, True, cross_triggered)

            sel_tau_indices  = ak.where(where_tautau, good_tau_indices, sel_tau_indices)

            tautau_presel_results, tautau_pairs, tautau_index_pairs = get_presel_pairs(events, 
                                                                                       sel_tau_indices,
                                                                                       sel_tau_indices,
                                                                                       ch_tautau)

            higgs_cand = ak.where(where_tautau, 
                                  self[select_higgs_cand](tautau_pairs,
                                                          tautau_index_pairs,
                                                          "tautau"), 
                                  higgs_cand)

            sel_mask["tautau_selection"] = is_tautau

            # hcand selection will be called here








        # some final type conversions
        channel_id           = ak.values_astype(channel_id, np.uint8)
        sel_electron_indices = ak.values_astype(sel_electron_indices, np.int32)
        sel_muon_indices     = ak.values_astype(sel_muon_indices, np.int32)
        sel_tau_indices      = ak.values_astype(sel_tau_indices, np.int32)

        # save new columns
        events = set_ak_column(events, "channel_id", channel_id)
        events = set_ak_column(events, "single_triggered", single_triggered)
        events = set_ak_column(events, "cross_triggered", cross_triggered)

        events = set_ak_column(events, "higgs_cand", higgs_cand)

        # Where should I apply the extra lepton veto? before hcand selection in the presence of more than one tau?
        return events, \
            sel_electron_indices, \
            sel_muon_indices, \
            sel_tau_indices, \
            Selectionresults(
                steps=sel_mask,
                objects={
                    
                }
            )
