# coding: utf-8

"""
Definition of triggers
"""

import order as od

from httcp.config.trigger_util import Trigger, TriggerLeg



# ----------------------------------------------- #
#                   Run2 2016 UL                  #
# ----------------------------------------------- #
def add_triggers_UL2016(config: od.Config, postfix: str) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    
    PreVFP:
    +----------------------------------------------------+----+---------------+---------------+---------------+---------------+---------------+
    |                        HLT                         | MC | /Tau_Run2016B | /Tau_Run2016C | /Tau_Run2016D | /Tau_Run2016E | /Tau_Run2016F |
    +----------------------------------------------------+----+---------------+---------------+---------------+---------------+---------------+
    |            HLT_Ele25_eta2p1_WPTight_Gsf            | √  |       √       |       √       |       √       |       √       |       √       |
    |                    HLT_IsoMu22                     | √  |       √       |       √       |       √       |       √       |       √       |
    |                   HLT_IsoTkMu22                    | √  |       √       |       √       |       √       |       √       |       √       |
    |                HLT_IsoTkMu22_eta2p1                | √  |       √       |       √       |       √       |       √       |       √       |
    |         HLT_IsoMu19_eta2p1_LooseIsoPFTau20         | √  |       √       |       √       |       √       |       √       |       √       |
    |    HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1     | √  |       √       |       √       |       √       |       √       |       √       |
    |     HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg     | √  |       √       |       √       |       √       |       √       |       √       |
    | HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg | √  |       X       |       X       |       X       |       X       |       X       |
    +----------------------------------------------------+----+---------------+---------------+---------------+---------------+---------------+
    PostVFP:
    +----------------------------------------------------+----+---------------+---------------+---------------+
    |                        HLT                         | MC | /Tau_Run2016F | /Tau_Run2016G | /Tau_Run2016H |
    +----------------------------------------------------+----+---------------+---------------+---------------+
    |            HLT_Ele25_eta2p1_WPTight_Gsf            | √  |       √       |       √       |       √       |
    |                    HLT_IsoMu22                     | √  |       √       |       √       |       √       |
    |                   HLT_IsoTkMu22                    | √  |       √       |       √       |       √       |
    |                HLT_IsoTkMu22_eta2p1                | √  |       √       |       √       |       √       |
    |         HLT_IsoMu19_eta2p1_LooseIsoPFTau20         | √  |       √       |       √       |       √       |
    |    HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1     | √  |       √       |       √       |       √       |
    |     HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg     | √  |       √       |       √       |       X       |
    | HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg | √  |       X       |       X       |       √       |
    +----------------------------------------------------+----+---------------+---------------+---------------+

    """
    config.x.triggers = od.UniqueObjectIndex(Trigger, [
        # ===>>> single electron
        Trigger(
            name="HLT_Ele25_eta2p1_WPTight_Gsf",
            id=111,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=26.0,
                    # filter names:
                    # hltEle25erWPTightGsfTrackIsoFilter
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        # ===>>> single muon
        Trigger(
            name="HLT_IsoMu22",
            id=131,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=23.0,
                    # filter names:
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu22_eta2p1",
            id=132,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=23.0,
                    # filter names:
                    # hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoTkMu22",
            id=133,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=23.0,
                    # filter names:
                    # hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoTkMu22_eta2p1",
            id=134,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=23.0,
                    # filter names:
                    # hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        # ===>>> e-tauh

        # ===>>> mu-tauh
        Trigger(
            name="HLT_IsoMu19_eta2p1_LooseIsoPFTau20",
            id=13151,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=20.0,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=21.0,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=1024 + 512,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1",
            id=13152,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=20.0,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=21.0,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=1024 + 512,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        # ===>>> tauh-tauh
        Trigger(
            name="HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg",
            id=15151,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=36.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=36.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.is_mc) if postfix == "preVFP" else (lambda dataset_inst: dataset_inst.x.era < "H"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg",
            id=15152,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=36.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=36.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc) if postfix == "preVFP" else (lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era == "H"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
    ])    



# ----------------------------------------------- #
#                   Run2 2017 UL                  #
# ----------------------------------------------- #
def add_triggers_UL2017(config: od.Config) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.

    +--------------------------------------------------------------------+----+---------------+---------------+---------------+---------------+---------------+
    |                                HLT                                 | MC | /Tau_Run2017B | /Tau_Run2017C | /Tau_Run2017D | /Tau_Run2017E | /Tau_Run2017F |
    +--------------------------------------------------------------------+----+---------------+---------------+---------------+---------------+---------------+
    |                       HLT_Ele27_WPTight_Gsf                        | √  |       √       |       √       |       √       |       √       |       √       |
    |                       HLT_Ele32_WPTight_Gsf                        | √  |       X       |       X       |       √       |       √       |       √       |
    |                            HLT_IsoMu24                             | √  |       √       |       √       |       √       |       √       |       √       |
    |                            HLT_IsoMu27                             | √  |       √       |       √       |       √       |       √       |       √       |
    | HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1 | √  |       √       |       √       |       √       |       √       |       √       |
    |      HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1      | √  |       √       |       √       |       √       |       √       |       √       |
    |          HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg          | √  |       √       |       √       |       √       |       √       |       √       |
    |     HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg      | √  |       √       |       √       |       √       |       √       |       √       |
    |      HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg      | √  |       √       |       √       |       √       |       √       |       √       |
    +--------------------------------------------------------------------+----+---------------+---------------+---------------+---------------+---------------+

    """
    config.x.triggers = od.UniqueObjectIndex(Trigger, [
        # ===>>> single electron
        Trigger(
            name="HLT_Ele27_WPTight_Gsf",
            id=111,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=28.0,
                    # filter names:
                    # hltEle32L1DoubleEGWPTightGsfTrackIsoFilter
                    # hltEGL1SingleEGOrFilter
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        Trigger(
            name="HLT_Ele32_WPTight_Gsf",
            id=112,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=28.0,
                    # filter names:
                    # hltEle32WPTightGsfTrackIsoFilter
                    trigger_bits=2 + 1024,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era >= "D"),
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        # ===>>> single muon
        Trigger(
            name="HLT_IsoMu24",
            id=131,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    # filter names:
                    # hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu27",
            id=132,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    # filter names:
                    # hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        # ===>>> e-tauh
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1",
            id=11151,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=25.0,
                    # filter names:
                    # hltEle24erWPTightGsfTrackIsoFilterForTau
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    # filter names:
                    # hltSelectedPFTau30LooseChargedIsolationL1HLTMatched
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=1024 + 256,
                ),
            ],
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
        ),
        # ===>>> mu-tauh
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1",
            id=13151,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=21.0,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=1024 + 512,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        # ===>>> tauh-tauh
        Trigger(
            name="HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg",
            id=15151,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
            id=15152,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            #applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
            id=15153,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=37.0,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=37.0,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            #applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
    ])    



# ----------------------------------------------- #
#                   Run2 2018 UL                  #
# ----------------------------------------------- #
def add_triggers_UL2018(config: od.Config) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.

    +-----------------------------------------------------------------------+----+---------------+---------------+---------------+---------------+
    |                                  HLT                                  | MC | /Tau_Run2018A | /Tau_Run2018B | /Tau_Run2018C | /Tau_Run2018D |
    +-----------------------------------------------------------------------+----+---------------+---------------+---------------+---------------+
    |                         HLT_Ele35_WPTight_Gsf                         | √  |       √       |       √       |       √       |       √       |
    |                         HLT_Ele32_WPTight_Gsf                         | √  |       √       |       √       |       √       |       √       |
    |                              HLT_IsoMu24                              | √  |       √       |       √       |       √       |       √       |
    |                              HLT_IsoMu27                              | √  |       √       |       √       |       √       |       √       |
    |   HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1  | X  |       √       |       √       |       X       |       X       |
    | HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1 | √  |       X       |       √       |       √       |       √       |
    |        HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1       | X  |       √       |       √       |       X       |       X       |
    |      HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1      | √  |       √       |       √       |       √       |       √       |
    |    HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1   | X  |       √       |       √       |       X       |       X       |
    |            HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg           | X  |       √       |       √       |       X       |       X       |
    |        HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg       | X  |       √       |       √       |       X       |       X       |
    |        HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg       | X  |       √       |       √       |       X       |       X       |
    |          HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg         | √  |       √       |       √       |       √       |       √       |
    +-----------------------------------------------------------------------+----+---------------+---------------+---------------+---------------+

    """
    config.x.triggers = od.UniqueObjectIndex(Trigger, [
        # ===>>> single electron
        Trigger(
            name="HLT_Ele35_WPTight_Gsf",
            id=111,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=36.0,
                    # filter names:
                    # hltEle32L1DoubleEGWPTightGsfTrackIsoFilter
                    # hltEGL1SingleEGOrFilter
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        Trigger(
            name="HLT_Ele32_WPTight_Gsf",
            id=112,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=28.0,
                    # filter names:
                    # hltEle32WPTightGsfTrackIsoFilter
                    trigger_bits=2 + 1024,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),

        # ===>>> single muon
        Trigger(
            name="HLT_IsoMu24",
            id=131,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    # filter names:
                    # hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu27",
            id=132,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=28.0,
                    # filter names:
                    # hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),

        # ===>>> e-tauh
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1",
            id=11151,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=25.0,
                    # filter names:
                    # hltEle24erWPTightGsfTrackIsoFilterForTau
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    # filter names:
                    # hltSelectedPFTau30LooseChargedIsolationL1HLTMatched
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=1024 + 256,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
        ),
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1",
            id=11152,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=25.0,
                    # filter names:
                    # hltEle24erWPTightGsfTrackIsoFilterForTau
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    # filter names:
                    # hltSelectedPFTau30LooseChargedIsolationL1HLTMatched
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=1024 + 256,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era >= "B"),
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
        ),
        # ===>>> mu-tauh
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1",
            id=13151,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=21.0,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=28.0,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=1024 + 512,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1",
            id=13152,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=21.0,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=1024 + 512,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1",
            id=13153,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=21.0,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=1024 + 512,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        # ===>>> tauh-tauh
        Trigger(
            name="HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg",
            id=15151,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
            id=15152,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
            id=15153,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=37.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=37.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg",
            id=15154,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=37.0,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=37.0,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
    ])    




# ----------------------------------------------- #
#                   Run3 2022                     #
# ----------------------------------------------- #
def add_triggers_run3_2022(config: od.Config, postfix: str) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    PreEE:
    +--------------------------------------------------------------------+----+---------------+---------------+
    |                                HLT                                 | MC | /Tau_Run2022C | /Tau_Run2022D |
    +--------------------------------------------------------------------+----+---------------+---------------+
    |                       HLT_Ele32_WPTight_Gsf                        | √  |       √       |       √       | 352494 - 362760
    |                       HLT_Ele27_WPTight_Gsf                        | √  |       √       |       √       | 352494 - 362760
    |                            HLT_IsoMu27                             | √  |       √       |       √       | 352494 - 362760
    |                            HLT_IsoMu24                             | √  |       √       |       √       | 352494 - 362760
    | HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 | √  |       √       |       √       | 352494 - 362760
    |      HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1      | √  |       √       |       √       | 352494 - 362760
    |          HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1          | √  |       √       |       X       | 352494 - 362760
    |          HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1           | √  |       √       |       X       | 352494 - 362760
    |           HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1            | √  |       √       |       √       | 355862 - 362760
    +--------------------------------------------------------------------+----+---------------+---------------+
    PostEE:
    +--------------------------------------------------------------------+----+---------------+---------------+---------------+
    |                                HLT                                 | MC | /Tau_Run2022E | /Tau_Run2022F | /Tau_Run2022G |
    +--------------------------------------------------------------------+----+---------------+---------------+---------------+
    |                       HLT_Ele32_WPTight_Gsf                        | √  |       √       |       √       |       √       |
    |                       HLT_Ele27_WPTight_Gsf                        | √  |       √       |       √       |       √       |
    |                            HLT_IsoMu27                             | √  |       √       |       √       |       √       |
    |                            HLT_IsoMu24                             | √  |       √       |       √       |       √       |
    | HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 | √  |       √       |       √       |       √       |
    |      HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1      | √  |       √       |       √       |       √       |
    |          HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1          | √  |       √       |       √       |       √       |
    |          HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1           | √  |       √       |       √       |       √       |
    |           HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1            | √  |       √       |       √       |       √       |
    +--------------------------------------------------------------------+----+---------------+---------------+---------------+
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger,[
        # ===>>> single electron
        Trigger(
            name="HLT_Ele32_WPTight_Gsf",
            id=111,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=33.0,
                    # filter names:
                    # hltEle32L1DoubleEGWPTightGsfTrackIsoFilter
                    # hltEGL1SingleEGOrFilter
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        Trigger(
            name="HLT_Ele27_WPTight_Gsf",
            id=112,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=28.0,
                    # filter names:
                    # hltEle32WPTightGsfTrackIsoFilter
                    trigger_bits=2 + 1024,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        # ===>>> single muon
        Trigger(
            name="HLT_IsoMu27",
            id=131,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=28.0,
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu24",
            id=132,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        # ===>>> e-tauh
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1",
            id=11151,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=25.0,
                    # filter names:
                    # hltEle24erWPTightGsfTrackIsoFilterForTau
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    # filter names:
                    # hltSelectedPFTau30LooseChargedIsolationL1HLTMatched
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=1024 + 256,
                ),
            ],
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
        ),
        # ===>>> mu-tauh
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1",
            id=13151,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=21.0,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=1024 + 512,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        # ===>>> tauh-tauh
        # https://cmshltinfo.app.cern.ch/summary?search=&year=2022&paths=true&prescaled=false&stream-types=Physics,Scouting,Parking
        # https://cmshltinfo.app.cern.ch/summary?search=HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1&year=2022&paths=true&prescaled=false&stream-types=Physics,Scouting,Parking
        #Trigger(
        #    name="HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1",
        #    id=15151,
        #    run_range=[352494,362760],
        #    legs=[
        #        TriggerLeg(
        #            pdg_id=15,
        #            min_pt=42.0,
        #            # filter names:
        #            # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
        #            trigger_bits=64,
        #        ),
        #        TriggerLeg(
        #            pdg_id=15,
        #            min_pt=42.0,
        #            # filter names:
        #            # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
        #            trigger_bits=64,
        #        ),
        #    ],
        #    #applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era != "D"),
        #    tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        #),
        # https://cmshltinfo.app.cern.ch/summary?search=HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1&year=2022&paths=true&prescaled=false&stream-types=Physics,Scouting,Parking
        #Trigger(
        #    name="HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1",
        #    id=15152,
        #    run_range=[352494,362760],
        #    legs=[
        #       TriggerLeg(
        #            pdg_id=15,
        #           min_pt=37.0,
        #           # filter names:
        #            # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
        #            trigger_bits=64,
        #        ),
        #        TriggerLeg(
        #            pdg_id=15,
        #            min_pt=37.0,
        #            # filter names:
        #            # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
        #            trigger_bits=64,
        #        ),
        #    ],
        #    #applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era != "D"),
        #    tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        #),
        # https://cmshltinfo.app.cern.ch/summary?search=HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1&year=2022&paths=true&prescaled=false&stream-types=Physics,Scouting,Parking
        Trigger(
            name="HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1",
            id=15153,
            run_range=[355862,362760],
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=37.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=37.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
    ])



# ----------------------------------------------- #
#                   Run3 2023                     #
# ----------------------------------------------- #    
def add_triggers_run3_2023(config: od.Config, postfix: str) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    PreBPix:
    +--------------------------------------------------------------------+----+------------------+
    |                                HLT                                 | MC | /Tau_Run2023C_v1 |
    +--------------------------------------------------------------------+----+------------------+
    |                       HLT_Ele32_WPTight_Gsf                        | √  |        √         |
    |                       HLT_Ele30_WPTight_Gsf                        | √  |        √         |
    |                            HLT_IsoMu27                             | √  |        √         |
    |                            HLT_IsoMu24                             | √  |        √         |
    | HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 | √  |        √         |
    |      HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1      | √  |        √         |
    |           HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1            | √  |        √         |
    +--------------------------------------------------------------------+----+------------------+
    PostBPix:
    +--------------------------------------------------------------------+----+------------------+------------------+
    |                                HLT                                 | MC | /Tau_Run2023D_v1 | /Tau_Run2023D_v2 |
    +--------------------------------------------------------------------+----+------------------+------------------+
    |                       HLT_Ele32_WPTight_Gsf                        | √  |        √         |        √         |
    |                       HLT_Ele30_WPTight_Gsf                        | √  |        √         |        √         |
    |                            HLT_IsoMu27                             | √  |        √         |        √         |
    |                            HLT_IsoMu24                             | √  |        √         |        √         |
    | HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 | √  |        √         |        √         |
    |      HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1      | √  |        √         |        √         |
    |           HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1            | √  |        √         |        √         |
    +--------------------------------------------------------------------+----+------------------+------------------+

    """
    config.x.triggers = od.UniqueObjectIndex(Trigger,[
        # ===>>> single electron
        Trigger(
            name="HLT_Ele32_WPTight_Gsf",
            id=111,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=33.0,
                    # filter names:
                    # hltEle32L1DoubleEGWPTightGsfTrackIsoFilter
                    # hltEGL1SingleEGOrFilter
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        Trigger(
            name="HLT_Ele30_WPTight_Gsf",
            id=112,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=31.0,
                    # filter names:
                    # hltEle32WPTightGsfTrackIsoFilter
                    trigger_bits=2 + 1024,
                ),
            ],
            #applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era >= "D"),
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        # ===>>> single muon
        Trigger(
            name="HLT_IsoMu27",
            id=131,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=28.0,
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu24",
            id=132,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        # ===>>> e-tauh
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1",
            id=11151,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=25.0,
                    # filter names:
                    # hltEle24erWPTightGsfTrackIsoFilterForTau
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    # filter names:
                    # hltSelectedPFTau30LooseChargedIsolationL1HLTMatched
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=1024 + 256,
                ),
            ],
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
        ),
        # ===>>> mu-tauh
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1",
            id=13151,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=21.0,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=1024 + 512,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        # ===>>> tauh-tauh
        Trigger(
            name="HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1",
            id=15151,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=37.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=37.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
    ])













# ----------------------------------------------- #
#                      EXTRAS                     #
# ----------------------------------------------- #
def add_triggers_run3_2022_tau_tau_postEE(config: od.Config) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.

    HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60
    HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75
    HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger,[
        #
        # cross-tau
        #
        Trigger(
            name="HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1",
            id=151,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),

        Trigger(
            name="HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1",
            id=152,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
    ])


def add_triggers_run2_UL2017_mu_tau(config: od.Config) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger,[
        
        #
        # single muon
        #
        Trigger(
            name="HLT_IsoMu24",
            id=131, #13 is for muon pdg_id, 1 because it's first muon trigger
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
    ])
