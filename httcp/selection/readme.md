# Structure of the selection in main.py 

- json_filter : golden lumi mask filter 
- trigger_selection
- met_filters
- jet_selection
- muon_selection
- electron_selection
- tau_selection
- selection of event with at least two leptons with at least one tau [before trigger obj matching]
- match_trigobj : trigger obj matching
- check if there are at least two leptons with at least one tau [after trigger obj matching]
- double_lepton_veto
- build hcand with etau_selection, mutau_selection, tautau_selection 
- get_categories : channel selection
- make sure events have at least one lepton pair in hcand (it is only applied on the events with one higgs candidate only)
- gen matching : gentau_selection
- add the mc weight :mc_weight
- rel-charge


## json_filter 


## trigger_selection :HLT trigger path selection.

For each trigger define in the config: 
    Check for each leg of the trigger that : 
        - trigger object pdg = trigger pdg 
        - pt cut
        - trigger bits match
        - all_legs_match ? at least one leg is matched 
    - Trigger fire and at leat one leg is matched : 
            any_fired_all_legs_match = any_fired_all_legs_match | fired_and_all_legs_match ? 

"trigger_ids": contains the trigger IDs of the fired and selected triggers


## met_filters

## jet_selection

## muon_selection
    ```
    good_selections = {
        "muon_pt_26"          : events.Muon.pt > 26,
        "muon_eta_2p4"        : abs(events.Muon.eta) < 2.4,
        "mediumID"            : events.Muon.mediumId == 1,
        "muon_dxy_0p045"      : abs(events.Muon.dxy) < 0.045,
        "muon_dz_0p2"         : abs(events.Muon.dz) < 0.2,
        "muon_iso_0p15"       : events.Muon.pfRelIso04_all < 0.15
    }
    single_veto_selections = {
        "muon_pt_10"          : events.Muon.pt > 10,
        "muon_eta_2p4"        : abs(events.Muon.eta) < 2.4,
        "mediumID"            : events.Muon.mediumId == 1,
        "muon_dxy_0p045"      : abs(events.Muon.dxy) < 0.045,
        "muon_dz_0p2"         : abs(events.Muon.dz) < 0.2,
        "muon_iso_0p3"        : events.Muon.pfRelIso04_all < 0.3
    }
    double_veto_selections = {
        "muon_pt_15"          : events.Muon.pt > 15,
        "muon_eta_2p4"        : abs(events.Muon.eta) < 2.4,
        "muon_isGlobal"       : events.Muon.isGlobal == True,
        "muon_isPF"           : events.Muon.isPFcand == True,
        #"muon_isTracker"      : events.Muon.isTracker ==True,
        "muon_dxy_0p045"      : abs(events.Muon.dxy) < 0.045,
        "muon_dz_0p2"         : abs(events.Muon.dz) < 0.2,
        "muon_iso_0p3"        : events.Muon.pfRelIso04_all < 0.3
    }
    ```

## electron_selection
    ```
    good_selections = {
        "electron_pt_25"          : events.Electron.pt > 25,
        "electron_eta_2p1"        : abs(events.Electron.eta) < 2.1,
        "electron_dxy_0p045"      : abs(events.Electron.dxy) < 0.045,
        "electron_dz_0p2"         : abs(events.Electron.dz) < 0.2,
        "electron_mva_iso_wp80"   : mva_iso_wp80 == 1,
    }
    single_veto_selections = {
        "electron_pt_10"          : events.Electron.pt > 10,
        "electron_eta_2p5"        : abs(events.Electron.eta) < 2.5,
        "electron_dxy_0p045"      : abs(events.Electron.dxy) < 0.045,
        "electron_dz_0p2"         : abs(events.Electron.dz) < 0.2,
        "electron_mva_noniso_wp90": mva_noniso_wp90 == 1,
        "electron_convVeto"       : events.Electron.convVeto == 1,
        #"electron_lostHits"       : events.Electron.lostHits <= 1,
        "electron_pfRelIso03_all" : events.Electron.pfRelIso03_all < 0.3,
    }
    double_veto_selections = {
        "electron_pt_15"          : events.Electron.pt > 15,
        "electron_eta_2p5"        : abs(events.Electron.eta) < 2.5,
        "electron_dxy_0p045"      : abs(events.Electron.dxy) < 0.045,
        "electron_dz_0p2"         : abs(events.Electron.dz) < 0.2,
        "electron_cutBased"       : events.Electron.cutBased == 1,
        "electron_pfRelIso03_all" : events.Electron.pfRelIso03_all < 0.3,
    }
    ```


## tau_selection

# Isolated Tau : pass Medium
    ```
    good_selections = {
        "tau_pt_20"     : events.Tau.pt > 20,
        "tau_eta_2p3"   : abs(events.Tau.eta) < 2.3,
        "tau_dz_0p2"    : abs(events.Tau.dz) < 0.2,
        "DeepTauVSjet"  : events.Tau.idDeepTau2018v2p5VSjet >= tau_tagger_wps.vs_j.Medium,
        "DeepTauVSe"    : events.Tau.idDeepTau2018v2p5VSe   >= tau_tagger_wps.vs_e.VVLoose,
        "DeepTauVSmu"   : events.Tau.idDeepTau2018v2p5VSmu  >= tau_tagger_wps.vs_m.Tight,
        "DecayMode"     : ((events.Tau.decayMode == 0) 
                        | (events.Tau.decayMode == 1)
                        | (events.Tau.decayMode == 2)
                        | (events.Tau.decayMode == 10)
                        | (events.Tau.decayMode == 11))
        #"CleanFromEle"  : ak.all(events.Tau.metric_table(events.Electron[electron_indices]) > 0.5, axis=2),
        #"CleanFromMu"   : ak.all(events.Tau.metric_table(events.Muon[muon_indices]) > 0.5, axis=2),
    }
    ```

## match_trigobj

For each trigger:
    If trigger is_single_mu or is_cross_mu:
        check the muon object
        if trigger is is_single_mu:
            check that there is one leg
            check tigger.legs[0] is muon
            match legs[0] 
            assign to single_muon_triggered
        if trigger is is_cross_mu:
            check that there are 2 legs
            check tigger.legs[0] is muon
            match legs[0] 
            assign to cross_muon_triggered
    if is_single_el or is_cross_el: 
        same as muon
    if is_cross_el or is_cross_mu or is_cross_tau:
        if is_cross_el or is_cross_mu:
            check that there are 2 legs
            check tigger.legs[1] is a tau
            match legs[0] by checking electron and events.TrigObj[leg_masks[0]] are in cone of dR O.5
            assign to is_cross_el or is_cross_mu
        if is_cross_tau: 
            check that there are 2 legs
            check tigger.legs[0] is a tau
            check tigger.legs[1] is a tau
            match both legs


## double_lepton_veto
inverse charge 
dR>0.15

## etau_selection
Sorting lep1 [Electron] by isolation [ascending]
Sorting lep2 [Tau] by DeepTau [descending]
Create pair of leps : probable higgs candidate -> leps_pair
e.g. [ [(e1,t1)], [(e1,t1),(e1,t2)], [(e1,t1),(e2,t1)], [], [(e1,t1),(e1,t2),(e2,t1),(e2,t2)] ]
and their indices -> lep_indices_pair 
```
 preselection = {
        "etau_is_os"         : (lep1.charge * lep2.charge) < 0,
        "etau_dr_0p5"        : (1*lep1).delta_r(1*lep2) > 0.5,  # deltaR(lep1, lep2) > 0.5,
        #"etau_mT_50"         : transverse_mass(lep1, events.MET) < 50
        "etau_mT_50"         : transverse_mass(lep1, met) < 50
    }
```
If multiple pairs, sort them and get first one.
Sorting:
- sorted by their isolation
- Sort the pairs based on pfRelIso03_all of the first object in each pair
- Check if the pfRelIso03_all values are the same for the first two objects in each pair
- Sort the pairs based on pt if pfRelIso03_all is the same for the first two objects
- Check if the pt values are the same for the first two objects in each pair
- if so, sort the pairs with tau rawDeepTau2017v2p1VSjet
- check if the first two pairs have taus with same rawDeepTau2017v2p1VSjet
- Sort the pairs based on pt if rawDeepTau2017v2p1VSjet is the same for the first two objects
- Extract the first object in each pair (lep1) and the second object (lep2)
- Concatenate lep1 and lep2 to create the final dtrpair

## mutau_selection
Similair to etau 
```
    preselection = {
        "is_os"         : (lep1.charge * lep2.charge) < 0,
        "dr_0p5"        : deltaR(lep1, lep2) > 0.5,
        "mT_50"         : transverse_mass(lep1, events.MET) < 50
    }
```
+ sorting

## tautau_selection
Similair to etau 
```
preselection = {
        "tautau_is_pt_40"      : (lep1.pt > 40) & (lep2.pt > 40),
        "tautau_is_eta_2p1"    : (np.abs(lep1.eta) < 2.1) & (np.abs(lep2.eta) < 2.1),
        "tautau_is_os"         : (lep1.charge * lep2.charge) < 0,
        "tautau_dr_0p5"        : (1*lep1).delta_r(1*lep2) > 0.5,  #deltaR(lep1, lep2) > 0.5,
    }
```
+ sorting

## get_categories

## gen matching 