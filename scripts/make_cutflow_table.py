import os
import sys
import pickle
import numpy as np
import awkward as ak
from tabulate import tabulate

lumi = 7.98
categories_to_check = [1,101,201,301]
isMC = True
campaign      = "run3_2022_preEE_nano_cp_tau_v12"
version       = "Run3_2022PreEE_full_20240825_203729"
cutflowpath   = "/afs/cern.ch/work/g/gsaha/public/IPHC/Work/ColumnFlowAnalyses/CPinHToTauTau/data/cf_store/analysis_httcp/cf.CreateCutflowHistograms"
mergehistpath = "/afs/cern.ch/work/g/gsaha/public/IPHC/Work/ColumnFlowAnalyses/CPinHToTauTau/data/cf_store/analysis_httcp/cf.MergeHistograms"


files = {
    "signal"     : [f"{cutflowpath}/{campaign}/h_ggf_tautau_prod_cp_even_sm/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/h_ggf_tautau_prod_cp_even_sm/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    3.274821,
                    175963.0],
    "wj_incl"    : [f"{cutflowpath}/{campaign}/wj_incl/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/wj_incl/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    63425.1,
                    10626959.0],
    "dy_lep_m50" : [f"{cutflowpath}/{campaign}/dy_lep_m50/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/dy_lep_m50/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    6282.6,
                    145286646.0],
    "ww"         : [f"{cutflowpath}/{campaign}/ww/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/ww/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    122.3,
                    15405496.0],
    "wz"         : [f"{cutflowpath}/{campaign}/wz/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/wz/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    41.1,
                    7527043.0],
    "zz"         : [f"{cutflowpath}/{campaign}/zz/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/zz/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    19.4,
                    1181750.0],
    "tt_dl"      : [f"{cutflowpath}/{campaign}/tt_dl/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/tt_dl/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    98.04,
                    21890633.0],
    "tt_sl"      : [f"{cutflowpath}/{campaign}/tt_sl/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/tt_sl/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    405.75,
                    133416781.0],
    "tt_fh"      : [f"{cutflowpath}/{campaign}/tt_fh/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/tt_fh/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    419.80,
                    42646240.0],
}

"""

files = {
    "data_e_C"       : [f"{cutflowpath}/{campaign}/data_e_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_e_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        1000, 175963.0],
    "data_e_D"       : [f"{cutflowpath}/{campaign}/data_e_D/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_e_D/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        63425.1, 59463375.0],
    "data_mu_C"      : [f"{cutflowpath}/{campaign}/data_mu_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_mu_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        1000, 175963.0],
    "data_mu_D"      : [f"{cutflowpath}/{campaign}/data_mu_D/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_mu_D/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        63425.1, 59463375.0],
    "data_single_mu_C" : [f"{cutflowpath}/{campaign}/data_single_mu_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                          f"{mergehistpath}/{campaign}/data_single_mu_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                          63425.1, 59463375.0],
    "data_tau_C"     : [f"{cutflowpath}/{campaign}/data_tau_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_tau_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        1000, 175963.0],
    "data_tau_D"     : [f"{cutflowpath}/{campaign}/data_tau_D/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_tau_D/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        63425.1, 59463375.0],
}
"""

selections = [
    'Initial',
    'json',
    'trigger',
    'met_filter',
    'b_veto',
    'has_2_or_more_leps_with_at_least_1_tau',
    'dilepton_veto',
    'has_at_least_1_pair_before_trigobj_matching',
    'has_at_least_1_pair_after_trigobj_matching',
    'extra_lepton_veto',
    'One_higgs_cand_per_event',
    'has_proper_tau_decay_products'
]
# temporary
if isMC: selections.remove('json')


file_names = list(files.keys())
files_count_dict = {}
for name, filelist in files.items():
    file = filelist[0]
    #print(f"File: {file}")
    # load pickle file
    file_obj = open(file, 'rb')
    data = pickle.load(file_obj)

    axes = data.axes
    category_axis  = axes['category']
    category_idxs  = [category_axis.index(catid) for catid in categories_to_check]
    selection_axis = axes['step']
    
    count_dict = {key:[] for key in categories_to_check}
    for i,catidx in enumerate(category_idxs):
        catid = categories_to_check[i]
        #print(f"categoryID : {catid}")
        for j in range(len(selections)):
            count = data[catidx, :, j, :, :].values()[0,0,0]
            count_dict[catid].append(count)
            #print(f"\t {selections[j]} : {count}")

    files_count_dict[name] = count_dict


#print(files_count_dict)

rev_dict = {key:{} for key in categories_to_check}
counts_from_merged_hist = []
for file_name, count_dict in files_count_dict.items():
    #print(f"file : {file_name}")
    mergehistfile = files[file_name][1]
    fptr = open(mergehistfile, 'rb')
    merge_data = pickle.load(fptr)
    
    xsec = files[file_name][2]
    nev  = files[file_name][-1]
    temp_merged = []
    for key, val in count_dict.items():
        cat_idx = merge_data.axes['category'].index(key)
        count_from_merged_hist = merge_data[cat_idx,:,:,:].values()[0,0,0]
        temp_merged.append(count_from_merged_hist)
        #print(f"cat id: {key}")
        counts = ak.Array(val)
        # abs eff
        abs_eff = ak.values_astype(counts/(val[0]*ak.ones_like(counts)), np.float32)
        # rel eff
        den = ak.Array([val[0]]+val[:-1])
        rel_eff = ak.values_astype(counts/den, np.float32)

        count_1 = (counts[0] * xsec * lumi * 1000)/nev if isMC else counts[0]
        count_2 = (counts[-1] * xsec * lumi * 1000)/nev if isMC else counts[-1]

        abs_eff = ak.to_list(abs_eff)
        ext_abs_eff = np.round(np.array([count_1] + abs_eff + [count_2]), decimals=3)[:,None]
        
        rev_dict[key][file_name] = ext_abs_eff
    counts_from_merged_hist.append(temp_merged)

#print(rev_dict)

#print(counts_from_merged_hist)
counts_from_merged_hist = np.round(np.array(counts_from_merged_hist), decimals=3)
#print(counts_from_merged_hist)
#from IPython import embed; embed()
counts_from_merged_hist = counts_from_merged_hist.T

headers = ["selections"] + file_names

for j, (key, val) in enumerate(rev_dict.items()):
    print(f"cat id: {key}")
    sels = ak.to_numpy(ak.Array(["Begins with (LumiScaled)"] + selections + ["Ends with (LumiScaled)"]))[:,None]
    for file, arr in val.items():
        sels = np.concatenate([sels, arr], axis=1)
    #from IPython import embed; embed()
    y = sels[:,1:][-1]
    y = y.astype(np.float64)
    total_y = round(np.sum(y), 2)
    sels = sels.tolist()
    sels.append(["From Merged Hist"]+counts_from_merged_hist[j].tolist())
    z = round(np.sum(counts_from_merged_hist[j]), 2)
    table = tabulate(sels, headers, tablefmt="pretty")
    #from IPython import embed; embed()
    print(table)
    print(f"Total Yield      : {total_y}")
    print(f"From Merged Hist : {z}\n")




"""
files = {
    "signal"     : [f"{cutflowpath}/{campaign}/h_ggf_tautau_prod_cp_even_sm/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/h_ggf_tautau_prod_cp_even_sm/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    3.274821,
                    175963.0],
    "wj_incl"    : [f"{cutflowpath}/{campaign}/wj_incl/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/wj_incl/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    63425.1,
                    10626959.0],
    "dy_lep_m50" : [f"{cutflowpath}/{campaign}/dy_lep_m50/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/dy_lep_m50/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    6282.6,
                    145286646.0],
    "ww"         : [f"{cutflowpath}/{campaign}/ww/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/ww/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    122.3,
                    15405496.0],
    "wz"         : [f"{cutflowpath}/{campaign}/wz/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/wz/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    41.1,
                    7527043.0],
    "zz"         : [f"{cutflowpath}/{campaign}/zz/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/zz/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    19.4,
                    1181750.0],
    "tt_dl"      : [f"{cutflowpath}/{campaign}/tt_dl/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/tt_dl/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    98.04,
                    21890633.0],
    "tt_sl"      : [f"{cutflowpath}/{campaign}/tt_sl/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/tt_sl/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    405.75,
                    133416781.0],
    "tt_fh"      : [f"{cutflowpath}/{campaign}/tt_fh/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_has_2_or_more_leps_with_at_least_1_tau_dilepton_veto_has_at_least_1_pair_before_trigobj_matching_has_at_least_1_pair_after_trigobj_matching_extra_lepton_veto_One_higgs_cand_per_event_has_proper_tau_decay_products/{version}/cutflow_hist__event.pickle",
                    f"{mergehistpath}/{campaign}/tt_fh/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                    419.80,
                    42646240.0],
}

"""
"""
files = {
    "data_e_C"       : [f"{cutflowpath}/{campaign}/data_e_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_e_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        1000, 175963.0],
    "data_e_D"       : [f"{cutflowpath}/{campaign}/data_e_D/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_e_D/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        63425.1, 59463375.0],
    "data_mu_C"      : [f"{cutflowpath}/{campaign}/data_mu_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_mu_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        1000, 175963.0],
    "data_mu_D"      : [f"{cutflowpath}/{campaign}/data_mu_D/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_mu_D/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        63425.1, 59463375.0],
    "data_single_mu_C" : [f"{cutflowpath}/{campaign}/data_single_mu_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                          f"{mergehistpath}/{campaign}/data_single_mu_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                          63425.1, 59463375.0],
    "data_tau_C"     : [f"{cutflowpath}/{campaign}/data_tau_C/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_tau_C/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        1000, 175963.0],
    "data_tau_D"     : [f"{cutflowpath}/{campaign}/data_tau_D/nominal/calib__main/sel__main/{version}/cutflow_hist__event.pickle",
                        f"{mergehistpath}/{campaign}/data_tau_D/nominal/calib__main/sel__main/prod__main/{version}/hist__event.pickle",
                        63425.1, 59463375.0],
}
"""
