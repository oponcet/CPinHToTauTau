import itertools
import json

# Define the categories and their IDs
channel = [{"name": "tautau__real_1", "id": 50000}]
abcd = [
    {"name": "hadA", "id": 100},
    {"name": "hadB", "id": 300},
    {"name": "hadA0", "id": 500},
    {"name": "hadB0", "id": 700},
    {"name": "hadC0", "id": 1100},
    {"name": "hadD0", "id": 900},
    {"name": "hadC", "id": 1500},
    {"name": "hadD", "id": 1300},
]
nj = [
    {"name": "has_0j", "id": 51},
    {"name": "has_1j", "id": 53},
    {"name": "has_2j", "id": 55},
]
dm = [
    {"name": "pi_1", "id": 1},
    {"name": "rho_1", "id": 3},
    {"name": "a1dm2_1", "id": 5},
    {"name": "a1dm10_1", "id": 7},
    {"name": "a1dm11_1", "id": 9},
]

# Generate all combinations
combinations = itertools.product(channel, abcd, nj, dm)

# Prepare a dictionary to organize combinations
categories = {}

for combo in combinations:
    combo_name = "__".join(cat["name"] for cat in combo)
    combo_id = sum(cat["id"] for cat in combo)
    
    # Extract DM name and jet category (has_0j, has_1j, has_2j) from the combination
    dm_name = combo[3]["name"]  # e.g., "pi_1", "rho_1", etc.
    jet_category = combo[2]["name"]  # e.g., "has_0j", "has_1j", "has_2j"
    abcd_category = combo[1]["name"]  # e.g., "hadA", "hadB", etc.
    
    # Create structure if not already in the dictionary
    if dm_name not in categories:
        categories[dm_name] = {}

    if jet_category not in categories[dm_name]:
        categories[dm_name][jet_category] = {
            "ABCD": {},
            "A0B0C0D0": {}
        }
    
    # Add combination to the appropriate category and type
    if abcd_category == "hadA0" or abcd_category == "hadB0" or abcd_category == "hadC0" or abcd_category == "hadD0":
        categories[dm_name][jet_category]["A0B0C0D0"][combo_name] = combo_id
    else:
        categories[dm_name][jet_category]["ABCD"][combo_name] = combo_id

# Define the JSON structure
data_config = {
    "paths": {
        "eos_path": "/eos/user/g/gsaha/CPinHToTauTauOutput/cf_store/analysis_httcp/cf.MergeHistograms/",
        "task": "run3_2022_preEE_nano_cp_tau_v14/",
        "hist_path": "nominal/calib__main/sel__main/prod__main/weight__main/Run3_2022PreEE_full_20241130_190734/"
    },
    "datasets": {
        "data": [
            "data_tau_C", "data_tau_D", "data_mu_C", "data_mu_D", "data_e_C", "data_e_D", "data_single_mu_C"
        ],
        "mc": [
            "wj_incl_madgraph", "wj_1j_madgraph", "wj_2j_madgraph", "wj_3j_madgraph", "wj_4j_madgraph",
            "dy_lep_m10to50_madgraph", "dy_lep_m50_madgraph", "dy_lep_m50_1j_madgraph", "dy_lep_m50_2j_madgraph", 
            "dy_lep_m50_3j_madgraph", "dy_lep_m50_4j_madgraph", "tt_dl", "tt_sl", "tt_fh",
            "st_tchannel_t", "st_tchannel_tbar", "st_tw_t_sl", "st_tw_t_dl", "st_tw_t_fh", "st_tw_tb_sl",
            "st_tw_tb_dl", "st_tw_tb_fh", "st_schannel_t", "st_schannel_tbar", "ww", "wz", "zz", "www", "wwz",
            "wzz", "zzz", "h_ggf_tautau_uncorrelatedDecay_CPodd_Filtered_ProdAndDecay", 
            "h_ggf_tautau_uncorrelatedDecay_MM_Filtered_ProdAndDecay", "h_ggf_tautau_uncorrelatedDecay_SM_Filtered_ProdAndDecay"
        ]
    },
    "categories": categories,
    "variables": [
        {
            "var1": "dphi_met_h1",
            "var2": "hcand_1_pt"
        },
        {
            "var1": "hcand_1_eta",
            "var2": "hcand_1_pt"
        },
        {
            "var1": "hcand_1_phi",
            "var2": "hcand_1_pt"
        },
        {
            "var1": "hcand_dr",
            "var2": "hcand_1_pt"
        },
        {
            "var1": "hcand_invm",
            "var2": "hcand_1_pt"
        },
        {
            "var1": "jet_1_pt",
            "var2": "hcand_1_pt"
        },
        {
            "var1": "jet_2_pt",
            "var2": "hcand_1_pt"
        },
        {
            "var1": "met_var_qcd_h1",
            "var2": "hcand_1_pt"
        },
        {
            "var1": "puppi_met_phi",
            "var2": "hcand_1_pt"
        },
        {
            "var1": "puppi_met_pt",
            "var2": "hcand_1_pt"
        }
    ]
}

# Specify the filename for the JSON file
filename = 'script_FF/fake_factors_v3dec.json'

# Save the data to a JSON file with pretty formatting
with open(filename, 'w') as json_file:
    json.dump(data_config, json_file, indent=4)

print(f"JSON file saved as {filename}")



''' 
tautau__real_1__hadA__has_0j__pi_1 => 50152
tautau__real_1__hadA__has_0j__rho_1 => 50154
tautau__real_1__hadA__has_0j__a1dm10_1 => 50158
tautau__real_1__hadA__has_0j__a1dm11_1 => 50160
tautau__real_1__hadA__has_1j__pi_1 => 50154
tautau__real_1__hadA__has_1j__rho_1 => 50156
tautau__real_1__hadA__has_1j__a1dm10_1 => 50160
tautau__real_1__hadA__has_1j__a1dm11_1 => 50162
tautau__real_1__hadA__has_2j__pi_1 => 50156
tautau__real_1__hadA__has_2j__rho_1 => 50158
tautau__real_1__hadA__has_2j__a1dm10_1 => 50162
tautau__real_1__hadA__has_2j__a1dm11_1 => 50164
tautau__real_1__hadB__has_0j__pi_1 => 50352
tautau__real_1__hadB__has_0j__rho_1 => 50354
tautau__real_1__hadB__has_0j__a1dm10_1 => 50358
tautau__real_1__hadB__has_0j__a1dm11_1 => 50360
tautau__real_1__hadB__has_1j__pi_1 => 50354
tautau__real_1__hadB__has_1j__rho_1 => 50356
tautau__real_1__hadB__has_1j__a1dm10_1 => 50360
tautau__real_1__hadB__has_1j__a1dm11_1 => 50362
tautau__real_1__hadB__has_2j__pi_1 => 50356
tautau__real_1__hadB__has_2j__rho_1 => 50358
tautau__real_1__hadB__has_2j__a1dm10_1 => 50362
tautau__real_1__hadB__has_2j__a1dm11_1 => 50364
tautau__real_1__hadA0__has_0j__pi_1 => 50552
tautau__real_1__hadA0__has_0j__rho_1 => 50554
tautau__real_1__hadA0__has_0j__a1dm10_1 => 50558
tautau__real_1__hadA0__has_0j__a1dm11_1 => 50560
tautau__real_1__hadA0__has_1j__pi_1 => 50554
tautau__real_1__hadA0__has_1j__rho_1 => 50556
tautau__real_1__hadA0__has_1j__a1dm10_1 => 50560
tautau__real_1__hadA0__has_1j__a1dm11_1 => 50562
tautau__real_1__hadA0__has_2j__pi_1 => 50556
tautau__real_1__hadA0__has_2j__rho_1 => 50558
tautau__real_1__hadA0__has_2j__a1dm10_1 => 50562
tautau__real_1__hadA0__has_2j__a1dm11_1 => 50564
tautau__real_1__hadB0__has_0j__pi_1 => 50752
tautau__real_1__hadB0__has_0j__rho_1 => 50754
tautau__real_1__hadB0__has_0j__a1dm10_1 => 50758
tautau__real_1__hadB0__has_0j__a1dm11_1 => 50760
tautau__real_1__hadB0__has_1j__pi_1 => 50754
tautau__real_1__hadB0__has_1j__rho_1 => 50756
tautau__real_1__hadB0__has_1j__a1dm10_1 => 50760
tautau__real_1__hadB0__has_1j__a1dm11_1 => 50762
tautau__real_1__hadB0__has_2j__pi_1 => 50756
tautau__real_1__hadB0__has_2j__rho_1 => 50758
tautau__real_1__hadB0__has_2j__a1dm10_1 => 50762
tautau__real_1__hadB0__has_2j__a1dm11_1 => 50764
tautau__real_1__hadC0__has_0j__pi_1 => 51152
tautau__real_1__hadC0__has_0j__rho_1 => 51154
tautau__real_1__hadC0__has_0j__a1dm10_1 => 51158
tautau__real_1__hadC0__has_0j__a1dm11_1 => 51160
tautau__real_1__hadC0__has_1j__pi_1 => 51154
tautau__real_1__hadC0__has_1j__rho_1 => 51156
tautau__real_1__hadC0__has_1j__a1dm10_1 => 51160
tautau__real_1__hadC0__has_1j__a1dm11_1 => 51162
tautau__real_1__hadC0__has_2j__pi_1 => 51156
tautau__real_1__hadC0__has_2j__rho_1 => 51158
tautau__real_1__hadC0__has_2j__a1dm10_1 => 51162
tautau__real_1__hadC0__has_2j__a1dm11_1 => 51164
tautau__real_1__hadD0__has_0j__pi_1 => 50952
tautau__real_1__hadD0__has_0j__rho_1 => 50954
tautau__real_1__hadD0__has_0j__a1dm10_1 => 50958
tautau__real_1__hadD0__has_0j__a1dm11_1 => 50960
tautau__real_1__hadD0__has_1j__pi_1 => 50954
tautau__real_1__hadD0__has_1j__rho_1 => 50956
tautau__real_1__hadD0__has_1j__a1dm10_1 => 50960
tautau__real_1__hadD0__has_1j__a1dm11_1 => 50962
tautau__real_1__hadD0__has_2j__pi_1 => 50956
tautau__real_1__hadD0__has_2j__rho_1 => 50958
tautau__real_1__hadD0__has_2j__a1dm10_1 => 50962
tautau__real_1__hadD0__has_2j__a1dm11_1 => 50964
tautau__real_1__hadC__has_0j__pi_1 => 51552
tautau__real_1__hadC__has_0j__rho_1 => 51554
tautau__real_1__hadC__has_0j__a1dm10_1 => 51558
tautau__real_1__hadC__has_0j__a1dm11_1 => 51560
tautau__real_1__hadC__has_1j__pi_1 => 51554
tautau__real_1__hadC__has_1j__rho_1 => 51556
tautau__real_1__hadC__has_1j__a1dm10_1 => 51560
tautau__real_1__hadC__has_1j__a1dm11_1 => 51562
tautau__real_1__hadC__has_2j__pi_1 => 51556
tautau__real_1__hadC__has_2j__rho_1 => 51558
tautau__real_1__hadC__has_2j__a1dm10_1 => 51562
tautau__real_1__hadC__has_2j__a1dm11_1 => 51564
tautau__real_1__hadD__has_0j__pi_1 => 51352
tautau__real_1__hadD__has_0j__rho_1 => 51354
tautau__real_1__hadD__has_0j__a1dm10_1 => 51358
tautau__real_1__hadD__has_0j__a1dm11_1 => 51360
tautau__real_1__hadD__has_1j__pi_1 => 51354
tautau__real_1__hadD__has_1j__rho_1 => 51356
tautau__real_1__hadD__has_1j__a1dm10_1 => 51360
tautau__real_1__hadD__has_1j__a1dm11_1 => 51362
tautau__real_1__hadD__has_2j__pi_1 => 51356
tautau__real_1__hadD__has_2j__rho_1 => 51358
tautau__real_1__hadD__has_2j__a1dm10_1 => 51362
tautau__real_1__hadD__has_2j__a1dm11_1 => 51364

'''
