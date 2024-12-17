'''
author: @oponcet
date: 06-12-2024
Description:
This script generates JSON configuration files for different decay modes (DMs) in the context of Fake Factor derivation. 

Key Features:
1. **Define Categories**:
   - `channel`: Specifies the type of event or process (e.g., `tautau__real_1`).
   - `abcd`: Represents control regions (e.g., `hadA`, `hadB`, `hadC0`, etc.).
   - `nj`: Indicates the number of jets (e.g., `has_0j`, `has_1j`, etc.).
   - `dm`: Specifies decay modes (e.g., `pi_1`, `rho_1`).

2. **Generate Combinations**:
   - Creates all possible combinations of the above categories.
   - Each combination is given a unique name and ID defined in httcp/main/category.py. 

3. **Organize Data by Decay Mode (DM)**:
   - Groups combinations into:
     - `ABCD`: General control categories.
     - `A0B0C0D0`: Special control categories (e.g., `hadA0`, `hadB0`).
   - Further organizes by decay mode (`dm_name`) and jet multiplicity.

4. **JSON Structure**:
   - Includes:
     - File paths for input data and histograms.
     - Datasets for data and Monte Carlo (MC) samples.
     - Variables used for analysis.

5. **Save Per-Decay Mode Files**:
   - Generates one JSON file per decay mode.
   - Each file contains only the relevant combinations and metadata.

Output:
- The JSON files are saved in the `script_FF` directory with filenames like `fake_factors_pi_1.json`, where `pi_1` is the decay mode.
- Each file includes categories, datasets, paths, and variables specific to the given decay mode.

Use Case:
- This script simplifies the configuration process for Fake Factor analysis by automating the generation of structured JSON files for each decay mode.
- It is designed to be modular, scalable, and reusable for high-energy physics workflows.
'''

import itertools
import json
import os

# Define the categories and their IDs
channel = [{"name": "tautau__real_1", "id": 41000000}] # 1000000 real_1 and 40000000 tautau
abcd = [
    {"name": "hadA", "id": 1000},
    {"name": "hadB", "id": 3000},
    {"name": "hadA0", "id": 5000},
    {"name": "hadB0", "id": 7000},
    {"name": "hadC0", "id": 11000},
    {"name": "hadD0", "id": 9000},
    {"name": "hadC", "id": 15000},
    {"name": "hadD", "id": 13000},
]
nj = [
    {"name": "has_0j", "id": 100000},
    {"name": "has_1j", "id": 200000},
    {"name": "has_2j", "id": 300000},
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
categories_by_dm = {}

for combo in combinations:
    combo_name = "__".join(cat["name"] for cat in combo)
    combo_id = sum(cat["id"] for cat in combo)
    
    # Extract DM name and jet category (has_0j, has_1j, has_2j) from the combination
    dm_name = combo[3]["name"]  # e.g., "pi_1", "rho_1", etc.
    jet_category = combo[2]["name"]  # e.g., "has_0j", "has_1j", "has_2j"
    abcd_category = combo[1]["name"]  # e.g., "hadA", "hadB", etc.
    
    # Initialize structure if not present for the current DM
    if dm_name not in categories_by_dm:
        categories_by_dm[dm_name] = {}

    if jet_category not in categories_by_dm[dm_name]:
        categories_by_dm[dm_name][jet_category] = {
            "ABCD": {},
            "A0B0C0D0": {}
        }
    
    # Add combination to the appropriate category and type
    if abcd_category == "hadA0" or abcd_category == "hadB0" or abcd_category == "hadC0" or abcd_category == "hadD0":
        categories_by_dm[dm_name][jet_category]["A0B0C0D0"][combo_name] = combo_id
    else:
        categories_by_dm[dm_name][jet_category]["ABCD"][combo_name] = combo_id

# Base JSON structure for each DM
base_data_config = {
    "paths": {
        "eos_path": "/eos/user/g/gsaha/CPinHToTauTauOutput/cf_store/analysis_httcp/cf.MergeHistograms/", 
        "task": "run3_2022_preEE_nano_cp_tau_v14/",
        "hist_path": "nominal/calib__main/sel__main/prod__main/weight__main/Run3_2022PreEE_full_20241206_193602/"
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
            "wzz", "zzz"
        ]
    },
    "variables": [
        {
            "var1": "dphi_met_h1",
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
            "var1": "hcand_2_pt",
            "var2": "hcand_1_pt"
        },
        {
            "var1": "hT",
            "var2": "hcand_1_pt"
        },
        {
            "var1": "n_jet",
            "var2": "hcand_1_pt"
        }
    ]
}

# Create output directory if it doesn't exist
output_dir = '/afs/cern.ch/user/o/oponcet/private/analysis/CPinHToTauTau/script_FF/fake_factor_derivation/inputs/inputs_json'
os.makedirs(output_dir, exist_ok=True)

# Save a separate JSON file for each DM
for dm_name, dm_categories in categories_by_dm.items():
    dm_data_config = base_data_config.copy()
    dm_data_config["categories"] = {dm_name: dm_categories}
    
    filename = f"{output_dir}/fake_factors_{dm_name}.json"
    with open(filename, 'w') as json_file:
        json.dump(dm_data_config, json_file, indent=4)
    
    print(f"JSON file for {dm_name} saved as {filename}")


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