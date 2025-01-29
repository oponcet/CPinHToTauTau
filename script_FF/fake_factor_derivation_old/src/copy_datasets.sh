#!/bin/bash

# Define the dataset lists
data=("data_tau_C" "data_tau_D")
mc=(
    "wj_incl_madgraph"
    "dy_lep_m10to50_madgraph" "dy_lep_m50_madgraph" "dy_lep_m50_1j_madgraph" "dy_lep_m50_2j_madgraph"
    "dy_lep_m50_3j_madgraph" "dy_lep_m50_4j_madgraph" "tt_dl" "tt_sl" "tt_fh"
    "st_tchannel_t" "st_tchannel_tbar" "st_tw_t_sl" "st_tw_t_dl" "st_tw_t_fh" "st_tw_tb_sl"
    "st_tw_tb_dl" "st_tw_tb_fh" "st_schannel_t" "st_schannel_tbar" "ww" "wz" "zz" "www" "wwz"
    "wzz" "zzz"
)

# Define base paths
source_base="/eos/user/g/gsaha/CPinHToTauTauOutput/cf_store/analysis_httcp/cf.MergeHistograms/run3_2022_preEE_nano_cp_tau_v14"
dest_base="/eos/user/o/oponcet/TauCP/FakeFactors/MergeHisto/IC_selection_Run3_2022PreEE_full_20250118_215651"

# Define the common subpath to copy
subpath="nominal/calib__main/sel__main/prod__main/weight__main/Run3_2022PreEE_full_20250118_215651"

# Function to copy directories
copy_directories() {
    local datasets=("$@")
    for dataset in "${datasets[@]}"; do
        source_dir="${source_base}/${dataset}/${subpath}"
        dest_dir="${dest_base}/${dataset}"
        
        echo "Copying from $source_dir to $dest_dir..."
        
        # Perform the copy
        cp -r "$source_dir" "$dest_dir"
        
        # Check if the copy was successful
        if [ $? -eq 0 ]; then
            echo "Successfully copied $dataset."
        else
            echo "Failed to copy $dataset. Check paths and permissions."
        fi
    done
}

# Copy data and MC datasets
echo "Starting to copy data datasets..."
copy_directories "${data[@]}"
echo "Starting to copy MC datasets..."
copy_directories "${mc[@]}"
echo "All datasets processed."
