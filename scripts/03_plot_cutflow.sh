#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes $processes
        --version $version
        --datasets $datasets
        --workflow local
        --selector-steps "trigger,b_veto,multiple_leptons,has_higgs_cand,extra_lepton_veto,dilepton_veto"  
    )
law run cf.PlotCutflow "${args[@]}"