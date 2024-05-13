#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
wrapper_args=(
    --configs $config
    --datasets $datasets
    --version $version
    --cf.ReduceEvents-workflow local
    --cf.ReduceEvents-branch -1
    "${@:2}"
    )
echo law run cf.ReduceEventsWrapper "${wrapper_args[@]}"
law run cf.ReduceEventsWrapper "${wrapper_args[@]}"