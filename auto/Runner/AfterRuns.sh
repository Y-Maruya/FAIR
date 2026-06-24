##!/bin/bash

config_file=$1
start_run_number=$2
num_veto_events=$3

cd /eos/user/y/ymaruya/FASER/AHCAL/FAIR/
# set +u
source ./setup.sh

./bin/fair_calib $config_file -r $start_run_number -n $num_veto_events

echo "Finished running $config_file for run $start_run_number with $num_veto_events veto events"