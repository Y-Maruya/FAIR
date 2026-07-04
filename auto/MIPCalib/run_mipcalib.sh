#!/bin/bash
set -e

start_run_number=$1
NUM_VETO_EVENTS=$2

cd ../../
source setup.sh
./bin/fair_calib /eos/user/y/ymaruya/FASER/AHCAL/FAIR/config/MIP_second_calib_later22789.yaml -r $start_run_number -n $NUM_VETO_EVENTS
cd auto/MIPCalib/

echo "MIP-calibration for run $start_run_number is done."
