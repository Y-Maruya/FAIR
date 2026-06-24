#!bin/bash
set -e

start_run_number=$1
num_veto_events=$2

cd ../../
source setup.sh
./bin/fair_calib config/InterCalib_muon_second_calib.yaml -r $start_run_number -n $num_veto_events
cd auto/InterCalib/

echo "Inter-calibration for run $start_run_number is done."