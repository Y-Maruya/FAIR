#!bin/bash
set -e

start_run_number=$1

cd ../../
source setup.sh
./bin/fair_calib config/SPE.yaml -r $start_run_number
cd auto/SPECalib/

echo "SPE-calibration for run $start_run_number is done."