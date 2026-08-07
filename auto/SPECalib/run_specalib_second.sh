#!bin/bash
set -e

start_run_number=$1
end_run_number=$2

cd ../../
source setup.sh
./bin/fair_calib config/SPE_second_calib.yaml -r $start_run_number -l $end_run_number
cd auto/SPECalib/

echo "SPE-calibration for run $start_run_number is done."