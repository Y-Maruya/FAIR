#!bin/bash

start_run_number=$1

cd ../../
source setup.sh
./bin/fair_calib config/InterCalib.yaml -r $start_run_number 
cd auto/InterCalib/

echo "Inter-calibration for run $start_run_number is done."