#!bin/bash

start_run_number=$1

source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh
cd ../../
./bin/fair_calib config/InterCalib.yaml -r $start_run_number
cd auto/InterCalib/

echo "Inter-calibration for run $start_run_number is done."