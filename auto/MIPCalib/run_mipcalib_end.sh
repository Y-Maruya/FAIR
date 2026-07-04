#!/bin/bash
set -e

start_run_number=$1
end_run_number=$2

cd ../../
source setup.sh
if [ $start_run_number -gt 22560 ]; then
    ./bin/fair_calib /eos/user/y/ymaruya/FASER/AHCAL/FAIR/config/MIP_second_calib_later22562.yaml -r $start_run_number -l $end_run_number
else
    ./bin/fair_calib /eos/user/y/ymaruya/FASER/AHCAL/FAIR/config/MIP_second_calib.yaml -r $start_run_number -l $end_run_number
fi
# ./bin/fair_calib /eos/user/y/ymaruya/FASER/AHCAL/FAIR/config/MIP_second_calib.yaml -r $start_run_number -l $end_run_number
cd auto/MIPCalib/

echo "MIP-calibration for run $start_run_number is done."
