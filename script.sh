#!/bin/bash
for run in 22334; do
    mkdir -p out_calibration/pedestal_calib/${run}
    ./bin/fair_single config/ped_${run}.yaml -i config/Input/ped_${run}.txt &
done