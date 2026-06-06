#!bin/bash

runnumber=$1


cp ped.yaml ped_${runnumber}.yaml
python3 changepath.py ped_${runnumber}.yaml --run-number $runnumber
python3 makerunlist_byrun.py --runnumber $runnumber
cd ../../
source setup.sh
mkdir -p out_calibration/pedestal_calib/${runnumber}
mkdir -p out_calibration/pedestal_calib/${runnumber}/json
./bin/fair_single auto/pedestal/ped_${runnumber}.yaml -i auto/pedestal/ped_${runnumber}.txt

cd auto/pedestal/
rm ped_${runnumber}.txt
rm ped_${runnumber}.yaml

echo "Pedestal calibration for run $runnumber is done."