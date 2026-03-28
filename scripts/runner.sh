#!bin/bash

for run in {20000..24000..1}
do
    echo "Processing run $run"
    # Run the calibration script for the current run
     python3 PrePostTriggersArchiver.py --run $run
done