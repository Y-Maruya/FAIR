#!/bin/bash
set -eo pipefail

cd /eos/user/y/ymaruya/FASER/AHCAL/FAIR/auto/SPECalib
set +u
source ../../setup.sh
echo "Running SPE-calibration for run $1 to $2"
bash run_specalib_second.sh "$1" "$2"
