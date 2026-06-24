#!/bin/bash
set -eo pipefail

cd /eos/user/y/ymaruya/FASER/AHCAL/FAIR/auto/MIPCalib
set +u
source ../../setup.sh

bash run_mipcalib.sh "$1" "$2"
