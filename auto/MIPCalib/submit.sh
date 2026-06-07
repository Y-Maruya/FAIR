#!/bin/bash
set -eo pipefail

cd /eos/user/y/ymaruya/FASER/AHCAL/FAIR/auto/MIPCalib
set +u
source ../../setup.sh

flock /tmp/fair_calib_exec.lock bash run_mipcalib.sh "$1"
