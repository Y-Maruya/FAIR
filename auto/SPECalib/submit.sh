#!/bin/bash
set -eo pipefail

cd /eos/user/y/ymaruya/FASER/AHCAL/FAIR/auto/SPECalib
set +u
source ../../setup.sh

flock /tmp/fair_calib_exec.lock bash run_specalib.sh "$1"
