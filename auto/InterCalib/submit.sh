#!/bin/bash
set -eo pipefail

cd /eos/user/y/ymaruya/FASER/AHCAL/FAIR/auto/InterCalib
set +u
source ../../setup.sh

flock /tmp/fair_calib_exec.lock bash run_intercalib_mip.sh "$1"
