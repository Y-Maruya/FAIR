#!/bin/bash
set -eo pipefail

cd /eos/user/y/ymaruya/FASER/AHCAL/FAIR/auto/SPECalib
set +u
source ../../setup.sh
set -u
# source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh

python3 -u /eos/user/y/ymaruya/FASER/AHCAL/FAIR/auto/SPECalib/multiple_runner.py -s "$1" -e "$2"