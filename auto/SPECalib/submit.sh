#!/bin/bash
set -eo pipefail

cd /eos/user/y/ymaruya/FASER/AHCAL/FAIR/auto/SPECalib
set +u
source ../../setup.sh

bash run_specalib.sh "$1"
