#!/bin/bash
set -eo pipefail

cd /eos/user/y/ymaruya/FASER/AHCAL/FAIR/auto/InterCalib
set +u
source ../../setup.sh

bash run_intercalib_mip.sh "$1" "$2"
