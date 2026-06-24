#!/bin/bash
# set -eo pipefail

cd /eos/user/y/ymaruya/FASER/AHCAL/FAIR/
# set +u
source ./setup.sh

./bin/fair_run "$1" -r "$2"

echo "Finished running $1 for run $2"