#!/usr/bin/env bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh

# SET UP FRONTIER SERVER and CORAL LOOKUP TO READ ATLAS CONDITIONS DATABASE
export FRONTIER_SERVER="(proxyurl=http://ca-proxy-atlas.cern.ch:3128)"
export CORAL_DBLOOKUP_PATH=/cvmfs/atlas.cern.ch/repo/sw/software/24.0/Athena/24.0.41/InstallArea/x86_64-el9-gcc13-opt/XML/AtlasAuthentication
export CORAL_AUTH_PATH="$CORAL_DBLOOKUP_PATH"
export COOL_ORA_ENABLE_ADAPTIVE_OPT=Y
echo "$FRONTIER_SERVER"
echo "$CORAL_DBLOOKUP_PATH"
echo "$COOL_ORA_ENABLE_ADAPTIVE_OPT"
ls "$CORAL_DBLOOKUP_PATH/dblookup.xml"
ls "$CORAL_AUTH_PATH/authentication.xml"


# Copy this file to setup_local.sh and edit it for your environment.
# Below is an example of how to set up the environment variables for FAIR. Adjust the paths as needed for your setup.

export FAIR_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Setting up FAIR environment"
echo "FAIR_BASE = ${FAIR_BASE}"

export FAIR_OUT="${FAIR_BASE}/out"
export FAIR_OUT_CALIB="${FAIR_BASE}/out_calibration"
export FAIR_DATA="/eos/experiment/faser/raw/2026/"
export FAIR_ROOT_OUT="/eos/project/f/faser-upgrade/AHCAL/FAIR/"

if [ ! -d "${FAIR_ROOT_OUT}" ]; then
    echo "Creating FAIR_ROOT_OUT: ${FAIR_ROOT_OUT}"
    mkdir -p "${FAIR_ROOT_OUT}"
fi

if [ ! -d "${FAIR_OUT}" ]; then
    echo "Creating FAIR_OUT: ${FAIR_OUT}"
    mkdir -p "${FAIR_OUT}"
fi

if [ ! -d "${FAIR_OUT_CALIB}" ]; then
    echo "Creating FAIR_OUT_CALIB: ${FAIR_OUT_CALIB}"
    mkdir -p "${FAIR_OUT_CALIB}"
fi

if [ ! -d "${FAIR_ROOT_OUT}" ]; then
    echo "Creating FAIR_ROOT_OUT: ${FAIR_ROOT_OUT}"
    mkdir -p "${FAIR_ROOT_OUT}"
fi

echo "FAIR_OUT = ${FAIR_OUT}"
echo "FAIR_OUT_CALIB = ${FAIR_OUT_CALIB}"
echo "FAIR_DATA = ${FAIR_DATA}"
echo "FAIR_ROOT_OUT = ${FAIR_ROOT_OUT}"
