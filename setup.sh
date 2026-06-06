#!/usr/bin/env bash

# Copy this file to setup_local.sh and edit it for your environment.
source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh

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
