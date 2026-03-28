#!/usr/bin/env bash
set -euo pipefail

LCG_SETUP="/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh"

if [[ ! -f "$LCG_SETUP" ]]; then
  echo "[ERROR] LCG setup not found: $LCG_SETUP"
  echo "- /cvmfs is likely not mounted in this environment."
  echo "- Please run on a host with CVMFS enabled (e.g. CERN lxplus or a CVMFS-enabled container)."
  exit 1
fi

# shellcheck disable=SC1090
source "$LCG_SETUP"

echo "[OK] Loaded LCG 105 environment"
echo "ROOTSYS=${ROOTSYS:-<unset>}"

if command -v geant4-config >/dev/null 2>&1; then
  echo "Geant4: $(geant4-config --version)"
else
  echo "[WARN] geant4-config not found in PATH"
fi

if command -v root-config >/dev/null 2>&1; then
  echo "ROOT: $(root-config --version)"
else
  echo "[WARN] root-config not found in PATH"
fi

if command -v pythia8-config >/dev/null 2>&1; then
  echo "Pythia8: $(pythia8-config --version 2>/dev/null || echo available)"
else
  echo "[WARN] pythia8-config not found in PATH"
fi
