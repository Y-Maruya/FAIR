#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fair_dir="$(cd "${script_dir}/../.." && pwd)"
decision_source="${script_dir}/mip_fit_decision_by_run_v2.C"
fit_source="${fair_dir}/calibration/module/MIP/fit_histograms_directfit_integral.C"
decision_exe="${script_dir}/mip_fit_decision_by_run_v2.exe"
fit_exe="${fair_dir}/calibration/module/MIP/fit_histograms_directfit_integral.exe"

base_dir="${script_dir}"
run_list=""
state_shifts=""
condition=""
decision_file="mip_neighborcheck_fitted_sher_direct_3.root"
input_file="mip_neighborcheck_nofit.root"
tree_name="mip_fit_results"
output_dir=""
efficiency_input="${fair_dir}/out/SimMuon/mip_neighborcheck_nofit.root"

usage() {
  cat <<EOF
Usage:
  $(basename "$0") --run RUN --state-shifts SHIFTS --condition EXPR [options]

Required:
  --run RUN                 Run directory, e.g. 21987-22163
  --state-shifts SHIFTS     Comma-separated bit shifts, e.g. 20 or 20,27
  --condition EXPR          TTreeFormula condition, e.g. 'direct_chi2_ndf>3'

Options:
  --base-dir DIR            Directory containing run directories [${base_dir}]
  --decision-file FILE      Fit-result ROOT file [${decision_file}]
  --input-file FILE         No-fit ROOT file used for plotting [${input_file}]
  --tree-name NAME          Decision tree name [${tree_name}]
  --output-dir DIR          Workflow output directory
  --efficiency-input FILE   MC efficiency input ROOT file [${efficiency_input}]
  -h, --help                Show this help

The selected state bits use AND matching. For example, --state-shifts 20,27
selects channels containing both (1<<20) and (1<<27).
EOF
}

while (($#)); do
  case "$1" in
    --run) run_list="$2"; shift 2 ;;
    --state-shifts) state_shifts="$2"; shift 2 ;;
    --condition) condition="$2"; shift 2 ;;
    --base-dir) base_dir="$2"; shift 2 ;;
    --decision-file) decision_file="$2"; shift 2 ;;
    --input-file) input_file="$2"; shift 2 ;;
    --tree-name) tree_name="$2"; shift 2 ;;
    --output-dir) output_dir="$2"; shift 2 ;;
    --efficiency-input) efficiency_input="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [[ -z "$run_list" || -z "$state_shifts" || -z "$condition" ]]; then
  usage >&2
  exit 1
fi

mask=0
IFS=',' read -ra shifts <<< "$state_shifts"
for shift_value in "${shifts[@]}"; do
  shift_value="${shift_value//[[:space:]]/}"
  if [[ ! "$shift_value" =~ ^[0-9]+$ ]] || ((shift_value < 0 || shift_value > 30)); then
    echo "Invalid state shift: $shift_value (expected 0-30)" >&2
    exit 1
  fi
  ((mask |= 1 << shift_value))
done

safe_condition="$(printf '%s' "$condition" | tr -cs '[:alnum:]' '_' | sed 's/^_//;s/_$//')"
safe_runs="$(printf '%s' "$run_list" | tr -cs '[:alnum:]-' '_')"
if [[ -z "$output_dir" ]]; then
  output_dir="${base_dir}/selected_direct_integral_${safe_runs}_mask_${mask}_${safe_condition}"
fi
selection_dir="${output_dir}/selection"
plot_dir="${output_dir}/plots"
fit_output="${output_dir}/selected_direct_integral_fitted.root"
mkdir -p "$selection_dir" "$plot_dir"

compile_if_needed() {
  local source="$1"
  local executable="$2"
  if [[ ! -x "$executable" || "$source" -nt "$executable" ]]; then
    echo "Compiling $(basename "$source")"
    g++ -O2 -std=c++17 "$source" $(root-config --cflags --libs) -o "$executable"
  fi
}

compile_if_needed "$decision_source" "$decision_exe"
compile_if_needed "$fit_source" "$fit_exe"

echo "State shifts : $state_shifts"
echo "State mask   : $mask"
echo "Condition    : $condition"
echo "Output       : $output_dir"

"$decision_exe" \
  --baseDir "$base_dir" \
  --runList "$run_list" \
  --fileName "$decision_file" \
  --treeName "$tree_name" \
  --outDir "$selection_dir" \
  --selectStateMask "$mask" \
  --selectCondition "$condition" \
  --exportOnly

cellid_file="$(find "$selection_dir" -maxdepth 1 -name '*_cellids.txt' -print -quit)"
if [[ -z "$cellid_file" ]]; then
  echo "Cellid output was not created" >&2
  exit 1
fi
if grep -q 'cellid={};' "$cellid_file"; then
  echo "No channels matched the state mask and condition. See: $cellid_file"
  exit 0
fi

if [[ "$run_list" == *,* ]]; then
  echo "Direct integral plotting requires one run directory; got: $run_list" >&2
  exit 1
fi
plot_input="${base_dir}/${run_list}/${input_file}"
"$fit_exe" "$plot_input" "$fit_output" "$cellid_file" "$plot_dir" "$efficiency_input"

echo
echo "Done."
echo "Cellids : $cellid_file"
echo "Details : ${cellid_file%_cellids.txt}_details.txt"
echo "Gallery : ${plot_dir}/index.html"
