#!/bin/bash
set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
    echo "Usage: $0 <input.root> <output_dir> <output_stem> <efficiency_input.root> [shard_count]" >&2
    exit 1
fi

input_file=$(realpath -m "$1")
output_dir=$(realpath -m "$2")
output_stem=$3
efficiency_input_file=$(realpath -m "$4")
shard_count=${5:-5}

if [[ ! "$shard_count" =~ ^[1-9][0-9]*$ ]]; then
    echo "shard_count must be a positive integer: $shard_count" >&2
    exit 1
fi
if [[ ! -r "$input_file" ]]; then
    echo "Cannot read input file: $input_file" >&2
    exit 1
fi
if [[ ! -r "$efficiency_input_file" ]]; then
    echo "Cannot read efficiency input file: $efficiency_input_file" >&2
    exit 1
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "${script_dir}/../../../../.." && pwd)
source_file="${repo_root}/calibration/module/MIP/fit_histograms_noplot_direct_sharded.C"

mkdir -p "$output_dir"
job_dir="${output_dir}/condor_${output_stem}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$job_dir"

executable="${job_dir}/fit_histograms_noplot_direct_sharded.exe"
echo "Building sharded fitter: $executable"
g++ -O3 -std=c++17 "$source_file" $(root-config --cflags --libs) -o "$executable"

fit_submit="${job_dir}/fit.sub"
merge_submit="${job_dir}/merge.sub"
dag_file="${job_dir}/fit.dag"

cat > "$fit_submit" <<EOF
universe = vanilla
executable = ${script_dir}/run_fit_shard.sh
arguments = "'${executable}' '${input_file}' '${output_dir}' '${output_stem}' '${efficiency_input_file}' \$(shard) ${shard_count}"
getenv = True
should_transfer_files = NO
request_cpus = 1
request_memory = 1 GB
+JobFlavour = "longlunch"
+AccountingGroup = "group_u_FASER.users"
output = ${job_dir}/fit_\$(shard).stdout
error = ${job_dir}/fit_\$(shard).stderr
log = ${job_dir}/fit.log
queue 1
EOF

cat > "$merge_submit" <<EOF
universe = vanilla
executable = ${script_dir}/merge_fit_shards.sh
arguments = "'${output_dir}' '${output_stem}' ${shard_count}"
getenv = True
should_transfer_files = NO
request_cpus = 1
request_memory = 1 GB
+JobFlavour = "espresso"
+AccountingGroup = "group_u_FASER.users"
output = ${job_dir}/merge.stdout
error = ${job_dir}/merge.stderr
log = ${job_dir}/merge.log
queue 1
EOF

fit_nodes=()
: > "$dag_file"
for ((shard = 0; shard < shard_count; ++shard)); do
    node="FIT_${shard}"
    fit_nodes+=("$node")
    echo "JOB ${node} ${fit_submit}" >> "$dag_file"
    echo "VARS ${node} shard=\"${shard}\"" >> "$dag_file"
    echo "RETRY ${node} 2" >> "$dag_file"
done
echo "JOB MERGE ${merge_submit}" >> "$dag_file"
echo "RETRY MERGE 1" >> "$dag_file"
echo "PARENT ${fit_nodes[*]} CHILD MERGE" >> "$dag_file"

echo "Submitting ${shard_count} fit jobs followed by merge"
echo "Job directory: $job_dir"
if [[ ${DRY_RUN:-0} == 1 ]]; then
    echo "DRY_RUN=1: DAG generated but not submitted"
else
    condor_submit_dag "$dag_file"
fi
