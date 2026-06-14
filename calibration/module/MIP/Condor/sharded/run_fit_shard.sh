#!/bin/bash
set -euo pipefail

executable=$1
input_file=$2
output_dir=$3
output_stem=$4
efficiency_input_file=$5
shard_index=$6
shard_count=$7

mkdir -p "$output_dir"
output_file="${output_dir}/${output_stem}_part_${shard_index}.root"

echo "Running shard ${shard_index}/${shard_count}"
"$executable" \
    "$input_file" \
    "$output_file" \
    "$efficiency_input_file" \
    "$shard_index" \
    "$shard_count"

