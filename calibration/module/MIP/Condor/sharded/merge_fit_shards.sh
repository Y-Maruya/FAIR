#!/bin/bash
set -euo pipefail

output_dir=$1
output_stem=$2
shard_count=$3

parts=()
for ((shard = 0; shard < shard_count; ++shard)); do
    part="${output_dir}/${output_stem}_part_${shard}.root"
    if [[ ! -s "$part" ]]; then
        echo "Missing or empty shard output: $part" >&2
        exit 1
    fi
    parts+=("$part")
done

merged="${output_dir}/${output_stem}.root"
temporary="${merged}.tmp"

rm -f "$temporary"
hadd -f "$temporary" "${parts[@]}"
mv "$temporary" "$merged"

echo "Merged output: $merged"

