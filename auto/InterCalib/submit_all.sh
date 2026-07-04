#!/bin/bash

input="run_starts.txt"

while read -r run n; do
    [[ -z "$run" ]] && continue
    [[ "$run" =~ ^# ]] && continue

    echo "run = $run, n = $n"

    bash submit.sh "$run" "$n"

done < "$input"