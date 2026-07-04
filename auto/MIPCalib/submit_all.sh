#!/bin/bash

input="runs.txt"

while read -r srun erun; do
    [[ -z "$srun" ]] && continue
    [[ "$srun" =~ ^# ]] && continue

    echo "start run = $srun, end run = $erun"

    bash submit.sh "$srun" "$erun"

done < "$input"