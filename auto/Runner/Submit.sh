#!/bin/bash

input="runs_2.txt"

while read -r config run; do
    [[ -z "$config" ]] && continue
    [[ "$config" =~ ^# ]] && continue

    echo "config = $config, run = $run"

    # ここに処理を書く
    bash Runs.sh "$config" "$run"

done < "$input"