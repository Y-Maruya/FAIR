#!/bin/bash
for run in 22074 22140 22160 22168 22206; do
    ./bin/fair_single config/ped_${run}.yaml -i config/Input/ped_${run}.txt &
done