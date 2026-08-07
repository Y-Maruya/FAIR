#!/bin/bash

cdir=$(pwd)
setupLCG_109
dirs=$(ls -d ./*)
echo "dirs: $dirs"
for dir in $dirs; do
    cd $dir
    root -l -b -q '/eos/user/y/ymaruya/FASER/AHCAL/FAIR/calibration/module/spe/refit_spe_from_hists.C("spe_analysis.root","spe_refit_gaus_min11.root",8192,11.0,40.0,500,0.01,4,true)' &
    cd $cdir
done