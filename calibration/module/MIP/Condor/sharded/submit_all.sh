#!/bin/bash

input="/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/auto/MIPCalib/runs.txt"
module load lxbatch/eossubmit
while read -r srun erun; do
    [[ -z "$srun" ]] && continue
    [[ "$srun" =~ ^# ]] && continue
    
    echo "start run = $srun, end run = $erun"
    cd /eos/user/y/ymaruya/FASER/AHCAL/FAIR/calibration/module/MIP/Condor/sharded/
    # echo "source submit_sharded.sh /eos/user/y/ymaruya/FASER/AHCAL/FAIR/out_calibration/mip_calib/${srun}-${erun}/mip_neighborcheck_nofit.root /eos/user/y/ymaruya/FASER/AHCAL/FAIR/out_calibration/mip_calib/${srun}-${erun}/ mip_neighborcheck_fitted_sher_direct_7 /eos/user/y/ymaruya/FASER/AHCAL/FAIR/out/SimMuon/mip_neighborcheck_nofit.root 5"
    source submit_sharded.sh /eos/user/y/ymaruya/FASER/AHCAL/FAIR/out_calibration/mip_calib/${srun}-${erun}/mip_neighborcheck_nofit_tmp.root /eos/user/y/ymaruya/FASER/AHCAL/FAIR/out_calibration/mip_calib/${srun}-${erun}/ mip_neighborcheck_fitted_sher_direct_8 /eos/user/y/ymaruya/FASER/AHCAL/FAIR/out/SimMuon/mip_neighborcheck_nofit.root 5
    # bash /afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/calibration/module/MIP/Condor/sharded/submit_sharded.sh "$srun" "$erun"

done < "$input"