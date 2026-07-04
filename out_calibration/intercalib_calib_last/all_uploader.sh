#!/bin/bash

# runNumber=$1
password=$1
for file in `ls ./intercalib_corrected/*json`
do
    echo "uploading $file"
    python3 ../upload2DB.py $file FASER $password
done