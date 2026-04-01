#!/bin/bash

runNumber=$1
password=$2
for file in `ls ${runNumber}/json/*json`
do
    echo "uploading $file"
    python3 ../upload2DB.py $file FASER $password
done