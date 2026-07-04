#!/bin/bash

# runNumber=$1
password=$1
for file in `ls ./mip_imputed_last/run22400*.json`
do
    echo "uploading $file"
    if [ -z "$password" ]; then
        echo "No password provided. Uploading without authentication."
        python3 ../upload2DB.py $file FASER
    else
        python3 ../upload2DB.py $file FASER $password
    fi
    #python3 ../upload2DB.py $file FASER $password
done