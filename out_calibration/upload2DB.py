import os
import sys
import json

import argparse
import numpy as np


args = sys.argv

json_open = open(args[1],'r') 
json_load = json.load(json_open)
user = args[2]
password = args[3]
print(json_load['RunNumber'])

if json_load['RunNumber'] !=int(0):
    print("upload")
    os.system('curl -s -d @'+str(args[1])+' -H "Content-Type: application/json" -X POST -u ' + user + ':' + password + ' https://ahcalib-calibrationdb.app.cern.ch/AddEntry')
    print('curl -s -d @'+str(args[1])+' -H "Content-Type: application/json" -X POST -u ' + user + ':' + password + ' https://ahcalib-calibrationdb.app.cern.ch/AddEntry')

else:
    print('upload failed')
