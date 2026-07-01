import json
import os
import sys
import requests

def read_json_from_web(runnumber):
    url = f"https://faser-runinfo.app.cern.ch/cgibin/getRunInfo.py?runno={runnumber}"
    try:        
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data from {url}: {e}")
        return None
    
def extract_run_type(run_info):
    if run_info and 'type' in run_info:
        return run_info['type']
    return None

def extract_configName(run_info):
    if run_info and 'configName' in run_info:
        return run_info['configName']
    return None

def last_run_number():
    with open("/eos/user/y/ymaruya/FASER/AHCAL/FAIR/auto/pedestal/ped_finished_runs.txt", "r") as f:
        finished_runs = f.read().splitlines()
    return finished_runs[-1] if finished_runs else "21985"

def isFinished(runnumber):
    with open("/eos/user/y/ymaruya/FASER/AHCAL/FAIR/auto/pedestal/ped_finished_runs.txt", "r") as f:
        finished_runs = f.read().splitlines()
    return str(runnumber) in finished_runs

def writeToFinishedRuns(runnumber):
    with open("/eos/user/y/ymaruya/FASER/AHCAL/FAIR/auto/pedestal/ped_finished_runs.txt", "a") as f:
        f.write(f"{runnumber}\n")

if __name__ == "__main__":
    cont_sum_notrun = 0
    for runnumber in range(int(last_run_number()) + 1, 30000):
        print(f"Checking run {runnumber}...")
        if isFinished(runnumber):
            continue
        run_info = read_json_from_web(runnumber)
        if not run_info:
            cont_sum_notrun += 1
            if cont_sum_notrun > 150:
                print("No more runs found, exiting.")
                break
            continue
        run_type = extract_run_type(run_info)
        config_name = extract_configName(run_info)
        if run_type != "AHCALCalib" or config_name != "AHCALCalibration_pedestal":
            continue
        print(f"Processing run {runnumber} with type {run_type} and config {config_name}")
        import subprocess
        subprocess.run(
            ["bash", "run_pedestalcalib.sh", str(runnumber)],
            cwd="/eos/user/y/ymaruya/FASER/AHCAL/FAIR/auto/pedestal",
            check=True
        )
        writeToFinishedRuns(runnumber)

