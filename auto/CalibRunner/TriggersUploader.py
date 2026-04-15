#count the event which has veto signal
from influxdb import InfluxDBClient
import os
import re
import json
from pathlib import Path
import argparse
from datetime import datetime,timedelta,timezone
# Read command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("password", help="Password for CalibrationDB")
parser.add_argument("--verbose", "-v", action="store_true",
                    help="Debug output")

args = parser.parse_args()
password = args.password
def replace_environment_variables(input_str):
    base_dir = Path(__file__).resolve().parent


    # get secret file path
    secret_file = "faser-secret.json"

    env_config = {}
    if Path(secret_file).exists():
        with open(secret_file) as f:
            env_config = json.load(f)
    else:
        print(f"Warning: Secret file {secret_file} not found. Will rely on environment variables only.")
    def replace(match):
        var_name = match.group(1)

        # priority 1: secret file
        if var_name in env_config:
            return str(env_config[var_name])

        # priority 2: environment variable
        env_value = os.getenv(var_name)
        if env_value is not None:
            return env_value

        # fallback
        return match.group(0)

    return re.sub(r"\$([A-Za-z0-9_]+)", replace, input_str)

def find_run_time(run):

    import json
    import requests

    # Make sure this is an int
    runnum = int(run)

    response = requests.get(f"https://faser-runinfo.app.cern.ch/cgibin/getRunInfo.py?runno={runnum}")
    run_json = response.json()

    if not run_json:
        print(f"Couldn't find run information for {runnum}!")
        return None, None
    run_type = run_json.get("type", "unknown")
    start_time = run_json.get("starttime", None)
    end_time = run_json.get("stoptime", None)
    if run_type != "AHCAL":
        print(f"Warning: Run {runnum} is of type {run_type}, not AHCAL. Results may be unreliable.")
        return None, None
    if not start_time:
        tstart = None
    else:
        tstart = datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc)

    if not end_time:
        tend = datetime.now(timezone.utc).timestamp()
    else:
        tend = datetime.strptime(end_time, "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc)

    return tstart.timestamp(), tend.timestamp()


client = InfluxDBClient(
    host='dbod-faser-influx-prod.cern.ch',  # Just the hostname
    port=8080,
    username=replace_environment_variables("$INFLUXUSER"),
    password=replace_environment_variables("$INFLUXPW"),
    database=replace_environment_variables("$INFLUXDB"),
    ssl=True,  # Add this to use HTTPS
    verify_ssl=True
)
def get_value_at_end_time(valuename, tend):
    query = f'SELECT first("value") FROM "{valuename}" WHERE time >= {int(tend*1000)}ms and time <= {int(tend*1000)+10000}ms GROUP BY time(1s) fill(null)'
    if args.verbose: print(f"Querying InfluxDB with: {query}")
    result = client.query(query)
    points = list(result.get_points())
    last_value = None
    for i in range(1, len(points)):
        if points[i-1]['first'] is None:
            points[i-1]['first'] = 0
        if points[i]['first'] is None:
            points[i]['first'] = 0
        if points[i]['first'] is not None and points[i-1]['first'] is not None:
            if points[i]['first'] - points[i-1]['first'] < 0:
                last_value = points[i-1]['first']
                break
        if i == len(points) - 1:
            last_value = points[i]['first']
    if not isinstance(last_value, (int, float)):
        print(f"Warning: Last value {last_value} for {valuename} is not a number. Setting to None.")
        last_value = None
        for point in points:
            print(f"{valuename}: {point['time']} -> {point['first']}")
        print(f"Last {valuename}: {last_value}")
    last_value = int(last_value) if last_value is not None else None

    if args.verbose:
        for point in points:
            print(f"{valuename}: {point['time']} -> {point['first']}")
        print(f"Last {valuename}: {last_value}")

    return last_value

def get_value_of_run(valuename, run):
    tstart, tend = find_run_time(run)
    if tstart is None or tend is None:
        print(f"Couldn't find valid start or end time for run {run}!")
        return -1
    tmp = get_value_at_end_time(valuename, tend)
    if tmp is None:
        return -1
    return tmp

def get_all_runs_from_runinfo():
    import json
    import requests
    try:
        with open("all_runs_checked.txt", "r") as f:
            checked_runs = set(int(line.strip()) for line in f if line.strip().isdigit())
    except FileNotFoundError:
        checked_runs = set([20000])  # Start from 20000 if no file exists
    all_runs = checked_runs.copy()
    start_run = max(checked_runs) if checked_runs else 20000
    for i in range(start_run + 1, 30000):
        if i % 100 == 0:
            print(f"Checking run {i}...")
        response = requests.get(f"https://faser-runinfo.app.cern.ch/cgibin/getRunInfo.py?runno={i}")
        run_json = response.json()
        if not run_json:
            break
        elif run_json["type"] == "AHCAL":
            all_runs.add(i)
            with open("all_runs_checked.txt", "a") as f:
                f.write(f"{i}\n")
    return all_runs

if __name__ == "__main__":
    try:
        with open("finished_runs.txt", "r") as f: # read and if not exist create finished_runs.txt
            finished_runs = set(int(line.strip()) for line in f if line.strip().isdigit())
    except FileNotFoundError:
        finished_runs = set()
    all_runs = get_all_runs_from_runinfo()
    not_finished_runs = [run for run in all_runs if run not in finished_runs]
    print(f"Found {len(all_runs)} total runs, {len(not_finished_runs)} of which are not finished.")
    for run in not_finished_runs:
        print(f"Processing run {run}...")
        any_vetoed = get_value_of_run("tlureceiver00-TriggerVetoed", run)
        non_vetoed = get_value_of_run("tlureceiver00-Trigger_i5_0_i6_0", run)
        non_physical = get_value_of_run("ahcaleventreceiver00-DuplicatedEventCount", run)
        physical = get_value_of_run("ahcaleventreceiver00-GoodEventsCount", run)
        PreVeto = get_value_of_run("tlureceiver00-PreVetoTriggers", run)
        PostVeto = get_value_of_run("tlureceiver00-PostVetoTriggers", run)
        Recorded = get_value_of_run("ahcaleventreceiver00-EventNumber", run)
        #make json file for each run
        tstart, tend = find_run_time(run)
        timestamp = int((tstart + tend) / 2)
        timestamp_str = datetime.fromtimestamp(timestamp, tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
        run_data = {
            "RunNumber": run,
            "TimeStamp": timestamp_str,
            "Layer": -1,
            "CalibrationType": "Triggers",
            "Status": 0,
            "Summary": {
                "AnyVetoed": any_vetoed,
                "NonVetoed": non_vetoed,
                "NonPhysical": non_physical,
                "Physical": physical,
                "PreVeto": PreVeto,
                "PostVeto": PostVeto,
                "Recorded": Recorded
            }
        }
        json_file = f"Triggers/run_{run}_triggers.json"
        with open(json_file, "w") as f:
            json.dump(run_data, f, indent=4)
        print("upload")
        os.system('curl -s -d @'+json_file+' -H "Content-Type: application/json" -X POST -u FASER:' + password + ' https://ahcalib-calibrationdb.app.cern.ch/AddEntry')
        print('curl -s -d @'+json_file+' -H "Content-Type: application/json" -X POST -u FASER:' + password + ' https://ahcalib-calibrationdb.app.cern.ch/AddEntry')
        with open("finished_runs.txt", "a") as f:
            f.write(f"{run}\n")
    print("All done!")