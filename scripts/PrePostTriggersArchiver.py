from influxdb import InfluxDBClient
import os
import re
import json
from pathlib import Path
import argparse
from datetime import datetime,timedelta,timezone
# Read command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--verbose", "-v", action="store_true",
                    help="Debug output")

parser.add_argument('--run', default="22000",
                    help="Specify run number to query (default: 22000)")

parser.add_argument("--output", "-o", default="/eos/user/y/ymaruya/prepost_triggers/",
                    help="Specify output directory (default: /eos/user/y/ymaruya/prepost_triggers/)")

args = parser.parse_args()

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

if __name__ == "__main__":
    tstart, tend = find_run_time(args.run)

    if tstart is None:
        print(f"Couldn't find start time for run {args.run}!")
        exit(1)

    if tend is None:
        print(f"Couldn't find end time for run {args.run}!")
        exit(1)

    query_triggerID = f'SELECT first("value") FROM "tlureceiver00-AHCALTLUTriggerID_32bit" WHERE time >= {int(tend*1000)}ms and time <= {int(tend*1000)+10000}ms GROUP BY time(1s) fill(null)'
    query_PreVetoTriggerID = f'SELECT first("value") FROM "tlureceiver00-PreVetoTriggers" WHERE time >= {int(tend*1000)}ms and time <= {int(tend*1000)+10000}ms GROUP BY time(1s) fill(null)'
    query_PostVetoTriggerID = f'SELECT first("value") FROM "tlureceiver00-PostVetoTriggers" WHERE time >= {int(tend*1000)}ms and time <= {int(tend*1000)+10000}ms GROUP BY time(1s) fill(null)'
    if args.verbose: print(f"Querying InfluxDB with: {query_triggerID}")
    result_triggerID = client.query(query_triggerID)
    result_PreVetoTriggerID = client.query(query_PreVetoTriggerID)
    result_PostVetoTriggerID = client.query(query_PostVetoTriggerID)
    points_triggerID = list(result_triggerID.get_points())
    points_PreVetoTriggerID = list(result_PreVetoTriggerID.get_points())
    points_PostVetoTriggerID = list(result_PostVetoTriggerID.get_points())
    first_triggerID = None
    first_PreVetoTriggerID = None
    first_PostVetoTriggerID = None
    last_triggerID = None
    last_PreVetoTriggerID = None
    last_PostVetoTriggerID = None
    for point in points_triggerID:
        if point['first'] is not None:
            if first_triggerID is None:
                first_triggerID = point['first']
            if first_triggerID is not None and point['first'] != first_triggerID:
                break
            last_triggerID = point['first']
    if not last_triggerID.is_integer():
        print(f"Warning: Last TriggerID {last_triggerID} is not an integer. Setting to None.")
        last_triggerID = None
    last_triggerID = int(last_triggerID) if last_triggerID is not None else None
    for point in points_PreVetoTriggerID:
        if point['first'] is not None:
            if first_PreVetoTriggerID is None:
                first_PreVetoTriggerID = point['first']
            if first_PreVetoTriggerID is not None and point['first'] != first_PreVetoTriggerID:
                break
            last_PreVetoTriggerID = point['first']
    if not last_PreVetoTriggerID.is_integer():
        print(f"Warning: Last PreVetoTriggerID {last_PreVetoTriggerID} is not an integer. Setting to None.")
        last_PreVetoTriggerID = None
    last_PreVetoTriggerID = int(last_PreVetoTriggerID) if last_PreVetoTriggerID is not None else None
    for point in points_PostVetoTriggerID:
        if point['first'] is not None:
            if first_PostVetoTriggerID is None:
                first_PostVetoTriggerID = point['first']
            if first_PostVetoTriggerID is not None and point['first'] != first_PostVetoTriggerID:
                break
            last_PostVetoTriggerID = point['first']
    if not last_PostVetoTriggerID.is_integer():
        print(f"Warning: Last PostVetoTriggerID {last_PostVetoTriggerID} is not an integer. Setting to None.")
        last_PostVetoTriggerID = None
    last_PostVetoTriggerID = int(last_PostVetoTriggerID) if last_PostVetoTriggerID is not None else None
    if args.verbose:
        for point in points_triggerID:
            print(f"TriggerID: {point['time']} -> {point['first']}")
        print(f"Last TriggerID: {last_triggerID}")
        for point in points_PreVetoTriggerID:
            print(f"PreVetoTriggerID: {point['time']} -> {point['first']}")
        print(f"Last PreVetoTriggerID: {last_PreVetoTriggerID}")
        for point in points_PostVetoTriggerID:
            print(f"PostVetoTriggerID: {point['time']} -> {point['first']}")
        print(f"Last PostVetoTriggerID: {last_PostVetoTriggerID}")
    
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    with open(args.output+"/run_"+args.run+"_prepost_triggers.txt", "w") as f:
        f.write(f"Run, Start Time, End Time, TriggerID, PreVetoTriggerID, PostVetoTriggerID\n")
        f.write(f"{args.run}, {datetime.fromtimestamp(tstart, tz=timezone.utc).isoformat()}, {datetime.fromtimestamp(tend, tz=timezone.utc).isoformat()}, {last_triggerID}, {last_PreVetoTriggerID}, {last_PostVetoTriggerID}\n")

