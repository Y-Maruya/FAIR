#!/usr/bin/env python3

import argparse
import json
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen


RUNINFO_URL = "https://faser-runinfo.app.cern.ch/cgibin/getRunInfo.py"
CALIB_DB_QUERY_URL = "https://ahcalib-calibrationdb.app.cern.ch/Query"

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent

def fetch_json(url, params, timeout):
    request_url = f"{url}?{urlencode(params)}"
    with urlopen(request_url, timeout=timeout) as response:
        return json.load(response)


def read_excluded_runs(path):
    excluded_runs = set()
    with path.open() as excluded_file:
        for line in excluded_file:
            try:
                excluded_runs.add(int(line.split()[0]))
            except (IndexError, ValueError):
                continue
    return excluded_runs


def is_ahcal_run(run_number, excluded_runs, timeout):
    if run_number in excluded_runs:
        print(f"Run {run_number} is excluded.")
        return False

    run_info = fetch_json(RUNINFO_URL, {"runno": run_number}, timeout)
    if run_info.get("runnumber") != run_number:
        raise RuntimeError(
            f"Runinfo mismatch for run {run_number}: "
            f"got runnumber={run_info.get('runnumber')!r}"
        )
    return run_info.get("type") == "AHCAL"

def write_config_and_runs_to_file(path, config, runs):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{config} {run}\n" for run in runs))
    print(f"Wrote {len(runs)} runs to {path}")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Check if runs are AHCAL runs and write them to a file."
    )
    parser.add_argument(
        "config",
        type=str,
        help="Configuration name to write to the output file.",
    )
    parser.add_argument(
        "--start-run_number",
        "-s",
        type=int,
        required=True,
        help="Starting run number to check.",
    )
    parser.add_argument(
        "--end-run_number",
        "-e",
        type=int,
        required=True,
        help="Ending run number to check.",
    )
    parser.add_argument(
        "--excluded-runs-file",
        type=Path,
        default=REPO_ROOT / "excluded_runs.txt",
        help="Excluded runs file used by ExecCalib.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=10,
        help="Timeout in seconds for HTTP requests.",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Run condor_submit after writing the start-run list.",
    )
    parser.add_argument(
        "--submit-file",
        type=Path,
        default=SCRIPT_DIR / "submit.sdf",
        help="Condor submit description used with --submit.",
    )
    return parser.parse_args()

def main():
    args = parse_args()
    excluded_runs = read_excluded_runs(args.excluded_runs_file)
    ahcal_runs = []

    for run_number in range(args.start_run_number, args.end_run_number + 1):
        if is_ahcal_run(run_number, excluded_runs, args.timeout):
            ahcal_runs.append(run_number)

    write_config_and_runs_to_file(
        Path(f"runs.txt"), args.config, ahcal_runs
    )
    mkdir_path = Path(f"condor_logs/out_{args.config}")
    mkdir_path.mkdir(parents=True, exist_ok=True)
    
    if args.submit:
        subprocess.run(
            ["condor_submit", str(args.submit_file)],
            cwd=args.submit_file.parent,
            check=True,
        )

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)