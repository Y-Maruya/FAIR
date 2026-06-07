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
DEFAULT_NUM_VETO_EVENTS = 7000 * 18 * 18
MAX_RUN_NUMBER = 30000
MAX_SEARCH_DISTANCE = 1000

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


def get_any_vetoed(run_number, timeout):
    response = fetch_json(
        CALIB_DB_QUERY_URL,
        {
            "CalibrationType": "Triggers",
            "RunNumber": run_number,
            "Layer": -1,
            "PerChannel": "False",
            "PerChip": "False",
        },
        timeout,
    )
    if isinstance(response, list):
        if not response:
            raise RuntimeError(f"No triggers data found for run {run_number}")
        response = response[0]

    try:
        return int(response["Summary"]["AnyVetoed"])
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError(
            f"Invalid triggers data for run {run_number}: missing Summary.AnyVetoed"
        ) from error


def select_run_group(start_run, num_veto_events, excluded_runs, timeout):
    selected_runs = []
    veto_events_count = 0

    for run_number in range(start_run, start_run + MAX_SEARCH_DISTANCE + 1):
        if is_ahcal_run(run_number, excluded_runs, timeout):
            any_vetoed = get_any_vetoed(run_number, timeout)
            veto_events_count += any_vetoed
            selected_runs.append(run_number)
            print(
                f"Run {run_number}: AnyVetoed={any_vetoed}, "
                f"Cumulative Veto Events={veto_events_count}"
            )
        if veto_events_count >= num_veto_events:
            break

    if not selected_runs:
        raise RuntimeError(
            f"No AHCAL runs selected within {MAX_SEARCH_DISTANCE} runs "
            f"after start run {start_run}"
        )
    return selected_runs


def find_start_runs(start_run, end_run, num_veto_events, excluded_runs, timeout):
    start_runs = []
    current_start = start_run

    while current_start <= end_run and current_start <= MAX_RUN_NUMBER:
        print(f"\nSelecting group starting from run {current_start}...")
        selected_runs = select_run_group(
            current_start, num_veto_events, excluded_runs, timeout
        )
        start_runs.append(current_start)

        next_start = selected_runs[-1] + 1
        print(
            f"Group start {current_start}: selected {selected_runs[0]}-"
            f"{selected_runs[-1]}; next start is {next_start}"
        )
        if next_start <= current_start:
            raise RuntimeError("Next start run did not advance")
        current_start = next_start

    return start_runs


def write_start_runs(path, start_runs):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{run_number}\n" for run_number in start_runs))
    print(f"\nWrote {len(start_runs)} Condor job start runs to {path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Select independent ExecCalib start runs for parallel Condor jobs."
    )
    parser.add_argument(
        "--start_run_number",
        "-s",
        type=int,
        required=True,
        help="First start run for SPE calibration.",
    )
    parser.add_argument(
        "--end_run_number",
        "-e",
        type=int,
        default=MAX_RUN_NUMBER,
        help="Last allowed job start run (inclusive).",
    )
    parser.add_argument(
        "--num-veto-events",
        "-n",
        type=int,
        default=DEFAULT_NUM_VETO_EVENTS,
        help="Target AnyVetoed count per ExecCalib job.",
    )
    parser.add_argument(
        "--excluded-runs-file",
        type=Path,
        default=REPO_ROOT / "excluded_runs.txt",
        help="Excluded runs file used by ExecCalib.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=SCRIPT_DIR / "run_starts.txt",
        help="Output file consumed by submit.sdf.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30,
        help="HTTP request timeout in seconds.",
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
    if args.start_run_number <= 0:
        raise ValueError("Start run number must be positive")
    if args.end_run_number < args.start_run_number:
        raise ValueError("End run number must not be smaller than start run number")
    if args.num_veto_events <= 0:
        raise ValueError("Number of veto events must be positive")

    excluded_runs = read_excluded_runs(args.excluded_runs_file)
    start_runs = find_start_runs(
        args.start_run_number,
        args.end_run_number,
        args.num_veto_events,
        excluded_runs,
        args.timeout,
    )
    write_start_runs(args.output, start_runs)

    if args.submit:
        subprocess.run(
            ["condor_submit", str(args.submit_file)],
            cwd=args.submit_file.parent,
            check=True,
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
