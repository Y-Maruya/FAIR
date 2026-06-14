#!/usr/bin/env python3

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


RUN_LIST_URL = "https://faser-runinfo.app.cern.ch/cgibin/getRunList.py"
RUN_INFO_URL = "https://faser-runinfo.app.cern.ch/cgibin/getRunInfo.py"
SNAPSHOT_VERSION = 1


def default_output_path():
    fair_base = os.environ.get("FAIR_BASE")
    if fair_base:
        return Path(fair_base) / "config/local/runinfo_snapshot.json"
    return Path(__file__).resolve().parent.parent / "config/local/runinfo_snapshot.json"


def fetch_json(url, params=None, timeout=30.0, attempts=5, retry_delay=1.0):
    request_url = url
    if params:
        request_url = f"{url}?{urlencode(params)}"

    for attempt in range(1, attempts + 1):
        try:
            request = Request(request_url, headers={"User-Agent": "FAIR-runinfo-snapshot/1"})
            with urlopen(request, timeout=timeout) as response:
                return json.load(response)
        except (HTTPError, URLError, TimeoutError, OSError, json.JSONDecodeError) as error:
            if attempt == attempts:
                raise RuntimeError(
                    f"Failed to fetch {request_url} after {attempts} attempts: {error}"
                ) from error
            delay = retry_delay * attempt
            print(
                f"Transient failure fetching {request_url} "
                f"(attempt {attempt}/{attempts}): {error}. Retrying in {delay:g}s.",
                file=sys.stderr,
            )
            time.sleep(delay)


def load_existing_details(path):
    if not path.exists():
        return {}
    try:
        snapshot = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        print(f"Ignoring unreadable existing snapshot {path}: {error}", file=sys.stderr)
        return {}

    details = snapshot.get("details", {})
    if not isinstance(details, dict):
        return {}
    valid_details = {}
    for run_number, detail in details.items():
        try:
            numeric_run_number = int(run_number)
        except (TypeError, ValueError):
            continue
        if isinstance(detail, dict) and detail.get("runnumber") == numeric_run_number:
            valid_details[str(numeric_run_number)] = detail
    return valid_details


def write_snapshot(path, run_types, details):
    snapshot = {
        "version": SNAPSHOT_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "run_types": run_types,
        "details": details,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_name(f".{path.name}.tmp")
    with temporary_path.open("w") as output:
        json.dump(snapshot, output, sort_keys=True, separators=(",", ":"))
        output.write("\n")
        output.flush()
        os.fsync(output.fileno())
    os.replace(temporary_path, path)


def build_snapshot(output, timeout, attempts, retry_delay):
    print(f"Fetching full run list from {RUN_LIST_URL}...")
    run_list = fetch_json(
        RUN_LIST_URL, timeout=timeout, attempts=attempts, retry_delay=retry_delay
    )
    if not isinstance(run_list, list):
        raise RuntimeError("Run list response is not a JSON array")

    run_types = {}
    ahcal_runs = []
    for run in run_list:
        if not isinstance(run, dict) or not isinstance(run.get("runnumber"), int):
            continue
        run_number = run["runnumber"]
        run_type = run.get("type")
        if isinstance(run_type, str):
            run_types[str(run_number)] = run_type
        if run_type == "AHCAL":
            ahcal_runs.append(run_number)

    existing_details = load_existing_details(output)
    details = {
        str(run_number): existing_details[str(run_number)]
        for run_number in ahcal_runs
        if str(run_number) in existing_details
    }
    write_snapshot(output, run_types, details)

    missing_runs = [run for run in ahcal_runs if str(run) not in details]
    print(
        f"Found {len(run_types)} runs and {len(ahcal_runs)} AHCAL runs; "
        f"reusing {len(details)} details and fetching {len(missing_runs)}."
    )

    for index, run_number in enumerate(missing_runs, start=1):
        detail = fetch_json(
            RUN_INFO_URL,
            {"runno": run_number},
            timeout=timeout,
            attempts=attempts,
            retry_delay=retry_delay,
        )
        if not isinstance(detail, dict) or detail.get("runnumber") != run_number:
            raise RuntimeError(
                f"RunInfo mismatch for run {run_number}: "
                f"got runnumber={getattr(detail, 'get', lambda *_: None)('runnumber')!r}"
            )
        details[str(run_number)] = detail
        write_snapshot(output, run_types, details)
        print(f"[{index}/{len(missing_runs)}] Saved RunInfo for run {run_number}")

    print(f"Wrote complete RunInfo snapshot to {output}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a resumable local RunInfo snapshot for all AHCAL runs."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=default_output_path(),
        help="Snapshot output path (default: ${FAIR_BASE}/config/local/runinfo_snapshot.json).",
    )
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--attempts", type=int, default=5)
    parser.add_argument("--retry-delay", type=float, default=1.0)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.timeout <= 0 or args.attempts <= 0 or args.retry_delay < 0:
        raise ValueError("timeout and attempts must be positive; retry-delay must be non-negative")
    build_snapshot(args.output, args.timeout, args.attempts, args.retry_delay)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
