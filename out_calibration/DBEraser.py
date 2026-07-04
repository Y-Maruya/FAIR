#!/usr/bin/env python3

import argparse
import json
import subprocess
import sys
from urllib.parse import urlencode


DELETE_URL = "https://ahcalib-calibrationdb.app.cern.ch/Delete"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Delete all calibration DB entries matching CalibrationType in a JSON file."
    )
    parser.add_argument("json_file", help="JSON file containing CalibrationType")
    parser.add_argument("user", help="Calibration DB user")
    parser.add_argument("password", help="Calibration DB password")
    parser.add_argument("--run", type=int, help="Run number to filter entries (optional)")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Query matching entries without deleting them",
    )
    return parser.parse_args()


def load_calibration_type(json_file):
    with open(json_file, "r") as f:
        data = json.load(f)

    calibration_type = data.get("CalibrationType")
    if not calibration_type:
        raise ValueError("CalibrationType is not found in the JSON file")

    return calibration_type


def main():
    args = parse_args()

    try:
        calibration_type = load_calibration_type(args.json_file)
    except (OSError, json.JSONDecodeError, ValueError) as e:
        print(f"failed to read calibration type: {e}", file=sys.stderr)
        return 1

    query = {"CalibrationType": calibration_type}
    if args.run is not None:
        query["RunNumber"] = args.run
    if not args.dry_run:
        query["Delete"] = "forSure"

    url = f"{DELETE_URL}?{urlencode(query)}"
    action = "query" if args.dry_run else "delete"

    print(f"{action}: CalibrationType={calibration_type}")
    print(f"curl -s -u {args.user}:******** '{url}'")

    result = subprocess.run(
        ["curl", "-s", "-u", f"{args.user}:{args.password}", url],
        check=False,
        text=True,
        capture_output=True,
    )

    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)

    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
