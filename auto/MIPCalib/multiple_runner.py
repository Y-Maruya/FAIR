import subprocess
import argparse
import sys

parser = argparse.ArgumentParser(description="Run mip-calibration for a range of runs.")
parser.add_argument("--start_run_number", "-s", type=int, required=True,
                    help="The starting run number for mip-calibration.")
parser.add_argument("--end_run_number", "-e", type=int, default=None,
                    help="The ending run number for mip-calibration inclusive.")

args = parser.parse_args()
start_run_number = args.start_run_number
end_run_number = args.end_run_number

workdir = "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/auto/MIPCalib"

while True:
    if end_run_number is not None and start_run_number > end_run_number:
        print("Reached the end run number, exiting.", flush=True)
        break

    if start_run_number > 30000:
        print("Reached run number 30000, exiting.", flush=True)
        break

    print(f"\n=== Processing run {start_run_number} ===", flush=True)

    result = subprocess.run(
        ["bash", "run_mipcalib.sh", str(start_run_number)],
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    # 子プロセスの cout/stdout を Condor stdout に流す
    if result.stdout:
        print(result.stdout, end="", flush=True)

    # 子プロセスの cerr/stderr を Condor stderr に流す
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr, flush=True)

    if result.returncode != 0:
        print(f"run_mipcalib.sh failed with return code {result.returncode}", file=sys.stderr, flush=True)
        sys.exit(result.returncode)

    next_run_number = None

    for line in result.stdout.splitlines():
        if "Next run after processing: " in line:
            print("Found:", line, flush=True)
            next_run_number = int(line.split("Next run after processing: ")[1])
            break

    if next_run_number is None:
        print("Did not find the next run number in the output, exiting.", file=sys.stderr, flush=True)
        sys.exit(1)

    if next_run_number <= start_run_number:
        print("Next run number is not greater than current run number, exiting to avoid infinite loop.",
              file=sys.stderr, flush=True)
        sys.exit(1)

    start_run_number = next_run_number