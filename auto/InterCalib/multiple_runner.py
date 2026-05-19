import subprocess
import argparse

# Read command-line arguments
parser = argparse.ArgumentParser(description="Run inter-calibration for a range of runs.")
parser.add_argument("--start_run_number", "-s", type=int, help="The starting run number for inter-calibration.")
parser.add_argument("--end_run_number", "-e", type=int, default=None, help="The ending run number for inter-calibration (inclusive). If not provided, it will produce until 30000.")

args = parser.parse_args()
start_run_number = args.start_run_number
end_run_number = args.end_run_number

while True:
    if end_run_number is not None and start_run_number > end_run_number:
        print("Reached the end run number, exiting.")
        break
    if start_run_number > 30000:
        print("Reached run number 30000, exiting.")
        break
    result = subprocess.run(
        ["bash","run_intercalib_muon.sh", f"{start_run_number}"],
        cwd="/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/auto/InterCalib",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    output = result.stdout
    error = result.stderr

    for line in output.splitlines():
        if "Next run after processing: " in line:
            print("Found:", line)
            next_run_number = int(line.split("Next run after processing: ")[1])
            if next_run_number <= start_run_number:
                print("Next run number is not greater than current run number, exiting to avoid infinite loop.")
                exit(1)
            start_run_number = next_run_number
            break
    else:
        print("Did not find the next run number in the output, exiting.")
        print("Output was:", output)
        print("Error was:", error)
        exit(1)
    
