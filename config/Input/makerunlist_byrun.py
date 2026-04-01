import glob
import os
import re

# output_file = "pedestal_run.txt"
# NOTE: Use tuples (run_number, output_file). Sets are unordered, and can swap values.
good_run_list = [
    (22509, "muon_run22509.txt")
]

# search path of raw data from /eos/experiment/faser/raw/2026/0"run_number"/Faser*-Physics-0"run_number"-*.raw
raw_data_path = "/eos/experiment/faser/raw/2026/0{}/Faser*-Physics-0{}-*.raw"

written_files = []
for run_number, output_file in good_run_list:
    written_files.append(output_file)
    with open(output_file, "w") as f:
        files = sorted(glob.glob(raw_data_path.format(run_number, run_number)))
        if not files:
            print(f"No raw data files found for run {run_number}")
            continue

        print(f"Raw data files for run {run_number}:")
        for file in files:
            base = os.path.basename(file)
            m = re.search(r"-(\d+)\.raw$", base)
            pool_number = int(m.group(1)) if m else -1
            f.write(f"{file} {run_number} {pool_number}\n")
            print(f"  {file} {run_number} {pool_number}")

print("Run lists written to:")
for name in written_files:
    print(f"  {name}")

