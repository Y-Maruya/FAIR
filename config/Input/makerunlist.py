import glob 
import os

output_file = "InternalTrigger.txt"
good_run_list = [
    21774,
    21773,
    21739,
    21729
]
# serch path of raw data from /eos/experiment/faser/raw/2026/Faser*-Physics-0"run_number"-{pool_number}.raw
raw_data_path = "/eos/experiment/faser/raw/2026/0{}/Faser*-Physics-0{}-*.raw"
with open(output_file, "w") as f:
    for run_number in good_run_list:
        # search raw data files for the run number
        files = glob.glob(raw_data_path.format(run_number, run_number))
        if len(files) == 0:
            print(f"No raw data files found for run {run_number}")
        else:
            print(f"Raw data files for run {run_number}:")
            for file in files:
                #extract the pool number from the file name
                pool_number = os.path.basename(file).split("-")[-1].replace(".raw", "")
                pool_number = int(pool_number)
                f.write(f"{file} {run_number} {pool_number}\n")
                print(f"  {file} {run_number} {pool_number}")

print(f"Full muon run list has been written to {output_file}")

