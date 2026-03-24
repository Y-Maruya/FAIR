import glob 
import os

output_file = "EHN1.txt"
good_run_list = glob.glob("/eos/project/f/faser-upgrade/AHCAL/AHCAL-data-EHN1/Faser*-Physics-0*-*.raw")
good_run_list = [int(os.path.basename(file).split("-")[2][1:]) for file in good_run_list]
good_run_list = sorted(set(good_run_list))
# serch path of raw data from /eos/project/f/faser-upgrade/AHCAL/AHCAL-data-EHN1/Faser*-Physics-0"run_number"-{pool_number}.raw
raw_data_path = "/eos/project/f/faser-upgrade/AHCAL/AHCAL-data-EHN1/Faser*-Physics-0{}-*.raw"
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

