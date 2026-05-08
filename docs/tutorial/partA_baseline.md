# Part A (45 min): Build and Baseline Run

## Objective
Confirm that you can build the project and run one known-good existing chain.

## Step 1: Build
Use your standard FAIR build flow.

Example (edit to your environment):
```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_105a_cuda/x86_64-el9-gcc11-opt/setup.sh
mkdir -p build
cd build
cmake ..
make -j
```

<details>
  <summary>Hint: if build fails</summary>

- Start from a clean build directory.
- Check dependency path variables.
- Compare with existing CI or teammate build command.
</details>

## Step 2: Run baseline executable
Run one known-good command prepared by the instructor.

Example template:
```bash
./bin/fair_single config/muon_run22106.yaml
```


## Step 3: Confirm output
Check:
- log contains normal end/finish message
- expected output file is created
- output file size is non-zero

## Deliverable
Record in your notebook/text file:
- build command
- baseline run command
- output path
