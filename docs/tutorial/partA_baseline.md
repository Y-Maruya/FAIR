# Part A (45 min): Build and Baseline Run

## Objective
Confirm that you can build the project and run one known-good existing chain.

## Step 1: Build
Use your standard FAIR build flow.

Example (edit to your environment):
```bash
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
./exe/<baseline_executable> <config_or_input>
```

<details>
  <summary>Answer (example run memo)</summary>

```text
Build: cmake .. && make -j
Run: ./exe/ExecCalib config/baseline.yaml
Output: output_baseline.root
```
</details>

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
