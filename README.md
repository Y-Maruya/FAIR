# FAIR

FAIR (FASER AHCAL Integrated Reconstruction framework) is the offline analysis and reconstruction framework for the FASER AHCAL project. The codebase integrates common HEP dependencies (ROOT, yaml-cpp, spdlog) to read RawHit/Hit data and run reconstruction algorithms in a modular pipeline.

## Introduction
FAIR is designed to be flexible and modular, allowing users to easily configure and run different reconstruction algorithms through YAML configuration files. Also, it supports writing all of the output from the module in ROOT format for further analysis, and reading all of the data from ROOT files or binary raw files.

## Repository overview

### Top-level layout

- `common/` -- Core framework utilities (algorithm interfaces, registry, event store, geometry).
- `IO/` -- I/O abstractions and concrete readers/writers.
- `calibration/` -- Calibration data and logic; reference files live under `reffiles/`.
- `adc_to_energy/` -- ADC-to-energy conversion algorithms (ROOT TTree input).
- `reco_alg/` -- Reconstruction algorithms and modules.
- `simulation/` -- Simulation and digitization.
- `exe/` -- Executable entry points (drivers such as `MultiInOut.cpp`).
    - `MultiInOut.cpp` -- Main driver that sets up the pipeline based on YAML config. Multi inputs files and multi outputs.
    - `MultiInOneOut.cpp` -- Similar to `MultiInOut`, but with a single output file from several input files.
    - `EventDisplay.cpp` -- Simple event display application from RecoHits.
- `config/` -- YAML configuration files and input lists.
    - `first.yaml` -- Example configuration file demonstrating a full reconstruction chain.
- `scripts/` -- Helper scripts (e.g., submission utilities).
- `faser-common/` -- For decoding the raw files.

## Setup

Example using an LCG view on CERN LXPLUS:
```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh
```

## Downdload Build

```bash
git clone --recursive <repository_url>
```

## Build
```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh 
cd FAIR
mkdir build && cd build
cmake ../
make -j16
make install
```

## Run
Multi input, single output example:
```bash
./bin/fair_single config/first.yaml -i <INPUT_FILE.txt> 
```
Multi input, multi output example:
```bash
./bin/fair_multi config/first.yaml -i <INPUT_FILE.txt> 
```

## Configuration
```yaml
run:
  input: <INPUT_FILE> # Input file this configuration is used only when not using -i option
  MC: false           # Whether the input is MC or real data
  log_file: <LOG_FILE> # Log output file path
  output: <OUTPUT_FILE> # Output ROOT file path
  nEvents: -1         # Number of events to process (-1 for all)
  log_level: INFO     # Log level (TRACE, DEBUG, INFO, WARN, ERROR, CRITICAL)
  runNumber: 0        # Run number to be used in the event, this config is used only when not using -i option
  poolIndex: 0       # Pool index to be used in the event, this config is used only when not using -i option
reader:                # Input module configuration
    type: <READER_TYPE>  # Type of the input module (e.g., RootRawHitReader, BinaryRawHitReader, RootInput)
    cfg:
      <PARAMETER>: <VALUE>
algs:                # List of algorithms to run in sequence
  - type: <ALGO_TYPE>    # Type of the algorithm
    cfg:                # Configuration for the algorithm
      <PARAMETER>: <VALUE>
  - type: <ALGO_TYPE>
    cfg:
      <PARAMETER>: <VALUE>
  - type: RootWriterAlg
    cfg: 
      outputlist:
        # - vector<AHCALRecoHit>
        - SimpleFittedTrack
        - Track
```

## Modules
### Input and output schemes
The framework supports the data flow between the modules through an event store(`common/EventStore.hpp`), which allows modules to read and write data products using string keys. The input modules read data from various formats (ROOT files, binary raw files) and populate the event store. The modules in the processing chain read data products from the event store and perform their respective computations, and the flow of the modules are defined in the YAML configuration file. The modules flow starts from the input module, followed by a sequence of processing algorithms as according to the order written in the configuration file.
Finally, while the output modules write specified data products from the event store to output files. You can specify which data products to write in the configuration of the output module.
### Input modules

- `RootRawHitReader` -- Reads RawHit and TLU data from ROOT files. 
    - Parameters:
        - `out_rawhits_key` -- Key for the output RawHits collection.
        - `out_tlu_key` -- Key for the output TLU data.
- `BinaryRawHitReader` -- Reads RawHit and TLU data from binary raw files.
    - Parameters:
        - `out_rawhits_key` -- Key for the output RawHits collection.
        - `out_tlu_key` -- Key for the output TLU data.
- `RootInput` -- Generic ROOT input module for reading the data from this framework's output files.
    - Parameters:
        - `inputlist` -- List of data products to read from the input ROOT file. 
        Add the following format for each data product:
        [`<C++_TYPE>`, `<KEY>`]
        ``` yaml
        inputlist: 
            - [vector<AHCALRawHit>, RawHits]
            - [AHCALTLURawData, TLURawData]
            - [vector<AHCALRecoHit>, RecoHits]
            - [SimpleFittedTrack, FittedTrack]
            - [Track, MuonKFTrack]
        ```

### Output modules
- `RootWriterAlg` -- Writes specified data products to ROOT files.
    - Parameters:
        - `outputlist` -- List of data products to write to the output ROOT file. 
        Add the following format for each data product:
        `<C++_TYPE>`
        ``` yaml
        outputlist:
            - vector<AHCALRecoHit>
            - SimpleFittedTrack
            - Track
        ```

### Algorithms
- `PedestalAlg` -- Pedestal calculation algorithm. implemented in ```calibration/module/pedestal/PedestalAlg.hpp```.
    - Parameters:
        - `in_rawhits_key` -- Key for the input RawHits collection.
        - `pedestal_to_file` -- Boolean flag to save the calculated pedestal to a file.
        - `pedestal_to_DB` -- Boolean flag to save the calculated pedestal to the database. _NOT IMPLMENTED YET_
        - `out_pedestal_filename` -- Output filename for the pedestal data.
        - `use_hittag` -- Boolean flag to use hittag information for pedestal calculation. (normally true for real data)
        - `selected_hittag` -- Hittag value to select for pedestal calculation. (Default: 0)
    
- `AdcToEnergyReadTTreeAlg` -- ADC-to-energy conversion algorithm using calibration data from ROOT TTrees. implemented in `adc_to_energy/AdcToEnergyReadTTreeAlg.hpp`.
    - Parameters:
        - `in_rawhits_key` -- Key for the input RawHits collection.
        - `out_recohits_key` -- Key for the output RecoHits collection.
        - `mip` 
           - `file` -- Path to the MIP calibration ROOT file.
           - `tree` -- Name of the TTree containing MIP calibration data.
           - `cellid_version` -- Layer position version for cell ID decoding. (Default: 1) 0 for old EHN1 style, 1 for UJ12 style.
        - `pedestal` 
           - `file` -- Path to the pedestal calibration ROOT file.
           - `tree` -- Name of the TTree containing pedestal calibration data.
           - `cellid_version` -- Layer position version for cell ID decoding. (Default: 1) 0 for old EHN1 style, 1 for UJ12 style.
        - `dac` 
           - `file` -- Path to the DAC calibration ROOT file.
           - `tree` -- Name of the TTree containing DAC calibration data.
           - `cellid_version` -- Layer position version for cell ID decoding. (Default: 1) 0 for old EHN1 style, 1 for UJ12 style.
- `TrackFitAlg` -- Simple track fitting algorithm using linear regression. implemented in `reco_alg/module/TrackFitAlg.hpp`.
    - Parameters:
        - `in_recohits_key` -- Key for the input RecoHits collection.
        - `out_track_key` -- Key for the output Track object.
        - `threshold_xy` -- The threshold in the XY plane to consider hits is in the track. (Default: 20.0 mm)
- `MuonKFAlg` -- Kalman filter-based muon track fitting algorithm. implemented in `reco_alg/module/MuonKFAlg.hpp`.
    - NOT CHECKED YET


    

