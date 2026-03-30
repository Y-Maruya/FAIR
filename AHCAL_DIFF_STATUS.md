# AHCAL-simulation parity status (FAIR)

Date: 2026-03-28

This note captures the current FAIR-side implementation status relative to the
AHCAL-simulation alignment work.

## Implemented in FAIR

- gFaser-derived primary event loading in Geant4 primary generator.
- Interaction/secondary extraction and SimData propagation.
- Birks constant configuration for active scintillator material.
- Trigger component construction in detector geometry.
- Pythia8 external decayer integration hooks and physics registration path.

## Remaining verification focus

- Full detector-geometry parity against AHCAL-simulation reference.
- Full label/category parity for interaction/secondary logic.
- Runtime verification of Pythia8 decayer in FAIR_HAS_PYTHIA8-enabled build.
- Physics-output parity checks (layer Edep, total Edep, timing, secondary spectra).

## Note

The corresponding AHCAL-simulation repository was not present in this execution
environment during this update, so direct file-by-file diff validation was not
performed here.
