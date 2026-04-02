# AHCAL-simulation parity status (FAIR)

Date: 2026-03-30

This document explicitly lists **known remaining differences** between FAIR's
current internal Geant4 implementation and the AHCAL-simulation (`2026`) model.

Reference checked:
- `AHCAL-simulation/src/DetectorConstruction.cc` (branch `2026`)

## Confirmed aligned (implemented in FAIR)

- `gFaser`-based primary reading path and `SimData` injection fields.
- Birks attenuation handling in SD hit accumulation.
- CaloUnit-like nested geometry concept (housing/passive/sensitive/attach).
- Optional trigger planes and optional Pythia8 decayer wiring.

## Known remaining differences (explicit)

1. **Cell-level geometry granularity**
   - AHCAL-simulation builds repeated `CaloUnitVolume` units at finer granularity.
   - FAIR currently builds one large unit per layer and derives `(chip,channel)` by
     inverse XY mapping in SD.
   - Impact: boundary/edge behavior and copy-number semantics can differ.

2. **Optical surface implementation**
   - AHCAL-simulation includes optical/visual infrastructure (`G4OpticalSurface`,
     `G4LogicalSkinSurface`, `G4LogicalBorderSurface`, material property tables).
   - FAIR currently does not define optical surfaces/material optical properties.
   - Impact: optical-photon-level response is not equivalent.

3. **Material selection path**
   - AHCAL-simulation uses material indices (`GetCaloMaterial(...)`) via detector
     parameter sets.
   - FAIR currently uses direct NIST materials (e.g. `G4_MYLAR`,
     `G4_POLYETHYLENE`, polystyrene density override).
   - Impact: if source indices map to custom compounds, dE/dx and shower response
     can differ.

4. **`Sensitive_dig_out(_ESR)` semantics**
   - FAIR added dig-out style air volumes with source-like dimensions.
   - Exact source material/index behavior and all placement semantics are not yet
     guaranteed bit-identical.

5. **Trigger volume naming/copy-number conventions**
   - AHCAL-simulation SD code distinguishes specific physical names
     (`HCALtriggerPhysicalFront`, `HCALDownstreamPhysical`).
   - FAIR currently keeps compatibility checks in SD, but detector-side naming and
     hierarchy are not fully source-identical.

6. **Output structure parity scope**
   - FAIR writes `AHCALSimHit` + `SimData` (framework EDM schema).
   - AHCAL-simulation writes additional ROOT vectors/truth structures in
     `EventAction/RunAction` not fully mirrored in FAIR EDM.

## Files where these differences currently live

- FAIR detector geometry: `simulation/geant4/DetectorConstruction.cpp`
- FAIR sensitive detector logic: `simulation/geant4/HcalSD.cpp`
- FAIR geometry/runtime config: `simulation/geant4/Geant4Config.hpp`
- FAIR reader mapping and config parse: `IO/reader/Geant4SimReader.cpp`
- FAIR algorithm mapping and config parse: `simulation/geant4/Geant4Alg.cpp`

## Next actions to close parity

1. Introduce true cell-unit placement hierarchy (source-style repeated units).
2. Port optical surfaces + material property tables from source definitions.
3. Replace hard-coded NIST mapping with source-equivalent material table/indices.
4. Reproduce trigger physical naming/hierarchy exactly.
5. Add parity validation: per-layer Edep/time and secondary distributions under
   same input/seed conditions.
