# Part C (75 min): Implement Cluster `AnaAlg`

## Objective
Analyze `std::vector<TutorialCluster>` and produce basic correlation plots.

## Input
- Cluster container produced by Part B.
- `TutorialCluster` type from `docs/tutorial/answer/files/TutorialCluster.hpp`
- `docs/tutorial/answer/files/SimpleClusteringRecoAlg.hpp`.

## 1) Required observables
- `Ncluster = clusters.size()`
- `totalEdep_event = sum(cluster.totalEdep)`
- `Ncluster vs totalEdep_event` correlation
- `NclusterAboveThreshold` where `cluster.totalEdep >= edep_threshold`
- `threshold vs Ncluster` scan plot (x: threshold, y: Ncluster above threshold)

## 2) Optional observables (recommended)
- `meanClusterSize = mean(cluster.nHitInCluster)`
- `meanClusterEdep = mean(cluster.avgEdep)`

## 3) Run Part C config (analysis inside same reco pipeline)
```bash
./exe/ExecCalib docs/tutorial/config/partC_ana.yaml
```

<details>
  <summary>Answer files (full example)</summary>

- `docs/tutorial/answer/files/ClusterBasicAnaAlg.hpp`
- `docs/tutorial/answer/files/ClusterBasicAnaAlg.cpp`
</details>

## Deliverable
- New cluster AnaAlg source/header
- `docs/tutorial/config/partC_ana.yaml`
- Output file containing required histograms


## Config note (threshold scan)
In `ClusterBasicAnaAlg` cfg, use:
- `threshold_scan_max`
- `threshold_scan_step`

Then scan threshold from `0` to `threshold_scan_max` with that step.


> `partC_ana.yaml` uses `###run-number###` placeholders; replace them before execution.

> This Part C config does **not** use `RootInput`; analysis is executed in the same chain right after reco/clustering.