# Part C Answer (Detailed)

## Answer files (use as-is)
- `docs/tutorial/answer/files/ClusterBasicAnaAlg.hpp`
- `docs/tutorial/answer/files/ClusterBasicAnaAlg.cpp`

## Explanation
- Prefer cluster summary variables (`totalEdep`, `nHitInCluster`, `avgEdep`) for speed/readability.
- Fill histograms including thresholded cluster multiplicity:
  - `hNcluster`
  - `hTotalEdep`
  - `hNclusterVsTotalEdep`
  - `hNclusterAboveThreshold` (`cluster.totalEdep >= edep_threshold`)
  - `hThresholdVsNcluster` (x: threshold, y: Ncluster above threshold)

## Validation points
- Hist entries > 0.
- Correlation plot not empty.
- Plot shape changes when `cluster_distance_mm` changes in Part B config.


## Config for threshold scan
Use:
- `threshold_scan_max`
- `threshold_scan_step`


## Output files
Set in cfg:
- `out_root_filename`
- `write_to_png`
- `out_png_dir`
