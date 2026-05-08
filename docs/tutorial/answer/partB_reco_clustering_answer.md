# Part B Answer (Detailed)

This section provides copy-ready answer files + explanation.

## Answer files (use as-is)
- `docs/tutorial/answer/files/SimpleClusteringRecoAlg.hpp`
- `docs/tutorial/answer/files/SimpleClusteringRecoAlg.cpp`
- `docs/tutorial/answer/files/partB_reco_answer.yaml`

## Explanation
- `cluster_distance_mm` is read from config, so behavior is tunable without recompilation.
- Clustering uses 3D distance from the last hit in current cluster.
- `SimpleClusteringRecoAlg` follows actual framework style (`initialize/execute/finalize/parse_cfg`).
- Each finalized cluster calls `updateSummary()` to fill:
  - `nHitInCluster`
  - `totalEdep`
  - `avgEdep`

## Validation points
- Cluster container exists.
- Cluster summaries are non-zero for non-empty events.
- Changing `cluster_distance_mm` changes cluster multiplicity trend.
