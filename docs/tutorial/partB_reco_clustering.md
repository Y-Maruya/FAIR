# Part B (90 min): Implement Clustering `RecoAlg`

## Objective
Create a clustering reconstruction algorithm and output `std::vector<TutorialCluster>`.

## Important design rule
For this tutorial package, define cluster struct in the tutorial answer files (`docs/tutorial/answer/files/`).
(Do not change production EDM headers for this tutorial step.)

---

## 1) Clustering concept (distance-based)
Use hit position (`Xpos`, `Ypos`, `Zpos`) and group hits when they are spatially close.

For each candidate hit `h_new` and current cluster seed/last hit `h_ref`, compute distance:

```text
d = sqrt((x_new-x_ref)^2 + (y_new-y_ref)^2 + (z_new-z_ref)^2)
```

If `d <= cluster_distance_mm`, include into same cluster; otherwise start new cluster.

### Configurable parameter (required)
`cluster_distance_mm` must come from YAML config, not hardcoded.

<details>
  <summary>Hint</summary>

Start with a simple value like `cluster_distance_mm: 45.0` and tune later.
</details>

---

## 2) Define tutorial cluster type in answer reco files
`Cluster` should contain:
- `std::vector<AHCALRecoHit> hits`
- `std::vector<int> index`
- summary vars: `nHitInCluster`, `totalEdep`, `avgEdep`

And have `updateSummary()`.

---

## 3) Implement RecoAlg
1. Copy a small existing RecoAlg.
2. Read input hits.
3. Cluster hits using `cluster_distance_mm`.
4. Fill `hits/index`.
5. Call `updateSummary()` for each cluster.
6. Output `std::vector<TutorialCluster>`.

<details>
  <summary>Answer files (full example)</summary>

- `docs/tutorial/answer/files/TutorialCluster.hpp`
- `docs/tutorial/answer/files/SimpleClusteringRecoAlg.hpp`
- `docs/tutorial/answer/files/SimpleClusteringRecoAlg.cpp`
- `docs/tutorial/answer/files/partB_reco_answer.yaml`

(These files are copy-ready answer examples.)
</details>

---

## 4) Run Part B config (`partB_reco.yaml`)
```bash
./bin/fair_single docs/tutorial/config/partB_reco.yaml
```

`partB_reco.yaml` should include `cluster_distance_mm`.

## Deliverable
- New clustering RecoAlg source/header (tutorial answer files)
- `docs/tutorial/config/partB_reco.yaml`
- Run command and output path


## CMake update (required)
Add your new RecoAlg `.cpp` into the corresponding module `CMakeLists.txt`, then rebuild.

## Execution patterns to try
- Single-chain execution (reco and later analysis in one alg sequence)
- Run-by-run mode (separate config/output per run)


> `partB_reco.yaml` uses `###run-number###` placeholders; replace them before execution.

## fair_single options (actual)
- Basic: `./bin/fair_single <config_yaml>`
- Multi-input text list: `./bin/fair_single <config_yaml> -i <input_list.txt>`
