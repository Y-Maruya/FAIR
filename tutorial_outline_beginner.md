# Half-Day FAIR Tutorial (Beginners)

## Goal (what students can do by the end)
In ~4 hours, students will be able to:
1. Build FAIR and run an existing executable.
2. Implement a simple clustering `RecoAlg`.
3. Output `std::vector<TutorialCluster>` from the reco stage.
4. Implement a simple `AnaAlg` to study cluster-level correlations.

---

## Scope and Style
- Simple and straight-forward: minimum theory, maximum hands-on.
- Copy-modify-run approach using existing code as templates.
- Finish with one working custom pipeline in half a day.

---

## 0. Before the tutorial (instructor prep)
Prepare these in advance:
- Working build environment (compiler, CMake, ROOT, dependencies).
- Small input sample that runs in a few minutes.
- One known-good baseline command and expected output.
- Starter branch/tag for students.

---

## 1. Half-day timeline (example: 4 hours)

### Part A (45 min): Build and baseline run
- Quick repository map (where algs/executables/configs are).
- Build command.
- Run one existing chain.
- Confirm output/log.

**Checkpoint A:** baseline workflow runs for everyone.

### Part B (90 min): Reco exercise — clustering
Students implement a minimal clustering algorithm in a new `RecoAlg`.

#### Required output data type
`RecoAlg` must output:
- `std::vector<TutorialCluster>`

Where `Cluster` contains:
- `std::vector<RecoHit> hits`
- `std::vector<int> index`

#### Work items
- Copy a small existing `RecoAlg` as template.
- Create/define `Cluster` structure (or use existing one if already available).
- Implement simple clustering logic (rule can be minimal and explicit).
- Fill `hits` and corresponding `index` for each cluster.
- Register/build the new algorithm and connect it into run chain.

**Checkpoint B:** run finishes and `std::vector<TutorialCluster>` is produced.

### Part C (75 min): Analysis exercise — cluster correlations
Students implement a new `AnaAlg` using the cluster output.

#### Minimum analysis items
- Number of clusters per event.
- Total edep per event (or per cluster set).
- Correlation between cluster count and total edep.

#### Work items
- Copy a simple existing `AnaAlg` as template.
- Read `std::vector<TutorialCluster>` from event data.
- Compute required observables.
- Fill histograms (e.g., 1D: `Ncluster`, `totalEdep`; 2D: `Ncluster` vs `totalEdep`).

**Checkpoint C:** analysis outputs are produced and readable.

### Part D (Optional, 15-30 min): Final validation
- Build passes.
- Chain runs end-to-end.
- Cluster output is present.
- Analysis histograms are present.
- Commands are documented for rerun.

**Final checkpoint (if Part D is done):** each student has a reproducible run with new reco+analysis.

---

## 2. Student tasks (explicit)

### Task 1: Baseline
- Build.
- Run baseline command.
- Record command and output path.

### Task 2: Clustering `RecoAlg`
- Add/confirm `Cluster` definition:
  - `std::vector<RecoHit>`
  - `std::vector<int> index`
- Implement clustering and output `std::vector<TutorialCluster>`.
- Rebuild and run without crash.

### Task 3: Cluster `AnaAlg`
- Read cluster container.
- Compute:
  - cluster count
  - total edep
  - count-edep correlation
- Write histogram outputs.

### Task 4: Reproducibility note
- Record:
  - commit hash
  - run command
  - output file names

---

## 3. Definition of done
Complete if all are true:
1. Code builds.
2. Baseline run works.
3. New clustering `RecoAlg` runs.
4. `std::vector<TutorialCluster>` is produced.
5. `AnaAlg` outputs include cluster count, total edep, and their correlation.
6. Workflow is documented for rerun.

---

## 4. Minimal instructor materials
- One-page cheat sheet:
  - build command
  - baseline run command
  - where to add Reco/Ana alg files
  - common errors (missing registration, CMake omission, wrong config path)
- Tiny reference implementation for fallback.

---

## 5. Optional extension (if time remains)
- Try two clustering parameter settings and compare plots.
- Add one more cluster observable (e.g., mean cluster size).
