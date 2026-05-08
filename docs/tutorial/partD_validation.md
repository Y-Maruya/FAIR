# Part D (Optional, 15-30 min): Validation and Reproducibility

## Objective
Ensure others can rerun your pipeline reliably.

## Validation checklist
- [ ] Build passes.
- [ ] Baseline chain runs.
- [ ] New clustering `RecoAlg` runs.
- [ ] `std::vector<TutorialCluster>` is produced.
- [ ] `AnaAlg` outputs include:
  - [ ] `Ncluster`
  - [ ] `totalEdep`
  - [ ] `Ncluster` vs `totalEdep`

<details>
  <summary>Hint: fast debug order</summary>

1. Confirm Part B output file has cluster branch.
2. Confirm Part C input points to Part B output.
3. Confirm `clusters_key` matches producer/consumer names.
</details>

<details>
  <summary>Answer (minimum reproducibility block)</summary>

```text
commit: <hash>
build: <exact command>
run partB: <exact command>
run partC: <exact command>
outputs: output_partB.root, output_partC.root
```
</details>

## Reproducibility note (required)
Create a short text file with:
- commit hash
- exact build command
- exact run command
- output file path(s)

## Definition of done
All checklist items are complete and another person can rerun your commands.
