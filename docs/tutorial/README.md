# FAIR Half-Day Tutorial (Master's Students)

This tutorial is a hands-on, straight-forward course to finish in about 4 hours.

## Learning goal
By the end, you will:
1. Build FAIR and run an existing pipeline.
2. Implement a clustering `RecoAlg` that outputs `std::vector<TutorialCluster>`.
3. Define/use tutorial `TutorialCluster` in `docs/tutorial/answer/files/`.
4. Implement an `AnaAlg` that studies cluster count, total edep, and their correlation.

## Structure
- `partA_baseline.md`: Build and baseline run.
- `partB_reco_clustering.md`: Implement clustering `RecoAlg` (with tutorial `TutorialCluster`).
- `partC_analysis_cluster.md`: Implement cluster analysis `AnaAlg`.
- `partD_validation.md`: Optional wrap-up (validation/reproducibility checklist).

## Answers
Reference answers are under:
- `docs/tutorial/answer/`

## Optional: hide/show answer blocks in Markdown
Yes — you can use HTML details tags in Markdown:

```markdown
<details>
  <summary>Show answer</summary>

  Answer text here.
</details>
```

Many renderers (GitHub, GitLab, VS Code preview) support this.


## Tutorial configs (YAML)
- `docs/tutorial/config/partB_reco.yaml`: run clustering RecoAlg flow
- `docs/tutorial/config/partC_ana.yaml`: run cluster analysis flow


## Build integration note (CMake)
When you add new tutorial alg files, also update the relevant module `CMakeLists.txt`:
- add new `.cpp` to target sources
- ensure include paths are visible
- rebuild (`cmake .. && make -j`)

## Two execution patterns to document
1. **Single-chain reco+analysis (recommended in this tutorial)**
   - Use one config where analysis alg runs after reco algs in the same `algs` chain.
   - Example:
     ```bash
     ./bin/fair_single docs/tutorial/config/partC_ana.yaml
     ```

2. **Run-by-run execution**
   - Use configs with `###run-number###` placeholders and replace per run.
   - Example:
     ```bash
     RUN=21897
     sed "s/###run-number###/${RUN}/g" docs/tutorial/config/partB_reco.yaml > /tmp/partB_${RUN}.yaml
     sed "s/###run-number###/${RUN}/g" docs/tutorial/config/partC_ana.yaml > /tmp/partC_${RUN}.yaml
     ./bin/fair_run docs/tutorial/config/partB_reco.yaml -r ${RUN}
     ./bin/fair_run docs/tutorial/config/partC_ana.yaml -r ${RUN}
     ```
