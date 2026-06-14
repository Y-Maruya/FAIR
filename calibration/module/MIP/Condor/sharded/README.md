# Sharded MIP Histogram Fitting

Channels are distributed using:

```text
cellid % shard_count == shard_index
```

Submit five fit jobs and automatically merge their result trees:

```bash
cd calibration/module/MIP/Condor/sharded
./submit_sharded.sh \
  /path/to/mip_neighborcheck_nofit.root \
  /path/to/output_directory \
  mip_neighborcheck_fitted \
  /path/to/simulation_mip_neighborcheck_nofit.root \
  5
```

The final output is `mip_neighborcheck_fitted.root`. It is created with `hadd`
after all five shard jobs finish successfully.

Generate the executable and DAG without submitting:

```bash
DRY_RUN=1 ./submit_sharded.sh input.root output_dir output_name simulation.root 5
```
