# MIPSimAlg

Use `MIPSimAlg` after digitization and track fitting. It reads `SimHits`, the
digitized `SimRawHits`, and `FittedTrack`. Its ADC output is produced by the
existing `MIPAlg` itself, with the same track, crossing, hit, pedestal, and
Landau–Gaussian fit settings as the data configuration. The additional ROOT
file contains visible energy spectra in MeV for the **same selected ADC hits**.

Example algorithm configuration (place after the digitizer and `TrackFitAlg`):

```yaml
- type: MIPSimAlg
  cfg:
    in_simhit_key: SimHits
    in_rawhit_key: SimRawHits
    in_track_key: FittedTrack
    string_track_struct: SimpleFittedTrack
    track_selection_string: (valid>0)*(chi2xperndf<2)*(chi2xperndf>0)*(chi2yperndf<2)*(abs(direction_x)<0.001)*(abs(direction_y)<0.001)
    out_mip_filename: mip_sim_adc.root
    out_simhit_filename: mip_sim_truth.root
    fit: true
    calculate_fwhm: true
    substrate_pedestal: true
    read_pedestal_from_ROOT: true
    in_pedestal_file: pedestal.root
    read_pedestal_from_DB: false
    xy_size_threshold: 20
    enable_neighbor_layer_crossing_selection: true
    mip_neighbor_upstream_layers: 2
    mip_neighbor_downstream_layers: 2
    mip_edge_layer_policy: true
    mip_reject_if_neighbor_cell_hit: true
    mip_neighbor_cell_search_radius: 1
    mip_skip_layers: [0, 2, 14, 28]
    mip_skip_chips: [[21, 1], [24, 3], [24, 4], [24, 5], [24, 6], [24, 7]]
```

Point `in_pedestal_file` at the pedestal calibration used by the digitizer, or
choose the DB settings used by the data MIP run. All ordinary `MIPAlg` options
can be copied unchanged from the data run. `mip_sim_adc.root` contains the
standard `mip` tree: `MPV` is the Landau component parameter and `max_x` is
the peak of the fitted convolution when `calculate_fwhm` is true.

Only digitized hits with `hittag == 1` which actually enter an ADC histogram
receive a truth entry. If the digitizer inserts a noise hit without a SimHit,
`selected_without_simhit` counts it. Multiple SimHits for one cell are summed.
The optional `edep_nbin` and `edep_max` parameters configure the truth
histograms (defaults: 200 bins from 0 to 2 MeV).
