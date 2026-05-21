# HG:LG Calibration Quality Metrics

This document describes the HG:LG calibration quality metrics feature implemented in `InterCalibAlg`.

## Overview

The `InterCalibAlg` algorithm now computes detailed quality metrics for each channel's HG:LG calibration fit. These metrics are stored in a ROOT TTree for later analysis, allowing users to identify problematic channels and understand fit quality without re-running the full calibration.

## Quick Start

Enable quality metrics in your YAML configuration:

```yaml
algs:
- type: InterCalibAlg
  cfg:
    # ... existing InterCalib parameters ...
    
    # Enable quality metrics output
    compute_quality_metrics: true
    quality_output_ttree: true
    quality_out_filename: "hglg_calib_quality.root"
```

The algorithm will automatically:
1. Compute fit quality metrics for each channel
2. Calculate layer-wise slope statistics
3. Assign quality flags based on configurable thresholds
4. Save results to a ROOT TTree for querying

## Output

### ROOT TTree: "HGLGCalibQuality"

The output file contains a TTree with 25 branches per channel:

**Channel Identifiers:**
- `run`: Run number
- `layer`, `chip`, `channel`: Geometric channel position
- `cellid`: Unique channel ID (cellid)

**Pedestal Values:**
- `hg_pedestal`: High Gain ADC pedestal
- `lg_pedestal`: Low Gain ADC pedestal

**Fit Parameters:**
- `slope`: HG/LG ratio (main calibration parameter)
- `intercept`: Linear fit intercept
- `slope_error`: Error on slope (currently placeholder)
- `intercept_error`: Error on intercept (currently placeholder)

**Fit Quality:**
- `chi2`: Chi-squared from fit
- `ndf`: Number of degrees of freedom
- `chi2_ndf`: Reduced chi-squared (chi2/ndf)
- `fit_status`: ROOT fit status (0 = success)

**Point Statistics:**
- `n_points_total`: Total calibration points before outlier rejection
- `n_points_used`: Points used in final fit
- `n_points_outlier`: Points rejected as outliers
- `outlier_fraction`: n_points_outlier / n_points_total

**Correlation & Layer Statistics:**
- `correlation`: Pearson correlation between HG_sub and LG_sub
- `slope_median_layer`: Median slope for this layer
- `slope_mad_layer`: Median absolute deviation of slopes in layer
- `slope_z_layer`: Layer-wise z-score: `(slope - median) / (1.4826 * MAD)`

**Quality Flags:**
- `quality_flag`: Bitmask of quality issues (see below)
- `is_pedestal_masked`: Channel failed pedestal calibration
- `is_hglg_calib_bad`: Channel failed HG:LG calibration quality checks
- `is_physics_mask_candidate`: Channel should be physics-masked (conservative)
- `manual_bad`: User marked as bad (reserved for interactive mode)
- `manual_good`: User marked as good (reserved for interactive mode)

## Quality Flags

The `quality_flag` field is a bitmask with the following bits:

| Flag | Value | Condition |
|------|-------|-----------|
| LOW_STAT | 1 << 0 | n_points_used < minimum_points |
| FIT_FAILED | 1 << 1 | fit_status != 0 or slope not finite |
| BAD_CORRELATION | 1 << 2 | correlation < bad_correlation_threshold |
| MANY_OUTLIERS | 1 << 3 | outlier_fraction > outlier_fraction_threshold |
| BAD_SLOPE_LAYER_Z | 1 << 4 | abs(slope_z_layer) > slope_layer_z_threshold |
| BAD_CHI2 | 1 << 5 | chi2_ndf > chi2_ndf_threshold |
| PEDESTAL_MASKED | 1 << 6 | Channel is pedestal-masked |
| MANUAL_BAD | 1 << 7 | User marked as bad |
| MANUAL_GOOD | 1 << 8 | User marked as good |

A channel is considered `is_hglg_calib_bad = true` if ANY of these flags are set:
LOW_STAT, FIT_FAILED, BAD_CORRELATION, MANY_OUTLIERS, BAD_SLOPE_LAYER_Z, BAD_CHI2

## Configuration Parameters

All parameters are optional with sensible defaults:

```yaml
# Enable/disable quality metrics computation
compute_quality_metrics: true

# Minimum number of points required for a channel to have valid fit
quality_minimum_points: 10

# Thresholds for quality flags
quality_bad_correlation_threshold: 0.8      # Min correlation to pass
quality_outlier_fraction_threshold: 0.2     # Max fraction of outliers
quality_slope_layer_z_threshold: 5.0        # Max |z-score| for slope
quality_chi2_ndf_threshold: 10.0            # Max chi2/ndf

# TTree output options
quality_output_ttree: true                       # Enable TTree output
quality_out_filename: "hglg_calib_quality.root"  # Output file name
quality_save_summary_plots: false                # Generate summary histograms (future)

# Interactive mode (Phase 2 - not yet implemented)
interactive_fit: false
interactive_mode: "bad_only"   # all, bad_only, layer, single_channel
interactive_layer: -1
interactive_chip: -1
interactive_channel: -1
```

## Usage Examples

### Query TTree for Bad Channels

```python
import ROOT
f = ROOT.TFile("hglg_calib_quality.root")
tree = f.Get("HGLGCalibQuality")

# Find all channels with bad correlation
tree.Draw("channel", "correlation < 0.8")

# Find all channels with high chi2
tree.Draw("layer:chip:channel", "chi2_ndf > 10")

# Count flags by type
tree.Project("hist_flags", "quality_flag", "")
```

### Analysis with C++

```cpp
TFile f("hglg_calib_quality.root");
TTree* tree = (TTree*)f.Get("HGLGCalibQuality");

int layer, chip, channel;
double slope, correlation, slope_z_layer;
uint32_t quality_flag;

tree->SetBranchAddress("layer", &layer);
tree->SetBranchAddress("chip", &chip);
tree->SetBranchAddress("channel", &channel);
tree->SetBranchAddress("slope", &slope);
tree->SetBranchAddress("correlation", &correlation);
tree->SetBranchAddress("slope_z_layer", &slope_z_layer);
tree->SetBranchAddress("quality_flag", &quality_flag);

for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    
    if (quality_flag & (1 << 2)) {  // BAD_CORRELATION
        std::cout << "Layer " << layer << " Chip " << chip 
                  << " Ch " << channel << ": correlation = " 
                  << correlation << std::endl;
    }
}
```

## Masking Policy

**Important**: A channel with bad HG:LG calibration is NOT automatically marked for physics masking.

- **Calibration-bad** (`is_hglg_calib_bad=true`): Channel's fit quality failed checks. Downstream processing should use layer median or chip median slope instead of this channel's fit.
- **Physics-mask candidate** (`is_physics_mask_candidate=true`): Channel should be excluded from physics analysis. Currently only set if:
  - Channel is pedestal-masked, OR
  - User manually marked as bad (interactive mode)

This conservative approach preserves data for analysis while flagging problematic calibrations.

## Future Enhancements (Phase 2)

The specification includes an interactive mode for visual inspection of fits:
- Display 2D histograms with fit lines for each channel
- Keyboard controls: Enter (next), q (quit), s (save canvas)
- Manual channel marking as good/bad
- Channel filtering (all, suspicious, layer, single channel)
- Canvas export with channel metadata

## Implementation Details

### Correlation Computation

Correlation is computed as:
```
r = cov(HG_sub, LG_sub) / (σ_HG × σ_LG)
```

Where HG_sub and LG_sub are pedestal-subtracted values from the fit points.

### Layer Statistics

For each layer:
1. Collect all valid channel slopes (exclude: pedestal-masked, failed fits, low stats, non-physical)
2. Compute median of slopes
3. Compute median absolute deviation (MAD)
4. For each channel: `z_score = (slope - median) / (1.4826 × MAD)`

The factor 1.4826 makes MAD consistent with standard deviation for Gaussian distributions.

### Quality Flag Logic

Flags are assigned independently based on thresholds, then combined to determine:
- `is_hglg_calib_bad`: True if any critical flags (LOW_STAT, FIT_FAILED, BAD_CORRELATION, MANY_OUTLIERS, BAD_SLOPE_LAYER_Z, BAD_CHI2) are set
- `is_physics_mask_candidate`: True only if channel is pedestal-masked or manually marked bad

This allows downstream processors to handle calibration-bad channels (e.g., by using layer statistics) while being conservative about physics masking.

## See Also

- [InterCalibAlg specification](../SPECIFICATION.md)
- ROOT TTree documentation
- AHCAL calibration framework
