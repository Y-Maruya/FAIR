from influxdb import InfluxDBClient
import argparse
import json
import math
import os
import re
from pathlib import Path
from datetime import datetime, timezone

import requests
import matplotlib.pyplot as plt

from Luminosity import (
    UB_TO_FB_INV,
    align_by_timestamp,
    build_ratio_lumi_denominator,
    correction_window_values,
    discover_lumi_measurement,
    find_run_time,
    influx_interval_str,
    make_client,
    mask_unphysical_for_plot,
    median_positive,
    moving_average_nan,
    query_counter_series,
    query_mean_series,
    replace_environment_variables,
    solve_true_rate_from_reco,
    to_iso_utc,
    to_rate_from_counter,
)


TRIGGER_MEASUREMENTS = {
    "AnyVetoed": "tlureceiver00-TriggerVetoed",
    "NonVetoed": "tlureceiver00-Trigger_i5_0_i6_0",
    "NonPhysical": "ahcaleventreceiver00-DuplicatedEventCount",
    "Physical": "ahcaleventreceiver00-GoodEventsCount",
    "PreVeto": "tlureceiver00-PreVetoTriggers",
    "PostVeto": "tlureceiver00-PostVetoTriggers",
    "Recorded": "ahcaleventreceiver00-EventNumber",
}

BEAM_INTENSITY_SCALE = 1.0e13
TITLE_FONT_SIZE = 20
LABEL_FONT_SIZE = 16
TICK_FONT_SIZE = 15
LEGEND_FONT_SIZE = 12
ANNOTATION_FONT_SIZE = 11
LUMINOSITY_DISPLAY_SCALE = 0.7


def scale_series_to_reference(values, reference_values):
    finite_values = [value for value in values if isinstance(value, (int, float)) and math.isfinite(value)]
    finite_reference = [value for value in reference_values if isinstance(value, (int, float)) and math.isfinite(value)]
    if not finite_values or not finite_reference:
        return [math.nan] * len(values), math.nan

    value_scale = median_positive(finite_values)
    reference_scale = median_positive(finite_reference)
    if not math.isfinite(value_scale) or not math.isfinite(reference_scale) or value_scale <= 0.0 or reference_scale <= 0.0:
        return [math.nan] * len(values), math.nan

    factor = (reference_scale / value_scale) * LUMINOSITY_DISPLAY_SCALE
    scaled = [value * factor if isinstance(value, (int, float)) and math.isfinite(value) else math.nan for value in values]
    return scaled, factor


def query_beam1_intensity_series(client, start_s, end_s, bin_seconds, verbose=False):
    interval = influx_interval_str(bin_seconds)
    query = (
        'SELECT mean("TotalBeamIntensity") / 10000000000000 AS v FROM "dip" '
        f"WHERE time >= '{to_iso_utc(start_s)}' and time <= '{to_iso_utc(end_s)}' and \"Beam\" = 'Beam1' "
        f"GROUP BY time({interval}) fill(null)"
    )
    if verbose:
        print(f"[beam] {query}")
    result = client.query(query)
    points = list(result.get_points())
    return [(datetime.fromisoformat(p["time"].replace("Z", "+00:00")).timestamp(), p.get("v")) for p in points]


def save_csv(path, ts_sorted, aligned, corrected_rate, corrected_rate_bisection, ratio_series, ratio_series_bisection, start_s):
    with open(path, "w", encoding="utf-8") as f:
        headers = [
            "timestamp_utc",
            "elapsed_s",
            "AnyVetoed_rate_Hz",
            "NonVetoed_rate_Hz",
            "NonPhysical_rate_Hz",
            "Physical_rate_Hz",
            "PreVeto_rate_Hz",
            "PostVeto_rate_Hz",
            "Recorded_rate_Hz",
            "beam1_intensity_scaled",
            "beam1_intensity_smoothed_scaled",
            "beam1_intensity_for_ratio_scaled",
            "luminosity_normalized",
            "luminosity_smoothed_normalized",
            "luminosity_scaled_for_display",
            "corrected_non_vetoed_rate_Hz",
            "corrected_non_vetoed_rate_bisection_Hz",
            "corrected_over_beam1_intensity",
            "corrected_over_beam1_intensity_bisection",
        ]
        f.write(",".join(headers) + "\n")
        for i, ts_s in enumerate(ts_sorted):
            stamp = datetime.fromtimestamp(ts_s, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            elapsed = ts_s - start_s

            def fmt(v):
                return "" if (v is None or (isinstance(v, float) and math.isnan(v))) else f"{v:.10g}"

            row = [
                stamp,
                f"{elapsed:.3f}",
                fmt(aligned["AnyVetoed_rate"][ts_s]),
                fmt(aligned["NonVetoed_rate"][ts_s]),
                fmt(aligned["NonPhysical_rate"][ts_s]),
                fmt(aligned["Physical_rate"][ts_s]),
                fmt(aligned["PreVeto_rate"][ts_s]),
                fmt(aligned["PostVeto_rate"][ts_s]),
                fmt(aligned["Recorded_rate"][ts_s]),
                fmt(aligned["Beam1Intensity"][ts_s]),
                fmt(aligned["Beam1Intensity_smoothed"][ts_s]),
                fmt(aligned["Beam1Intensity_for_ratio"][ts_s]),
                fmt(aligned["Luminosity"][ts_s]),
                fmt(aligned["Luminosity_smoothed"][ts_s]),
                fmt(aligned["Luminosity_scaled"][ts_s]),
                fmt(corrected_rate[i]),
                fmt(corrected_rate_bisection[i]),
                fmt(ratio_series[i]),
                fmt(ratio_series_bisection[i]),
            ]
            f.write(",".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Fetch trigger/beam-intensity time series for one run and make comparison plots.")
    parser.add_argument("--run", required=True, type=int, help="AHCAL run number")
    parser.add_argument("--bin-seconds", type=float, default=0.2, help="Bin width in seconds for beam intensity InfluxDB time binning")
    parser.add_argument("--trigger-bin-seconds", type=float, default=1.0, help="Bin width in seconds for trigger counter time binning")
    parser.add_argument("--correction-window-seconds", type=float, default=30.0, help="Window size n [s] used in corrected-rate calculation")
    parser.add_argument("--correction-smooth-method", choices=["gaussian", "moving-average"], default="gaussian", help="Smoothing method over n-second window before correction")
    parser.add_argument("--beam-plot-smooth-bins", type=int, default=11, help="Moving-average window (bins) for displayed beam intensity")
    parser.add_argument("--beam-ratio-smooth-bins", type=int, default=5, help="Moving-average window (bins) for beam intensity denominator in ratio")
    parser.add_argument("--beam-ratio-min-frac", type=float, default=0.05, help="Minimum denominator floor as fraction of median smoothed beam intensity")
    parser.add_argument("--beam-min-abs", type=float, default=1.0e-6, help="Absolute near-zero beam intensity threshold; ratio is masked below this value")
    parser.add_argument("--assumed-deadtime-ms", type=float, default=1.3, help="Assumed deadtime D in ms for bisection correction")
    parser.add_argument("--assumed-auto-triggered-ms", type=float, default=3.9, help="Assumed auto-triggered window L in ms for bisection correction")
    parser.add_argument("--plot-start-seconds", type=float, default=0.0, help="Plot start time in seconds from run start")
    parser.add_argument("--plot-end-seconds", type=float, default=None, help="Plot end time in seconds from run start; defaults to run end")
    parser.add_argument("--plot-max-rate-hz", type=float, default=1.0e4, help="Upper plot cap for corrected rate")
    parser.add_argument("--plot-max-ratio", type=float, default=1.0e12, help="Upper plot cap for ratio")
    parser.add_argument("--lumi-measurement", default=None, help="Measurement name in DcsData for luminosity")
    parser.add_argument("--lumi-field", default="Luminosity", help="Field name in DcsData luminosity measurement")
    parser.add_argument("--output-dir", default=".", help="Output directory for plots and CSV/JSON")
    parser.add_argument("--verbose", "-v", action="store_true", help="Debug output")
    args = parser.parse_args()

    if args.bin_seconds <= 0:
        raise RuntimeError("--bin-seconds must be > 0")
    if args.trigger_bin_seconds <= 0:
        raise RuntimeError("--trigger-bin-seconds must be > 0")
    if args.correction_window_seconds <= 0:
        raise RuntimeError("--correction-window-seconds must be > 0")
    if args.beam_plot_smooth_bins <= 0:
        raise RuntimeError("--beam-plot-smooth-bins must be > 0")
    if args.beam_ratio_smooth_bins <= 0:
        raise RuntimeError("--beam-ratio-smooth-bins must be > 0")
    if args.beam_ratio_min_frac < 0:
        raise RuntimeError("--beam-ratio-min-frac must be >= 0")
    if args.beam_min_abs < 0:
        raise RuntimeError("--beam-min-abs must be >= 0")
    if args.assumed_deadtime_ms < 0:
        raise RuntimeError("--assumed-deadtime-ms must be >= 0")
    if args.assumed_auto_triggered_ms <= 0:
        raise RuntimeError("--assumed-auto-triggered-ms must be > 0")
    if args.plot_start_seconds < 0:
        raise RuntimeError("--plot-start-seconds must be >= 0")
    if args.plot_max_rate_hz <= 0:
        raise RuntimeError("--plot-max-rate-hz must be > 0")
    if args.plot_max_ratio <= 0:
        raise RuntimeError("--plot-max-ratio must be > 0")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    start_s, end_s = find_run_time(args.run)
    if end_s <= start_s:
        raise RuntimeError(f"Invalid run window: start={start_s}, end={end_s}")

    trigger_db_name = replace_environment_variables("$INFLUXDB")
    if trigger_db_name.startswith("$"):
        raise RuntimeError("Could not resolve $INFLUXDB for trigger database")

    trigger_client = make_client(trigger_db_name)
    dcs_client = make_client("DcsData")

    print(f"Run {args.run}: start={to_iso_utc(start_s)} end={to_iso_utc(end_s)}")
    print(
        f"Beam1 intensity bin={args.bin_seconds}s, trigger bin={args.trigger_bin_seconds}s, "
        f"correction window n={args.correction_window_seconds}s, method={args.correction_smooth_method}"
    )
    print(f"Trigger DB={trigger_db_name}, Beam DB=DcsData, Luminosity DB=DcsData")

    counter_rate_series = {}
    missing_measurements = []
    for label, measurement in TRIGGER_MEASUREMENTS.items():
        counter_series = query_counter_series(
            trigger_client,
            measurement,
            start_s,
            end_s,
            args.trigger_bin_seconds,
            field="value",
            verbose=args.verbose,
        )
        if not counter_series:
            missing_measurements.append(measurement)
            counter_rate_series[f"{label}_rate"] = []
            continue
        counter_rate_series[f"{label}_rate"] = to_rate_from_counter(counter_series, args.trigger_bin_seconds)

    if missing_measurements:
        print(f"Warning: missing trigger measurements in run window: {missing_measurements}")

    beam_series = query_beam1_intensity_series(
        dcs_client,
        start_s,
        end_s,
        args.bin_seconds,
        verbose=args.verbose,
    )
    if not beam_series:
        raise RuntimeError("No Beam1 intensity points found in run window")

    lumi_measurement = args.lumi_measurement or discover_lumi_measurement(dcs_client, verbose=args.verbose)
    if not lumi_measurement:
        raise RuntimeError("Could not discover luminosity measurement")

    lumi_field_candidates = [args.lumi_field]
    for fallback_field in ["Luminosity", "value"]:
        if fallback_field not in lumi_field_candidates:
            lumi_field_candidates.append(fallback_field)

    lumi_series = []
    used_lumi_field = None
    for lumi_field in lumi_field_candidates:
        lumi_series = query_mean_series(
            dcs_client,
            lumi_measurement,
            start_s,
            end_s,
            args.bin_seconds,
            field=lumi_field,
            verbose=args.verbose,
        )
        if lumi_series:
            used_lumi_field = lumi_field
            break

    if not lumi_series:
        raise RuntimeError(
            "No luminosity points found in run window "
            f"(measurement={lumi_measurement}, tried fields={lumi_field_candidates})"
        )
    print(f"Using luminosity measurement={lumi_measurement}, field={used_lumi_field}")

    series_dict = dict(counter_rate_series)
    series_dict["Beam1Intensity"] = beam_series
    series_dict["Luminosity"] = lumi_series
    ts_sorted, aligned = align_by_timestamp(series_dict)

    for key in series_dict:
        data = [math.nan if aligned[key].get(ts) is None else aligned[key].get(ts) for ts in ts_sorted]
        for i, ts_s in enumerate(ts_sorted):
            aligned[key][ts_s] = data[i]

    beam_values = [aligned["Beam1Intensity"][ts_s] for ts_s in ts_sorted]
    beam_plot_smoothed = moving_average_nan(beam_values, args.beam_plot_smooth_bins)
    for i, ts_s in enumerate(ts_sorted):
        aligned.setdefault("Beam1Intensity_smoothed", {})[ts_s] = beam_plot_smoothed[i]

    beam_ratio_denom, beam_ratio_floor = build_ratio_lumi_denominator(
        beam_values,
        args.beam_ratio_smooth_bins,
        args.beam_ratio_min_frac,
    )
    for i, ts_s in enumerate(ts_sorted):
        aligned.setdefault("Beam1Intensity_for_ratio", {})[ts_s] = beam_ratio_denom[i]

    lumi_values = [aligned["Luminosity"][ts_s] for ts_s in ts_sorted]
    lumi_plot_smoothed = moving_average_nan(lumi_values, args.beam_plot_smooth_bins)
    for i, ts_s in enumerate(ts_sorted):
        aligned.setdefault("Luminosity_smoothed", {})[ts_s] = lumi_plot_smoothed[i]

    lumi_scaled_for_display, lumi_display_scale = scale_series_to_reference(
        lumi_plot_smoothed,
        beam_plot_smoothed,
    )
    for i, ts_s in enumerate(ts_sorted):
        aligned.setdefault("Luminosity_scaled", {})[ts_s] = lumi_scaled_for_display[i]

    correction_window_bins = max(1, int(round(args.correction_window_seconds / args.trigger_bin_seconds)))
    if correction_window_bins > 1 and correction_window_bins % 2 == 0:
        correction_window_bins += 1

    any_for_corr = correction_window_values(
        [aligned["AnyVetoed_rate"][ts_s] for ts_s in ts_sorted],
        correction_window_bins,
        method=args.correction_smooth_method,
    )
    non_for_corr = correction_window_values(
        [aligned["NonVetoed_rate"][ts_s] for ts_s in ts_sorted],
        correction_window_bins,
        method=args.correction_smooth_method,
    )
    np_for_corr = correction_window_values(
        [aligned["NonPhysical_rate"][ts_s] for ts_s in ts_sorted],
        correction_window_bins,
        method=args.correction_smooth_method,
    )
    phys_for_corr = correction_window_values(
        [aligned["Physical_rate"][ts_s] for ts_s in ts_sorted],
        correction_window_bins,
        method=args.correction_smooth_method,
    )

    corrected_non_veto_rate = []
    corrected_non_veto_rate_bisection = []
    ratio_series = []
    ratio_series_bisection = []
    linear_upper_limit_bins = []
    bisection_upper_limit_bins = []
    invalid_ratio_bins = 0
    near_zero_beam_bins = 0

    for i, ts_s in enumerate(ts_sorted):
        any_vetoed = any_for_corr[i]
        non_vetoed = non_for_corr[i]
        non_physical = np_for_corr[i]
        physical = phys_for_corr[i]
        beam_raw = aligned["Beam1Intensity"][ts_s]
        beam_for_ratio = aligned["Beam1Intensity_for_ratio"][ts_s]

        linear_hit_upper_limit = False
        if any(math.isnan(x) for x in [non_vetoed, non_physical, physical]):
            corrected_linear = math.nan
        else:
            denom = 1.0 - (args.assumed_deadtime_ms * non_physical + 1.4 * physical) / 1000.0
            if denom <= 0.0:
                corrected_linear = math.nan
            else:
                corrected_linear = non_vetoed / denom
                if corrected_linear > args.plot_max_rate_hz:
                    linear_hit_upper_limit = True
                    corrected_linear = math.nan
        corrected_non_veto_rate.append(corrected_linear)
        linear_upper_limit_bins.append(linear_hit_upper_limit)

        bisection_hit_upper_limit = False
        if any(math.isnan(x) for x in [any_vetoed, non_vetoed]):
            corrected_bisection = math.nan
        else:
            reco_total = any_vetoed + non_vetoed
            true_total = solve_true_rate_from_reco(
                reco_total,
                args.assumed_deadtime_ms / 1000.0,
                args.assumed_auto_triggered_ms / 1000.0,
            )
            if true_total > args.plot_max_rate_hz:
                bisection_hit_upper_limit = True
                true_total = args.plot_max_rate_hz
            corrected_bisection = math.nan if reco_total <= 0.0 else non_vetoed * (true_total / reco_total)
        corrected_non_veto_rate_bisection.append(corrected_bisection)
        bisection_upper_limit_bins.append(bisection_hit_upper_limit)

        beam_near_zero = isinstance(beam_raw, (int, float)) and not math.isnan(beam_raw) and beam_raw < args.beam_min_abs
        if beam_near_zero:
            ratio_series.append(math.nan)
            ratio_series_bisection.append(math.nan)
            invalid_ratio_bins += 1
            near_zero_beam_bins += 1
        else:
            if math.isnan(corrected_linear) or math.isnan(beam_for_ratio) or beam_for_ratio == 0:
                ratio_series.append(math.nan)
            else:
                ratio_series.append((corrected_linear / beam_for_ratio) * UB_TO_FB_INV)

            if math.isnan(corrected_bisection) or math.isnan(beam_for_ratio) or beam_for_ratio == 0:
                ratio_series_bisection.append(math.nan)
            else:
                ratio_series_bisection.append((corrected_bisection / beam_for_ratio) * UB_TO_FB_INV)

    elapsed = [ts - start_s for ts in ts_sorted]
    plot_end_seconds = args.plot_end_seconds if args.plot_end_seconds is not None else (elapsed[-1] if elapsed else 0.0)
    if plot_end_seconds <= args.plot_start_seconds:
        raise RuntimeError("--plot-end-seconds must be greater than --plot-start-seconds")

    plot_indices = [i for i, e in enumerate(elapsed) if args.plot_start_seconds <= e <= plot_end_seconds]
    if not plot_indices:
        raise RuntimeError("No points fall inside the requested plot time range")

    plot_elapsed = [elapsed[i] for i in plot_indices]
    plot_linear = [corrected_non_veto_rate[i] for i in plot_indices]
    plot_bisection = [corrected_non_veto_rate_bisection[i] for i in plot_indices]
    plot_beam_raw = [aligned["Beam1Intensity"][ts_sorted[i]] for i in plot_indices]
    plot_beam_smooth = [aligned["Beam1Intensity_smoothed"][ts_sorted[i]] for i in plot_indices]
    plot_lumi_raw = [aligned["Luminosity"][ts_sorted[i]] for i in plot_indices]
    plot_lumi_smooth = [aligned["Luminosity_smoothed"][ts_sorted[i]] for i in plot_indices]
    plot_ratio = [ratio_series[i] for i in plot_indices]
    plot_ratio_bisection = [ratio_series_bisection[i] for i in plot_indices]
    plot_linear_hit_upper_limit = [linear_upper_limit_bins[i] for i in plot_indices]
    plot_bisection_hit_upper_limit = [bisection_upper_limit_bins[i] for i in plot_indices]

    plot_linear = mask_unphysical_for_plot(plot_linear, args.plot_max_rate_hz)
    plot_bisection = mask_unphysical_for_plot(plot_bisection, args.plot_max_rate_hz)
    plot_ratio = mask_unphysical_for_plot(plot_ratio, args.plot_max_ratio)
    plot_ratio_bisection = mask_unphysical_for_plot(plot_ratio_bisection, args.plot_max_ratio)
    plot_lumi_smooth_scaled, lumi_plot_scale = scale_series_to_reference(plot_lumi_smooth, plot_beam_smooth)

    has_linear = any((not math.isnan(v)) for v in plot_linear)
    has_bisection = any((not math.isnan(v)) for v in plot_bisection)
    has_ratio = any((not math.isnan(v)) for v in plot_ratio)
    has_ratio_bisection = any((not math.isnan(v)) for v in plot_ratio_bisection)
    has_linear_cap = any(plot_linear_hit_upper_limit)
    has_bisection_cap = any(plot_bisection_hit_upper_limit)

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(12, 10),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 2]},
    )

    ax_top_right = ax_top.twinx()
    if has_linear:
        ax_top.plot(plot_elapsed, plot_linear, color="tab:blue", linewidth=2.2, label="Corrected NonVetoed rate")
    if has_linear_cap:
        cap_times = [t for t, hit in zip(plot_elapsed, plot_linear_hit_upper_limit) if hit]
        ax_top.plot(
            cap_times,
            [0.97] * len(cap_times),
            color="tab:orange",
            marker="o",
            linestyle="None",
            label="linear correction hit upper limit",
            zorder=5,
        )
    if has_bisection:
        ax_top.plot(plot_elapsed, plot_bisection, color="tab:purple", linewidth=2.0, linestyle="--", label="Corrected NonVetoed rate (bisection)")
    if has_bisection_cap:
        cap_times = [t for t, hit in zip(plot_elapsed, plot_bisection_hit_upper_limit) if hit]
        ax_top.plot(
            cap_times,
            [0.97] * len(cap_times),
            color="tab:orange",
            marker="v",
            linestyle="None",
            label="bisection hit upper limit",
            zorder=5,
        )
        ax_top.text(
            0.99,
            0.94,
            f"orange markers: corrected rate clipped at >= {args.plot_max_rate_hz:.0f} Hz",
            transform=ax_top.transAxes,
            ha="right",
            va="top",
            color="tab:orange",
            fontsize=ANNOTATION_FONT_SIZE,
        )
    if plot_beam_smooth and any((not math.isnan(v)) for v in plot_beam_smooth):
        ax_top_right.plot(plot_elapsed, plot_beam_raw, color="tab:red", linewidth=1.0, alpha=0.22, label="Beam1 intensity (raw)")
        ax_top_right.plot(plot_elapsed, plot_beam_smooth, color="tab:red", linewidth=2.1, alpha=0.95, label="Beam1 intensity (smoothed)")
    if plot_lumi_smooth and any((not math.isnan(v)) for v in plot_lumi_smooth):
        ax_top_right.plot(plot_elapsed, plot_lumi_raw, color="tab:green", linewidth=1.0, alpha=0.22, label="Luminosity (raw)")
    if plot_lumi_smooth_scaled and any((not math.isnan(v)) for v in plot_lumi_smooth_scaled):
        ax_top_right.plot(plot_elapsed, plot_lumi_smooth_scaled, color="tab:green", linewidth=2.1, alpha=0.95, label=f"Luminosity (scaled x{lumi_plot_scale:.3g})")
    ax_top.set_ylabel("Corrected NonVetoed Rate [Hz]", color="tab:blue", fontsize=LABEL_FONT_SIZE, labelpad=10)
    ax_top_right.set_ylabel("Beam1 intensity / luminosity", color="tab:red", fontsize=LABEL_FONT_SIZE, labelpad=16)
    ax_top.set_ylim(bottom=0.0)
    ax_top_right.set_ylim(bottom=0.0)
    ax_top.grid(True, alpha=0.25)
    ax_top.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax_top_right.tick_params(axis="y", labelsize=TICK_FONT_SIZE, pad=4)
    lines_left, labels_left = ax_top.get_legend_handles_labels()
    lines_right, labels_right = ax_top_right.get_legend_handles_labels()
    ax_top.legend(lines_left + lines_right, labels_left + labels_right, loc="upper right", fontsize=LEGEND_FONT_SIZE)
    ax_top.set_title(f"Run {args.run}: corrected NonVetoed rate, Beam1 intensity, and luminosity", fontsize=TITLE_FONT_SIZE)

    if has_ratio:
        ax_bottom.plot(plot_elapsed, plot_ratio, color="tab:green", linewidth=2.0, label="Corrected NonVetoed rate / Beam1 intensity")
    if has_linear_cap:
        cap_times = [t for t, hit in zip(plot_elapsed, plot_linear_hit_upper_limit) if hit]
        ax_bottom.plot(
            cap_times,
            [0.96] * len(cap_times),
            color="tab:orange",
            marker="o",
            linestyle="None",
            transform=ax_bottom.get_xaxis_transform(),
            label="linear correction hit upper limit",
            zorder=5,
        )
    if has_ratio_bisection:
        ax_bottom.plot(plot_elapsed, plot_ratio_bisection, color="tab:purple", linewidth=1.8, linestyle="--", label="Corrected NonVetoed rate / Beam1 intensity (bisection)")
    if has_bisection_cap:
        cap_times = [t for t, hit in zip(plot_elapsed, plot_bisection_hit_upper_limit) if hit]
        ax_bottom.plot(
            cap_times,
            [0.96] * len(cap_times),
            color="tab:orange",
            marker="v",
            linestyle="None",
            transform=ax_bottom.get_xaxis_transform(),
            label="bisection hit upper limit",
            zorder=5,
        )
        ax_bottom.text(
            0.99,
            0.94,
            "orange markers: corrected rate clipped at upper limit",
            transform=ax_bottom.transAxes,
            ha="right",
            va="top",
            color="tab:orange",
            fontsize=ANNOTATION_FONT_SIZE,
        )
    ax_bottom.set_xlabel("Elapsed time from run start [s]", fontsize=LABEL_FONT_SIZE)
    ax_bottom.set_ylabel("Rate / Beam1 intensity", fontsize=LABEL_FONT_SIZE)
    ax_bottom.set_ylim(bottom=0.0)
    ax_bottom.grid(True, alpha=0.25)
    ax_bottom.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax_bottom.legend(loc="upper right", fontsize=LEGEND_FONT_SIZE)

    fig.tight_layout()
    plot_path = output_dir / f"run_{args.run}_beamintensity_plots.png"
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)

    csv_path = output_dir / f"run_{args.run}_beamintensity_timeseries.csv"
    save_csv(
        csv_path,
        ts_sorted,
        aligned,
        corrected_non_veto_rate,
        corrected_non_veto_rate_bisection,
        ratio_series,
        ratio_series_bisection,
        start_s,
    )

    meta = {
        "RunNumber": args.run,
        "StartUTC": to_iso_utc(start_s),
        "EndUTC": to_iso_utc(end_s),
        "BinSeconds": args.bin_seconds,
        "TriggerBinSeconds": args.trigger_bin_seconds,
        "CorrectionWindowSeconds": args.correction_window_seconds,
        "CorrectionWindowBins": correction_window_bins,
        "CorrectionSmoothMethod": args.correction_smooth_method,
        "AssumedDeadtimeMs": args.assumed_deadtime_ms,
        "AssumedAutoTriggeredMs": args.assumed_auto_triggered_ms,
        "TriggerDB": trigger_db_name,
        "BeamDB": "DcsData",
        "LuminosityDB": "DcsData",
        "BeamMeasurement": "dip",
        "BeamField": "TotalBeamIntensity",
        "BeamTag": "Beam1",
        "BeamScaleApplied": "x1/1e13",
        "LuminosityMeasurement": lumi_measurement,
        "LuminosityField": used_lumi_field,
        "BeamPlotSmoothBins": args.beam_plot_smooth_bins,
        "BeamRatioSmoothBins": args.beam_ratio_smooth_bins,
        "BeamRatioMinFrac": args.beam_ratio_min_frac,
        "BeamMinAbs": args.beam_min_abs,
        "LuminosityDisplayScale": lumi_display_scale,
        "BeamRatioFloor": beam_ratio_floor,
        "RatioScaleApplied": "x1e9 from rate/intensity",
        "MissingTriggerMeasurements": missing_measurements,
        "LinearUpperLimitBins": sum(1 for hit in linear_upper_limit_bins if hit),
        "BisectionUpperLimitBins": sum(1 for hit in bisection_upper_limit_bins if hit),
        "PlotStartSeconds": args.plot_start_seconds,
        "PlotEndSeconds": plot_end_seconds,
        "PlottedPoints": len(plot_indices),
        "Plots": [str(plot_path)],
        "CSV": str(csv_path),
    }
    json_path = output_dir / f"run_{args.run}_beamintensity_summary.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    print(f"Saved: {plot_path}")
    print(f"Saved: {csv_path}")
    print(f"Saved: {json_path}")
    print(f"Beam ratio floor: {beam_ratio_floor}")
    if invalid_ratio_bins > 0:
        print(f"Warning: {invalid_ratio_bins} bins were dropped in ratio due to invalid/too-small beam intensity denominator.")
    if near_zero_beam_bins > 0:
        print(f"Warning: {near_zero_beam_bins} bins were dropped in ratio because beam intensity < {args.beam_min_abs}.")


if __name__ == "__main__":
    main()
