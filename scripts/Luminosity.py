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


TRIGGER_MEASUREMENTS = {
    "AnyVetoed": "tlureceiver00-TriggerVetoed",
    "NonVetoed": "tlureceiver00-Trigger_i5_0_i6_0",
    "NonPhysical": "ahcaleventreceiver00-DuplicatedEventCount",
    "Physical": "ahcaleventreceiver00-GoodEventsCount",
    "PreVeto": "tlureceiver00-PreVetoTriggers",
    "PostVeto": "tlureceiver00-PostVetoTriggers",
    "Recorded": "ahcaleventreceiver00-EventNumber",
}

UB_TO_FB_INV = 1.0e9
MAX_PLOT_BISECTION_RATE_HZ = 1.0e4
MAX_PLOT_BISECTION_RATIO_EVENTS_PER_FB = 1.0e12
TITLE_FONT_SIZE = 20
LABEL_FONT_SIZE = 16
TICK_FONT_SIZE = 15
LEGEND_FONT_SIZE = 12
ANNOTATION_FONT_SIZE = 11




def replace_environment_variables(input_str):
    secret_file = "faser-secret.json"
    env_config = {}
    if Path(secret_file).exists():
        with open(secret_file, encoding="utf-8") as f:
            env_config = json.load(f)

    def replace(match):
        var_name = match.group(1)
        if var_name in env_config:
            return str(env_config[var_name])
        env_value = os.getenv(var_name)
        if env_value is not None:
            return env_value
        return match.group(0)

    return re.sub(r"\$([A-Za-z0-9_]+)", replace, input_str)


def find_run_time(run):
    runnum = int(run)
    response = requests.get(
        f"https://faser-runinfo.app.cern.ch/cgibin/getRunInfo.py?runno={runnum}",
        timeout=20,
    )
    response.raise_for_status()
    run_json = response.json()

    if not run_json:
        raise RuntimeError(f"Could not find run information for run {runnum}")

    run_type = run_json.get("type", "unknown")
    if run_type != "AHCAL":
        raise RuntimeError(f"Run {runnum} is of type {run_type}, not AHCAL")

    start_time = run_json.get("starttime")
    end_time = run_json.get("stoptime")
    if not start_time:
        raise RuntimeError(f"Run {runnum} has no starttime")

    tstart = datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc)
    if end_time:
        tend = datetime.strptime(end_time, "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc)
    else:
        tend = datetime.now(timezone.utc)
    return tstart.timestamp(), tend.timestamp()


def make_client(database):
    return InfluxDBClient(
        host="dbod-faser-influx-prod.cern.ch",
        port=8080,
        username=replace_environment_variables("$INFLUXUSER"),
        password=replace_environment_variables("$INFLUXPW"),
        database=database,
        ssl=True,
        verify_ssl=True,
    )


def to_iso_utc(ts_s):
    return datetime.fromtimestamp(ts_s, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_time_to_ts_s(value):
    if value.endswith("Z"):
        value = value[:-1] + "+00:00"
    return datetime.fromisoformat(value).timestamp()


def influx_interval_str(seconds):
    if seconds <= 0:
        raise RuntimeError("Interval seconds must be > 0")
    if seconds >= 1.0 and abs(seconds - round(seconds)) < 1e-9:
        return f"{int(round(seconds))}s"
    ms = int(round(seconds * 1000.0))
    if ms <= 0:
        raise RuntimeError("Interval too small to represent in ms")
    return f"{ms}ms"


def query_counter_series(client, measurement, start_s, end_s, bin_seconds, field="value", verbose=False):
    interval = influx_interval_str(bin_seconds)
    query = (
        f'SELECT max("{field}") AS v FROM "{measurement}" '
        f"WHERE time >= '{to_iso_utc(start_s)}' and time <= '{to_iso_utc(end_s)}' "
        f"GROUP BY time({interval}) fill(previous)"
    )
    if verbose:
        print(f"[counter] {measurement}: {query}")
    result = client.query(query)
    points = list(result.get_points())
    return [(parse_time_to_ts_s(p["time"]), p.get("v")) for p in points]


def query_mean_series(client, measurement, start_s, end_s, bin_seconds, field="value", verbose=False):
    interval = influx_interval_str(bin_seconds)
    query = (
        f'SELECT mean("{field}") AS v FROM "{measurement}" '
        f"WHERE time >= '{to_iso_utc(start_s)}' and time <= '{to_iso_utc(end_s)}' "
        f"GROUP BY time({interval}) fill(null)"
    )
    if verbose:
        print(f"[mean] {measurement}: {query}")
    result = client.query(query)
    points = list(result.get_points())
    return [(parse_time_to_ts_s(p["time"]), p.get("v")) for p in points]


def to_rate_from_counter(counter_series, bin_seconds):
    rate = []
    prev = None
    for ts_s, value in counter_series:
        if value is None or not isinstance(value, (int, float)):
            rate.append((ts_s, math.nan))
            continue
        if prev is None or prev is None:
            rate.append((ts_s, math.nan))
            prev = value
            continue
        delta = value - prev
        if delta < 0:
            # Counter reset/rollover: drop this bin to avoid non-physical negative rate.
            rate.append((ts_s, math.nan))
        else:
            rate.append((ts_s, delta / float(bin_seconds)))
        prev = value
    return rate


def solve_true_rate_from_reco(r_reco, d_s, l_s):
    # Solve R_reco = R_true(1-exp(-L*R_true)) / (1-exp(-L*R_true)+D*R_true)
    # for R_true with bisection. D and L are in seconds.
    if not math.isfinite(r_reco) or r_reco <= 0.0:
        return 0.0
    if not math.isfinite(d_s) or d_s < 0.0 or not math.isfinite(l_s) or l_s <= 0.0:
        return r_reco

    def reco_from_true(r_true):
        one_minus_exp = 1.0 - math.exp(-l_s * r_true)
        denom = one_minus_exp + d_s * r_true
        if denom <= 0.0:
            return 0.0
        return (r_true * one_minus_exp) / denom

    lo = 0.0
    hi = max(1.0, r_reco * 2.0 + 1.0)
    for _ in range(100):
        if reco_from_true(hi) >= r_reco:
            break
        hi *= 2.0

    for _ in range(120):
        mid = 0.5 * (lo + hi)
        reco_mid = reco_from_true(mid)
        if reco_mid < r_reco:
            lo = mid
        else:
            hi = mid

    return 0.5 * (lo + hi)


def median_positive(values):
    pos = sorted(v for v in values if isinstance(v, (int, float)) and not math.isnan(v) and v > 0)
    if not pos:
        return math.nan
    n = len(pos)
    mid = n // 2
    if n % 2 == 1:
        return pos[mid]
    return 0.5 * (pos[mid - 1] + pos[mid])


def moving_average_nan(values, window_bins):
    if window_bins <= 1:
        return list(values)
    out = []
    half = window_bins // 2
    n = len(values)
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        vals = [v for v in values[lo:hi] if isinstance(v, (int, float)) and not math.isnan(v)]
        if vals:
            out.append(sum(vals) / len(vals))
        else:
            out.append(math.nan)
    return out


def gaussian_smooth_nan(values, window_bins):
    if window_bins <= 1:
        return list(values)
    radius = max(1, window_bins // 2)
    sigma = max(1.0, window_bins / 3.0)
    kernel = []
    for k in range(-radius, radius + 1):
        kernel.append(math.exp(-0.5 * (k / sigma) ** 2))

    n = len(values)
    out = []
    for i in range(n):
        wsum = 0.0
        vsum = 0.0
        for idx, k in enumerate(range(-radius, radius + 1)):
            j = i + k
            if j < 0 or j >= n:
                continue
            v = values[j]
            if not isinstance(v, (int, float)) or math.isnan(v):
                continue
            w = kernel[idx]
            wsum += w
            vsum += w * v
        out.append(vsum / wsum if wsum > 0.0 else math.nan)
    return out


def correction_window_values(values, window_bins, method="gaussian"):
    if method == "gaussian":
        return gaussian_smooth_nan(values, window_bins)
    return moving_average_nan(values, window_bins)


def build_ratio_lumi_denominator(lumi_rate_values, smooth_bins, min_frac):
    smoothed = moving_average_nan(lumi_rate_values, smooth_bins)
    med = median_positive(smoothed)
    if math.isnan(med):
        return [math.nan] * len(smoothed), math.nan
    floor = min_frac * med
    denom = []
    for v in smoothed:
        if not isinstance(v, (int, float)) or math.isnan(v) or v < floor:
            denom.append(math.nan)
        else:
            denom.append(v)
    return denom, floor


def mask_unphysical_for_plot(values, max_abs_value):
    masked = []
    for value in values:
        if not isinstance(value, (int, float)) or not math.isfinite(value):
            masked.append(math.nan)
        elif abs(value) > max_abs_value:
            masked.append(math.nan)
        else:
            masked.append(value)
    return masked


def clip_unphysical_for_plot(values, max_abs_value):
    clipped = []
    for value in values:
        if not isinstance(value, (int, float)) or not math.isfinite(value):
            clipped.append(math.nan)
        elif value > max_abs_value:
            clipped.append(max_abs_value)
        elif value < -max_abs_value:
            clipped.append(-max_abs_value)
        else:
            clipped.append(value)
    return clipped


def convert_lumi_to_ub(series):
    converted = []
    for ts_s, value in series:
        if isinstance(value, (int, float)) and not math.isnan(value):
            converted.append((ts_s, value))
        else:
            converted.append((ts_s, math.nan))
    return converted


def align_by_timestamp(series_dict):
    all_ts = set()
    for series in series_dict.values():
        for ts_s, _ in series:
            all_ts.add(ts_s)
    ts_sorted = sorted(all_ts)
    aligned = {key: {} for key in series_dict}
    for key, series in series_dict.items():
        for ts_s, value in series:
            aligned[key][ts_s] = value
    return ts_sorted, aligned


def discover_lumi_measurement(client, verbose=False):
    query = 'SHOW MEASUREMENTS WITH MEASUREMENT =~ /(?i)dip/'
    result = client.query(query)
    measurements = []
    for item in result.raw.get("series", []):
        for row in item.get("values", []):
            measurements.append(row[0])
    if verbose:
        print(f"DIP candidates: {measurements}")
    if not measurements:
        raise RuntimeError("No dip-like measurements found in DcsData")

    preferred = [m for m in measurements if "lumi" in m.lower() or "luminosity" in m.lower()]
    if len(preferred) == 1:
        return preferred[0]
    if len(preferred) > 1:
        raise RuntimeError(f"Multiple lumi-like dip measurements found. Use --lumi-measurement. Candidates: {preferred}")
    if len(measurements) == 1:
        return measurements[0]
    raise RuntimeError(f"Multiple dip measurements found. Use --lumi-measurement. Candidates: {measurements}")


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
            "luminosity_ub_inv_s",
            "luminosity_smoothed_ub_inv_s",
            "luminosity_for_ratio_ub_inv_s",
            "corrected_any_vetoed_rate_Hz",
            "corrected_any_vetoed_rate_bisection_Hz",
            "corrected_over_luminosity_events_per_fb_inv",
            "corrected_over_luminosity_bisection_events_per_fb_inv",
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
                fmt(aligned["Luminosity"][ts_s]),
                fmt(aligned["Luminosity_smoothed"][ts_s]),
                fmt(aligned["Luminosity_for_ratio"][ts_s]),
                fmt(corrected_rate[i]),
                fmt(corrected_rate_bisection[i]),
                fmt(ratio_series[i]),
                fmt(ratio_series_bisection[i]),
            ]
            f.write(",".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Fetch trigger/lumi time series for one run and make rate plots.")
    parser.add_argument("--run", required=True, type=int, help="AHCAL run number")
    parser.add_argument("--bin-seconds", type=float, default=0.2, help="Bin width in seconds for luminosity InfluxDB time binning")
    parser.add_argument("--trigger-bin-seconds", type=float, default=1.0, help="Bin width in seconds for trigger counter time binning")
    parser.add_argument("--correction-window-seconds", type=float, default=10.0, help="Window size n [s] used in corrected-rate calculation")
    parser.add_argument("--correction-smooth-method", choices=["gaussian", "moving-average"], default="gaussian", help="Smoothing method over n-second window before correction")
    parser.add_argument("--lumi-measurement", default=None, help="Measurement name in DcsData for luminosity")
    parser.add_argument("--lumi-field", default="Luminosity", help="Field name in DcsData luminosity measurement")
    parser.add_argument("--lumi-plot-smooth-bins", type=int, default=11, help="Moving-average window (bins) for displayed luminosity")
    parser.add_argument("--ratio-lumi-smooth-bins", type=int, default=5, help="Moving-average window (bins) for luminosity denominator in ratio")
    parser.add_argument("--ratio-lumi-min-frac", type=float, default=0.05, help="Minimum denominator floor as fraction of median smoothed luminosity")
    parser.add_argument("--ratio-lumi-min-abs", type=float, default=1.0e6, help="Absolute near-zero luminosity threshold in ub^-1/s; ratio is masked when luminosity is below this value")
    parser.add_argument("--assumed-deadtime-ms", type=float, default=1.22, help="Assumed deadtime D in ms for bisection correction")
    parser.add_argument("--assumed-auto-triggered-ms", type=float, default=3.9, help="Assumed auto-triggered window L in ms for bisection correction")
    parser.add_argument("--assumed-deadtime-ms-first", type=float, default=1.348, help="Assumed deadtime D in ms for first correction; if 0, uses --assumed-deadtime-ms for first correction as well")
    parser.add_argument("--plot-start-seconds", type=float, default=0.0, help="Plot start time in seconds from run start")
    parser.add_argument("--plot-end-seconds", type=float, default=None, help="Plot end time in seconds from run start; defaults to run end")
    parser.add_argument("--output-dir", default=".", help="Output directory for plots and CSV/JSON")
    parser.add_argument("--verbose", "-v", action="store_true", help="Debug output")
    args = parser.parse_args()

    if args.bin_seconds <= 0:
        raise RuntimeError("--bin-seconds must be > 0")
    if args.trigger_bin_seconds <= 0:
        raise RuntimeError("--trigger-bin-seconds must be > 0")
    if args.correction_window_seconds <= 0:
        raise RuntimeError("--correction-window-seconds must be > 0")
    if args.lumi_plot_smooth_bins <= 0:
        raise RuntimeError("--lumi-plot-smooth-bins must be > 0")
    if args.ratio_lumi_smooth_bins <= 0:
        raise RuntimeError("--ratio-lumi-smooth-bins must be > 0")
    if args.ratio_lumi_min_frac < 0:
        raise RuntimeError("--ratio-lumi-min-frac must be >= 0")
    if args.ratio_lumi_min_abs < 0:
        raise RuntimeError("--ratio-lumi-min-abs must be >= 0")
    if args.assumed_deadtime_ms < 0:
        raise RuntimeError("--assumed-deadtime-ms must be >= 0")
    if args.assumed_auto_triggered_ms <= 0:
        raise RuntimeError("--assumed-auto-triggered-ms must be > 0")
    if args.plot_start_seconds < 0:
        raise RuntimeError("--plot-start-seconds must be >= 0")

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

    lumi_measurement = args.lumi_measurement or discover_lumi_measurement(dcs_client, verbose=args.verbose)
    print(f"Run {args.run}: start={to_iso_utc(start_s)} end={to_iso_utc(end_s)}")
    print(
        f"Luminosity bin={args.bin_seconds}s, trigger bin={args.trigger_bin_seconds}s, "
        f"correction window n={args.correction_window_seconds}s, method={args.correction_smooth_method}"
    )
    print(f"Trigger DB={trigger_db_name}, Lumi DB=DcsData, Lumi measurement={lumi_measurement}")

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

    lumi_mean_series = query_mean_series(
        dcs_client,
        lumi_measurement,
        start_s,
        end_s,
        args.bin_seconds,
        field=args.lumi_field,
        verbose=args.verbose,
    )
    if not lumi_mean_series:
        raise RuntimeError(f"No luminosity points found for {lumi_measurement} in run window")
    lumi_mean_series = convert_lumi_to_ub(lumi_mean_series)

    series_dict = dict(counter_rate_series)
    series_dict["Luminosity"] = lumi_mean_series
    ts_sorted, aligned = align_by_timestamp(series_dict)

    for key in series_dict:
        data = [math.nan if aligned[key].get(ts) is None else aligned[key].get(ts) for ts in ts_sorted]
        for i, ts_s in enumerate(ts_sorted):
            aligned[key][ts_s] = data[i]

    lumi_values = [aligned["Luminosity"][ts_s] for ts_s in ts_sorted]
    lumi_plot_smoothed = moving_average_nan(lumi_values, args.lumi_plot_smooth_bins)
    for i, ts_s in enumerate(ts_sorted):
        aligned.setdefault("Luminosity_smoothed", {})[ts_s] = lumi_plot_smoothed[i]

    lumi_ratio_denom, lumi_ratio_floor = build_ratio_lumi_denominator(
        lumi_values,
        args.ratio_lumi_smooth_bins,
        args.ratio_lumi_min_frac,
    )
    for i, ts_s in enumerate(ts_sorted):
        aligned.setdefault("Luminosity_for_ratio", {})[ts_s] = lumi_ratio_denom[i]

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

    corrected_any_rate = []
    corrected_any_rate_bisection = []
    ratio_series_events_per_fb = []
    ratio_series_bisection_events_per_fb = []
    linear_upper_limit_bins = []
    bisection_upper_limit_bins = []
    invalid_bins = 0
    invalid_ratio_bins = 0
    near_zero_lumi_bins = 0
    for i, ts_s in enumerate(ts_sorted):
        any_vetoed = any_for_corr[i]
        non_vetoed = non_for_corr[i]
        non_physical = np_for_corr[i]
        physical = phys_for_corr[i]
        raw_luminosity = aligned["Luminosity"][ts_s]
        luminosity_for_ratio = aligned["Luminosity_for_ratio"][ts_s]

        if any(math.isnan(x) for x in [any_vetoed, non_physical, physical]):
            corrected = math.nan
            linear_hit_upper_limit = False
        else:
            denom = 1.0 - ( (args.assumed_deadtime_ms *(non_physical) + args.assumed_deadtime_ms_first*(physical) )/ 1000.0 )
            if denom <= 0.0:
                corrected = math.nan
                invalid_bins += 1
                linear_hit_upper_limit = False
            else:
                corrected = any_vetoed / denom
                linear_hit_upper_limit = False
                if corrected > MAX_PLOT_BISECTION_RATE_HZ:
                    linear_hit_upper_limit = True
                    corrected = math.nan
        corrected_any_rate.append(corrected)
        linear_upper_limit_bins.append(linear_hit_upper_limit)

        if any(math.isnan(x) for x in [any_vetoed, non_vetoed]):
            corrected_bisection = math.nan
            bisection_hit_upper_limit = False
        else:
            # Match VetoAnaAlg usage: infer total reconstructed trigger rate and
            # scale AnyVetoed by true/reco total from bisection solution.
            reco_non_vetoed = non_for_corr[i]
            if math.isnan(reco_non_vetoed):
                corrected_bisection = math.nan
                bisection_hit_upper_limit = False
                corrected_any_rate_bisection.append(corrected_bisection)
                is_near_zero_raw_lumi = (
                    isinstance(raw_luminosity, (int, float)) and not math.isnan(raw_luminosity) and raw_luminosity < args.ratio_lumi_min_abs
                )
                if is_near_zero_raw_lumi:
                    ratio_series_events_per_fb.append(math.nan)
                    ratio_series_bisection_events_per_fb.append(math.nan)
                    invalid_ratio_bins += 1
                    near_zero_lumi_bins += 1
                elif math.isnan(corrected) or math.isnan(luminosity_for_ratio) or luminosity_for_ratio == 0:
                    ratio_series_events_per_fb.append(math.nan)
                    ratio_series_bisection_events_per_fb.append(math.nan)
                    invalid_ratio_bins += 1
                else:
                    ratio_events_per_ub = corrected / luminosity_for_ratio
                    ratio_series_events_per_fb.append(ratio_events_per_ub * UB_TO_FB_INV)
                    ratio_series_bisection_events_per_fb.append(math.nan)
                continue

            reco_total = any_vetoed + reco_non_vetoed
            true_total = solve_true_rate_from_reco(
                reco_total,
                args.assumed_deadtime_ms / 1000.0,
                args.assumed_auto_triggered_ms / 1000.0,
            )
            bisection_hit_upper_limit = False
            if true_total > MAX_PLOT_BISECTION_RATE_HZ:
                # Sanity check to avoid unphysically large correction from bisection solution.
                bisection_hit_upper_limit = True
                true_total = MAX_PLOT_BISECTION_RATE_HZ
                print(f"Warning: true total rate from bisection is very large ({true_total} Hz) , reco_total={reco_total} Hz")
            if reco_total <= 0.0:
                corrected_bisection = math.nan
            else:
                corrected_bisection = any_vetoed * (true_total / reco_total)
        corrected_any_rate_bisection.append(corrected_bisection)
        bisection_upper_limit_bins.append(bisection_hit_upper_limit)

        is_near_zero_raw_lumi = (
            isinstance(raw_luminosity, (int, float)) and not math.isnan(raw_luminosity) and raw_luminosity < args.ratio_lumi_min_abs
        )
        if is_near_zero_raw_lumi:
            ratio_series_events_per_fb.append(math.nan)
            ratio_series_bisection_events_per_fb.append(math.nan)
            invalid_ratio_bins += 1
            near_zero_lumi_bins += 1
        else:
            if math.isnan(corrected) or math.isnan(luminosity_for_ratio) or luminosity_for_ratio == 0:
                ratio_series_events_per_fb.append(math.nan)
                invalid_ratio_bins += 1
            else:
                ratio_events_per_ub = corrected / luminosity_for_ratio
                ratio_series_events_per_fb.append(ratio_events_per_ub * UB_TO_FB_INV)

            if math.isnan(corrected_bisection) or math.isnan(luminosity_for_ratio) or luminosity_for_ratio == 0:
                ratio_series_bisection_events_per_fb.append(math.nan)
            else:
                ratio_bisection_events_per_ub = corrected_bisection / luminosity_for_ratio
                ratio_series_bisection_events_per_fb.append(ratio_bisection_events_per_ub * UB_TO_FB_INV)

    elapsed = [ts - start_s for ts in ts_sorted]

    plot_end_seconds = args.plot_end_seconds
    if plot_end_seconds is None:
        plot_end_seconds = elapsed[-1] if elapsed else 0.0
    if plot_end_seconds <= args.plot_start_seconds:
        raise RuntimeError("--plot-end-seconds must be greater than --plot-start-seconds")

    plot_indices = [i for i, e in enumerate(elapsed) if args.plot_start_seconds <= e <= plot_end_seconds]
    if not plot_indices:
        raise RuntimeError("No points fall inside the requested plot time range")

    plot_elapsed = [elapsed[i] for i in plot_indices]
    plot_corrected_any_rate = [corrected_any_rate[i] for i in plot_indices]
    plot_corrected_any_rate_bisection = [corrected_any_rate_bisection[i] for i in plot_indices]
    plot_any_veto_raw = [aligned["AnyVetoed_rate"][ts_sorted[i]] for i in plot_indices]
    plot_any_veto_smoothed = [any_for_corr[i] for i in plot_indices]
    plot_lumi_for_ratio = [aligned["Luminosity_for_ratio"][ts_sorted[i]] for i in plot_indices]
    plot_pre_ratio_raw = []
    plot_pre_ratio_smoothed = []
    for i in range(len(plot_indices)):
        denom = plot_lumi_for_ratio[i]
        raw_rate = plot_any_veto_raw[i]
        smooth_rate = plot_any_veto_smoothed[i]
        if not isinstance(denom, (int, float)) or math.isnan(denom) or denom == 0.0:
            plot_pre_ratio_raw.append(math.nan)
            plot_pre_ratio_smoothed.append(math.nan)
        else:
            if isinstance(raw_rate, (int, float)) and not math.isnan(raw_rate):
                plot_pre_ratio_raw.append((raw_rate / denom) * UB_TO_FB_INV)
            else:
                plot_pre_ratio_raw.append(math.nan)
            if isinstance(smooth_rate, (int, float)) and not math.isnan(smooth_rate):
                plot_pre_ratio_smoothed.append((smooth_rate / denom) * UB_TO_FB_INV)
            else:
                plot_pre_ratio_smoothed.append(math.nan)
    plot_lumi_raw = [aligned["Luminosity"][ts_sorted[i]] for i in plot_indices]
    plot_lumi_smooth = [aligned["Luminosity_smoothed"][ts_sorted[i]] for i in plot_indices]
    plot_ratio = [ratio_series_events_per_fb[i] for i in plot_indices]
    plot_ratio_bisection = [ratio_series_bisection_events_per_fb[i] for i in plot_indices]
    plot_linear_hit_upper_limit = [linear_upper_limit_bins[i] for i in plot_indices]
    plot_bisection_hit_upper_limit = [bisection_upper_limit_bins[i] for i in plot_indices]
    plot_corrected_any_rate = mask_unphysical_for_plot(plot_corrected_any_rate, MAX_PLOT_BISECTION_RATE_HZ)
    plot_corrected_any_rate_bisection = mask_unphysical_for_plot(plot_corrected_any_rate_bisection, MAX_PLOT_BISECTION_RATE_HZ)
    plot_any_veto_raw = clip_unphysical_for_plot(plot_any_veto_raw, MAX_PLOT_BISECTION_RATE_HZ)
    plot_any_veto_smoothed = clip_unphysical_for_plot(plot_any_veto_smoothed, MAX_PLOT_BISECTION_RATE_HZ)
    plot_pre_ratio_raw = mask_unphysical_for_plot(plot_pre_ratio_raw, MAX_PLOT_BISECTION_RATIO_EVENTS_PER_FB)
    plot_pre_ratio_smoothed = mask_unphysical_for_plot(plot_pre_ratio_smoothed, MAX_PLOT_BISECTION_RATIO_EVENTS_PER_FB)
    plot_ratio = mask_unphysical_for_plot(plot_ratio, MAX_PLOT_BISECTION_RATIO_EVENTS_PER_FB)
    plot_ratio_bisection = mask_unphysical_for_plot(plot_ratio_bisection, MAX_PLOT_BISECTION_RATIO_EVENTS_PER_FB)
    has_linear_correction = any((not math.isnan(v)) for v in plot_corrected_any_rate)
    has_linear_ratio = any((not math.isnan(v)) for v in plot_ratio)
    has_linear_upper_limit_hits = any(plot_linear_hit_upper_limit)
    has_plottable_bisection_correction = any((not math.isnan(v)) for v in plot_corrected_any_rate_bisection)
    has_plottable_bisection_ratio = any((not math.isnan(v)) for v in plot_ratio_bisection)
    has_bisection_upper_limit_hits = any(plot_bisection_hit_upper_limit)
    has_any_veto_raw = any((not math.isnan(v)) for v in plot_any_veto_raw)
    has_any_veto_smoothed = any((not math.isnan(v)) for v in plot_any_veto_smoothed)
    has_pre_ratio_raw = any((not math.isnan(v)) for v in plot_pre_ratio_raw)
    has_pre_ratio_smoothed = any((not math.isnan(v)) for v in plot_pre_ratio_smoothed)

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(12, 10),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 2]},
    )

    ax_top_right = ax_top.twinx()
    if has_linear_correction:
        ax_top.plot(plot_elapsed, plot_corrected_any_rate, color="tab:blue", linewidth=2.2, label="Corrected AnyVetoed rate")
    if has_linear_upper_limit_hits:
        cap_times = [t for t, hit in zip(plot_elapsed, plot_linear_hit_upper_limit) if hit]
        ax_top.plot(
            cap_times,
            [0.97] * len(cap_times),
            color="tab:orange",
            marker="o",
            linestyle="None",
            transform=ax_top.get_xaxis_transform(),
            label="linear correction hit upper limit",
            zorder=5,
        )
    # if has_plottable_bisection_correction:
    #     ax_top.plot(plot_elapsed, plot_corrected_any_rate_bisection, color="tab:purple", linewidth=2.0, linestyle="--", label="Corrected AnyVetoed rate (bisection)")
    # if has_bisection_upper_limit_hits:
    #     cap_times = [t for t, hit in zip(plot_elapsed, plot_bisection_hit_upper_limit) if hit]
    #     ax_top.plot(
    #         cap_times,
    #         [0.97] * len(cap_times),
    #         color="tab:orange",
    #         marker="v",
    #         linestyle="None",
    #         transform=ax_top.get_xaxis_transform(),
    #         label="bisection hit upper limit",
    #         zorder=5,
    #     )
    #     ax_top.text(
    #         0.99,
    #         0.94,
    #         f"orange markers: bisection clipped at >={MAX_PLOT_BISECTION_RATE_HZ:.0f} Hz",
    #         transform=ax_top.transAxes,
    #         ha="right",
    #         va="top",
    #         color="tab:orange",
    #         fontsize=ANNOTATION_FONT_SIZE,
    #     )
    ax_top_right.plot(plot_elapsed, plot_lumi_raw, color="tab:red", linewidth=1.0, alpha=0.22, label="Luminosity (raw)")
    ax_top_right.plot(plot_elapsed, plot_lumi_smooth, color="tab:red", linewidth=2.1, alpha=0.95, label="Luminosity (smoothed)")
    ax_top.set_ylabel("Corrected AnyVetoed Rate [Hz]", color="tab:blue", fontsize=LABEL_FONT_SIZE, labelpad=10)
    ax_top_right.set_ylabel("Luminosity [ub^-1/s]", color="tab:red", fontsize=LABEL_FONT_SIZE, labelpad=16)
    # ax_top.set_ylim(bottom=0.0)
    # ax_top_right.set_ylim(bottom=0.0)
    ax_top.set_yscale("log")
    ax_top.set_ylim(bottom=100)
    ax_top_right.set_yscale("log")
    ax_top_right.set_ylim(bottom=100)
    ax_top.grid(True, alpha=0.25)
    ax_top.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax_top_right.tick_params(axis="y", labelsize=TICK_FONT_SIZE, pad=4)
    lines_left, labels_left = ax_top.get_legend_handles_labels()
    lines_right, labels_right = ax_top_right.get_legend_handles_labels()
    ax_top.legend(lines_left + lines_right, labels_left + labels_right, loc="upper right", fontsize=LEGEND_FONT_SIZE)
    ax_top.set_title(f"Run {args.run}: corrected AnyVetoed rate and luminosity", fontsize=TITLE_FONT_SIZE)

    if has_linear_ratio:
        ax_bottom.plot(plot_elapsed, plot_ratio, color="tab:green", linewidth=2.0, label="Corrected AnyVetoed rate / luminosity")
    if has_linear_upper_limit_hits:
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
    # if has_plottable_bisection_ratio:
    #     ax_bottom.plot(plot_elapsed, plot_ratio_bisection, color="tab:purple", linewidth=1.8, linestyle="--", label="Corrected AnyVetoed rate / luminosity (bisection)")
    # if has_bisection_upper_limit_hits:
    #     cap_times = [t for t, hit in zip(plot_elapsed, plot_bisection_hit_upper_limit) if hit]
    #     ax_bottom.plot(
    #         cap_times,
    #         [0.96] * len(cap_times),
    #         color="tab:orange",
    #         marker="v",
    #         linestyle="None",
    #         transform=ax_bottom.get_xaxis_transform(),
    #         label="bisection hit upper limit",
    #         zorder=5,
    #     )
    #     ax_bottom.text(
    #         0.99,
    #         0.94,
    #         "orange markers: bisection clipped at upper limit",
    #         transform=ax_bottom.transAxes,
    #         ha="right",
    #         va="top",
    #         color="tab:orange",
    #         fontsize=ANNOTATION_FONT_SIZE,
    #     )
    ax_bottom.set_xlabel("Elapsed time from run start [s]", fontsize=LABEL_FONT_SIZE)
    ax_bottom.set_ylabel("Rate / Luminosity [events/fb^-1]", fontsize=LABEL_FONT_SIZE)
    ax_bottom.set_ylim(bottom=0.0)
    ax_bottom.grid(True, alpha=0.25)
    ax_bottom.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax_bottom.legend(loc="upper right", fontsize=LEGEND_FONT_SIZE)

    fig.tight_layout()
    plot_path = output_dir / f"run_{args.run}_stacked_plots.png"
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)

    fig_overlay, (ax_top_overlay, ax_bottom_overlay) = plt.subplots(
        2,
        1,
        figsize=(12, 10),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 2]},
    )

    ax_top_overlay_right = ax_top_overlay.twinx()
    if has_linear_correction:
        ax_top_overlay.plot(plot_elapsed, plot_corrected_any_rate, color="tab:blue", linewidth=2.2, label="Corrected AnyVetoed rate")
    if has_any_veto_raw:
        ax_top_overlay.plot(
            plot_elapsed,
            plot_any_veto_raw,
            color="0.45",
            linewidth=1.1,
            alpha=0.28,
            label="AnyVetoed rate (raw, pre-correction)",
        )
    if has_any_veto_smoothed:
        ax_top_overlay.plot(
            plot_elapsed,
            plot_any_veto_smoothed,
            color="0.30",
            linewidth=1.9,
            alpha=0.90,
            linestyle=":",
            label="AnyVetoed rate (smoothed, pre-correction)",
        )
    if has_linear_upper_limit_hits:
        cap_times = [t for t, hit in zip(plot_elapsed, plot_linear_hit_upper_limit) if hit]
        ax_top_overlay.plot(
            cap_times,
            [0.97] * len(cap_times),
            color="tab:orange",
            marker="o",
            linestyle="None",
            transform=ax_top_overlay.get_xaxis_transform(),
            label="linear correction hit upper limit",
            zorder=5,
        )
    # if has_plottable_bisection_correction:
    #     ax_top_overlay.plot(
    #         plot_elapsed,
    #         plot_corrected_any_rate_bisection,
    #         color="tab:purple",
    #         linewidth=2.0,
    #         linestyle="--",
    #         label="Corrected AnyVetoed rate (bisection)",
    #     )
    # if has_bisection_upper_limit_hits:
    #     cap_times = [t for t, hit in zip(plot_elapsed, plot_bisection_hit_upper_limit) if hit]
    #     ax_top_overlay.plot(
    #         cap_times,
    #         [0.97] * len(cap_times),
    #         color="tab:orange",
    #         marker="v",
    #         linestyle="None",
    #         transform=ax_top_overlay.get_xaxis_transform(),
    #         label="bisection hit upper limit",
    #         zorder=5,
    #     )
    #     ax_top_overlay.text(
    #         0.99,
    #         0.94,
    #         f"orange markers: bisection clipped at >={MAX_PLOT_BISECTION_RATE_HZ:.0f} Hz",
    #         transform=ax_top_overlay.transAxes,
    #         ha="right",
    #         va="top",
    #         color="tab:orange",
    #         fontsize=ANNOTATION_FONT_SIZE,
    #     )
    ax_top_overlay_right.plot(plot_elapsed, plot_lumi_raw, color="tab:red", linewidth=1.0, alpha=0.22, label="Luminosity (raw)")
    ax_top_overlay_right.plot(plot_elapsed, plot_lumi_smooth, color="tab:red", linewidth=2.1, alpha=0.95, label="Luminosity (smoothed)")
    ax_top_overlay.set_ylabel("Veto Rate [Hz]", color="tab:blue", fontsize=LABEL_FONT_SIZE, labelpad=10)
    ax_top_overlay_right.set_ylabel("Luminosity [ub^-1/s]", color="tab:red", fontsize=LABEL_FONT_SIZE, labelpad=16)
    ax_top_overlay.set_yscale("log")
    ax_top_overlay.set_ylim(bottom=100)
    ax_top_overlay_right.set_yscale("log")
    ax_top_overlay_right.set_ylim(bottom=100)
    ax_top_overlay.grid(True, alpha=0.25)
    ax_top_overlay.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax_top_overlay_right.tick_params(axis="y", labelsize=TICK_FONT_SIZE, pad=4)
    lines_left_overlay, labels_left_overlay = ax_top_overlay.get_legend_handles_labels()
    lines_right_overlay, labels_right_overlay = ax_top_overlay_right.get_legend_handles_labels()
    ax_top_overlay.legend(
        lines_left_overlay + lines_right_overlay,
        labels_left_overlay + labels_right_overlay,
        loc="upper right",
        fontsize=LEGEND_FONT_SIZE,
    )
    ax_top_overlay.set_title(
        f"Run {args.run}: corrected + pre-correction AnyVetoed rate and luminosity",
        fontsize=TITLE_FONT_SIZE,
    )

    if has_linear_ratio:
        ax_bottom_overlay.plot(plot_elapsed, plot_ratio, color="tab:green", linewidth=2.0, label="Corrected AnyVetoed rate / luminosity")
    if has_pre_ratio_raw:
        ax_bottom_overlay.plot(
            plot_elapsed,
            plot_pre_ratio_raw,
            color="0.45",
            linewidth=1.1,
            alpha=0.30,
            label="AnyVetoed/Luminosity (raw, pre-correction)",
        )
    if has_pre_ratio_smoothed:
        ax_bottom_overlay.plot(
            plot_elapsed,
            plot_pre_ratio_smoothed,
            color="0.30",
            linewidth=1.9,
            alpha=0.90,
            linestyle=":",
            label="AnyVetoed/Luminosity (smoothed, pre-correction)",
        )
    if has_linear_upper_limit_hits:
        cap_times = [t for t, hit in zip(plot_elapsed, plot_linear_hit_upper_limit) if hit]
        ax_bottom_overlay.plot(
            cap_times,
            [0.96] * len(cap_times),
            color="tab:orange",
            marker="o",
            linestyle="None",
            transform=ax_bottom_overlay.get_xaxis_transform(),
            label="linear correction hit upper limit",
            zorder=5,
        )
    # if has_plottable_bisection_ratio:
    #     ax_bottom_overlay.plot(
    #         plot_elapsed,
    #         plot_ratio_bisection,
    #         color="tab:purple",
    #         linewidth=1.8,
    #         linestyle="--",
    #         label="Corrected AnyVetoed rate / luminosity (bisection)",
    #     )
    # if has_bisection_upper_limit_hits:
    #     cap_times = [t for t, hit in zip(plot_elapsed, plot_bisection_hit_upper_limit) if hit]
    #     ax_bottom_overlay.plot(
    #         cap_times,
    #         [0.96] * len(cap_times),
    #         color="tab:orange",
    #         marker="v",
    #         linestyle="None",
    #         transform=ax_bottom_overlay.get_xaxis_transform(),
    #         label="bisection hit upper limit",
    #         zorder=5,
    #     )
    #     ax_bottom_overlay.text(
    #         0.99,
    #         0.94,
    #         "orange markers: bisection clipped at upper limit",
    #         transform=ax_bottom_overlay.transAxes,
    #         ha="right",
    #         va="top",
    #         color="tab:orange",
    #         fontsize=ANNOTATION_FONT_SIZE,
    #     )
    ax_bottom_overlay.set_xlabel("Elapsed time from run start [s]", fontsize=LABEL_FONT_SIZE)
    ax_bottom_overlay.set_ylabel("Rate / Luminosity [events/fb^-1]", fontsize=LABEL_FONT_SIZE)
    ax_bottom_overlay.set_ylim(bottom=0.0)
    ax_bottom_overlay.grid(True, alpha=0.25)
    ax_bottom_overlay.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax_bottom_overlay.legend(loc="upper right", fontsize=LEGEND_FONT_SIZE)

    fig_overlay.tight_layout()
    corrected_overlay_plot_path = output_dir / f"run_{args.run}_stacked_plots_with_precorrection.png"
    fig_overlay.savefig(corrected_overlay_plot_path, dpi=180)
    plt.close(fig_overlay)

    plot_compare_rate = plot_corrected_any_rate
    plot_compare_bisection = plot_corrected_any_rate_bisection
    plot_compare_pre = plot_any_veto_smoothed
    plot_compare_ratio_bisection_over_rate = []
    plot_compare_ratio_rate_over_pre = []
    plot_compare_ratio_bisection_over_pre = []
    for rate_value, bisection_value, pre_value in zip(plot_compare_rate, plot_compare_bisection, plot_compare_pre):
        if not isinstance(rate_value, (int, float)) or not math.isfinite(rate_value) or rate_value == 0.0:
            plot_compare_ratio_bisection_over_rate.append(math.nan)
        else:
            if not isinstance(bisection_value, (int, float)) or not math.isfinite(bisection_value):
                plot_compare_ratio_bisection_over_rate.append(math.nan)
            else:
                plot_compare_ratio_bisection_over_rate.append(bisection_value / rate_value)

        if not isinstance(pre_value, (int, float)) or not math.isfinite(pre_value) or pre_value == 0.0:
            plot_compare_ratio_rate_over_pre.append(math.nan)
            plot_compare_ratio_bisection_over_pre.append(math.nan)
        else:
            if not isinstance(rate_value, (int, float)) or not math.isfinite(rate_value):
                plot_compare_ratio_rate_over_pre.append(math.nan)
            else:
                plot_compare_ratio_rate_over_pre.append(rate_value / pre_value)

            if not isinstance(bisection_value, (int, float)) or not math.isfinite(bisection_value):
                plot_compare_ratio_bisection_over_pre.append(math.nan)
            else:
                plot_compare_ratio_bisection_over_pre.append(bisection_value / pre_value)

    has_compare_rate = any((not math.isnan(v)) for v in plot_compare_rate)
    has_compare_bisection = any((not math.isnan(v)) for v in plot_compare_bisection)
    has_compare_pre = any((not math.isnan(v)) for v in plot_compare_pre)
    has_compare_ratio_bisection_over_rate = any((not math.isnan(v)) for v in plot_compare_ratio_bisection_over_rate)
    has_compare_ratio_rate_over_pre = any((not math.isnan(v)) for v in plot_compare_ratio_rate_over_pre)
    has_compare_ratio_bisection_over_pre = any((not math.isnan(v)) for v in plot_compare_ratio_bisection_over_pre)

    fig_compare, (ax_compare_top, ax_compare_bottom) = plt.subplots(
        2,
        1,
        figsize=(12, 10),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 2]},
    )

    if has_compare_rate:
        ax_compare_top.plot(
            plot_elapsed,
            plot_compare_rate,
            color="tab:blue",
            linewidth=2.2,
            label="rate-based correction",
        )
    if has_compare_bisection:
        ax_compare_top.plot(
            plot_elapsed,
            plot_compare_bisection,
            color="tab:purple",
            linewidth=2.0,
            linestyle="--",
            label="bisection correction",
        )
    if has_compare_pre:
        ax_compare_top.plot(
            plot_elapsed,
            plot_compare_pre,
            color="0.30",
            linewidth=1.9,
            alpha=0.90,
            linestyle=":",
            label="pre-correction (smoothed)",
        )
    ax_compare_top.set_ylabel("Corrected Rate [Hz]", fontsize=LABEL_FONT_SIZE)
    ax_compare_top.set_title(f"Run {args.run}: bisection vs rate correction", fontsize=TITLE_FONT_SIZE)
    ax_compare_top.set_ylim(bottom=0.0)
    ax_compare_top.grid(True, alpha=0.25)
    ax_compare_top.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    if has_compare_rate or has_compare_bisection:
        ax_compare_top.legend(loc="upper right", fontsize=LEGEND_FONT_SIZE)

    if has_compare_ratio_bisection_over_rate:
        ax_compare_bottom.plot(
            plot_elapsed,
            plot_compare_ratio_bisection_over_rate,
            color="tab:green",
            linewidth=2.0,
            label="bisection / rate",
        )
    if has_compare_ratio_rate_over_pre:
        ax_compare_bottom.plot(
            plot_elapsed,
            plot_compare_ratio_rate_over_pre,
            color="tab:blue",
            linewidth=1.8,
            alpha=0.95,
            label="rate / pre-correction",
        )
    if has_compare_ratio_bisection_over_pre:
        ax_compare_bottom.plot(
            plot_elapsed,
            plot_compare_ratio_bisection_over_pre,
            color="tab:purple",
            linewidth=1.8,
            alpha=0.95,
            linestyle="--",
            label="bisection / pre-correction",
        )
    ax_compare_bottom.axhline(1.0, color="0.5", linestyle="--", linewidth=1.0, alpha=0.8)
    ax_compare_bottom.set_xlabel("Elapsed time from run start [s]", fontsize=LABEL_FONT_SIZE)
    ax_compare_bottom.set_ylabel("Ratio", fontsize=LABEL_FONT_SIZE)
    ax_compare_bottom.set_ylim(bottom=0.0)
    ax_compare_bottom.grid(True, alpha=0.25)
    ax_compare_bottom.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    if has_compare_ratio_bisection_over_rate or has_compare_ratio_rate_over_pre or has_compare_ratio_bisection_over_pre:
        ax_compare_bottom.legend(loc="upper right", fontsize=LEGEND_FONT_SIZE)

    fig_compare.tight_layout()
    compare_plot_path = output_dir / f"run_{args.run}_bisection_vs_rate_comparison.png"
    fig_compare.savefig(compare_plot_path, dpi=180)
    plt.close(fig_compare)

    fig_any, (ax_any, ax_any_ratio) = plt.subplots(
        2,
        1,
        figsize=(12, 10),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 2]},
    )
    if has_any_veto_raw:
        ax_any.plot(plot_elapsed, plot_any_veto_raw, color="tab:blue", linewidth=1.2, alpha=0.28, label="AnyVetoed rate (raw, pre-correction)")
    if has_any_veto_smoothed:
        ax_any.plot(
            plot_elapsed,
            plot_any_veto_smoothed,
            color="tab:blue",
            linewidth=2.4,
            alpha=0.95,
            label=(
                f"AnyVetoed rate (window-smoothed, n={args.correction_window_seconds:g}s, "
                f"{args.correction_smooth_method})"
            ),
        )
    ax_any.set_ylabel("AnyVetoed Rate [Hz]", fontsize=LABEL_FONT_SIZE)
    ax_any.set_title(f"Run {args.run}: pre-correction AnyVetoed rate", fontsize=TITLE_FONT_SIZE)
    ax_any.set_ylim(bottom=0.0)
    ax_any.grid(True, alpha=0.25)
    ax_any.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    if has_any_veto_raw or has_any_veto_smoothed:
        ax_any.legend(loc="upper right", fontsize=LEGEND_FONT_SIZE)

    if has_pre_ratio_raw:
        ax_any_ratio.plot(
            plot_elapsed,
            plot_pre_ratio_raw,
            color="tab:green",
            linewidth=1.2,
            alpha=0.30,
            label="AnyVetoed/Luminosity (raw, pre-correction)",
        )
    if has_pre_ratio_smoothed:
        ax_any_ratio.plot(
            plot_elapsed,
            plot_pre_ratio_smoothed,
            color="tab:green",
            linewidth=2.2,
            alpha=0.95,
            label="AnyVetoed/Luminosity (smoothed, pre-correction)",
        )
    ax_any_ratio.set_xlabel("Elapsed time from run start [s]", fontsize=LABEL_FONT_SIZE)
    ax_any_ratio.set_ylabel("Rate / Luminosity [events/fb^-1]", fontsize=LABEL_FONT_SIZE)
    ax_any_ratio.set_ylim(bottom=0.0)
    ax_any_ratio.grid(True, alpha=0.25)
    ax_any_ratio.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    if has_pre_ratio_raw or has_pre_ratio_smoothed:
        ax_any_ratio.legend(loc="upper right", fontsize=LEGEND_FONT_SIZE)

    fig_any.tight_layout()
    anyveto_plot_path = output_dir / f"run_{args.run}_anyveto_precorrection.png"
    fig_any.savefig(anyveto_plot_path, dpi=180)
    plt.close(fig_any)

    csv_path = output_dir / f"run_{args.run}_timeseries.csv"
    save_csv(
        csv_path,
        ts_sorted,
        aligned,
        corrected_any_rate,
        corrected_any_rate_bisection,
        ratio_series_events_per_fb,
        ratio_series_bisection_events_per_fb,
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
        "LumiDB": "DcsData",
        "LumiMeasurement": lumi_measurement,
        "LumiField": args.lumi_field,
        "LuminosityUnit": "ub^-1/s",
        "LuminosityInputUnit": "Hz/ub (= ub^-1 s^-1)",
        "LuminosityScaleApplied": "x1",
        "LuminosityQuery": "mean",
        "LumiPlotSmoothBins": args.lumi_plot_smooth_bins,
        "RatioLumiSmoothBins": args.ratio_lumi_smooth_bins,
        "RatioLumiMinFrac": args.ratio_lumi_min_frac,
        "RatioLumiMinAbs": args.ratio_lumi_min_abs,
        "RatioLumiFloor": lumi_ratio_floor,
        "RatioUnit": "events/fb^-1",
        "RatioScaleApplied": "x1e9 from events/ub",
        "MissingTriggerMeasurements": missing_measurements,
        "InvalidDenominatorBins": invalid_bins,
        "InvalidRatioBins": invalid_ratio_bins,
        "NearZeroLuminosityRatioBins": near_zero_lumi_bins,
        "LinearUpperLimitBins": sum(1 for hit in linear_upper_limit_bins if hit),
        "BisectionUpperLimitBins": sum(1 for hit in bisection_upper_limit_bins if hit),
        "PlotStartSeconds": args.plot_start_seconds,
        "PlotEndSeconds": plot_end_seconds,
        "PlottedPoints": len(plot_indices),
        "BisectionVsRateComparisonPlot": str(compare_plot_path),
        "Plots": [str(plot_path), str(corrected_overlay_plot_path), str(compare_plot_path), str(anyveto_plot_path)],
        "CorrectedWithPreCorrectionPlot": str(corrected_overlay_plot_path),
        "PreCorrectionAnyVetoPlot": str(anyveto_plot_path),
        "CSV": str(csv_path),
    }
    json_path = output_dir / f"run_{args.run}_summary.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    print(f"Saved: {plot_path}")
    print(f"Saved: {corrected_overlay_plot_path}")
    print(f"Saved: {compare_plot_path}")
    print(f"Saved: {anyveto_plot_path}")
    print(f"Saved: {csv_path}")
    print(f"Saved: {json_path}")
    print(f"Ratio luminosity floor: {lumi_ratio_floor}")
    if invalid_bins > 0:
        print(f"Warning: {invalid_bins} bins had non-positive denominator in correction and were dropped.")
    if invalid_ratio_bins > 0:
        print(f"Warning: {invalid_ratio_bins} bins were dropped in ratio due to invalid/too-small luminosity denominator.")
    if near_zero_lumi_bins > 0:
        print(f"Warning: {near_zero_lumi_bins} bins were dropped in ratio because raw luminosity < {args.ratio_lumi_min_abs}.")


if __name__ == "__main__":
    main()