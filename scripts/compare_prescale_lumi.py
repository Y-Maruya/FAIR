import argparse
import csv
import math
import sys
from array import array
from dataclasses import dataclass
from pathlib import Path

import ROOT

from Luminosity import discover_lumi_measurement, find_run_time, make_client, query_mean_series, to_iso_utc


@dataclass
class RunLumiSummary:
    run_number: int
    start_s: float
    end_s: float
    n_points: int
    average: float
    rms: float
    std_error: float


@dataclass
class GraphPoint:
    x: float
    y: float
    ex_low: float
    ex_high: float
    ey_low: float
    ey_high: float
    run_number: int


def parse_run_numbers(run_list):
    run_numbers = []
    if run_list is None:
        return run_numbers

    for token in run_list.split(","):
        token = token.strip()
        if not token:
            continue
        run_numbers.append(int(token))
    return run_numbers


def summarize_lumi_series(series):
    values = [value for _, value in series if isinstance(value, (int, float)) and math.isfinite(value)]
    if not values:
        return math.nan, math.nan, math.nan, 0

    average = sum(values) / float(len(values))
    variance = sum((value - average) ** 2 for value in values) / float(len(values))
    rms = math.sqrt(variance)
    std_error = rms / math.sqrt(len(values)) if len(values) > 0 else math.nan
    return average, rms, std_error, len(values)


def build_run_summaries(run_numbers, lumi_client, lumi_measurement, lumi_field, bin_seconds, verbose=False):
    summaries = []
    for run_number in run_numbers:
        start_s, end_s = find_run_time(run_number)
        lumi_series = query_mean_series(
            lumi_client,
            lumi_measurement,
            start_s,
            end_s,
            bin_seconds,
            field=lumi_field,
            verbose=verbose,
        )
        average, rms, std_error, n_points = summarize_lumi_series(lumi_series)
        if not math.isfinite(average):
            raise RuntimeError(f"No valid luminosity values found for run {run_number}")
        summaries.append(
            RunLumiSummary(
                run_number=run_number,
                start_s=start_s,
                end_s=end_s,
                n_points=n_points,
                average=average,
                rms=rms,
                std_error=std_error,
            )
        )
    return summaries


def read_veto_points(run_numbers, rate_dir, lumi_by_run):
    points = []

    for run_number in run_numbers:
        lumi_scale = lumi_by_run.get(run_number)
        if lumi_scale is None or not math.isfinite(lumi_scale) or lumi_scale <= 0.0:
            print(f"Warning: skipping run {run_number} because average luminosity is invalid", file=sys.stderr)
            continue

        file_name = Path(rate_dir) / f"run{run_number}" / "rate_ana.root"
        input_file = ROOT.TFile.Open(str(file_name), "READ")
        if input_file is None or input_file.IsZombie():
            print(f"Warning: cannot open {file_name}", file=sys.stderr)
            continue

        input_graph = input_file.Get("g_veto")
        if input_graph is None:
            print(f"Warning: graph g_veto is missing in {file_name}", file=sys.stderr)
            input_file.Close()
            continue

        for point_index in range(input_graph.GetN()):
            x = array("d", [0.0])
            y = array("d", [0.0])
            input_graph.GetPoint(point_index, x, y)
            x_value = x[0]
            y_value = y[0]
            points.append(
                GraphPoint(
                    x=x_value,
                    y=y_value / lumi_scale,
                    ex_low=input_graph.GetErrorXlow(point_index),
                    ex_high=input_graph.GetErrorXhigh(point_index),
                    ey_low=input_graph.GetErrorYlow(point_index) / lumi_scale,
                    ey_high=input_graph.GetErrorYhigh(point_index) / lumi_scale,
                    run_number=run_number,
                )
            )

        input_file.Close()

    points.sort(key=lambda point: (point.x, point.run_number))
    return points


def build_combined_graph(points, graph_name, graph_title):
    graph = ROOT.TGraphAsymmErrors()
    graph.SetName(graph_name)
    graph.SetTitle(graph_title)

    for index, point in enumerate(points):
        graph.SetPoint(index, point.x, point.y)
        graph.SetPointError(index, point.ex_low, point.ex_high, point.ey_low, point.ey_high)

    return graph


def set_graph_style(graph, color, marker_style):
    if graph is None:
        return
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    graph.SetMarkerStyle(marker_style)
    graph.SetMarkerSize(1.2)
    graph.SetLineWidth(2)


def estimate_amplitude(graph):
    if graph is None or graph.GetN() == 0:
        return 1.0

    sum_ax = 0.0
    valid_points = 0
    for point_index in range(graph.GetN()):
        x = array("d", [0.0])
        y = array("d", [0.0])
        graph.GetPoint(point_index, x, y)
        x_value = x[0]
        y_value = y[0]
        if x_value <= 0.0:
            continue
        sum_ax += x_value * y_value
        valid_points += 1

    if valid_points == 0:
        return 1.0
    return sum_ax / float(valid_points)


def get_average_y_error(graph, point_index):
    if graph is None:
        return 0.0
    return 0.5 * (graph.GetErrorYlow(point_index) + graph.GetErrorYhigh(point_index))


def build_relative_residuals(graph, fit_function, graph_name, graph_title):
    residuals = ROOT.TGraphErrors()
    residuals.SetName(graph_name)
    residuals.SetTitle(graph_title)

    if graph is None or fit_function is None:
        return residuals

    output_index = 0
    for point_index in range(graph.GetN()):
        x = array("d", [0.0])
        y = array("d", [0.0])
        graph.GetPoint(point_index, x, y)
        x_value = x[0]
        y_value = y[0]
        fit_value = fit_function.Eval(x_value)
        if fit_value == 0.0:
            continue

        residual_percent = 100.0 * (y_value - fit_value) / fit_value
        residual_error = 100.0 * get_average_y_error(graph, point_index) / abs(fit_value)
        x_error = 0.5 * (graph.GetErrorXlow(point_index) + graph.GetErrorXhigh(point_index))
        residuals.SetPoint(output_index, x_value, residual_percent)
        residuals.SetPointError(output_index, x_error, residual_error)
        output_index += 1

    return residuals


def build_run_lumi_histogram(summaries):
    histogram = ROOT.TH1D(
        "h_run_avg_lumi",
        "Run average luminosity;Run number;Average luminosity [ub^{-1} s^{-1}]",
        len(summaries),
        0.5,
        len(summaries) + 0.5,
    )
    histogram.SetDirectory(0)

    for index, summary in enumerate(summaries, start=1):
        histogram.GetXaxis().SetBinLabel(index, str(summary.run_number))
        histogram.SetBinContent(index, summary.average)
        histogram.SetBinError(index, summary.std_error if math.isfinite(summary.std_error) else 0.0)

    return histogram


def write_run_summary_csv(path, summaries):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "run",
            "start_utc",
            "end_utc",
            "n_lumi_points",
            "avg_lumi_ub_inv_s",
            "rms_lumi_ub_inv_s",
            "stderr_lumi_ub_inv_s",
        ])
        for summary in summaries:
            writer.writerow([
                summary.run_number,
                to_iso_utc(summary.start_s),
                to_iso_utc(summary.end_s),
                summary.n_points,
                f"{summary.average:.12g}",
                f"{summary.rms:.12g}" if math.isfinite(summary.rms) else "",
                f"{summary.std_error:.12g}" if math.isfinite(summary.std_error) else "",
            ])


def get_x_range(graph):
    x_min = math.inf
    x_max = -math.inf
    for point_index in range(graph.GetN()):
        x = array("d", [0.0])
        y = array("d", [0.0])
        graph.GetPoint(point_index, x, y)
        x_value = x[0]
        y_value = y[0]
        if x_value <= 0.0 or y_value <= 0.0:
            continue
        x_min = min(x_min, x_value)
        x_max = max(x_max, x_value)

    if not math.isfinite(x_min) or not math.isfinite(x_max):
        raise RuntimeError("No valid positive x values found in normalized g_veto")
    return x_min, x_max


def get_y_range(graph, fit_functions):
    y_min = math.inf
    y_max = -math.inf

    for point_index in range(graph.GetN()):
        x = array("d", [0.0])
        y = array("d", [0.0])
        graph.GetPoint(point_index, x, y)
        y_value = y[0]
        if not math.isfinite(y_value) or y_value <= 0.0:
            continue
        y_min = min(y_min, y_value)
        y_max = max(y_max, y_value)

    for fit_function, x_values in fit_functions:
        if fit_function is None:
            continue
        for x_value in x_values:
            fit_value = fit_function.Eval(x_value)
            if math.isfinite(fit_value) and fit_value > 0.0:
                y_min = min(y_min, fit_value)
                y_max = max(y_max, fit_value)

    if not math.isfinite(y_min) or y_min <= 0.0:
        y_min = 1e-6
    if not math.isfinite(y_max) or y_max <= y_min:
        y_max = y_min * 10.0

    return y_min, y_max


def residual_range(graph):
    y_min = math.inf
    y_max = -math.inf
    for point_index in range(graph.GetN()):
        x = array("d", [0.0])
        y = array("d", [0.0])
        graph.GetPoint(point_index, x, y)
        y_value = y[0]
        if not math.isfinite(y_value):
            continue
        y_min = min(y_min, y_value)
        y_max = max(y_max, y_value)
    if not math.isfinite(y_min) or not math.isfinite(y_max):
        return -10.0, 10.0
    if y_min == y_max:
        return y_min - 1.0, y_max + 1.0
    pad = 0.25 * max(abs(y_min), abs(y_max))
    return y_min - pad, y_max + pad


def main():
    parser = argparse.ArgumentParser(
        description="Fetch run-average luminosity, normalize g_veto by it, and fit the normalized veto rate."
    )
    parser.add_argument(
        "--run-list",
        default="22711,22712,22713,22714,22715,22716,22717,22718",
        help="Comma-separated AHCAL run numbers",
    )
    parser.add_argument(
        "--rate-dir",
        default="/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out/RateAna",
        help="Directory containing run*/rate_ana.root files",
    )
    parser.add_argument("--bin-seconds", type=float, default=0.2, help="Luminosity InfluxDB bin width in seconds")
    parser.add_argument("--lumi-measurement", default=None, help="Luminosity measurement name in DcsData")
    parser.add_argument("--lumi-field", default="Luminosity", help="Field name of the luminosity measurement")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    parser.add_argument("--output-root", default="combined_rate_ana_lumi_norm.root", help="Output ROOT file name")
    parser.add_argument("--output-prefix", default="compare_prescale_lumi", help="Output file prefix for images")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose debug output")
    args = parser.parse_args()

    if args.bin_seconds <= 0.0:
        raise RuntimeError("--bin-seconds must be > 0")

    run_numbers = parse_run_numbers(args.run_list)
    if not run_numbers:
        raise RuntimeError("No valid run numbers were provided")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    lumi_client = make_client("DcsData")
    lumi_measurement = args.lumi_measurement or discover_lumi_measurement(lumi_client, verbose=args.verbose)

    summaries = build_run_summaries(
        run_numbers,
        lumi_client,
        lumi_measurement,
        args.lumi_field,
        args.bin_seconds,
        verbose=args.verbose,
    )
    lumi_by_run = {summary.run_number: summary.average for summary in summaries}

    print("Run average luminosity summary:")
    for summary in summaries:
        print(
            f"  run {summary.run_number}: avg={summary.average:.6g} ub^-1 s^-1, "
            f"rms={summary.rms:.6g}, n={summary.n_points}"
        )

    points = read_veto_points(run_numbers, args.rate_dir, lumi_by_run)
    if not points:
        raise RuntimeError("No veto points could be loaded")

    normalized_graph = build_combined_graph(
        points,
        "g_veto_lumi_norm",
        "Veto Trigger Rate / <Luminosity>;CalibRate;Rate / <Luminosity>",
    )
    set_graph_style(normalized_graph, ROOT.kRed + 1, 20)

    x_min, x_max = get_x_range(normalized_graph)

    veto_fit = ROOT.TF1("f_veto_prescale_lumi", "[0]/(2*x+[1])", x_min, x_max)
    veto_fit.SetParName(0, "A")
    veto_fit.SetParameter(0, estimate_amplitude(normalized_graph))
    veto_fit.SetParName(1, "B")
    veto_fit.SetParameter(1, 1.0)
    veto_fit.SetLineColor(ROOT.kRed + 2)
    veto_fit.SetLineWidth(2)

    veto_fit_a2x = ROOT.TF1("f_veto_prescale_lumi_a2x", "[0]/(2*x)", x_min, x_max)
    veto_fit_a2x.SetParName(0, "A")
    veto_fit_a2x.SetParameter(0, estimate_amplitude(normalized_graph))
    veto_fit_a2x.SetLineColor(ROOT.kBlue + 2)
    veto_fit_a2x.SetLineStyle(2)
    veto_fit_a2x.SetLineWidth(2)

    fit_result = normalized_graph.Fit(veto_fit, "Q0S")
    fit_result_a2x = normalized_graph.Fit(veto_fit_a2x, "Q0S+")

    residuals = build_relative_residuals(
        normalized_graph,
        veto_fit,
        "g_veto_lumi_norm_residuals",
        "Veto Trigger Rate / <Luminosity> residuals;CalibRate; (data-fit)/fit [%]",
    )
    residuals_a2x = build_relative_residuals(
        normalized_graph,
        veto_fit_a2x,
        "g_veto_lumi_norm_residuals_a2x",
        "Veto Trigger Rate / <Luminosity> residuals (A/(2*x));CalibRate; (data-fit)/fit [%]",
    )

    if fit_result.Get() is not None:
        print(
            f"A/(2*x+B): A={veto_fit.GetParameter(0):.6g}, B={veto_fit.GetParameter(1):.6g}, "
            f"chi2/ndf={fit_result.Chi2():.3f}/{fit_result.Ndf()}={fit_result.Chi2() / (fit_result.Ndf() if fit_result.Ndf() > 0 else 1):.3f}, "
            f"prob={fit_result.Prob():.3g}"
        )
    if fit_result_a2x.Get() is not None:
        print(
            f"A/(2*x): A={veto_fit_a2x.GetParameter(0):.6g}, "
            f"chi2/ndf={fit_result_a2x.Chi2():.3f}/{fit_result_a2x.Ndf()}={fit_result_a2x.Chi2() / (fit_result_a2x.Ndf() if fit_result_a2x.Ndf() > 0 else 1):.3f}, "
            f"prob={fit_result_a2x.Prob():.3g}"
        )

    lumi_hist = build_run_lumi_histogram(summaries)

    y_min, y_max = get_y_range(normalized_graph, [(veto_fit, [x_min, x_max]), (veto_fit_a2x, [x_min, x_max])])
    y_min = max(y_min * 0.5, 1e-9)
    y_max = y_max * 2.0

    plot_canvas = ROOT.TCanvas("c_compare_prescale_lumi", "Normalized trigger rate", 900, 700)
    plot_canvas.SetLogx()
    plot_canvas.SetLogy()
    normalized_graph.SetMinimum(y_min)
    normalized_graph.SetMaximum(y_max)
    normalized_graph.Draw("AP")
    veto_fit.Draw("SAME")
    veto_fit_a2x.Draw("SAME")
    plot_legend = ROOT.TLegend(0.58, 0.70, 0.88, 0.88)
    plot_legend.SetBorderSize(0)
    plot_legend.SetFillStyle(0)
    plot_legend.AddEntry(normalized_graph, "g_veto / <L>", "lp")
    plot_legend.AddEntry(
        veto_fit,
        f"A/(2*x+B): A={veto_fit.GetParameter(0):.3g}, B={veto_fit.GetParameter(1):.3g}",
        "l",
    )
    plot_legend.AddEntry(veto_fit_a2x, f"A/(2*x): A={veto_fit_a2x.GetParameter(0):.3g}", "l")
    plot_legend.Draw()
    plot_latex = ROOT.TLatex()
    plot_latex.SetNDC()
    plot_latex.SetTextSize(0.035)
    if fit_result.Get() is not None:
        plot_latex.DrawLatex(
            0.15,
            0.84,
            f"A/(2*x+B): #chi^{{2}}/ndf = {fit_result.Chi2():.2f}/{fit_result.Ndf()} = {fit_result.Chi2() / (fit_result.Ndf() if fit_result.Ndf() > 0 else 1):.2f}",
        )
    if fit_result_a2x.Get() is not None:
        plot_latex.DrawLatex(
            0.15,
            0.78,
            f"A/(2*x): #chi^{{2}}/ndf = {fit_result_a2x.Chi2():.2f}/{fit_result_a2x.Ndf()} = {fit_result_a2x.Chi2() / (fit_result_a2x.Ndf() if fit_result_a2x.Ndf() > 0 else 1):.2f}",
        )

    fit_canvas = ROOT.TCanvas("c_veto_prescale_lumi_fit", "Veto trigger rate fit", 900, 900)
    top_pad = ROOT.TPad("top_pad", "top_pad", 0.0, 0.30, 1.0, 1.0)
    bottom_pad = ROOT.TPad("bottom_pad", "bottom_pad", 0.0, 0.0, 1.0, 0.30)
    top_pad.SetBottomMargin(0.02)
    bottom_pad.SetTopMargin(0.05)
    bottom_pad.SetBottomMargin(0.25)
    top_pad.Draw()
    bottom_pad.Draw()
    top_pad.SetLogx()
    top_pad.SetLogy()
    bottom_pad.SetLogx()

    top_pad.cd()
    normalized_graph.Draw("AP")
    veto_fit.Draw("SAME")
    veto_fit_a2x.Draw("SAME")
    plot_legend.Draw()
    fit_latex = ROOT.TLatex()
    fit_latex.SetNDC()
    fit_latex.SetTextSize(0.035)
    if fit_result.Get() is not None:
        fit_latex.DrawLatex(
            0.15,
            0.84,
            f"A/(2*x+B): #chi^{{2}}/ndf = {fit_result.Chi2():.2f}/{fit_result.Ndf()} = {fit_result.Chi2() / (fit_result.Ndf() if fit_result.Ndf() > 0 else 1):.2f}",
        )
    if fit_result_a2x.Get() is not None:
        fit_latex.DrawLatex(
            0.15,
            0.78,
            f"A/(2*x): #chi^{{2}}/ndf = {fit_result_a2x.Chi2():.2f}/{fit_result_a2x.Ndf()} = {fit_result_a2x.Chi2() / (fit_result_a2x.Ndf() if fit_result_a2x.Ndf() > 0 else 1):.2f}",
        )

    bottom_pad.cd()
    residuals.SetTitle(";CalibRate; (data-fit)/fit [%]")
    residuals.SetMarkerStyle(20)
    residuals.SetMarkerSize(1.0)
    residuals.SetMarkerColor(ROOT.kRed + 1)
    residuals.SetLineColor(ROOT.kRed + 1)
    residuals_a2x.SetMarkerStyle(21)
    residuals_a2x.SetMarkerSize(1.0)
    residuals_a2x.SetMarkerColor(ROOT.kBlue + 1)
    residuals_a2x.SetLineColor(ROOT.kBlue + 1)
    res_min, res_max = residual_range(residuals)
    res2_min, res2_max = residual_range(residuals_a2x)
    residuals.SetMinimum(min(res_min, res2_min))
    residuals.SetMaximum(max(res_max, res2_max))
    residuals.Draw("AP")
    residuals_a2x.Draw("P SAME")
    zero_line = ROOT.TLine(x_min, 0.0, x_max, 0.0)
    zero_line.SetLineStyle(2)
    zero_line.SetLineColor(ROOT.kGray + 2)
    zero_line.Draw("SAME")
    residual_legend = ROOT.TLegend(0.14, 0.70, 0.42, 0.88)
    residual_legend.SetBorderSize(0)
    residual_legend.SetFillStyle(0)
    residual_legend.AddEntry(residuals, "Residuals: A/(2*x+B)", "lp")
    residual_legend.AddEntry(residuals_a2x, "Residuals: A/(2*x)", "lp")
    residual_legend.Draw()

    lumi_canvas = ROOT.TCanvas("c_run_average_lumi", "Run average luminosity", 900, 600)
    lumi_hist.SetLineColor(ROOT.kBlue + 1)
    lumi_hist.SetMarkerColor(ROOT.kBlue + 1)
    lumi_hist.SetMarkerStyle(20)
    lumi_hist.Draw("E1")

    plot_canvas.Modified()
    plot_canvas.Update()
    fit_canvas.Modified()
    fit_canvas.Update()
    lumi_canvas.Modified()
    lumi_canvas.Update()

    output_root = output_dir / args.output_root
    root_file = ROOT.TFile.Open(str(output_root), "RECREATE")
    lumi_hist.Write()
    normalized_graph.Write()
    veto_fit.Write()
    veto_fit_a2x.Write()
    residuals.Write()
    residuals_a2x.Write()
    plot_canvas.Write()
    fit_canvas.Write()
    lumi_canvas.Write()
    root_file.Close()

    plot_canvas.SaveAs(str(output_dir / f"{args.output_prefix}.png"))
    plot_canvas.SaveAs(str(output_dir / f"{args.output_prefix}.pdf"))
    fit_canvas.SaveAs(str(output_dir / f"{args.output_prefix}_fit.png"))
    fit_canvas.SaveAs(str(output_dir / f"{args.output_prefix}_fit.pdf"))
    lumi_canvas.SaveAs(str(output_dir / f"{args.output_prefix}_run_average_lumi.png"))
    lumi_canvas.SaveAs(str(output_dir / f"{args.output_prefix}_run_average_lumi.pdf"))

    write_run_summary_csv(output_dir / f"{args.output_prefix}_run_summary.csv", summaries)
    print(f"Wrote {output_root}")


if __name__ == "__main__":
    main()
