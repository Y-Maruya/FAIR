#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLine.h"
#include "TSystem.h"
#include "TTree.h"

// ============================================================================
// Data Structures
// ============================================================================

struct InterCalibEntry {
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int ch = -1;
  int channel_index = -1;
  long long n_points = 0;
  double intercept = -999.0;  // p0
  double slope = -999.0;       // p1
  double chi2 = -1.0;
  int ndf = -1;
  double chi2_ndf = -1.0;
  int fit_status = 999;
  int fit_ok = 0;
  double x_mm = -999.0;
  double y_mm = -999.0;
};

struct RunData {
  int run = 0;
  std::string label;
  std::string path;
  std::vector<InterCalibEntry> entries;
  std::map<int, InterCalibEntry> byCellId;
  bool hasFitFlags = false;
  bool hasCoordinates = false;
};

struct SummaryStats {
  double mean = 0.0;
  double rms = 0.0;
  double min = 0.0;
  double max = 0.0;
  int n = 0;
};

struct CellShiftInfo {
  int cellid = -1;
  double spread = 0.0;
  double refValue = 0.0;
  double meanValue = 0.0;
  double x_mm = 0.0;
  double y_mm = 0.0;
};

// ============================================================================
// Utility Functions
// ============================================================================

std::vector<int> parseRunList(const char *csv) {
  std::vector<int> runs;
  if (!csv) return runs;

  if (std::string(csv) == std::string("all")) {
    gSystem->Exec("ls -d */ | sed 's#/##' > run_list.txt");
    std::ifstream ifs("run_list.txt");
    std::string line;
    while (std::getline(ifs, line)) {
      line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
      if (line.empty()) continue;
      runs.push_back(std::atoi(line.c_str()));
    }
    return runs;
  }

  std::stringstream ss(csv);
  std::string token;
  while (std::getline(ss, token, ',')) {
    token.erase(std::remove_if(token.begin(), token.end(), ::isspace), token.end());
    if (token.empty()) continue;
    runs.push_back(std::atoi(token.c_str()));
  }
  return runs;
}

bool setBranchIfExists(TTree *tree, const char *name, void *addr) {
  if (!tree || !tree->GetBranch(name)) return false;
  tree->SetBranchAddress(name, addr);
  return true;
}

bool loadRunFromFile(const std::string &filePath,
                     int run,
                     const std::string &label,
                     RunData &out) {
  out = RunData();
  out.run = run;
  out.label = label;
  out.path = filePath;

  TFile *file = TFile::Open(out.path.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Failed to open file: " << out.path << std::endl;
    return false;
  }

  TTree *tree = dynamic_cast<TTree *>(file->Get("intercalib"));
  if (!tree) {
    std::cerr << "Missing 'intercalib' tree in: " << out.path << std::endl;
    file->Close();
    return false;
  }

  InterCalibEntry e;
  if (!setBranchIfExists(tree, "cellid", &e.cellid) ||
      !setBranchIfExists(tree, "layer", &e.layer) ||
      !setBranchIfExists(tree, "chip", &e.chip) ||
      !setBranchIfExists(tree, "ch", &e.ch) ||
      !setBranchIfExists(tree, "n_points", &e.n_points) ||
      !setBranchIfExists(tree, "intercept", &e.intercept) ||
      !setBranchIfExists(tree, "slope", &e.slope) ||
      !setBranchIfExists(tree, "chi2_ndf", &e.chi2_ndf)) {
    std::cerr << "Missing required branches in intercalib tree: " << out.path << std::endl;
    file->Close();
    return false;
  }

  // Optional branches
  setBranchIfExists(tree, "channel_index", &e.channel_index);
  setBranchIfExists(tree, "chi2", &e.chi2);
  setBranchIfExists(tree, "ndf", &e.ndf);
  setBranchIfExists(tree, "fit_status", &e.fit_status);
  out.hasFitFlags = setBranchIfExists(tree, "fit_ok", &e.fit_ok);
  bool hasX = setBranchIfExists(tree, "x_mm", &e.x_mm);
  bool hasY = setBranchIfExists(tree, "y_mm", &e.y_mm);
  out.hasCoordinates = hasX && hasY;

  const Long64_t nEntries = tree->GetEntries();
  out.entries.reserve(nEntries);

  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    out.entries.push_back(e);
    out.byCellId[e.cellid] = e;
  }

  if (!out.hasFitFlags) {
    std::cout << "Warning: fit_ok branch not found" << std::endl;
  }
  if (!out.hasCoordinates) {
    std::cout << "Warning: x_mm/y_mm branches not found" << std::endl;
  }

  file->Close();
  return true;
}

bool loadRun(const std::string &baseDir, int run, RunData &out) {
  const std::string path = baseDir + "/" + std::to_string(run) + "/intercalib.root";
  return loadRunFromFile(path, run, std::to_string(run), out);
}

void setGraphStyle(TGraph *g, int color, int markerStyle) {
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(markerStyle);
  g->SetMarkerSize(1.2);
  g->SetLineWidth(2);
}

bool getBranchValue(const InterCalibEntry &e, const std::string &branch, bool onlyFitOk, double &value) {
  if (onlyFitOk && !e.fit_ok) return false;

  if (branch == "slope") { value = e.slope; return true; }
  if (branch == "intercept") { value = e.intercept; return true; }
  if (branch == "chi2_ndf") { value = e.chi2_ndf; return true; }
  if (branch == "n_points") { value = (double)e.n_points; return true; }
  if (branch == "fit_ok") { value = e.fit_ok ? 1.0 : 0.0; return true; }

  return false;
}

SummaryStats computeStats(const std::vector<double> &values) {
  SummaryStats s;
  s.n = static_cast<int>(values.size());
  if (values.empty()) return s;

  s.min = *std::min_element(values.begin(), values.end());
  s.max = *std::max_element(values.begin(), values.end());
  s.mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();

  double sum2 = 0.0;
  for (double v : values) {
    sum2 += (v - s.mean) * (v - s.mean);
  }
  s.rms = std::sqrt(sum2 / values.size());
  return s;
}

SummaryStats extractStats(const RunData &runData, const std::string &what, bool onlyFitOk) {
  std::vector<double> values;
  values.reserve(runData.entries.size());

  for (const auto &e : runData.entries) {
    double val;
    if (getBranchValue(e, what, onlyFitOk, val)) {
      values.push_back(val);
    }
  }

  return computeStats(values);
}

std::vector<CellShiftInfo> rankCellsBySpread(const std::vector<RunData> &runs,
                                             const std::string &branch,
                                             bool onlyFitOk) {
  std::vector<CellShiftInfo> ranked;
  if (runs.size() < 2) return ranked;

  const RunData &ref = runs.front();
  for (const auto &kv : ref.byCellId) {
    const int cellid = kv.first;
    const InterCalibEntry &ref_e = kv.second;

    double ref_val;
    if (!getBranchValue(ref_e, branch, onlyFitOk, ref_val)) continue;

    std::vector<double> values;
    values.push_back(ref_val);

    for (size_t i = 1; i < runs.size(); ++i) {
      auto it = runs[i].byCellId.find(cellid);
      if (it == runs[i].byCellId.end()) continue;

      double val;
      if (!getBranchValue(it->second, branch, onlyFitOk, val)) continue;

      values.push_back(val);
    }

    if (values.size() < 2) continue;

    SummaryStats stats = computeStats(values);
    CellShiftInfo info;
    info.cellid = cellid;
    info.spread = stats.rms;
    info.refValue = ref_val;
    info.meanValue = stats.mean;
    info.x_mm = ref_e.x_mm;
    info.y_mm = ref_e.y_mm;
    ranked.push_back(info);
  }

  std::sort(ranked.begin(), ranked.end(),
            [](const CellShiftInfo &a, const CellShiftInfo &b) {
              if (a.spread != b.spread) return a.spread > b.spread;
              return a.cellid < b.cellid;
            });

  return ranked;
}

// ============================================================================
// Visualization Functions
// ============================================================================

TH1D *makeDeltaHist(const char *name, const char *title,
                    const std::vector<RunData> &runs,
                    const std::string &branch, bool onlyFitOk) {
  if (runs.size() < 2) return nullptr;

  std::vector<double> deltas;
  const RunData &ref = runs.front();

  for (size_t i = 1; i < runs.size(); ++i) {
    for (const auto &kv : ref.byCellId) {
      const int cellid = kv.first;
      const InterCalibEntry &ref_e = kv.second;

      double ref_val;
      if (!getBranchValue(ref_e, branch, onlyFitOk, ref_val)) continue;

      auto it = runs[i].byCellId.find(cellid);
      if (it == runs[i].byCellId.end()) continue;

      double val;
      if (!getBranchValue(it->second, branch, onlyFitOk, val)) continue;

      deltas.push_back(val - ref_val);
    }
  }

  if (deltas.empty()) return nullptr;

  auto mm = std::minmax_element(deltas.begin(), deltas.end());
  double xmin = *mm.first;
  double xmax = *mm.second;
  if (xmin == xmax) {
    xmin -= 0.1;
    xmax += 0.1;
  }

  TH1D *h = new TH1D(name, title, 80, xmin - 0.1 * std::abs(xmin), xmax + 0.1 * std::abs(xmax));
  for (double d : deltas) h->Fill(d);
  h->SetLineWidth(2);
  h->SetStats(1);
  return h;
}

std::vector<double> collectDeltaValues(const RunData &ref,
                                       const RunData &target,
                                       const std::string &branch, bool onlyFitOk) {
  std::vector<double> deltas;
  deltas.reserve(ref.byCellId.size());

  for (const auto &kv : ref.byCellId) {
    const int cellid = kv.first;
    const InterCalibEntry &ref_e = kv.second;

    double ref_val;
    if (!getBranchValue(ref_e, branch, onlyFitOk, ref_val)) continue;

    auto it = target.byCellId.find(cellid);
    if (it == target.byCellId.end()) continue;

    double val;
    if (!getBranchValue(it->second, branch, onlyFitOk, val)) continue;

    deltas.push_back(val - ref_val);
  }

  return deltas;
}

void drawDeltaFromReferenceByRun(const std::vector<RunData> &runs,
                                 const std::string &branch,
                                 const std::string &branchLabel,
                                 const std::string &outTag,
                                 const std::string &outDir,
                                 bool onlyFitOk = true) {
  if (runs.size() < 2) return;

  const RunData &ref = runs.front();
  std::vector<std::vector<double>> allDeltas(runs.size());
  double globalMin = 1e30;
  double globalMax = -1e30;

  for (size_t i = 1; i < runs.size(); ++i) {
    allDeltas[i] = collectDeltaValues(ref, runs[i], branch, onlyFitOk);
    if (!allDeltas[i].empty()) {
      double mn = *std::min_element(allDeltas[i].begin(), allDeltas[i].end());
      double mx = *std::max_element(allDeltas[i].begin(), allDeltas[i].end());
      globalMin = std::min(globalMin, mn);
      globalMax = std::max(globalMax, mx);
    }
  }

  if (globalMin > globalMax) return;
  if (globalMin == globalMax) {
    globalMin -= 0.1;
    globalMax += 0.1;
  }

  const double span = std::max(1e-6, globalMax - globalMin);
  const double xmin = globalMin - 0.1 * span;
  const double xmax = globalMax + 0.1 * span;
  const int nbins = 80;

  TCanvas *c = new TCanvas(Form("c_%s_%s", branch.c_str(), outTag.c_str()),
                           Form("Delta %s vs run", branchLabel.c_str()),
                           1200, 600);
  c->Divide(2, 1);

  // Left: overlaid histograms
  c->cd(1);
  TH1D *hFirst = nullptr;
  int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2};
  for (size_t i = 1; i < runs.size(); ++i) {
    TH1D *h = new TH1D(Form("hDelta_%zu", i), Form("%s vs run %d", branchLabel.c_str(), runs[i].run),
                       nbins, xmin, xmax);
    for (double d : allDeltas[i]) h->Fill(d);
    h->SetLineColor(colors[(i - 1) % 6]);
    h->SetLineWidth(2);
    if (i == 1) {
      hFirst = h;
      hFirst->Draw();
    } else {
      h->Draw("same");
    }
  }
  gPad->SetLogy(0);

  // Right: summary bar chart
  c->cd(2);
  std::vector<double> means;
  std::vector<double> rmss;
  std::vector<std::string> labels;

  means.push_back(0.0);
  rmss.push_back(0.0);
  labels.push_back(ref.label);

  for (size_t i = 1; i < runs.size(); ++i) {
    SummaryStats stats = computeStats(allDeltas[i]);
    means.push_back(stats.mean);
    rmss.push_back(stats.rms);
    labels.push_back(runs[i].label);
  }

  TH1D *hBar = new TH1D("hBar", Form("Mean #pm RMS of %s deltas", branchLabel.c_str()),
                        runs.size(), 0, runs.size());
  for (size_t i = 0; i < means.size(); ++i) {
    hBar->SetBinContent(i + 1, means[i]);
    hBar->SetBinError(i + 1, rmss[i]);
    hBar->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
  }
  hBar->SetLineWidth(2);
  hBar->SetMarkerStyle(20);
  hBar->SetMarkerSize(1.0);
  hBar->Draw("E");
  gPad->SetBottomMargin(0.3);
  gPad->SetGridy();

  c->SaveAs((outDir + "/Delta_" + branch + "_" + outTag + ".pdf").c_str());
}

void drawRepresentativeTrendSet(const std::vector<RunData> &runs,
                                const std::vector<CellShiftInfo> &cells,
                                const std::string &branch,
                                const std::string &branchLabel,
                                const std::string &label,
                                const std::string &outDir,
                                const std::string &canvasTag,
                                bool onlyFitOk) {
  if (cells.empty()) return;

  const int n = static_cast<int>(cells.size());
  const int nCols = 2;
  const int nRows = (n + nCols - 1) / nCols;
  TCanvas *c = new TCanvas(Form("c_%s_%s", branch.c_str(), canvasTag.c_str()),
                           Form("%s %s", branchLabel.c_str(), label.c_str()),
                           1400, 420 * nRows);
  c->Divide(nCols, nRows);

  for (int i = 0; i < n; ++i) {
    c->cd(i + 1);

    std::vector<double> runIndices, values;
    int colors[] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7};

    int runIndex = 0;
    for (const auto &run : runs) {
      auto it = run.byCellId.find(cells[i].cellid);
      if (it == run.byCellId.end()) continue;

      double val;
      if (!getBranchValue(it->second, branch, onlyFitOk, val)) continue;

      runIndices.push_back(runIndex);
      values.push_back(val);
      runIndex++;
    }

    if (!runIndices.empty()) {
      TGraph *g = new TGraph(runIndices.size(), runIndices.data(), values.data());
      g->SetMarkerColor(kBlack);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(1.2);
      g->SetLineColor(kBlack);
      g->SetLineWidth(2);
      g->SetTitle(Form("CellID %d", cells[i].cellid));
      g->GetXaxis()->SetLimits(-0.5, runIndex - 0.5);
      g->GetXaxis()->SetTitle("Run index");
      g->GetYaxis()->SetTitle(branchLabel.c_str());
      g->Draw("ALP");
      gPad->SetGridy();

      TLatex tl;
      tl.SetNDC(true);
      tl.SetTextSize(0.08);
      tl.DrawLatex(0.15, 0.92, Form("CellID %d (L%d C%d Ch%d)", cells[i].cellid,
                                     (cells[i].cellid / 100000),
                                     ((cells[i].cellid / 10000) % 10),
                                     (cells[i].cellid % 10000)));
    }
  }

  c->SaveAs((outDir + "/" + branch + "_" + canvasTag + ".pdf").c_str());
}

void drawSummaryGraphs(const std::vector<RunData> &runs, const std::string &outDir) {
  std::vector<std::string> branches = {"slope", "intercept", "chi2_ndf", "n_points"};
  std::vector<std::string> labels = {"Slope (p1)", "Intercept (p0)", "Chi2/ndf", "N Points"};

  for (size_t b = 0; b < branches.size(); ++b) {
    std::vector<double> runNumbers;
    std::vector<double> means;
    std::vector<double> rmss;

    for (size_t i = 0; i < runs.size(); ++i) {
      SummaryStats stats = extractStats(runs[i], branches[b], true);
      if (stats.n == 0) continue;

      runNumbers.push_back(i);
      means.push_back(stats.mean);
      rmss.push_back(stats.rms);
    }

    if (means.empty()) continue;

    TCanvas *c = new TCanvas(Form("c_summary_%s", branches[b].c_str()),
                             Form("Summary: %s", labels[b].c_str()), 800, 600);

    TGraph *gMean = new TGraph(means.size(), runNumbers.data(), means.data());
    gMean->SetTitle(Form("%s vs run", labels[b].c_str()));
    gMean->SetMarkerColor(kBlack);
    gMean->SetMarkerStyle(20);
    gMean->SetMarkerSize(1.2);
    gMean->SetLineColor(kBlack);
    gMean->SetLineWidth(2);
    gMean->GetXaxis()->SetTitle("Run index");
    gMean->GetYaxis()->SetTitle(labels[b].c_str());
    gMean->Draw("ALP");

    c->SaveAs((outDir + "/Summary_" + branches[b] + ".pdf").c_str());
  }
}

void printSummaryTable(const std::vector<RunData> &runs) {
  std::vector<std::string> branches = {"slope", "intercept", "chi2_ndf", "n_points"};
  std::vector<std::string> labels = {"Slope (p1)", "Intercept (p0)", "Chi2/ndf", "N Points"};

  std::cout << "\n";
  std::cout << "================================================================================\n";
  std::cout << "INTERCALIB COMPARISON SUMMARY\n";
  std::cout << "================================================================================\n\n";

  for (size_t b = 0; b < branches.size(); ++b) {
    std::cout << labels[b] << ":\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << "Run          ";
    std::cout << "  Mean (±RMS)       ";
    std::cout << "  Min        Max        Entries\n";
    std::cout << std::string(80, '-') << "\n";

    for (const auto &run : runs) {
      SummaryStats stats = extractStats(run, branches[b], false);
      printf("%-12s  %8.4f ± %-8.4f  %8.4f  %8.4f  %7d\n",
             run.label.c_str(), stats.mean, stats.rms, stats.min, stats.max, stats.n);
    }
    std::cout << "\n";
  }

  // Fit quality summary
  std::cout << "Fit Quality (fit_ok):\n";
  std::cout << std::string(80, '-') << "\n";
  std::cout << "Run          ";
  std::cout << "  fit_ok=1   fit_ok=0   Total      Fraction OK\n";
  std::cout << std::string(80, '-') << "\n";

  for (const auto &run : runs) {
    int ok_count = 0;
    int total = run.entries.size();
    for (const auto &e : run.entries) {
      if (e.fit_ok) ok_count++;
    }
    printf("%-12s  %8d   %8d   %8d   %.2f%%\n",
           run.label.c_str(), ok_count, total - ok_count, total,
           100.0 * ok_count / (total > 0 ? total : 1));
  }

  std::cout << "\n================================================================================\n\n";
}

void writeRepresentativeSummary(const std::string &fileName,
                                const std::vector<CellShiftInfo> &shifted,
                                const std::vector<CellShiftInfo> &stable,
                                const std::string &branch) {
  std::ofstream ofs(fileName.c_str());
  if (!ofs) return;

  ofs << "# branch = " << branch << "\n";
  ofs << "# shifted channels (top 12)\n";
  ofs << "# cellid spread refValue meanValue x_mm y_mm\n";
  for (size_t i = 0; i < std::min(size_t(12), shifted.size()); ++i) {
    const auto &info = shifted[i];
    ofs << info.cellid << " " << info.spread << " " << info.refValue << " "
        << info.meanValue << " " << info.x_mm << " " << info.y_mm << "\n";
  }

  ofs << "\n# stable channels (bottom 12)\n";
  ofs << "# cellid spread refValue meanValue x_mm y_mm\n";
  int start = std::max(int(stable.size()) - 12, 0);
  for (int i = start; i < (int)stable.size(); ++i) {
    const auto &info = stable[i];
    ofs << info.cellid << " " << info.spread << " " << info.refValue << " "
        << info.meanValue << " " << info.x_mm << " " << info.y_mm << "\n";
  }
}

// ============================================================================
// Main Comparison Function
// ============================================================================

void compare(const char *baseDir = ".",
             const char *runList = "22163,22296",
             const char *outDir = "comparison_out") {
  // Create output directory
  gSystem->Exec(Form("mkdir -p %s", outDir));

  // Parse run list
  std::vector<int> runNumbers = parseRunList(runList);
  if (runNumbers.empty()) {
    std::cerr << "No runs to process" << std::endl;
    return;
  }

  std::cout << "Loading " << runNumbers.size() << " runs..." << std::endl;

  // Load all runs
  std::vector<RunData> allRuns;
  for (int run : runNumbers) {
    RunData data;
    if (loadRun(baseDir, run, data)) {
      std::cout << "  Run " << run << ": " << data.entries.size() << " entries" << std::endl;
      allRuns.push_back(data);
    } else {
      std::cerr << "  Failed to load run " << run << std::endl;
    }
  }

  if (allRuns.size() < 2) {
    std::cerr << "Need at least 2 runs for comparison" << std::endl;
    return;
  }

  std::cout << "Successfully loaded " << allRuns.size() << " runs\n" << std::endl;

  // Print summary
  printSummaryTable(allRuns);

  // Compare metrics
  std::cout << "Generating comparison plots..." << std::endl;

  std::vector<std::pair<std::string, std::string>> metrics = {
      {"slope", "Slope (p1)"},
      {"intercept", "Intercept (p0)"},
      {"chi2_ndf", "Chi2/ndf"},
      {"n_points", "N Points"}};

  for (const auto &[branch, label] : metrics) {
    // Delta distributions
    drawDeltaFromReferenceByRun(allRuns, branch, label, "all_runs", outDir, true);

    // Representative channels
    auto shifted = rankCellsBySpread(allRuns, branch, true);
    if (!shifted.empty()) {
      drawRepresentativeTrendSet(allRuns, shifted, branch, label, "shifted channels",
                                 outDir, "shifted", true);
      writeRepresentativeSummary(
          std::string(outDir) + "/shifted_" + branch + ".txt", shifted, shifted, branch);
    }
  }

  // Summary graphs
  drawSummaryGraphs(allRuns, outDir);

  std::cout << "\nComparison complete. Output in: " << outDir << std::endl;
}
