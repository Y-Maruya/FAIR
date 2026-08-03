#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLine.h"
#include "TSystem.h"
#include "TLatex.h"

#include "TTree.h"

// ============================================================================
// Data Structures
// ============================================================================

const bool kDrawLegacyTH2DPlots = false;

struct InterCalibEntry {
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int ch = -1;
  int channel_index = -1;
  long long n_points = 0;
  int n_fit_points_over500 = -1;
  double intercept = -999.0;  // p0
  double slope = -999.0;       // p1
  double slope_error = -999.0;
  double chi2 = -1.0;
  int ndf = -1;
  double chi2_ndf = -1.0;
  int fit_status = 999;
  int fit_ok = 0;
  double x_mm = -999.0;
  double y_mm = -999.0;
  double hg_adc_saturation = -999.0;
};

struct RunData {
  std::string dirName = "";
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
// Quality Flag Definitions
// ============================================================================

enum HGLGQualityFlag : uint32_t {
    HGLG_GOOD                = 0,
    HGLG_LOW_STAT            = 1 << 0,  // Bit 0
    HGLG_FIT_FAILED          = 1 << 1,  // Bit 1
    HGLG_BAD_RMSPERSIGMA     = 1 << 2,  // Bit 2
    HGLG_LOW_GAIN            = 1 << 3,  // Bit 3
    HGLG_PEDESTAL_MASKED     = 1 << 6,  // Bit 6
    HGLG_MANUAL_BAD          = 1 << 7,  // Bit 7
    HGLG_MANUAL_GOOD         = 1 << 8   // Bit 8
};

struct QualityFlagEntry {
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;
  int fit_ok = 0;
  double x_mm = -999.0;
  double y_mm = -999.0;
  std::vector<uint32_t> flagsByRun;  // quality_flag for each run
  int n_fit_points_over500 = -1;
  std::vector<double> distance_rms_vals;
  std::vector<double> distance_sigma_vals;
};

struct QualityFlagRunData {
  std::string dirName = "";
  std::string label;
  std::string path;
  std::vector<QualityFlagEntry> entries;
  std::map<int, QualityFlagEntry> byCellId;
};

// ============================================================================
// Utility Functions
// ============================================================================

std::vector<std::string> parseDirectoryRanges(const char *csv) {
  std::vector<std::string> dirs;
  if (!csv) return dirs;

  if (std::string(csv) == std::string("all")) {
    gSystem->Exec("ls -d */ | sed 's#/##' > dir_list.txt");
    std::ifstream ifs("dir_list.txt");
    std::string line;
    while (std::getline(ifs, line)) {
      line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
      if (line.empty()) continue;
      dirs.push_back(line);
    }
    return dirs;
  }

  std::stringstream ss(csv);
  std::string token;
  while (std::getline(ss, token, ',')) {
    token.erase(std::remove_if(token.begin(), token.end(), ::isspace), token.end());
    if (token.empty()) continue;
    dirs.push_back(token);
  }
  return dirs;
}

bool setBranchIfExists(TTree *tree, const char *name, void *addr) {
  if (!tree || !tree->GetBranch(name)) return false;
  tree->SetBranchAddress(name, addr);
  return true;
}

bool loadRunFromFile(const std::string &filePath,
                     const std::string &dirName,
                     const std::string &label,
                     RunData &out) {
  out = RunData();
  out.dirName = dirName;
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
  setBranchIfExists(tree, "slope_error", &e.slope_error);
  setBranchIfExists(tree, "fit_status", &e.fit_status);
  out.hasFitFlags = setBranchIfExists(tree, "fit_ok", &e.fit_ok);
  bool hasX = setBranchIfExists(tree, "x_mm", &e.x_mm);
  bool hasY = setBranchIfExists(tree, "y_mm", &e.y_mm);
  out.hasCoordinates = hasX && hasY;
  setBranchIfExists(tree, "hg_adc_saturation", &e.hg_adc_saturation);

  const Long64_t nEntries = tree->GetEntries();
  out.entries.reserve(nEntries);

  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (e.layer == 24 && (e.chip == 4 || e.chip == 6 || e.chip == 8)) {
      // Skip known bad channels in layer 24
      continue;
    }
    if (e.fit_status == 999) {
      // Skip entries with invalid fit status
      continue;
    }
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

bool loadRun(const std::string &baseDir, const std::string &dirName, RunData &out) {
  const std::string path = baseDir + "/" + dirName + "/intercalib_adj_muon.root";
  return loadRunFromFile(path, dirName, dirName, out);
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
  if (branch == "hg_adc_saturation") { value = e.hg_adc_saturation; return true; }

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

SummaryStats extractStats(const RunData &runData, const std::string &what, bool onlyFitOk, bool mask_L21C1) {
  std::vector<double> values;
  values.reserve(runData.entries.size());

  for (const auto &e : runData.entries) {
    double val;
    if (mask_L21C1 && e.layer == 21 && e.chip == 1) {
      // Skip known bad channel L21C1
      continue;
    }
    if (getBranchValue(e, what, onlyFitOk, val)) {
      values.push_back(val);
    }
  }

  return computeStats(values);
}

std::vector<CellShiftInfo> rankCellsBySpread(const std::vector<RunData> &runs,
                                             const std::string &branch,
                                             bool onlyFitOk,
                                             long long minNPoints = -1) {
  std::vector<CellShiftInfo> ranked;
  if (runs.size() < 2) return ranked;

  const RunData &ref = runs.front();
  for (const auto &kv : ref.byCellId) {
    const int cellid = kv.first;
    const InterCalibEntry &ref_e = kv.second;
    if (minNPoints >= 0 && ref_e.n_points <= minNPoints) continue;

    double ref_val;
    if (!getBranchValue(ref_e, branch, onlyFitOk, ref_val)) continue;

    std::vector<double> values;
    values.push_back(ref_val);

    for (size_t i = 1; i < runs.size(); ++i) {
      auto it = runs[i].byCellId.find(cellid);
      if (it == runs[i].byCellId.end()) continue;
      if (minNPoints >= 0 && it->second.n_points <= minNPoints) continue;

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

std::vector<CellShiftInfo> excludeLowPointsAtShortFitRange(
    const std::vector<RunData> &runs,
    const std::vector<CellShiftInfo> &ranked,
    double maxFitRange = 2000.0,
    long long maxNPoints = 5000) {
  std::vector<CellShiftInfo> filtered;
  filtered.reserve(ranked.size());

  for (const auto &info : ranked) {
    bool exclude = false;
    for (const auto &run : runs) {
      auto it = run.byCellId.find(info.cellid);
      if (it == run.byCellId.end()) continue;

      const double fitRange = it->second.hg_adc_saturation - 600.0;
      if (fitRange <= maxFitRange && it->second.n_points <= maxNPoints) {
        exclude = true;
        break;
      }
    }
    if (!exclude) filtered.push_back(info);
  }

  return filtered;
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
                                       const std::string &branch, bool onlyFitOk,
                                       bool mask_L21C1 = true,
                                       bool mask_L3 = false,
                                       bool npoint_cut = false) {
  std::vector<double> deltas;
  deltas.reserve(ref.byCellId.size());

  for (const auto &kv : ref.byCellId) {
    const int cellid = kv.first;
    const InterCalibEntry &ref_e = kv.second;
    if (mask_L21C1 && int(cellid/100000) == 21 && int(cellid/10000)%10 == 1) continue;
    if (int(cellid/100000) == 24 && int(cellid/10000)%10 > 3) continue;
    if (mask_L3 && int(cellid/100000) == 3) continue;
    // if (int(cellid/100000) < 20) continue;
    double ref_val;
    if (!getBranchValue(ref_e, branch, onlyFitOk, ref_val)) continue;
    if (npoint_cut && ref_e.n_points < 10000) continue;

    auto it = target.byCellId.find(cellid);
    if (it == target.byCellId.end()) continue;

    double val;
    if (!getBranchValue(it->second, branch, onlyFitOk, val)) continue;

    if (npoint_cut && it->second.n_points < 10000) continue;
    deltas.push_back(val - ref_val);
    if (branch == "slope" && std::abs(ref_val) > 1e-6 && std::abs(val - ref_val) > 3) {
      // For slope, also collect relative difference
      std::cout << "CellID " << cellid << ": ref=" << ref_val << ", target=" << val
                << ", delta=" << (val - ref_val) << ", rel_delta=" << ((val - ref_val) / ref_val) * 100.0 << " %\n";
      std::cout << "  Ref entry: layer=" << ref_e.layer << ", chip=" << ref_e.chip << ", ch=" << ref_e.ch
                << ", n_points=" << ref_e.n_points << ", chi2_ndf=" << ref_e.chi2_ndf
                << ", fit_ok=" << (ref_e.fit_ok ? "yes" : "no") << "\n";
      std::cout << "  Target entry: layer=" << it->second.layer << ", chip=" << it->second.chip << ", ch=" << it->second.ch
                << ", n_points=" << it->second.n_points << ", chi2_ndf=" << it->second.chi2_ndf
                << ", fit_ok=" << (it->second.fit_ok ? "yes" : "no") << "\n";
    }
  }

  return deltas;
}

void drawDeltaFromReferenceByRun(const std::vector<RunData> &runs,
                                 const std::string &branch,
                                 const std::string &branchLabel,
                                 const std::string &outTag,
                                 const std::string &outDir,
                                 bool onlyFitOk = true,
                                 bool mask_L21C1 = true,
                                 bool mask_L3 = true) {
  if (runs.size() < 2) return;

  const RunData &ref = runs.front();
  std::vector<std::vector<double>> allDeltas(runs.size());
  double globalMin = 1e30;
  double globalMax = -1e30;

  for (size_t i = 1; i < runs.size(); ++i) {
    allDeltas[i] = collectDeltaValues(ref, runs[i], branch, onlyFitOk, mask_L21C1, mask_L3);
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
  c->cd(1)->SetLogy(true);
  TH1D *hFirst = nullptr;
  int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2};
  for (size_t i = 1; i < runs.size(); ++i) {
    TH1D *h = new TH1D(Form("hDelta_%zu", i), Form("%s vs %s", branchLabel.c_str(), runs[i].label.c_str()),
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
  // gPad->SetLogy(0);

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
                                bool onlyFitOk,
                                size_t maxCells = 12,
                                bool showNPoints = false,
                                long long minNPoints = -1) {
  if (cells.empty()) return;

  const int n = static_cast<int>(std::min(maxCells, cells.size()));
  const int nCols = showNPoints ? 3 : 2;
  const int nRows = (n + nCols - 1) / nCols;
  TCanvas *c = new TCanvas(Form("c_%s_%s", branch.c_str(), canvasTag.c_str()),
                           Form("%s %s", branchLabel.c_str(), label.c_str()),
                           showNPoints ? 2100 : 1400, 420 * nRows);
  c->Divide(nCols, nRows);

  for (int i = 0; i < n; ++i) {
    c->cd(i + 1);

    std::vector<double> runIndices, values;
    std::vector<long long> nPoints;

    for (size_t runIndex = 0; runIndex < runs.size(); ++runIndex) {
      const auto &run = runs[runIndex];
      auto it = run.byCellId.find(cells[i].cellid);
      if (it == run.byCellId.end()) continue;
      if (minNPoints >= 0 && it->second.n_points <= minNPoints) continue;

      double val;
      if (!getBranchValue(it->second, branch, onlyFitOk, val)) continue;

      runIndices.push_back(static_cast<double>(runIndex));
      values.push_back(val);
      nPoints.push_back(it->second.n_points);
    }

    if (!runIndices.empty()) {
      TGraph *g = new TGraph(runIndices.size(), runIndices.data(), values.data());
      g->SetMarkerColor(kBlack);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(1.2);
      g->SetLineColor(kBlack);
      g->SetLineWidth(2);
      g->SetTitle(Form("CellID %d", cells[i].cellid));
      g->GetXaxis()->SetLimits(-0.5, runs.size() - 0.5);
      g->GetXaxis()->SetTitle("Run index");
      g->GetYaxis()->SetTitle(branchLabel.c_str());
      if (showNPoints) {
        const auto mm = std::minmax_element(values.begin(), values.end());
        const double span = std::max(1e-6, *mm.second - *mm.first);
        g->SetMinimum(*mm.first - 0.12 * span);
        g->SetMaximum(*mm.second + 0.25 * span);
      }
      g->Draw("ALP");
      gPad->SetGridy();

      TLatex tl;
      tl.SetNDC(true);
      tl.SetTextSize(0.055);
      tl.DrawLatex(0.12, 0.93, Form("#%d CellID %d (L%d C%d Ch%d), RMS=%.3g",
                                    i + 1, cells[i].cellid,
                                    (cells[i].cellid / 100000),
                                    ((cells[i].cellid / 10000) % 10),
                                    (cells[i].cellid % 10000), cells[i].spread));

      if (showNPoints) {
        tl.SetNDC(false);
        tl.SetTextSize(0.035);
        tl.SetTextAlign(21);
        const auto mm = std::minmax_element(values.begin(), values.end());
        const double offset = 0.055 * std::max(1e-6, *mm.second - *mm.first);
        for (size_t j = 0; j < runIndices.size(); ++j) {
          tl.DrawLatex(runIndices[j], values[j] + offset, Form("n=%lld", nPoints[j]));
        }
      }
    }

  }

  c->SaveAs((outDir + "/" + branch + "_" + canvasTag + ".pdf").c_str());
}

double computePearsonCorrelation(const std::vector<double> &x,
                                 const std::vector<double> &y) {
  if (x.size() != y.size() || x.size() < 2) return 0.0;

  const double meanX = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
  const double meanY = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
  double cov = 0.0;
  double varX = 0.0;
  double varY = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    const double dx = x[i] - meanX;
    const double dy = y[i] - meanY;
    cov += dx * dy;
    varX += dx * dx;
    varY += dy * dy;
  }
  if (varX <= 0.0 || varY <= 0.0) return 0.0;
  return cov / std::sqrt(varX * varY);
}

void drawNPointsSlopeCorrelations(const std::vector<RunData> &runs,
                                  const std::vector<CellShiftInfo> &ranked,
                                  const std::string &outDir,
                                  size_t maxCells = 30,
                                  const std::string &subDirName =
                                      "slope_vs_n_points_shifted_channels") {
  if (runs.empty() || ranked.empty()) return;

  const std::string corrOutDir = outDir + "/" + subDirName;
  gSystem->mkdir(corrOutDir.c_str(), true);

  const size_t nCells = std::min(maxCells, ranked.size());
  int nWritten = 0;

  for (size_t rank = 0; rank < nCells; ++rank) {
    const auto &info = ranked[rank];
    std::vector<double> nPoints;
    std::vector<double> slopes;
    std::vector<std::string> labels;
    const InterCalibEntry *firstEntry = nullptr;

    for (const auto &run : runs) {
      auto it = run.byCellId.find(info.cellid);
      if (it == run.byCellId.end() || !it->second.fit_ok) continue;
      if (!std::isfinite(it->second.slope) || it->second.n_points <= 0) continue;
      if (!firstEntry) firstEntry = &it->second;
      nPoints.push_back(static_cast<double>(it->second.n_points));
      slopes.push_back(it->second.slope);
      labels.push_back(run.label);
    }

    if (nPoints.size() < 2 || !firstEntry) continue;

    const auto xRange = std::minmax_element(nPoints.begin(), nPoints.end());
    const auto yRange = std::minmax_element(slopes.begin(), slopes.end());
    double xMin = *xRange.first;
    double xMax = *xRange.second;
    double yMin = *yRange.first;
    double yMax = *yRange.second;
    const double xSpan = std::max(1.0, xMax - xMin);
    const double ySpan = std::max(1e-6, yMax - yMin);
    xMin = std::max(0.0, xMin - 0.12 * xSpan);
    xMax += 0.12 * xSpan;
    yMin -= 0.18 * ySpan;
    yMax += 0.28 * ySpan;

    TCanvas *c = new TCanvas(Form("c_slope_vs_npoints_cell%d", info.cellid),
                             Form("Slope vs n_points cell %d", info.cellid),
                             1100, 850);
    TGraph *g = new TGraph(nPoints.size(), nPoints.data(), slopes.data());
    g->SetTitle(Form("Rank %zu, CellID %d (L%d C%d Ch%d);n_points;Slope (p1)",
                     rank + 1, info.cellid, firstEntry->layer, firstEntry->chip,
                     firstEntry->ch));
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.4);
    g->SetMarkerColor(kBlue + 1);
    g->SetLineColor(kBlue + 1);
    g->GetXaxis()->SetLimits(xMin, xMax);
    g->SetMinimum(yMin);
    g->SetMaximum(yMax);
    g->Draw("AP");
    gPad->SetGrid();

    TLatex tl;
    tl.SetNDC(true);
    tl.SetTextSize(0.035);
    tl.DrawLatex(0.15, 0.86, Form("slope RMS across runs: %.4g", info.spread));
    tl.DrawLatex(0.15, 0.81, Form("Pearson r: %.3f", computePearsonCorrelation(nPoints, slopes)));
    tl.DrawLatex(0.15, 0.76, Form("Points: %zu", nPoints.size()));

    tl.SetNDC(false);
    tl.SetTextSize(0.026);
    tl.SetTextAlign(12);
    const double labelDx = 0.018 * (xMax - xMin);
    const double labelDy = 0.018 * (yMax - yMin);
    for (size_t i = 0; i < nPoints.size(); ++i) {
      tl.DrawLatex(nPoints[i] + labelDx, slopes[i] + labelDy, labels[i].c_str());
    }

    const std::string pdfPath =
        corrOutDir + Form("/rank%02zu_cell%d_L%02d_C%02d_ch%02d_slope_vs_n_points.pdf",
                          rank + 1, info.cellid, firstEntry->layer, firstEntry->chip,
                          firstEntry->ch);
    c->SaveAs(pdfPath.c_str());
    ++nWritten;
    delete c;
  }

  if (nWritten > 0) {
    std::cout << "  Wrote slope-vs-n_points shifted-channel plots: "
              << nWritten << " files in " << corrOutDir << std::endl;
  }
}

void writeTopSlopeTrendSummary(const std::string &fileName,
                               const std::vector<RunData> &runs,
                               const std::vector<CellShiftInfo> &ranked,
                               size_t maxCells = 30) {
  std::ofstream ofs(fileName.c_str());
  if (!ofs) return;

  ofs << "# Slope-shift ranking by RMS across runs\n";
  ofs << "# rank cellid layer chip channel slope_rms run slope n_points fit_range\n";
  const size_t n = std::min(maxCells, ranked.size());
  for (size_t i = 0; i < n; ++i) {
    const auto &info = ranked[i];
    for (const auto &run : runs) {
      auto it = run.byCellId.find(info.cellid);
      if (it == run.byCellId.end() || !it->second.fit_ok) continue;
      ofs << i + 1 << " " << info.cellid << " "
          << it->second.layer << " " << it->second.chip << " " << it->second.ch << " "
          << info.spread << " " << run.label << " " << it->second.slope << " "
          << it->second.n_points << " " << it->second.hg_adc_saturation - 600.0 << "\n";
    }
    ofs << "\n";
  }
}

void drawSlopeRmsDistribution(const std::vector<CellShiftInfo> &ranked,
                              const std::string &outDir,
                              const std::string &tag = "") {
  if (ranked.empty()) return;

  std::vector<double> rmsValues;
  rmsValues.reserve(ranked.size());
  for (const auto &info : ranked) {
    if (std::isfinite(info.spread) && info.spread >= 0.0) {
      rmsValues.push_back(info.spread);
    }
  }
  if (rmsValues.empty()) return;

  const SummaryStats stats = computeStats(rmsValues);
  const double xMax = std::max(1e-6, stats.max * 1.05);
  const std::string suffix = tag.empty() ? "" : "_" + tag;
  const std::string titleSuffix = tag.empty() ? "" : " (" + tag + ")";
  TH1D *h = new TH1D(("h_slope_rms_across_runs" + suffix).c_str(),
                     ("Slope RMS across runs" + titleSuffix +
                      ";RMS of slope across runs;Channels").c_str(),
                     100, 0.0, xMax);
  for (double rms : rmsValues) h->Fill(rms);
  h->SetLineColor(kBlue + 1);
  h->SetFillColorAlpha(kBlue + 1, 0.25);
  h->SetLineWidth(2);

  TCanvas *c = new TCanvas(("c_slope_rms_across_runs" + suffix).c_str(),
                           ("Slope RMS across runs" + titleSuffix).c_str(), 1400, 600);
  c->Divide(2, 1);

  c->cd(1);
  h->Draw("HIST");
  gPad->SetGridy();

  TLatex tl;
  tl.SetNDC(true);
  tl.SetTextSize(0.035);
  tl.DrawLatex(0.58, 0.84, Form("Channels: %d", stats.n));
  tl.DrawLatex(0.58, 0.79, Form("Mean RMS: %.4g", stats.mean));
  tl.DrawLatex(0.58, 0.74, Form("Max RMS: %.4g", stats.max));

  c->cd(2);
  TH1D *hLog = dynamic_cast<TH1D *>(h->Clone(
      ("h_slope_rms_across_runs_log" + suffix).c_str()));
  hLog->Draw("HIST");
  gPad->SetLogy();
  gPad->SetGridy();

  c->SaveAs((outDir + "/Slope_RMS_across_runs" + suffix + ".pdf").c_str());
}

void drawMeanNPointsVsHgSaturationRms(const std::vector<RunData> &runs,
                                      const std::string &outDir,
                                      bool onlyFitOk = true,
                                      const std::string &tag = "") {
  if (runs.size() < 2) return;

  std::set<int> allCellIds;
  for (const auto &run : runs) {
    for (const auto &[cellid, entry] : run.byCellId) {
      allCellIds.insert(cellid);
    }
  }

  std::vector<double> meanNPoints;
  std::vector<double> hgSaturationRms;
  std::vector<int> cellids;
  for (int cellid : allCellIds) {
    std::vector<double> nPointsVals;
    std::vector<double> hgVals;

    for (const auto &run : runs) {
      auto it = run.byCellId.find(cellid);
      if (it == run.byCellId.end()) continue;
      if (onlyFitOk && !it->second.fit_ok) continue;
      if (it->second.n_points <= 0) continue;
      if (!std::isfinite(it->second.hg_adc_saturation) ||
          it->second.hg_adc_saturation <= 0.0) continue;

      nPointsVals.push_back(static_cast<double>(it->second.n_points));
      hgVals.push_back(it->second.hg_adc_saturation);
    }

    if (hgVals.size() < 2 || nPointsVals.empty()) continue;
    SummaryStats nStats = computeStats(nPointsVals);
    SummaryStats hgStats = computeStats(hgVals);
    meanNPoints.push_back(nStats.mean);
    hgSaturationRms.push_back(hgStats.rms);
    cellids.push_back(cellid);
  }

  if (meanNPoints.empty()) return;

  const double corr = computePearsonCorrelation(meanNPoints, hgSaturationRms);
  const auto xRange = std::minmax_element(meanNPoints.begin(), meanNPoints.end());
  const auto yRange = std::minmax_element(hgSaturationRms.begin(), hgSaturationRms.end());
  double xMin = *xRange.first;
  double xMax = *xRange.second;
  double yMin = *yRange.first;
  double yMax = *yRange.second;
  const double xSpan = std::max(1.0, xMax - xMin);
  const double ySpan = std::max(1e-6, yMax - yMin);
  xMin = std::max(0.0, xMin - 0.08 * xSpan);
  xMax += 0.08 * xSpan;
  yMin = std::max(0.0, yMin - 0.08 * ySpan);
  yMax += 0.12 * ySpan;

  const std::string suffix = tag.empty() ? "" : "_" + tag;
  TCanvas *c = new TCanvas(("c_mean_npoints_vs_hg_saturation_rms" + suffix).c_str(),
                           "Mean n_points vs HG saturation RMS",
                           1100, 850);
  TGraph *g = new TGraph(meanNPoints.size(), meanNPoints.data(), hgSaturationRms.data());
  g->SetTitle("All channels;Mean n_points across runs;RMS of HG ADC saturation across runs");
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.7);
  g->SetMarkerColor(kBlue + 1);
  g->GetXaxis()->SetLimits(xMin, xMax);
  g->SetMinimum(yMin);
  g->SetMaximum(yMax);
  g->Draw("AP");
  gPad->SetGrid();

  TLatex tl;
  tl.SetNDC(true);
  tl.SetTextSize(0.035);
  tl.DrawLatex(0.15, 0.86, Form("Channels: %zu", meanNPoints.size()));
  tl.DrawLatex(0.15, 0.81, Form("Pearson r: %.3f", corr));

  c->SaveAs((outDir + "/MeanNPoints_vs_HGADCSaturation_RMS" + suffix + ".pdf").c_str());

  std::ofstream ofs((outDir + "/MeanNPoints_vs_HGADCSaturation_RMS" + suffix + ".txt").c_str());
  if (ofs) {
    ofs << "# cellid mean_n_points hg_adc_saturation_rms\n";
    for (size_t i = 0; i < meanNPoints.size(); ++i) {
      ofs << cellids[i] << " " << meanNPoints[i] << " " << hgSaturationRms[i] << "\n";
    }
  }
}

std::vector<CellShiftInfo> collectHighMeanNPointsHighHgSaturationRmsCells(
    const std::vector<RunData> &runs,
    double minMeanNPoints = 8000.0,
    double minHgSaturationRms = 70.0,
    bool onlyFitOk = true) {
  std::vector<CellShiftInfo> selected;
  if (runs.size() < 2) return selected;

  std::set<int> allCellIds;
  for (const auto &run : runs) {
    for (const auto &[cellid, entry] : run.byCellId) {
      allCellIds.insert(cellid);
    }
  }

  for (int cellid : allCellIds) {
    std::vector<double> nPointsVals;
    std::vector<double> hgVals;
    const InterCalibEntry *firstEntry = nullptr;

    for (const auto &run : runs) {
      auto it = run.byCellId.find(cellid);
      if (it == run.byCellId.end()) continue;
      if (onlyFitOk && !it->second.fit_ok) continue;
      if (it->second.n_points <= 0) continue;
      if (!std::isfinite(it->second.hg_adc_saturation) ||
          it->second.hg_adc_saturation <= 0.0) continue;

      if (!firstEntry) firstEntry = &it->second;
      nPointsVals.push_back(static_cast<double>(it->second.n_points));
      hgVals.push_back(it->second.hg_adc_saturation);
    }

    if (!firstEntry || hgVals.size() < 2 || nPointsVals.empty()) continue;

    const SummaryStats nStats = computeStats(nPointsVals);
    const SummaryStats hgStats = computeStats(hgVals);
    if (nStats.mean <= minMeanNPoints || hgStats.rms <= minHgSaturationRms) {
      continue;
    }

    CellShiftInfo info;
    info.cellid = cellid;
    info.spread = hgStats.rms;
    info.refValue = nStats.mean;
    info.meanValue = hgStats.mean;
    info.x_mm = firstEntry->x_mm;
    info.y_mm = firstEntry->y_mm;
    selected.push_back(info);
  }

  std::sort(selected.begin(), selected.end(),
            [](const CellShiftInfo &a, const CellShiftInfo &b) {
              if (a.spread != b.spread) return a.spread > b.spread;
              return a.cellid < b.cellid;
            });

  return selected;
}

void writeHighMeanNPointsHighHgSaturationRmsSummary(
    const std::string &fileName,
    const std::vector<CellShiftInfo> &cells) {
  std::ofstream ofs(fileName.c_str());
  if (!ofs) return;

  ofs << "# Channels with mean n_points > 8000 and HG ADC saturation RMS > 70\n";
  ofs << "# rank cellid hg_adc_saturation_rms mean_n_points mean_hg_adc_saturation x_mm y_mm\n";
  for (size_t i = 0; i < cells.size(); ++i) {
    const auto &info = cells[i];
    ofs << i + 1 << " " << info.cellid << " "
        << info.spread << " " << info.refValue << " "
        << info.meanValue << " " << info.x_mm << " " << info.y_mm << "\n";
  }
}

std::set<int> collectChannelsForHgSaturationStability(
    const std::vector<RunData> &runs,
    double minMeanNPoints = 10000.0,
    bool onlyFitOk = true) {
  std::set<int> selected;
  if (runs.size() < 2) return selected;

  std::set<int> allCellIds;
  for (const auto &run : runs) {
    for (const auto &[cellid, entry] : run.byCellId) {
      allCellIds.insert(cellid);
    }
  }

  for (int cellid : allCellIds) {
    std::vector<double> nPointsVals;
    bool excludedLayer = false;

    for (const auto &run : runs) {
      auto it = run.byCellId.find(cellid);
      if (it == run.byCellId.end()) continue;
      if (it->second.layer == 1 || it->second.layer == 3 ||
          (it->second.layer == 13 && it->second.chip == 5)) {
        excludedLayer = true;
        break;
      }
      if (onlyFitOk && !it->second.fit_ok) continue;
      if (it->second.n_points > 0) {
        nPointsVals.push_back(static_cast<double>(it->second.n_points));
      }
    }

    if (excludedLayer || nPointsVals.empty()) continue;
    if (computeStats(nPointsVals).mean >= minMeanNPoints) {
      selected.insert(cellid);
    }
  }

  return selected;
}

std::vector<RunData> filterRunsKeepingCellIds(const std::vector<RunData> &runs,
                                             const std::set<int> &keptCellIds) {
  std::vector<RunData> filteredRuns;
  filteredRuns.reserve(runs.size());

  for (const auto &run : runs) {
    RunData filtered = run;
    filtered.entries.clear();
    filtered.byCellId.clear();
    filtered.entries.reserve(run.entries.size());

    for (const auto &entry : run.entries) {
      if (!keptCellIds.count(entry.cellid)) continue;
      filtered.entries.push_back(entry);
      filtered.byCellId[entry.cellid] = entry;
    }

    filteredRuns.push_back(filtered);
  }

  return filteredRuns;
}

void drawHgSaturationDistributionByRun(const std::vector<RunData> &runs,
                                       const std::string &outDir,
                                       const std::string &tag) {
  if (runs.empty()) return;

  std::vector<std::vector<double>> valuesByRun(runs.size());
  double globalMin = 1e30;
  double globalMax = -1e30;

  for (size_t i = 0; i < runs.size(); ++i) {
    for (const auto &entry : runs[i].entries) {
      if (!entry.fit_ok) continue;
      if (!std::isfinite(entry.hg_adc_saturation) ||
          entry.hg_adc_saturation <= 0.0) continue;
      valuesByRun[i].push_back(entry.hg_adc_saturation);
    }
    if (!valuesByRun[i].empty()) {
      auto range = std::minmax_element(valuesByRun[i].begin(), valuesByRun[i].end());
      globalMin = std::min(globalMin, *range.first);
      globalMax = std::max(globalMax, *range.second);
    }
  }

  if (globalMin > globalMax) return;
  if (globalMin == globalMax) {
    globalMin -= 1.0;
    globalMax += 1.0;
  }

  const double span = std::max(1.0, globalMax - globalMin);
  const double xMin = globalMin - 0.05 * span;
  const double xMax = globalMax + 0.05 * span;
  const int nbins = 80;
  int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
                  kOrange + 7, kCyan + 2, kViolet + 1, kGray + 2};

  TCanvas *c = new TCanvas(("c_hg_saturation_distribution_" + tag).c_str(),
                           "HG ADC saturation distribution", 1200, 800);
  TLegend *leg = new TLegend(0.68, 0.18, 0.94, 0.88);
  leg->SetBorderSize(1);
  leg->SetFillStyle(0);

  bool first = true;
  for (size_t i = 0; i < runs.size(); ++i) {
    if (valuesByRun[i].empty()) continue;
    TH1D *h = new TH1D(Form("h_hg_saturation_%s_%zu", tag.c_str(), i),
                       "HG ADC saturation distribution;HG ADC saturation;Channels",
                       nbins, xMin, xMax);
    for (double v : valuesByRun[i]) h->Fill(v);
    h->SetLineColor(colors[i % 8]);
    h->SetLineWidth(2);
    h->SetStats(0);
    h->Draw(first ? "HIST" : "HIST SAME");
    first = false;
    leg->AddEntry(h, Form("%s (N=%zu)", runs[i].label.c_str(), valuesByRun[i].size()), "l");
  }

  leg->Draw();
  gPad->SetGridy();
  c->SaveAs((outDir + "/HGADCSaturation_Distribution_" + tag + ".pdf").c_str());
}

void drawHgSaturationRmsDistribution(const std::vector<RunData> &runs,
                                     const std::string &outDir,
                                     const std::string &tag) {
  auto ranked = rankCellsBySpread(runs, "hg_adc_saturation", true);
  if (ranked.empty()) return;

  std::vector<double> rmsValues;
  rmsValues.reserve(ranked.size());
  for (const auto &info : ranked) {
    if (std::isfinite(info.spread) && info.spread >= 0.0) {
      rmsValues.push_back(info.spread);
    }
  }
  if (rmsValues.empty()) return;

  const SummaryStats stats = computeStats(rmsValues);
  const double xMax = std::max(1e-6, stats.max * 1.05);
  TH1D *h = new TH1D(("h_hg_saturation_rms_" + tag).c_str(),
                     "HG ADC saturation RMS across runs;RMS of HG ADC saturation;Channels",
                     80, 0.0, xMax);
  for (double v : rmsValues) h->Fill(v);
  h->SetLineColor(kBlue + 1);
  h->SetFillColorAlpha(kBlue + 1, 0.25);
  h->SetLineWidth(2);

  TCanvas *c = new TCanvas(("c_hg_saturation_rms_" + tag).c_str(),
                           "HG ADC saturation RMS across runs", 1100, 800);
  h->Draw("HIST");
  gPad->SetGridy();

  TLatex tl;
  tl.SetNDC(true);
  tl.SetTextSize(0.035);
  tl.DrawLatex(0.58, 0.84, Form("Channels: %d", stats.n));
  tl.DrawLatex(0.58, 0.79, Form("Mean RMS: %.4g", stats.mean));
  tl.DrawLatex(0.58, 0.74, Form("Max RMS: %.4g", stats.max));

  c->SaveAs((outDir + "/HGADCSaturation_RMS_" + tag + ".pdf").c_str());

  std::ofstream ofs((outDir + "/HGADCSaturation_RMS_" + tag + ".txt").c_str());
  if (ofs) {
    ofs << "# cellid hg_adc_saturation_rms mean_hg_adc_saturation\n";
    for (const auto &info : ranked) {
      ofs << info.cellid << " " << info.spread << " " << info.meanValue << "\n";
    }
  }
}

void drawShiftedChannelAccumulatedHistograms(const std::vector<RunData> &runs,
                                             const std::vector<CellShiftInfo> &ranked,
                                             const std::string &outDir,
                                             size_t maxCells = 30,
                                             const std::string &subDirName =
                                                 "shifted_channel_accumulated",
                                             const std::string &scoreLabel =
                                                 "Slope RMS across runs") {
  if (runs.empty() || ranked.empty()) return;

  const std::string histOutDir = outDir + "/" + subDirName;
  gSystem->mkdir(histOutDir.c_str(), true);

  const size_t nCells = std::min(maxCells, ranked.size());
  for (size_t rank = 0; rank < nCells; ++rank) {
    const auto &info = ranked[rank];
    const InterCalibEntry *refEntry = nullptr;
    for (const auto &run : runs) {
      auto it = run.byCellId.find(info.cellid);
      if (it != run.byCellId.end()) {
        refEntry = &it->second;
        break;
      }
    }
    if (!refEntry) continue;

    const int layer = refEntry->layer;
    const int chip = refEntry->chip;
    const int channel = refEntry->ch;
    const std::string histPath =
        Form("AccumulatedHistograms/h_accum_L%02d_C%02d_ch%02d", layer, chip, channel);
    const std::string pdfPath =
        histOutDir + Form("/rank%02zu_cell%d_L%02d_C%02d_ch%02d.pdf",
                          rank + 1, info.cellid, layer, chip, channel);

    TCanvas *c = new TCanvas(Form("c_accum_cell%d", info.cellid),
                             Form("Accumulated histogram cell %d", info.cellid),
                             1000, 850);
    bool pdfOpened = false;
    int nPages = 0;

    for (const auto &run : runs) {
      TFile *file = TFile::Open(run.path.c_str(), "READ");
      if (!file || file->IsZombie()) {
        if (file) file->Close();
        continue;
      }

      TH2D *h = dynamic_cast<TH2D *>(file->Get(histPath.c_str()));
      if (!h) {
        file->Close();
        continue;
      }

      c->Clear();
      c->SetRightMargin(0.14);
      c->SetLogz();
      h->SetStats(0);
      h->GetXaxis()->SetRangeUser(-100.0, 300.0);
      h->SetTitle(Form("Rank %zu, CellID %d (L%02d C%02d ch%02d), run %s",
                       rank + 1, info.cellid, layer, chip, channel, run.label.c_str()));
      h->Draw("COLZ");

      auto entryIt = run.byCellId.find(info.cellid);
      TLine *fitMinLine = nullptr;
      TLine *fitMaxLine = nullptr;
      TLatex tl;
      tl.SetNDC(true);
      tl.SetTextSize(0.03);
      tl.DrawLatex(0.12, 0.92, Form("%s: %.4g", scoreLabel.c_str(), info.spread));
      if (entryIt != run.byCellId.end()) {
        const double hgFitMin = 50.0;
        const double hgFitMax = entryIt->second.hg_adc_saturation - 600.0;
        const double xMin = -100.0;
        const double xMax = 300.0;

        fitMinLine = new TLine(xMin, hgFitMin, xMax, hgFitMin);
        fitMinLine->SetLineColor(kRed + 1);
        fitMinLine->SetLineStyle(2);
        fitMinLine->SetLineWidth(3);
        fitMinLine->Draw();

        fitMaxLine = new TLine(xMin, hgFitMax, xMax, hgFitMax);
        fitMaxLine->SetLineColor(kRed + 1);
        fitMaxLine->SetLineStyle(2);
        fitMaxLine->SetLineWidth(3);
        fitMaxLine->Draw();

        tl.DrawLatex(0.12, 0.88,
                     Form("slope: %.6g, intercept: %.6g",
                          entryIt->second.slope, entryIt->second.intercept));
        tl.DrawLatex(0.12, 0.84,
                     Form("chi2/ndf: %.6g, n_points: %lld, n_fit_points_over500: %d",
                          entryIt->second.chi2_ndf, entryIt->second.n_points,
                          entryIt->second.n_fit_points_over500));
        tl.SetTextColor(kRed + 1);
        tl.DrawLatex(0.12, 0.80,
                     Form("HG fit range: %.0f < HG < %.0f",
                          hgFitMin, hgFitMax));
        tl.SetTextColor(kBlack);
      }

      if (!pdfOpened) {
        c->Print((pdfPath + "[").c_str());
        pdfOpened = true;
      }
      c->Print(pdfPath.c_str());
      ++nPages;
      c->Clear();
      delete fitMinLine;
      delete fitMaxLine;
      file->Close();
    }

    if (pdfOpened) {
      c->Print((pdfPath + "]").c_str());
      std::cout << "  Wrote " << nPages << " accumulated-histogram pages: "
                << pdfPath << std::endl;
    }
    delete c;
  }
}

std::vector<CellShiftInfo> collectChannelsByLayerChip(const std::vector<RunData> &runs,
                                                      int targetLayer,
                                                      int targetChip) {
  std::map<int, CellShiftInfo> cellsById;

  for (const auto &run : runs) {
    for (const auto &[cellid, entry] : run.byCellId) {
      if (entry.layer != targetLayer || entry.chip != targetChip) continue;
      if (cellsById.count(cellid)) continue;

      CellShiftInfo info;
      info.cellid = cellid;
      info.spread = static_cast<double>(entry.ch);
      info.refValue = entry.hg_adc_saturation;
      info.meanValue = static_cast<double>(entry.n_points);
      info.x_mm = entry.x_mm;
      info.y_mm = entry.y_mm;
      cellsById[cellid] = info;
    }
  }

  std::vector<CellShiftInfo> cells;
  cells.reserve(cellsById.size());
  for (const auto &[cellid, info] : cellsById) {
    cells.push_back(info);
  }

  std::sort(cells.begin(), cells.end(),
            [](const CellShiftInfo &a, const CellShiftInfo &b) {
              const int chA = a.cellid % 10000;
              const int chB = b.cellid % 10000;
              if (chA != chB) return chA < chB;
              return a.cellid < b.cellid;
            });

  return cells;
}

void drawSummaryGraphs(const std::vector<RunData> &runs, const std::string &outDir) {
  std::vector<std::string> branches = {"slope", "intercept", "chi2_ndf", "n_points", "hg_adc_saturation"};
  std::vector<std::string> labels = {"Slope (p1)", "Intercept (p0)", "Chi2/ndf", "N Points", "HG ADC Saturation"};

  for (size_t b = 0; b < branches.size(); ++b) {
    std::vector<double> runNumbers;
    std::vector<double> means;
    std::vector<double> rmss;
    std::vector<std::string> xLabels;

    for (size_t i = 0; i < runs.size(); ++i) {
      SummaryStats stats = extractStats(runs[i], branches[b], true,true);
      if (stats.n == 0) continue;

      runNumbers.push_back(static_cast<double>(xLabels.size() + 1));
      means.push_back(stats.mean);
      rmss.push_back(stats.rms);
      xLabels.push_back(runs[i].label);
    }

    if (means.empty()) continue;

    TCanvas *c = new TCanvas(Form("c_summary_%s", branches[b].c_str()),
                             Form("Summary: %s", labels[b].c_str()), 800, 600);

    double yMin = means.front() - rmss.front();
    double yMax = means.front() + rmss.front();
    for (size_t i = 0; i < means.size(); ++i) {
      yMin = std::min(yMin, means[i] - rmss[i]);
      yMax = std::max(yMax, means[i] + rmss[i]);
    }
    const double ySpan = std::max(1e-6, yMax - yMin);

    TH1D *frame = new TH1D(Form("hFrame_%s", branches[b].c_str()),
                           Form("%s;RunNumber-RunNumber;%s", labels[b].c_str(), labels[b].c_str()),
                           static_cast<int>(xLabels.size()), 0.5, xLabels.size() + 0.5);
    frame->SetMinimum(yMin - 0.15 * ySpan);
    frame->SetMaximum(yMax + 0.15 * ySpan);
    frame->SetStats(0);
    for (size_t i = 0; i < xLabels.size(); ++i) {
      frame->GetXaxis()->SetBinLabel(static_cast<int>(i) + 1, xLabels[i].c_str());
    }
    frame->Draw();

    std::vector<double> xErrors(means.size(), 0.0);
    TGraphErrors *gMean = new TGraphErrors(means.size(), runNumbers.data(), means.data(),
                                           xErrors.data(), rmss.data());
    gMean->SetTitle(Form("%s vs run", labels[b].c_str()));
    gMean->SetMarkerColor(kBlack);
    gMean->SetMarkerStyle(20);
    gMean->SetMarkerSize(1.2);
    gMean->SetLineColor(kBlack);
    gMean->SetLineWidth(2);
    gMean->GetXaxis()->SetTitle("RunNumber-RunNumber");
    gMean->GetYaxis()->SetTitle(labels[b].c_str());
    gMean->Draw("P SAME");

    c->SaveAs((outDir + "/Summary_" + branches[b] + ".pdf").c_str());
  }
}

void drawParameterDistributions(const std::vector<RunData> &runs, const std::string &outDir) {
  std::vector<std::string> branches = {"slope", "intercept", "chi2_ndf", "n_points", "hg_adc_saturation"};
  std::vector<std::string> labels = {"Slope (p1)", "Intercept (p0)", "Chi2/ndf", "N Points", "HG ADC Saturation"};

  int colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7, 
                  kCyan + 2, kViolet + 2, kPink, kSpring, kTeal, kAzure, kRose};

  for (size_t b = 0; b < branches.size(); ++b) {
    // Collect global min/max for consistent axis
    double global_min = 1e30;
    double global_max = -1e30;
    std::vector<std::vector<double>> all_values(runs.size());

    for (size_t i = 0; i < runs.size(); ++i) {
      for (const auto &e : runs[i].entries) {
        double val;
        if (getBranchValue(e, branches[b], true, val)) {
          all_values[i].push_back(val);
          global_min = std::min(global_min, val);
          global_max = std::max(global_max, val);
        }
      }
    }

    if (global_min > global_max) continue;
    if (global_min == global_max) {
      global_min -= 0.1;
      global_max += 0.1;
    }

    const double range_span = std::max(1e-6, global_max - global_min);
    double x_min = global_min - 0.05 * range_span;
    double x_max = global_max + 0.05 * range_span;
    const int nbins = 100;

    TCanvas *c = new TCanvas(Form("c_dist_%s", branches[b].c_str()),
                             Form("%s Distribution", labels[b].c_str()),
                             1200, 600);
    c->Divide(2, 1);

    // Left: overlaid histograms
    c->cd(1);
    TH1D *hFirst = nullptr;
    TLegend *leg = new TLegend(0.6, 0.55, 0.90, 0.90);
    leg->SetFillStyle(0);
    leg->SetBorderSize(1);

    for (size_t i = 0; i < runs.size(); ++i) {
      if (all_values[i].empty()) continue;
      if (branches[b] == "slope"){
        x_min = std::max(20.0, x_min);
        x_max = std::min(40.0, x_max);
      }
      TH1D *h = new TH1D(Form("h_dist_%s_run%zu", branches[b].c_str(), i),
                         Form("%s - Run %s", labels[b].c_str(), runs[i].label.c_str()),
                         nbins, x_min, x_max);
      h->SetLineColor(colors[i % 12]);
      h->SetLineWidth(2);
      h->SetFillStyle(0);

      for (double val : all_values[i]) {
        h->Fill(val);
      }

      if (i == 0) {
        hFirst = h;
        hFirst->Draw();
      } else {
        h->Draw("same");
      }

      leg->AddEntry(h, runs[i].label.c_str(), "l");
    }

    leg->Draw();
    // gPad->SetGridy();
    gStyle->SetOptStat(0);

    // Right: statistics table
    c->cd(2);
    gPad->SetMargin(0.1, 0.1, 0.1, 0.1);
    
    TH1D *hStats = new TH1D("hStats", "", 1, 0, 1);
    hStats->SetStats(0);
    hStats->Draw();

    TLatex tl;
    tl.SetNDC(true);
    tl.SetTextSize(0.035);
    tl.SetTextFont(42);

    double y_pos = 0.95;
    tl.DrawLatex(0.15, y_pos, Form("%s Statistics", labels[b].c_str()));
    y_pos -= 0.08;

    for (size_t i = 0; i < runs.size(); ++i) {
      if (all_values[i].empty()) continue;

      SummaryStats stats = computeStats(all_values[i]);
      tl.SetTextColor(colors[i % 12]);
      tl.DrawLatex(0.15, y_pos, Form("%s:", runs[i].label.c_str()));
      y_pos -= 0.05;
      tl.DrawLatex(0.20, y_pos, Form("  Mean: %.4f", stats.mean));
      y_pos -= 0.04;
      tl.DrawLatex(0.20, y_pos, Form("  RMS:  %.4f", stats.rms));
      y_pos -= 0.04;
      tl.DrawLatex(0.20, y_pos, Form("  Min:  %.4f", stats.min));
      y_pos -= 0.04;
      tl.DrawLatex(0.20, y_pos, Form("  Max:  %.4f", stats.max));
      y_pos -= 0.04;
      tl.DrawLatex(0.20, y_pos, Form("  N:    %d", stats.n));
      y_pos -= 0.08;
      tl.SetTextColor(kBlack);
    }

    c->SaveAs((outDir + "/Distribution_" + branches[b] + ".pdf").c_str());
  }
}

void printSummaryTable(const std::vector<RunData> &runs) {
  std::vector<std::string> branches = {"slope", "intercept", "chi2_ndf", "n_points", "hg_adc_saturation"};
  std::vector<std::string> labels = {"Slope (p1)", "Intercept (p0)", "Chi2/ndf", "N Points", "HG ADC Saturation"};

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
      SummaryStats stats = extractStats(run, branches[b], false, false);
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
// Quality Flag Comparison Functions
// ============================================================================

std::string flagToString(uint32_t flag) {
  std::string result;
  if (flag == 0) return "GOOD";
  if (flag & HGLG_LOW_STAT) result += "LOW_STAT|";
  if (flag & HGLG_FIT_FAILED) result += "FIT_FAILED|";
  if (flag & HGLG_BAD_RMSPERSIGMA) result += "BAD_RMS|";
  if (flag & HGLG_LOW_GAIN) result += "LOW_GAIN|";
  if (flag & HGLG_PEDESTAL_MASKED) result += "PED_MASKED|";
  if (flag & HGLG_MANUAL_BAD) result += "MANUAL_BAD|";
  if (flag & HGLG_MANUAL_GOOD) result += "MANUAL_GOOD|";
  if (!result.empty()) result.pop_back();  // Remove trailing '|'
  return result;
}

void addLowGainFlags(std::vector<QualityFlagRunData> &runs,
                     int maxNFitPointsOver500 = 50) {
  for (auto &run : runs) {
    int count = 0;
    for (auto &[cellid, entry] : run.byCellId) {
      if (entry.n_fit_points_over500 < 0) continue;
      if (entry.n_fit_points_over500 <= maxNFitPointsOver500) {
        if (!entry.flagsByRun.empty()) {
          entry.flagsByRun[0] |= HGLG_LOW_GAIN;
          ++count;
        }
      }
    }
    std::cout << "  Run " << run.dirName << ": LOW_GAIN flagged "
              << count << " channels (n_fit_points_over500 <= "
              << maxNFitPointsOver500 << ")" << std::endl;
  }
}

void addRMSPerSigmaFlags(std::vector<QualityFlagRunData> &runs, bool onlyFitOk = true, bool mask_L21C1 = true, bool mask_L13C5 = true) {
  // For each run, calculate threshold based on that run's RMS/Sigma ratios
  for (auto &run : runs) {
    // First pass: collect RMS/Sigma ratios for this run
    std::vector<double> run_ratios;
    for (auto &[cellid, entry] : run.byCellId) {
      if (entry.distance_rms_vals.empty() || entry.distance_sigma_vals.empty()) {
        continue;
      }
      if (mask_L21C1 && int(cellid/100000) == 21 && int(cellid/10000)%10 == 1) {
        // Skip known bad channel L21C1
        continue;
      }
      if (mask_L13C5 && int(cellid/100000) == 13 && int(cellid/10000)%10 == 5) {
        // Skip known bad channel L13C5
        continue;
      }
      if (onlyFitOk && !entry.fit_ok) {
        continue;
      }
      
      double dist_rms = entry.distance_rms_vals[0];
      double dist_sigma = entry.distance_sigma_vals[0];
      
      if (dist_sigma > 0 && dist_rms >= 0) {
        double ratio = dist_rms / dist_sigma;
        run_ratios.push_back(ratio);
      }
    }
    
    // Calculate threshold for this run as RMS of ratios * 10
    double mean_ratio = 0.0;
    if (!run_ratios.empty()) {
      mean_ratio = std::accumulate(run_ratios.begin(), run_ratios.end(), 0.0) / run_ratios.size();
    }
    
    double rms_ratio = 0.0;
    if (!run_ratios.empty()) {
      double sum_sq_diff = 0.0;
      for (double ratio : run_ratios) {
        sum_sq_diff += (ratio - mean_ratio) * (ratio - mean_ratio);
      }
      rms_ratio = std::sqrt(sum_sq_diff / run_ratios.size());
    }
    
    double rms_per_sigma_threshold = mean_ratio + rms_ratio * 8.0;
    std::cout << "  Run " << run.dirName << ": Mean = " << mean_ratio << ", RMS = " << rms_ratio 
              << ", Threshold = " << rms_per_sigma_threshold << std::endl;
    
    // Second pass: apply flags based on this run's threshold
    for (auto &[cellid, entry] : run.byCellId) {
      if (entry.distance_rms_vals.empty() || entry.distance_sigma_vals.empty()) {
        continue;
      }
      
      double dist_rms = entry.distance_rms_vals[0];
      double dist_sigma = entry.distance_sigma_vals[0];
      
      if (dist_sigma > 0 && dist_rms >= 0) {
        double ratio = dist_rms / dist_sigma;
        if (ratio > rms_per_sigma_threshold) {
          // Add BAD_RMSPERSIGMA flag to the existing quality_flag
          if (!entry.flagsByRun.empty()) {
            entry.flagsByRun[0] |= HGLG_BAD_RMSPERSIGMA;
          }
        }
      }
    }
  }
}

bool loadQualityFlagRun(const std::string &baseDir, const std::string &dirName, QualityFlagRunData &out) {
  out = QualityFlagRunData();
  out.dirName = dirName;
  out.label = dirName;
  
  const std::string path = baseDir + "/" + dirName + "/hglg_calib_quality_muon.root";
  out.path = path;

  TFile *file = TFile::Open(out.path.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Failed to open quality file: " << out.path << std::endl;
    return false;
  }

  TTree *tree = dynamic_cast<TTree *>(file->Get("HGLGCalibQuality"));
  if (!tree) {
    std::cerr << "Missing 'HGLGCalibQuality' tree in: " << out.path << std::endl;
    file->Close();
    return false;
  }

  int p_cellid = -1, p_layer = -1, p_chip = -1, p_channel = -1;
  double p_x_mm = -999.0, p_y_mm = -999.0;
  uint32_t p_quality_flag = 0;
  double p_distance_rms = -1.0, p_distance_sigma = -1.0;
  int p_n_fit_points_over500 = -1;

  tree->SetBranchAddress("cellid", &p_cellid);
  tree->SetBranchAddress("layer", &p_layer);
  tree->SetBranchAddress("chip", &p_chip);
  tree->SetBranchAddress("channel", &p_channel);
  tree->SetBranchAddress("quality_flag", &p_quality_flag);

  bool hasCoordinates = false;
  if (tree->GetBranch("x_mm") && tree->GetBranch("y_mm")) {
    tree->SetBranchAddress("x_mm", &p_x_mm);
    tree->SetBranchAddress("y_mm", &p_y_mm);
    hasCoordinates = true;
  }

  bool hasDistanceRMS = false;
  bool hasDistanceSigma = false;
  if (tree->GetBranch("distance_rms")) {
    tree->SetBranchAddress("distance_rms", &p_distance_rms);
    hasDistanceRMS = true;
  }
  if (tree->GetBranch("distance_sigma")) {
    tree->SetBranchAddress("distance_sigma", &p_distance_sigma);
    hasDistanceSigma = true;
  }
  bool hasNFitPointsOver500 = false;
  if (tree->GetBranch("n_fit_points_over500")) {
    tree->SetBranchAddress("n_fit_points_over500", &p_n_fit_points_over500);
    hasNFitPointsOver500 = true;
  }
  if (!hasNFitPointsOver500) {
    std::cout << "Warning: n_fit_points_over500 branch not found in "
              << out.path << std::endl;
  }

  const Long64_t nEntries = tree->GetEntries();
  out.entries.reserve(nEntries);

  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    QualityFlagEntry e;
    e.cellid = p_cellid;
    e.layer = p_layer;
    e.chip = p_chip;
    e.channel = p_channel;
    e.fit_ok = 1;
    if (hasCoordinates) {
      e.x_mm = p_x_mm;
      e.y_mm = p_y_mm;
    }
    e.flagsByRun.push_back(p_quality_flag);
    if (hasDistanceRMS) {
      e.distance_rms_vals.push_back(p_distance_rms);
    }
    if (hasDistanceSigma) {
      e.distance_sigma_vals.push_back(p_distance_sigma);
    }
    if (hasNFitPointsOver500) {
      e.n_fit_points_over500 = p_n_fit_points_over500;
      // std::cout << "  CellID " << e.cellid << ": n_fit_points_over500 = " << e.n_fit_points_over500 << std::endl;
    } else {
      e.n_fit_points_over500 = -1; // Default value if not present
      // std::cout << "  CellID " << e.cellid << ": n_fit_points_over500 not present in tree, set to -1" << std::endl;
    }
    if (p_quality_flag == 1<<0 || p_quality_flag == 1<<1) {
      e.fit_ok = 0;
    }
    out.entries.push_back(e);
    out.byCellId[e.cellid] = e;
  }

  file->Close();
  return true;
}

void copyNFitPointsOver500FromQuality(RunData &run,
                                      const QualityFlagRunData &qualityRun) {
  for (auto &entry : run.entries) {
    auto qIt = qualityRun.byCellId.find(entry.cellid);
    if (qIt != qualityRun.byCellId.end()) {
      entry.n_fit_points_over500 = qIt->second.n_fit_points_over500;
    }
  }

  for (auto &[cellid, entry] : run.byCellId) {
    auto qIt = qualityRun.byCellId.find(cellid);
    if (qIt != qualityRun.byCellId.end()) {
      entry.n_fit_points_over500 = qIt->second.n_fit_points_over500;
    }
  }
}

uint32_t badQualityFlagMask() {
  return HGLG_LOW_STAT | HGLG_FIT_FAILED | HGLG_BAD_RMSPERSIGMA |
         HGLG_LOW_GAIN | HGLG_PEDESTAL_MASKED | HGLG_MANUAL_BAD;
}

bool hasBadQualityFlag(const QualityFlagEntry &entry) {
  return !entry.flagsByRun.empty() && (entry.flagsByRun[0] & badQualityFlagMask());
}

bool isBadInQualityRun(const QualityFlagRunData *qualityRun, int cellid) {
  if (!qualityRun) return false;
  auto it = qualityRun->byCellId.find(cellid);
  if (it == qualityRun->byCellId.end()) return false;
  return hasBadQualityFlag(it->second);
}

std::map<std::string, const QualityFlagRunData *> mapQualityRunsByDir(
    const std::vector<QualityFlagRunData> &runs) {
  std::map<std::string, const QualityFlagRunData *> byDir;
  for (const auto &run : runs) {
    byDir[run.dirName] = &run;
  }
  return byDir;
}

void drawSlopeErrorOverSlopeGoodOverlay(
    const std::vector<RunData> &runs,
    const std::vector<QualityFlagRunData> &qualityRuns,
    const std::string &outDir) {
  if (runs.empty() || qualityRuns.empty()) return;

  const auto qualityByDir = mapQualityRunsByDir(qualityRuns);
  std::vector<std::vector<double>> ratiosByRun(runs.size());
  double globalMin = 1e30;
  double globalMax = -1e30;

  for (size_t i = 0; i < runs.size(); ++i) {
    const QualityFlagRunData *qRun = nullptr;
    auto qIt = qualityByDir.find(runs[i].dirName);
    if (qIt != qualityByDir.end()) qRun = qIt->second;

    for (const auto &entry : runs[i].entries) {
      if (!entry.fit_ok) continue;
      if (isBadInQualityRun(qRun, entry.cellid)) continue;
      if (!std::isfinite(entry.slope) || !std::isfinite(entry.slope_error)) continue;
      if (std::abs(entry.slope) <= 1e-12 || entry.slope_error < 0.0) continue;

      const double ratio = entry.slope_error / entry.slope;
      if (!std::isfinite(ratio)) continue;

      ratiosByRun[i].push_back(ratio);
      globalMin = std::min(globalMin, ratio);
      globalMax = std::max(globalMax, ratio);
    }
  }

  if (globalMin > globalMax) {
    std::cerr << "No valid slope_error/slope entries found for non-BAD channels"
              << std::endl;
    return;
  }
  if (globalMin == globalMax) {
    globalMin -= 0.1 * std::abs(globalMin == 0.0 ? 1.0 : globalMin);
    globalMax += 0.1 * std::abs(globalMax == 0.0 ? 1.0 : globalMax);
  }

  const double span = std::max(1e-12, globalMax - globalMin);
  const double xMin = globalMin - 0.05 * span;
  const double xMax = globalMax + 0.05 * span;
  const int nbins = 100;
  int colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7,
                  kCyan + 2, kViolet + 2, kPink, kSpring, kTeal, kAzure, kRose};

  TCanvas *c = new TCanvas("c_slope_error_over_slope_good",
                           "slope_error/slope for non-BAD channels",
                           1200, 700);
  TLegend *leg = new TLegend(0.58, 0.52, 0.90, 0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(1);

  TH1D *hFirst = nullptr;
  for (size_t i = 0; i < runs.size(); ++i) {
    if (ratiosByRun[i].empty()) continue;

    TH1D *h = new TH1D(Form("h_slope_error_over_slope_good_run%zu", i),
                       "slope_error/slope for non-BAD channels;slope_error / slope;Channels",
                       nbins, xMin, xMax);
    h->SetLineColor(colors[i % 12]);
    h->SetLineWidth(2);
    h->SetFillStyle(0);

    for (double ratio : ratiosByRun[i]) {
      h->Fill(ratio);
    }

    SummaryStats stats = computeStats(ratiosByRun[i]);
    leg->AddEntry(
        h,
        Form("%s (N=%d, mean=%.3g, RMS=%.3g)",
             runs[i].label.c_str(), stats.n, stats.mean, stats.rms),
        "l");

    if (!hFirst) {
      hFirst = h;
      hFirst->Draw("hist");
    } else {
      h->Draw("hist same");
    }
  }

  if (!hFirst) {
    delete c;
    return;
  }

  leg->Draw();
  gPad->SetGridy();
  gStyle->SetOptStat(0);
  c->SaveAs((outDir + "/Distribution_slope_error_over_slope_good.pdf").c_str());
}

std::set<int> collectHighBadRateCellIds(const std::vector<QualityFlagRunData> &runs,
                                        double maxKeptBadRate = 0.50) {
  std::set<int> badCellIds;
  if (runs.empty()) return badCellIds;

  std::set<int> allCellIds;
  for (const auto &run : runs) {
    for (const auto &[cellid, entry] : run.byCellId) {
      allCellIds.insert(cellid);
    }
  }

  for (int cellid : allCellIds) {
    int badCount = 0;
    for (const auto &run : runs) {
      auto it = run.byCellId.find(cellid);
      if (it != run.byCellId.end() && hasBadQualityFlag(it->second)) {
        ++badCount;
      }
    }
    if (static_cast<double>(badCount) / runs.size() > maxKeptBadRate) {
        badCellIds.insert(cellid);
      }
  }
  return badCellIds;
}

std::vector<RunData> buildRunsWithLowRateBadImputed(
    const std::vector<RunData> &runs,
    const std::vector<QualityFlagRunData> &qualityRuns,
    const std::set<int> &excludedCellIds) {
  struct AverageValues {
    double slope = 0.0;
    double intercept = 0.0;
    double chi2_ndf = 0.0;
    double hg_adc_saturation = 0.0;
    double n_points = 0.0;
    double n_fit_points_over500 = 0.0;
    int count = 0;
  };

  const auto qualityByDir = mapQualityRunsByDir(qualityRuns);
  std::map<int, AverageValues> averages;

  for (const auto &run : runs) {
    const QualityFlagRunData *qRun = nullptr;
    auto qIt = qualityByDir.find(run.dirName);
    if (qIt != qualityByDir.end()) qRun = qIt->second;

    for (const auto &entry : run.entries) {
      if (excludedCellIds.count(entry.cellid)) continue;
      if (isBadInQualityRun(qRun, entry.cellid)) continue;
      auto &avg = averages[entry.cellid];
      avg.slope += entry.slope;
      avg.intercept += entry.intercept;
      avg.chi2_ndf += entry.chi2_ndf;
      avg.hg_adc_saturation += entry.hg_adc_saturation;
      avg.n_points += static_cast<double>(entry.n_points);
      avg.n_fit_points_over500 += static_cast<double>(entry.n_fit_points_over500);
      ++avg.count;
    }
  }

  std::vector<RunData> handledRuns;
  handledRuns.reserve(runs.size());
  for (const auto &run : runs) {
    const QualityFlagRunData *qRun = nullptr;
    auto qIt = qualityByDir.find(run.dirName);
    if (qIt != qualityByDir.end()) qRun = qIt->second;

    RunData handled = run;
    handled.entries.clear();
    handled.byCellId.clear();
    handled.entries.reserve(run.entries.size());

    for (auto entry : run.entries) {
      if (excludedCellIds.count(entry.cellid)) continue;

      if (isBadInQualityRun(qRun, entry.cellid)) {
        auto avgIt = averages.find(entry.cellid);
        if (avgIt != averages.end() && avgIt->second.count > 0) {
          const auto &avg = avgIt->second;
          const double n = static_cast<double>(avg.count);
          entry.slope = avg.slope / n;
          entry.intercept = avg.intercept / n;
          entry.chi2_ndf = avg.chi2_ndf / n;
          entry.hg_adc_saturation = avg.hg_adc_saturation / n;
          entry.n_points = static_cast<long long>(std::llround(avg.n_points / n));
          entry.n_fit_points_over500 =
              static_cast<int>(std::lround(avg.n_fit_points_over500 / n));
          entry.fit_ok = 1;
          entry.fit_status = 0;
        }
      }

      handled.entries.push_back(entry);
      handled.byCellId[entry.cellid] = entry;
    }

    handledRuns.push_back(handled);
  }

  return handledRuns;
}

std::vector<RunData> filterRunsExcludingCellIds(
    const std::vector<RunData> &runs,
    const std::set<int> &excludedCellIds) {
  std::vector<RunData> filteredRuns;
  filteredRuns.reserve(runs.size());

  for (const auto &run : runs) {
    RunData filtered = run;
    filtered.entries.clear();
    filtered.byCellId.clear();
    filtered.entries.reserve(run.entries.size());

    for (const auto &entry : run.entries) {
      if (excludedCellIds.count(entry.cellid)) continue;
      filtered.entries.push_back(entry);
      filtered.byCellId[entry.cellid] = entry;
    }

    filteredRuns.push_back(filtered);
  }

  return filteredRuns;
}

void writeBadExcludedCellList(const std::string &fileName,
                              const std::set<int> &badCellIds) {
  std::ofstream ofs(fileName.c_str());
  if (!ofs) return;

  ofs << "# CellIDs excluded from bad-excluded comparison plots\n";
  ofs << "# Excluded only when bad_rate > 50%; bad_rate <= 50% is imputed from good-run average\n";
  ofs << "# Bad bits: LOW_STAT, FIT_FAILED, BAD_RMSPERSIGMA, LOW_GAIN, PED_MASKED, MANUAL_BAD\n";
  ofs << "# cellid\n";
  for (int cellid : badCellIds) {
    ofs << cellid << "\n";
  }
}

void printQualityFlagTable(const std::vector<QualityFlagRunData> &runs, const std::string &outDir) {
  if (runs.empty()) return;

  std::ofstream ofs((outDir + "/QualityFlags_table.txt").c_str());
  if (!ofs) return;

  ofs << "===============================================================================\n";
  ofs << "QUALITY FLAG COMPARISON TABLE\n";
  ofs << "===============================================================================\n\n";
  ofs << "Flags: GOOD=0, LOW_STAT(b0), FIT_FAILED(b1), BAD_RMS(b2), LOW_GAIN(b3), PED_MASKED(b6), ";
  ofs << "MANUAL_BAD(b7), MANUAL_GOOD(b8)\n\n";

  // Collect all unique cellids
  std::set<int> all_cellids;
  for (const auto &run : runs) {
    for (const auto &[cid, entry] : run.byCellId) {
      all_cellids.insert(cid);
    }
  }

  // Print header
  ofs << std::left << std::setw(10) << "CellID";
  for (const auto &run : runs) {
    ofs << std::setw(20) << run.label;
  }
  ofs << "\n";
  ofs << std::string(10 + 20 * runs.size(), '-') << "\n";

  // Print data for each cellid
  for (int cid : all_cellids) {
    ofs << std::left << std::setw(10) << cid;
    for (const auto &run : runs) {
      auto it = run.byCellId.find(cid);
      if (it != run.byCellId.end() && !it->second.flagsByRun.empty()) {
        uint32_t flag = it->second.flagsByRun[0];
        ofs << std::left << std::setw(20) << flagToString(flag);
      } else {
        ofs << std::left << std::setw(20) << "N/A";
      }
    }
    ofs << "\n";
  }

  ofs << "\n===============================================================================\n";

  // Print detailed flag bit status
  ofs << "\nBIT-WISE FLAG STATUS:\n";
  ofs << "===============================================================================\n";

  std::vector<std::pair<uint32_t, std::string>> flags = {
      {HGLG_LOW_STAT, "LOW_STAT (Bit 0)"},
      {HGLG_FIT_FAILED, "FIT_FAILED (Bit 1)"},
      {HGLG_BAD_RMSPERSIGMA, "BAD_RMSPERSIGMA (Bit 2)"},
      {HGLG_LOW_GAIN, "LOW_GAIN (Bit 3)"},
      {HGLG_PEDESTAL_MASKED, "PEDESTAL_MASKED (Bit 6)"},
      {HGLG_MANUAL_BAD, "MANUAL_BAD (Bit 7)"},
      {HGLG_MANUAL_GOOD, "MANUAL_GOOD (Bit 8)"}
  };

  for (const auto &[flagbit, flagname] : flags) {
    ofs << "\n" << flagname << ":\n";
    ofs << std::string(80, '-') << "\n";
    ofs << std::left << std::setw(10) << "CellID";
    for (const auto &run : runs) {
      ofs << std::setw(15) << run.label;
    }
    ofs << "\n" << std::string(10 + 15 * runs.size(), '-') << "\n";

    for (int cid : all_cellids) {
      ofs << std::left << std::setw(10) << cid;
      for (const auto &run : runs) {
        auto it = run.byCellId.find(cid);
        if (it != run.byCellId.end() && !it->second.flagsByRun.empty()) {
          uint32_t flag = it->second.flagsByRun[0];
          bool set = (flag & flagbit) != 0;
          ofs << std::left << std::setw(15) << (set ? "SET" : "NOT_SET");
        } else {
          ofs << std::left << std::setw(15) << "N/A";
        }
      }
      ofs << "\n";
    }
  }

  // Print RMS/Sigma ratio info
  ofs << "\n===============================================================================\n";
  ofs << "\nRMS/SIGMA RATIO (for BAD_RMSPERSIGMA flag calculation):\n";
  ofs << "===============================================================================\n";
  ofs << std::left << std::setw(10) << "CellID";
  for (const auto &run : runs) {
    ofs << std::setw(20) << run.label;
  }
  ofs << "\n";
  ofs << std::string(10 + 20 * runs.size(), '-') << "\n";

  for (int cid : all_cellids) {
    ofs << std::left << std::setw(10) << cid;
    for (const auto &run : runs) {
      auto it = run.byCellId.find(cid);
      if (it != run.byCellId.end() && !it->second.distance_rms_vals.empty() && !it->second.distance_sigma_vals.empty()) {
        double dist_rms = it->second.distance_rms_vals[0];
        double dist_sigma = it->second.distance_sigma_vals[0];
        if (dist_sigma > 0 && dist_rms >= 0) {
          double ratio = dist_rms / dist_sigma;
          char buf[20];
          snprintf(buf, sizeof(buf), "%.3f", ratio);
          ofs << std::left << std::setw(20) << buf;
        } else {
          ofs << std::left << std::setw(20) << "invalid";
        }
      } else {
        ofs << std::left << std::setw(20) << "N/A";
      }
    }
    ofs << "\n";
  }

  ofs << "\n===============================================================================\n\n";
}

void printQualityFlagStatistics(const std::vector<QualityFlagRunData> &runs, const std::string &outDir) {
  if (runs.empty()) return;

  std::ofstream ofs((outDir + "/QualityFlags_statistics.txt").c_str());
  if (!ofs) return;

  ofs << "===============================================================================\n";
  ofs << "QUALITY FLAG STATISTICS (Count per Run)\n";
  ofs << "===============================================================================\n\n";

  std::vector<std::pair<uint32_t, std::string>> flags = {
      {HGLG_LOW_STAT, "LOW_STAT (Bit 0)"},
      {HGLG_FIT_FAILED, "FIT_FAILED (Bit 1)"},
      {HGLG_BAD_RMSPERSIGMA, "BAD_RMSPERSIGMA (Bit 2)"},
      {HGLG_LOW_GAIN, "LOW_GAIN (Bit 3)"},
      {HGLG_PEDESTAL_MASKED, "PEDESTAL_MASKED (Bit 6)"},
      {HGLG_MANUAL_BAD, "MANUAL_BAD (Bit 7)"},
      {HGLG_MANUAL_GOOD, "MANUAL_GOOD (Bit 8)"}
  };

  for (const auto &[flagbit, flagname] : flags) {
    ofs << "\n" << flagname << ":\n";
    ofs << std::string(100, '=') << "\n";
    ofs << std::left << std::setw(20) << "Run";
    ofs << std::setw(20) << "SET";
    ofs << std::setw(20) << "NOT_SET";
    ofs << std::setw(20) << "N/A";
    ofs << std::setw(20) << "Total\n";
    ofs << std::string(100, '-') << "\n";

    for (const auto &run : runs) {
      int count_set = 0;
      int count_not_set = 0;
      int count_na = 0;

      for (const auto &[cid, entry] : run.byCellId) {
        if (entry.flagsByRun.empty()) {
          count_na++;
        } else {
          uint32_t flag = entry.flagsByRun[0];
          if (flag & flagbit) {
            count_set++;
          } else {
            count_not_set++;
          }
        }
      }

      int total = count_set + count_not_set + count_na;
      ofs << std::left << std::setw(20) << run.label;
      ofs << std::setw(20) << count_set;
      ofs << std::setw(20) << count_not_set;
      ofs << std::setw(20) << count_na;
      ofs << std::setw(20) << total << "\n";
    }
  }

  ofs << "\n===============================================================================\n\n";
}

void drawQualityFlagTrends(const std::vector<QualityFlagRunData> &runs, const std::string &outDir) {
  if (runs.size() < 2) return;

  // Collect all unique cellids
  std::set<int> all_cellids;
  for (const auto &run : runs) {
    for (const auto &[cid, entry] : run.byCellId) {
      all_cellids.insert(cid);
    }
  }

  // Select representative cellids with most changes
  std::vector<std::pair<int, int>> cellid_changes;
  for (int cid : all_cellids) {
    std::set<uint32_t> unique_flags;
    for (const auto &run : runs) {
      auto it = run.byCellId.find(cid);
      if (it != run.byCellId.end() && !it->second.flagsByRun.empty()) {
        unique_flags.insert(it->second.flagsByRun[0]);
      }
    }
    int n_changes = unique_flags.size();
    cellid_changes.push_back({cid, n_changes});
  }

  // Sort by number of changes descending
  std::sort(cellid_changes.begin(), cellid_changes.end(),
            [](const auto &a, const auto &b) { return a.second > b.second; });

  // Draw trends for top changed cellids
  int n_to_plot = std::min(20, (int)cellid_changes.size());
  const int nCols = 4;
  const int nRows = (n_to_plot + nCols - 1) / nCols;

  TCanvas *c = new TCanvas("c_quality_flags_trends", "Quality Flag Trends", 1600, 400 * nRows);
  c->Divide(nCols, nRows);

  for (int i = 0; i < n_to_plot; ++i) {
    c->cd(i + 1);
    int cellid = cellid_changes[i].first;

    std::vector<double> run_idx;
    std::vector<double> flag_vals;
    int run_count = 0;

    for (const auto &run : runs) {
      auto it = run.byCellId.find(cellid);
      if (it != run.byCellId.end() && !it->second.flagsByRun.empty()) {
        run_idx.push_back(run_count);
        flag_vals.push_back((double)it->second.flagsByRun[0]);
      }
      run_count++;
    }

    if (!run_idx.empty()) {
      TGraph *g = new TGraph(run_idx.size(), run_idx.data(), flag_vals.data());
      g->SetMarkerStyle(20);
      g->SetMarkerSize(1.0);
      g->SetMarkerColor(kBlack);
      g->SetLineColor(kBlack);
      g->SetLineWidth(2);
      g->SetTitle(Form("CellID %d", cellid));
      g->GetXaxis()->SetTitle("Run Index");
      g->GetYaxis()->SetTitle("Quality Flag");
      g->Draw("ALP");
      gPad->SetGridy();
    }
  }

  c->SaveAs((outDir + "/QualityFlags_trends.pdf").c_str());
}

void drawQualityFlagHeatmap(const std::vector<QualityFlagRunData> &runs, const std::string &outDir) {
  if (runs.empty()) return;

  // Create 2D histogram using coordinate: 36*9*layer+36*chip+channel vs run
  int n_runs = runs.size();
  
  // Calculate coordinate range
  int max_coord = 0;
  
  for (const auto &run : runs) {
    for (const auto &[cid, entry] : run.byCellId) {
      int layer = entry.layer;
      int chip = entry.chip;
      int channel = entry.channel;
      
      // Skip if any coordinate is invalid
      if (layer < 0 || chip < 0 || channel < 0) continue;
      
      int coord = 36 * 9 * layer + 36 * chip + channel;
      if (coord > max_coord) max_coord = coord;
    }
  }

  TH2D *h = new TH2D("h_quality_flags_heatmap", 
                     "Quality Flag Heatmap (36*9*layer+36*chip+channel); Run; Coordinate",
                     n_runs, 0, n_runs, max_coord + 1, 0, max_coord + 1);
  h->SetStats(0);

  for (int col = 0; col < (int)runs.size(); ++col) {
    for (const auto &[cid, entry] : runs[col].byCellId) {
      int layer = entry.layer;
      int chip = entry.chip;
      int channel = entry.channel;
      
      if (layer < 0 || chip < 0 || channel < 0) continue;
      
      int coord = 36 * 9 * layer + 36 * chip + channel;
      if (!entry.flagsByRun.empty()) {
        double flag_val = (double)entry.flagsByRun[0];
        h->SetBinContent(col + 1, coord + 1, flag_val);
      }
    }
  }

  // Set x-axis labels (run names)
  for (int i = 0; i < (int)runs.size(); ++i) {
    h->GetXaxis()->SetBinLabel(i + 1, runs[i].label.c_str());
  }

  TCanvas *c = new TCanvas("c_quality_flags_heatmap", "Quality Flag Heatmap", 1200, 800);
  gStyle->SetOptStat(0);
  h->Draw("COLZ");
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.1);
  c->SaveAs((outDir + "/QualityFlags_heatmap.pdf").c_str());
}

void printBadRatePerChannel(const std::vector<QualityFlagRunData> &runs, const std::string &outDir) {
  if (runs.empty()) return;

  struct ChannelBadRate {
    int cellid = -1;
    int layer = -1;
    int chip = -1;
    int channel = -1;
    double x_mm = -999.0;
    double y_mm = -999.0;
    int bad_count = 0;
    int total_runs = 0;
    double bad_rate = 0.0;
  };

  // Collect all unique cellids
  std::set<int> all_cellids;
  for (const auto &run : runs) {
    for (const auto &[cid, entry] : run.byCellId) {
      all_cellids.insert(cid);
    }
  }

  // Process each flag type separately
  std::vector<std::pair<uint32_t, std::string>> flag_types = {
      {HGLG_BAD_RMSPERSIGMA, "BAD_RMSPERSIGMA"},
      {HGLG_LOW_GAIN, "LOW_GAIN"},
      {HGLG_FIT_FAILED, "FIT_FAILED"},
      {HGLG_LOW_STAT, "LOW_STAT"}
  };

  for (const auto &[flagbit, flagname] : flag_types) {
    std::vector<ChannelBadRate> bad_rates;

    for (int cid : all_cellids) {
      ChannelBadRate cbr;
      cbr.cellid = cid;
      cbr.total_runs = runs.size();
      cbr.bad_count = 0;

      for (const auto &run : runs) {
        auto it = run.byCellId.find(cid);
        if (it != run.byCellId.end()) {
          if (cbr.layer < 0) {
            cbr.layer = it->second.layer;
            cbr.chip = it->second.chip;
            cbr.channel = it->second.channel;
            cbr.x_mm = it->second.x_mm;
            cbr.y_mm = it->second.y_mm;
          }
          if (!it->second.flagsByRun.empty() && (it->second.flagsByRun[0] & flagbit)) {
            cbr.bad_count++;
          }
        }
      }

      if (cbr.bad_count > 0) {
        cbr.bad_rate = 100.0 * cbr.bad_count / cbr.total_runs;
        bad_rates.push_back(cbr);
      }
    }

    // Sort by bad_rate descending
    std::sort(bad_rates.begin(), bad_rates.end(),
              [](const ChannelBadRate &a, const ChannelBadRate &b) {
                if (a.bad_rate != b.bad_rate) return a.bad_rate > b.bad_rate;
                return a.cellid < b.cellid;
              });

    // Write output file for this flag type
    std::ofstream ofs((outDir + "/BadRate_" + flagname + ".txt").c_str());
    if (!ofs) continue;

    ofs << "===============================================================================\n";
    ofs << "BAD RATE PER CHANNEL (" << flagname << ")\n";
    ofs << "===============================================================================\n\n";
    ofs << "Total Runs: " << runs.size() << "\n";
    ofs << "Channels with at least one " << flagname << " flag: " << bad_rates.size() << "\n\n";

    ofs << std::left << std::setw(10) << "CellID";
    ofs << std::setw(8) << "Layer";
    ofs << std::setw(8) << "Chip";
    ofs << std::setw(10) << "Channel";
    ofs << std::setw(10) << "X (mm)";
    ofs << std::setw(10) << "Y (mm)";
    ofs << std::setw(10) << "Bad Count";
    ofs << std::setw(10) << "Total";
    ofs << std::setw(12) << "Bad Rate\n";
    ofs << std::string(108, '-') << "\n";

    for (const auto &cbr : bad_rates) {
      ofs << std::left << std::setw(10) << cbr.cellid;
      ofs << std::setw(8) << cbr.layer;
      ofs << std::setw(8) << cbr.chip;
      ofs << std::setw(10) << cbr.channel;
      char buf_x[15], buf_y[15], buf_rate[15];
      snprintf(buf_x, sizeof(buf_x), "%.2f", cbr.x_mm);
      snprintf(buf_y, sizeof(buf_y), "%.2f", cbr.y_mm);
      snprintf(buf_rate, sizeof(buf_rate), "%.1f%%", cbr.bad_rate);
      ofs << std::left << std::setw(10) << buf_x;
      ofs << std::left << std::setw(10) << buf_y;
      ofs << std::left << std::setw(10) << cbr.bad_count;
      ofs << std::left << std::setw(10) << cbr.total_runs;
      ofs << std::left << std::setw(12) << buf_rate << "\n";
    }

    ofs << "\n===============================================================================\n";

    // Summary statistics
    ofs << "\nSUMMARY STATISTICS:\n";
    ofs << "===============================================================================\n";

    int count_100 = 0, count_50_100 = 0, count_1_50 = 0;
    for (const auto &cbr : bad_rates) {
      if (cbr.bad_rate == 100.0) count_100++;
      else if (cbr.bad_rate >= 50.0) count_50_100++;
      else if (cbr.bad_rate > 0) count_1_50++;
    }

    ofs << "Channels always bad (100%): " << count_100 << "\n";
    ofs << "Channels mostly bad (50-100%): " << count_50_100 << "\n";
    ofs << "Channels sometimes bad (1-50%): " << count_1_50 << "\n";
    ofs << "Total problematic channels: " << bad_rates.size() << "\n";

    ofs << "\n===============================================================================\n\n";
  }
}

std::vector<CellShiftInfo> collectFlaggedCells(
    const std::vector<QualityFlagRunData> &runs,
    uint32_t flagbit,
    bool mask_L21C1 = true) {
  std::vector<CellShiftInfo> flaggedCells;
  if (runs.empty()) return flaggedCells;

  std::set<int> allCellIds;
  for (const auto &run : runs) {
    for (const auto &[cid, entry] : run.byCellId) {
      allCellIds.insert(cid);
    }
  }

  for (int cid : allCellIds) {
    int badCount = 0;
    const QualityFlagEntry *firstEntry = nullptr;
    for (const auto &run : runs) {
      auto it = run.byCellId.find(cid);
      if (it == run.byCellId.end()) continue;
      if (!firstEntry) firstEntry = &it->second;
      if (!it->second.flagsByRun.empty() &&
          (it->second.flagsByRun[0] & flagbit)) {
        ++badCount;
      }
    }

    if (badCount == 0 || !firstEntry) continue;
    if (mask_L21C1 && firstEntry->layer == 21 && firstEntry->chip == 1) {
      continue;
    }

    CellShiftInfo info;
    info.cellid = cid;
    info.spread = static_cast<double>(badCount);
    info.refValue = 100.0 * badCount / runs.size();
    info.meanValue = static_cast<double>(runs.size());
    info.x_mm = firstEntry->x_mm;
    info.y_mm = firstEntry->y_mm;
    flaggedCells.push_back(info);
  }

  std::sort(flaggedCells.begin(), flaggedCells.end(),
            [](const CellShiftInfo &a, const CellShiftInfo &b) {
              if (a.spread != b.spread) return a.spread > b.spread;
              return a.cellid < b.cellid;
            });

  return flaggedCells;
}

std::vector<CellShiftInfo> collectBadRmsPerSigmaCells(
    const std::vector<QualityFlagRunData> &runs,
    bool mask_L21C1 = true) {
  return collectFlaggedCells(runs, HGLG_BAD_RMSPERSIGMA, mask_L21C1);
}

void writeFlaggedAccumulatedSummary(
    const std::string &fileName,
    const std::vector<CellShiftInfo> &badCells,
    const std::string &flagName) {
  std::ofstream ofs(fileName.c_str());
  if (!ofs) return;

  ofs << "# " << flagName << " channels used for accumulated 2D plots\n";
  ofs << "# Layer21Chip1 is excluded\n";
  ofs << "# rank cellid bad_run_count bad_rate_percent total_runs x_mm y_mm\n";
  for (size_t i = 0; i < badCells.size(); ++i) {
    const auto &info = badCells[i];
    ofs << i + 1 << " " << info.cellid << " "
        << static_cast<int>(info.spread) << " "
        << info.refValue << " " << static_cast<int>(info.meanValue) << " "
        << info.x_mm << " " << info.y_mm << "\n";
  }
}

void writeBadRmsPerSigmaAccumulatedSummary(
    const std::string &fileName,
    const std::vector<CellShiftInfo> &badCells) {
  writeFlaggedAccumulatedSummary(fileName, badCells, "BAD_RMSPERSIGMA");
}

void drawDistanceRmsVsSigmaDistributions(const std::vector<QualityFlagRunData> &runs, const std::string &outDir) {
  if (runs.empty()) return;

  // Collect min/max of RMS/Sigma ratio across all runs for common range
  double global_min = 1e10, global_max = -1e10;

  for (const auto &run : runs) {
    for (const auto &[cid, entry] : run.byCellId) {
      // RMS and Sigma must be paired by index
      size_t min_size = std::min(entry.distance_rms_vals.size(), entry.distance_sigma_vals.size());
      if (int(cid/100000)==21 && int(cid/10000)%10==1) {
        continue;
      }
      if (int(cid/100000)==13 && int(cid/10000)%10==5) {
        continue;
      }
      for (size_t i = 0; i < min_size; ++i) {
        double rms = entry.distance_rms_vals[i];
        double sigma = entry.distance_sigma_vals[i];
        if (rms >= 0 && sigma > 0) {
          double ratio = rms / sigma;
          global_min = std::min(global_min, ratio);
          global_max = std::max(global_max, ratio);
        }
      }
    }
  }

  if (global_min > global_max) {
    std::cerr << "No valid RMS/Sigma ratio data found" << std::endl;
    return;
  }

  if (global_min == global_max) {
    global_min -= 0.1;
    global_max += 0.1;
  }

  // Color palette for different runs
  std::vector<int> colors = {kBlue, kRed, kGreen, kMagenta, kCyan, kYellow+2, kOrange, 
                             kGray+1, kViolet, kPink, kSpring, kTeal, kAzure, kRose};

  TCanvas *c = new TCanvas("c_distance_distributions", "RMS/Sigma Ratio Distributions per Run", 1000, 700);
  // c->SetGridy();

  // Create legend
  TLegend *leg = new TLegend(0.65, 0.2, 0.95, 0.95);
  leg->SetFillStyle(0);
  leg->SetBorderSize(1);
  // Draw histograms per run with "same" option
  bool first_draw = true;
  for (size_t run_idx = 0; run_idx < runs.size(); ++run_idx) {
    const auto &run = runs[run_idx];

    // Collect RMS/Sigma ratio values for this run
    std::vector<double> ratio_vals;
    for (const auto &[cid, entry] : run.byCellId) {
      size_t min_size = std::min(entry.distance_rms_vals.size(), entry.distance_sigma_vals.size());
      if (int(cid/100000)==21 && int(cid/10000)%10==1) {
        continue;
      }
      if (int(cid/100000)==13 && int(cid/10000)%10==5) {
        continue;
      }
      for (size_t i = 0; i < min_size; ++i) {
        double rms = entry.distance_rms_vals[i];
        double sigma = entry.distance_sigma_vals[i];
        if (rms >= 0 && sigma > 0) {
          ratio_vals.push_back(rms / sigma);
        }
      }
    }

    if (ratio_vals.empty()) continue;

    // Create histogram for this run
    TH1D *h = new TH1D(Form("h_ratio_run%zu", run_idx), 
                       Form("RMS/Sigma Ratio - Run %s; Ratio; Count", run.dirName.c_str()),
                       100, global_min, global_max);
    h->SetLineColor(colors[run_idx % colors.size()]);
    h->SetLineWidth(2);
    h->SetFillStyle(0);  // No fill for overlaid

    for (double val : ratio_vals) {
      h->Fill(val);
    }

    // Draw with "same" option for subsequent runs
    if (first_draw) {
      h->Draw();
      first_draw = false;
    } else {
      h->Draw("same");
    }
    double avg = std::accumulate(ratio_vals.begin(), ratio_vals.end(), 0.0) / ratio_vals.size();
    double rms = 0.0;
    for (double val : ratio_vals) {      
        rms += (val - avg) * (val - avg);
    }
    rms = std::sqrt(rms / ratio_vals.size());
    
    leg->AddEntry(h, Form("%s (N=%zu) Avg-1: %.4f RMS: %.4f", run.dirName.c_str(), ratio_vals.size(), avg-1, rms), "l");
  }

  leg->Draw();
  c->SetLogy();  // Optional: set log scale for better visibility
  c->SaveAs((outDir + "/Distance_RMS_vs_Sigma_distributions.pdf").c_str());
}

// ============================================================================
// Main Quality Flag Comparison Function
// ============================================================================

void compareQualityFlags(const char *baseDir = ".",
                         const char *runList = "all",
                         const char *outDir = "quality_flags_out") {
  // Create output directory
  gSystem->Exec(Form("mkdir -p %s", outDir));

  // Parse directory ranges
  std::vector<std::string> dirNames = parseDirectoryRanges(runList);
  if (dirNames.empty()) {
    std::cerr << "No directories to process" << std::endl;
    return;
  }

  std::cout << "Loading quality flag data from " << dirNames.size() << " directories..." << std::endl;

  // Load all quality flag data
  std::vector<QualityFlagRunData> allQualityFlags;
  for (const auto &dirName : dirNames) {
    QualityFlagRunData qfdata;
    if (loadQualityFlagRun(baseDir, dirName, qfdata)) {
      std::cout << "  " << dirName << ": " << qfdata.entries.size() << " entries" << std::endl;
      allQualityFlags.push_back(qfdata);
    } else {
      std::cerr << "  Warning: Failed to load quality file from " << dirName << std::endl;
    }
  }

  if (allQualityFlags.empty()) {
    std::cerr << "No quality flag data could be loaded" << std::endl;
    return;
  }

  std::cout << "Successfully loaded " << allQualityFlags.size() << " quality datasets\n" << std::endl;

  // Add RMS/Sigma flags based on distance metrics
  std::cout << "Computing BAD_RMSPERSIGMA flags..." << std::endl;
  addRMSPerSigmaFlags(allQualityFlags);
  std::cout << "  Done\n" << std::endl;

  std::cout << "Computing LOW_GAIN flags..." << std::endl;
  addLowGainFlags(allQualityFlags, 50);
  std::cout << "  Done\n" << std::endl;

  // Generate outputs
  std::cout << "Generating quality flag comparison outputs..." << std::endl;

  printQualityFlagTable(allQualityFlags, outDir);
  std::cout << "  Printed quality flag table" << std::endl;

  printQualityFlagStatistics(allQualityFlags, outDir);
  std::cout << "  Printed quality flag statistics" << std::endl;

  printBadRatePerChannel(allQualityFlags, outDir);
  std::cout << "  Printed bad rate per channel" << std::endl;

  drawQualityFlagTrends(allQualityFlags, outDir);
  std::cout << "  Drew quality flag trends" << std::endl;

  drawQualityFlagHeatmap(allQualityFlags, outDir);
  std::cout << "  Drew quality flag heatmap" << std::endl;

  drawDistanceRmsVsSigmaDistributions(allQualityFlags, outDir);
  std::cout << "  Drew distance RMS vs Sigma distributions" << std::endl;

  auto badRmsCells = collectBadRmsPerSigmaCells(allQualityFlags, true);
  writeBadRmsPerSigmaAccumulatedSummary(
      std::string(outDir) + "/BAD_RMSPERSIGMA_accumulated_channels.txt",
      badRmsCells);
  std::cout << "  Found " << badRmsCells.size()
            << " BAD_RMSPERSIGMA channels for accumulated 2D plots"
            << " (excluding Layer21Chip1)" << std::endl;

  auto badLowGainCells = collectFlaggedCells(allQualityFlags, HGLG_LOW_GAIN, true);
  writeFlaggedAccumulatedSummary(
      std::string(outDir) + "/BAD_LOWGAIN_accumulated_channels.txt",
      badLowGainCells, "BAD_LOWGAIN");
  std::cout << "  Found " << badLowGainCells.size()
            << " BAD_LOWGAIN channels for accumulated 2D plots"
            << " (excluding Layer21Chip1)" << std::endl;

  if (kDrawLegacyTH2DPlots && (!badRmsCells.empty() || !badLowGainCells.empty())) {
    std::vector<RunData> interCalibRuns;
    interCalibRuns.reserve(allQualityFlags.size());
    for (const auto &qRun : allQualityFlags) {
      RunData data;
      if (loadRun(baseDir, qRun.dirName, data)) {
        copyNFitPointsOver500FromQuality(data, qRun);
        interCalibRuns.push_back(data);
      } else {
        std::cerr << "  Warning: Failed to load intercalib file for 2D plots from "
                  << qRun.dirName << std::endl;
      }
    }

    if (!interCalibRuns.empty() && !badRmsCells.empty()) {
      drawShiftedChannelAccumulatedHistograms(
          interCalibRuns, badRmsCells, outDir, badRmsCells.size(),
          "BAD_RMSPERSIGMA_accumulated",
          "BAD_RMSPERSIGMA run count");
      std::cout << "  Drew BAD_RMSPERSIGMA accumulated 2D plots" << std::endl;
    }
    if (!interCalibRuns.empty() && !badLowGainCells.empty()) {
      drawShiftedChannelAccumulatedHistograms(
          interCalibRuns, badLowGainCells, outDir, badLowGainCells.size(),
          "BAD_LOWGAIN_accumulated",
          "BAD_LOWGAIN run count");
      std::cout << "  Drew BAD_LOWGAIN accumulated 2D plots" << std::endl;
    }
  }

  std::cout << "\nQuality flag comparison complete. Output in: " << outDir << std::endl;
}

// ============================================================================
// Main Comparison Function
// ============================================================================

void compare(const char *baseDir = ".",
             const char *runList = "all",
             const char *outDir = "comparison_out_directory") {
  // Create output directory
  gSystem->Exec(Form("mkdir -p %s", outDir));

  // Parse directory ranges
  std::vector<std::string> dirNames = parseDirectoryRanges(runList);
  if (dirNames.empty()) {
    std::cerr << "No directories to process" << std::endl;
    return;
  }

  std::cout << "Loading " << dirNames.size() << " directories..." << std::endl;

  // Load all directories
  std::vector<RunData> allRuns;
  std::vector<QualityFlagRunData> qualityRunsForCompare;
  for (const auto &dirName : dirNames) {
    if (gSystem->AccessPathName((std::string(baseDir) + "/" + dirName+"/hglg_calib_quality_muon.root").c_str()) != 0) {
      std::cerr << "  Warning: Directory not accessible: " << dirName << std::endl;
      continue;
    }
    RunData data;
    if (loadRun(baseDir, dirName, data)) {
      QualityFlagRunData qfdata;
      if (loadQualityFlagRun(baseDir, dirName, qfdata)) {
        copyNFitPointsOver500FromQuality(data, qfdata);
        qualityRunsForCompare.push_back(qfdata);
      } else {
        std::cerr << "  Warning: Failed to load quality file for n_fit_points_over500 from "
                  << dirName << std::endl;
      }
      std::cout << "  " << dirName << ": " << data.entries.size() << " entries" << std::endl;
      allRuns.push_back(data);
    } else {
      std::cerr << "  Failed to load directory " << dirName << std::endl;
    }
  }

  if (allRuns.size() < 2) {
    std::cerr << "Need at least 2 runs for comparison" << std::endl;
    return;
  }

  std::cout << "Successfully loaded " << allRuns.size() << " runs\n" << std::endl;

  // Print summary
  printSummaryTable(allRuns);

  std::vector<RunData> badExcludedRuns;
  std::set<int> badExcludedCellIds;
  if (!qualityRunsForCompare.empty()) {
    std::cout << "Computing BAD flags for bad-excluded comparison plots..." << std::endl;
    addRMSPerSigmaFlags(qualityRunsForCompare);
    addLowGainFlags(qualityRunsForCompare, 50);
    badExcludedCellIds = collectHighBadRateCellIds(qualityRunsForCompare, 0.50);
    writeBadExcludedCellList(
        std::string(outDir) + "/bad_excluded_cellids.txt",
        badExcludedCellIds);
    badExcludedRuns = buildRunsWithLowRateBadImputed(
        allRuns, qualityRunsForCompare, badExcludedCellIds);
    std::cout << "  Excluding " << badExcludedCellIds.size()
              << " channels with bad_rate > 50%; bad_rate <= 50% is imputed"
              << " from good-run average" << std::endl;
    drawSlopeErrorOverSlopeGoodOverlay(
        allRuns, qualityRunsForCompare, std::string(outDir));
    std::cout << "  Drew slope_error/slope overlay for non-BAD channels"
              << std::endl;
  } else {
    std::cout << "Warning: no quality flag data loaded; bad-excluded plots will be skipped"
              << std::endl;
  }

  // Compare metrics
  std::cout << "Generating comparison plots..." << std::endl;

  std::vector<std::pair<std::string, std::string>> metrics = {
      {"slope", "Slope (p1)"},
      {"intercept", "Intercept (p0)"},
      {"chi2_ndf", "Chi2/ndf"},
      {"n_points", "N Points"},
      {"hg_adc_saturation", "HG ADC Saturation"}};

  for (const auto &[branch, label] : metrics) {
    // Delta distributions
    drawDeltaFromReferenceByRun(allRuns, branch, label, "all_runs", outDir, true);

    // Representative channels
    auto shifted = rankCellsBySpread(allRuns, branch, true);
    if (branch == "slope" && !shifted.empty()) {
      drawNPointsSlopeCorrelations(allRuns, shifted, outDir, 30);
    }
    if (!shifted.empty()) {
      const bool isSlope = branch == "slope";
      drawRepresentativeTrendSet(allRuns, shifted, branch, label, "shifted channels",
                                 outDir, "shifted", true, isSlope ? 30 : 12, isSlope);
      writeRepresentativeSummary(
          std::string(outDir) + "/shifted_" + branch + ".txt", shifted, shifted, branch);
      if (isSlope) {
        writeTopSlopeTrendSummary(
            std::string(outDir) + "/top30_slope_trends.txt", allRuns, shifted, 30);
        drawSlopeRmsDistribution(shifted, outDir);
        if (kDrawLegacyTH2DPlots) {
          drawShiftedChannelAccumulatedHistograms(allRuns, shifted, outDir, 30);
        }

        std::vector<CellShiftInfo> shiftedNoLayer3;
        shiftedNoLayer3.reserve(shifted.size());
        for (const auto &info : shifted) {
          auto it = allRuns.front().byCellId.find(info.cellid);
          if (it != allRuns.front().byCellId.end() && it->second.layer != 3) {
            shiftedNoLayer3.push_back(info);
          }
        }
        if (!shiftedNoLayer3.empty()) {
          drawRepresentativeTrendSet(allRuns, shiftedNoLayer3, branch, label,
                                     "shifted channels excluding Layer 3",
                                     outDir, "shifted_no_layer3", true, 30, true);
          writeTopSlopeTrendSummary(
              std::string(outDir) + "/top30_slope_trends_no_layer3.txt",
              allRuns, shiftedNoLayer3, 30);
          drawSlopeRmsDistribution(shiftedNoLayer3, outDir, "no_layer3");
          if (kDrawLegacyTH2DPlots) {
            drawShiftedChannelAccumulatedHistograms(
                allRuns, shiftedNoLayer3, outDir, 30,
                "shifted_channel_accumulated_no_layer3");
          }
          if (branch == "slope" && !shiftedNoLayer3.empty()) {
            drawNPointsSlopeCorrelations(allRuns, shiftedNoLayer3, outDir, 30, "no_layer3");
          }
        }

        const auto shiftedFitRangeNPointsCut =
            excludeLowPointsAtShortFitRange(allRuns, shifted, 2000.0, 5000);
        if (!shiftedFitRangeNPointsCut.empty()) {
          drawRepresentativeTrendSet(
              allRuns, shiftedFitRangeNPointsCut, branch, label,
              "shifted channels after fitRange/n_points cut",
              outDir, "shifted_fitrange_npoints_cut", true, 30, true);
          writeTopSlopeTrendSummary(
              std::string(outDir) + "/top30_slope_trends_fitrange_npoints_cut.txt",
              allRuns, shiftedFitRangeNPointsCut, 30);
          drawSlopeRmsDistribution(
              shiftedFitRangeNPointsCut, outDir, "fitrange_npoints_cut");
          if (kDrawLegacyTH2DPlots) {
            drawShiftedChannelAccumulatedHistograms(
                allRuns, shiftedFitRangeNPointsCut, outDir, 30,
                "shifted_channel_accumulated_fitrange_npoints_cut");
          }
        }

        std::vector<CellShiftInfo> shiftedFitRangeNPointsCutNoLayer3;
        shiftedFitRangeNPointsCutNoLayer3.reserve(shiftedFitRangeNPointsCut.size());
        for (const auto &info : shiftedFitRangeNPointsCut) {
          auto it = allRuns.front().byCellId.find(info.cellid);
          if (it != allRuns.front().byCellId.end() && it->second.layer != 3) {
            shiftedFitRangeNPointsCutNoLayer3.push_back(info);
          }
        }
        if (!shiftedFitRangeNPointsCutNoLayer3.empty()) {
          drawRepresentativeTrendSet(
              allRuns, shiftedFitRangeNPointsCutNoLayer3, branch, label,
              "shifted channels after fitRange/n_points cut, excluding Layer 3",
              outDir, "shifted_fitrange_npoints_cut_no_layer3", true, 30, true);
          writeTopSlopeTrendSummary(
              std::string(outDir) +
                  "/top30_slope_trends_fitrange_npoints_cut_no_layer3.txt",
              allRuns, shiftedFitRangeNPointsCutNoLayer3, 30);
          drawSlopeRmsDistribution(
              shiftedFitRangeNPointsCutNoLayer3, outDir,
              "fitrange_npoints_cut_no_layer3");
          if (kDrawLegacyTH2DPlots) {
            drawShiftedChannelAccumulatedHistograms(
                allRuns, shiftedFitRangeNPointsCutNoLayer3, outDir, 30,
                "shifted_channel_accumulated_fitrange_npoints_cut_no_layer3");
          }
        }
      }
    }

    if (branch == "hg_adc_saturation") {
      auto shiftedHighNPoints = rankCellsBySpread(allRuns, branch, true, 10000);
      if (!shiftedHighNPoints.empty()) {
        drawRepresentativeTrendSet(
            allRuns, shiftedHighNPoints, branch, label,
            "shifted channels with n_points > 10000",
            outDir, "shifted_npoints_gt10000", true, 12, true, 10000);
        writeRepresentativeSummary(
            std::string(outDir) + "/shifted_hg_adc_saturation_npoints_gt10000.txt",
            shiftedHighNPoints, shiftedHighNPoints, branch);
      }
    }

    if (badExcludedRuns.size() >= 2 && !qualityRunsForCompare.empty()) {
      drawDeltaFromReferenceByRun(
          badExcludedRuns, branch, label, "bad_excluded", outDir, true);

      auto shiftedBadExcluded = rankCellsBySpread(badExcludedRuns, branch, true);
      if (!shiftedBadExcluded.empty()) {
        const bool isSlope = branch == "slope";
        drawRepresentativeTrendSet(
            badExcludedRuns, shiftedBadExcluded, branch, label,
            "shifted channels excluding BAD channels",
            outDir, "shifted_bad_excluded", true, isSlope ? 30 : 12, isSlope);
        writeRepresentativeSummary(
            std::string(outDir) + "/shifted_bad_excluded_" + branch + ".txt",
            shiftedBadExcluded, shiftedBadExcluded, branch);

        if (isSlope) {
          writeTopSlopeTrendSummary(
              std::string(outDir) + "/top30_slope_trends_bad_excluded.txt",
              badExcludedRuns, shiftedBadExcluded, 30);
          drawSlopeRmsDistribution(shiftedBadExcluded, outDir, "bad_excluded");
          if (kDrawLegacyTH2DPlots) {
            drawShiftedChannelAccumulatedHistograms(
                badExcludedRuns, shiftedBadExcluded, outDir, 30,
                "shifted_channel_accumulated_bad_excluded");
          }
          drawNPointsSlopeCorrelations(
              badExcludedRuns, shiftedBadExcluded, outDir, 30,
              "slope_vs_n_points_shifted_channels_bad_excluded");
        }
      }

      if (branch == "hg_adc_saturation") {
        auto shiftedBadExcludedHighNPoints =
            rankCellsBySpread(badExcludedRuns, branch, true, 10000);
        if (!shiftedBadExcludedHighNPoints.empty()) {
          drawRepresentativeTrendSet(
              badExcludedRuns, shiftedBadExcludedHighNPoints, branch, label,
              "shifted channels excluding BAD channels with n_points > 10000",
              outDir, "shifted_bad_excluded_npoints_gt10000",
              true, 12, true, 10000);
          writeRepresentativeSummary(
              std::string(outDir) +
                  "/shifted_bad_excluded_hg_adc_saturation_npoints_gt10000.txt",
              shiftedBadExcludedHighNPoints, shiftedBadExcludedHighNPoints, branch);
        }
      }
    }
  }

  // Summary graphs
  drawSummaryGraphs(allRuns, outDir);

  // Parameter distributions
  drawParameterDistributions(allRuns, outDir);

  // Correlation between statistics across runs
  drawMeanNPointsVsHgSaturationRms(allRuns, outDir, true);

  const std::string hgStabilityTag = "excludeL1L3L13C5_meanNPoints_ge10000";
  auto hgStabilityCellIds =
      collectChannelsForHgSaturationStability(allRuns, 10000.0, true);
  std::ofstream hgStabilityList(
      (std::string(outDir) + "/HGADCSaturation_stability_" +
       hgStabilityTag + "_channels.txt").c_str());
  if (hgStabilityList) {
    hgStabilityList << "# Channels used for HG ADC saturation stability plots\n";
    hgStabilityList << "# Selection: layer != 1, layer != 3, not L13C5, mean n_points >= 10000\n";
    hgStabilityList << "# cellid\n";
    for (int cellid : hgStabilityCellIds) {
      hgStabilityList << cellid << "\n";
    }
  }
  if (!hgStabilityCellIds.empty()) {
    auto hgStabilityRuns = filterRunsKeepingCellIds(allRuns, hgStabilityCellIds);
    drawHgSaturationDistributionByRun(hgStabilityRuns, outDir, hgStabilityTag);
    drawDeltaFromReferenceByRun(
        hgStabilityRuns, "hg_adc_saturation", "HG ADC Saturation",
        hgStabilityTag, outDir, true, false, false);
    drawHgSaturationRmsDistribution(hgStabilityRuns, outDir, hgStabilityTag);
    std::cout << "  Drew HG saturation stability plots for "
              << hgStabilityCellIds.size()
              << " channels (excluding L1/L3/L13C5 and mean n_points < 10000)"
              << std::endl;
  }

  auto highMeanNPointsHighHgRmsCells =
      collectHighMeanNPointsHighHgSaturationRmsCells(allRuns, 8000.0, 70.0, true);
  writeHighMeanNPointsHighHgSaturationRmsSummary(
      std::string(outDir) +
          "/MeanNPoints_gt8000_HGADCSaturationRMS_gt70_channels.txt",
      highMeanNPointsHighHgRmsCells);
  if (!highMeanNPointsHighHgRmsCells.empty()) {
    drawShiftedChannelAccumulatedHistograms(
        allRuns, highMeanNPointsHighHgRmsCells, outDir,
        highMeanNPointsHighHgRmsCells.size(),
        "MeanNPoints_gt8000_HGADCSaturationRMS_gt70_accumulated",
        "HG ADC saturation RMS");
    std::cout << "  Drew accumulated 2D plots for "
              << highMeanNPointsHighHgRmsCells.size()
              << " channels with mean n_points > 8000 and HG saturation RMS > 70"
              << std::endl;
  }

  // Run-by-run accumulated 2D plots for all channels in Layer 13, Chip 5
  auto layer13Chip5Cells = collectChannelsByLayerChip(allRuns, 13, 5);
  if (kDrawLegacyTH2DPlots && !layer13Chip5Cells.empty()) {
    drawShiftedChannelAccumulatedHistograms(
        allRuns, layer13Chip5Cells, outDir, layer13Chip5Cells.size(),
        "L13C5_accumulated_all_channels",
        "Channel");
    std::cout << "  Drew Layer 13 Chip 5 accumulated 2D plots for "
              << layer13Chip5Cells.size() << " channels" << std::endl;
  }

  std::cout << "\nComparison complete. Output in: " << outDir << std::endl;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " [compare|quality_flags|both] [options]\n";
    std::cout << "Options:\n";
    std::cout << "  --baseDir <path>       Base directory containing run subdirectories (default: .)\n";
    std::cout << "  --runList <list>      Comma-separated list of runs or ranges (default: all)\n";
    std::cout << "  --outDir <path>       Output directory for results (default: comparison_out_directory)\n";
    return 1;
  }

  std::string mode = argv[1];
  std::string baseDir = ".";
  std::string runList = "all";
  std::string outDir = "comparison_out_directory";

  for (int i = 2; i < argc; ++i) {
    if (std::string(argv[i]) == "--baseDir" && i + 1 < argc) {
      baseDir = argv[++i];
    } else if (std::string(argv[i]) == "--runList" && i + 1 < argc) {
      runList = argv[++i];
    } else if (std::string(argv[i]) == "--outDir" && i + 1 < argc) {
      outDir = argv[++i];
    }
  }

  if (mode == "compare") {
    compare(baseDir.c_str(), runList.c_str(), outDir.c_str());
  } else if (mode == "quality_flags") {
    compareQualityFlags(baseDir.c_str(), runList.c_str(), outDir.c_str());
  } else if (mode == "both") {
    compare(baseDir.c_str(), runList.c_str(), outDir.c_str());
    compareQualityFlags(baseDir.c_str(), runList.c_str(), outDir.c_str());
  } else {
    std::cerr << "Unknown mode: " << mode << "\n";
    return 1;
  }

  return 0;
}
