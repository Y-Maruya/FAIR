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
    HGLG_PEDESTAL_MASKED     = 1 << 6,  // Bit 6
    HGLG_MANUAL_BAD          = 1 << 7,  // Bit 7
    HGLG_MANUAL_GOOD         = 1 << 8   // Bit 8
};

struct QualityFlagEntry {
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;
  double x_mm = -999.0;
  double y_mm = -999.0;
  std::vector<uint32_t> flagsByRun;  // quality_flag for each run
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
                                       const std::string &branch, bool onlyFitOk,
                                       bool mask_L21C1 = true,
                                       bool mask_L3 = true,
                                       bool npoint_cut = true) {
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
    const double x_min = global_min - 0.05 * range_span;
    const double x_max = global_max + 0.05 * range_span;
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
  if (flag & HGLG_PEDESTAL_MASKED) result += "PED_MASKED|";
  if (flag & HGLG_MANUAL_BAD) result += "MANUAL_BAD|";
  if (flag & HGLG_MANUAL_GOOD) result += "MANUAL_GOOD|";
  if (!result.empty()) result.pop_back();  // Remove trailing '|'
  return result;
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
      // if (onlyFitOk && !entry.fit_ok) {
      //   continue;
      // }
      
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

  const Long64_t nEntries = tree->GetEntries();
  out.entries.reserve(nEntries);

  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    QualityFlagEntry e;
    e.cellid = p_cellid;
    e.layer = p_layer;
    e.chip = p_chip;
    e.channel = p_channel;
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
    out.entries.push_back(e);
    out.byCellId[e.cellid] = e;
  }

  file->Close();
  return true;
}

void printQualityFlagTable(const std::vector<QualityFlagRunData> &runs, const std::string &outDir) {
  if (runs.empty()) return;

  std::ofstream ofs((outDir + "/QualityFlags_table.txt").c_str());
  if (!ofs) return;

  ofs << "===============================================================================\n";
  ofs << "QUALITY FLAG COMPARISON TABLE\n";
  ofs << "===============================================================================\n\n";
  ofs << "Flags: GOOD=0, LOW_STAT(b0), FIT_FAILED(b1), BAD_RMS(b2), PED_MASKED(b6), ";
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
  for (const auto &dirName : dirNames) {
    if (gSystem->AccessPathName((std::string(baseDir) + "/" + dirName+"/hglg_calib_quality_muon.root").c_str()) != 0) {
      std::cerr << "  Warning: Directory not accessible: " << dirName << std::endl;
      continue;
    }
    RunData data;
    if (loadRun(baseDir, dirName, data)) {
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
    if (!shifted.empty()) {
      drawRepresentativeTrendSet(allRuns, shifted, branch, label, "shifted channels",
                                 outDir, "shifted", true);
      writeRepresentativeSummary(
          std::string(outDir) + "/shifted_" + branch + ".txt", shifted, shifted, branch);
    }
  }

  // Summary graphs
  drawSummaryGraphs(allRuns, outDir);

  // Parameter distributions
  drawParameterDistributions(allRuns, outDir);

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