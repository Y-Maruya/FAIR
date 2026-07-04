#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
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
#include "TSystem.h"
#include "TTree.h"

struct MipEntry {
  int cellid = -1;
  double MPV = 0.0;
  double width = 0.0;
  double TotalArea = 0.0;
  int entries = 0;
  double gaus_sigma = 0.0;
  double max_x = 0.0;
  double FWHM = 0.0;
  double chi2 = 0.0;
  int ndf = 0;
  double chi2perndf = 0.0;
  int fit_status = 0;
  int fit_ok = 1;
  double x_mm = 0.0;
  double y_mm = 0.0;
};

struct RunData {
  int run = 0;
  std::string label;
  std::string path;
  std::vector<MipEntry> entries;
  std::map<int, MipEntry> byCellId;
  bool hasFitFlag = false;
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

std::vector<int> parseRunList(const char *csv) {
  std::vector<int> runs;
  if (!csv) return runs;

  std::stringstream ss(csv);
  std::string token;
  while (std::getline(ss, token, ',')) {
    token.erase(std::remove_if(token.begin(), token.end(), ::isspace), token.end());
    if (token.empty()) continue;
    // support single run numbers and ranges like 22133-22136
    const auto dashPos = token.find('-');
    if (dashPos == std::string::npos) {
      char *endptr = nullptr;
      long val = std::strtol(token.c_str(), &endptr, 10);
      if (endptr != token.c_str() && *endptr == '\0') {
        runs.push_back(static_cast<int>(val));
      } else {
        std::cerr << "[compare] warning: invalid run token '" << token << "' skipped\n";
      }
    } else {
      std::string a = token.substr(0, dashPos);
      std::string b = token.substr(dashPos + 1);
      if (a.empty() || b.empty()) {
        std::cerr << "[compare] warning: invalid range '" << token << "' skipped\n";
        continue;
      }
      char *endptr = nullptr;
      long start = std::strtol(a.c_str(), &endptr, 10);
      if (endptr == a.c_str() || *endptr != '\0') {
        std::cerr << "[compare] warning: invalid range start '" << a << "' skipped\n";
        continue;
      }
      endptr = nullptr;
      long end = std::strtol(b.c_str(), &endptr, 10);
      if (endptr == b.c_str() || *endptr != '\0') {
        std::cerr << "[compare] warning: invalid range end '" << b << "' skipped\n";
        continue;
      }
      if (start > end) std::swap(start, end);
      for (long r = start; r <= end; ++r) runs.push_back(static_cast<int>(r));
    }
  }
  return runs;
}

std::vector<std::string> parseTokenList(const char *csv) {
  std::vector<std::string> tokens;
  if (!csv) return tokens;
  std::stringstream ss(csv);
  std::string token;
  while (std::getline(ss, token, ',')) {
    token.erase(std::remove_if(token.begin(), token.end(), ::isspace), token.end());
    if (token.empty()) continue;
    tokens.push_back(token);
  }
  return tokens;
}

std::string sanitizeLabel(const std::string &label) {
  std::string out = label;
  for (char &c : out) {
    if (!std::isalnum(static_cast<unsigned char>(c))) c = '_';
  }
  return out;
}

bool setBranchIfExists(TTree *tree, const char *name, void *addr) {
  if (!tree || !tree->GetBranch(name)) return false;
  tree->SetBranchAddress(name, addr);
  return true;
}

inline int HBUPositionOrder[40] = {
    39, 38, 37, 27, 14, 6,  7,  9,  12, 0, 2,  3,  5,  8,
    10, 11, 13, 15, 16, 1,  17, 18, 19, 20, 21, 22, 23, 24,
    25, 4,  26, 28, 29, 30, 31, 32, 33, 35, 34, 36};

inline int PosToLayerID(int position_order) {
  for (int i = 0; i < 40; ++i) {
    if (HBUPositionOrder[i] == position_order) {
      return i;
    }
  }
  return -1;
}

int mapCellIdForComparison(int rawCellId) {
  const int layer = rawCellId / 100000;
  const int new_layer = PosToLayerID(layer);
  if (new_layer < 0) return rawCellId;
  return new_layer * 100000 + (rawCellId % 100000);
}

int mapCellIdForComparisonInverse(int mappedCellId) {
  const int layer = mappedCellId / 100000;
  if (layer < 0 || layer >= 40) return mappedCellId;
  const int new_layer = HBUPositionOrder[layer];
  return new_layer * 100000 + (mappedCellId % 100000);
}

void setGraphStyle(TGraph *g, int color, int markerStyle) {
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(markerStyle);
  g->SetMarkerSize(1.1);
  g->SetLineWidth(2);
}

bool getBranchValue(const MipEntry &e, const std::string &branch, bool onlyFitOk,
                    double &value) {
  if (onlyFitOk && !e.fit_ok) return false;

  if (branch == "MPV") {
    value = e.MPV;
    return true;
  }
  if (branch == "width") {
    value = e.width;
    return true;
  }
  if (branch == "TotalArea") {
    value = e.TotalArea;
    return true;
  }
  if (branch == "gaus_sigma") {
    value = e.gaus_sigma;
    return true;
  }
  if (branch == "max_x") {
    value = e.max_x;
    return true;
  }
  if (branch == "FWHM") {
    value = e.FWHM;
    return true;
  }
  if (branch == "chi2") {
    value = e.chi2;
    return true;
  }
  if (branch == "chi2perndf") {
    value = e.chi2perndf;
    return true;
  }
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

SummaryStats extractStats(const RunData &runData, const std::string &what,
                          bool onlyFitOk) {
  std::vector<double> values;
  values.reserve(runData.entries.size());

  for (const auto &e : runData.entries) {
    double v = 0.0;
    if (getBranchValue(e, what, onlyFitOk, v)) values.push_back(v);
  }

  return computeStats(values);
}

std::vector<double> collectMetricValues(const RunData &run,
                                        const std::string &branch,
                                        bool onlyFitOk) {
  std::vector<double> values;
  values.reserve(run.entries.size());

  for (const auto &e : run.entries) {
    double v = 0.0;
    if (getBranchValue(e, branch, onlyFitOk, v)) {
      values.push_back(v);
    }
  }
  return values;
}

bool loadRunFromFile(const std::string &filePath,
                     int run,
                     const std::string &label,
                     bool applyCellIdMapping,
                     RunData &out) {
  out = RunData();
  out.run = run;
  out.label = label;
  out.path = filePath;

  TFile *file = TFile::Open(out.path.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "[compare] cannot open " << out.path << std::endl;
    if (file) file->Close();
    return false;
  }

  TTree *tree = dynamic_cast<TTree *>(file->Get("mip"));
  if (!tree) {
    std::cerr << "[compare] TTree 'mip' not found in " << out.path << std::endl;
    file->Close();
    return false;
  }

  MipEntry e;
  if (!setBranchIfExists(tree, "cellid", &e.cellid) ||
      !setBranchIfExists(tree, "MPV", &e.MPV) ||
      !setBranchIfExists(tree, "width", &e.width) ||
      !setBranchIfExists(tree, "TotalArea", &e.TotalArea) ||
      !setBranchIfExists(tree, "entries", &e.entries) ||
      !setBranchIfExists(tree, "gaus_sigma", &e.gaus_sigma) ||
      !setBranchIfExists(tree, "max_x", &e.max_x) ||
      !setBranchIfExists(tree, "FWHM", &e.FWHM) ||
      !setBranchIfExists(tree, "chi2", &e.chi2) ||
      !setBranchIfExists(tree, "ndf", &e.ndf) ||
      !setBranchIfExists(tree, "chi2perndf", &e.chi2perndf)) {
    std::cerr << "[compare] missing mandatory MIP branches in " << out.path
              << std::endl;
    file->Close();
    return false;
  }

  setBranchIfExists(tree, "fit_status", &e.fit_status);
  out.hasFitFlag = setBranchIfExists(tree, "fit_ok", &e.fit_ok);
  bool hasX = setBranchIfExists(tree, "x_mm", &e.x_mm);
  bool hasY = setBranchIfExists(tree, "y_mm", &e.y_mm);
  out.hasCoordinates = hasX && hasY;

  const Long64_t nEntries = tree->GetEntries();
  out.entries.reserve(nEntries);

  for (Long64_t i = 0; i < nEntries; ++i) {
    e.fit_ok = 1;
    e.x_mm = 0.0;
    e.y_mm = 0.0;
    tree->GetEntry(i);
    MipEntry mapped = e;
    mapped.cellid = applyCellIdMapping ? mapCellIdForComparison(e.cellid) : e.cellid;
    out.entries.push_back(mapped);
    out.byCellId[mapped.cellid] = mapped;
  }

  if (!out.hasFitFlag) {
    std::cout << "[compare] run " << run
              << " has no fit_ok branch, using all cells as valid" << std::endl;
  }
  if (!out.hasCoordinates) {
    std::cout << "[compare] run " << run
              << " has no x_mm/y_mm coordinates" << std::endl;
  }

  file->Close();
  return true;
}

bool loadRun(const std::string &baseDir, int run, RunData &out) {
  const std::string path = baseDir + "/" + std::to_string(run) + "/mip.root";
  return loadRunFromFile(path, run, std::to_string(run), false, out);
}

bool loadTestBeamRun(const std::string &baseDir,
                     const std::string &testBeamDir,
                     int pseudoRun,
                     RunData &out) {
  const std::string path = baseDir + "/" + testBeamDir + "/mip.root";
  return loadRunFromFile(path, pseudoRun, testBeamDir, true, out);
}

std::vector<CellShiftInfo> rankCellsBySpread(const std::vector<RunData> &runs,
                                             const std::string &branch,
                                             bool onlyFitOk) {
  std::vector<CellShiftInfo> ranked;
  if (runs.size() < 2) return ranked;

  const RunData &ref = runs.front();
  for (const auto &kv : ref.byCellId) {
    const int cellid = kv.first;
    std::vector<double> values;
    values.reserve(runs.size());

    double refValue = 0.0;
    bool haveRef = false;
    double x_mm = kv.second.x_mm;
    double y_mm = kv.second.y_mm;

    for (const auto &run : runs) {
      auto it = run.byCellId.find(cellid);
      if (it == run.byCellId.end()) {
        values.clear();
        break;
      }

      double value = 0.0;
      if (!getBranchValue(it->second, branch, onlyFitOk, value)) {
        values.clear();
        break;
      }

      if (!haveRef) {
        refValue = value;
        haveRef = true;
      }
      values.push_back(value);
    }

    if (values.size() != runs.size()) continue;

    auto mm = std::minmax_element(values.begin(), values.end());
    CellShiftInfo info;
    info.cellid = cellid;
    info.spread = *mm.second - *mm.first;
    info.refValue = refValue;
    info.meanValue = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    info.x_mm = x_mm;
    info.y_mm = y_mm;
    ranked.push_back(info);
  }

  std::sort(ranked.begin(), ranked.end(),
            [](const CellShiftInfo &a, const CellShiftInfo &b) {
              if (a.spread != b.spread) return a.spread > b.spread;
              return a.cellid < b.cellid;
            });
  return ranked;
}

void writeRepresentativeSummary(const std::string &fileName,
                                const std::vector<CellShiftInfo> &shifted,
                                const std::vector<CellShiftInfo> &stable,
                                const std::string &branch) {
  std::ofstream ofs(fileName.c_str());
  if (!ofs) return;

  ofs << "# branch = " << branch << "\n";
  ofs << "# shifted channels\n";
  ofs << "# cellid spread refValue meanValue x_mm y_mm\n";
  for (const auto &info : shifted) {
    ofs << info.cellid << " " << info.spread << " " << info.refValue << " "
        << info.meanValue << " " << info.x_mm << " " << info.y_mm << "\n";
  }

  ofs << "\n# stable channels\n";
  ofs << "# cellid spread refValue meanValue x_mm y_mm\n";
  for (const auto &info : stable) {
    ofs << info.cellid << " " << info.spread << " " << info.refValue << " "
        << info.meanValue << " " << info.x_mm << " " << info.y_mm << "\n";
  }
}

bool isTestBeamLikePath(const std::string &path) {
  std::string lower = path;
  std::transform(lower.begin(), lower.end(), lower.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return lower.find("ehn1") != std::string::npos ||
         lower.find("testbeam") != std::string::npos;
}

TH1 *loadMipHistogram(const std::string &filePath,
                      int mappedCellId,
                      bool inverseMapForLookup) {
  int lookupCellId = mappedCellId;
  if (inverseMapForLookup) {
    lookupCellId = mapCellIdForComparisonInverse(mappedCellId);
  }

  const int layer = lookupCellId / 100000;
  const int chip = (lookupCellId / 10000) % 10;

  TFile *file = TFile::Open(filePath.c_str(), "READ");
  if (!file || file->IsZombie()) {
    if (file) file->Close();
    return nullptr;
  }

  std::vector<std::string> candidates;
  candidates.push_back(Form("MIP/Layer%02d/Chip%d/hMIP_%d", layer, chip, lookupCellId));
  candidates.push_back(Form("MIP/Layer%d/Chip%d/hMIP_%d", layer, chip, lookupCellId));
  candidates.push_back(Form("MIP/Layer%02d/Chip%02d/hMIP_%d", layer, chip, lookupCellId));
  candidates.push_back(Form("hMIP_%d", lookupCellId));

  TH1 *hist = nullptr;
  for (const auto &name : candidates) {
    hist = dynamic_cast<TH1 *>(file->Get(name.c_str()));
    if (hist) break;
  }

  TH1 *cloned = nullptr;
  if (hist) {
    static long long cloneCounter = 0;
    ++cloneCounter;
    cloned = dynamic_cast<TH1 *>(
        hist->Clone(Form("hMIP_%d_clone_%lld", mappedCellId, cloneCounter)));
    if (cloned) cloned->SetDirectory(nullptr);
  }

  file->Close();
  return cloned;
}

void drawRepresentativeHistogramSet(const std::vector<RunData> &runs,
                                    const std::vector<CellShiftInfo> &cells,
                                    const std::string &label,
                                    const std::string &outDir,
                                    const std::string &canvasTag) {
  if (cells.empty()) return;

  const int n = static_cast<int>(cells.size());
  const int nCols = 2;
  const int nRows = (n + nCols - 1) / nCols;
  TCanvas *c = new TCanvas(Form("c_mip_%s", canvasTag.c_str()),
                           Form("MIP %s", label.c_str()), 1400, 450 * nRows);
  c->Divide(nCols, nRows);

  int colors[] = {kRed + 1,     kBlue + 1,    kGreen + 2,  kMagenta + 1,
                  kOrange + 7,  kCyan + 2,    kRed + 3,    kBlue + 3,
                  kGreen + 4,   kMagenta + 3, kOrange + 9, kCyan + 4};

  for (int i = 0; i < n; ++i) {
    c->cd(i + 1);
    gPad->SetGrid();

    TLegend *leg = new TLegend(0.54, 0.60, 0.88, 0.88);
    bool first = true;
    double maxY = 0.0;
    std::vector<TH1 *> drawn;

    for (size_t irun = 0; irun < runs.size(); ++irun) {
      const bool isTb = isTestBeamLikePath(runs[irun].path);
      TH1 *h = loadMipHistogram(runs[irun].path, cells[i].cellid, isTb);
      if (!h) continue;

      if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
    //   h->GetXaxis()->SetRangeUser(0, 3000);
      h->Rebin(16);
      h->SetLineColor(colors[irun % 12]);
      h->SetLineWidth(2);
      h->SetStats(0);
      h->SetTitle(Form("MPV cell %d (spread = %.3f);ADC counts;Normalized entries",
                       cells[i].cellid, cells[i].spread));

      maxY = std::max(maxY, h->GetMaximum());
      if (first) {
        h->Draw("HIST");
        first = false;
      } else {
        h->Draw("HIST SAME");
      }
      leg->AddEntry(h, runs[irun].label.c_str(), "l");
      drawn.push_back(h);
    }

    if (!drawn.empty()) {
      drawn.front()->SetMaximum(maxY * 1.25);
      bool haveCoords = std::abs(cells[i].x_mm) > 1e-12 || std::abs(cells[i].y_mm) > 1e-12;
      if (haveCoords) {
        leg->SetHeader(Form("cellid=%d, (x,y)=(%.1f, %.1f)", cells[i].cellid,
                            cells[i].x_mm, cells[i].y_mm),
                       "C");
      } else {
        leg->SetHeader(Form("cellid=%d, coords=N/A", cells[i].cellid), "C");
      }
      leg->Draw();
    }
  }

  c->SaveAs((outDir + "/MPV_" + canvasTag + ".pdf").c_str());
  c->SaveAs((outDir + "/MPV_" + canvasTag + ".png").c_str());
}

void drawRepresentativeHistogramZoomSet(const std::vector<RunData> &runs,
                                    const std::vector<CellShiftInfo> &cells,
                                    const std::string &label,
                                    const std::string &outDir,
                                    const std::string &canvasTag) {
  if (cells.empty()) return;

  const int n = static_cast<int>(cells.size());
  const int nCols = 2;
  const int nRows = (n + nCols - 1) / nCols;
  TCanvas *c = new TCanvas(Form("c_mip_%s", canvasTag.c_str()),
                           Form("MIP %s", label.c_str()), 1400, 450 * nRows);
  c->Divide(nCols, nRows);

  int colors[] = {kRed + 1,     kBlue + 1,    kGreen + 2,  kMagenta + 1,
                  kOrange + 7,  kCyan + 2,    kRed + 3,    kBlue + 3,
                  kGreen + 4,   kMagenta + 3, kOrange + 9, kCyan + 4};

  for (int i = 0; i < n; ++i) {
    c->cd(i + 1);
    gPad->SetGrid();

    TLegend *leg = new TLegend(0.54, 0.60, 0.88, 0.88);
    bool first = true;
    double maxY = 0.0;
    std::vector<TH1 *> drawn;

    for (size_t irun = 0; irun < runs.size(); ++irun) {
      const bool isTb = isTestBeamLikePath(runs[irun].path);
      TH1 *h = loadMipHistogram(runs[irun].path, cells[i].cellid, isTb);
      if (!h) continue;

      if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
      h->Rebin(16);
      h->GetXaxis()->SetRangeUser(240, 480+240);
      h->SetLineColor(colors[irun % 12]);
      h->SetLineWidth(2);
      h->SetStats(0);
      h->SetTitle(Form("MPV cell %d (spread = %.3f);ADC counts;Normalized entries",
                       cells[i].cellid, cells[i].spread));

      maxY = std::max(maxY, h->GetMaximum());
      if (first) {
        h->Draw("HIST");
        first = false;
      } else {
        h->Draw("HIST SAME");
      }
      leg->AddEntry(h, runs[irun].label.c_str(), "l");
      drawn.push_back(h);
    }

    if (!drawn.empty()) {
      drawn.front()->SetMaximum(maxY * 1.25);
      bool haveCoords = std::abs(cells[i].x_mm) > 1e-12 || std::abs(cells[i].y_mm) > 1e-12;
      if (haveCoords) {
        leg->SetHeader(Form("cellid=%d, (x,y)=(%.1f, %.1f)", cells[i].cellid,
                            cells[i].x_mm, cells[i].y_mm),
                       "C");
      } else {
        leg->SetHeader(Form("cellid=%d, coords=N/A", cells[i].cellid), "C");
      }
      leg->Draw();
    }
  }

  c->SaveAs((outDir + "/MPV_zoom_" + canvasTag + ".pdf").c_str());
  c->SaveAs((outDir + "/MPV_zoom_" + canvasTag + ".png").c_str());
}

void drawRepresentativeTrendSet(const std::vector<RunData> &runs,
                                const std::vector<CellShiftInfo> &cells,
                                const std::string &branch,
                                const std::string &label,
                                const std::string &outDir,
                                const std::string &canvasTag,
                                bool onlyFitOk) {
  if (cells.empty()) return;

  const int n = static_cast<int>(cells.size());
  const int nCols = 2;
  const int nRows = (n + nCols - 1) / nCols;
  TCanvas *c = new TCanvas(Form("c_%s_%s", branch.c_str(), canvasTag.c_str()),
                           Form("%s %s", branch.c_str(), label.c_str()), 1400,
                           420 * nRows);
  c->Divide(nCols, nRows);

  for (int i = 0; i < n; ++i) {
    c->cd(i + 1);
    gPad->SetGrid();

    std::vector<double> x;
    std::vector<double> y;
    x.reserve(runs.size());
    y.reserve(runs.size());

    for (const auto &run : runs) {
      auto it = run.byCellId.find(cells[i].cellid);
      if (it == run.byCellId.end()) continue;
      double v = 0.0;
      if (!getBranchValue(it->second, branch, onlyFitOk, v)) continue;
      x.push_back(run.run);
      y.push_back(v);
    }

    if (x.empty()) continue;

    TGraph *g = new TGraph(static_cast<int>(x.size()), &x[0], &y[0]);
    setGraphStyle(g, kBlue + 1, 20);
    g->SetTitle(Form("%s trend, cell %d (spread=%.3f);Run;%s",
                     branch.c_str(), cells[i].cellid, cells[i].spread,
                     branch.c_str()));
    g->Draw("ALP");
  }

  c->SaveAs((outDir + "/" + branch + "_" + canvasTag + ".pdf").c_str());
  c->SaveAs((outDir + "/" + branch + "_" + canvasTag + ".png").c_str());
}

void drawRepresentativeHistograms(const std::vector<RunData> &runs,
                                  const std::string &outDir,
                                  int nRepresentative = 8) {
  const std::string branch = "MPV";
  std::vector<CellShiftInfo> ranked = rankCellsBySpread(runs, branch, true);
  if (ranked.empty()) return;

  const int nTake = std::min(nRepresentative, static_cast<int>(ranked.size()));
  std::vector<CellShiftInfo> shifted(ranked.begin(), ranked.begin() + nTake);

  std::vector<CellShiftInfo> stable;
  for (auto it = ranked.rbegin(); it != ranked.rend() &&
                               static_cast<int>(stable.size()) < nTake;
       ++it) {
    stable.push_back(*it);
  }
  std::reverse(stable.begin(), stable.end());

  writeRepresentativeSummary(outDir + "/" + branch +
                                 "_representative_channels.txt",
                             shifted, stable, branch);

  drawRepresentativeHistogramSet(runs, shifted, "shifted channels", outDir,
                                 "shifted_channels");
  drawRepresentativeHistogramSet(runs, stable, "stable channels", outDir,
                                 "stable_channels");

  drawRepresentativeHistogramZoomSet(runs, shifted, "shifted channels (zoom)", outDir,
                                 "shifted_channels_zoom");
  drawRepresentativeHistogramZoomSet(runs, stable, "stable channels (zoom)", outDir,
                                 "stable_channels_zoom");
//   drawRepresentativeTrendSet(runs, shifted, branch, "shifted channels", outDir,
//                              "shifted_channels", true);
//   drawRepresentativeTrendSet(runs, stable, branch, "stable channels", outDir,
//                              "stable_channels", true);
}

TH1D *makeDeltaHist(const char *name,
                    const char *title,
                    const std::vector<RunData> &runs,
                    const std::string &branch,
                    bool onlyFitOk) {
  if (runs.size() < 2) return nullptr;

  std::vector<double> deltas;
  const RunData &ref = runs.front();

  for (size_t i = 1; i < runs.size(); ++i) {
    for (const auto &kv : ref.byCellId) {
      const int cellid = kv.first;
      const auto it = runs[i].byCellId.find(cellid);
      if (it == runs[i].byCellId.end()) continue;

      double va = 0.0;
      double vb = 0.0;
      if (!getBranchValue(kv.second, branch, onlyFitOk, va)) continue;
      if (!getBranchValue(it->second, branch, onlyFitOk, vb)) continue;
      deltas.push_back(vb - va);
    }
  }

  if (deltas.empty()) return nullptr;

  auto mm = std::minmax_element(deltas.begin(), deltas.end());
  double xmin = *mm.first;
  double xmax = *mm.second;
  if (xmin == xmax) {
    xmin -= 1.0;
    xmax += 1.0;
  }

  const double span = std::max(1e-6, xmax - xmin);
  TH1D *h = new TH1D(name, title, 80, xmin - 0.1 * span, xmax + 0.1 * span);
  for (double d : deltas) h->Fill(d);
  h->SetLineWidth(2);
  return h;
}

std::vector<double> collectDeltaValues(const RunData &ref,
                                       const RunData &target,
                                       const std::string &branch,
                                       bool onlyFitOk) {
  std::vector<double> deltas;
  deltas.reserve(ref.byCellId.size());

  for (const auto &kv : ref.byCellId) {
    const int cellid = kv.first;
    const auto it = target.byCellId.find(cellid);
    if (it == target.byCellId.end()) continue;

    double va = 0.0;
    double vb = 0.0;
    if (!getBranchValue(kv.second, branch, onlyFitOk, va)) continue;
    if (!getBranchValue(it->second, branch, onlyFitOk, vb)) continue;
    deltas.push_back(vb - va);
  }

  return deltas;
}

void drawDeltaFromReferenceByRun(const std::vector<RunData> &runs,
                                 const std::string &branch,
                                 const std::string &title,
                                 const std::string &outTag,
                                 const std::string &outDir,
                                 bool onlyFitOk = true) {
  if (runs.size() < 2) return;

  const RunData &ref = runs.front();
  std::vector<std::vector<double> > allDeltas(runs.size());
  double globalMin = 1e30;
  double globalMax = -1e30;

  for (size_t i = 1; i < runs.size(); ++i) {
    allDeltas[i] = collectDeltaValues(ref, runs[i], branch, onlyFitOk);
    if (allDeltas[i].empty()) continue;

    auto mm = std::minmax_element(allDeltas[i].begin(), allDeltas[i].end());
    globalMin = std::min(globalMin, *mm.first);
    globalMax = std::max(globalMax, *mm.second);
  }

  if (globalMin > globalMax) return;
  if (globalMin == globalMax) {
    globalMin -= 1.0;
    globalMax += 1.0;
  }

  const double span = std::max(1e-6, globalMax - globalMin);
  const double xmin = globalMin - 0.1 * span;
  const double xmax = globalMax + 0.1 * span;
  double ymax = 0.0;

  TCanvas *c =
      new TCanvas(Form("c_%s_delta_ref", outTag.c_str()), title.c_str(), 1000, 700);
  c->SetGrid();

  TLegend *leg = new TLegend(0.62, 0.64, 0.88, 0.88);
  int colors[] = {kRed + 1,     kBlue + 1,    kGreen + 2,  kMagenta + 1,
                  kOrange + 7,  kCyan + 2,    kRed + 3,    kBlue + 3,
                  kGreen + 4,   kMagenta + 3, kOrange + 9, kCyan + 4};
  bool first = true;
  std::vector<TH1D *> drawn;

  for (size_t i = 1; i < runs.size(); ++i) {
    if (allDeltas[i].empty()) continue;

    const std::string san = sanitizeLabel(runs[i].label);
    const std::string sanRef = sanitizeLabel(ref.label);
    TH1D *h = new TH1D(
      Form("h_%s_delta_ref_%s", outTag.c_str(), san.c_str()),
      Form("%s (ref %s);#Delta %s;Normalized entries", title.c_str(),
         ref.label.c_str(), branch.c_str()),
      80, xmin, xmax);
    for (double d : allDeltas[i]) h->Fill(d);
    if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
    if (ymax < h->GetMaximum()) ymax = h->GetMaximum() * 1.25;
    h->SetLineColor(colors[(i - 1) % 12]);
    h->SetLineWidth(2);
    h->SetStats(0);

    if (first) {
      h->Draw("HIST");
      first = false;
    } else {
      h->Draw("HIST SAME");
    }
    drawn.push_back(h);
    leg->AddEntry(h, Form("%s - %s", runs[i].label.c_str(), ref.label.c_str()), "l");
  }

  for (TH1D *hist : drawn) hist->SetMaximum(ymax);

  if (!first) {
    leg->Draw();
    c->SaveAs((outDir + "/" + outTag + "_delta_from_first_run_overlay.pdf").c_str());
    c->SaveAs((outDir + "/" + outTag + "_delta_from_first_run_overlay.png").c_str());
  }
}

void drawDeltaByLayer(const std::vector<RunData> &runs,
                      const std::string &branch,
                      const std::string &outDir,
                      bool onlyFitOk = true) {
  if (runs.size() < 2) return;

  const RunData &ref = runs.front();
  int nLayers = 40;
  int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2};

  for (size_t irun = 1; irun < runs.size(); ++irun) {
    std::vector<double> sum(nLayers, 0.0);
    std::vector<int> count(nLayers, 0);

    for (const auto &kv : ref.byCellId) {
      const int cellid = kv.first;
      const int layer = cellid / 100000;
      if (layer < 0 || layer >= nLayers) continue;

      auto it = runs[irun].byCellId.find(cellid);
      if (it == runs[irun].byCellId.end()) continue;

      double va = 0.0, vb = 0.0;
      if (!getBranchValue(kv.second, branch, onlyFitOk, va)) continue;
      if (!getBranchValue(it->second, branch, onlyFitOk, vb)) continue;

      sum[layer] += (vb - va);
      ++count[layer];
    }

    std::vector<int> layers;
    std::vector<double> meanDelta;
    for (int L = 0; L < nLayers; ++L) {
      if (count[L] > 0) {
        layers.push_back(L);
        meanDelta.push_back(sum[L] / count[L]);
      }
    }

    if (layers.empty()) continue;

    std::vector<double> xd(layers.size()), yd(meanDelta.size());
    for (size_t i = 0; i < layers.size(); ++i) {
      xd[i] = static_cast<double>(layers[i]);
      yd[i] = meanDelta[i];
    }

    TCanvas *c = new TCanvas(Form("c_delta_layer_%s_%s", branch.c_str(), runs[irun].label.c_str()),
                             Form("%s mean #Delta per layer (%s - %s)", branch.c_str(), runs[irun].label.c_str(), ref.label.c_str()),
                             1200, 600);
    c->SetGrid();

    TGraph *g = new TGraph(static_cast<int>(xd.size()), &xd[0], &yd[0]);
    setGraphStyle(g, colors[(irun - 1) % 6], 20);
    g->SetTitle(Form("%s mean #Delta per layer;Layer;#Delta %s", branch.c_str(), branch.c_str()));
    g->Draw("ALP");

    // annotate counts per layer in the legend
    TLegend *leg = new TLegend(0.65, 0.70, 0.88, 0.88);
    leg->AddEntry(g, Form("%s - %s (mean per layer)", runs[irun].label.c_str(), ref.label.c_str()), "lp");
    leg->Draw();

    const std::string san = sanitizeLabel(runs[irun].label);
    c->SaveAs((outDir + Form("/delta_%s_by_layer_%s.pdf", branch.c_str(), san.c_str())).c_str());
    c->SaveAs((outDir + Form("/delta_%s_by_layer_%s.png", branch.c_str(), san.c_str())).c_str());
  }
}

void drawLayerOverlays(const std::vector<RunData> &runs,
                       const std::string &outDir,
                       bool onlyFitOk = true) {
  if (runs.size() < 2) return;

  const RunData &ref = runs.front();
  const int nLayers = 40;

  for (int layer = 0; layer < nLayers; ++layer) {
    std::vector<double> x_runs;
    std::vector<double> delta_mpv;
    std::vector<double> delta_fwhm;

    // compute reference means for this layer
    auto computeMeanForLayer = [&](const RunData &r, const std::string &branch) -> std::pair<bool, double> {
      double sum = 0.0;
      int cnt = 0;
      for (const auto &kv : r.byCellId) {
        const int cid = kv.first;
        const int lay = cid / 100000;
        if (lay != layer) continue;
        double v = 0.0;
        if (!getBranchValue(kv.second, branch, onlyFitOk, v)) continue;
        sum += v;
        ++cnt;
      }
      if (cnt == 0) return {false, 0.0};
      return {true, sum / cnt};
    };

    auto ref_mpv_pair = computeMeanForLayer(ref, "MPV");
    auto ref_fwhm_pair = computeMeanForLayer(ref, "FWHM");
    if (!ref_mpv_pair.first && !ref_fwhm_pair.first) continue;

    for (size_t ir = 0; ir < runs.size(); ++ir) {
      const RunData &r = runs[ir];
      auto mv_mpv = computeMeanForLayer(r, "MPV");
      auto mv_fwhm = computeMeanForLayer(r, "FWHM");

      // require that at least one of MPV/FWHM is present to include this run
      if (!mv_mpv.first && !mv_fwhm.first) continue;

      x_runs.push_back(static_cast<double>(r.run));
      if (mv_mpv.first && ref_mpv_pair.first)
        delta_mpv.push_back(mv_mpv.second - ref_mpv_pair.second);
      else
        delta_mpv.push_back(0.0/0.0); // NaN placeholder for missing

      if (mv_fwhm.first && ref_fwhm_pair.first)
        delta_fwhm.push_back(mv_fwhm.second - ref_fwhm_pair.second);
      else
        delta_fwhm.push_back(0.0/0.0);
    }

    if (x_runs.empty()) continue;

    // Filter out NaN entries so graphs have matching points
    std::vector<double> x_mpv, y_mpv, x_fwhm, y_fwhm;
    for (size_t i = 0; i < x_runs.size(); ++i) {
      if (std::isfinite(delta_mpv[i])) {
        x_mpv.push_back(x_runs[i]);
        y_mpv.push_back(delta_mpv[i]);
      }
      if (std::isfinite(delta_fwhm[i])) {
        x_fwhm.push_back(x_runs[i]);
        y_fwhm.push_back(delta_fwhm[i]);
      }
    }

    if (x_mpv.empty() && x_fwhm.empty()) continue;

    TCanvas *c = new TCanvas(Form("c_layer_overlay_%02d", layer),
                             Form("Layer %02d: delta MPV/FWHM overlay", layer), 1200, 600);
    c->SetGrid();

    TLegend *leg = new TLegend(0.62, 0.70, 0.88, 0.88);

    if (!x_mpv.empty()) {
      TGraph *gmpv = new TGraph(static_cast<int>(x_mpv.size()), &x_mpv[0], &y_mpv[0]);
      setGraphStyle(gmpv, kRed + 1, 20);
      gmpv->SetTitle(Form("Layer %02d;Run;#Delta MPV/FWHM", layer));
      gmpv->Draw("ALP");
      leg->AddEntry(gmpv, "#Delta MPV", "lp");
    }

    if (!x_fwhm.empty()) {
      TGraph *gfwhm = new TGraph(static_cast<int>(x_fwhm.size()), &x_fwhm[0], &y_fwhm[0]);
      setGraphStyle(gfwhm, kBlue + 1, 21);
      if (!x_mpv.empty()) gfwhm->Draw("LP SAME");
      else gfwhm->Draw("ALP");
      leg->AddEntry(gfwhm, "#Delta FWHM", "lp");
    }

    leg->Draw();

    c->SaveAs((outDir + Form("/layer_%02d_delta_overlay.pdf", layer)).c_str());
    c->SaveAs((outDir + Form("/layer_%02d_delta_overlay.png", layer)).c_str());
  }
}

void drawDeltaByLayerOverlay(const std::vector<RunData> &runs,
                             const std::string &branch,
                             const std::string &outDir,
                             bool onlyFitOk = true) {
  if (runs.size() < 2) return;

  const RunData &ref = runs.front();
  const int nLayers = 40;
  int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2,
                  kRed + 3, kBlue + 3, kGreen + 4, kMagenta + 3, kOrange + 9, kCyan + 4};

  for (int layer = 0; layer < nLayers; ++layer) {
    std::vector<std::vector<double>> layerDeltas(runs.size());
    double globalMin = 1e30;
    double globalMax = -1e30;

    // compute deltas for each run (relative to ref) in this layer only
    for (size_t ir = 1; ir < runs.size(); ++ir) {
      for (const auto &kv : ref.byCellId) {
        const int cellid = kv.first;
        const int cLayer = cellid / 100000;
        if (cLayer != layer) continue;

        auto it = runs[ir].byCellId.find(cellid);
        if (it == runs[ir].byCellId.end()) continue;

        double va = 0.0, vb = 0.0;
        if (!getBranchValue(kv.second, branch, onlyFitOk, va)) continue;
        if (!getBranchValue(it->second, branch, onlyFitOk, vb)) continue;

        double delta = vb - va;
        layerDeltas[ir].push_back(delta);
        globalMin = std::min(globalMin, delta);
        globalMax = std::max(globalMax, delta);
      }
    }

    if (globalMin > globalMax) continue;
    if (globalMin == globalMax) {
      globalMin -= 1.0;
      globalMax += 1.0;
    }

    const double span = std::max(1e-6, globalMax - globalMin);
    const double xmin = globalMin - 0.1 * span;
    const double xmax = globalMax + 0.1 * span;
    double ymax = 0.0;

    TCanvas *c = new TCanvas(Form("c_layer_%02d_delta_%s_overlay", layer, branch.c_str()),
                             Form("Layer %02d #Delta %s overlay", layer, branch.c_str()),
                             1000, 700);
    c->SetGrid();

    TLegend *leg = new TLegend(0.62, 0.64, 0.88, 0.88);
    bool first = true;
    std::vector<TH1D *> drawn;

    for (size_t ir = 1; ir < runs.size(); ++ir) {
      if (layerDeltas[ir].empty()) continue;

      const std::string san = sanitizeLabel(runs[ir].label);
      TH1D *h = new TH1D(
        Form("h_layer_%02d_delta_%s_%s", layer, branch.c_str(), san.c_str()),
        Form("Layer %02d #Delta %s;#Delta %s;Normalized entries", layer, branch.c_str(), branch.c_str()),
        80, xmin, xmax);

      for (double d : layerDeltas[ir]) h->Fill(d);
      if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
      if (ymax < h->GetMaximum()) ymax = h->GetMaximum() * 1.25;
      h->SetLineColor(colors[(ir - 1) % 12]);
      h->SetLineWidth(2);
      h->SetStats(0);

      if (first) {
        h->Draw("HIST");
        first = false;
      } else {
        h->Draw("HIST SAME");
      }
      drawn.push_back(h);
      leg->AddEntry(h, Form("%s - %s", runs[ir].label.c_str(), ref.label.c_str()), "l");
    }

    for (TH1D *hist : drawn) hist->SetMaximum(ymax);

    if (!first) {
      leg->Draw();
      c->SaveAs((outDir + Form("/layer_%02d_delta_%s_overlay.pdf", layer, branch.c_str())).c_str());
      c->SaveAs((outDir + Form("/layer_%02d_delta_%s_overlay.png", layer, branch.c_str())).c_str());
    }
  }
}

void drawSummaryGraphs(const std::vector<RunData> &runs, const std::string &outDir) {
  const int n = static_cast<int>(runs.size());
  if (n == 0) return;
  std::vector<double> x(n), mpv(n), width(n), fwhm(n), chi2ndf(n), fitok(n), gausSigma(n);
  for (int i = 0; i < n; ++i) {
    x[i] = runs[i].run;

    mpv[i] = extractStats(runs[i], "MPV", true).mean;
    width[i] = extractStats(runs[i], "width", true).mean;
    fwhm[i] = extractStats(runs[i], "FWHM", true).mean;
    chi2ndf[i] = extractStats(runs[i], "chi2perndf", true).mean;
    gausSigma[i] = extractStats(runs[i], "gaus_sigma", true).mean;

    int ok = 0;
    for (const auto &e : runs[i].entries) {
      if (e.fit_ok) ++ok;
    }
    fitok[i] = 100.0 * ok / std::max(1, static_cast<int>(runs[i].entries.size()));
  }

  // Produce separate canvases (one metric per canvas) to avoid overlaying
  // and to produce 8 output files as requested.

  // 1) MPV vs run
  {
    TCanvas *c = new TCanvas("c_mpv", "Mean MPV vs run", 1000, 700);
    TGraph *g = new TGraph(n, &x[0], &mpv[0]);
    setGraphStyle(g, kRed + 1, 20);
    g->SetTitle("Mean MPV vs run;Run;ADC counts");
    g->Draw("ALP");
    c->SetGrid();
    c->SaveAs((outDir + "/summary_MPV.pdf").c_str());
    c->SaveAs((outDir + "/summary_MPV.png").c_str());
  }

  // 2) width vs run
  {
    TCanvas *c = new TCanvas("c_width", "Mean width vs run", 1000, 700);
    TGraph *g = new TGraph(n, &x[0], &width[0]);
    setGraphStyle(g, kBlue + 1, 21);
    g->SetTitle("Mean width vs run;Run;ADC counts");
    g->Draw("ALP");
    c->SetGrid();
    c->SaveAs((outDir + "/summary_width.pdf").c_str());
    c->SaveAs((outDir + "/summary_width.png").c_str());
  }

  // 3) FWHM vs run
  {
    TCanvas *c = new TCanvas("c_fwhm", "Mean FWHM vs run", 1000, 700);
    TGraph *g = new TGraph(n, &x[0], &fwhm[0]);
    setGraphStyle(g, kMagenta + 1, 20);
    g->SetTitle("Mean FWHM vs run;Run;ADC counts");
    g->Draw("ALP");
    c->SetGrid();
    c->SaveAs((outDir + "/summary_FWHM.pdf").c_str());
    c->SaveAs((outDir + "/summary_FWHM.png").c_str());
  }

  // 4) gaus_sigma vs run
  {
    TCanvas *c = new TCanvas("c_gaus", "Mean gaus_sigma vs run", 1000, 700);
    TGraph *g = new TGraph(n, &x[0], &gausSigma[0]);
    setGraphStyle(g, kCyan + 2, 21);
    g->SetTitle("Mean gaus_sigma vs run;Run;ADC counts");
    g->Draw("ALP");
    c->SetGrid();
    c->SaveAs((outDir + "/summary_gaus_sigma.pdf").c_str());
    c->SaveAs((outDir + "/summary_gaus_sigma.png").c_str());
  }

  // 5) chi2perndf vs run
  {
    TCanvas *c = new TCanvas("c_chi2ndf", "Mean chi2/ndf vs run", 1000, 700);
    TGraph *g = new TGraph(n, &x[0], &chi2ndf[0]);
    setGraphStyle(g, kGreen + 2, 20);
    g->SetTitle("Mean chi2/ndf vs run;Run;chi2/ndf");
    g->Draw("ALP");
    c->SetGrid();
    c->SaveAs((outDir + "/summary_chi2perndf.pdf").c_str());
    c->SaveAs((outDir + "/summary_chi2perndf.png").c_str());
  }

  // 6) fit_ok percentage vs run
  {
    TCanvas *c = new TCanvas("c_fitok", "fit_ok [%] vs run", 1000, 700);
    TGraph *g = new TGraph(n, &x[0], &fitok[0]);
    setGraphStyle(g, kOrange + 7, 21);
    g->SetTitle("fit_ok [%] vs run;Run;fit_ok [%]");
    g->Draw("ALP");
    c->SetGrid();
    c->SaveAs((outDir + "/summary_fit_ok.pdf").c_str());
    c->SaveAs((outDir + "/summary_fit_ok.png").c_str());
  }

  // 7) Delta MPV histogram (relative to first run)
  {
    TH1D *hDeltaMpv = makeDeltaHist("hDeltaMpv",
                                    Form("#Delta MPV relative to run %d;#Delta MPV;Cells",
                                         runs.front().run),
                                    runs, "MPV", true);
    if (hDeltaMpv) {
      TCanvas *c = new TCanvas("c_deltampv", "#Delta MPV", 1000, 700);
      hDeltaMpv->SetLineColor(kRed + 1);
      hDeltaMpv->SetLineWidth(2);
      hDeltaMpv->Draw("HIST");
      c->SetGrid();
      c->SaveAs((outDir + "/delta_MPV_from_first_run.pdf").c_str());
      c->SaveAs((outDir + "/delta_MPV_from_first_run.png").c_str());
    }
  }

  // 8) Delta FWHM histogram (relative to first run)
  {
    TH1D *hDeltaFwhm = makeDeltaHist("hDeltaFwhm",
                                     Form("#Delta FWHM relative to run %d;#Delta FWHM;Cells",
                                          runs.front().run),
                                     runs, "FWHM", true);
    if (hDeltaFwhm) {
      TCanvas *c = new TCanvas("c_deltafwhm", "#Delta FWHM", 1000, 700);
      hDeltaFwhm->SetLineColor(kBlue + 1);
      hDeltaFwhm->SetLineWidth(2);
      hDeltaFwhm->Draw("HIST");
      c->SetGrid();
      c->SaveAs((outDir + "/delta_FWHM_from_first_run.pdf").c_str());
      c->SaveAs((outDir + "/delta_FWHM_from_first_run.png").c_str());
    }
  }
}

void drawOverlayHistograms(const std::vector<RunData> &runs, const std::string &outDir) {
  struct PlotDef {
    std::string key;
    std::string title;
    std::string xTitle;
    bool onlyFitOk;
  };

  std::vector<PlotDef> defs = {
      {"MPV", "MIP MPV", "ADC counts", true},
      {"width", "MIP width", "ADC counts", true},
      {"FWHM", "MIP FWHM", "ADC counts", true},
      {"chi2perndf", "MIP chi2/ndf", "chi2/ndf", true},
      {"gaus_sigma", "MIP gaus sigma", "ADC counts", true},
      {"TotalArea", "MIP TotalArea", "a.u.", true},
  };

  int colors[] = {kRed + 1,     kBlue + 1,    kGreen + 2,  kMagenta + 1,
                  kOrange + 7,  kCyan + 2,    kRed + 3,    kBlue + 3,
                  kGreen + 4,   kMagenta + 3, kOrange + 9, kCyan + 4};

  for (const auto &def : defs) {
    double globalMin = 1e30;
    double globalMax = -1e30;

    for (const auto &run : runs) {
      std::vector<double> values = collectMetricValues(run, def.key, def.onlyFitOk);
      if (values.empty()) continue;
      SummaryStats s = computeStats(values);
      globalMin = std::min(globalMin, s.min);
      globalMax = std::max(globalMax, s.max);
    }

    if (globalMin > globalMax) continue;
    if (globalMin == globalMax) {
      globalMin -= 1.0;
      globalMax += 1.0;
    }

    const double span = std::max(1e-6, globalMax - globalMin);
    const double xmin = globalMin - 0.05 * span;
    const double xmax = globalMax + 0.05 * span;
    const int bins = 120;

    TCanvas *c =
        new TCanvas(Form("c_%s", def.key.c_str()), def.title.c_str(), 1000, 700);
    TLegend *leg = new TLegend(0.65, 0.70, 0.88, 0.88);
    bool first = true;
    double ymax = 0.0;
    std::vector<TH1D *> drawn;

    for (size_t i = 0; i < runs.size(); ++i) {
      const std::string san = sanitizeLabel(runs[i].label);
      TH1D *h = new TH1D(
          Form("h_%s_%s", def.key.c_str(), san.c_str()),
          Form("%s;%s;Normalized entries", def.title.c_str(), def.xTitle.c_str()), bins,
          xmin, xmax);

      for (const auto &e : runs[i].entries) {
        double v = 0.0;
        if (getBranchValue(e, def.key, def.onlyFitOk, v)) h->Fill(v);
      }

      if (h->GetEntries() <= 0) continue;
      if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
      if (ymax < h->GetMaximum()) ymax = h->GetMaximum() * 1.25;
      h->SetLineColor(colors[i % 12]);
      h->SetLineWidth(2);
      h->SetStats(0);

      if (first) {
        h->Draw("HIST");
        first = false;
      } else {
        h->Draw("HIST SAME");
      }
      leg->AddEntry(h, runs[i].label.c_str(), "l");
      drawn.push_back(h);
    }

    for (TH1D *hist : drawn) hist->SetMaximum(ymax);

    if (!first) {
      leg->Draw();
      c->SetGrid();
      c->SaveAs((outDir + "/" + def.key + "_overlay.pdf").c_str());
      c->SaveAs((outDir + "/" + def.key + "_overlay.png").c_str());
    }
  }
}

void drawCellTrend(const std::vector<RunData> &runs, int cellid,
                   const std::string &outDir) {
  if (runs.empty() || cellid < 0) return;

  struct TrendDef {
    std::string key;
    std::string title;
    int color;
    int marker;
  };
  std::vector<TrendDef> defs = {
      {"MPV", "MPV", kRed + 1, 20},
      {"width", "width", kBlue + 1, 21},
      {"FWHM", "FWHM", kMagenta + 1, 20},
      {"chi2perndf", "chi2perndf", kGreen + 2, 21},
      {"gaus_sigma", "gaus_sigma", kCyan + 2, 20},
      {"TotalArea", "TotalArea", kOrange + 7, 21},
  };

  TCanvas *c = new TCanvas(Form("c_cell_%d", cellid), Form("cell %d", cellid),
                           1400, 1100);
  c->Divide(3, 2);

  bool foundAny = false;
  for (size_t idef = 0; idef < defs.size(); ++idef) {
    std::vector<double> x;
    std::vector<double> y;
    x.reserve(runs.size());
    y.reserve(runs.size());

    for (const auto &run : runs) {
      auto it = run.byCellId.find(cellid);
      if (it == run.byCellId.end()) continue;
      double v = 0.0;
      if (!getBranchValue(it->second, defs[idef].key, true, v)) continue;
      x.push_back(run.run);
      y.push_back(v);
    }

    c->cd(static_cast<int>(idef) + 1);
    gPad->SetGrid();

    if (x.empty()) continue;
    foundAny = true;
    TGraph *g = new TGraph(static_cast<int>(x.size()), &x[0], &y[0]);
    setGraphStyle(g, defs[idef].color, defs[idef].marker);
    g->SetTitle(Form("Cell %d %s;Run;%s", cellid, defs[idef].title.c_str(),
                     defs[idef].title.c_str()));
    g->Draw("ALP");
  }

  if (!foundAny) {
    std::cerr << "[compare] cellid " << cellid
              << " not found with fit_ok=1 in loaded runs" << std::endl;
    return;
  }

  c->SaveAs((outDir + Form("/cell_%d_compare_mip.pdf", cellid)).c_str());
  c->SaveAs((outDir + Form("/cell_%d_compare_mip.png", cellid)).c_str());
}

void printSummaryTable(const std::vector<RunData> &runs) {
  std::cout << "\n=== MIP summary by run ===" << std::endl;
  std::cout << "run\tnCell\tMPV(mean)\twidth(mean)\tFWHM(mean)\tchi2perndf(mean)"
               "\tgaus_sigma(mean)\tTotalArea(mean)\tfit_ok[%]"
            << std::endl;

  for (const auto &run : runs) {
    SummaryStats sMpv = extractStats(run, "MPV", true);
    SummaryStats sWidth = extractStats(run, "width", true);
    SummaryStats sFwhm = extractStats(run, "FWHM", true);
    SummaryStats sChi2Ndf = extractStats(run, "chi2perndf", true);
    SummaryStats sGaus = extractStats(run, "gaus_sigma", true);
    SummaryStats sArea = extractStats(run, "TotalArea", true);

    int ok = 0;
    for (const auto &e : run.entries) {
      if (e.fit_ok) ++ok;
    }

    const double frac = 100.0 * ok / std::max(1, static_cast<int>(run.entries.size()));

    std::cout << run.label << "\t" << run.entries.size() << "\t" << sMpv.mean << "\t"
              << sWidth.mean << "\t" << sFwhm.mean << "\t" << sChi2Ndf.mean << "\t"
              << sGaus.mean << "\t" << sArea.mean << "\t" << frac << std::endl;
  }
  std::cout << std::endl;
}

void compare(const char *baseDir = ".",
             const char *runListCsv =
                 "all",
             int focusCellId = -1,
             const char *outDirName = "compare_plots_mip",
             const char *testBeamDir = "EHN1") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  // Accept either numeric run lists (e.g. "22133,22135") or directory-like
  // tokens (e.g. "21987-22133,22135-22135") where each token is a folder
  // name under `baseDir` containing a `mip.root` file. Detect which mode to
  // use by checking for existing directories under `baseDir`.
  std::vector<std::string> tokens = parseTokenList(runListCsv);
  if (tokens.empty()) {
    std::cerr << "[compare] empty run list" << std::endl;
    return;
  }

  bool tokensAreDirs = false;
  for (const auto &t : tokens) {
    const std::string p = std::string(baseDir) + "/" + t + "/mip.root";
    if (gSystem->AccessPathName(p.c_str()) == 0) {
      tokensAreDirs = true;
      break;
    }
  }

  std::vector<int> runNumbers;
  std::vector<RunData> runs;

  if (tokensAreDirs) {
    int idx = 1;
    for (const auto &t : tokens) {
      const std::string path = std::string(baseDir) + "/" + t + "/mip.root";
      RunData data;
      if (loadRunFromFile(path, idx, t, false, data)) {
        std::cout << "[compare] loaded run " << data.label << " from " << data.path
                  << " with " << data.entries.size() << " cells" << std::endl;
        runs.push_back(data);
        ++idx;
      } else {
        std::cout << "[compare] skipping missing/invalid token path " << path << std::endl;
      }
    }
  } else {
    runNumbers = parseRunList(runListCsv);
    if (runNumbers.empty()) {
      std::cerr << "[compare] empty run list" << std::endl;
      return;
    }

    for (int run : runNumbers) {
      RunData data;
      if (loadRun(baseDir, run, data)) {
        std::cout << "[compare] loaded run " << data.label << " from " << data.path
                  << " with " << data.entries.size() << " cells" << std::endl;
        runs.push_back(data);
      }
    }
  }
  std::string outDir = std::string(baseDir) + "/" + outDirName;
  gSystem->mkdir(outDir.c_str(), true);

  if (testBeamDir && std::string(testBeamDir).size() > 0) {
    int pseudoRun = 1;
    if (!runNumbers.empty()) {
      pseudoRun = *std::max_element(runNumbers.begin(), runNumbers.end()) + 1;
    } else if (!runs.empty()) {
      // use next index after last loaded run
      pseudoRun = runs.back().run + 1;
    }
    RunData tbData;
    if (loadTestBeamRun(baseDir, testBeamDir, pseudoRun, tbData)) {
      std::cout << "[compare] loaded TestBeam run " << tbData.label << " from "
                << tbData.path << " with " << tbData.entries.size()
                << " cells (cellid mapped)" << std::endl;
      runs.push_back(tbData);
    } else {
      std::cout << "[compare] TestBeam data not loaded from " << testBeamDir << std::endl;
    }
  }

  if (runs.empty()) {
    std::cerr << "[compare] no runs loaded" << std::endl;
    return;
  }

  printSummaryTable(runs);
  drawSummaryGraphs(runs, outDir);
  drawOverlayHistograms(runs, outDir);
  drawRepresentativeHistograms(runs, outDir);
  drawDeltaFromReferenceByRun(runs, "MPV", "MPV shift from first run", "MPV", outDir,
                              true);
  drawDeltaFromReferenceByRun(runs, "FWHM", "FWHM shift from first run", "FWHM",
                              outDir, true);
  drawDeltaFromReferenceByRun(runs, "chi2perndf", "chi2perndf shift from first run",
                              "chi2perndf", outDir, true);

  // layer-by-layer delta plots for MPV and FWHM
  drawDeltaByLayer(runs, "MPV", outDir, true);
  drawDeltaByLayer(runs, "FWHM", outDir, true);

  // per-layer overlay plots (MPV and FWHM) vs run
  drawLayerOverlays(runs, outDir, true);

  // per-layer overlay plots (delta MPV and delta FWHM from first run)
  drawDeltaByLayerOverlay(runs, "MPV", outDir, true);
  drawDeltaByLayerOverlay(runs, "FWHM", outDir, true);

  if (focusCellId >= 0) {
    drawCellTrend(runs, focusCellId, outDir);
  }

  std::cout << "[compare] output directory: " << outDir << std::endl;
  std::cout << "[compare] done" << std::endl;
}
