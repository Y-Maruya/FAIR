#include <algorithm>
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
#include "TLine.h"
#include "TSystem.h"
#include "TTree.h"

struct PedestalEntry {
  int cellid = -1;
  double highgain_peak = 0.0;
  double lowgain_peak = 0.0;
  double highgain_sigma = 0.0;
  double lowgain_sigma = 0.0;
  int entries_hg = 0;
  int entries_lg = 0;
  int fitStatus_hg = 0;
  int fitStatus_lg = 0;
  int fitOk_hg = 0;
  int fitOk_lg = 0;
  double x_mm = 0.0;
  double y_mm = 0.0;
};

struct RunData {
  int run = 0;
  std::string path;
  std::vector<PedestalEntry> entries;
  std::map<int, PedestalEntry> byCellId;
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
    runs.push_back(std::atoi(token.c_str()));
  }
  return runs;
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
    if (what == "highgain_peak") {
      if (!onlyFitOk || e.fitOk_hg) values.push_back(e.highgain_peak);
    } else if (what == "lowgain_peak") {
      if (!onlyFitOk || e.fitOk_lg) values.push_back(e.lowgain_peak);
    } else if (what == "highgain_sigma") {
      if (!onlyFitOk || e.fitOk_hg) values.push_back(e.highgain_sigma);
    } else if (what == "lowgain_sigma") {
      if (!onlyFitOk || e.fitOk_lg) values.push_back(e.lowgain_sigma);
    }
  }

  return computeStats(values);
}

bool loadRun(const std::string &baseDir, int run, RunData &out) {
  out = RunData();
  out.run = run;
  out.path = baseDir + "/" + std::to_string(run) + "/pedestal.root";

  TFile *file = TFile::Open(out.path.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "[compare] cannot open " << out.path << std::endl;
    if (file) file->Close();
    return false;
  }

  TTree *tree = dynamic_cast<TTree *>(file->Get("pedestal"));
  if (!tree) {
    std::cerr << "[compare] TTree 'pedestal' not found in " << out.path << std::endl;
    file->Close();
    return false;
  }

  PedestalEntry e;
  tree->SetBranchAddress("cellid", &e.cellid);
  tree->SetBranchAddress("highgain_peak", &e.highgain_peak);
  tree->SetBranchAddress("lowgain_peak", &e.lowgain_peak);
  tree->SetBranchAddress("highgain_sigma", &e.highgain_sigma);
  tree->SetBranchAddress("lowgain_sigma", &e.lowgain_sigma);
  tree->SetBranchAddress("entries_hg", &e.entries_hg);
  tree->SetBranchAddress("entries_lg", &e.entries_lg);
  tree->SetBranchAddress("fitStatus_hg", &e.fitStatus_hg);
  tree->SetBranchAddress("fitStatus_lg", &e.fitStatus_lg);
  tree->SetBranchAddress("fitOk_hg", &e.fitOk_hg);
  tree->SetBranchAddress("fitOk_lg", &e.fitOk_lg);
  tree->SetBranchAddress("x_mm", &e.x_mm);
  tree->SetBranchAddress("y_mm", &e.y_mm);

  const Long64_t nEntries = tree->GetEntries();
  out.entries.reserve(nEntries);

  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    out.entries.push_back(e);
    out.byCellId[e.cellid] = e;
  }

  file->Close();
  return true;
}

void setGraphStyle(TGraph *g, int color, int markerStyle) {
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(markerStyle);
  g->SetMarkerSize(1.2);
  g->SetLineWidth(2);
}

bool getBranchValue(const PedestalEntry &e, const std::string &branch, bool onlyFitOk, double &value) {
  if (branch == "highgain_peak") {
    if (onlyFitOk && !e.fitOk_hg) return false;
    value = e.highgain_peak;
    return true;
  }
  if (branch == "lowgain_peak") {
    if (onlyFitOk && !e.fitOk_lg) return false;
    value = e.lowgain_peak;
    return true;
  }
  if (branch == "highgain_sigma") {
    if (onlyFitOk && !e.fitOk_hg) return false;
    value = e.highgain_sigma;
    return true;
  }
  if (branch == "lowgain_sigma") {
    if (onlyFitOk && !e.fitOk_lg) return false;
    value = e.lowgain_sigma;
    return true;
  }
  return false;
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

  std::sort(ranked.begin(), ranked.end(), [](const CellShiftInfo &a, const CellShiftInfo &b) {
    if (a.spread != b.spread) return a.spread > b.spread;
    return a.cellid < b.cellid;
  });
  return ranked;
}

TH1 *loadPedestalHistogram(const std::string &filePath, const std::string &gain, int cellid) {
  TFile *file = TFile::Open(filePath.c_str(), "READ");
  if (!file || file->IsZombie()) {
    if (file) file->Close();
    return nullptr;
  }

  std::vector<std::string> candidates = {
      Form("Pedestal/%s/hPed%s_%d", gain.c_str(), gain.c_str(), cellid),
      Form("%s/hPed%s_%d", gain.c_str(), gain.c_str(), cellid),
      Form("hPed%s_%d", gain.c_str(), cellid),
  };

  TH1 *hist = nullptr;
  for (const auto &name : candidates) {
    hist = dynamic_cast<TH1 *>(file->Get(name.c_str()));
    if (hist) break;
  }

  TH1 *cloned = nullptr;
  if (hist) {
    static long long cloneCounter = 0;
    ++cloneCounter;
    cloned = dynamic_cast<TH1 *>(hist->Clone(Form("%s_clone_%d_%lld", gain.c_str(), cellid, cloneCounter)));
    if (cloned) cloned->SetDirectory(nullptr);
  }

  file->Close();
  return cloned;
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
    ofs << info.cellid << " " << info.spread << " " << info.refValue << " " << info.meanValue
        << " " << info.x_mm << " " << info.y_mm << "\n";
  }

  ofs << "\n# stable channels\n";
  ofs << "# cellid spread refValue meanValue x_mm y_mm\n";
  for (const auto &info : stable) {
    ofs << info.cellid << " " << info.spread << " " << info.refValue << " " << info.meanValue
        << " " << info.x_mm << " " << info.y_mm << "\n";
  }
}

void drawRepresentativeHistogramSet(const std::vector<RunData> &runs,
                                    const std::vector<CellShiftInfo> &cells,
                                    const std::string &gain,
                                    const std::string &label,
                                    const std::string &outDir,
                                    const std::string &canvasTag) {
  if (cells.empty()) return;

  const int n = static_cast<int>(cells.size());
  const int nCols = 2;
  const int nRows = (n + nCols - 1) / nCols;
  TCanvas *c = new TCanvas(Form("c_%s_%s", gain.c_str(), canvasTag.c_str()),
                           Form("%s %s", gain.c_str(), label.c_str()), 1400, 450 * nRows);
  c->Divide(nCols, nRows);

  int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kRed + 3, kBlue + 3, kGreen + 4, kMagenta + 3, kOrange + 9, kCyan + 4};

  for (int i = 0; i < n; ++i) {
    c->cd(i + 1);
    gPad->SetGrid();

    TLegend *leg = new TLegend(0.58, 0.62, 0.88, 0.88);
    bool first = true;
    double maxY = 0.0;
    std::vector<TH1 *> drawn;

    for (size_t irun = 0; irun < runs.size(); ++irun) {
      TH1 *h = loadPedestalHistogram(runs[irun].path, gain, cells[i].cellid);
      if (!h) continue;

      if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
      h->GetXaxis()->SetRangeUser(cells[i].meanValue - 50, cells[i].meanValue + 50);
      h->SetLineColor(colors[irun % 12]);
      h->SetLineWidth(2);
      h->SetStats(0);
      h->SetTitle(Form("%s cell %d (spread = %.3f);ADC counts;Normalized entries",
                       gain.c_str(), cells[i].cellid, cells[i].spread));

      maxY = std::max(maxY, h->GetMaximum());
      if (first) {
        h->Draw("HIST");
        first = false;
      } else {
        h->Draw("HIST SAME");
      }
      leg->AddEntry(h, Form("Run %d", runs[irun].run), "l");
      drawn.push_back(h);
    }

    if (!drawn.empty()) {
      drawn.front()->SetMaximum(maxY * 1.25);
      leg->SetHeader(Form("cellid=%d, (x,y)=(%.1f, %.1f)", cells[i].cellid, cells[i].x_mm, cells[i].y_mm), "C");
      leg->Draw();
    }
  }

  c->SaveAs((outDir + "/" + gain + "_" + canvasTag + ".pdf").c_str());
  c->SaveAs((outDir + "/" + gain + "_" + canvasTag + ".png").c_str());
}

void drawRepresentativeHistograms(const std::vector<RunData> &runs,
                                  const std::string &outDir,
                                  int nRepresentative = 8) {
  struct GainDef {
    std::string branch;
    std::string gain;
  };

  std::vector<GainDef> defs = {
      {"highgain_peak", "HG"},
      {"lowgain_peak", "LG"},
  };

  for (const auto &def : defs) {
    std::vector<CellShiftInfo> ranked = rankCellsBySpread(runs, def.branch, true);
    if (ranked.empty()) continue;

    const int nTake = std::min(nRepresentative, static_cast<int>(ranked.size()));
    std::vector<CellShiftInfo> shifted(ranked.begin(), ranked.begin() + nTake);

    std::vector<CellShiftInfo> stable;
    for (auto it = ranked.rbegin(); it != ranked.rend() && static_cast<int>(stable.size()) < nTake; ++it) {
      stable.push_back(*it);
    }
    std::reverse(stable.begin(), stable.end());

    writeRepresentativeSummary(outDir + "/" + def.branch + "_representative_channels.txt",
                               shifted, stable, def.branch);
    drawRepresentativeHistogramSet(runs, shifted, def.gain, "shifted channels", outDir,
                                   def.branch + "_shifted_channels");
    drawRepresentativeHistogramSet(runs, stable, def.gain, "stable channels", outDir,
                                   def.branch + "_stable_channels");
  }
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

      const PedestalEntry &a = kv.second;
      const PedestalEntry &b = it->second;
      double va = 0.0;
      double vb = 0.0;
      bool ok = true;

      if (branch == "highgain_peak") {
        ok = !onlyFitOk || (a.fitOk_hg && b.fitOk_hg);
        va = a.highgain_peak;
        vb = b.highgain_peak;
      } else if (branch == "lowgain_peak") {
        ok = !onlyFitOk || (a.fitOk_lg && b.fitOk_lg);
        va = a.lowgain_peak;
        vb = b.lowgain_peak;
      } else if (branch == "highgain_sigma") {
        ok = !onlyFitOk || (a.fitOk_hg && b.fitOk_hg);
        va = a.highgain_sigma;
        vb = b.highgain_sigma;
      } else if (branch == "lowgain_sigma") {
        ok = !onlyFitOk || (a.fitOk_lg && b.fitOk_lg);
        va = a.lowgain_sigma;
        vb = b.lowgain_sigma;
      }

      if (ok) deltas.push_back(vb - va);
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

  TH1D *h = new TH1D(name, title, 80, xmin - 0.1 * std::abs(xmin), xmax + 0.1 * std::abs(xmax));
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

    const PedestalEntry &a = kv.second;
    const PedestalEntry &b = it->second;
    double va = 0.0;
    double vb = 0.0;
    bool ok = true;

    if (branch == "highgain_peak") {
      ok = !onlyFitOk || (a.fitOk_hg && b.fitOk_hg);
      va = a.highgain_peak;
      vb = b.highgain_peak;
    } else if (branch == "lowgain_peak") {
      ok = !onlyFitOk || (a.fitOk_lg && b.fitOk_lg);
      va = a.lowgain_peak;
      vb = b.lowgain_peak;
    } else if (branch == "highgain_sigma") {
      ok = !onlyFitOk || (a.fitOk_hg && b.fitOk_hg);
      va = a.highgain_sigma;
      vb = b.highgain_sigma;
    } else if (branch == "lowgain_sigma") {
      ok = !onlyFitOk || (a.fitOk_lg && b.fitOk_lg);
      va = a.lowgain_sigma;
      vb = b.lowgain_sigma;
    }

    if (ok) deltas.push_back(vb - va);
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
  std::vector<std::vector<double>> allDeltas(runs.size());
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
  TCanvas *c = new TCanvas(Form("c_%s_delta_ref", outTag.c_str()), title.c_str(), 1000, 700);
  c->SetGrid();

  TLegend *leg = new TLegend(0.62, 0.64, 0.88, 0.88);
  int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kRed + 3, kBlue + 3, kGreen + 4, kMagenta + 3, kOrange + 9, kCyan + 4};
  bool first = true;
  std::vector<TH1D *> drawn;
  for (size_t i = 1; i < runs.size(); ++i) {
    if (allDeltas[i].empty()) continue;

    TH1D *h = new TH1D(Form("h_%s_delta_ref_%d", outTag.c_str(), runs[i].run),
                       Form("%s (ref run %d);#Delta ADC;Normalized entries", title.c_str(), ref.run),
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
    leg->AddEntry(h, Form("Run %d - %d", runs[i].run, ref.run), "l");
  }
  for (TH1D *hist : drawn) {
    hist->SetMaximum(ymax);
  }

  if (!first) {
    leg->Draw();
    c->SaveAs((outDir + "/" + outTag + "_delta_from_first_run_overlay.pdf").c_str());
    c->SaveAs((outDir + "/" + outTag + "_delta_from_first_run_overlay.png").c_str());
  }

}

void drawSummaryGraphs(const std::vector<RunData> &runs, const std::string &outDir) {
  const int n = static_cast<int>(runs.size());
  if (n == 0) return;

  std::vector<double> x(n), hgPeak(n), lgPeak(n), hgSigma(n), lgSigma(n), fitHg(n), fitLg(n);
  for (int i = 0; i < n; ++i) {
    x[i] = runs[i].run;

    SummaryStats sHgPeak = extractStats(runs[i], "highgain_peak", true);
    SummaryStats sLgPeak = extractStats(runs[i], "lowgain_peak", true);
    SummaryStats sHgSigma = extractStats(runs[i], "highgain_sigma", true);
    SummaryStats sLgSigma = extractStats(runs[i], "lowgain_sigma", true);

    hgPeak[i] = sHgPeak.mean;
    lgPeak[i] = sLgPeak.mean;
    hgSigma[i] = sHgSigma.mean;
    lgSigma[i] = sLgSigma.mean;

    int okHg = 0;
    int okLg = 0;
    for (const auto &e : runs[i].entries) {
      if (e.fitOk_hg) ++okHg;
      if (e.fitOk_lg) ++okLg;
    }
    fitHg[i] = 100.0 * okHg / std::max(1, static_cast<int>(runs[i].entries.size()));
    fitLg[i] = 100.0 * okLg / std::max(1, static_cast<int>(runs[i].entries.size()));
  }

  TCanvas *c = new TCanvas("c_summary", "pedestal summary", 1200, 900);
  c->Divide(2, 2);

  c->cd(1);
  TMultiGraph *mgPeak = new TMultiGraph();
  TGraph *gHgPeak = new TGraph(n, &x[0], &hgPeak[0]);
  TGraph *gLgPeak = new TGraph(n, &x[0], &lgPeak[0]);
  gHgPeak->SetTitle("HG peak");
  gLgPeak->SetTitle("LG peak");
  setGraphStyle(gHgPeak, kRed + 1, 20);
  setGraphStyle(gLgPeak, kBlue + 1, 21);
  mgPeak->Add(gHgPeak, "LP");
  mgPeak->Add(gLgPeak, "LP");
  mgPeak->SetTitle("Mean pedestal peak vs run;Run;ADC counts");
  mgPeak->Draw("A");
  {
    TLegend *leg = new TLegend(0.22, 0.42, 0.48, 0.68);
    leg->AddEntry(gHgPeak, "HG peak", "lp");
    leg->AddEntry(gLgPeak, "LG peak", "lp");
    leg->Draw();
  }

  c->cd(2);
  TMultiGraph *mgSigma = new TMultiGraph();
  TGraph *gHgSigma = new TGraph(n, &x[0], &hgSigma[0]);
  TGraph *gLgSigma = new TGraph(n, &x[0], &lgSigma[0]);
  gHgSigma->SetTitle("HG sigma");
  gLgSigma->SetTitle("LG sigma");
  setGraphStyle(gHgSigma, kMagenta + 1, 20);
  setGraphStyle(gLgSigma, kCyan + 2, 21);
  mgSigma->Add(gHgSigma, "LP");
  mgSigma->Add(gLgSigma, "LP");
  mgSigma->SetTitle("Mean pedestal sigma vs run;Run;ADC counts");
  mgSigma->Draw("A");
  {
    TLegend *leg = new TLegend(0.22, 0.32, 0.48, 0.58);
    leg->AddEntry(gHgSigma, "HG sigma", "lp");
    leg->AddEntry(gLgSigma, "LG sigma", "lp");
    leg->Draw();
  }

  c->cd(3);
  TMultiGraph *mgFit = new TMultiGraph();
  TGraph *gFitHg = new TGraph(n, &x[0], &fitHg[0]);
  TGraph *gFitLg = new TGraph(n, &x[0], &fitLg[0]);
  gFitHg->SetTitle("HG fit OK");
  gFitLg->SetTitle("LG fit OK");
  setGraphStyle(gFitHg, kOrange + 7, 20);
  setGraphStyle(gFitLg, kGreen + 2, 21);
  mgFit->Add(gFitHg, "LP");
  mgFit->Add(gFitLg, "LP");
  mgFit->SetTitle("Fit success rate vs run;Run;Fit OK [%]");
  mgFit->Draw("A");
  {
    TLegend *leg = new TLegend(0.12, 0.40, 0.38, 0.66);
    leg->AddEntry(gFitHg, "HG fit OK", "lp");
    leg->AddEntry(gFitLg, "LG fit OK", "lp");
    leg->Draw();
  }

  c->cd(4);
  gPad->SetGrid();
  TH1D *hDeltaHg = makeDeltaHist("hDeltaHg", Form("#Delta HG peak relative to run %d;#Delta ADC;Cells", runs.front().run), runs, "highgain_peak", true);
  TH1D *hDeltaLg = makeDeltaHist("hDeltaLg", Form("#Delta LG peak relative to run %d;#Delta ADC;Cells", runs.front().run), runs, "lowgain_peak", true);
  if (hDeltaHg) {
    hDeltaHg->SetLineColor(kRed + 1);
    hDeltaHg->Draw("HIST");
  }
  if (hDeltaLg) {
    hDeltaLg->SetLineColor(kBlue + 1);
    if (hDeltaHg) hDeltaLg->Draw("HIST SAME");
    else hDeltaLg->Draw("HIST");
  }
  if (hDeltaHg || hDeltaLg) gPad->BuildLegend();

  c->SaveAs((outDir + "/summary_compare.pdf").c_str());
  c->SaveAs((outDir + "/summary_compare.png").c_str());
}

void drawOverlayHistograms(const std::vector<RunData> &runs, const std::string &outDir) {
  struct PlotDef {
    std::string key;
    std::string title;
    std::string xTitle;
    bool onlyFitOk;
  };

  std::vector<PlotDef> defs = {
      {"highgain_peak", "High-gain pedestal peak", "ADC counts", true},
      {"lowgain_peak", "Low-gain pedestal peak", "ADC counts", true},
      {"highgain_sigma", "High-gain pedestal sigma", "ADC counts", true},
      {"lowgain_sigma", "Low-gain pedestal sigma", "ADC counts", true},
  };

  int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kRed + 3, kBlue + 3, kGreen + 4, kMagenta + 3, kOrange + 9, kCyan + 4};

  for (const auto &def : defs) {
    std::vector<TH1D *> hists;
    double globalMin = 1e30;
    double globalMax = -1e30;

    for (size_t i = 0; i < runs.size(); ++i) {
      std::vector<double> values;
      values.reserve(runs[i].entries.size());
      for (const auto &e : runs[i].entries) {
        if (def.key == "highgain_peak") {
          if (!def.onlyFitOk || e.fitOk_hg) values.push_back(e.highgain_peak);
        } else if (def.key == "lowgain_peak") {
          if (!def.onlyFitOk || e.fitOk_lg) values.push_back(e.lowgain_peak);
        } else if (def.key == "highgain_sigma") {
          if (!def.onlyFitOk || e.fitOk_hg) values.push_back(e.highgain_sigma);
        } else if (def.key == "lowgain_sigma") {
          if (!def.onlyFitOk || e.fitOk_lg) values.push_back(e.lowgain_sigma);
        }
      }

      if (values.empty()) continue;

      SummaryStats s = computeStats(values);
      globalMin = std::min(globalMin, s.min);
      globalMax = std::max(globalMax, s.max);
    }

    if (globalMin >= globalMax) {
      globalMin -= 1.0;
      globalMax += 1.0;
    }
    int globalBins = globalMax - globalMin + 1;

    TCanvas *c = new TCanvas(Form("c_%s", def.key.c_str()), def.title.c_str(), 1000, 700);
    TLegend *leg = new TLegend(0.65, 0.70, 0.88, 0.88);
    bool first = true;

    for (size_t i = 0; i < runs.size(); ++i) {
      TH1D *h = new TH1D(Form("h_%s_%d", def.key.c_str(), runs[i].run),
                         Form("%s; %s; Normalized entries", def.title.c_str(), def.xTitle.c_str()),
                         globalBins, globalMin, globalMax);

      for (const auto &e : runs[i].entries) {
        if (def.key == "highgain_peak") {
          if (!def.onlyFitOk || e.fitOk_hg) h->Fill(e.highgain_peak);
        } else if (def.key == "lowgain_peak") {
          if (!def.onlyFitOk || e.fitOk_lg) h->Fill(e.lowgain_peak);
        } else if (def.key == "highgain_sigma") {
          if (!def.onlyFitOk || e.fitOk_hg) h->Fill(e.highgain_sigma);
        } else if (def.key == "lowgain_sigma") {
          if (!def.onlyFitOk || e.fitOk_lg) h->Fill(e.lowgain_sigma);
        }
      }

      if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
      h->SetLineColor(colors[i % 12]);
      h->SetLineWidth(2);
      h->SetStats(0);
      if (first) {
        h->Draw("HIST");
        first = false;
      } else {
        h->Draw("HIST SAME");
      }
      leg->AddEntry(h, Form("Run %d", runs[i].run), "l");
      hists.push_back(h);
    }

    leg->Draw();
    c->SetGrid();
    c->SaveAs((outDir + "/" + def.key + "_overlay.pdf").c_str());
    c->SaveAs((outDir + "/" + def.key + "_overlay.png").c_str());
  }
}

void drawCellTrend(const std::vector<RunData> &runs, int cellid, const std::string &outDir) {
  if (runs.empty() || cellid < 0) return;

  std::vector<double> x, hgPeak, lgPeak, hgSigma, lgSigma;
  for (const auto &run : runs) {
    auto it = run.byCellId.find(cellid);
    if (it == run.byCellId.end()) continue;
    x.push_back(run.run);
    hgPeak.push_back(it->second.highgain_peak);
    lgPeak.push_back(it->second.lowgain_peak);
    hgSigma.push_back(it->second.highgain_sigma);
    lgSigma.push_back(it->second.lowgain_sigma);
  }

  if (x.empty()) {
    std::cerr << "[compare] cellid " << cellid << " not found in loaded runs" << std::endl;
    return;
  }

  TCanvas *c = new TCanvas(Form("c_cell_%d", cellid), Form("cell %d", cellid), 1200, 900);
  c->Divide(2, 2);

  c->cd(1);
  TGraph *gHgPeak = new TGraph(x.size(), &x[0], &hgPeak[0]);
  setGraphStyle(gHgPeak, kRed + 1, 20);
  gHgPeak->SetTitle(Form("Cell %d HG peak;Run;ADC counts", cellid));
  gHgPeak->Draw("ALP");

  c->cd(2);
  TGraph *gLgPeak = new TGraph(x.size(), &x[0], &lgPeak[0]);
  setGraphStyle(gLgPeak, kBlue + 1, 21);
  gLgPeak->SetTitle(Form("Cell %d LG peak;Run;ADC counts", cellid));
  gLgPeak->Draw("ALP");

  c->cd(3);
  TGraph *gHgSigma = new TGraph(x.size(), &x[0], &hgSigma[0]);
  setGraphStyle(gHgSigma, kMagenta + 1, 20);
  gHgSigma->SetTitle(Form("Cell %d HG sigma;Run;ADC counts", cellid));
  gHgSigma->Draw("ALP");

  c->cd(4);
  TGraph *gLgSigma = new TGraph(x.size(), &x[0], &lgSigma[0]);
  setGraphStyle(gLgSigma, kCyan + 2, 21);
  gLgSigma->SetTitle(Form("Cell %d LG sigma;Run;ADC counts", cellid));
  gLgSigma->Draw("ALP");

  c->SaveAs((outDir + Form("/cell_%d_compare.pdf", cellid)).c_str());
  c->SaveAs((outDir + Form("/cell_%d_compare.png", cellid)).c_str());
}

void printSummaryTable(const std::vector<RunData> &runs) {
  std::cout << "\n=== Pedestal summary by run ===" << std::endl;
  std::cout << "run\tnCell\tHGpeak(mean)\tLGpeak(mean)\tHGsigma(mean)\tLGsigma(mean)\tFitOK_HG[%]\tFitOK_LG[%]" << std::endl;

  for (const auto &run : runs) {
    SummaryStats sHgPeak = extractStats(run, "highgain_peak", true);
    SummaryStats sLgPeak = extractStats(run, "lowgain_peak", true);
    SummaryStats sHgSigma = extractStats(run, "highgain_sigma", true);
    SummaryStats sLgSigma = extractStats(run, "lowgain_sigma", true);

    int okHg = 0;
    int okLg = 0;
    for (const auto &e : run.entries) {
      if (e.fitOk_hg) ++okHg;
      if (e.fitOk_lg) ++okLg;
    }

    double fracHg = 100.0 * okHg / std::max(1, static_cast<int>(run.entries.size()));
    double fracLg = 100.0 * okLg / std::max(1, static_cast<int>(run.entries.size()));

    std::cout << run.run << "\t"
              << run.entries.size() << "\t"
              << sHgPeak.mean << "\t"
              << sLgPeak.mean << "\t"
              << sHgSigma.mean << "\t"
              << sLgSigma.mean << "\t"
              << fracHg << "\t"
              << fracLg << std::endl;
  }
  std::cout << std::endl;
}

void compare(const char *baseDir = ".",
             const char *runListCsv = "22074,22140,22160,22168,22206,22249,22287,22324,22334",
             int focusCellId = -1,
             const char *outDirName = "compare_plots") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  std::vector<int> runNumbers = parseRunList(runListCsv);
  if (runNumbers.empty()) {
    std::cerr << "[compare] empty run list" << std::endl;
    return;
  }

  std::string outDir = std::string(baseDir) + "/" + outDirName;
  gSystem->mkdir(outDir.c_str(), true);

  std::vector<RunData> runs;
  for (int run : runNumbers) {
    RunData data;
    if (loadRun(baseDir, run, data)) {
      std::cout << "[compare] loaded run " << run << " from " << data.path
                << " with " << data.entries.size() << " cells" << std::endl;
      runs.push_back(data);
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
  drawDeltaFromReferenceByRun(runs, "highgain_peak", "HG peak shift from first run",
                              "highgain_peak", outDir, true);
  drawDeltaFromReferenceByRun(runs, "lowgain_peak", "LG peak shift from first run",
                              "lowgain_peak", outDir, true);

  if (focusCellId >= 0) {
    drawCellTrend(runs, focusCellId, outDir);
  }

  std::cout << "[compare] output directory: " << outDir << std::endl;
  std::cout << "[compare] done" << std::endl;
}
