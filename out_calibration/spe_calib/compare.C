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

struct SpeEntry {
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;
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
  int npad = 0;
  double gain = 0.0;
  double gainErr = 0.0;
  double snr = 0.0;
  int kp = 0;
  double k0 = 0.0;
  double peakAmp = 0.0;
  double noiseMed = 0.0;
};

struct RunData {
  int run = 0;
  std::string label;
  std::string path;
  std::vector<SpeEntry> entries;
  std::map<int, SpeEntry> byCellId;
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

void setGraphStyle(TGraph *g, int color, int markerStyle) {
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(markerStyle);
  g->SetMarkerSize(1.1);
  g->SetLineWidth(2);
}

bool getBranchValue(const SpeEntry &e, const std::string &branch, bool onlyFitOk,
                    double &value) {
  if (onlyFitOk && !e.fit_ok) return false;

  if (branch == "gain") {
    value = e.gain;
    return true;
  }
  if (branch == "gainErr") {
    value = e.gainErr;
    return true;
  }
  if (branch == "snr") {
    value = e.snr;
    return true;
  }
  if (branch == "noiseMed") {
    value = e.noiseMed;
    return true;
  }
  if (branch == "peakAmp") {
    value = e.peakAmp;
    return true;
  }
  if (branch == "k0") {
    value = e.k0;
    return true;
  }
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

  TTree *tree = dynamic_cast<TTree *>(file->Get("spe"));
  if (!tree) {
    std::cerr << "[compare] TTree 'spe' not found in " << out.path << std::endl;
    file->Close();
    return false;
  }

  SpeEntry e;
  if (!setBranchIfExists(tree, "cellid", &e.cellid) ||
      !setBranchIfExists(tree, "gain", &e.gain) ||
      !setBranchIfExists(tree, "snr", &e.snr)) {
    std::cerr << "[compare] missing mandatory SPE branches (cellid, gain, snr) in " << out.path
              << std::endl;
    file->Close();
    return false;
  }

  // Optional branches
  setBranchIfExists(tree, "layer", &e.layer);
  setBranchIfExists(tree, "chip", &e.chip);
  setBranchIfExists(tree, "channel", &e.channel);
  setBranchIfExists(tree, "MPV", &e.MPV);
  setBranchIfExists(tree, "width", &e.width);
  setBranchIfExists(tree, "TotalArea", &e.TotalArea);
  setBranchIfExists(tree, "entries", &e.entries);
  setBranchIfExists(tree, "gaus_sigma", &e.gaus_sigma);
  setBranchIfExists(tree, "max_x", &e.max_x);
  setBranchIfExists(tree, "FWHM", &e.FWHM);
  setBranchIfExists(tree, "chi2", &e.chi2);
  setBranchIfExists(tree, "ndf", &e.ndf);
  setBranchIfExists(tree, "chi2perndf", &e.chi2perndf);
  setBranchIfExists(tree, "fit_status", &e.fit_status);
  out.hasFitFlag = setBranchIfExists(tree, "fit_ok", &e.fit_ok);
  bool hasX = setBranchIfExists(tree, "x_mm", &e.x_mm);
  bool hasY = setBranchIfExists(tree, "y_mm", &e.y_mm);
  out.hasCoordinates = hasX && hasY;
  setBranchIfExists(tree, "npad", &e.npad);
  setBranchIfExists(tree, "gainErr", &e.gainErr);
  setBranchIfExists(tree, "kp", &e.kp);
  setBranchIfExists(tree, "k0", &e.k0);
  setBranchIfExists(tree, "peakAmp", &e.peakAmp);
  setBranchIfExists(tree, "noiseMed", &e.noiseMed);

  const Long64_t nEntries = tree->GetEntries();
  out.entries.reserve(nEntries);

  for (Long64_t i = 0; i < nEntries; ++i) {
    e.fit_ok = 1;
    e.x_mm = 0.0;
    e.y_mm = 0.0;
    tree->GetEntry(i);
    SpeEntry entry = e;
    out.entries.push_back(entry);
    out.byCellId[entry.cellid] = entry;
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
  const std::string path = baseDir + "/" + std::to_string(run) + "/spe_analysis.root";
  return loadRunFromFile(path, run, std::to_string(run), false, out);
}

bool loadTestBeamRun(const std::string &baseDir,
                     const std::string &testBeamDir,
                     int pseudoRun,
                     RunData &out) {
  const std::string path = baseDir + "/" + testBeamDir + "/spe_analysis.root";
  return loadRunFromFile(path, pseudoRun, testBeamDir, false, out);
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
  if (outTag == "gain") {
    c->SetLogy();
  }
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

void drawSummaryGraphs(const std::vector<RunData> &runs, const std::string &outDir) {
  const int n = static_cast<int>(runs.size());
  if (n == 0) return;

  std::vector<double> x(n), gain(n), gainErr(n), snr(n), noiseMed(n), fitok(n), peakAmp(n);
  for (int i = 0; i < n; ++i) {
    x[i] = runs[i].run;

    gain[i] = extractStats(runs[i], "gain", true).mean;
    gainErr[i] = extractStats(runs[i], "gainErr", true).mean;
    snr[i] = extractStats(runs[i], "snr", true).mean;
    noiseMed[i] = extractStats(runs[i], "noiseMed", true).mean;
    peakAmp[i] = extractStats(runs[i], "peakAmp", true).mean;

    int ok = 0;
    for (const auto &e : runs[i].entries) {
      if (e.fit_ok) ++ok;
    }
    fitok[i] = 100.0 * ok / std::max(1, static_cast<int>(runs[i].entries.size()));
  }

  TCanvas *c = new TCanvas("c_summary", "SPE summary", 1200, 900);
  c->Divide(2, 2);

  c->cd(1);
  TMultiGraph *mgGain = new TMultiGraph();
  TGraph *gGain = new TGraph(n, &x[0], &gain[0]);
  TGraph *gGainErr = new TGraph(n, &x[0], &gainErr[0]);
  setGraphStyle(gGain, kRed + 1, 20);
  setGraphStyle(gGainErr, kBlue + 1, 21);
  mgGain->Add(gGain, "LP");
  mgGain->Add(gGainErr, "LP");
  mgGain->SetTitle("Mean gain/gainErr vs run;Run;Gain");
  mgGain->Draw("A");
  {
    TLegend *leg = new TLegend(0.16, 0.42, 0.42, 0.68);
    leg->AddEntry(gGain, "gain", "lp");
    leg->AddEntry(gGainErr, "gainErr", "lp");
    leg->Draw();
  }

  c->cd(2);
  TMultiGraph *mgNoise = new TMultiGraph();
  TGraph *gSnr = new TGraph(n, &x[0], &snr[0]);
  TGraph *gNoiseMed = new TGraph(n, &x[0], &noiseMed[0]);
  setGraphStyle(gSnr, kMagenta + 1, 20);
  setGraphStyle(gNoiseMed, kCyan + 2, 21);
  mgNoise->Add(gSnr, "LP");
  mgNoise->Add(gNoiseMed, "LP");
  mgNoise->SetTitle("Mean SNR/noiseMed vs run;Run;Value");
  mgNoise->Draw("A");
  {
    TLegend *leg = new TLegend(0.16, 0.42, 0.42, 0.68);
    leg->AddEntry(gSnr, "SNR", "lp");
    leg->AddEntry(gNoiseMed, "noiseMed", "lp");
    leg->Draw();
  }
    {
      // Calculate average entries for SNR legend
      SummaryStats sEntries = extractStats(runs[0], "entries", false);
      TLegend *leg2 = new TLegend(0.16, 0.12, 0.42, 0.28);
      leg2->AddEntry((TObject*)0, Form("Avg entries: %.0f", sEntries.mean), "");
      leg2->Draw();
    }
  c->cd(3);
  TMultiGraph *mgQuality = new TMultiGraph();
  TGraph *gPeakAmp = new TGraph(n, &x[0], &peakAmp[0]);
  TGraph *gFitOk = new TGraph(n, &x[0], &fitok[0]);
  setGraphStyle(gPeakAmp, kGreen + 2, 20);
  setGraphStyle(gFitOk, kOrange + 7, 21);
  mgQuality->Add(gPeakAmp, "LP");
  mgQuality->Add(gFitOk, "LP");
  mgQuality->SetTitle("Peak amplitude and fit quality vs run;Run;Value");
  mgQuality->Draw("A");
  {
    TLegend *leg = new TLegend(0.16, 0.42, 0.42, 0.68);
    leg->AddEntry(gPeakAmp, "peakAmp (mean)", "lp");
    leg->AddEntry(gFitOk, "fit_ok [%]", "lp");
    leg->Draw();
  }

  c->cd(4);
  gPad->SetGrid();
  TH1D *hDeltaGain = makeDeltaHist("hDeltaGain",
                                   Form("#Delta gain relative to run %d;#Delta gain;Cells",
                                        runs.front().run),
                                   runs, "gain", true);
  TH1D *hDeltaSnr = makeDeltaHist("hDeltaSnr",
                                  Form("#Delta SNR relative to run %d;#Delta SNR;Cells",
                                       runs.front().run),
                                  runs, "snr", true);

  if (hDeltaGain) {
    hDeltaGain->SetLineColor(kRed + 1);
    hDeltaGain->Draw("HIST");
  }
  if (hDeltaSnr) {
    hDeltaSnr->SetLineColor(kBlue + 1);
    if (hDeltaGain)
      hDeltaSnr->Draw("HIST SAME");
    else
      hDeltaSnr->Draw("HIST");
  }
  if (hDeltaGain || hDeltaSnr) gPad->BuildLegend();

  c->SaveAs((outDir + "/summary_compare_spe.pdf").c_str());
  c->SaveAs((outDir + "/summary_compare_spe.png").c_str());
}

void drawOverlayHistograms(const std::vector<RunData> &runs, const std::string &outDir) {
  struct PlotDef {
    std::string key;
    std::string title;
    std::string xTitle;
    bool onlyFitOk;
  };

  std::vector<PlotDef> defs = {
      {"gain", "SPE gain", "Gain (mV/PE)", true},
      {"gainErr", "SPE gain error", "Gain error (mV/PE)", true},
      {"snr", "SPE SNR", "SNR", true},
      {"noiseMed", "SPE noise median", "Noise (mV)", true},
      {"peakAmp", "SPE peak amplitude", "Peak amplitude (ADC)", true},
      {"k0", "SPE k0", "k0", true},
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
      {"gain", "gain", kRed + 1, 20},
      {"gainErr", "gainErr", kBlue + 1, 21},
      {"snr", "SNR", kMagenta + 1, 20},
      {"noiseMed", "noiseMed", kGreen + 2, 21},
      {"peakAmp", "peakAmp", kCyan + 2, 20},
      {"k0", "k0", kOrange + 7, 21},
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

  c->SaveAs((outDir + Form("/cell_%d_compare_spe.pdf", cellid)).c_str());
  c->SaveAs((outDir + Form("/cell_%d_compare_spe.png", cellid)).c_str());
}

void printSummaryTable(const std::vector<RunData> &runs) {
  std::cout << "\n=== SPE summary by run ===" << std::endl;
  std::cout << "run\tnCell\tgain(mean)\tgainErr(mean)\tSNR(mean)\tnoiseMed(mean)"
               "\tpeakAmp(mean)\tk0(mean)\tfit_ok[%]"
            << std::endl;

  for (const auto &run : runs) {
    SummaryStats sGain = extractStats(run, "gain", true);
    SummaryStats sGainErr = extractStats(run, "gainErr", true);
    SummaryStats sSnr = extractStats(run, "snr", true);
    SummaryStats sNoiseMed = extractStats(run, "noiseMed", true);
    SummaryStats sPeakAmp = extractStats(run, "peakAmp", true);
    SummaryStats sK0 = extractStats(run, "k0", true);

    int ok = 0;
    for (const auto &e : run.entries) {
      if (e.fit_ok) ++ok;
    }

    const double frac = 100.0 * ok / std::max(1, static_cast<int>(run.entries.size()));

    std::cout << run.label << "\t" << run.entries.size() << "\t" << sGain.mean << "\t"
              << sGainErr.mean << "\t" << sSnr.mean << "\t" << sNoiseMed.mean << "\t"
              << sPeakAmp.mean << "\t" << sK0.mean << "\t" << frac << std::endl;
  }
  std::cout << std::endl;
}

void compare(const char *baseDir = ".",
             const char *runListCsv =
                 "21987-22133,22135-22135,22137-22142,22159-22163,22165-22198,22205-22251,22268-22268",
             int focusCellId = -1,
             const char *outDirName = "compare_plots_spe",
             const char *testBeamDir = "") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  // Accept either numeric run lists (e.g. "22133,22135") or directory-like
  // tokens (e.g. "21987-22133,22135-22135") where each token is a folder
  // name under `baseDir` containing a `spe_analysis.root` file. Detect which mode to
  // use by checking for existing directories under `baseDir`.
  std::vector<std::string> tokens = parseTokenList(runListCsv);
  if (tokens.empty()) {
    std::cerr << "[compare] empty run list" << std::endl;
    return;
  }

  bool tokensAreDirs = false;
  for (const auto &t : tokens) {
    const std::string p = std::string(baseDir) + "/" + t + "/spe_analysis.root";
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
      const std::string path = std::string(baseDir) + "/" + t + "/spe_analysis.root";
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

  if (runs.empty()) {
    std::cerr << "[compare] no runs loaded" << std::endl;
    return;
  }

  printSummaryTable(runs);
  drawSummaryGraphs(runs, outDir);
  drawOverlayHistograms(runs, outDir);
  drawDeltaFromReferenceByRun(runs, "gain", "gain shift from first run", "gain", outDir,
                              true);
  drawDeltaFromReferenceByRun(runs, "snr", "SNR shift from first run", "snr",
                              outDir, true);
  drawDeltaFromReferenceByRun(runs, "noiseMed", "noiseMed shift from first run",
                              "noiseMed", outDir, true);

  if (focusCellId >= 0) {
    drawCellTrend(runs, focusCellId, outDir);
  }

  std::cout << "[compare] output directory: " << outDir << std::endl;
  std::cout << "[compare] done" << std::endl;
}
