#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
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
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
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
  bool fit_mis_excluded = false;
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

struct DriftCandidate {
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;
  double snrThreshold = 0.0;
  std::string type;
  int direction = 0;
  int validRuns = 0;
  double referenceGain = 0.0;
  double changePercent = 0.0;
  double slopePercentPerRun = 0.0;
  double spearman = 0.0;
  double residualMadPercent = 0.0;
  double lineResidualMadPercent = 0.0;
  double stepResidualMadPercent = 0.0;
  double lineInterceptPercent = 0.0;
  double stepBeforePercent = 0.0;
  double stepAfterPercent = 0.0;
  double stepX = -1.0;
  int stepRunIndex = -1;
  double score = 0.0;
};

struct ConsensusDriftCandidate {
  DriftCandidate representative;
  std::vector<DriftCandidate> matches;
  int consensusCount = 0;
  double score = 0.0;
};

struct StableChannelInfo {
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;
  int validRuns = 0;
  double referenceGain = 0.0;
  double deltaMadPercent = 0.0;
  double deltaRmsPercent = 0.0;
  double maxAbsDeltaPercent = 0.0;
  double slopePercentPerRun = 0.0;
  double spearman = 0.0;
  std::vector<double> deltasPercent;
};

struct GainRunRmsInfo {
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;
  int validRuns = 0;
  double meanGain = 0.0;
  double gainRms = 0.0;
  double minGain = 0.0;
  double maxGain = 0.0;
  std::vector<double> gainValues;
};

const double kDriftSnrThresholds[] = {3.0, 4.0, 5.0};
const int kDriftConsensusMin = 2;
const double kDriftMinAbsSpearman = 0.7;
const double kDriftMaxResidualMadPercent = 2.0;
const double kFitMisGainDeltaThreshold = 4.0;
const double kFitMisHighGainThreshold = 25.0;
const double kFitMisHighGainThresholdForSingle = 38.0;

std::vector<int> parseRunList(const char *csv) {
  std::vector<int> runs;
  if (!csv) return runs;

  std::stringstream ss(csv);
  std::string token;
  while (std::getline(ss, token, ',')) {
    token.erase(std::remove_if(token.begin(), token.end(),
                               [](unsigned char c) { return std::isspace(c); }),
                token.end());
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

std::string normalizeRunToken(std::string token) {
  token.erase(std::remove_if(token.begin(), token.end(),
                             [](unsigned char c) { return std::isspace(c); }),
              token.end());
  if (token.size() >= 2 &&
      ((token.front() == '"' && token.back() == '"') ||
       (token.front() == '\'' && token.back() == '\''))) {
    token = token.substr(1, token.size() - 2);
  }
  return token;
}

bool parseIntegerToken(const std::string &token, long &value) {
  if (token.empty()) return false;
  char *endptr = nullptr;
  value = std::strtol(token.c_str(), &endptr, 10);
  return endptr != token.c_str() && *endptr == '\0';
}

bool isNumericRunRangeToken(const std::string &token) {
  const auto dashPos = token.find('-');
  if (dashPos == std::string::npos || dashPos != token.rfind('-')) return false;
  const std::string a = token.substr(0, dashPos);
  const std::string b = token.substr(dashPos + 1);
  long start = 0;
  long end = 0;
  return parseIntegerToken(a, start) && parseIntegerToken(b, end);
}

bool parseRunTokenRange(const std::string &token, long &start, long &end) {
  if (parseIntegerToken(token, start)) {
    end = start;
    return true;
  }

  const auto dashPos = token.find('-');
  if (dashPos == std::string::npos || dashPos != token.rfind('-')) return false;
  const std::string a = token.substr(0, dashPos);
  const std::string b = token.substr(dashPos + 1);
  if (!parseIntegerToken(a, start) || !parseIntegerToken(b, end)) return false;
  if (start > end) std::swap(start, end);
  return true;
}

bool isAllRunToken(const std::string &token) {
  std::string lower = token;
  std::transform(lower.begin(), lower.end(), lower.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return lower == "all";
}

std::vector<std::string> parseTokenList(const char *csv) {
  std::vector<std::string> tokens;
  if (!csv) return tokens;
  std::stringstream ss(csv);
  std::string token;
  while (std::getline(ss, token, ',')) {
    token = normalizeRunToken(token);
    if (token.empty()) continue;
    tokens.push_back(token);
  }
  return tokens;
}

std::vector<std::string> listAllRunTokens(const char *baseDir) {
  std::vector<std::string> tokens;
  const std::string base = (baseDir && std::string(baseDir).size() > 0) ? baseDir : ".";
  void *dir = gSystem->OpenDirectory(base.c_str());
  if (!dir) {
    std::cerr << "[compare] warning: cannot open baseDir '" << base << "'\n";
    return tokens;
  }

  const char *entry = nullptr;
  while ((entry = gSystem->GetDirEntry(dir))) {
    const std::string token = entry;
    if (token == "." || token == "..") continue;
    const std::string path = base + "/" + token + "/spe_analysis.root";
    if (gSystem->AccessPathName(path.c_str()) == 0) tokens.push_back(token);
  }
  gSystem->FreeDirectory(dir);

  std::sort(tokens.begin(), tokens.end(),
            [](const std::string &a, const std::string &b) {
              long ai = 0;
              long bi = 0;
              const bool aIsInt = parseIntegerToken(a, ai);
              const bool bIsInt = parseIntegerToken(b, bi);
              if (aIsInt && bIsInt) return ai < bi;
              if (aIsInt != bIsInt) return aIsInt;
              return a < b;
            });
  return tokens;
}

std::vector<std::string> listRunTokensInRange(const char *baseDir,
                                              long start,
                                              long end) {
  std::vector<std::string> matches;
  const std::vector<std::string> allTokens = listAllRunTokens(baseDir);
  for (const auto &token : allTokens) {
    long tokenStart = 0;
    long tokenEnd = 0;
    if (!parseRunTokenRange(token, tokenStart, tokenEnd)) continue;
    if (tokenEnd < start || tokenStart > end) continue;
    matches.push_back(token);
  }
  return matches;
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
  if (onlyFitOk && e.fit_mis_excluded) return false;
  if (onlyFitOk && e.snr < 3.0) return false; // skip very low SNR cells

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

double computeMedian(std::vector<double> values) {
  if (values.empty()) return 0.0;

  const size_t mid = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + mid, values.end());
  const double upper = values[mid];
  if (values.size() % 2 != 0) return upper;

  std::nth_element(values.begin(), values.begin() + mid - 1, values.end());
  return 0.5 * (values[mid - 1] + upper);
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
    e.fit_mis_excluded = false;
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

int entryLayer(const SpeEntry &entry) {
  if (entry.layer >= 0) return entry.layer;
  if (entry.cellid >= 0) return entry.cellid / 100000;
  return -1;
}

int entryChip(const SpeEntry &entry) {
  if (entry.chip >= 0) return entry.chip;
  if (entry.cellid >= 0) return (entry.cellid / 10000) % 10;
  return -1;
}

int entryChannel(const SpeEntry &entry) {
  if (entry.channel >= 0) return entry.channel;
  if (entry.cellid >= 0) return entry.cellid % 10000;
  return -1;
}

bool isAcceptedSpePoint(const SpeEntry &entry) {
  return entry.fit_ok && !entry.fit_mis_excluded && entry.snr >= 5.0 &&
         entry.gain > 0.0;
}

bool isNonExcludedFitPoint(const SpeEntry &entry) {
  return entry.fit_ok && !entry.fit_mis_excluded && entry.gain > 0.0 && entry.snr >= 5.0;
}

bool isSurvivingSpePlotPoint(const SpeEntry &entry) {
  return entry.fit_ok && !entry.fit_mis_excluded && entry.snr > 5.0 &&
         entry.gain > 0.0 && std::isfinite(entry.gain);
}

std::set<int> selectCellsBySurvivalRate(const std::vector<RunData> &runs,
                                        double minRate = 0.50) {
  std::map<int, int> survivingCounts;
  for (const auto &run : runs) {
    std::set<int> survivedInRun;
    for (const auto &entry : run.entries) {
      if (entry.cellid < 0 || !isSurvivingSpePlotPoint(entry)) continue;
      survivedInRun.insert(entry.cellid);
    }
    for (int cellid : survivedInRun) {
      ++survivingCounts[cellid];
    }
  }

  std::set<int> selected;
  const double denom = static_cast<double>(runs.size());
  if (denom <= 0.0) return selected;
  for (const auto &kv : survivingCounts) {
    if (static_cast<double>(kv.second) / denom > minRate) selected.insert(kv.first);
  }
  return selected;
}

int markHighGainFittingMis(std::vector<RunData> &runs,
                           const std::string &reportFile,
                           double minGainDelta = kFitMisGainDeltaThreshold,
                           double highGainThreshold = kFitMisHighGainThreshold) {
  if (runs.size() < 2) return 0;
  minGainDelta = std::max(0.0, minGainDelta);
  highGainThreshold = std::max(0.0, highGainThreshold);

  struct GainObservation {
    size_t irun = 0;
    int cellid = -1;
    int layer = -1;
    int chip = -1;
    int channel = -1;
    double gain = 0.0;
    double snr = 0.0;
  };

  struct FitMisRecord {
    GainObservation high;
    GainObservation reference;
    double delta = 0.0;
  };

  std::map<int, std::vector<GainObservation> > observationsByCell;
  for (size_t irun = 0; irun < runs.size(); ++irun) {
    for (const auto &entry : runs[irun].entries) {
      if (!isAcceptedSpePoint(entry)) continue;
      GainObservation obs;
      obs.irun = irun;
      obs.cellid = entry.cellid;
      obs.layer = entryLayer(entry);
      obs.chip = entry.chip;
      obs.channel = entry.channel;
      obs.gain = entry.gain;
      obs.snr = entry.snr;
      observationsByCell[entry.cellid].push_back(obs);
    }
  }

  std::vector<std::set<int> > excludedByRun(runs.size());
  std::map<std::pair<size_t, int>, FitMisRecord> records;

  for (const auto &kv : observationsByCell) {
    const std::vector<GainObservation> &observations = kv.second;
    if (observations.size() < 2) continue;

    for (const auto &high : observations) {
      if (high.gain > kFitMisHighGainThresholdForSingle) {
        excludedByRun[high.irun].insert(high.cellid);
        continue;
      }
      if (high.gain < highGainThreshold) continue;

      bool foundLowerReference = false;
      GainObservation bestReference;
      double bestDelta = 0.0;
      for (const auto &reference : observations) {
        if (reference.irun == high.irun) continue;
        const double delta = high.gain - reference.gain;
        if (delta >= minGainDelta && delta > bestDelta) {
          foundLowerReference = true;
          bestReference = reference;
          bestDelta = delta;
        }
      }

      if (!foundLowerReference) continue;
      excludedByRun[high.irun].insert(high.cellid);
      FitMisRecord record;
      record.high = high;
      record.reference = bestReference;
      record.delta = bestDelta;
      records[std::make_pair(high.irun, high.cellid)] = record;
    }
  }

  for (size_t irun = 0; irun < runs.size(); ++irun) {
    if (excludedByRun[irun].empty()) continue;

    for (auto &entry : runs[irun].entries) {
      if (excludedByRun[irun].count(entry.cellid) == 0) continue;
      entry.fit_mis_excluded = true;
      entry.fit_ok = 0;
    }
    for (int cellid : excludedByRun[irun]) {
      auto it = runs[irun].byCellId.find(cellid);
      if (it == runs[irun].byCellId.end()) continue;
      it->second.fit_mis_excluded = true;
      it->second.fit_ok = 0;
    }
  }

  std::ofstream ofs(reportFile.c_str());
  if (!ofs) {
    std::cerr << "[compare] cannot write fitting-mis report " << reportFile
              << std::endl;
  } else {
    ofs << "# High-gain fitting-mis exclusions\n"
        << "# Excluded point: gain >= " << highGainThreshold
        << " and gain - another valid gain for the same cell >= " << minGainDelta
        << "\n"
        << "# Valid comparison point before exclusion: fit_ok=1, SNR>=3, gain>0\n\n";
    ofs << std::left << std::setw(20) << "run"
        << std::setw(10) << "cellid"
        << std::setw(7) << "layer"
        << std::setw(7) << "chip"
        << std::setw(9) << "channel"
        << std::right << std::setw(14) << "gain"
        << std::setw(14) << "SNR"
        << std::setw(20) << "reference_run"
        << std::setw(18) << "reference_gain"
        << std::setw(14) << "delta_gain" << "\n";
    ofs << std::string(133, '-') << "\n";
    ofs << std::fixed << std::setprecision(4);
    for (const auto &kv : records) {
      const FitMisRecord &record = kv.second;
      ofs << std::left << std::setw(20) << runs[record.high.irun].label
          << std::setw(10) << record.high.cellid
          << std::setw(7) << record.high.layer
          << std::setw(7) << record.high.chip
          << std::setw(9) << record.high.channel
          << std::right << std::setw(14) << record.high.gain
          << std::setw(14) << record.high.snr
          << std::setw(20) << runs[record.reference.irun].label
          << std::setw(18) << record.reference.gain
          << std::setw(14) << record.delta << "\n";
    }
  }

  int nExcluded = 0;
  for (const auto &cells : excludedByRun) nExcluded += static_cast<int>(cells.size());
  std::cout << "[compare] excluded " << nExcluded
            << " high-gain fitting-mis points (gain delta >= " << minGainDelta
            << ", high gain >= " << highGainThreshold << ")" << std::endl;
  return nExcluded;
}

void writeAndDrawLayerEfficiencyByRun(const std::vector<RunData> &runs,
                                      const std::string &outDir) {
  if (runs.empty()) return;

  struct LayerRunCounts {
    int total = 0;
    int excluded = 0;
    int accepted = 0;
  };

  std::set<int> layers;
  std::vector<std::map<int, LayerRunCounts> > counts(runs.size());
  for (size_t irun = 0; irun < runs.size(); ++irun) {
    for (const auto &entry : runs[irun].entries) {
      const int layer = entryLayer(entry);
      if (layer < 0) continue;
      layers.insert(layer);
      LayerRunCounts &count = counts[irun][layer];
      ++count.total;
      if (entry.fit_mis_excluded) ++count.excluded;
      if (isAcceptedSpePoint(entry)) ++count.accepted;
    }
  }
  if (layers.empty()) return;

  const std::string textFile = outDir + "/layer_efficiency_after_fitmis.txt";
  std::ofstream ofs(textFile.c_str());
  if (!ofs) {
    std::cerr << "[compare] cannot write layer efficiency file " << textFile
              << std::endl;
  } else {
    ofs << "# Layer efficiency by run after high-gain fitting-mis exclusion\n"
        << "# accepted: fit_ok=1, not fitting-mis-excluded, SNR>=3, gain>0\n"
        << "# denominator_after_exclusion = total cells in layer - fitting_mis_excluded\n\n";
    ofs << std::left << std::setw(20) << "run"
        << std::right << std::setw(8) << "layer"
        << std::setw(10) << "total"
        << std::setw(11) << "excluded"
        << std::setw(14) << "denom_after"
        << std::setw(11) << "accepted"
        << std::setw(16) << "eff_after[%]"
        << std::setw(14) << "eff_all[%]" << "\n";
    ofs << std::string(104, '-') << "\n";
    ofs << std::fixed << std::setprecision(4);
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      for (int layer : layers) {
        const auto it = counts[irun].find(layer);
        LayerRunCounts count;
        if (it != counts[irun].end()) count = it->second;
        const int denomAfter = std::max(0, count.total);
        const double effAfter =
            denomAfter > 0 ? 100.0 * count.accepted / denomAfter : 0.0;
        const double effAll =
            count.total > 0 ? 100.0 * count.accepted / count.total : 0.0;
        ofs << std::left << std::setw(20) << runs[irun].label
            << std::right << std::setw(8) << layer
            << std::setw(10) << count.total
            << std::setw(11) << count.excluded
            << std::setw(14) << denomAfter
            << std::setw(11) << count.accepted
            << std::setw(16) << effAfter
            << std::setw(14) << effAll << "\n";
      }
    }
  }

  const int colors[] = {kRed + 1,     kBlue + 1,    kGreen + 2,  kMagenta + 1,
                        kOrange + 7,  kCyan + 2,    kViolet + 1, kSpring + 5,
                        kAzure + 1,   kPink + 1,    kTeal + 3,   kGray + 2};
  const int markers[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33};

  TCanvas *c = new TCanvas("c_layer_efficiency_after_fitmis",
                           "Layer efficiency after fitting-mis exclusion", 1200, 800);
  c->SetGrid();
  TMultiGraph *mg = new TMultiGraph();
  TLegend *leg = new TLegend(0.72, 0.12, 0.90, 0.88);
  if (layers.size() > 24) {
    leg->SetNColumns(3);
  } else if (layers.size() > 12) {
    leg->SetNColumns(2);
  }

  int ilayer = 0;
  for (int layer : layers) {
    std::vector<double> x;
    std::vector<double> y;
    x.reserve(runs.size());
    y.reserve(runs.size());
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      const auto it = counts[irun].find(layer);
      LayerRunCounts count;
      if (it != counts[irun].end()) count = it->second;
      const int denomAfter = std::max(0, count.total);
      if (denomAfter <= 0) continue;
      x.push_back(runs[irun].run);
      y.push_back(100.0 * count.accepted / denomAfter);
    }
    if (x.empty()) continue;
    TGraph *g = new TGraph(static_cast<int>(x.size()), &x[0], &y[0]);
    setGraphStyle(g, colors[ilayer % 12], markers[ilayer % 12]);
    mg->Add(g, "LP");
    leg->AddEntry(g, Form("layer %d", layer), "lp");
    ++ilayer;
  }

  if (ilayer == 0) {
    delete leg;
    delete mg;
    delete c;
    return;
  }

  mg->SetTitle("Layer efficiency after fitting-mis exclusion;Run;Efficiency [%]");
  mg->Draw("A");
  mg->GetYaxis()->SetRangeUser(0.0, 105.0);
  leg->Draw();
  c->SaveAs((outDir + "/layer_efficiency_after_fitmis.pdf").c_str());
  c->SaveAs((outDir + "/layer_efficiency_after_fitmis.png").c_str());

  TCanvas *cByRun = new TCanvas("c_layer_efficiency_by_run_after_fitmis",
                                "Layer efficiency by run after fitting-mis exclusion",
                                1600, 800);
  cByRun->SetGrid();
  cByRun->SetRightMargin(0.30);
  TMultiGraph *mgByRun = new TMultiGraph();
  TLegend *legByRun = new TLegend(0.70, 0.12, 1.0, 0.88);
  if (runs.size() > 24) {
    legByRun->SetNColumns(3);
  } else if (runs.size() > 12) {
    legByRun->SetNColumns(2);
  }

  int drawnRuns = 0;
  for (size_t irun = 0; irun < runs.size(); ++irun) {
    std::vector<double> x;
    std::vector<double> y;
    x.reserve(layers.size());
    y.reserve(layers.size());
    double totalEff = 0.0;
    double layerover5Eff = 0.0;
    int totalDenom = 0;
    int layerover5Denom = 0;
    for (int layer : layers) {
      const auto it = counts[irun].find(layer);
      LayerRunCounts count;
      if (it != counts[irun].end()) count = it->second;
      const int denomAfter = std::max(0, count.total);
      totalDenom += denomAfter;
      totalEff += static_cast<double>(count.accepted);
      if (layer >= 5) {
        layerover5Denom += denomAfter;
        layerover5Eff += static_cast<double>(count.accepted);
      }
      if (denomAfter <= 0) continue;
      x.push_back(layer);
      y.push_back(100.0 * count.accepted / denomAfter);
    }
    if (x.empty()) continue;
    TGraph *g = new TGraph(static_cast<int>(x.size()), &x[0], &y[0]);
    setGraphStyle(g, colors[drawnRuns % 12], markers[drawnRuns % 12]);
    mgByRun->Add(g, "LP");
    legByRun->AddEntry(
      g,
      (runs[irun].label + " eff: " +
       std::to_string(totalDenom > 0 ? totalEff / totalDenom * 100.0 : 0.0).substr(0, 4) + "%" +
        " (L5+: " + std::to_string(layerover5Denom > 0 ? layerover5Eff / layerover5Denom * 100.0 : 0.0).substr(0, 4) + "%)")
        .c_str(),
      "lp");
    ++drawnRuns;
  }

  if (drawnRuns > 0) {
    mgByRun->SetTitle(
        "Layer efficiency by run after fitting-mis exclusion;Layer;Efficiency [%]");
    mgByRun->Draw("A");
    mgByRun->GetYaxis()->SetRangeUser(0.0, 105.0);
    legByRun->Draw();
    cByRun->SaveAs((outDir + "/layer_efficiency_by_run_after_fitmis.pdf").c_str());
    cByRun->SaveAs((outDir + "/layer_efficiency_by_run_after_fitmis.png").c_str());
  }

  std::cout << "[compare] wrote layer efficiency after fitting-mis exclusion to "
            << textFile << std::endl;
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

void writeShiftedGainReport(const std::vector<RunData> &runs,
                            const std::string &fileName,
                            double relativeThreshold = 0.10,
                            int minValidRuns = 3) {
  if (runs.empty()) return;
  relativeThreshold = std::max(0.0, relativeThreshold);
  minValidRuns = std::max(1, minValidRuns);

  struct ChannelReference {
    double medianGain = 0.0;
    double gainMad = 0.0;
    int validRuns = 0;
  };

  std::set<int> allCellIds;
  for (const auto &run : runs) {
    for (const auto &kv : run.byCellId) allCellIds.insert(kv.first);
  }

  std::map<int, ChannelReference> references;
  std::set<int> shiftedCellIds;
  std::vector<std::set<int> > shiftedByRun(runs.size());

  for (int cellid : allCellIds) {
    std::vector<double> gains;
    for (const auto &run : runs) {
      const auto it = run.byCellId.find(cellid);
      if (it == run.byCellId.end()) continue;
      double gain = 0.0;
      if (getBranchValue(it->second, "gain", true, gain) && gain > 0.0) {
        gains.push_back(gain);
      }
    }
    if (static_cast<int>(gains.size()) < minValidRuns) continue;

    ChannelReference ref;
    ref.medianGain = computeMedian(gains);
    ref.validRuns = static_cast<int>(gains.size());
    std::vector<double> deviations;
    deviations.reserve(gains.size());
    for (double gain : gains) deviations.push_back(std::abs(gain - ref.medianGain));
    ref.gainMad = computeMedian(deviations);
    references[cellid] = ref;

    if (ref.medianGain == 0.0) continue;
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      const auto it = runs[irun].byCellId.find(cellid);
      if (it == runs[irun].byCellId.end()) continue;
      double gain = 0.0;
      if (!getBranchValue(it->second, "gain", true, gain) || gain <= 0.0) continue;
      const double relativeShift = std::abs(gain - ref.medianGain) / std::abs(ref.medianGain);
      if (relativeShift >= relativeThreshold) {
        shiftedCellIds.insert(cellid);
        shiftedByRun[irun].insert(cellid);
      }
    }
  }

  std::ofstream ofs(fileName.c_str());
  if (!ofs) {
    std::cerr << "[compare] cannot write shifted gain report " << fileName << std::endl;
    return;
  }

  ofs << "# Shifted SPE gain channel report\n"
      << "# Valid measurement: fit_ok=1, not fitting-mis-excluded, SNR>=3, and gain>0\n"
      << "# Per-channel reference: median gain over valid runs\n"
      << "# SHIFTED: abs(gain - reference_gain) / abs(reference_gain) >= "
      << 100.0 * relativeThreshold << "%\n"
      << "# A channel is listed when at least one run is SHIFTED and it has at least "
      << minValidRuns << " valid runs.\n\n";

  ofs << "=== RUN SUMMARY ===\n";
  ofs << std::left << std::setw(20) << "run"
      << std::right << std::setw(9) << "nCell"
      << std::setw(9) << "nValid"
      << std::setw(14) << "gain_mean"
      << std::setw(14) << "gain_rms"
      << std::setw(14) << "SNR_mean"
      << std::setw(14) << "SNR_rms"
      << std::setw(16) << "entries_mean"
      << std::setw(14) << "nShifted" << "\n";
  ofs << std::string(124, '-') << "\n";

  ofs << std::fixed << std::setprecision(4);
  for (size_t irun = 0; irun < runs.size(); ++irun) {
    std::vector<double> gainValues;
    std::vector<double> snrValues;
    std::vector<double> entryValues;
    gainValues.reserve(runs[irun].entries.size());
    snrValues.reserve(runs[irun].entries.size());
    entryValues.reserve(runs[irun].entries.size());
    for (const auto &entry : runs[irun].entries) {
      double validGain = 0.0;
      if (getBranchValue(entry, "gain", true, validGain) && validGain > 0.0) {
        gainValues.push_back(validGain);
        snrValues.push_back(entry.snr);
        entryValues.push_back(entry.entries);
      }
    }
    const SummaryStats gain = computeStats(gainValues);
    const SummaryStats snr = computeStats(snrValues);
    const SummaryStats entries = computeStats(entryValues);

    ofs << std::left << std::setw(20) << runs[irun].label
        << std::right << std::setw(9) << runs[irun].entries.size()
        << std::setw(9) << gain.n
        << std::setw(14) << gain.mean
        << std::setw(14) << gain.rms
        << std::setw(14) << snr.mean
        << std::setw(14) << snr.rms
        << std::setw(16) << entries.mean
        << std::setw(14) << shiftedByRun[irun].size() << "\n";
  }

  ofs << "\n=== SHIFTED CHANNELS: GAIN AND SNR BY RUN ===\n";
  ofs << "# Two rows per channel: gain followed by SNR. Only channels shifted in at least "
         "one run are listed.\n";
  ofs << "# Values with fit_ok=0, fitting-mis exclusion, SNR<3, or gain<=0 are "
         "written as NA and are not used for the reference or shift decision.\n";
  ofs << "# A '*' after gain marks a run satisfying the gain-shift criterion.\n";
  ofs << std::left << std::setw(10) << "cellid"
      << std::setw(7) << "layer"
      << std::setw(7) << "chip"
      << std::setw(9) << "channel"
      << std::setw(8) << "metric"
      << std::right << std::setw(12) << "ref_gain"
      << std::setw(12) << "gain_MAD"
      << std::setw(11) << "ref_nRun"
      << std::setw(13) << "max_abs_d[%]"
      << std::setw(20) << "max_shift_run";
  for (const auto &run : runs) {
    ofs << std::setw(18) << run.label;
  }
  ofs << "\n";

  for (int cellid : shiftedCellIds) {
    const ChannelReference &ref = references[cellid];
    const SpeEntry *identity = nullptr;
    double maxAbsDeltaPercent = 0.0;
    std::string maxShiftRun = "NA";

    for (const auto &run : runs) {
      const auto it = run.byCellId.find(cellid);
      if (it == run.byCellId.end()) continue;
      if (!identity) identity = &it->second;
      const SpeEntry &entry = it->second;
      const bool valid = isAcceptedSpePoint(entry);
      if (!valid || ref.medianGain == 0.0) continue;
      const double deltaPercent =
          100.0 * (entry.gain - ref.medianGain) / std::abs(ref.medianGain);
      if (std::abs(deltaPercent) > maxAbsDeltaPercent) {
        maxAbsDeltaPercent = std::abs(deltaPercent);
        maxShiftRun = run.label;
      }
    }

    if (!identity) continue;
    std::vector<std::string> gainTexts;
    std::vector<std::string> snrTexts;
    gainTexts.reserve(runs.size());
    snrTexts.reserve(runs.size());
    for (const auto &run : runs) {
      const auto it = run.byCellId.find(cellid);
      if (it == run.byCellId.end()) {
        gainTexts.push_back("NA");
        snrTexts.push_back("NA");
        continue;
      }
      const SpeEntry &entry = it->second;
      const bool valid = isAcceptedSpePoint(entry);
      if (!valid) {
        gainTexts.push_back("NA");
        snrTexts.push_back("NA");
        continue;
      }

      const double deltaPercent =
          ref.medianGain != 0.0
              ? 100.0 * (entry.gain - ref.medianGain) / std::abs(ref.medianGain)
              : 0.0;
      const bool shifted = std::abs(deltaPercent) >= 100.0 * relativeThreshold;
      std::ostringstream gainText;
      gainText << std::fixed << std::setprecision(4) << entry.gain;
      if (shifted) gainText << "*";
      std::ostringstream snrText;
      snrText << std::fixed << std::setprecision(4) << entry.snr;
      gainTexts.push_back(gainText.str());
      snrTexts.push_back(snrText.str());
    }

    ofs << std::left << std::setw(10) << cellid
        << std::setw(7) << identity->layer
        << std::setw(7) << identity->chip
        << std::setw(9) << identity->channel
        << std::setw(8) << "gain"
        << std::right << std::setw(12) << ref.medianGain
        << std::setw(12) << ref.gainMad
        << std::setw(11) << ref.validRuns
        << std::setw(13) << maxAbsDeltaPercent
        << std::setw(20) << maxShiftRun;
    for (const auto &value : gainTexts) ofs << std::setw(18) << value;
    ofs << "\n";

    ofs << std::left << std::setw(10) << ""
        << std::setw(7) << ""
        << std::setw(7) << ""
        << std::setw(9) << ""
        << std::setw(8) << "SNR"
        << std::right << std::setw(12) << ""
        << std::setw(12) << ""
        << std::setw(11) << ""
        << std::setw(13) << ""
        << std::setw(20) << "";
    for (const auto &value : snrTexts) ofs << std::setw(18) << value;
    ofs << "\n\n";
  }

  ofs << "# shifted channels: " << shiftedCellIds.size() << " / "
      << references.size() << " channels with sufficient valid runs\n";
  std::cout << "[compare] wrote " << shiftedCellIds.size()
            << " shifted gain channels to " << fileName << std::endl;
}

struct GainDeltaStudy {
  double snrThreshold = 0.0;
  std::vector<double> allDeltasPercent;
  std::vector<std::vector<double> > stableDeltasByRun;
  int referenceChannels = 0;
  int stableChannels = 0;
};

GainDeltaStudy buildGainDeltaStudy(const std::vector<RunData> &runs,
                                   double snrThreshold,
                                   double stableLimitPercent = 10.0,
                                   int minValidRuns = 3) {
  GainDeltaStudy study;
  study.snrThreshold = snrThreshold;
  study.stableDeltasByRun.resize(runs.size());

  std::set<int> allCellIds;
  for (const auto &run : runs) {
    for (const auto &kv : run.byCellId) allCellIds.insert(kv.first);
  }

  for (int cellid : allCellIds) {
    std::vector<double> gains;
    for (const auto &run : runs) {
      const auto it = run.byCellId.find(cellid);
      if (it == run.byCellId.end()) continue;
      const SpeEntry &entry = it->second;
      if (isAcceptedSpePoint(entry) && entry.snr > snrThreshold) {
        gains.push_back(entry.gain);
      }
    }
    if (static_cast<int>(gains.size()) < minValidRuns) continue;

    const double referenceGain = computeMedian(gains);
    if (referenceGain == 0.0) continue;
    ++study.referenceChannels;

    bool stable = true;
    std::vector<std::pair<size_t, double> > deltasByRun;
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      const auto it = runs[irun].byCellId.find(cellid);
      if (it == runs[irun].byCellId.end()) continue;
      const SpeEntry &entry = it->second;
      if (!isAcceptedSpePoint(entry) || entry.snr <= snrThreshold) continue;

      const double deltaPercent =
          100.0 * (entry.gain - referenceGain) / std::abs(referenceGain);
      study.allDeltasPercent.push_back(deltaPercent);
      deltasByRun.push_back(std::make_pair(irun, deltaPercent));
      if (std::abs(deltaPercent) > stableLimitPercent) stable = false;
    }

    if (!stable) continue;
    ++study.stableChannels;
    for (const auto &item : deltasByRun) {
      study.stableDeltasByRun[item.first].push_back(item.second);
    }
  }

  return study;
}

void drawGainDeltaSnrStudies(const std::vector<RunData> &runs,
                             const std::string &outDir,
                             double stableLimitPercent = 10.0,
                             int minValidRuns = 3) {
  if (runs.empty()) return;

  std::vector<GainDeltaStudy> studies;
  for (double threshold : kDriftSnrThresholds) {
    studies.push_back(
        buildGainDeltaStudy(runs, threshold, stableLimitPercent, minValidRuns));
  }

  const int colors[] = {kRed + 1, kBlue + 1, kGreen + 2};
  const int markers[] = {20, 21, 22};

  {
    TCanvas *c = new TCanvas("c_gain_delta_snr_thresholds",
                             "Gain delta by SNR threshold", 1000, 700);
    c->SetGrid();
    c->SetLogy();
    TLegend *leg = new TLegend(0.60, 0.67, 0.89, 0.89);
    bool first = true;
    double ymax = 0.0;
    std::vector<TH1D *> histograms;

    for (size_t i = 0; i < studies.size(); ++i) {
      TH1D *h = new TH1D(
          Form("h_gain_delta_snr_gt%d", static_cast<int>(studies[i].snrThreshold)),
          "Gain relative delta by SNR threshold;#Delta gain / median gain [%];"
          "Normalized entries",
          200, -100.0, 100.0);
      for (double delta : studies[i].allDeltasPercent) h->Fill(delta);
      if (h->Integral() > 0.0) h->Scale(1.0 / h->Integral());
      h->SetLineColor(colors[i]);
      h->SetLineWidth(2);
      h->SetStats(0);
      ymax = std::max(ymax, h->GetMaximum() * 1.3);

      if (first) {
        h->Draw("HIST");
        first = false;
      } else {
        h->Draw("HIST SAME");
      }
      histograms.push_back(h);
      leg->AddEntry(
          h,
          Form("SNR > %.0f (N=%zu)", studies[i].snrThreshold,
               studies[i].allDeltasPercent.size()),
          "l");
    }
    for (TH1D *h : histograms) {
      h->SetMaximum(ymax);
      h->SetMinimum(1e-7);
    }
    leg->Draw();
    c->SaveAs((outDir + "/gain_delta_snr_thresholds_overlay.pdf").c_str());
    c->SaveAs((outDir + "/gain_delta_snr_thresholds_overlay.png").c_str());
  }

  int runColors[] = {kRed + 1,     kBlue + 1,    kGreen + 2,  kMagenta + 1,
                     kOrange + 7,  kCyan + 2,    kRed + 3,    kBlue + 3,
                     kGreen + 4,   kMagenta + 3, kOrange + 9, kCyan + 4,
                     kGray + 2};
  for (const auto &study : studies) {
    const int threshold = static_cast<int>(study.snrThreshold);
    TCanvas *c = new TCanvas(Form("c_stable_gain_delta_snr_gt%d", threshold),
                             "Stable gain delta by run", 1200, 800);
    c->SetGrid();
    TLegend *leg = new TLegend(0.68, 0.50, 0.90, 0.89);
    bool first = true;
    double ymax = 0.0;
    std::vector<TH1D *> histograms;

    for (size_t irun = 0; irun < runs.size(); ++irun) {
      if (study.stableDeltasByRun[irun].empty()) continue;
      TH1D *h = new TH1D(
          Form("h_stable_delta_snr_gt%d_run_%s", threshold,
               sanitizeLabel(runs[irun].label).c_str()),
          Form("Stable channels (SNR > %d, |#Delta| #leq %.0f%%);"
               "#Delta gain / median gain [%%];Normalized entries",
               threshold, stableLimitPercent),
          100, -stableLimitPercent, stableLimitPercent);
      for (double delta : study.stableDeltasByRun[irun]) h->Fill(delta);
      if (h->Integral() > 0.0) h->Scale(1.0 / h->Integral());
      h->SetLineColor(runColors[irun % 13]);
      h->SetLineWidth(2);
      h->SetStats(0);
      ymax = std::max(ymax, h->GetMaximum() * 1.25);

      if (first) {
        h->Draw("HIST");
        first = false;
      } else {
        h->Draw("HIST SAME");
      }
      histograms.push_back(h);
      leg->AddEntry(
          h,
          Form("%s (mean=%.3f%%)", runs[irun].label.c_str(),
               computeStats(study.stableDeltasByRun[irun]).mean),
          "l");
    }
    for (TH1D *h : histograms) h->SetMaximum(ymax);
    leg->Draw();
    c->SaveAs(
        (outDir + Form("/stable_gain_delta_by_run_snr_gt%d.pdf", threshold)).c_str());
    c->SaveAs(
        (outDir + Form("/stable_gain_delta_by_run_snr_gt%d.png", threshold)).c_str());
  }

  {
    TCanvas *c =
        new TCanvas("c_stable_gain_delta_mean_trend", "Stable delta mean trend",
                    1400, 800);
    c->SetGrid();
    c->SetBottomMargin(0.24);
    c->SetRightMargin(0.04);
    TLegend *leg = new TLegend(0.70, 0.70, 0.90, 0.89);
    double maxAbsMean = 0.0;
    for (const auto &study : studies) {
      for (const auto &values : study.stableDeltasByRun) {
        if (values.empty()) continue;
        maxAbsMean = std::max(maxAbsMean, std::abs(computeStats(values).mean));
      }
    }
    const double meanRange = std::max(0.5, 1.25 * maxAbsMean);
    TH1D *frame = new TH1D(
        "h_stable_gain_delta_mean_frame",
        Form("Mean gain delta of channels always within %.0f%%;;"
             "Mean #Delta gain / median gain [%%]",
             stableLimitPercent),
        static_cast<int>(runs.size()), 0.5, runs.size() + 0.5);
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      frame->GetXaxis()->SetBinLabel(static_cast<int>(irun) + 1,
                                    runs[irun].label.c_str());
    }
    frame->GetXaxis()->LabelsOption("v");
    frame->SetMinimum(-meanRange);
    frame->SetMaximum(meanRange);
    frame->SetStats(0);
    frame->Draw();

    for (size_t istudy = 0; istudy < studies.size(); ++istudy) {
      std::vector<double> x;
      std::vector<double> y;
      for (size_t irun = 0; irun < runs.size(); ++irun) {
        if (studies[istudy].stableDeltasByRun[irun].empty()) continue;
        x.push_back(irun + 1.0);
        y.push_back(computeStats(studies[istudy].stableDeltasByRun[irun]).mean);
      }
      if (x.empty()) continue;
      TGraph *g = new TGraph(static_cast<int>(x.size()), &x[0], &y[0]);
      setGraphStyle(g, colors[istudy], markers[istudy]);
      g->Draw("LP SAME");
      leg->AddEntry(
          g,
          Form("SNR > %.0f (%d stable channels)", studies[istudy].snrThreshold,
               studies[istudy].stableChannels),
          "lp");
    }
    TLine *zero = new TLine(0.5, 0.0, runs.size() + 0.5, 0.0);
    zero->SetLineStyle(2);
    zero->Draw();
    leg->Draw();
    c->SaveAs((outDir + "/stable_gain_delta_mean_trend.pdf").c_str());
    c->SaveAs((outDir + "/stable_gain_delta_mean_trend.png").c_str());
  }

  std::ofstream ofs((outDir + "/stable_gain_delta_mean_trend.txt").c_str());
  if (ofs) {
    ofs << "# delta[%] = 100 * (gain - per-channel median gain) / median gain\n"
        << "# Stable channels never exceed |delta| > " << stableLimitPercent
        << "% among valid measurements for each SNR threshold.\n"
        << "# snr_threshold run stable_channels entries mean_delta_percent "
           "rms_delta_percent\n";
    ofs << std::fixed << std::setprecision(6);
    for (const auto &study : studies) {
      for (size_t irun = 0; irun < runs.size(); ++irun) {
        const SummaryStats stats = computeStats(study.stableDeltasByRun[irun]);
        ofs << study.snrThreshold << " " << runs[irun].label << " "
            << study.stableChannels << " " << stats.n << " " << stats.mean << " "
            << stats.rms << "\n";
      }
    }
  }
}

double medianAbsoluteDeviation(const std::vector<double> &values, double center) {
  std::vector<double> deviations;
  deviations.reserve(values.size());
  for (double value : values) deviations.push_back(std::abs(value - center));
  return computeMedian(deviations);
}

std::vector<double> computeRanks(const std::vector<double> &values) {
  std::vector<size_t> order(values.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&values](size_t a, size_t b) {
    if (values[a] != values[b]) return values[a] < values[b];
    return a < b;
  });

  std::vector<double> ranks(values.size(), 0.0);
  size_t begin = 0;
  while (begin < order.size()) {
    size_t end = begin + 1;
    while (end < order.size() && values[order[end]] == values[order[begin]]) ++end;
    const double rank = 0.5 * (begin + end - 1) + 1.0;
    for (size_t i = begin; i < end; ++i) ranks[order[i]] = rank;
    begin = end;
  }
  return ranks;
}

double computePearsonCorrelation(const std::vector<double> &a,
                                 const std::vector<double> &b) {
  if (a.size() != b.size() || a.size() < 2) return 0.0;
  const double meanA = std::accumulate(a.begin(), a.end(), 0.0) / a.size();
  const double meanB = std::accumulate(b.begin(), b.end(), 0.0) / b.size();
  double covariance = 0.0;
  double varianceA = 0.0;
  double varianceB = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    const double da = a[i] - meanA;
    const double db = b[i] - meanB;
    covariance += da * db;
    varianceA += da * da;
    varianceB += db * db;
  }
  if (varianceA <= 0.0 || varianceB <= 0.0) return 0.0;
  return covariance / std::sqrt(varianceA * varianceB);
}

double computeSpearmanCorrelation(const std::vector<double> &x,
                                  const std::vector<double> &y) {
  return computePearsonCorrelation(computeRanks(x), computeRanks(y));
}

bool fitTheilSen(const std::vector<double> &x,
                 const std::vector<double> &y,
                 double &slope,
                 double &intercept,
                 double &residualMad) {
  if (x.size() != y.size() || x.size() < 2) return false;
  std::vector<double> slopes;
  for (size_t i = 0; i < x.size(); ++i) {
    for (size_t j = i + 1; j < x.size(); ++j) {
      if (x[j] != x[i]) slopes.push_back((y[j] - y[i]) / (x[j] - x[i]));
    }
  }
  if (slopes.empty()) return false;
  slope = computeMedian(slopes);

  std::vector<double> intercepts;
  intercepts.reserve(x.size());
  for (size_t i = 0; i < x.size(); ++i) intercepts.push_back(y[i] - slope * x[i]);
  intercept = computeMedian(intercepts);

  std::vector<double> residuals;
  residuals.reserve(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    residuals.push_back(y[i] - (intercept + slope * x[i]));
  }
  residualMad = medianAbsoluteDeviation(residuals, 0.0);
  return true;
}

bool findPersistentStep(const std::vector<double> &x,
                        const std::vector<double> &y,
                        double minChangePercent,
                        double maxSegmentMadPercent,
                        double &beforeMedian,
                        double &afterMedian,
                        double &stepX,
                        int &stepRunIndex,
                        double &residualMad) {
  if (x.size() != y.size() || x.size() < 6) return false;

  bool found = false;
  double bestResidual = 1e30;
  for (size_t split = 3; split + 3 <= y.size(); ++split) {
    std::vector<double> before(y.begin(), y.begin() + split);
    std::vector<double> after(y.begin() + split, y.end());
    const double beforeMed = computeMedian(before);
    const double afterMed = computeMedian(after);
    const double beforeMad = medianAbsoluteDeviation(before, beforeMed);
    const double afterMad = medianAbsoluteDeviation(after, afterMed);
    const double step = afterMed - beforeMed;
    const double segmentNoise = std::max(beforeMad, afterMad);
    if (std::abs(step) < minChangePercent) continue;
    if (beforeMad > maxSegmentMadPercent || afterMad > maxSegmentMadPercent) continue;
    if (std::abs(step) < 3.0 * segmentNoise) continue;

    int beforeConsistent = 0;
    int afterConsistent = 0;
    for (double value : before) {
      if (std::abs(value - beforeMed) < std::abs(value - afterMed)) ++beforeConsistent;
    }
    for (double value : after) {
      if (std::abs(value - afterMed) < std::abs(value - beforeMed)) ++afterConsistent;
    }
    if (3 * beforeConsistent < 2 * static_cast<int>(before.size()) ||
        3 * afterConsistent < 2 * static_cast<int>(after.size())) {
      continue;
    }

    std::vector<double> residuals;
    residuals.reserve(y.size());
    for (size_t i = 0; i < split; ++i) residuals.push_back(y[i] - beforeMed);
    for (size_t i = split; i < y.size(); ++i) residuals.push_back(y[i] - afterMed);
    const double candidateResidual = medianAbsoluteDeviation(residuals, 0.0);
    if (candidateResidual >= bestResidual) continue;

    found = true;
    bestResidual = candidateResidual;
    beforeMedian = beforeMed;
    afterMedian = afterMed;
    stepX = 0.5 * (x[split - 1] + x[split]);
    stepRunIndex = static_cast<int>(std::lround(x[split]));
    residualMad = candidateResidual;
  }
  return found;
}

std::map<int, DriftCandidate> findDriftCandidates(
    const std::vector<RunData> &runs,
    double snrThreshold,
    double minChangePercent = 3.0,
    int minValidRuns = 8,
    double minAbsSpearman = kDriftMinAbsSpearman,
    double maxResidualMadPercent = kDriftMaxResidualMadPercent) {
  std::map<int, DriftCandidate> candidates;
  std::set<int> allCellIds;
  for (const auto &run : runs) {
    for (const auto &kv : run.byCellId) allCellIds.insert(kv.first);
  }

  for (int cellid : allCellIds) {
    std::vector<double> x;
    std::vector<double> gains;
    const SpeEntry *identity = nullptr;
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      const auto it = runs[irun].byCellId.find(cellid);
      if (it == runs[irun].byCellId.end()) continue;
      if (!identity) identity = &it->second;
      const SpeEntry &entry = it->second;
      if (!isAcceptedSpePoint(entry) || entry.snr <= snrThreshold) continue;
      x.push_back(static_cast<double>(irun));
      gains.push_back(entry.gain);
    }
    if (!identity || static_cast<int>(gains.size()) < minValidRuns) continue;

    const double referenceGain = computeMedian(gains);
    if (referenceGain <= 0.0) continue;
    std::vector<double> deltaPercent;
    deltaPercent.reserve(gains.size());
    for (double gain : gains) {
      deltaPercent.push_back(100.0 * (gain - referenceGain) / referenceGain);
    }

    double slope = 0.0;
    double intercept = 0.0;
    double lineResidual = 0.0;
    const bool haveLine = fitTheilSen(x, deltaPercent, slope, intercept, lineResidual);
    const double spearman = computeSpearmanCorrelation(x, deltaPercent);
    const double predictedChange =
        haveLine ? slope * (x.back() - x.front()) : 0.0;
    const bool driftPass =
        haveLine && std::abs(predictedChange) >= minChangePercent &&
        std::abs(spearman) >= minAbsSpearman &&
        lineResidual <= maxResidualMadPercent;

    double beforeMedian = 0.0;
    double afterMedian = 0.0;
    double stepX = -1.0;
    int stepRunIndex = -1;
    double stepResidual = 0.0;
    const bool stepPass =
        findPersistentStep(x, deltaPercent, minChangePercent, maxResidualMadPercent,
                           beforeMedian, afterMedian, stepX, stepRunIndex, stepResidual);

    if (!driftPass && !stepPass) continue;

    DriftCandidate candidate;
    candidate.cellid = cellid;
    candidate.layer = identity->layer;
    candidate.chip = identity->chip;
    candidate.channel = identity->channel;
    candidate.snrThreshold = snrThreshold;
    candidate.validRuns = static_cast<int>(gains.size());
    candidate.referenceGain = referenceGain;
    candidate.slopePercentPerRun = slope;
    candidate.spearman = spearman;
    candidate.lineResidualMadPercent = lineResidual;
    candidate.stepResidualMadPercent = stepPass ? stepResidual : -1.0;
    candidate.lineInterceptPercent = intercept;
    candidate.stepBeforePercent = beforeMedian;
    candidate.stepAfterPercent = afterMedian;
    candidate.stepX = stepX;
    candidate.stepRunIndex = stepRunIndex;

    const bool chooseStep =
        stepPass && (!driftPass || stepResidual <= 0.70 * lineResidual);
    if (chooseStep) {
      candidate.type = "step";
      candidate.changePercent = afterMedian - beforeMedian;
      candidate.residualMadPercent = stepResidual;
    } else {
      candidate.type = "drift";
      candidate.changePercent = predictedChange;
      candidate.residualMadPercent = lineResidual;
    }
    candidate.direction = candidate.changePercent >= 0.0 ? 1 : -1;
    candidate.score =
        std::abs(candidate.changePercent) / std::max(0.10, candidate.residualMadPercent);
    candidates[cellid] = candidate;
  }
  return candidates;
}

std::vector<ConsensusDriftCandidate> buildDriftConsensus(
    const std::vector<std::map<int, DriftCandidate> > &candidateSets,
    int minConsensus = 2) {
  std::set<int> allCellIds;
  for (const auto &candidateSet : candidateSets) {
    for (const auto &kv : candidateSet) allCellIds.insert(kv.first);
  }

  std::vector<ConsensusDriftCandidate> consensus;
  for (int cellid : allCellIds) {
    std::vector<DriftCandidate> available;
    for (const auto &candidateSet : candidateSets) {
      const auto it = candidateSet.find(cellid);
      if (it == candidateSet.end()) continue;
      available.push_back(it->second);
    }

    std::vector<DriftCandidate> matches;
    double bestGroupScore = -1.0;
    for (const auto &anchor : available) {
      std::vector<DriftCandidate> group;
      double groupScore = 0.0;
      for (const auto &candidate : available) {
        if (candidate.direction != anchor.direction ||
            candidate.type != anchor.type) {
          continue;
        }
        if (candidate.type == "step" &&
            std::abs(candidate.stepRunIndex - anchor.stepRunIndex) > 1) {
          continue;
        }
        group.push_back(candidate);
        groupScore += candidate.score;
      }
      if (group.size() > matches.size() ||
          (group.size() == matches.size() && groupScore > bestGroupScore)) {
        matches = group;
        bestGroupScore = groupScore;
      }
    }
    if (static_cast<int>(matches.size()) < minConsensus) continue;
    std::sort(matches.begin(), matches.end(),
              [](const DriftCandidate &a, const DriftCandidate &b) {
                if (a.score != b.score) return a.score > b.score;
                return a.snrThreshold > b.snrThreshold;
              });

    ConsensusDriftCandidate item;
    item.representative = matches.front();
    item.matches = matches;
    item.consensusCount = static_cast<int>(matches.size());
    item.score = matches.front().score * item.consensusCount;
    consensus.push_back(item);
  }

  std::sort(consensus.begin(), consensus.end(),
            [](const ConsensusDriftCandidate &a, const ConsensusDriftCandidate &b) {
              if (a.consensusCount != b.consensusCount)
                return a.consensusCount > b.consensusCount;
              if (a.score != b.score) return a.score > b.score;
              return a.representative.cellid < b.representative.cellid;
            });
  return consensus;
}

std::string candidateThresholdsText(const std::vector<DriftCandidate> &matches) {
  std::ostringstream os;
  for (size_t i = 0; i < matches.size(); ++i) {
    if (i > 0) os << ",";
    os << ">" << static_cast<int>(matches[i].snrThreshold);
  }
  return os.str();
}

void writeDriftCandidateRows(std::ofstream &ofs,
                             const DriftCandidate &candidate,
                             const std::vector<RunData> &runs,
                             const std::string &thresholds,
                             int consensusCount,
                             int rank,
                             double displaySnrThreshold = -1.0) {
  const double displayThreshold =
      displaySnrThreshold >= 0.0 ? displaySnrThreshold : candidate.snrThreshold;
  ofs << std::left << std::setw(7) << rank
      << std::setw(10) << candidate.cellid
      << std::setw(7) << candidate.layer
      << std::setw(7) << candidate.chip
      << std::setw(9) << candidate.channel
      << std::setw(8) << "gain"
      << std::setw(8) << candidate.type
      << std::setw(10) << (candidate.direction > 0 ? "up" : "down")
      << std::setw(12) << thresholds
      << std::right << std::setw(8) << consensusCount
      << std::setw(11) << candidate.validRuns
      << std::setw(14) << candidate.referenceGain
      << std::setw(14) << candidate.changePercent
      << std::setw(14) << candidate.slopePercentPerRun
      << std::setw(12) << candidate.spearman
      << std::setw(14) << candidate.residualMadPercent
      << std::setw(18)
      << (candidate.stepRunIndex >= 0 ? runs[candidate.stepRunIndex].label : "NA");

  for (const auto &run : runs) {
    const auto it = run.byCellId.find(candidate.cellid);
    if (it == run.byCellId.end() || !isAcceptedSpePoint(it->second) ||
        it->second.snr <= displayThreshold) {
      ofs << std::setw(18) << "NA";
    } else {
      ofs << std::setw(18) << it->second.gain;
    }
  }
  ofs << "\n";

  ofs << std::left << std::setw(7) << ""
      << std::setw(10) << ""
      << std::setw(7) << ""
      << std::setw(7) << ""
      << std::setw(9) << ""
      << std::setw(8) << "SNR"
      << std::setw(8) << ""
      << std::setw(10) << ""
      << std::setw(12) << ""
      << std::right << std::setw(8) << ""
      << std::setw(11) << ""
      << std::setw(14) << ""
      << std::setw(14) << ""
      << std::setw(14) << ""
      << std::setw(12) << ""
      << std::setw(14) << ""
      << std::setw(18) << "";
  for (const auto &run : runs) {
    const auto it = run.byCellId.find(candidate.cellid);
    if (it == run.byCellId.end() || !isAcceptedSpePoint(it->second) ||
        it->second.snr <= displayThreshold) {
      ofs << std::setw(18) << "NA";
    } else {
      ofs << std::setw(18) << it->second.snr;
    }
  }
  ofs << "\n\n";
}

void writeDriftCandidateFile(const std::string &fileName,
                             const std::vector<RunData> &runs,
                             const std::vector<DriftCandidate> &candidates,
                             double minChangePercent,
                             int minValidRuns) {
  std::ofstream ofs(fileName.c_str());
  if (!ofs) {
    std::cerr << "[compare] cannot write drift candidate file " << fileName << std::endl;
    return;
  }
  ofs << "# Good-fit SPE gain drift/step candidates\n"
      << "# Valid point: fit_ok=1, gain>0, and SNR above the candidate threshold\n"
      << "# Minimum predicted drift or persistent step change: " << minChangePercent
      << "%; minimum valid runs: " << minValidRuns << "\n"
      << "# Drift: Theil-Sen change >= threshold, |Spearman|>=0.7, residual MAD<=2%\n"
      << "# Step: >=3 valid points before/after, change >= threshold, segment MAD<=2%\n"
      << "# Two rows per candidate: gain followed by SNR\n\n";
  ofs << std::left << std::setw(7) << "rank"
      << std::setw(10) << "cellid"
      << std::setw(7) << "layer"
      << std::setw(7) << "chip"
      << std::setw(9) << "channel"
      << std::setw(8) << "metric"
      << std::setw(8) << "type"
      << std::setw(10) << "direction"
      << std::setw(12) << "SNR_cuts"
      << std::right << std::setw(8) << "nAgree"
      << std::setw(11) << "validRuns"
      << std::setw(14) << "ref_gain"
      << std::setw(14) << "change[%]"
      << std::setw(14) << "slope[%/run]"
      << std::setw(12) << "spearman"
      << std::setw(14) << "resid_MAD[%]"
      << std::setw(18) << "step_after_run";
  for (const auto &run : runs) ofs << std::setw(18) << run.label;
  ofs << "\n";
  ofs << std::fixed << std::setprecision(4);
  for (size_t i = 0; i < candidates.size(); ++i) {
    const DriftCandidate &candidate = candidates[i];
    writeDriftCandidateRows(
        ofs, candidate, runs,
        ">" + std::to_string(static_cast<int>(candidate.snrThreshold)), 1,
        static_cast<int>(i) + 1);
  }
  ofs << "# candidate count: " << candidates.size() << "\n";
}

void writeConsensusDriftFile(const std::string &fileName,
                             const std::vector<RunData> &runs,
                             const std::vector<ConsensusDriftCandidate> &consensus,
                             double minChangePercent,
                             int minValidRuns) {
  std::ofstream ofs(fileName.c_str());
  if (!ofs) {
    std::cerr << "[compare] cannot write consensus drift file " << fileName << std::endl;
    return;
  }
  ofs << "# Consensus good-fit SPE gain drift/step candidates\n"
      << "# Candidate must be detected with the same type and direction in at least two "
         "of SNR>3, SNR>4, SNR>5\n"
      << "# Step candidates must also agree on the transition position within one run "
         "group.\n"
      << "# Minimum change: " << minChangePercent
      << "%; minimum valid runs: " << minValidRuns << "\n"
      << "# Values use the highest-scoring reproduced candidate model; run values show "
         "all good points with SNR>3.\n\n";
  ofs << std::left << std::setw(7) << "rank"
      << std::setw(10) << "cellid"
      << std::setw(7) << "layer"
      << std::setw(7) << "chip"
      << std::setw(9) << "channel"
      << std::setw(8) << "metric"
      << std::setw(8) << "type"
      << std::setw(10) << "direction"
      << std::setw(12) << "SNR_cuts"
      << std::right << std::setw(8) << "nAgree"
      << std::setw(11) << "validRuns"
      << std::setw(14) << "ref_gain"
      << std::setw(14) << "change[%]"
      << std::setw(14) << "slope[%/run]"
      << std::setw(12) << "spearman"
      << std::setw(14) << "resid_MAD[%]"
      << std::setw(18) << "step_after_run";
  for (const auto &run : runs) ofs << std::setw(18) << run.label;
  ofs << "\n";
  ofs << std::fixed << std::setprecision(4);
  for (size_t i = 0; i < consensus.size(); ++i) {
    const ConsensusDriftCandidate &item = consensus[i];
    writeDriftCandidateRows(ofs, item.representative, runs,
                            candidateThresholdsText(item.matches),
                            item.consensusCount, static_cast<int>(i) + 1, 3.0);
  }
  ofs << "# consensus candidate count: " << consensus.size() << "\n";
}

void drawDriftCandidateTrend(const std::vector<RunData> &runs,
                             const ConsensusDriftCandidate &item,
                             const std::string &outDir,
                             int rank) {
  const DriftCandidate &candidate = item.representative;
  TCanvas *c = new TCanvas(Form("c_drift_%d", candidate.cellid),
                           Form("Drift candidate %d", candidate.cellid), 1300, 900);
  c->Divide(1, 2);

  std::vector<double> xGain;
  std::vector<double> yGain;
  std::vector<double> xSnr;
  std::vector<double> ySnr;
  for (size_t irun = 0; irun < runs.size(); ++irun) {
    const auto it = runs[irun].byCellId.find(candidate.cellid);
    if (it == runs[irun].byCellId.end() || !isNonExcludedFitPoint(it->second)) {
      continue;
    }
    xSnr.push_back(irun + 1.0);
    ySnr.push_back(it->second.snr);
    if (it->second.snr > 3.0) {
      xGain.push_back(irun + 1.0);
      yGain.push_back(it->second.gain);
    }
  }
  if (xGain.empty() || xSnr.empty()) return;

  c->cd(1);
  gPad->SetGrid();
  gPad->SetBottomMargin(0.18);
  const SummaryStats gainStats = computeStats(yGain);
  const double gainSpan = std::max(0.5, gainStats.max - gainStats.min);
  TH1D *gainFrame = new TH1D(
      Form("h_drift_gain_frame_%d", candidate.cellid),
      Form("Rank %d cell %d: %s %s %.2f%% (SNR%s);;Gain",
           rank, candidate.cellid, candidate.type.c_str(),
           candidate.direction > 0 ? "up" : "down", candidate.changePercent,
           candidateThresholdsText(item.matches).c_str()),
      static_cast<int>(runs.size()), 0.5, runs.size() + 0.5);
  for (size_t irun = 0; irun < runs.size(); ++irun) {
    gainFrame->GetXaxis()->SetBinLabel(irun + 1, runs[irun].label.c_str());
  }
  gainFrame->GetXaxis()->LabelsOption("v");
  gainFrame->SetMinimum(gainStats.min - 0.2 * gainSpan);
  gainFrame->SetMaximum(gainStats.max + 0.2 * gainSpan);
  gainFrame->SetStats(0);
  gainFrame->Draw();
  TGraph *gainGraph = new TGraph(xGain.size(), &xGain[0], &yGain[0]);
  setGraphStyle(gainGraph, kBlack, 20);
  gainGraph->Draw("LP SAME");

  std::vector<double> modelX;
  std::vector<double> modelY;
  if (candidate.type == "step") {
    modelX = {1.0, candidate.stepX + 1.0, candidate.stepX + 1.0,
              static_cast<double>(runs.size())};
    modelY = {
        candidate.referenceGain * (1.0 + candidate.stepBeforePercent / 100.0),
        candidate.referenceGain * (1.0 + candidate.stepBeforePercent / 100.0),
        candidate.referenceGain * (1.0 + candidate.stepAfterPercent / 100.0),
        candidate.referenceGain * (1.0 + candidate.stepAfterPercent / 100.0)};
  } else {
    modelX = {1.0, static_cast<double>(runs.size())};
    modelY = {
        candidate.referenceGain *
            (1.0 + candidate.lineInterceptPercent / 100.0),
        candidate.referenceGain *
            (1.0 + (candidate.lineInterceptPercent +
                    candidate.slopePercentPerRun * (runs.size() - 1)) /
                       100.0)};
  }
  TGraph *model = new TGraph(modelX.size(), &modelX[0], &modelY[0]);
  model->SetLineColor(kRed + 1);
  model->SetLineWidth(3);
  model->Draw("L SAME");

  c->cd(2);
  gPad->SetGrid();
  gPad->SetBottomMargin(0.18);
  const SummaryStats snrStats = computeStats(ySnr);
  TH1D *snrFrame = new TH1D(
      Form("h_drift_snr_frame_%d", candidate.cellid),
      "SNR used to validate gain calculation;;SNR",
      static_cast<int>(runs.size()), 0.5, runs.size() + 0.5);
  for (size_t irun = 0; irun < runs.size(); ++irun) {
    snrFrame->GetXaxis()->SetBinLabel(irun + 1, runs[irun].label.c_str());
  }
  snrFrame->GetXaxis()->LabelsOption("v");
  snrFrame->SetMinimum(0.0);
  snrFrame->SetMaximum(std::max(6.0, snrStats.max * 1.2));
  snrFrame->SetStats(0);
  snrFrame->Draw();
  TGraph *snrGraph = new TGraph(xSnr.size(), &xSnr[0], &ySnr[0]);
  setGraphStyle(snrGraph, kBlue + 1, 21);
  snrGraph->Draw("LP SAME");
  for (int threshold = 3; threshold <= 5; ++threshold) {
    TLine *line = new TLine(0.5, threshold, runs.size() + 0.5, threshold);
    line->SetLineColor(threshold == 3 ? kRed + 1
                                     : (threshold == 4 ? kOrange + 7 : kGreen + 2));
    line->SetLineStyle(2);
    line->Draw();
  }

  const std::string prefix =
      outDir + Form("/rank%03d_cell%d_%s", rank, candidate.cellid,
                    candidate.type.c_str());
  c->SaveAs((prefix + ".pdf").c_str());
  c->SaveAs((prefix + ".png").c_str());
}

void runDriftAnalysis(const std::vector<RunData> &runs,
                      const std::string &outDir,
                      double minChangePercent = 3.0,
                      int minValidRuns = 8,
                      int topPlots = 50) {
  std::vector<std::map<int, DriftCandidate> > candidateSets;
  for (double threshold : kDriftSnrThresholds) {
    std::map<int, DriftCandidate> candidateMap =
        findDriftCandidates(runs, threshold, minChangePercent, minValidRuns);
    std::vector<DriftCandidate> candidates;
    for (const auto &kv : candidateMap) candidates.push_back(kv.second);
    std::sort(candidates.begin(), candidates.end(),
              [](const DriftCandidate &a, const DriftCandidate &b) {
                if (a.score != b.score) return a.score > b.score;
                return a.cellid < b.cellid;
              });
    writeDriftCandidateFile(
        outDir + Form("/drift_candidates_snr_gt%d.txt", static_cast<int>(threshold)),
        runs, candidates, minChangePercent, minValidRuns);
    std::cout << "[compare] SNR > " << threshold << ": " << candidates.size()
              << " drift/step candidates" << std::endl;
    candidateSets.push_back(candidateMap);
  }

  const std::vector<ConsensusDriftCandidate> consensus =
      buildDriftConsensus(candidateSets, kDriftConsensusMin);
  writeConsensusDriftFile(outDir + "/drift_candidates_consensus.txt", runs,
                          consensus, minChangePercent, minValidRuns);
  std::cout << "[compare] consensus drift/step candidates: " << consensus.size()
            << std::endl;

  const std::string plotDir = outDir + "/drift_top50";
  gSystem->mkdir(plotDir.c_str(), true);
  const int nPlots =
      std::min(std::max(0, topPlots), static_cast<int>(consensus.size()));
  for (int i = 0; i < nPlots; ++i) {
    drawDriftCandidateTrend(runs, consensus[i], plotDir, i + 1);
  }
}

std::vector<StableChannelInfo> rankStableChannels(const std::vector<RunData> &runs,
                                                  double snrThreshold = 5.0,
                                                  int minValidRuns = 8) {
  std::set<int> allCellIds;
  for (const auto &run : runs) {
    for (const auto &kv : run.byCellId) allCellIds.insert(kv.first);
  }

  std::vector<StableChannelInfo> ranked;
  for (int cellid : allCellIds) {
    const SpeEntry *identity = nullptr;
    std::vector<double> x;
    std::vector<double> gains;
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      const auto it = runs[irun].byCellId.find(cellid);
      if (it == runs[irun].byCellId.end()) continue;
      if (!identity) identity = &it->second;
      const SpeEntry &entry = it->second;
      if (!isAcceptedSpePoint(entry) || entry.snr <= snrThreshold) continue;
      x.push_back(irun);
      gains.push_back(entry.gain);
    }
    if (!identity || static_cast<int>(gains.size()) < minValidRuns) continue;

    StableChannelInfo info;
    info.cellid = cellid;
    info.layer = identity->layer;
    info.chip = identity->chip;
    info.channel = identity->channel;
    info.validRuns = static_cast<int>(gains.size());
    info.referenceGain = computeMedian(gains);
    if (info.referenceGain <= 0.0) continue;
    for (double gain : gains) {
      info.deltasPercent.push_back(
          100.0 * (gain - info.referenceGain) / info.referenceGain);
    }
    info.deltaMadPercent = medianAbsoluteDeviation(info.deltasPercent, 0.0);
    info.deltaRmsPercent = computeStats(info.deltasPercent).rms;
    for (double delta : info.deltasPercent) {
      info.maxAbsDeltaPercent = std::max(info.maxAbsDeltaPercent, std::abs(delta));
    }
    double intercept = 0.0;
    double residualMad = 0.0;
    fitTheilSen(x, info.deltasPercent, info.slopePercentPerRun, intercept, residualMad);
    info.spearman = computeSpearmanCorrelation(x, info.deltasPercent);
    ranked.push_back(info);
  }

  std::sort(ranked.begin(), ranked.end(),
            [](const StableChannelInfo &a, const StableChannelInfo &b) {
              if (a.validRuns != b.validRuns) return a.validRuns > b.validRuns;
              if (a.deltaRmsPercent != b.deltaRmsPercent)
                return a.deltaRmsPercent < b.deltaRmsPercent;
              if (a.maxAbsDeltaPercent != b.maxAbsDeltaPercent)
                return a.maxAbsDeltaPercent < b.maxAbsDeltaPercent;
              if (a.deltaMadPercent != b.deltaMadPercent)
                return a.deltaMadPercent < b.deltaMadPercent;
              return a.cellid < b.cellid;
            });
  return ranked;
}

void writeStableChannelRanking(const std::string &fileName,
                               const std::vector<RunData> &runs,
                               const std::vector<StableChannelInfo> &ranked,
                               double snrThreshold) {
  std::ofstream ofs(fileName.c_str());
  if (!ofs) return;
  ofs << "# Stable SPE gain channels ranked by valid run count, delta RMS, maximum "
         "absolute delta, and MAD\n"
      << "# Valid point: fit_ok=1, gain>0, SNR>" << snrThreshold << "\n"
      << "# delta[%] = 100 * (gain - per-channel median gain) / median gain\n\n";
  ofs << std::left << std::setw(7) << "rank"
      << std::setw(10) << "cellid"
      << std::setw(7) << "layer"
      << std::setw(7) << "chip"
      << std::setw(9) << "channel"
      << std::right << std::setw(11) << "validRuns"
      << std::setw(14) << "ref_gain"
      << std::setw(14) << "MAD[%]"
      << std::setw(14) << "RMS[%]"
      << std::setw(16) << "max_abs[%]"
      << std::setw(16) << "slope[%/run]"
      << std::setw(12) << "spearman";
  for (const auto &run : runs) ofs << std::setw(18) << run.label;
  ofs << "\n" << std::fixed << std::setprecision(4);

  for (size_t i = 0; i < ranked.size(); ++i) {
    const StableChannelInfo &info = ranked[i];
    ofs << std::left << std::setw(7) << i + 1
        << std::setw(10) << info.cellid
        << std::setw(7) << info.layer
        << std::setw(7) << info.chip
        << std::setw(9) << info.channel
        << std::right << std::setw(11) << info.validRuns
        << std::setw(14) << info.referenceGain
        << std::setw(14) << info.deltaMadPercent
        << std::setw(14) << info.deltaRmsPercent
        << std::setw(16) << info.maxAbsDeltaPercent
        << std::setw(16) << info.slopePercentPerRun
        << std::setw(12) << info.spearman;
    for (const auto &run : runs) {
      const auto it = run.byCellId.find(info.cellid);
      if (it == run.byCellId.end() || !isAcceptedSpePoint(it->second) ||
          it->second.snr <= snrThreshold) {
        ofs << std::setw(18) << "NA";
      } else {
        ofs << std::setw(18) << it->second.gain;
      }
    }
    ofs << "\n";
  }
}

void drawStableChannelDistributions(const std::vector<RunData> &runs,
                                    const std::string &outDir,
                                    double snrThreshold = 5.0,
                                    int minValidRuns = 8,
                                    int topChannels = 10) {
  const std::vector<StableChannelInfo> ranked =
      rankStableChannels(runs, snrThreshold, minValidRuns);
  if (ranked.empty()) return;

  writeStableChannelRanking(outDir + "/stable_channels_ranked.txt", runs, ranked,
                            snrThreshold);

  {
    TCanvas *c =
        new TCanvas("c_stable_channel_metric_distributions",
                    "Stable channel metric distributions", 1400, 900);
    c->Divide(2, 2);
    TH1D *hMad =
        new TH1D("h_stable_channel_mad", "Channel gain stability;Delta MAD [%];Channels",
                 120, 0.0, 6.0);
    TH1D *hRms =
        new TH1D("h_stable_channel_rms", "Channel gain stability;Delta RMS [%];Channels",
                 120, 0.0, 12.0);
    TH1D *hMax = new TH1D("h_stable_channel_max_abs",
                          "Channel gain stability;Maximum |delta| [%];Channels",
                          120, 0.0, 30.0);
    TH1D *hValid =
        new TH1D("h_stable_channel_valid_runs",
                 "Channel gain stability;Valid run groups;Channels",
                 runs.size(), 0.5, runs.size() + 0.5);
    for (const auto &info : ranked) {
      hMad->Fill(info.deltaMadPercent);
      hRms->Fill(info.deltaRmsPercent);
      hMax->Fill(info.maxAbsDeltaPercent);
      hValid->Fill(info.validRuns);
    }
    c->cd(1); gPad->SetGrid(); hMad->Draw("HIST");
    c->cd(2); gPad->SetGrid(); hRms->Draw("HIST");
    c->cd(3); gPad->SetGrid(); hMax->Draw("HIST");
    c->cd(4); gPad->SetGrid(); hValid->Draw("HIST");
    c->SaveAs((outDir + "/stable_channel_metric_distributions.pdf").c_str());
    c->SaveAs((outDir + "/stable_channel_metric_distributions.png").c_str());
  }

  const int nTop = std::min(std::max(1, topChannels), static_cast<int>(ranked.size()));
  {
    TCanvas *c =
        new TCanvas("c_most_stable_channel_delta_distributions",
                    "Most stable channel delta distributions", 1100, 750);
    c->SetGrid();
    TLegend *leg = new TLegend(0.62, 0.55, 0.90, 0.89);
    int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7,
                    kCyan + 2, kRed + 3, kBlue + 3, kGreen + 4, kGray + 2};
    double ymax = 0.0;
    std::vector<TH1D *> histograms;
    for (int i = 0; i < nTop; ++i) {
      TH1D *h = new TH1D(
          Form("h_stable_delta_rank%d_cell%d", i + 1, ranked[i].cellid),
          Form("Most stable channels (SNR > %.0f);Delta gain / median gain [%%];"
               "Normalized entries",
               snrThreshold),
          80, -5.0, 5.0);
      for (double delta : ranked[i].deltasPercent) h->Fill(delta);
      if (h->Integral() > 0.0) h->Scale(1.0 / h->Integral());
      h->SetLineColor(colors[i % 10]);
      h->SetLineWidth(2);
      h->SetStats(0);
      ymax = std::max(ymax, h->GetMaximum() * 1.25);
      h->Draw(i == 0 ? "HIST" : "HIST SAME");
      histograms.push_back(h);
      leg->AddEntry(h, Form("#%d cell %d (MAD=%.4f%%)", i + 1, ranked[i].cellid,
                            ranked[i].deltaMadPercent), "l");
    }
    for (TH1D *h : histograms) h->SetMaximum(ymax);
    leg->Draw();
    c->SaveAs((outDir + "/most_stable_channels_delta_distribution.pdf").c_str());
    c->SaveAs((outDir + "/most_stable_channels_delta_distribution.png").c_str());
  }

  const StableChannelInfo &best = ranked.front();
  TCanvas *c =
      new TCanvas("c_most_stable_channel_trend", "Most stable channel", 1300, 850);
  c->Divide(1, 2);
  std::vector<double> x;
  std::vector<double> gain;
  std::vector<double> snr;
  for (size_t irun = 0; irun < runs.size(); ++irun) {
    const auto it = runs[irun].byCellId.find(best.cellid);
    if (it == runs[irun].byCellId.end() || !isAcceptedSpePoint(it->second) ||
        it->second.snr <= snrThreshold) continue;
    x.push_back(irun + 1.0);
    gain.push_back(it->second.gain);
    snr.push_back(it->second.snr);
  }
  c->cd(1); gPad->SetGrid();
  TGraph *gGain = new TGraph(x.size(), &x[0], &gain[0]);
  setGraphStyle(gGain, kBlack, 20);
  gGain->SetTitle(Form("Most stable cell %d: MAD=%.4f%% RMS=%.4f%%;Run group;Gain",
                       best.cellid, best.deltaMadPercent, best.deltaRmsPercent));
  gGain->Draw("ALP");
  c->cd(2); gPad->SetGrid();
  TGraph *gSnr = new TGraph(x.size(), &x[0], &snr[0]);
  setGraphStyle(gSnr, kBlue + 1, 21);
  gSnr->SetTitle(Form("Most stable cell %d;Run group;SNR", best.cellid));
  gSnr->Draw("ALP");
  c->SaveAs((outDir + "/most_stable_channel_trend.pdf").c_str());
  c->SaveAs((outDir + "/most_stable_channel_trend.png").c_str());

  std::cout << "[compare] most stable channel: " << best.cellid
            << " validRuns=" << best.validRuns << " MAD=" << best.deltaMadPercent
            << "% RMS=" << best.deltaRmsPercent << "%" << std::endl;
}

std::vector<GainRunRmsInfo> buildGainRunRmsInfo(const std::vector<RunData> &runs,
                                                int minValidRuns = 2) {
  const std::set<int> selectedCells = selectCellsBySurvivalRate(runs, 0.50);

  std::vector<GainRunRmsInfo> infos;
  for (int cellid : selectedCells) {
    const SpeEntry *identity = nullptr;
    std::vector<double> gains;
    gains.reserve(runs.size());

    for (const auto &run : runs) {
      const auto it = run.byCellId.find(cellid);
      if (it == run.byCellId.end()) continue;
      if (!identity) identity = &it->second;
      const SpeEntry &entry = it->second;
      if (!isSurvivingSpePlotPoint(entry)) continue;
      gains.push_back(entry.gain);
    }

    if (!identity || static_cast<int>(gains.size()) < minValidRuns) continue;
    const SummaryStats gainStats = computeStats(gains);

    GainRunRmsInfo info;
    info.cellid = cellid;
    info.layer = entryLayer(*identity);
    info.chip = entryChip(*identity);
    info.channel = entryChannel(*identity);
    info.validRuns = static_cast<int>(gains.size());
    info.meanGain = gainStats.mean;
    info.gainRms = gainStats.rms;
    info.minGain = gainStats.min;
    info.maxGain = gainStats.max;
    info.gainValues = gains;
    infos.push_back(info);
  }

  std::sort(infos.begin(), infos.end(),
            [](const GainRunRmsInfo &a, const GainRunRmsInfo &b) {
              if (a.gainRms != b.gainRms) return a.gainRms > b.gainRms;
              if (a.validRuns != b.validRuns) return a.validRuns > b.validRuns;
              return a.cellid < b.cellid;
            });
  return infos;
}

void writeGainRunRmsReport(const std::string &fileName,
                           const std::vector<RunData> &runs,
                           const std::vector<GainRunRmsInfo> &infos,
                           int minValidRuns) {
  std::ofstream ofs(fileName.c_str());
  if (!ofs) return;

  ofs << "# Simple RMS of the gain distribution for each channel\n"
      << "# Selected channel: survival rate > 50% across loaded runs\n"
      << "# Survival/valid point: fit_ok=1, not fitting-mis-excluded, SNR>5, gain>0, finite gain\n"
      << "# gain_RMS is the standard RMS/spread of gain values across runs for one channel\n"
      << "# Per-run columns contain gain values\n"
      << "# Minimum valid runs: " << minValidRuns << "\n\n";

  ofs << std::left << std::setw(7) << "rank"
      << std::setw(10) << "cellid"
      << std::setw(7) << "layer"
      << std::setw(7) << "chip"
      << std::setw(9) << "channel"
      << std::right << std::setw(11) << "validRuns"
      << std::setw(14) << "mean_gain"
      << std::setw(14) << "gain_RMS"
      << std::setw(14) << "min_gain"
      << std::setw(14) << "max_gain";
  for (const auto &run : runs) ofs << std::setw(18) << run.label;
  ofs << "\n" << std::fixed << std::setprecision(4);

  for (size_t i = 0; i < infos.size(); ++i) {
    const GainRunRmsInfo &info = infos[i];
    ofs << std::left << std::setw(7) << i + 1
        << std::setw(10) << info.cellid
        << std::setw(7) << info.layer
        << std::setw(7) << info.chip
        << std::setw(9) << info.channel
        << std::right << std::setw(11) << info.validRuns
        << std::setw(14) << info.meanGain
        << std::setw(14) << info.gainRms
        << std::setw(14) << info.minGain
        << std::setw(14) << info.maxGain;
    for (const auto &run : runs) {
      const auto it = run.byCellId.find(info.cellid);
      if (it == run.byCellId.end() || !isSurvivingSpePlotPoint(it->second)) {
        ofs << std::setw(18) << "NA";
      } else {
        ofs << std::setw(18) << it->second.gain;
      }
    }
    ofs << "\n";
  }
}

void drawGainRunRmsDistribution(const std::vector<RunData> &runs,
                                const std::string &outDir,
                                int minValidRuns = 2) {
  if (runs.size() < 2) return;

  const std::vector<GainRunRmsInfo> infos =
      buildGainRunRmsInfo(runs, minValidRuns);
  if (infos.empty()) {
    std::cout << "[compare] no channels with survival rate > 50% for gain RMS"
              << std::endl;
    return;
  }

  writeGainRunRmsReport(outDir + "/gain_rms_by_channel.txt", runs, infos,
                        minValidRuns);

  std::vector<double> gainRmsValues;
  gainRmsValues.reserve(infos.size());
  for (const auto &info : infos) {
    gainRmsValues.push_back(info.gainRms);
  }

  const SummaryStats rmsStats = computeStats(gainRmsValues);
  const double rmsMax = rmsStats.max > 0.0 ? rmsStats.max * 1.1 : 1.0;

  TCanvas *c =
      new TCanvas("c_gain_rms_distribution",
                  "Channel-by-channel gain RMS distribution", 900, 700);
  c->SetGrid();

  TH1D *h = new TH1D("h_gain_rms_by_channel",
                     "Gain RMS by channel;RMS of gain across runs;Channels",
                     120, 0.0, rmsMax);

  for (const auto &info : infos) {
    h->Fill(info.gainRms);
  }

  h->SetLineColor(kRed + 1);
  h->SetLineWidth(2);
  h->Draw("HIST");

  c->SaveAs((outDir + "/gain_rms_distribution.pdf").c_str());
  c->SaveAs((outDir + "/gain_rms_distribution.png").c_str());

  std::cout << "[compare] wrote channel gain RMS distribution for "
            << infos.size() << " channels" << std::endl;
}

void drawHighGainRmsChannelTrends(const std::vector<RunData> &runs,
                                  const std::string &outDir,
                                  int topPlots = 50,
                                  int minValidRuns = 2) {
  if (runs.size() < 2 || topPlots <= 0) return;

  const std::vector<GainRunRmsInfo> infos =
      buildGainRunRmsInfo(runs, minValidRuns);
  if (infos.empty()) {
    std::cout << "[compare] no channels with survival rate > 50% for gain RMS trends"
              << std::endl;
    return;
  }

  const std::string plotDir = outDir + "/gain_rms_trends";
  gSystem->mkdir(plotDir.c_str(), true);
  const int nPlots = std::min(topPlots, static_cast<int>(infos.size()));

  for (int i = 0; i < nPlots; ++i) {
    const GainRunRmsInfo &info = infos[i];
    std::vector<double> xGain;
    std::vector<double> yGain;
    std::vector<double> xSnr;
    std::vector<double> ySnr;
    xGain.reserve(runs.size());
    yGain.reserve(runs.size());
    xSnr.reserve(runs.size());
    ySnr.reserve(runs.size());

    for (size_t irun = 0; irun < runs.size(); ++irun) {
      const auto it = runs[irun].byCellId.find(info.cellid);
      if (it == runs[irun].byCellId.end() || !isSurvivingSpePlotPoint(it->second)) {
        continue;
      }
      xGain.push_back(irun + 1.0);
      yGain.push_back(it->second.gain);
      xSnr.push_back(irun + 1.0);
      ySnr.push_back(it->second.snr);
    }
    if (xGain.empty() || xSnr.empty()) continue;

    TCanvas *c =
        new TCanvas(Form("c_gain_rms_trend_%d", info.cellid),
                    Form("High gain RMS channel %d", info.cellid), 1300, 900);
    c->Divide(1, 2);

    c->cd(1);
    gPad->SetGrid();
    gPad->SetBottomMargin(0.18);
    const SummaryStats gainStats = computeStats(yGain);
    const double gainSpan = std::max(0.5, gainStats.max - gainStats.min);
    TH1D *gainFrame = new TH1D(
        Form("h_gain_rms_frame_%d", info.cellid),
        Form("Rank %d cell %d: gain RMS %.4f, mean %.4f;;Gain",
             i + 1, info.cellid, info.gainRms, info.meanGain),
        static_cast<int>(runs.size()), 0.5, runs.size() + 0.5);
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      gainFrame->GetXaxis()->SetBinLabel(irun + 1, runs[irun].label.c_str());
    }
    gainFrame->GetXaxis()->LabelsOption("v");
    gainFrame->SetMinimum(gainStats.min - 0.2 * gainSpan);
    gainFrame->SetMaximum(gainStats.max + 0.2 * gainSpan);
    gainFrame->SetStats(0);
    gainFrame->Draw();
    TGraph *gainGraph = new TGraph(xGain.size(), &xGain[0], &yGain[0]);
    setGraphStyle(gainGraph, kBlack, 20);
    gainGraph->Draw("LP SAME");
    TLine *meanLine = new TLine(0.5, info.meanGain, runs.size() + 0.5, info.meanGain);
    meanLine->SetLineColor(kRed + 1);
    meanLine->SetLineStyle(2);
    meanLine->SetLineWidth(2);
    meanLine->Draw();

    c->cd(2);
    gPad->SetGrid();
    gPad->SetBottomMargin(0.18);
    const SummaryStats snrStats = computeStats(ySnr);
    TH1D *snrFrame = new TH1D(
        Form("h_gain_rms_snr_frame_%d", info.cellid),
        "SNR used for gain RMS selection;;SNR",
        static_cast<int>(runs.size()), 0.5, runs.size() + 0.5);
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      snrFrame->GetXaxis()->SetBinLabel(irun + 1, runs[irun].label.c_str());
    }
    snrFrame->GetXaxis()->LabelsOption("v");
    snrFrame->SetMinimum(0.0);
    snrFrame->SetMaximum(std::max(6.0, snrStats.max * 1.2));
    snrFrame->SetStats(0);
    snrFrame->Draw();
    TGraph *snrGraph = new TGraph(xSnr.size(), &xSnr[0], &ySnr[0]);
    setGraphStyle(snrGraph, kBlue + 1, 21);
    snrGraph->Draw("LP SAME");
    TLine *snrLine = new TLine(0.5, 5.0, runs.size() + 0.5, 5.0);
    snrLine->SetLineColor(kGreen + 2);
    snrLine->SetLineStyle(2);
    snrLine->Draw();

    const std::string prefix =
        plotDir + Form("/rank%03d_cell%d_gain_rms", i + 1, info.cellid);
    c->SaveAs((prefix + ".pdf").c_str());
    c->SaveAs((prefix + ".png").c_str());
    delete c;
  }

  std::cout << "[compare] wrote top " << nPlots
            << " high gain-RMS channel trends to " << plotDir << std::endl;
}

void writeAndDrawLayerMeanGainShift(const std::vector<RunData> &runs,
                                    const std::string &outDir) {
  if (runs.size() < 2) return;

  const std::set<int> selectedCells = selectCellsBySurvivalRate(runs, 0.50);
  if (selectedCells.empty()) {
    std::cout << "[compare] no channels with survival rate > 50% for layer mean gain shift"
              << std::endl;
    return;
  }

  std::vector<std::map<int, std::vector<double> > > gainsByRunLayer(runs.size());
  std::set<int> layers;
  for (size_t irun = 0; irun < runs.size(); ++irun) {
    for (const auto &entry : runs[irun].entries) {
      if (selectedCells.count(entry.cellid) == 0) continue;
      if (!isSurvivingSpePlotPoint(entry)) continue;
      const int layer = entryLayer(entry);
      if (layer < 0) continue;
      gainsByRunLayer[irun][layer].push_back(entry.gain);
      layers.insert(layer);
    }
  }
  if (layers.empty()) return;

  std::map<int, std::vector<double> > meanByLayer;
  std::map<int, std::vector<int> > countByLayer;
  std::map<int, double> referenceByLayer;
  for (int layer : layers) {
    meanByLayer[layer].assign(runs.size(), 0.0);
    countByLayer[layer].assign(runs.size(), 0);
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      const auto it = gainsByRunLayer[irun].find(layer);
      if (it == gainsByRunLayer[irun].end() || it->second.empty()) continue;
      const SummaryStats stats = computeStats(it->second);
      meanByLayer[layer][irun] = stats.mean;
      countByLayer[layer][irun] = stats.n;
      if (referenceByLayer.count(layer) == 0 && stats.mean > 0.0) {
        referenceByLayer[layer] = stats.mean;
      }
    }
  }

  std::ofstream ofs((outDir + "/layer_mean_gain_shift.txt").c_str());
  if (ofs) {
    ofs << "# Layer mean gain shift\n"
        << "# Selected channel: survival rate > 50% across loaded runs\n"
        << "# Survival/valid point: fit_ok=1, not fitting-mis-excluded, SNR>5, gain>0, finite gain\n"
        << "# shift[%] = 100 * (layer_mean_gain - first_valid_layer_mean_gain) / first_valid_layer_mean_gain\n\n";
    ofs << std::left << std::setw(20) << "run"
        << std::setw(8) << "layer"
        << std::right << std::setw(12) << "nCell"
        << std::setw(16) << "mean_gain"
        << std::setw(16) << "ref_gain"
        << std::setw(16) << "shift[%]" << "\n";
    ofs << std::fixed << std::setprecision(4);
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      for (int layer : layers) {
        const int nCell = countByLayer[layer][irun];
        if (nCell <= 0 || referenceByLayer.count(layer) == 0 ||
            referenceByLayer[layer] <= 0.0) {
          continue;
        }
        const double meanGain = meanByLayer[layer][irun];
        const double refGain = referenceByLayer[layer];
        const double shift = 100.0 * (meanGain - refGain) / refGain;
        ofs << std::left << std::setw(20) << runs[irun].label
            << std::setw(8) << layer
            << std::right << std::setw(12) << nCell
            << std::setw(16) << meanGain
            << std::setw(16) << refGain
            << std::setw(16) << shift << "\n";
      }
    }
  }

  TCanvas *c = new TCanvas("c_layer_mean_gain_shift",
                           "Layer mean gain shift", 1400, 900);
  c->SetGrid();
  c->SetBottomMargin(0.22);
  c->SetRightMargin(0.24);
  TLegend *leg = new TLegend(0.80, 0.12, 0.98, 0.88);
  leg->SetTextSize(0.022);
  if (layers.size() > 24) {
    leg->SetNColumns(3);
  } else if (layers.size() > 12) {
    leg->SetNColumns(2);
  }
  int colors[] = {kRed + 1,     kBlue + 1,    kGreen + 2,  kMagenta + 1,
                  kOrange + 7,  kCyan + 2,    kRed + 3,    kBlue + 3,
                  kGreen + 4,   kMagenta + 3, kOrange + 9, kCyan + 4};
  int markers[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33};

  double ymin = 1e30;
  double ymax = -1e30;
  int drawn = 0;
  std::vector<TGraph *> graphs;
  struct LayerShiftLabel {
    int layer = -1;
    double x = 0.0;
    double y = 0.0;
    double maxAbsShift = 0.0;
    int color = kBlack;
  };
  std::vector<LayerShiftLabel> labelCandidates;
  for (int layer : layers) {
    if (referenceByLayer.count(layer) == 0 || referenceByLayer[layer] <= 0.0) continue;
    std::vector<double> x;
    std::vector<double> y;
    LayerShiftLabel label;
    label.layer = layer;
    label.color = colors[drawn % 12];
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      if (countByLayer[layer][irun] <= 0) continue;
      const double shift =
          100.0 * (meanByLayer[layer][irun] - referenceByLayer[layer]) /
          referenceByLayer[layer];
      x.push_back(irun + 1.0);
      y.push_back(shift);
      ymin = std::min(ymin, shift);
      ymax = std::max(ymax, shift);
      if (std::abs(shift) > label.maxAbsShift) {
        label.maxAbsShift = std::abs(shift);
        label.x = irun + 1.0;
        label.y = shift;
      }
    }
    if (x.empty()) continue;
    TGraph *g = new TGraph(x.size(), &x[0], &y[0]);
    setGraphStyle(g, colors[drawn % 12], markers[drawn % 12]);
    graphs.push_back(g);
    leg->AddEntry(g, Form("Layer %d", layer), "lp");
    labelCandidates.push_back(label);
    ++drawn;
  }

  if (drawn > 0) {
    const double span = std::max(0.5, ymax - ymin);
    TH1D *frame = new TH1D("h_layer_mean_gain_shift_frame",
                           "Layer mean gain shift;;#Delta layer mean gain [%]",
                           static_cast<int>(runs.size()), 0.5, runs.size() + 0.5);
    const int labelStride =
        runs.size() > 20 ? std::max(1, static_cast<int>(runs.size() / 18)) : 1;
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      const bool showLabel = irun == 0 || irun + 1 == runs.size() ||
                             (static_cast<int>(irun) % labelStride == 0);
      frame->GetXaxis()->SetBinLabel(irun + 1,
                                     showLabel ? runs[irun].label.c_str() : "");
    }
    frame->GetXaxis()->LabelsOption("v");
    frame->GetXaxis()->SetLabelSize(0.026);
    frame->GetYaxis()->SetLabelSize(0.032);
    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->SetMinimum(ymin - 0.2 * span);
    frame->SetMaximum(ymax + 0.2 * span);
    frame->SetStats(0);
    frame->Draw();
    for (TGraph *g : graphs) g->Draw("LP SAME");
    TLine *zero = new TLine(0.5, 0.0, runs.size() + 0.5, 0.0);
    zero->SetLineColor(kBlack);
    zero->SetLineStyle(2);
    zero->Draw();
    std::sort(labelCandidates.begin(), labelCandidates.end(),
              [](const LayerShiftLabel &a, const LayerShiftLabel &b) {
                return a.maxAbsShift > b.maxAbsShift;
              });
    const int nLabels =
        std::min(7, static_cast<int>(labelCandidates.size()));
    for (int i = 0; i < nLabels; ++i) {
      const LayerShiftLabel &label = labelCandidates[i];
      const double dx = runs.size() > 1 ? 0.35 : 0.0;
      const bool rightSide = label.x > 0.72 * runs.size();
      const double textX = rightSide ? label.x - dx : label.x + dx;
      const double textY = label.y + (label.y >= 0.0 ? 0.06 : -0.06) * span;
      TLatex *latex = new TLatex(textX, textY, Form("L%d", label.layer));
      latex->SetTextSize(0.026);
      latex->SetTextColor(label.color);
      latex->SetTextAlign(rightSide ? 32 : 12);
      latex->Draw();
    }
    leg->Draw();
    c->SaveAs((outDir + "/layer_mean_gain_shift.pdf").c_str());
    c->SaveAs((outDir + "/layer_mean_gain_shift.png").c_str());
  }
  delete c;

  std::cout << "[compare] wrote layer mean gain shift" << std::endl;
}

void writeAndDrawLayerMeanGainRmsTrends(const std::vector<RunData> &runs,
                                        const std::string &outDir) {
  if (runs.empty()) return;

  const std::set<int> selectedCells = selectCellsBySurvivalRate(runs, 0.50);
  if (selectedCells.empty()) {
    std::cout << "[compare] no channels with survival rate > 50% for layer mean gain RMS trends"
              << std::endl;
    return;
  }

  std::vector<std::map<int, std::vector<double> > > gainsByRunLayer(runs.size());
  std::set<int> layers;
  for (size_t irun = 0; irun < runs.size(); ++irun) {
    for (const auto &entry : runs[irun].entries) {
      if (selectedCells.count(entry.cellid) == 0) continue;
      if (!isSurvivingSpePlotPoint(entry)) continue;
      const int layer = entryLayer(entry);
      if (layer < 0) continue;
      gainsByRunLayer[irun][layer].push_back(entry.gain);
      layers.insert(layer);
    }
  }
  if (layers.empty()) return;

  const std::string plotDir = outDir + "/layer_mean_gain_rms_trends";
  gSystem->mkdir(plotDir.c_str(), true);

  std::ofstream ofs((plotDir + "/layer_mean_gain_rms_trends.txt").c_str());
  if (ofs) {
    ofs << "# Layer mean gain trend with layer RMS as y error\n"
        << "# Selected channel: survival rate > 50% across loaded runs\n"
        << "# Survival/valid point: fit_ok=1, not fitting-mis-excluded, SNR>5, gain>0, finite gain\n\n";
    ofs << std::left << std::setw(20) << "run"
        << std::setw(8) << "layer"
        << std::right << std::setw(12) << "nCell"
        << std::setw(16) << "mean_gain"
        << std::setw(16) << "gain_RMS" << "\n";
    ofs << std::fixed << std::setprecision(4);
  }

  for (int layer : layers) {
    std::vector<double> x;
    std::vector<double> ex;
    std::vector<double> y;
    std::vector<double> ey;
    x.reserve(runs.size());
    ex.reserve(runs.size());
    y.reserve(runs.size());
    ey.reserve(runs.size());

    double ymin = 1e30;
    double ymax = -1e30;
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      const auto it = gainsByRunLayer[irun].find(layer);
      if (it == gainsByRunLayer[irun].end() || it->second.empty()) continue;
      const SummaryStats stats = computeStats(it->second);
      x.push_back(irun + 1.0);
      ex.push_back(0.0);
      y.push_back(stats.mean);
      ey.push_back(stats.rms);
      ymin = std::min(ymin, stats.mean - stats.rms);
      ymax = std::max(ymax, stats.mean + stats.rms);
      if (ofs) {
        ofs << std::left << std::setw(20) << runs[irun].label
            << std::setw(8) << layer
            << std::right << std::setw(12) << stats.n
            << std::setw(16) << stats.mean
            << std::setw(16) << stats.rms << "\n";
      }
    }
    if (x.empty()) continue;

    TCanvas *c = new TCanvas(Form("c_layer_mean_gain_rms_%d", layer),
                             Form("Layer %d mean gain RMS trend", layer),
                             1200, 760);
    c->SetGrid();
    c->SetBottomMargin(0.24);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.06);

    const double span = std::max(0.5, ymax - ymin);
    TH1D *frame = new TH1D(
        Form("h_layer_mean_gain_rms_frame_%d", layer),
        Form("Gain layer mean trend: layer %d;;Gain", layer),
        static_cast<int>(runs.size()), 0.5, runs.size() + 0.5);
    const int labelStride =
        runs.size() > 20 ? std::max(1, static_cast<int>(runs.size() / 18)) : 1;
    for (size_t irun = 0; irun < runs.size(); ++irun) {
      const bool showLabel = irun == 0 || irun + 1 == runs.size() ||
                             (static_cast<int>(irun) % labelStride == 0);
      frame->GetXaxis()->SetBinLabel(irun + 1,
                                     showLabel ? runs[irun].label.c_str() : "");
    }
    frame->GetXaxis()->LabelsOption("v");
    frame->GetXaxis()->SetLabelSize(0.032);
    frame->GetYaxis()->SetLabelSize(0.032);
    frame->GetYaxis()->SetTitleSize(0.040);
    frame->SetMinimum(ymin - 0.15 * span);
    frame->SetMaximum(ymax + 0.15 * span);
    frame->SetStats(0);
    frame->Draw();

    TGraphErrors *g =
        new TGraphErrors(x.size(), &x[0], &y[0], &ex[0], &ey[0]);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.0);
    g->SetMarkerColor(kBlue + 1);
    g->SetLineColor(kBlue + 1);
    g->SetLineWidth(2);
    g->Draw("PE1 SAME");

    TLatex *xTitle = new TLatex(runs.size() + 0.15, ymin - 0.08 * span, "Run");
    xTitle->SetTextSize(0.035);
    xTitle->SetTextAlign(12);
    xTitle->Draw();

    const std::string prefix =
        plotDir + Form("/layer%02d_gain_mean_rms", layer);
    c->SaveAs((prefix + ".pdf").c_str());
    c->SaveAs((prefix + ".png").c_str());
    delete c;
  }

  std::cout << "[compare] wrote layer mean gain RMS trends to " << plotDir
            << std::endl;
}

void drawGain2DMapsByRun(const std::vector<RunData> &runs,
                         const std::string &outDir) {
  if (runs.empty()) return;

  const std::set<int> selectedCells = selectCellsBySurvivalRate(runs, 0.50);
  if (selectedCells.empty()) {
    std::cout << "[compare] no channels with survival rate > 50% for gain 2D maps"
              << std::endl;
    return;
  }

  int maxX = -1;
  int maxChannel = -1;
  double maxGain = 0.0;
  int nCellsWithCoordinates = 0;
  for (const auto &run : runs) {
    for (const auto &entry : run.entries) {
      if (selectedCells.count(entry.cellid) == 0) continue;
      const int layer = entryLayer(entry);
      const int chip = entryChip(entry);
      const int channel = entryChannel(entry);
      if (layer < 0 || chip < 0 || channel < 0) continue;

      const int x = layer * 9 + chip;
      if (x < 0) continue;
      maxX = std::max(maxX, x);
      maxChannel = std::max(maxChannel, channel);
      if (isSurvivingSpePlotPoint(entry)) {
        maxGain = std::max(maxGain, entry.gain);
      }
      ++nCellsWithCoordinates;
    }
  }

  if (maxX < 0 || maxChannel < 0 || nCellsWithCoordinates == 0) {
    std::cout << "[compare] no valid layer/chip/channel coordinates for gain 2D maps"
              << std::endl;
    return;
  }

  const std::string mapDir = outDir + "/gain_2d_maps";
  gSystem->mkdir(mapDir.c_str(), true);

  TFile *rootOut = TFile::Open((mapDir + "/gain_2d_maps.root").c_str(), "RECREATE");
  const int nX = maxX + 1;
  const int nY = maxChannel + 1;
  const double zMax = maxGain > 0.0 ? maxGain * 1.05 : 1.0;

  for (const auto &run : runs) {
    const std::string san = sanitizeLabel(run.label);
    TH2D *h = new TH2D(
        Form("h_gain_2d_map_run_%s", san.c_str()),
        Form("Gain 2D map run %s;layer*9+chip;channel;gain", run.label.c_str()),
        nX, -0.5, maxX + 0.5, nY, -0.5, maxChannel + 0.5);
    h->SetDirectory(nullptr);
    h->SetStats(0);
    h->SetMinimum(0.0);
    h->SetMaximum(zMax);

    for (const auto &entry : run.entries) {
      if (selectedCells.count(entry.cellid) == 0) continue;
      const int layer = entryLayer(entry);
      const int chip = entryChip(entry);
      const int channel = entryChannel(entry);
      if (layer < 0 || chip < 0 || channel < 0) continue;

      const int x = layer * 9 + chip;
      if (x < 0) continue;
      double gain = std::isfinite(entry.gain) ? entry.gain : 0.0;
      if (!isSurvivingSpePlotPoint(entry)) gain = 0.0;
      h->SetBinContent(h->GetXaxis()->FindBin(x), h->GetYaxis()->FindBin(channel),
                       gain);
    }

    TCanvas *c = new TCanvas(Form("c_gain_2d_map_run_%s", san.c_str()),
                             "Gain 2D map", 1200, 800);
    c->SetRightMargin(0.16);
    c->SetLeftMargin(0.10);
    c->SetBottomMargin(0.12);
    c->SetGrid();
    h->Draw("COLZ");
    c->SaveAs((mapDir + "/gain_2d_map_run_" + san + ".pdf").c_str());
    c->SaveAs((mapDir + "/gain_2d_map_run_" + san + ".png").c_str());

    if (rootOut && !rootOut->IsZombie()) {
      rootOut->cd();
      h->Write();
    }
    delete c;
    delete h;
  }

  if (rootOut) {
    rootOut->Close();
    delete rootOut;
  }

  std::cout << "[compare] wrote gain 2D maps to " << mapDir << std::endl;
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
      if (!isAcceptedSpePoint(kv.second) || !isAcceptedSpePoint(it->second)) continue;
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
    if (!isAcceptedSpePoint(kv.second) || !isAcceptedSpePoint(it->second)) continue;
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
      if (e.fit_ok && !e.fit_mis_excluded) ++ok;
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
  const std::set<int> selectedCells = selectCellsBySurvivalRate(runs, 0.50);
  if (selectedCells.empty()) {
    std::cout << "[compare] no channels with survival rate > 50% for overlay plots"
              << std::endl;
    return;
  }

  struct PlotDef {
    std::string key;
    std::string title;
    std::string xTitle;
    bool onlyFitOk;
  };

  std::vector<PlotDef> defs = {
      {"gain", "SPE gain", "Gain (ADC/PE)", true},
      {"gainErr", "SPE gain error", "Gain error (ADC/PE)", true},
      {"snr", "SPE SNR", "SNR", true},
      {"noiseMed", "SPE noise median", "Noise (ADC)", true},
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
      std::vector<double> values;
      values.reserve(run.entries.size());
      for (const auto &e : run.entries) {
        if (selectedCells.count(e.cellid) == 0) continue;
        if (!isSurvivingSpePlotPoint(e)) continue;
        double v = 0.0;
        if (getBranchValue(e, def.key, def.onlyFitOk, v)) values.push_back(v);
      }
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

    const int bins = 50;

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
        if (selectedCells.count(e.cellid) == 0) continue;
        if (!isSurvivingSpePlotPoint(e)) continue;
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
      if (e.fit_ok && !e.fit_mis_excluded) ++ok;
    }

    const double frac = 100.0 * ok / std::max(1, static_cast<int>(run.entries.size()));

    std::cout << run.label << "\t" << run.entries.size() << "\t" << sGain.mean << "\t"
              << sGainErr.mean << "\t" << sSnr.mean << "\t" << sNoiseMed.mean << "\t"
              << sPeakAmp.mean << "\t" << sK0.mean << "\t" << frac << std::endl;
  }
  std::cout << std::endl;
}

bool testGainDriftDetectionSynthetic() {
  std::vector<RunData> runs(13);
  for (size_t irun = 0; irun < runs.size(); ++irun) {
    runs[irun].run = static_cast<int>(irun + 1);
    runs[irun].label = "synthetic_" + std::to_string(irun + 1);

    auto add = [&](int cellid, double gain, double snr, int fitOk = 1) {
      SpeEntry entry;
      entry.cellid = cellid;
      entry.layer = cellid / 100000;
      entry.chip = (cellid / 10000) % 10;
      entry.channel = cellid % 10000;
      entry.gain = gain;
      entry.snr = snr;
      entry.fit_ok = fitOk;
      runs[irun].entries.push_back(entry);
      runs[irun].byCellId[cellid] = entry;
    };

    add(100001, 10.0 * (1.0 + 0.005 * irun), 7.0);              // drift
    add(100002, irun < 6 ? 10.0 : 10.5, 7.0);                   // step
    add(100003, irun == 6 ? 11.0 : 10.0, 7.0);                  // one bad point
    add(100004, 10.0 + 0.01 * std::sin(irun), 7.0);             // stable
    add(100005, 10.0 * (1.0 + 0.005 * irun), irun < 7 ? 7 : 2); // too few
    add(100006, 10.0 * (1.0 + 0.005 * irun), 2.0);              // low SNR
    add(100007, 10.0 * (1.0 + 0.005 * irun), 7.0, 0);           // fit failed
  }

  const std::map<int, DriftCandidate> candidates =
      findDriftCandidates(runs, 5.0, 3.0, 8);
  const std::vector<StableChannelInfo> stable = rankStableChannels(runs, 5.0, 8);
  const bool ok = candidates.count(100001) == 1 &&
                  candidates.at(100001).type == "drift" &&
                  candidates.count(100002) == 1 &&
                  candidates.at(100002).type == "step" &&
                  candidates.count(100003) == 0 &&
                  candidates.count(100004) == 0 &&
                  candidates.count(100005) == 0 &&
                  candidates.count(100006) == 0 &&
                  candidates.count(100007) == 0 &&
                  !stable.empty() && stable.front().cellid == 100004;
  std::cout << "[compare:test] synthetic drift detection " << (ok ? "PASS" : "FAIL")
            << " (candidates=" << candidates.size() << ")" << std::endl;
  return ok;
}

void compare(const char *baseDir = ".",
             const char *runListCsv =
                 "21987-22163,22165-22296,22298-22560,22562-22797,22819-22856,22868-23181,23184-23563,23583-24064",
             int focusCellId = -1,
             const char *outDirName = "compare_plots_spe",
             const char *testBeamDir = "",
             double shiftedGainRelativeThreshold = 0.10,
             int shiftedGainMinValidRuns = 3,
             double driftMinChangePercent = 3.0,
             int driftMinValidRuns = 8,
             int driftTopPlots = 50) {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  // Accept numeric run lists (e.g. "22133,22135"), numeric ranges
  // (e.g. "21987-24064"), "all", or non-numeric directory-like tokens where each
  // token is a folder under `baseDir` containing a `spe_analysis.root` file.
  std::vector<std::string> tokens = parseTokenList(runListCsv);
  if (tokens.empty()) {
    std::cerr << "[compare] empty run list" << std::endl;
    return;
  }

  std::vector<std::string> expandedTokens;
  for (const auto &t : tokens) {
    if (isAllRunToken(t)) {
      const std::vector<std::string> allTokens = listAllRunTokens(baseDir);
      expandedTokens.insert(expandedTokens.end(), allTokens.begin(), allTokens.end());
    } else {
      long start = 0;
      long end = 0;
      const std::vector<std::string> rangeTokens =
          parseRunTokenRange(t, start, end) ? listRunTokensInRange(baseDir, start, end)
                                           : std::vector<std::string>();
      if (!rangeTokens.empty()) {
        expandedTokens.insert(expandedTokens.end(), rangeTokens.begin(), rangeTokens.end());
      } else {
        expandedTokens.push_back(t);
      }
    }
  }

  std::vector<std::string> uniqueTokens;
  std::set<std::string> seenTokens;
  for (const auto &t : expandedTokens) {
    if (seenTokens.insert(t).second) uniqueTokens.push_back(t);
  }
  tokens.swap(uniqueTokens);
  if (tokens.empty()) {
    std::cerr << "[compare] empty run list" << std::endl;
    return;
  }

  std::vector<RunData> runs;
  int pseudoRun = 1;
  for (const auto &t : tokens) {
    long runNumber = 0;
    if (parseIntegerToken(t, runNumber)) {
      RunData data;
      if (loadRun(baseDir, static_cast<int>(runNumber), data)) {
        std::cout << "[compare] loaded run " << data.label << " from " << data.path
                  << " with " << data.entries.size() << " cells" << std::endl;
        runs.push_back(data);
      }
      continue;
    }

    const std::string path = std::string(baseDir) + "/" + t + "/spe_refit_gaus_min11.root";
    if (gSystem->AccessPathName(path.c_str()) == 0) {
      RunData data;
      if (loadRunFromFile(path, pseudoRun, t, false, data)) {
        std::cout << "[compare] loaded run " << data.label << " from " << data.path
                  << " with " << data.entries.size() << " cells" << std::endl;
        runs.push_back(data);
        ++pseudoRun;
      } else {
        std::cout << "[compare] skipping missing/invalid token path " << path << std::endl;
      }
    } else {
      const std::vector<int> runNumbers = parseRunList(t.c_str());
      if (runNumbers.empty()) {
        std::cout << "[compare] skipping invalid run token " << t << std::endl;
        continue;
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
  }
  std::string outDir = std::string(baseDir) + "/" + outDirName;
  gSystem->mkdir(outDir.c_str(), true);

  if (runs.empty()) {
    std::cerr << "[compare] no runs loaded" << std::endl;
    return;
  }

  markHighGainFittingMis(runs, outDir + "/fitting_mis_gain_outliers.txt");
  writeAndDrawLayerEfficiencyByRun(runs, outDir);

  printSummaryTable(runs);
  writeShiftedGainReport(runs, outDir + "/shifted_gain_channels.txt",
                         shiftedGainRelativeThreshold, shiftedGainMinValidRuns);
  drawGainDeltaSnrStudies(runs, outDir, 10.0, shiftedGainMinValidRuns);
  runDriftAnalysis(runs, outDir, driftMinChangePercent, driftMinValidRuns,
                   driftTopPlots);
  drawStableChannelDistributions(runs, outDir, 5.0, driftMinValidRuns, 10);
  drawGainRunRmsDistribution(runs, outDir, 2);
  drawHighGainRmsChannelTrends(runs, outDir, driftTopPlots, 2);
  writeAndDrawLayerMeanGainShift(runs, outDir);
  writeAndDrawLayerMeanGainRmsTrends(runs, outDir);
  drawGain2DMapsByRun(runs, outDir);
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
