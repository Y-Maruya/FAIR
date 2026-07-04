#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace mip_imputation {

const int kChannelsPerLayer = 36 * 9;

struct Record {
  int runStart = -1;
  int runEnd = -1;
  std::string runLabel;
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;
  int entries = 0;
  int decision = -1;
  std::string decisionName;
  int state = 0;
  std::string stateNames;
  double mpv = -1.0;
  double width = -1.0;
  double gausSigma = -1.0;
  double directThreshold = -1.0;
  double directWidth = -1.0;
  double mpvError = -1.0;
  double widthError = -1.0;
  double gausSigmaError = -1.0;
  double directThresholdError = -1.0;
  double directWidthError = -1.0;
};

struct RunData {
  int runStart = -1;
  int runEnd = -1;
  std::string label;
  std::vector<Record> records;
  std::map<int, size_t> indexByCell;
  std::map<std::pair<int, int>, size_t> indexByLayerChannel;
};

struct ParameterInfo {
  const char *name;
  const char *jsonName;
  const char *jsonErrorName;
  const char *title;
  double Record::*valueMember;
  double Record::*errorMember;
};

const std::vector<ParameterInfo> &parameters() {
  static const std::vector<ParameterInfo> params = {
      {"mpv", "MPV", "MPVError", "MPV", &Record::mpv, &Record::mpvError},
      {"width", "Width", "WidthError", "Width", &Record::width, &Record::widthError},
      {"gaus_sigma", "GausSigma", "GausSigmaError", "Gaussian sigma",
       &Record::gausSigma, &Record::gausSigmaError},
      {"direct_threshold", "Threshold", "ThresholdError", "Threshold",
       &Record::directThreshold, &Record::directThresholdError},
      {"direct_width", "ThresholdWidth", "ThresholdWidthError", "Threshold width",
       &Record::directWidth, &Record::directWidthError}};
  return params;
}

struct Summary {
  size_t n = 0;
  double mean = 0.0;
  double rms = 0.0;
  double median = 0.0;
  double q16 = 0.0;
  double q84 = 0.0;
  double minimum = 0.0;
  double maximum = 0.0;
};

struct ChannelKey {
  int layer = -1;
  int index = -1;

  bool operator<(const ChannelKey &other) const {
    if (layer != other.layer) return layer < other.layer;
    return index < other.index;
  }
};

struct ImputeStats {
  int validRuns = 0;
  bool canImpute = false;
  double mean = -999.0;
  double rms = 0.0;
  double meanError = 0.0;
  double imputedError = -999.0;
};

struct StoredPoint {
  double value = -999.0;
  double error = -999.0;
  bool imputed = false;
  bool reference = false;
  bool fullFit = false;
};

typedef std::vector<std::map<ChannelKey, StoredPoint> > FullFitValueMap;

std::string directoryName(const std::string &path) {
  const size_t slash = path.find_last_of('/');
  if (slash == std::string::npos) return ".";
  if (slash == 0) return "/";
  return path.substr(0, slash);
}

std::string joinPath(const std::string &left, const std::string &right) {
  if (left.empty() || left == ".") return right;
  if (right.empty()) return left;
  if (right[0] == '/') return right;
  return left.back() == '/' ? left + right : left + "/" + right;
}

bool fileExists(const std::string &path) {
  return gSystem->AccessPathName(path.c_str()) == false;
}

std::string resolveFullFitPath(const std::string &inputPath,
                               const std::string &fullFitPath) {
  if (fullFitPath.empty()) return "";
  if (fullFitPath[0] == '/' || fileExists(fullFitPath)) return fullFitPath;

  const std::string inputDir = directoryName(inputPath);
  const std::string inputBaseDir = directoryName(inputDir);
  const std::string siblingPath = joinPath(inputBaseDir, fullFitPath);
  return fileExists(siblingPath) ? siblingPath : fullFitPath;
}

bool validValue(double value) {
  return std::isfinite(value) && value >= 0.0;
}

bool validError(double value) {
  return std::isfinite(value) && value >= 0.0;
}

int channelIndex(int chip, int channel) {
  if (chip < 0 || chip >= 9 || channel < 0 || channel >= 36) return -1;
  return chip * 36 + channel;
}

double quantile(const std::vector<double> &sorted, double probability) {
  if (sorted.empty()) return 0.0;
  const double position = probability * (sorted.size() - 1);
  const size_t low = static_cast<size_t>(std::floor(position));
  const size_t high = static_cast<size_t>(std::ceil(position));
  if (low == high) return sorted[low];
  return sorted[low] + (position - low) * (sorted[high] - sorted[low]);
}

Summary summarize(std::vector<double> values) {
  Summary result;
  result.n = values.size();
  if (values.empty()) return result;
  std::sort(values.begin(), values.end());
  double sum = 0.0;
  double sum2 = 0.0;
  for (double value : values) {
    sum += value;
    sum2 += value * value;
  }
  result.mean = sum / values.size();
  result.rms = std::sqrt(std::max(0.0, sum2 / values.size() - result.mean * result.mean));
  result.median = quantile(values, 0.50);
  result.q16 = quantile(values, 0.16);
  result.q84 = quantile(values, 0.84);
  result.minimum = values.front();
  result.maximum = values.back();
  return result;
}

std::pair<double, double> valueRange(const std::vector<double> &values) {
  const double low = *std::min_element(values.begin(), values.end());
  const double high = *std::max_element(values.begin(), values.end());
  const double padding = high > low ? 0.10 * (high - low)
                                    : std::max(1.0, 0.10 * std::abs(low));
  return std::make_pair(low - padding, high + padding);
}

int plotColor(size_t index) {
  static const int colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1,
                               kOrange + 7, kCyan + 2, kViolet + 1, kTeal + 3,
                               kPink + 7, kAzure + 7, kSpring + 5, kGray + 2};
  return colors[index % (sizeof(colors) / sizeof(colors[0]))];
}

std::string fileSafeString(const std::string &text) {
  std::string safe;
  for (char c : text) {
    safe += std::isalnum(static_cast<unsigned char>(c)) ? c : '_';
  }
  while (safe.find("__") != std::string::npos) safe.replace(safe.find("__"), 2, "_");
  if (!safe.empty() && safe.back() == '_') safe.pop_back();
  return safe.empty() ? "run" : safe;
}

std::string runKey(int start, int end, const std::string &label) {
  return std::to_string(start) + ":" + std::to_string(end) + ":" + label;
}

bool setBranchIfExists(TTree *tree, const char *name, void *address) {
  if (!tree->GetBranch(name)) return false;
  tree->SetBranchAddress(name, address);
  return true;
}

bool loadRuns(const char *inputPath, const char *treeName, std::vector<RunData> &runs) {
  TFile *file = TFile::Open(inputPath, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open " << inputPath << std::endl;
    delete file;
    return false;
  }

  TTree *tree = dynamic_cast<TTree *>(file->Get(treeName));
  if (!tree) {
    std::cerr << "Cannot find TTree '" << treeName << "' in " << inputPath << std::endl;
    file->Close();
    delete file;
    return false;
  }

  const char *required[] = {"run_start", "run_end", "cellid", "layer", "chip",
                            "channel", "entries", "decision", "state", "mpv",
                            "width", "gaus_sigma", "direct_threshold",
                            "direct_width"};
  for (const char *branch : required) {
    if (!tree->GetBranch(branch)) {
      std::cerr << "Missing required branch: " << branch << std::endl;
      file->Close();
      delete file;
      return false;
    }
  }

  Record value;
  std::string *runLabel = nullptr;
  std::string *decisionName = nullptr;
  std::string *stateNames = nullptr;

  tree->SetBranchAddress("run_start", &value.runStart);
  tree->SetBranchAddress("run_end", &value.runEnd);
  tree->SetBranchAddress("cellid", &value.cellid);
  tree->SetBranchAddress("layer", &value.layer);
  tree->SetBranchAddress("chip", &value.chip);
  tree->SetBranchAddress("channel", &value.channel);
  tree->SetBranchAddress("entries", &value.entries);
  tree->SetBranchAddress("decision", &value.decision);
  tree->SetBranchAddress("state", &value.state);
  tree->SetBranchAddress("mpv", &value.mpv);
  tree->SetBranchAddress("width", &value.width);
  tree->SetBranchAddress("gaus_sigma", &value.gausSigma);
  tree->SetBranchAddress("direct_threshold", &value.directThreshold);
  tree->SetBranchAddress("direct_width", &value.directWidth);
  setBranchIfExists(tree, "mpv_error", &value.mpvError);
  setBranchIfExists(tree, "width_error", &value.widthError);
  setBranchIfExists(tree, "gaus_sigma_error", &value.gausSigmaError);
  setBranchIfExists(tree, "direct_threshold_error", &value.directThresholdError);
  setBranchIfExists(tree, "direct_width_error", &value.directWidthError);
  if (tree->GetBranch("run_label")) tree->SetBranchAddress("run_label", &runLabel);
  if (tree->GetBranch("decision_name")) tree->SetBranchAddress("decision_name", &decisionName);
  if (tree->GetBranch("state_names")) tree->SetBranchAddress("state_names", &stateNames);

  std::map<std::string, size_t> runIndex;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    value = Record();
    tree->GetEntry(entry);
    value.runLabel = runLabel ? *runLabel
                              : std::to_string(value.runStart) + "-" + std::to_string(value.runEnd);
    value.decisionName = decisionName ? *decisionName : std::to_string(value.decision);
    value.stateNames = stateNames ? *stateNames : std::to_string(value.state);
    const std::string key = runKey(value.runStart, value.runEnd, value.runLabel);
    if (!runIndex.count(key)) {
      RunData run;
      run.runStart = value.runStart;
      run.runEnd = value.runEnd;
      run.label = value.runLabel;
      runIndex[key] = runs.size();
      runs.push_back(run);
    }
    runs[runIndex[key]].records.push_back(value);
  }

  file->Close();
  delete file;

  std::sort(runs.begin(), runs.end(), [](const RunData &a, const RunData &b) {
    if (a.runStart != b.runStart) return a.runStart < b.runStart;
    if (a.runEnd != b.runEnd) return a.runEnd < b.runEnd;
    return a.label < b.label;
  });

  for (RunData &run : runs) {
    for (size_t i = 0; i < run.records.size(); ++i) {
      const Record &record = run.records[i];
      run.indexByCell[record.cellid] = i;
      const int idx = channelIndex(record.chip, record.channel);
      if (record.layer >= 0 && idx >= 0) {
        run.indexByLayerChannel[std::make_pair(record.layer, idx)] = i;
      }
    }
  }
  return !runs.empty();
}

std::set<int> collectLayers(const std::vector<RunData> &runs) {
  std::set<int> layers;
  for (const RunData &run : runs) {
    for (const Record &record : run.records) {
      if (record.layer >= 0) layers.insert(record.layer);
    }
  }
  return layers;
}

std::set<ChannelKey> collectChannelKeys(const std::vector<RunData> &runs) {
  std::set<ChannelKey> keys;
  for (const RunData &run : runs) {
    for (const Record &record : run.records) {
      const int idx = channelIndex(record.chip, record.channel);
      if (record.layer < 0 || idx < 0) continue;
      ChannelKey key;
      key.layer = record.layer;
      key.index = idx;
      keys.insert(key);
    }
  }
  return keys;
}

bool isValidForParameter(const Record &record, const ParameterInfo &parameter, int minEntries) {
  return record.entries >= minEntries && validValue(record.*(parameter.valueMember));
}

bool useLayerReferenceFallback(size_t parameterIndex) {
  return parameterIndex < 3; // MPV, Width, GausSigma only.
}

bool containsStateName(const Record &record, const std::string &stateName) {
  return record.stateNames.find(stateName) != std::string::npos;
}

bool isThresholdStorageParameter(const ParameterInfo &parameter) {
  return parameter.valueMember == &Record::directThreshold ||
         parameter.valueMember == &Record::directWidth;
}

bool hasSuppressedThresholdState(const Record &record) {
  const int kEfficiencyRatioOver1p0 = 1 << 25;
  const int kEfficiencyRatioOver0p98 = 1 << 26;
  return (record.state & kEfficiencyRatioOver1p0) ||
         (record.state & kEfficiencyRatioOver0p98) ||
         containsStateName(record, "efficiency_ratio_over_1p0") ||
         containsStateName(record, "efficiency_ratio_0p98_to_1p0") ||
         containsStateName(record, "efficiency_ratio_over1p0") ||
         containsStateName(record, "efficiency_ratio_0p98_to1p0");
}

bool fullFitSucceededForParameter(const Record &record, size_t parameterIndex) {
  const int kFirstFitSuccessDecision = 1;
  const int kThresholdFitSuccessDecision = 2;
  const int kFirstFitSuccessState = 1 << 7;
  const int kSecondFitSuccessState = 1 << 16;
  const int kThresholdFitSuccessState = 1 << 24;

  const bool firstFitSuccess =
      record.decision == kFirstFitSuccessDecision ||
      record.decisionName == "first_fit_success" ||
      (record.state & kFirstFitSuccessState) ||
      containsStateName(record, "first_fit_success");
  const bool secondFitSuccess =
      (record.state & kSecondFitSuccessState) ||
      containsStateName(record, "second_fit_success");
  const bool thresholdFitSuccess =
      record.decision == kThresholdFitSuccessDecision ||
      record.decisionName == "threshold_fit_success" ||
      (record.state & kThresholdFitSuccessState) ||
      containsStateName(record, "threshold_fit_success");

  if (parameterIndex < 3) return firstFitSuccess || secondFitSuccess;
  return thresholdFitSuccess;
}

FullFitValueMap loadFullFitValues(const char *fullFitPath, const char *treeName) {
  FullFitValueMap values(parameters().size());
  if (!fullFitPath || !std::string(fullFitPath).size()) return values;

  std::vector<RunData> fullFitRuns;
  if (!loadRuns(fullFitPath, treeName, fullFitRuns)) {
    std::cerr << "Full-fit fallback disabled: could not load " << fullFitPath
              << std::endl;
    return values;
  }

  for (const RunData &run : fullFitRuns) {
    for (const Record &record : run.records) {
      const int idx = channelIndex(record.chip, record.channel);
      if (record.layer < 0 || idx < 0) continue;
      ChannelKey key;
      key.layer = record.layer;
      key.index = idx;

      for (size_t p = 0; p < parameters().size(); ++p) {
        const ParameterInfo &parameter = parameters()[p];
        if (!fullFitSucceededForParameter(record, p)) continue;
        const double value = record.*(parameter.valueMember);
        if (!validValue(value)) continue;

        StoredPoint point;
        point.value = value;
        const double error = record.*(parameter.errorMember);
        point.error = validError(error) ? error : -999.0;
        point.imputed = true;
        point.fullFit = true;
        values[p][key] = point;
      }
    }
  }

  std::cout << "Loaded full-fit fallback values from " << fullFitPath << ":";
  for (size_t p = 0; p < parameters().size(); ++p) {
    std::cout << " " << parameters()[p].jsonName << "=" << values[p].size();
  }
  std::cout << std::endl;
  return values;
}

std::map<int, Summary> layerSuccessSummaries(const RunData &run,
                                             const ParameterInfo &parameter,
                                             int minEntries) {
  std::map<int, std::vector<double> > valuesByLayer;
  for (const Record &record : run.records) {
    if (record.layer < 0) continue;
    if (!isValidForParameter(record, parameter, minEntries)) continue;
    valuesByLayer[record.layer].push_back(record.*(parameter.valueMember));
  }

  std::map<int, Summary> summaries;
  for (const auto &item : valuesByLayer) {
    summaries[item.first] = summarize(item.second);
  }
  return summaries;
}

std::vector<std::map<ChannelKey, ImputeStats> > computeImputeStats(
    const std::vector<RunData> &runs,
    const std::set<ChannelKey> &channelKeys,
    int minEntries,
    double minValidRunFraction) {
  const double clampedFraction = std::max(0.0, std::min(1.0, minValidRunFraction));
  const int requiredRuns =
      static_cast<int>(std::ceil(clampedFraction * static_cast<double>(runs.size())));
  std::vector<std::map<ChannelKey, ImputeStats> > stats(parameters().size());

  for (size_t p = 0; p < parameters().size(); ++p) {
    const ParameterInfo &parameter = parameters()[p];
    for (const ChannelKey &key : channelKeys) {
      std::vector<double> values;
      std::vector<double> errors;
      for (const RunData &run : runs) {
        const auto it = run.indexByLayerChannel.find(std::make_pair(key.layer, key.index));
        if (it == run.indexByLayerChannel.end()) continue;
        const Record &record = run.records[it->second];
        if (!isValidForParameter(record, parameter, minEntries)) continue;
        values.push_back(record.*(parameter.valueMember));
        const double error = record.*(parameter.errorMember);
        if (validError(error)) errors.push_back(error);
      }

      ImputeStats item;
      item.validRuns = values.size();
      if (!values.empty()) {
        const Summary valueSummary = summarize(values);
        item.mean = valueSummary.mean;
        item.rms = valueSummary.rms;
        if (!errors.empty()) {
          item.meanError = summarize(errors).mean;
        }
        item.imputedError =
            std::sqrt(item.rms * item.rms + item.meanError * item.meanError);
      }
      item.canImpute = item.validRuns >= requiredRuns && !values.empty();
      stats[p][key] = item;
    }
  }

  return stats;
}

std::vector<std::vector<std::map<ChannelKey, StoredPoint> > > buildStoredValues(
    const std::vector<RunData> &runs,
    const std::set<ChannelKey> &channelKeys,
    const std::vector<std::map<ChannelKey, ImputeStats> > &stats,
    const FullFitValueMap &fullFitValues,
    int minEntries,
    double minValidRunFraction) {
  std::vector<std::vector<std::map<ChannelKey, StoredPoint> > > stored(
      parameters().size(), std::vector<std::map<ChannelKey, StoredPoint> >(runs.size()));

  for (size_t p = 0; p < parameters().size(); ++p) {
    const ParameterInfo &parameter = parameters()[p];
    for (size_t r = 0; r < runs.size(); ++r) {
      const RunData &run = runs[r];
      const std::map<int, Summary> layerSummaries =
          layerSuccessSummaries(run, parameter, minEntries);
      for (const ChannelKey &key : channelKeys) {
        StoredPoint point;
        const auto recIt = run.indexByLayerChannel.find(std::make_pair(key.layer, key.index));
        if (recIt != run.indexByLayerChannel.end()) {
          const Record &record = run.records[recIt->second];
          if (isThresholdStorageParameter(parameter) &&
              hasSuppressedThresholdState(record)) {
            stored[p][r][key] = point;
            continue;
          }
          if (isValidForParameter(record, parameter, minEntries)) {
            point.value = record.*(parameter.valueMember);
            const double error = record.*(parameter.errorMember);
            point.error = validError(error) ? error : -999.0;
            stored[p][r][key] = point;
            continue;
          }
        }

        const auto statIt = stats[p].find(key);
        const double successFraction =
            runs.empty() || statIt == stats[p].end()
                ? 0.0
                : static_cast<double>(statIt->second.validRuns) / runs.size();
        bool filledByFullFit = false;
        if (p < fullFitValues.size()) {
          const auto fullFitIt = fullFitValues[p].find(key);
          if (fullFitIt != fullFitValues[p].end() &&
              validValue(fullFitIt->second.value)) {
            point = fullFitIt->second;
            const double fitError = point.error;
            const bool hasFitError = validError(fitError);
            const bool hasRunRms =
                statIt != stats[p].end() && statIt->second.validRuns > 0 &&
                std::isfinite(statIt->second.rms) && statIt->second.rms >= 0.0;
            if (hasFitError && hasRunRms) {
              const double runRms = statIt->second.rms;
              point.error = std::sqrt(fitError * fitError + runRms * runRms);
            } else if (hasRunRms) {
              point.error = statIt->second.rms;
            }
            point.imputed = true;
            point.fullFit = true;
            filledByFullFit = true;
          }
        }
        bool filledByChannelAverage = false;
        if (!filledByFullFit &&
            statIt != stats[p].end() && statIt->second.canImpute) {
          point.value = statIt->second.mean;
          point.error = statIt->second.imputedError;
          point.imputed = true;
          filledByChannelAverage = true;
        }
        if (!filledByFullFit && !filledByChannelAverage &&
            useLayerReferenceFallback(p) && successFraction <= minValidRunFraction) {
          const auto layerIt = layerSummaries.find(key.layer);
          if (layerIt != layerSummaries.end() && layerIt->second.n) {
            point.value = layerIt->second.mean;
            point.error = layerIt->second.rms;
            point.imputed = true;
            point.reference = true;
          }
        }
        stored[p][r][key] = point;
      }
    }
  }

  return stored;
}

json readJsonFile(const std::string &path) {
  std::ifstream input(path.c_str());
  if (!input.is_open()) return json();
  try {
    json content = json::parse(input);
    return content;
  } catch (const std::exception &error) {
    std::cerr << "Failed to parse JSON " << path << ": " << error.what() << std::endl;
  }
  return json();
}

std::string referenceJsonRunLabel(const std::string &runLabel, int runNumber) {
  if (runLabel == "23184-23563") {
    if (23184 <= runNumber && runNumber <= 23433) return "23184-23433";
    if (23436 <= runNumber && runNumber <= 23563) return "23436-23563";
  }
  return runLabel;
}

std::string referenceJsonDirectory(const std::string &baseDir,
                                   const std::string &runLabel,
                                   const std::string &jsonSubdir,
                                   int runNumber) {
  return joinPath(joinPath(baseDir, referenceJsonRunLabel(runLabel, runNumber)),
                  jsonSubdir);
}

std::string referenceJsonPath(const std::string &baseDir,
                              const std::string &runLabel,
                              const std::string &jsonSubdir,
                              int runNumber,
                              int layer) {
  return joinPath(referenceJsonDirectory(baseDir, runLabel, jsonSubdir, runNumber),
                  Form("run%d_mip_Layer%d.json", runNumber, layer));
}

bool parseLayerFromReferenceJsonName(const std::string &name, int &layer) {
  const std::string marker = "_mip_Layer";
  const size_t markerPos = name.rfind(marker);
  if (markerPos == std::string::npos) return false;

  const size_t digitStart = markerPos + marker.size();
  size_t digitEnd = digitStart;
  while (digitEnd < name.size() &&
         std::isdigit(static_cast<unsigned char>(name[digitEnd]))) {
    ++digitEnd;
  }
  if (digitEnd == digitStart || name.substr(digitEnd) != ".json") return false;

  layer = std::atoi(name.substr(digitStart, digitEnd - digitStart).c_str());
  return layer >= 0;
}

void collectLayersFromReferenceJsonDirectory(const std::string &directory,
                                             std::set<int> &layers) {
  void *dir = gSystem->OpenDirectory(directory.c_str());
  if (!dir) return;

  const char *entry = nullptr;
  while ((entry = gSystem->GetDirEntry(dir))) {
    int layer = -1;
    if (parseLayerFromReferenceJsonName(entry, layer)) layers.insert(layer);
  }
  gSystem->FreeDirectory(dir);
}

std::set<int> collectOutputLayers(const std::vector<RunData> &runs,
                                  const std::string &jsonBaseDir,
                                  const std::string &jsonSubdir) {
  std::set<int> layers = collectLayers(runs);
  std::set<std::string> referenceDirs;
  for (const RunData &run : runs) {
    referenceDirs.insert(referenceJsonDirectory(jsonBaseDir, run.label, jsonSubdir,
                                                run.runStart));
    if (run.label == "23184-23563") {
      referenceDirs.insert(referenceJsonDirectory(jsonBaseDir, run.label, jsonSubdir,
                                                  23184));
      referenceDirs.insert(referenceJsonDirectory(jsonBaseDir, run.label, jsonSubdir,
                                                  23436));
    }
  }

  for (const std::string &directory : referenceDirs) {
    collectLayersFromReferenceJsonDirectory(directory, layers);
  }
  return layers;
}

std::vector<int> sameDataRunsFromJson(const json &content, int fallbackRun) {
  if (content.is_object() && content.contains("Summary") &&
      content["Summary"].is_object() &&
      content["Summary"].contains("SameDataRuns") &&
      content["Summary"]["SameDataRuns"].is_array()) {
    try {
      return content["Summary"]["SameDataRuns"].get<std::vector<int> >();
    } catch (const std::exception &) {
    }
  }
  return std::vector<int>(1, fallbackRun);
}

std::vector<int> mergedSameDataRunsForReference(const std::string &baseDir,
                                                const std::string &runLabel,
                                                const std::string &jsonSubdir,
                                                int runStart,
                                                int layer,
                                                const json &primaryJson) {
  std::vector<int> runs = sameDataRunsFromJson(primaryJson, runStart);
  if (runLabel != "23184-23563") return runs;

  const std::string secondPath =
      referenceJsonPath(baseDir, runLabel, jsonSubdir, 23436, layer);
  const json secondJson = readJsonFile(secondPath);
  if (!secondJson.is_object()) {
    std::cerr << "Reference JSON not found or invalid for split run group: "
              << secondPath << std::endl;
    return runs;
  }

  const std::vector<int> secondRuns = sameDataRunsFromJson(secondJson, 23436);
  runs.insert(runs.end(), secondRuns.begin(), secondRuns.end());
  std::sort(runs.begin(), runs.end());
  runs.erase(std::unique(runs.begin(), runs.end()), runs.end());
  return runs;
}

template <typename T>
std::vector<T> sizedArray(int size, const T &value) {
  return std::vector<T>(size, value);
}

int imputedBit(size_t parameterIndex) {
  return 1 << static_cast<int>(parameterIndex);
}

json imputedBitDefinitions() {
  json definitions = json::object();
  for (size_t p = 0; p < parameters().size(); ++p) {
    definitions[parameters()[p].jsonName] = imputedBit(p);
  }
  return definitions;
}

void pruneUnwantedPerChannelKeys(json &output) {
  if (!output.contains("PerChannel") || !output["PerChannel"].is_object()) return;
  static const char *removeKeys[] = {
      "Chi2", "chi2", "EntriesADCLe50", "FWHM", "FitStatus", "MaxX",
      "NDF", "RatioADCLe50", "StateName", "StateNames", "TotalArea",
      "MPVImputed", "WidthImputed", "GausSigmaImputed", "ThresholdImputed",
      "ThresholdWidthImputed"};
  for (const char *key : removeKeys) {
    output["PerChannel"].erase(key);
  }
}

void fillDecisionArrays(const RunData &run,
                        int layer,
                        std::vector<int> &decision,
                        std::vector<int> &state,
                        std::vector<std::string> &decisionName) {
  for (int idx = 0; idx < kChannelsPerLayer; ++idx) {
    const auto it = run.indexByLayerChannel.find(std::make_pair(layer, idx));
    if (it == run.indexByLayerChannel.end()) continue;
    const Record &record = run.records[it->second];
    decision[idx] = record.decision;
    state[idx] = record.state;
    decisionName[idx] = record.decisionName;
  }
}

void writeJsonOutputs(
    const std::vector<RunData> &runs,
    const std::set<int> &layers,
    const std::vector<std::vector<std::map<ChannelKey, StoredPoint> > > &stored,
    const std::string &outDir,
    const std::string &jsonBaseDir,
    const std::string &jsonSubdir,
    int minEntries,
    double minValidRunFraction) {
  gSystem->mkdir(outDir.c_str(), true);

  for (size_t r = 0; r < runs.size(); ++r) {
    const RunData &run = runs[r];
    for (int layer : layers) {
      std::vector<std::vector<double> > values(parameters().size());
      std::vector<std::vector<double> > errors(parameters().size());
      std::vector<int> imputedMask = sizedArray(kChannelsPerLayer, 0);
      std::vector<int> refMask = sizedArray(kChannelsPerLayer, 0);
      std::vector<int> fullFitMask = sizedArray(kChannelsPerLayer, 0);
      std::vector<int> imputedCounts(parameters().size(), 0);
      std::vector<int> refCounts(parameters().size(), 0);
      std::vector<int> fullFitCounts(parameters().size(), 0);
      std::vector<int> unfilledCounts(parameters().size(), 0);

      for (size_t p = 0; p < parameters().size(); ++p) {
        values[p] = sizedArray(kChannelsPerLayer, -999.0);
        errors[p] = sizedArray(kChannelsPerLayer, -999.0);
        for (int idx = 0; idx < kChannelsPerLayer; ++idx) {
          ChannelKey key;
          key.layer = layer;
          key.index = idx;
          const auto it = stored[p][r].find(key);
          if (it == stored[p][r].end()) {
            ++unfilledCounts[p];
            continue;
          }
          values[p][idx] = it->second.value;
          errors[p][idx] = it->second.error;
          if (it->second.imputed) {
            imputedMask[idx] |= imputedBit(p);
            ++imputedCounts[p];
          }
          if (it->second.reference) {
            refMask[idx] |= imputedBit(p);
            ++refCounts[p];
          }
          if (it->second.fullFit) {
            fullFitMask[idx] |= imputedBit(p);
            ++fullFitCounts[p];
          }
          if (!validValue(it->second.value)) ++unfilledCounts[p];
        }
      }

      std::vector<int> decision = sizedArray(kChannelsPerLayer, -1);
      std::vector<int> state = sizedArray(kChannelsPerLayer, 0);
      std::vector<std::string> decisionName(kChannelsPerLayer, "");
      fillDecisionArrays(run, layer, decision, state, decisionName);

      const std::string refPath =
          referenceJsonPath(jsonBaseDir, run.label, jsonSubdir, run.runStart, layer);
      json refJson = readJsonFile(refPath);
      if (!refJson.is_object()) {
        std::cerr << "Reference JSON not found or invalid: " << refPath
                  << " (metadata will be limited)" << std::endl;
      }

      const std::vector<int> sameDataRuns = mergedSameDataRunsForReference(
          jsonBaseDir, run.label, jsonSubdir, run.runStart, layer, refJson);
      for (int sameRun : sameDataRuns) {
        const std::string samePath =
            referenceJsonPath(jsonBaseDir, run.label, jsonSubdir, sameRun, layer);
        const std::string sameReferenceRunLabel = referenceJsonRunLabel(run.label, sameRun);
        json sameJson = readJsonFile(samePath);
        json output = sameJson.is_object() ? sameJson
                                           : (refJson.is_object() ? refJson : json::object());
        pruneUnwantedPerChannelKeys(output);
        if (!sameJson.is_object() && fileExists(samePath)) {
          std::cerr << "Using fallback metadata for " << samePath << std::endl;
        }

        output["RunNumber"] = sameRun;
        output["Layer"] = layer;
        output["CalibrationType"] = "MIPImputed";
        if (!output.contains("Status")) output["Status"] = 0;
        if (!output.contains("Summary") || !output["Summary"].is_object()) {
          output["Summary"] = json::object();
        }
        output["Summary"]["MIPImputationMinEntries"] = minEntries;
        output["Summary"]["MIPImputationMinValidRunFraction"] = minValidRunFraction;
        output["Summary"]["MIPImputationReferenceJsonSubdir"] = jsonSubdir;
        output["Summary"]["MIPImputationReferenceRunLabel"] = run.label;
        output["Summary"]["MIPImputationReferenceJsonRunLabel"] = sameReferenceRunLabel;
        output["Summary"]["SameDataRuns"] = sameDataRuns;
        output["Summary"]["NumUsedRun"] = static_cast<int>(sameDataRuns.size());
        output["Summary"]["ImputedBitDefinitions"] = imputedBitDefinitions();
        output["Summary"]["ImputedAnyChannels"] =
            std::count_if(imputedMask.begin(), imputedMask.end(),
                          [](int value) { return value != 0; });
        output["Summary"]["RefBitDefinitions"] = imputedBitDefinitions();
        output["Summary"]["RefAnyChannels"] =
            std::count_if(refMask.begin(), refMask.end(),
                          [](int value) { return value != 0; });
        output["Summary"]["FullFitBitDefinitions"] = imputedBitDefinitions();
        output["Summary"]["FullFitAnyChannels"] =
            std::count_if(fullFitMask.begin(), fullFitMask.end(),
                          [](int value) { return value != 0; });

        if (!output.contains("PerChannel") || !output["PerChannel"].is_object()) {
          output["PerChannel"] = json::object();
        }

        for (size_t p = 0; p < parameters().size(); ++p) {
          const ParameterInfo &parameter = parameters()[p];
          output["PerChannel"][parameter.jsonName] = values[p];
          output["PerChannel"][parameter.jsonErrorName] = errors[p];
          output["Summary"][std::string(parameter.jsonName) + "ImputedChannels"] =
              imputedCounts[p];
          output["Summary"][std::string(parameter.jsonName) + "RefChannels"] =
              refCounts[p];
          output["Summary"][std::string(parameter.jsonName) + "FullFitChannels"] =
              fullFitCounts[p];
          output["Summary"][std::string(parameter.jsonName) + "UnfilledChannels"] =
              unfilledCounts[p];
        }

        output["PerChannel"]["Imputed"] = imputedMask;
        output["PerChannel"]["Ref"] = refMask;
        output["PerChannel"]["FullFit"] = fullFitMask;
        output["PerChannel"]["Decision"] = decision;
        output["PerChannel"]["State"] = state;
        output["PerChannel"]["DecisionName"] = decisionName;

        const std::string outPath =
            joinPath(outDir, Form("run%d_mip_imputed_Layer%d.json", sameRun, layer));
        std::ofstream out(outPath.c_str());
        if (!out.is_open()) {
          std::cerr << "Cannot write JSON: " << outPath << std::endl;
          continue;
        }
        out << output.dump(2) << std::endl;
        std::cout << "Wrote " << outPath << std::endl;
      }
    }
  }
}

std::vector<double> collectStoredLayerValues(
    const std::map<ChannelKey, StoredPoint> &storedRun,
    int layer) {
  std::vector<double> values;
  for (int idx = 0; idx < kChannelsPerLayer; ++idx) {
    ChannelKey key;
    key.layer = layer;
    key.index = idx;
    const auto it = storedRun.find(key);
    if (it == storedRun.end()) continue;
    if (validValue(it->second.value)) values.push_back(it->second.value);
  }
  return values;
}

int storageModeCode(const StoredPoint &point) {
  if (!validValue(point.value)) return 0;
  if (point.fullFit) return 3;
  if (point.reference) return 4;
  if (point.imputed) return 2;
  return 1;
}

void writeStorageModeLegend(const std::string &mapDir) {
  std::ofstream csv(joinPath(mapDir, "storage_mode_codes.csv").c_str());
  csv << "code,mode\n";
  csv << "0,missing\n";
  csv << "1,measured\n";
  csv << "2,channel_average_imputed\n";
  csv << "3,full_fit_fallback\n";
  csv << "4,layer_reference\n";
}

void drawStorageModeMaps(
    const std::vector<RunData> &runs,
    const std::set<int> &layers,
    const std::vector<std::vector<std::map<ChannelKey, StoredPoint> > > &stored,
    const std::string &plotDir,
    TFile &rootFile) {
  if (runs.empty() || layers.empty()) return;

  const std::string mapDir = joinPath(plotDir, "storage_mode_maps");
  gSystem->mkdir(mapDir.c_str(), true);
  writeStorageModeLegend(mapDir);

  const int minX = (*layers.begin()) * 9;
  const int maxX = (*layers.rbegin()) * 9 + 8;
  const int nXBins = std::max(1, maxX - minX + 1);

  for (size_t p = 0; p < parameters().size(); ++p) {
    const ParameterInfo &parameter = parameters()[p];
    const std::string parameterDir = joinPath(mapDir, parameter.jsonName);
    gSystem->mkdir(parameterDir.c_str(), true);

    for (size_t r = 0; r < runs.size(); ++r) {
      const std::string runTag = fileSafeString(runs[r].label);
      TCanvas canvas(Form("c_storage_mode_%s_%s",
                          parameter.jsonName, runTag.c_str()),
                     "", 1600, 850);
      canvas.SetLeftMargin(0.08);
      canvas.SetRightMargin(0.18);
      canvas.SetBottomMargin(0.12);
      canvas.SetGrid();

      TH2D hist(Form("h_storage_mode_%s_%s",
                     parameter.jsonName, runTag.c_str()),
                Form("%s storage mode: %s;layer*9 + chip;channel;storage mode",
                     parameter.title, runs[r].label.c_str()),
                nXBins, minX - 0.5, maxX + 0.5,
                36, -0.5, 35.5);
      hist.SetStats(false);
      hist.SetMinimum(-0.5);
      hist.SetMaximum(4.5);
      hist.SetContour(5);
      hist.GetZaxis()->SetTitle(
          "0 missing, 1 measured, 2 avg, 3 full-fit, 4 layer-ref");
      hist.GetZaxis()->SetTitleOffset(1.25);

      for (int layer : layers) {
        for (int chip = 0; chip < 9; ++chip) {
          const int layerChip = layer * 9 + chip;
          const int xBin = layerChip - minX + 1;
          for (int channel = 0; channel < 36; ++channel) {
            ChannelKey key;
            key.layer = layer;
            key.index = chip * 36 + channel;
            int mode = 0;
            const auto it = stored[p][r].find(key);
            if (it != stored[p][r].end()) mode = storageModeCode(it->second);
            hist.SetBinContent(xBin, channel + 1, mode);
          }
        }
      }

      hist.Draw("COLZ");

      const std::string baseName =
          joinPath(parameterDir,
                   Form("%s_%s_storage_mode_map",
                        parameter.jsonName, runTag.c_str()));
      canvas.SaveAs((baseName + ".pdf").c_str());
      canvas.SaveAs((baseName + ".png").c_str());

      rootFile.cd();
      canvas.Write();
      hist.Write();
    }
  }
}

void drawLayerMedianTrends(
    const std::vector<RunData> &runs,
    const std::set<int> &layers,
    const std::vector<std::vector<std::map<ChannelKey, StoredPoint> > > &stored,
    const std::string &plotDir,
    TFile &rootFile) {
  for (size_t p = 0; p < parameters().size(); ++p) {
    const ParameterInfo &parameter = parameters()[p];
    const std::string parameterDir = joinPath(plotDir, parameter.jsonName);
    gSystem->mkdir(parameterDir.c_str(), true);

    for (int layer : layers) {
      std::vector<double> x;
      std::vector<double> y;
      std::vector<double> exLow;
      std::vector<double> exHigh;
      std::vector<double> eyLow;
      std::vector<double> eyHigh;
      for (size_t r = 0; r < runs.size(); ++r) {
        const Summary summary = summarize(collectStoredLayerValues(stored[p][r], layer));
        if (!summary.n) continue;
        x.push_back(r + 1.0);
        y.push_back(summary.median);
        exLow.push_back(0.0);
        exHigh.push_back(0.0);
        eyLow.push_back(summary.median - summary.q16);
        eyHigh.push_back(summary.q84 - summary.median);
      }
      if (x.empty()) continue;

      double low = y.front() - eyLow.front();
      double high = y.front() + eyHigh.front();
      for (size_t i = 1; i < y.size(); ++i) {
        low = std::min(low, y[i] - eyLow[i]);
        high = std::max(high, y[i] + eyHigh[i]);
      }
      const double padding = high > low ? 0.10 * (high - low)
                                        : std::max(1.0, 0.10 * std::abs(low));

      TCanvas canvas(Form("c_stored_%s_layer%02d_median_trend",
                          parameter.jsonName, layer),
                     "", 1200, 750);
      canvas.SetBottomMargin(0.22);
      canvas.SetLeftMargin(0.10);
      canvas.SetGrid();
      TH1D frame(Form("frame_stored_%s_layer%02d_median_trend",
                      parameter.jsonName, layer),
                 Form("Stored %s median trend: layer %02d;Run;%s",
                      parameter.title, layer, parameter.jsonName),
                 runs.size(), 0.5, runs.size() + 0.5);
      for (size_t r = 0; r < runs.size(); ++r) {
        frame.GetXaxis()->SetBinLabel(r + 1, runs[r].label.c_str());
      }
      frame.GetXaxis()->LabelsOption("v");
      frame.SetMinimum(low - padding);
      frame.SetMaximum(high + padding);
      frame.SetStats(false);
      frame.Draw();

      TGraphAsymmErrors graph(x.size(), &x[0], &y[0],
                              &exLow[0], &exHigh[0], &eyLow[0], &eyHigh[0]);
      graph.SetName(Form("g_stored_%s_layer%02d_median_trend",
                         parameter.jsonName, layer));
      graph.SetMarkerStyle(20);
      graph.SetLineWidth(2);
      graph.SetMarkerColor(kBlue + 1);
      graph.SetLineColor(kBlue + 1);
      graph.Draw("P SAME");

      const std::string baseName =
          joinPath(parameterDir, Form("%s_layer%02d_median_trend",
                                      parameter.jsonName, layer));
      canvas.SaveAs((baseName + ".pdf").c_str());
      canvas.SaveAs((baseName + ".png").c_str());
      rootFile.cd();
      canvas.Write();
      frame.Write();
      graph.Write();
    }
  }
}

std::set<ChannelKey> imputedChannelKeys(
    const std::vector<std::map<ChannelKey, StoredPoint> > &storedParam,
    bool referenceOnly = false) {
  std::set<ChannelKey> result;
  for (const auto &runMap : storedParam) {
    for (const auto &item : runMap) {
      if (referenceOnly) {
        if (item.second.reference) result.insert(item.first);
      } else if (item.second.imputed) {
        result.insert(item.first);
      }
    }
  }
  return result;
}

std::string channelKeyName(const ChannelKey &key) {
  const int chip = key.index / 36;
  const int channel = key.index % 36;
  return Form("layer%02d_chip%02d_channel%02d", key.layer, chip, channel);
}

void drawImputedChannelTrends(
    const std::vector<RunData> &runs,
    const std::vector<std::vector<std::map<ChannelKey, StoredPoint> > > &stored,
    const std::string &plotDir,
    TFile &rootFile,
    const std::string &channelSubdir = "imputed_channels",
    bool referenceOnly = false) {
  const std::string channelDir = joinPath(plotDir, channelSubdir);
  gSystem->mkdir(channelDir.c_str(), true);
  const std::string objectTag = referenceOnly ? "reference_channel" : "stored";

  for (size_t p = 0; p < parameters().size(); ++p) {
    const ParameterInfo &parameter = parameters()[p];
    const std::set<ChannelKey> keys = imputedChannelKeys(stored[p], referenceOnly);
    if (keys.empty()) continue;

    const std::string parameterDir = joinPath(channelDir, parameter.jsonName);
    gSystem->mkdir(parameterDir.c_str(), true);

    for (const ChannelKey &key : keys) {
      std::vector<double> realX;
      std::vector<double> realY;
      std::vector<double> realEX;
      std::vector<double> realEY;
      std::vector<double> imputedX;
      std::vector<double> imputedY;
      std::vector<double> imputedEX;
      std::vector<double> imputedEY;
      std::vector<double> fullFitX;
      std::vector<double> fullFitY;
      std::vector<double> fullFitEX;
      std::vector<double> fullFitEY;
      std::vector<double> referenceX;
      std::vector<double> referenceY;
      std::vector<double> referenceEX;
      std::vector<double> referenceEY;
      std::vector<double> missingX;
      std::vector<double> rangeValues;

      for (size_t r = 0; r < runs.size(); ++r) {
        const auto it = stored[p][r].find(key);
        if (it == stored[p][r].end() || !validValue(it->second.value)) {
          missingX.push_back(r + 1.0);
          continue;
        }
        rangeValues.push_back(it->second.value);
        if (validError(it->second.error)) {
          rangeValues.push_back(it->second.value - it->second.error);
          rangeValues.push_back(it->second.value + it->second.error);
        }
        if (it->second.fullFit) {
          fullFitX.push_back(r + 1.0);
          fullFitY.push_back(it->second.value);
          fullFitEX.push_back(0.0);
          fullFitEY.push_back(validError(it->second.error) ? it->second.error : 0.0);
        } else if (it->second.reference) {
          referenceX.push_back(r + 1.0);
          referenceY.push_back(it->second.value);
          referenceEX.push_back(0.0);
          referenceEY.push_back(validError(it->second.error) ? it->second.error : 0.0);
        } else if (it->second.imputed) {
          imputedX.push_back(r + 1.0);
          imputedY.push_back(it->second.value);
          imputedEX.push_back(0.0);
          imputedEY.push_back(validError(it->second.error) ? it->second.error : 0.0);
        } else {
          realX.push_back(r + 1.0);
          realY.push_back(it->second.value);
          realEX.push_back(0.0);
          realEY.push_back(validError(it->second.error) ? it->second.error : 0.0);
        }
      }
      if (rangeValues.empty()) continue;

      const std::pair<double, double> range = valueRange(rangeValues);
      const double span = std::max(1.0, range.second - range.first);
      const double missingY = range.first - 0.08 * span;
      const double frameYMin = range.first - 0.18 * span;
      const double frameYMax = range.second + 0.12 * span;
      std::vector<double> missingYValues(missingX.size(), missingY);

      const std::string keyName = channelKeyName(key);
      TCanvas canvas(Form("c_%s_%s_%s_trend",
                          objectTag.c_str(), parameter.jsonName, keyName.c_str()),
                     "", 1200, 750);
      canvas.SetBottomMargin(0.22);
      canvas.SetLeftMargin(0.10);
      canvas.SetGrid();

      TH1D frame(Form("frame_%s_%s_%s_trend",
                      objectTag.c_str(), parameter.jsonName, keyName.c_str()),
                 Form("%s %s trend: %s;Run;%s",
                      referenceOnly ? "Layer-reference" : "Stored",
                      parameter.title, keyName.c_str(), parameter.jsonName),
                 runs.size(), 0.5, runs.size() + 0.5);
      for (size_t r = 0; r < runs.size(); ++r) {
        frame.GetXaxis()->SetBinLabel(r + 1, runs[r].label.c_str());
      }
      frame.GetXaxis()->LabelsOption("v");
      frame.SetMinimum(frameYMin);
      frame.SetMaximum(frameYMax);
      frame.SetStats(false);
      frame.Draw();

      TGraphErrors realGraph;
      if (!realX.empty()) {
        realGraph = TGraphErrors(realX.size(), &realX[0], &realY[0], &realEX[0], &realEY[0]);
        realGraph.SetName(Form("g_%s_%s_%s_real",
                               objectTag.c_str(), parameter.jsonName, keyName.c_str()));
        realGraph.SetMarkerStyle(20);
        realGraph.SetMarkerColor(kBlue + 1);
        realGraph.SetLineColor(kBlue + 1);
        realGraph.SetLineWidth(2);
        realGraph.Draw("PL SAME");
      }

      TGraphErrors imputedGraph;
      if (!imputedX.empty()) {
        imputedGraph = TGraphErrors(imputedX.size(), &imputedX[0], &imputedY[0],
                                    &imputedEX[0], &imputedEY[0]);
        imputedGraph.SetName(Form("g_%s_%s_%s_imputed",
                                  objectTag.c_str(), parameter.jsonName, keyName.c_str()));
        imputedGraph.SetMarkerStyle(21);
        imputedGraph.SetMarkerColor(kRed + 1);
        imputedGraph.SetLineColor(kRed + 1);
        imputedGraph.SetLineWidth(2);
        imputedGraph.Draw("P SAME");
      }

      TGraphErrors fullFitGraph;
      if (!fullFitX.empty()) {
        fullFitGraph = TGraphErrors(fullFitX.size(), &fullFitX[0], &fullFitY[0],
                                    &fullFitEX[0], &fullFitEY[0]);
        fullFitGraph.SetName(Form("g_%s_%s_%s_full_fit",
                                  objectTag.c_str(), parameter.jsonName, keyName.c_str()));
        fullFitGraph.SetMarkerStyle(33);
        fullFitGraph.SetMarkerColor(kOrange + 7);
        fullFitGraph.SetLineColor(kOrange + 7);
        fullFitGraph.SetLineWidth(2);
        fullFitGraph.Draw("P SAME");
      }

      TGraphErrors referenceGraph;
      if (!referenceX.empty()) {
        referenceGraph = TGraphErrors(referenceX.size(), &referenceX[0], &referenceY[0],
                                      &referenceEX[0], &referenceEY[0]);
        referenceGraph.SetName(Form("g_%s_%s_%s_reference",
                                    objectTag.c_str(), parameter.jsonName, keyName.c_str()));
        referenceGraph.SetMarkerStyle(22);
        referenceGraph.SetMarkerColor(kGreen + 2);
        referenceGraph.SetLineColor(kGreen + 2);
        referenceGraph.SetLineWidth(2);
        referenceGraph.Draw("P SAME");
      }

      TGraph missingGraph;
      if (!missingX.empty()) {
        missingGraph.SetName(Form("g_%s_%s_%s_missing",
                                  objectTag.c_str(), parameter.jsonName, keyName.c_str()));
        missingGraph.SetMarkerStyle(5);
        missingGraph.SetMarkerSize(1.3);
        missingGraph.SetMarkerColor(kGray + 2);
        missingGraph.SetLineColor(kGray + 2);
        for (size_t i = 0; i < missingX.size(); ++i) {
          missingGraph.SetPoint(i, missingX[i], missingYValues[i]);
        }
        missingGraph.Draw("P SAME");
      }

      TLegend legend(0.56, 0.68, 0.94, 0.89);
      legend.SetBorderSize(0);
      if (!realX.empty()) legend.AddEntry(&realGraph, "measured", "pl");
      if (!imputedX.empty()) legend.AddEntry(&imputedGraph, "channel-average imputed", "p");
      if (!fullFitX.empty()) legend.AddEntry(&fullFitGraph, "full-fit fallback", "p");
      if (!referenceX.empty()) legend.AddEntry(&referenceGraph, "layer-reference", "p");
      if (!missingX.empty()) legend.AddEntry(&missingGraph, "unfilled", "p");
      legend.Draw();

      const std::string baseName =
          joinPath(parameterDir, std::string(parameter.jsonName) + "_" + keyName + "_trend");
      canvas.SaveAs((baseName + ".pdf").c_str());
      canvas.SaveAs((baseName + ".png").c_str());
      rootFile.cd();
      canvas.Write();
      frame.Write();
      if (!realX.empty()) realGraph.Write();
      if (!imputedX.empty()) imputedGraph.Write();
      if (!fullFitX.empty()) fullFitGraph.Write();
      if (!referenceX.empty()) referenceGraph.Write();
      if (!missingX.empty()) missingGraph.Write();
    }
  }
}

size_t parameterIndexByJsonName(const std::string &jsonName) {
  for (size_t p = 0; p < parameters().size(); ++p) {
    if (jsonName == parameters()[p].jsonName) return p;
  }
  return parameters().size();
}

void drawMpvFillFractionByLayer(
    const std::vector<RunData> &runs,
    const std::set<int> &layers,
    const std::vector<std::vector<std::map<ChannelKey, StoredPoint> > > &stored,
    const std::string &plotDir,
    TFile &rootFile) {
  const size_t mpvIndex = parameterIndexByJsonName("MPV");
  if (mpvIndex >= parameters().size() || layers.empty()) return;

  const std::string fractionDir = joinPath(plotDir, "MPVFillFraction");
  gSystem->mkdir(fractionDir.c_str(), true);

  std::ofstream csv(joinPath(fractionDir, "mpv_fill_fraction_by_layer.csv").c_str());
  csv << "run_start,run_end,run_label,layer,filled_channels,total_channels,filled_fraction\n";

  TCanvas canvas("c_mpv_fill_fraction_by_layer", "", 1200, 750);
  canvas.SetBottomMargin(0.12);
  canvas.SetLeftMargin(0.10);
  canvas.SetGrid();

  const int minLayer = *layers.begin();
  const int maxLayer = *layers.rbegin();
  TH1D frame("frame_mpv_fill_fraction_by_layer",
             "Final MPV filled channel fraction by layer;Layer;Fraction of channels with MPV",
             std::max(1, maxLayer - minLayer + 1),
             minLayer - 0.5, maxLayer + 0.5);
  frame.SetMinimum(0.0);
  frame.SetMaximum(1.05);
  frame.SetStats(false);
  frame.Draw();

  TLegend legend(0.55, 0.18, 0.92, 0.45);
  legend.SetBorderSize(0);
  if (runs.size() > 6) legend.SetNColumns(2);

  std::vector<TGraph *> graphs;
  for (size_t r = 0; r < runs.size(); ++r) {
    std::vector<double> x;
    std::vector<double> y;
    for (int layer : layers) {
      int filled = 0;
      for (int idx = 0; idx < kChannelsPerLayer; ++idx) {
        ChannelKey key;
        key.layer = layer;
        key.index = idx;
        const auto it = stored[mpvIndex][r].find(key);
        if (it != stored[mpvIndex][r].end() && validValue(it->second.value)) {
          ++filled;
        }
      }
      const double fraction = static_cast<double>(filled) / kChannelsPerLayer;
      x.push_back(layer);
      y.push_back(fraction);
      csv << runs[r].runStart << "," << runs[r].runEnd << ","
          << runs[r].label << "," << layer << "," << filled << ","
          << kChannelsPerLayer << "," << std::setprecision(12)
          << fraction << "\n";
    }
    if (x.empty()) continue;

    TGraph *graph = new TGraph(x.size(), &x[0], &y[0]);
    graph->SetName(Form("g_mpv_fill_fraction_by_layer_run%zu", r));
    graph->SetMarkerStyle(20 + (r % 10));
    graph->SetMarkerSize(0.9);
    graph->SetMarkerColor(plotColor(r));
    graph->SetLineColor(plotColor(r));
    graph->SetLineWidth(2);
    graph->Draw("P SAME");
    legend.AddEntry(graph, runs[r].label.c_str(), "p");
    graphs.push_back(graph);
  }

  legend.Draw();
  const std::string baseName = joinPath(fractionDir, "mpv_fill_fraction_by_layer");
  canvas.SaveAs((baseName + ".pdf").c_str());
  canvas.SaveAs((baseName + ".png").c_str());

  rootFile.cd();
  canvas.Write();
  frame.Write();
  for (TGraph *graph : graphs) {
    graph->Write();
    delete graph;
  }
}

void drawStoredTrends(
    const std::vector<RunData> &runs,
    const std::set<int> &layers,
    const std::vector<std::vector<std::map<ChannelKey, StoredPoint> > > &stored,
    const std::string &outDir) {
  const std::string plotDir = joinPath(outDir, "StoredMIPTrends");
  gSystem->mkdir(plotDir.c_str(), true);
  TFile rootFile(joinPath(plotDir, "stored_mip_trends.root").c_str(), "RECREATE");
  drawStorageModeMaps(runs, layers, stored, plotDir, rootFile);
  drawLayerMedianTrends(runs, layers, stored, plotDir, rootFile);
  drawMpvFillFractionByLayer(runs, layers, stored, plotDir, rootFile);
  drawImputedChannelTrends(runs, stored, plotDir, rootFile);
  drawImputedChannelTrends(runs, stored, plotDir, rootFile,
                           "reference_channels", true);
  rootFile.Close();
}

void drawStorageModeMapsOnly(
    const std::vector<RunData> &runs,
    const std::set<int> &layers,
    const std::vector<std::vector<std::map<ChannelKey, StoredPoint> > > &stored,
    const std::string &outDir) {
  const std::string plotDir = joinPath(outDir, "StoredMIPTrends");
  gSystem->mkdir(plotDir.c_str(), true);
  TFile rootFile(joinPath(plotDir, "storage_mode_maps.root").c_str(), "RECREATE");
  drawStorageModeMaps(runs, layers, stored, plotDir, rootFile);
  rootFile.Close();
}

} // namespace mip_imputation

void apply_mip_imputation(
    const char *inputFile = "mip_fit_decision_by_run_out_8/FinalMIPCalibration.root",
    const char *outDir = "mip_imputed",
    const char *treeName = "mip_calibration",
    const char *jsonBaseDir = ".",
    const char *jsonSubdir = "json_neighborcheck_nofit",
    bool writeJson = true,
    int minEntries = 200,
    double minValidRunFraction = 0.30,
    const char *fullFitFile = "full_fit/FinalMIPCalibration.root",
    bool storageModeMapsOnly = false) {
  using namespace mip_imputation;
  gStyle->SetOptStat(0);

  std::vector<RunData> runs;
  if (!loadRuns(inputFile, treeName, runs)) return;
  if (runs.empty()) {
    std::cerr << "No runs loaded" << std::endl;
    return;
  }

  const std::string resolvedJsonBaseDir =
      jsonBaseDir && std::string(jsonBaseDir).size() ? jsonBaseDir : ".";
  const std::string resolvedJsonSubdir =
      jsonSubdir && std::string(jsonSubdir).size() ? jsonSubdir : "json_neighborcheck_nofit";
  const std::set<int> layers =
      collectOutputLayers(runs, resolvedJsonBaseDir, resolvedJsonSubdir);
  const std::set<ChannelKey> channelKeys = collectChannelKeys(runs);
  const std::string resolvedFullFitFile =
      resolveFullFitPath(inputFile ? inputFile : "", fullFitFile ? fullFitFile : "");
  const FullFitValueMap fullFitValues =
      loadFullFitValues(resolvedFullFitFile.c_str(), treeName);
  std::cout << "Loaded " << runs.size() << " run groups, "
            << channelKeys.size() << " channels, " << layers.size()
            << " layers from " << inputFile << std::endl;

  const std::vector<std::map<ChannelKey, ImputeStats> > stats =
      computeImputeStats(runs, channelKeys, minEntries, minValidRunFraction);
  const std::vector<std::vector<std::map<ChannelKey, StoredPoint> > > stored =
      buildStoredValues(runs, channelKeys, stats, fullFitValues,
                        minEntries, minValidRunFraction);
  if (storageModeMapsOnly) {
    drawStorageModeMapsOnly(runs, layers, stored, outDir);
    std::cout << "Storage mode maps written to "
              << joinPath(joinPath(outDir, "StoredMIPTrends"),
                          "storage_mode_maps")
              << std::endl;
    return;
  }
  if (writeJson) {
    writeJsonOutputs(runs, layers, stored, outDir,
                     resolvedJsonBaseDir,
                     resolvedJsonSubdir,
                     minEntries, minValidRunFraction);
    drawStoredTrends(runs, layers, stored, outDir);
    std::cout << "JSON and trend outputs written to " << outDir << std::endl;
  } else {
    std::cout << "JSON output skipped. Pass writeJson=true or --write-json to write files."
              << std::endl;
  }
}

#ifndef __CLING__
int main(int argc, char **argv) {
  const char *inputFile = "mip_fit_decision_by_run_out_8/FinalMIPCalibration.root";
  const char *outDir = "mip_imputed_last";
  const char *treeName = "mip_calibration";
  const char *jsonBaseDir = ".";
  const char *jsonSubdir = "json_neighborcheck_nofit";
  const char *fullFitFile = "full_fit/FinalMIPCalibration.root";
  bool writeJson = true;
  bool storageModeMapsOnly = false;
  int minEntries = 200;
  double minValidRunFraction = 0.30;

  int positional = 0;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--write-json" || arg == "-j") {
      writeJson = true;
    } else if (arg == "--no-write-json") {
      writeJson = false;
    } else if (arg == "--storage-mode-maps-only") {
      storageModeMapsOnly = true;
      writeJson = false;
    } else if (arg == "--tree-name" && i + 1 < argc) {
      treeName = argv[++i];
    } else if (arg == "--json-base-dir" && i + 1 < argc) {
      jsonBaseDir = argv[++i];
    } else if (arg == "--json-subdir" && i + 1 < argc) {
      jsonSubdir = argv[++i];
    } else if (arg == "--full-fit-file" && i + 1 < argc) {
      fullFitFile = argv[++i];
    } else if (arg == "--min-entries" && i + 1 < argc) {
      minEntries = std::atoi(argv[++i]);
    } else if (arg == "--min-valid-run-fraction" && i + 1 < argc) {
      minValidRunFraction = std::atof(argv[++i]);
    } else if (!arg.empty() && arg[0] == '-') {
      std::cerr << "Unknown option: " << arg << std::endl;
      std::cerr << "Usage: " << argv[0]
                << " [inputFile] [outDir] [jsonBaseDir]"
                << " [--tree-name mip_calibration]"
                << " [--json-subdir json_neighborcheck_nofit]"
                << " [--full-fit-file full_fit/FinalMIPCalibration.root]"
                << " [--min-entries 100]"
                << " [--min-valid-run-fraction 0.25]"
                << " [--storage-mode-maps-only]"
                << " [--no-write-json]" << std::endl;
      return 1;
    } else {
      if (positional == 0) {
        inputFile = argv[i];
      } else if (positional == 1) {
        outDir = argv[i];
      } else if (positional == 2) {
        jsonBaseDir = argv[i];
      } else {
        std::cerr << "Unexpected positional argument: " << arg << std::endl;
        return 1;
      }
      ++positional;
    }
  }

  apply_mip_imputation(inputFile, outDir, treeName, jsonBaseDir, jsonSubdir,
                       writeJson, minEntries, minValidRunFraction, fullFitFile,
                       storageModeMapsOnly);
  return 0;
}
#endif
