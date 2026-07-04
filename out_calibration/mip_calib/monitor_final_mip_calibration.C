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
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TKey.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

namespace final_mip_monitor {

const int kFirstFitSuccessDecision = 1;
const int kThresholdFitSuccessDecision = 2;
const int kSecondFitSuccessState = 1 << 16;
const int kEfficiencyRatioUnder0p98 = 1 << 27;

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
  double mpvError = 0.0;
  double width = -1.0;
  double gausSigma = -1.0;
  double directThreshold = -1.0;
  double directWidth = -1.0;
};

struct RunData {
  int runStart = -1;
  int runEnd = -1;
  std::string label;
  std::vector<Record> records;
  std::map<int, size_t> indexByCell;
};

std::string directoryName(const std::string &path);
std::string parentDirectory(const std::string &path);
std::string joinPath(const std::string &left, const std::string &right);
bool fileExists(const std::string &path);

std::string resolveHistBaseDir(const char *inputFile, const char *requestedBaseDir,
                               const std::vector<RunData> &runs,
                               const std::string &histFileName) {
  if (requestedBaseDir && std::string(requestedBaseDir).size()) return requestedBaseDir;
  std::vector<std::string> candidates;
  candidates.push_back(directoryName(inputFile));
  candidates.push_back(parentDirectory(inputFile));
  candidates.push_back(".");
  for (const std::string &candidate : candidates) {
    if (runs.empty()) return candidate;
    const std::string probe = joinPath(joinPath(candidate, runs[0].label), histFileName);
    if (fileExists(probe)) return candidate;
  }
  return candidates.front();
}

struct VariableInfo {
  const char *name;
  const char *title;
  double Record::*member;
  bool direct;
};

const std::vector<VariableInfo> &variables() {
  static const std::vector<VariableInfo> vars = {
      {"mpv", "Adopted MPV", &Record::mpv, false},
      {"width", "Adopted width", &Record::width, false},
      {"gaus_sigma", "Adopted Gaussian sigma", &Record::gausSigma, false},
      {"direct_threshold", "Adopted direct threshold", &Record::directThreshold, true},
      {"direct_width", "Adopted direct width", &Record::directWidth, true}};
  return vars;
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
  double sum = 0.0, sum2 = 0.0;
  for (double value : values) { sum += value; sum2 += value * value; }
  result.mean = sum / values.size();
  result.rms = std::sqrt(std::max(0.0, sum2 / values.size() - result.mean * result.mean));
  result.median = quantile(values, 0.50);
  result.q16 = quantile(values, 0.16);
  result.q84 = quantile(values, 0.84);
  result.minimum = values.front();
  result.maximum = values.back();
  return result;
}

bool validValue(double value) { return std::isfinite(value) && value >= 0.0; }

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

std::string parentDirectory(const std::string &path) {
  const std::string dir = directoryName(path);
  if (dir == path || dir == "." || dir == "/") return dir;
  return directoryName(dir);
}

bool fileExists(const std::string &path) {
  return gSystem->AccessPathName(path.c_str()) == false;
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

bool belongsToGroup(const Record &record, const std::string &group,
                    const VariableInfo &variable) {
  const double value = record.*(variable.member);
  if (!validValue(value)) return false;
  if (group == "all_adopted") return !variable.direct;
  if (group == "first_fit") return !variable.direct &&
      record.decision == kFirstFitSuccessDecision;
  if (group == "second_fit") return !variable.direct &&
      (record.state & kSecondFitSuccessState);
  if (group == "threshold_fit") return variable.direct &&
      record.decision == kThresholdFitSuccessDecision &&
      (record.state & kEfficiencyRatioUnder0p98);
  return false;
}

bool isSuccessRecord(const Record &record) {
  return record.decision == kFirstFitSuccessDecision ||
         (record.state & kSecondFitSuccessState) ||
         (record.decision == kThresholdFitSuccessDecision &&
          (record.state & kEfficiencyRatioUnder0p98));
}

bool isFirstFitSuccessRecord(const Record &record) {
  return record.decision == kFirstFitSuccessDecision;
}

std::vector<std::string> groupsFor(const VariableInfo &variable) {
  if (variable.direct) return {"threshold_fit"};
  return {"all_adopted", "first_fit", "second_fit"};
}

std::vector<double> collect(const RunData &run, const VariableInfo &variable,
                            const std::string &group) {
  std::vector<double> values;
  for (const Record &record : run.records)
    if (belongsToGroup(record, group, variable)) values.push_back(record.*(variable.member));
  return values;
}

std::vector<double> collectLayerMpv(const RunData &run, int layer) {
  std::vector<double> values;
  for (const Record &record : run.records) {
    if (record.layer == layer && validValue(record.mpv)) values.push_back(record.mpv);
  }
  return values;
}

std::string runKey(int start, int end, const std::string &label) {
  return std::to_string(start) + ":" + std::to_string(end) + ":" + label;
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
    file->Close(); delete file;
    return false;
  }

  const char *required[] = {"run_start", "run_end", "cellid", "layer", "chip", "channel",
                            "entries", "decision", "state", "mpv", "width", "gaus_sigma",
                            "direct_threshold", "direct_width"};
  for (const char *branch : required) {
    if (!tree->GetBranch(branch)) {
      std::cerr << "Missing required branch: " << branch << std::endl;
      file->Close(); delete file;
      return false;
    }
  }

  Record value;
  std::string *runLabel = nullptr, *decisionName = nullptr, *stateNames = nullptr;
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
  if (tree->GetBranch("mpv_error")) tree->SetBranchAddress("mpv_error", &value.mpvError);
  tree->SetBranchAddress("width", &value.width);
  tree->SetBranchAddress("gaus_sigma", &value.gausSigma);
  tree->SetBranchAddress("direct_threshold", &value.directThreshold);
  tree->SetBranchAddress("direct_width", &value.directWidth);
  if (tree->GetBranch("run_label")) tree->SetBranchAddress("run_label", &runLabel);
  if (tree->GetBranch("decision_name")) tree->SetBranchAddress("decision_name", &decisionName);
  if (tree->GetBranch("state_names")) tree->SetBranchAddress("state_names", &stateNames);

  std::map<std::string, size_t> runIndex;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    tree->GetEntry(entry);
    value.runLabel = runLabel ? *runLabel : std::to_string(value.runStart) + "-" + std::to_string(value.runEnd);
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
  file->Close(); delete file;

  std::sort(runs.begin(), runs.end(), [](const RunData &a, const RunData &b) {
    if (a.runStart != b.runStart) return a.runStart < b.runStart;
    if (a.runEnd != b.runEnd) return a.runEnd < b.runEnd;
    return a.label < b.label;
  });
  for (RunData &run : runs)
    for (size_t i = 0; i < run.records.size(); ++i) run.indexByCell[run.records[i].cellid] = i;
  return !runs.empty();
}

std::pair<double, double> histogramRange(const std::vector<double> &values) {
  const double low = *std::min_element(values.begin(), values.end());
  const double high = *std::max_element(values.begin(), values.end());
  const double padding = high > low ? 0.03 * (high - low) : std::max(1.0, 0.05 * std::abs(low));
  return std::make_pair(low - padding, high + padding);
}

int runColor(size_t index) {
  static const int colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1,
                               kOrange + 7, kCyan + 2, kViolet + 1, kTeal + 3,
                               kPink + 7, kAzure + 7, kSpring + 5, kGray + 2};
  return colors[index % (sizeof(colors) / sizeof(colors[0]))];
}

void writeSummary(const std::vector<RunData> &runs, const std::string &outDir) {
  std::ofstream csv((outDir + "/FinalValueSummary_byRun.csv").c_str());
  csv << "run_start,run_end,run_label,total_channels,group,variable,valid_channels,valid_fraction,mean,rms,median,q16,q84,min,max\n";
  for (const RunData &run : runs)
    for (const VariableInfo &variable : variables())
      for (const std::string &group : groupsFor(variable)) {
        const Summary summary = summarize(collect(run, variable, group));
        if (!summary.n) continue;
        csv << run.runStart << "," << run.runEnd << "," << run.label << ","
            << run.records.size() << "," << group << "," << variable.name << "," << summary.n
            << "," << (run.records.empty() ? 0.0 : double(summary.n) / run.records.size()) << ","
            << std::setprecision(12)
            << summary.mean << "," << summary.rms << "," << summary.median << ","
            << summary.q16 << "," << summary.q84 << "," << summary.minimum << ","
            << summary.maximum << "\n";
      }
}

std::vector<std::string> splitStates(const std::string &states) {
  std::vector<std::string> result;
  std::string token;
  std::stringstream stream(states);
  while (std::getline(stream, token, '|')) {
    const size_t first = token.find_first_not_of(" \t");
    const size_t last = token.find_last_not_of(" \t");
    if (first != std::string::npos) result.push_back(token.substr(first, last - first + 1));
  }
  return result;
}

void writeStateSummary(const std::vector<RunData> &runs, const std::string &outDir) {
  std::ofstream csv((outDir + "/FinalValueSummary_byRunState.csv").c_str());
  csv << "run_start,run_end,run_label,state,state_channels,variable,valid_channels,valid_fraction_in_state,mean,rms,median,q16,q84\n";
  for (const RunData &run : runs) {
    std::map<std::string, std::map<std::string, std::vector<double> > > values;
    std::map<std::string, size_t> stateCounts;
    for (const Record &record : run.records) {
      for (const std::string &state : splitStates(record.stateNames)) {
        stateCounts[state]++;
        for (const VariableInfo &variable : variables()) {
          const double value = record.*(variable.member);
          if (validValue(value)) values[state][variable.name].push_back(value);
        }
      }
    }
    for (const auto &state : values)
      for (const auto &variable : state.second) {
        const Summary summary = summarize(variable.second);
        csv << run.runStart << "," << run.runEnd << "," << run.label << "," << state.first
            << "," << stateCounts[state.first] << "," << variable.first << "," << summary.n
            << "," << (stateCounts[state.first] ? double(summary.n) / stateCounts[state.first] : 0.0)
            << "," << std::setprecision(12)
            << summary.mean << "," << summary.rms << "," << summary.median << ","
            << summary.q16 << "," << summary.q84 << "\n";
      }
  }
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

std::map<int, std::string> collectDecisionNames(const std::vector<RunData> &runs) {
  std::map<int, std::string> names;
  for (const RunData &run : runs) {
    for (const Record &record : run.records) {
      if (!names.count(record.decision) || names[record.decision].empty()) {
        names[record.decision] = record.decisionName.empty()
            ? std::to_string(record.decision) : record.decisionName;
      }
    }
  }
  return names;
}

std::map<int, size_t> countDecisionsInLayer(const RunData &run, int layer,
                                            size_t &total) {
  std::map<int, size_t> counts;
  total = 0;
  for (const Record &record : run.records) {
    if (record.layer != layer) continue;
    ++total;
    counts[record.decision]++;
  }
  return counts;
}

std::string csvQuote(const std::string &text) {
  if (text.find_first_of(",\"\n") == std::string::npos) return text;
  std::string escaped = "\"";
  for (char c : text) {
    if (c == '"') escaped += "\"\"";
    else escaped += c;
  }
  escaped += "\"";
  return escaped;
}

void writeDecisionLayerSummary(const std::vector<RunData> &runs,
                               const std::string &outDir) {
  const std::set<int> layers = collectLayers(runs);
  const std::map<int, std::string> decisionNames = collectDecisionNames(runs);
  std::map<int, int> layerBins;
  int layerBin = 1;
  for (int layer : layers) layerBins[layer] = layerBin++;
  std::map<int, int> decisionBins;
  int decisionBin = 1;
  for (const auto &decision : decisionNames) decisionBins[decision.first] = decisionBin++;

  std::ofstream csv((outDir + "/FinalDecisionFraction_byRunLayer.csv").c_str());
  TFile rootFile((outDir + "/FinalDecisionFraction_byRunLayer.root").c_str(), "RECREATE");
  const std::string stackDir = outDir + "/DecisionLayerSummaryStacks";
  const std::string successDir = outDir + "/RunLayerSuccessRate";
  gSystem->mkdir(stackDir.c_str(), true);
  gSystem->mkdir(successDir.c_str(), true);
  csv << "run_start,run_end,run_label,layer,total_channels,decision,decision_name,channels,fraction\n";
  if (layers.empty() || decisionNames.empty()) {
    rootFile.Close();
    return;
  }
  for (size_t r = 0; r < runs.size(); ++r) {
    const RunData &run = runs[r];
    TH2D decisionLayerSummary(Form("DecisionLayerSummary_run%zu", r),
                              Form("Decision layer summary: %s;Layer;Decision;Fraction",
                                   run.label.c_str()),
                              layers.size(), 0.5, layers.size() + 0.5,
                              decisionNames.size(), 0.5, decisionNames.size() + 0.5);
    for (const auto &layer : layerBins) {
      decisionLayerSummary.GetXaxis()->SetBinLabel(layer.second, Form("%d", layer.first));
    }
    for (const auto &decision : decisionNames) {
      decisionLayerSummary.GetYaxis()->SetBinLabel(
          decisionBins[decision.first],
          Form("%s (%d)", decision.second.c_str(), decision.first));
    }
    for (int layer : layers) {
      size_t total = 0;
      const std::map<int, size_t> counts = countDecisionsInLayer(run, layer, total);
      if (!total) continue;
      for (const auto &decision : decisionNames) {
        const auto count = counts.find(decision.first);
        const size_t channels = count == counts.end() ? 0 : count->second;
        const double fraction = double(channels) / total;
        csv << run.runStart << "," << run.runEnd << "," << csvQuote(run.label) << ","
            << layer << "," << total << "," << decision.first << ","
            << csvQuote(decision.second) << "," << channels << ","
            << std::setprecision(12) << fraction << "\n";
        decisionLayerSummary.SetBinContent(layerBins[layer], decisionBins[decision.first], fraction);
      }
    }
    rootFile.cd();
    decisionLayerSummary.Write();
  }

  rootFile.mkdir("RunLayerSuccessRate");
  rootFile.cd("RunLayerSuccessRate");
  for (size_t r = 0; r < runs.size(); ++r) {
    const RunData &run = runs[r];
    TCanvas canvas(Form("c_success_rate_by_layer_run%zu", r), "", 1200, 750);
    TH1D hist(Form("h_success_rate_by_layer_run%zu", r),
              Form("Success rate by layer: %s;Layer;Success rate", run.label.c_str()),
              layers.size(), 0.5, layers.size() + 0.5);
    TH1D ineffHist(Form("h_inefficiency_by_layer_run%zu", r),
                   Form("Inefficiency by layer: %s;Layer;Inefficiency", run.label.c_str()),
                   layers.size(), 0.5, layers.size() + 0.5);
    TH1D firstFitHist(Form("h_first_fit_success_rate_by_layer_run%zu", r),
                      Form("First-fit success rate by layer: %s;Layer;First-fit success rate",
                           run.label.c_str()),
                      layers.size(), 0.5, layers.size() + 0.5);
    hist.SetMinimum(0.0);
    hist.SetMaximum(1.08);
    hist.SetMarkerStyle(20);
    hist.SetMarkerSize(1.0);
    hist.SetMarkerColor(kBlue + 1);
    hist.SetLineColor(kBlue + 1);
    hist.SetLineWidth(2);
    ineffHist.SetMinimum(1.0e-4);
    ineffHist.SetMaximum(1.0);
    ineffHist.SetMarkerStyle(20);
    ineffHist.SetMarkerSize(1.0);
    ineffHist.SetMarkerColor(kRed + 1);
    ineffHist.SetLineColor(kRed + 1);
    ineffHist.SetLineWidth(2);
    firstFitHist.SetMinimum(0.0);
    firstFitHist.SetMaximum(1.08);
    firstFitHist.SetMarkerStyle(21);
    firstFitHist.SetMarkerSize(1.0);
    firstFitHist.SetMarkerColor(kGreen + 2);
    firstFitHist.SetLineColor(kGreen + 2);
    firstFitHist.SetLineWidth(2);

    double successRateSum = 0.0;
    double firstFitSuccessRateSum = 0.0;
    size_t layersWithChannels = 0;
    for (int layer : layers) {
      size_t total = 0;
      size_t success = 0;
      size_t firstFitSuccess = 0;
      for (const Record &record : run.records) {
        if (record.layer != layer) continue;
        ++total;
        if (isSuccessRecord(record)) ++success;
        if (isFirstFitSuccessRecord(record)) ++firstFitSuccess;
      }
      hist.GetXaxis()->SetBinLabel(layerBins[layer], Form("%d", layer));
      ineffHist.GetXaxis()->SetBinLabel(layerBins[layer], Form("%d", layer));
      firstFitHist.GetXaxis()->SetBinLabel(layerBins[layer], Form("%d", layer));
      if (total) {
        const double successRate = double(success) / total;
        const double firstFitSuccessRate = double(firstFitSuccess) / total;
        hist.SetBinContent(layerBins[layer], successRate);
        ineffHist.SetBinContent(layerBins[layer], std::max(0.0, 1.0 - successRate));
        firstFitHist.SetBinContent(layerBins[layer], firstFitSuccessRate);
        successRateSum += successRate;
        firstFitSuccessRateSum += firstFitSuccessRate;
        ++layersWithChannels;
      }
    }
    const double meanSuccessRate = layersWithChannels ? successRateSum / layersWithChannels : 0.0;
    const double meanFirstFitSuccessRate =
        layersWithChannels ? firstFitSuccessRateSum / layersWithChannels : 0.0;

    canvas.SetBottomMargin(0.18);
    canvas.SetLeftMargin(0.10);
    canvas.SetGrid();
    hist.Draw("P HIST");
    firstFitHist.Draw("P HIST SAME");
    hist.GetXaxis()->LabelsOption("v");
    TLegend successLegend(0.58, 0.25, 0.89, 0.39);
    successLegend.SetBorderSize(0);
    successLegend.AddEntry(&hist, Form("all success (mean %.3f)", meanSuccessRate), "pl");
    successLegend.AddEntry(&firstFitHist,
                           Form("first-fit success (mean %.3f)", meanFirstFitSuccessRate), "pl");
    successLegend.Draw();
    const std::string baseName = successDir + "/success_rate_" + fileSafeString(run.label);
    canvas.SaveAs((baseName + ".pdf").c_str());
    canvas.SaveAs((baseName + ".png").c_str());
    canvas.Write();
    hist.Write();

    TCanvas ineffCanvas(Form("c_inefficiency_by_layer_run%zu", r), "", 1200, 750);
    ineffCanvas.SetBottomMargin(0.18);
    ineffCanvas.SetLeftMargin(0.10);
    ineffCanvas.SetGrid();
    ineffCanvas.SetLogy();
    ineffHist.Draw("P HIST");
    ineffHist.GetXaxis()->LabelsOption("v");
    const std::string ineffBaseName = successDir + "/inefficiency_" + fileSafeString(run.label);
    ineffCanvas.SaveAs((ineffBaseName + ".pdf").c_str());
    ineffCanvas.SaveAs((ineffBaseName + ".png").c_str());
    ineffCanvas.Write();
    ineffHist.Write();

    TCanvas firstFitCanvas(Form("c_first_fit_success_rate_by_layer_run%zu", r), "", 1200, 750);
    firstFitCanvas.SetBottomMargin(0.18);
    firstFitCanvas.SetLeftMargin(0.10);
    firstFitCanvas.SetGrid();
    firstFitHist.Draw("P HIST");
    firstFitHist.GetXaxis()->LabelsOption("v");
    const std::string firstFitBaseName =
        successDir + "/first_fit_success_rate_" + fileSafeString(run.label);
    firstFitCanvas.SaveAs((firstFitBaseName + ".pdf").c_str());
    firstFitCanvas.SaveAs((firstFitBaseName + ".png").c_str());
    firstFitCanvas.Write();
    firstFitHist.Write();
  }

  rootFile.mkdir("DecisionLayerSummaryStacks");
  rootFile.cd("DecisionLayerSummaryStacks");
  for (int layer : layers) {
    TCanvas canvas(Form("c_DecisionLayerSummary_layer%02d", layer), "", 1200, 750);
    THStack stack(Form("DecisionLayerSummary_layer%02d", layer),
                  Form("Decision layer summary: layer %02d;Run;Fraction", layer));
    TLegend legend(0.58, 0.55, 0.92, 0.89);
    legend.SetBorderSize(0);
    if (decisionNames.size() > 6) legend.SetNColumns(2);

    std::vector<TH1D *> histograms;
    size_t decisionIndex = 0;
    for (const auto &decision : decisionNames) {
      TH1D *hist = new TH1D(Form("h_DecisionLayerSummary_layer%02d_decision_%d",
                                  layer, decision.first),
                            decision.second.c_str(),
                            runs.size(), 0.5, runs.size() + 0.5);
      for (size_t r = 0; r < runs.size(); ++r) {
        hist->GetXaxis()->SetBinLabel(r + 1, runs[r].label.c_str());
        size_t total = 0;
        const std::map<int, size_t> counts = countDecisionsInLayer(runs[r], layer, total);
        if (!total) continue;
        const auto count = counts.find(decision.first);
        const size_t channels = count == counts.end() ? 0 : count->second;
        hist->SetBinContent(r + 1, double(channels) / total);
      }
      hist->SetFillColor(runColor(decisionIndex));
      hist->SetLineColor(kBlack);
      hist->SetLineWidth(1);
      stack.Add(hist);
      legend.AddEntry(hist, Form("%s (%d)", decision.second.c_str(), decision.first), "f");
      histograms.push_back(hist);
      ++decisionIndex;
    }

    canvas.SetBottomMargin(0.22);
    canvas.SetLeftMargin(0.10);
    canvas.SetGrid();
    stack.SetMinimum(0.0);
    stack.SetMaximum(1.08);
    stack.Draw("hist");
    stack.GetXaxis()->LabelsOption("v");
    stack.GetYaxis()->SetTitleOffset(1.15);
    legend.Draw();
    canvas.SaveAs((stackDir + Form("/DecisionLayerSummary_layer%02d.png", layer)).c_str());
    canvas.Write();
    stack.Write();
    for (TH1D *hist : histograms) {
      hist->Write();
      delete hist;
    }
  }
  rootFile.Close();
}

void drawTrends(const std::vector<RunData> &runs, const std::string &outDir) {
  const std::string plotDir = outDir + "/Trends";
  gSystem->mkdir(plotDir.c_str(), true);
  TFile rootFile((plotDir + "/final_value_trends.root").c_str(), "RECREATE");
  for (const VariableInfo &variable : variables()) {
    for (const std::string &group : groupsFor(variable)) {
      std::vector<double> x, y, exLow, exHigh, eyLow, eyHigh;
      for (size_t r = 0; r < runs.size(); ++r) {
        const Summary summary = summarize(collect(runs[r], variable, group));
        if (!summary.n) continue;
        x.push_back(r + 1.0); y.push_back(summary.median);
        exLow.push_back(0.0); exHigh.push_back(0.0);
        eyLow.push_back(summary.median - summary.q16);
        eyHigh.push_back(summary.q84 - summary.median);
      }
      if (x.empty()) continue;
      const std::string groupDir = plotDir + "/" + group;
      gSystem->mkdir(groupDir.c_str(), true);
      TCanvas canvas(("c_trend_" + group + "_" + variable.name).c_str(), "", 1200, 750);
      TGraphAsymmErrors graph(x.size(), &x[0], &y[0], &exLow[0], &exHigh[0], &eyLow[0], &eyHigh[0]);
      graph.SetName(("g_trend_" + group + "_" + variable.name).c_str());
      graph.SetMarkerStyle(20); graph.SetLineWidth(2);
      double low = y[0] - eyLow[0], high = y[0] + eyHigh[0];
      for (size_t i = 1; i < y.size(); ++i) {
        low = std::min(low, y[i] - eyLow[i]); high = std::max(high, y[i] + eyHigh[i]);
      }
      const double padding = high > low ? 0.10 * (high - low) : 1.0;
      TH1D frame(("frame_trend_" + group + "_" + variable.name).c_str(),
                 (std::string(variable.title) + " trend: " + group +
                  " (median, 16-84%);Run;" + variable.name).c_str(),
                 runs.size(), 0.5, runs.size() + 0.5);
      for (size_t r = 0; r < runs.size(); ++r) frame.GetXaxis()->SetBinLabel(r + 1, runs[r].label.c_str());
      frame.GetXaxis()->LabelsOption("v"); frame.SetMinimum(low - padding); frame.SetMaximum(high + padding);
      frame.Draw(); graph.Draw("P SAME");
      canvas.SetBottomMargin(0.22); canvas.SetGrid();
      canvas.SaveAs((groupDir + "/" + variable.name + "_trend.pdf").c_str());
      canvas.SaveAs((groupDir + "/" + variable.name + "_trend.png").c_str());
      rootFile.cd(); graph.Write();
    }
  }
  std::set<int> layers;
  for (const RunData &run : runs) {
    for (const Record &record : run.records) {
      if (record.layer >= 0 && validValue(record.mpv)) layers.insert(record.layer);
    }
  }
  if (!layers.empty()) {
    const std::string layerDir = plotDir + "/by_layer";
    gSystem->mkdir(layerDir.c_str(), true);
    for (int layer : layers) {
      std::vector<double> x, y, exLow, exHigh, eyLow, eyHigh;
      for (size_t r = 0; r < runs.size(); ++r) {
        const Summary summary = summarize(collectLayerMpv(runs[r], layer));
        if (!summary.n) continue;
        x.push_back(r + 1.0);
        y.push_back(summary.median);
        exLow.push_back(0.0);
        exHigh.push_back(0.0);
        eyLow.push_back(summary.median - summary.q16);
        eyHigh.push_back(summary.q84 - summary.median);
      }
      if (x.empty()) continue;
      TCanvas canvas(Form("c_trend_layer_%02d_mpv", layer), "", 1200, 750);
      TGraphAsymmErrors graph(x.size(), &x[0], &y[0], &exLow[0], &exHigh[0],
                              &eyLow[0], &eyHigh[0]);
      graph.SetName(Form("g_trend_layer_%02d_mpv", layer));
      graph.SetMarkerStyle(20);
      graph.SetLineWidth(2);
      double low = y[0] - eyLow[0], high = y[0] + eyHigh[0];
      for (size_t i = 1; i < y.size(); ++i) {
        low = std::min(low, y[i] - eyLow[i]);
        high = std::max(high, y[i] + eyHigh[i]);
      }
      const double padding = high > low ? 0.10 * (high - low) : 1.0;
      TH1D frame(Form("frame_trend_layer_%02d_mpv", layer),
                 Form("MPV trend: layer %02d (median, 16-84%%);Run;mpv", layer),
                 runs.size(), 0.5, runs.size() + 0.5);
      for (size_t r = 0; r < runs.size(); ++r) {
        frame.GetXaxis()->SetBinLabel(r + 1, runs[r].label.c_str());
      }
      frame.GetXaxis()->LabelsOption("v");
      frame.SetMinimum(low - padding);
      frame.SetMaximum(high + padding);
      frame.Draw();
      graph.Draw("P SAME");
      canvas.SetBottomMargin(0.22);
      canvas.SetGrid();
      canvas.SaveAs((layerDir + Form("/mpv_layer%02d_trend.pdf", layer)).c_str());
      canvas.SaveAs((layerDir + Form("/mpv_layer%02d_trend.png", layer)).c_str());
      rootFile.cd();
      graph.Write();
    }
  }
  const std::set<int> decisionLayers = collectLayers(runs);
  const std::map<int, std::string> decisionNames = collectDecisionNames(runs);
  if (!decisionLayers.empty() && !decisionNames.empty()) {
    const std::string decisionDir = plotDir + "/decision_fraction_by_layer";
    gSystem->mkdir(decisionDir.c_str(), true);
    for (int layer : decisionLayers) {
      TCanvas canvas(Form("c_decision_fraction_layer_%02d", layer), "", 1200, 750);
      TH1D frame(Form("frame_decision_fraction_layer_%02d", layer),
                 Form("Final decision fraction trend: layer %02d;Run;Fraction", layer),
                 runs.size(), 0.5, runs.size() + 0.5);
      for (size_t r = 0; r < runs.size(); ++r) {
        frame.GetXaxis()->SetBinLabel(r + 1, runs[r].label.c_str());
      }
      frame.GetXaxis()->LabelsOption("v");
      frame.SetMinimum(0.0);
      frame.SetMaximum(1.08);
      frame.Draw();

      TLegend legend(0.58, 0.55, 0.92, 0.89);
      legend.SetBorderSize(0);
      if (decisionNames.size() > 6) legend.SetNColumns(2);
      std::vector<TGraphAsymmErrors *> graphs;
      size_t decisionIndex = 0;
      for (const auto &decision : decisionNames) {
        std::vector<double> x, y, exLow, exHigh, eyLow, eyHigh;
        for (size_t r = 0; r < runs.size(); ++r) {
          size_t total = 0;
          const std::map<int, size_t> counts = countDecisionsInLayer(runs[r], layer, total);
          if (!total) continue;
          const auto count = counts.find(decision.first);
          const size_t channels = count == counts.end() ? 0 : count->second;
          x.push_back(r + 1.0);
          y.push_back(double(channels) / total);
          exLow.push_back(0.0);
          exHigh.push_back(0.0);
          eyLow.push_back(0.0);
          eyHigh.push_back(0.0);
        }
        if (x.empty()) {
          ++decisionIndex;
          continue;
        }
        TGraphAsymmErrors *graph = new TGraphAsymmErrors(
            x.size(), &x[0], &y[0], &exLow[0], &exHigh[0], &eyLow[0], &eyHigh[0]);
        graph->SetName(Form("g_decision_fraction_layer_%02d_decision_%d",
                            layer, decision.first));
        graph->SetMarkerStyle(20 + (decisionIndex % 10));
        graph->SetMarkerSize(0.9);
        graph->SetLineWidth(2);
        graph->SetLineColor(runColor(decisionIndex));
        graph->SetMarkerColor(runColor(decisionIndex));
        graph->Draw("PL SAME");
        legend.AddEntry(graph,
                        Form("%s (%d)", decision.second.c_str(), decision.first),
                        "pl");
        graphs.push_back(graph);
        ++decisionIndex;
      }
      legend.Draw();
      canvas.SetBottomMargin(0.22);
      canvas.SetGrid();
      canvas.SaveAs((decisionDir + Form("/final_decision_fraction_layer%02d_trend.pdf", layer)).c_str());
      canvas.SaveAs((decisionDir + Form("/final_decision_fraction_layer%02d_trend.png", layer)).c_str());
      rootFile.cd();
      for (TGraphAsymmErrors *graph : graphs) {
        graph->Write();
        delete graph;
      }
    }
  }
  rootFile.Close();
}

void drawOverlays(const std::vector<RunData> &runs, const std::string &outDir) {
  const std::string plotDir = outDir + "/OverlayDistributions";
  const std::string layerPlotDir = outDir + "/LayerOverlayDistributions";
  gSystem->mkdir(plotDir.c_str(), true);
  gSystem->mkdir(layerPlotDir.c_str(), true);
  TFile rootFile((plotDir + "/final_value_distributions.root").c_str(), "RECREATE");
  TFile layerRootFile((layerPlotDir + "/final_value_layer_distributions.root").c_str(), "RECREATE");
  const std::set<int> layers = collectLayers(runs);
  for (const VariableInfo &variable : variables()) {
    for (const std::string &group : groupsFor(variable)) {
      std::vector<std::vector<double> > byRun(runs.size());
      std::vector<double> all;
      for (size_t r = 0; r < runs.size(); ++r) {
        byRun[r] = collect(runs[r], variable, group);
        all.insert(all.end(), byRun[r].begin(), byRun[r].end());
      }
      if (all.empty()) continue;
      const std::pair<double, double> range = histogramRange(all);
      const std::string groupDir = plotDir + "/" + group;
      gSystem->mkdir(groupDir.c_str(), true);
      TCanvas canvas(("c_overlay_" + group + "_" + variable.name).c_str(), "", 1000, 800);
      TLegend legend(0.67, 0.60, 0.92, 0.89); legend.SetBorderSize(0);
      if (runs.size() > 6) legend.SetNColumns(2);
      std::vector<TH1D *> histograms;
      double maximum = 0.0;
      for (size_t r = 0; r < runs.size(); ++r) {
        if (byRun[r].empty()) continue;
        TH1D *hist = new TH1D(("h_overlay_" + group + "_" + variable.name + "_" + std::to_string(r)).c_str(),
                              "", 60, range.first, range.second);
        hist->SetDirectory(nullptr); hist->SetStats(false);
        for (double value : byRun[r]) hist->Fill(value);
        hist->Scale(1.0 / hist->Integral(), "width"); hist->SetLineColor(runColor(r)); hist->SetLineWidth(2);
        maximum = std::max(maximum, hist->GetMaximum()); histograms.push_back(hist);
        legend.AddEntry(hist, (runs[r].label + " (N=" + std::to_string(byRun[r].size()) + ")").c_str(), "l");
      }
      histograms[0]->SetTitle((std::string(variable.title) + " distributions: " + group +
                               ";" + variable.name + ";Density").c_str());
      histograms[0]->SetMaximum(1.18 * maximum); histograms[0]->Draw("HIST");
      for (size_t i = 1; i < histograms.size(); ++i) histograms[i]->Draw("HIST SAME");
      legend.Draw();
      canvas.SaveAs((groupDir + "/" + variable.name + "_overlay.pdf").c_str());
      canvas.SaveAs((groupDir + "/" + variable.name + "_overlay.png").c_str());
      rootFile.cd(); for (TH1D *hist : histograms) { hist->Write(); delete hist; }

      const std::string layerGroupDir = layerPlotDir + "/" + group;
      gSystem->mkdir(layerGroupDir.c_str(), true);
      for (int layer : layers) {
        std::vector<std::vector<double> > byRunInLayer(runs.size());
        std::vector<double> allInLayer;
        for (size_t r = 0; r < runs.size(); ++r) {
          for (const Record &record : runs[r].records) {
            if (record.layer != layer || !belongsToGroup(record, group, variable)) continue;
            byRunInLayer[r].push_back(record.*(variable.member));
            allInLayer.push_back(record.*(variable.member));
          }
        }
        if (allInLayer.empty()) continue;
        const std::pair<double, double> layerRange = histogramRange(allInLayer);

        TCanvas layerCanvas(Form("c_layer_run_overlay_%s_%s_layer%02d",
                                 group.c_str(), variable.name, layer),
                            "", 1200, 850);
        TLegend layerLegend(0.48, 0.50, 0.94, 0.89);
        layerLegend.SetBorderSize(0);
        layerLegend.SetTextSize(0.025);
        std::vector<TH1D *> layerHistograms;
        double layerMaximum = 0.0;
        for (size_t r = 0; r < runs.size(); ++r) {
          if (byRunInLayer[r].empty()) continue;
          TH1D *hist = new TH1D(
              Form("h_layer_run_overlay_%s_%s_layer%02d_run%zu",
                   group.c_str(), variable.name, layer, r),
              "", 60, layerRange.first, layerRange.second);
          hist->SetDirectory(nullptr);
          hist->SetStats(false);
          for (double value : byRunInLayer[r]) hist->Fill(value);
          if (hist->Integral() > 0.0) hist->Scale(1.0 / hist->Integral(), "width");
          hist->SetLineColor(runColor(r));
          hist->SetLineWidth(2);
          layerMaximum = std::max(layerMaximum, hist->GetMaximum());
          layerHistograms.push_back(hist);
          const Summary summary = summarize(byRunInLayer[r]);
          layerLegend.AddEntry(hist,
                               Form("avg=%.3g RMS=%.3g N=%zu  %s",
                                    summary.mean, summary.rms,
                                    byRunInLayer[r].size(), runs[r].label.c_str()),
                               "l");
        }
        if (layerHistograms.empty()) continue;
        layerHistograms[0]->SetTitle(
            (std::string(variable.title) + " run distributions: " + group +
             Form(", layer %02d;", layer) + variable.name + ";Density").c_str());
        layerHistograms[0]->SetMaximum(1.18 * layerMaximum);
        layerHistograms[0]->Draw("HIST");
        for (size_t i = 1; i < layerHistograms.size(); ++i) {
          layerHistograms[i]->Draw("HIST SAME");
        }
        layerLegend.Draw();
        const std::string layerBaseName =
            layerGroupDir + "/" + variable.name + Form("_layer%02d_run_overlay", layer);
        layerCanvas.SaveAs((layerBaseName + ".pdf").c_str());
        layerCanvas.SaveAs((layerBaseName + ".png").c_str());
        layerRootFile.cd();
        layerCanvas.Write();
        for (TH1D *hist : layerHistograms) {
          hist->Write();
          delete hist;
        }
      }
    }
  }
  rootFile.Close();
  layerRootFile.Close();
}

struct DeltaRecord {
  const Record *current = nullptr;
  double delta = 0.0;
  double relative = 0.0;
};

std::vector<DeltaRecord> matchedDeltas(const RunData &reference, const RunData &current,
                                       const VariableInfo &variable, const std::string &group) {
  std::vector<DeltaRecord> result;
  for (const Record &record : current.records) {
    if (!belongsToGroup(record, group, variable)) continue;
    const auto referenceIndex = reference.indexByCell.find(record.cellid);
    if (referenceIndex == reference.indexByCell.end()) continue;
    const Record &old = reference.records[referenceIndex->second];
    if (!belongsToGroup(old, group, variable)) continue;
    const double oldValue = old.*(variable.member);
    if (std::abs(oldValue) < 1e-12) continue;
    const double delta = record.*(variable.member) - oldValue;
    DeltaRecord item; item.current = &record; item.delta = delta; item.relative = delta / oldValue;
    result.push_back(item);
  }
  return result;
}

void writeAndDrawDeltas(const std::vector<RunData> &runs, const std::string &outDir,
                        double relativeThreshold, double robustZThreshold) {
  if (runs.size() < 2) return;
  const std::string plotDir = outDir + "/DeltaFromFirstRun";
  gSystem->mkdir(plotDir.c_str(), true);
  TFile rootFile((plotDir + "/delta_from_first_run.root").c_str(), "RECREATE");
  std::ofstream summary((outDir + "/DeltaSummary_fromFirstRun.csv").c_str());
  std::ofstream outliers((outDir + "/LargeShiftChannels.csv").c_str());
  summary << "reference_run,run,group,variable,matched_channels,mean_delta,rms_delta,median_delta,q16_delta,q84_delta,median_relative_delta,median_relative_shift_flag\n";
  outliers << "reference_run,run,group,variable,cellid,layer,chip,channel,decision,decision_name,state,state_names,value,delta,relative_delta,robust_z,relative_flag,robust_flag\n";

  for (const VariableInfo &variable : variables()) {
    for (const std::string &group : groupsFor(variable)) {
      std::vector<std::vector<DeltaRecord> > byRun(runs.size());
      std::vector<double> allDeltas;
      std::vector<double> allRelatives;
      for (size_t r = 1; r < runs.size(); ++r) {
        byRun[r] = matchedDeltas(runs[0], runs[r], variable, group);
        std::vector<double> deltas, relatives;
        for (const DeltaRecord &item : byRun[r]) { deltas.push_back(item.delta); relatives.push_back(item.relative); }
        allDeltas.insert(allDeltas.end(), deltas.begin(), deltas.end());
        allRelatives.insert(allRelatives.end(), relatives.begin(), relatives.end());
        if (deltas.empty()) continue;
        const Summary deltaSummary = summarize(deltas);
        const Summary relativeSummary = summarize(relatives);
        summary << runs[0].label << "," << runs[r].label << "," << group << "," << variable.name
                << "," << deltas.size() << "," << std::setprecision(12) << deltaSummary.mean << ","
                << deltaSummary.rms << "," << deltaSummary.median << "," << deltaSummary.q16 << ","
                << deltaSummary.q84 << "," << relativeSummary.median << ","
                << (std::abs(relativeSummary.median) >= relativeThreshold) << "\n";

        std::vector<double> absoluteDeviations;
        for (double delta : deltas) absoluteDeviations.push_back(std::abs(delta - deltaSummary.median));
        const double mad = summarize(absoluteDeviations).median;
        for (const DeltaRecord &item : byRun[r]) {
          const double robustZ = mad > 1e-12 ? 0.67448975 * (item.delta - deltaSummary.median) / mad : 0.0;
          const bool relativeFlag = std::abs(item.relative) >= relativeThreshold;
          const bool robustFlag = mad > 1e-12 && std::abs(robustZ) >= robustZThreshold;
          if (!relativeFlag && !robustFlag) continue;
          const Record &record = *item.current;
          outliers << runs[0].label << "," << runs[r].label << "," << group << "," << variable.name
                   << "," << record.cellid << "," << record.layer << "," << record.chip << ","
                   << record.channel << "," << record.decision << "," << record.decisionName << ","
                   << record.state << ",\"" << record.stateNames << "\"," << record.*(variable.member)
                   << "," << item.delta << "," << item.relative << "," << robustZ << ","
                   << relativeFlag << "," << robustFlag << "\n";
        }
      }
      if (allDeltas.empty()) continue;
      std::pair<double, double> range = histogramRange(allDeltas);
      range.first = std::min(range.first, 0.0); range.second = std::max(range.second, 0.0);
      const std::string groupDir = plotDir + "/" + group;
      gSystem->mkdir(groupDir.c_str(), true);
      TCanvas canvas(("c_delta_" + group + "_" + variable.name).c_str(), "", 1000, 800);
      TLegend legend(0.65, 0.60, 0.92, 0.89);
      legend.SetHeader(("Reference: " + runs[0].label + " (same source)").c_str()); legend.SetBorderSize(0);
      if (runs.size() > 7) legend.SetNColumns(2);
      std::vector<TH1D *> histograms; double maximum = 0.0;
      for (size_t r = 1; r < runs.size(); ++r) {
        if (byRun[r].empty()) continue;
        TH1D *hist = new TH1D(("h_delta_" + group + "_" + variable.name + "_" + std::to_string(r)).c_str(),
                              "", 60, range.first, range.second);
        hist->SetDirectory(nullptr); hist->SetStats(false);
        for (const DeltaRecord &item : byRun[r]) hist->Fill(item.delta);
        hist->Scale(1.0 / hist->Integral(), "width"); hist->SetLineColor(runColor(r)); hist->SetLineWidth(2);
        maximum = std::max(maximum, hist->GetMaximum()); histograms.push_back(hist);
        legend.AddEntry(hist, (runs[r].label + " (N=" + std::to_string(byRun[r].size()) + ")").c_str(), "l");
      }
      histograms[0]->SetTitle((std::string(variable.title) + " change from first run: " + group +
                               ";#Delta " + variable.name + ";Density").c_str());
      histograms[0]->SetMaximum(1.18 * maximum); histograms[0]->Draw("HIST");
      for (size_t i = 1; i < histograms.size(); ++i) histograms[i]->Draw("HIST SAME");
      TLine zero(0.0, 0.0, 0.0, 1.18 * maximum); zero.SetLineStyle(2); zero.Draw(); legend.Draw();
      canvas.SaveAs((groupDir + "/" + variable.name + "_delta_overlay.pdf").c_str());
      canvas.SaveAs((groupDir + "/" + variable.name + "_delta_overlay.png").c_str());
      rootFile.cd(); for (TH1D *hist : histograms) { hist->Write(); delete hist; }
      if (std::string(variable.name) == "mpv" && !allRelatives.empty()) {
        std::pair<double, double> relativeRange = histogramRange(allRelatives);
        relativeRange.first = std::min(relativeRange.first, 0.0);
        relativeRange.second = std::max(relativeRange.second, 0.0);
        TCanvas relativeCanvas(("c_relative_delta_" + group + "_" + variable.name).c_str(), "", 1000, 800);
        TLegend relativeLegend(0.65, 0.60, 0.92, 0.89);
        relativeLegend.SetHeader(("Reference: " + runs[0].label + " (same source)").c_str());
        relativeLegend.SetBorderSize(0);
        if (runs.size() > 7) relativeLegend.SetNColumns(2);
        std::vector<TH1D *> relativeHistograms;
        double relativeMaximum = 0.0;
        for (size_t r = 1; r < runs.size(); ++r) {
          if (byRun[r].empty()) continue;
          TH1D *hist = new TH1D(("h_relative_delta_" + group + "_" + variable.name + "_" + std::to_string(r)).c_str(),
                                "", 60, relativeRange.first, relativeRange.second);
          hist->SetDirectory(nullptr);
          hist->SetStats(false);
          for (const DeltaRecord &item : byRun[r]) hist->Fill(item.relative);
          hist->Scale(1.0 / hist->Integral(), "width");
          hist->SetLineColor(runColor(r));
          hist->SetLineWidth(2);
          relativeMaximum = std::max(relativeMaximum, hist->GetMaximum());
          relativeHistograms.push_back(hist);
          relativeLegend.AddEntry(hist, (runs[r].label + " (N=" + std::to_string(byRun[r].size()) + ")").c_str(), "l");
        }
        relativeHistograms[0]->SetTitle((std::string(variable.title) + " relative change from first run: " + group +
                                         ";#Delta MPV / MPV;Density").c_str());
        relativeHistograms[0]->SetMaximum(1.18 * relativeMaximum);
        relativeHistograms[0]->Draw("HIST");
        for (size_t i = 1; i < relativeHistograms.size(); ++i) relativeHistograms[i]->Draw("HIST SAME");
        TLine relativeZero(0.0, 0.0, 0.0, 1.18 * relativeMaximum);
        relativeZero.SetLineStyle(2);
        relativeZero.Draw();
        relativeLegend.Draw();
        relativeCanvas.SaveAs((groupDir + "/" + variable.name + "_relative_delta_overlay.pdf").c_str());
        relativeCanvas.SaveAs((groupDir + "/" + variable.name + "_relative_delta_overlay.png").c_str());
        rootFile.cd();
        for (TH1D *hist : relativeHistograms) { hist->Write(); delete hist; }
      }
    }
  }
  rootFile.Close();
}

std::set<int> writeLayerMpvThreeSigmaOutliers(const std::vector<RunData> &runs,
                                              const std::string &outDir,
                                              int minEntries = 100,
                                              size_t minChannelsInLayer = 3) {
  const std::set<int> layers = collectLayers(runs);
  std::set<int> outlierCells;
  std::ofstream csv((outDir + "/LayerMPVThreeSigmaOutliers.csv").c_str());
  csv << "run_start,run_end,run_label,layer,cellid,chip,channel,entries,mpv,"
      << "layer_mean,layer_rms,sigma_distance,decision,decision_name,state,state_names\n";

  for (const RunData &run : runs) {
    for (int layer : layers) {
      std::vector<const Record *> records;
      std::vector<double> mpvs;
      for (const Record &record : run.records) {
        if (record.layer != layer) continue;
        if (record.entries < minEntries || !validValue(record.mpv)) continue;
        records.push_back(&record);
        mpvs.push_back(record.mpv);
      }
      if (mpvs.size() < minChannelsInLayer) continue;
      const Summary summary = summarize(mpvs);
      if (summary.rms <= 1e-12) continue;

      for (const Record *record : records) {
        const double sigmaDistance = (record->mpv - summary.mean) / summary.rms;
        if (std::abs(sigmaDistance) < 3.0) continue;
        outlierCells.insert(record->cellid);
        csv << run.runStart << "," << run.runEnd << "," << csvQuote(run.label) << ","
            << layer << "," << record->cellid << "," << record->chip << ","
            << record->channel << "," << record->entries << ","
            << std::setprecision(12) << record->mpv << "," << summary.mean << ","
            << summary.rms << "," << sigmaDistance << "," << record->decision << ","
            << csvQuote(record->decisionName) << "," << record->state << ","
            << csvQuote(record->stateNames) << "\n";
      }
    }
  }
  return outlierCells;
}

void drawLayerMpvThreeSigmaOutlierTrends(const std::vector<RunData> &runs,
                                         const std::string &outDir,
                                         const std::set<int> &outlierCells,
                                         int minEntries = 100) {
  if (outlierCells.empty()) return;
  const std::string plotDir = outDir + "/LayerMPVThreeSigmaOutlierTrends";
  gSystem->mkdir(plotDir.c_str(), true);
  TFile rootFile((plotDir + "/layer_mpv_three_sigma_outlier_trends.root").c_str(), "RECREATE");

  for (int cellid : outlierCells) {
    const Record *metadata = nullptr;
    for (const RunData &run : runs) {
      const auto index = run.indexByCell.find(cellid);
      if (index != run.indexByCell.end()) {
        metadata = &run.records[index->second];
        break;
      }
    }
    if (!metadata) continue;

    std::vector<double> x, y, ex, ey;
    std::vector<double> missingX, missingYPlaceholder;
    for (size_t r = 0; r < runs.size(); ++r) {
      const auto index = runs[r].indexByCell.find(cellid);
      if (index == runs[r].indexByCell.end()) {
        missingX.push_back(r + 1.0);
        continue;
      }
      const Record &record = runs[r].records[index->second];
      if (!validValue(record.mpv) || record.entries < minEntries) {
        missingX.push_back(r + 1.0);
        continue;
      }
      x.push_back(r + 1.0);
      y.push_back(record.mpv);
      ex.push_back(0.0);
      ey.push_back(validValue(record.mpvError) ? record.mpvError : 0.0);
    }
    if (y.empty()) continue;

    std::vector<double> yRangeValues = y;
    for (size_t i = 0; i < y.size(); ++i) {
      yRangeValues.push_back(y[i] - ey[i]);
      yRangeValues.push_back(y[i] + ey[i]);
    }
    std::pair<double, double> range = histogramRange(yRangeValues);
    const double span = std::max(1.0, range.second - range.first);
    const double missingY = range.first - 0.08 * span;
    const double frameYMin = range.first - 0.18 * span;
    const double frameYMax = range.second + 0.12 * span;
    missingYPlaceholder.assign(missingX.size(), missingY);

    const std::string channelKey =
        Form("cell%d_layer%02d_chip%02d_channel%02d",
             cellid, metadata->layer, metadata->chip, metadata->channel);
    TCanvas canvas(("c_layer_mpv_three_sigma_outlier_trend_" + channelKey).c_str(), "", 1200, 750);
    canvas.SetBottomMargin(0.22);
    canvas.SetLeftMargin(0.10);
    canvas.SetGrid();

    TH1D frame(("frame_layer_mpv_three_sigma_outlier_trend_" + channelKey).c_str(),
               Form("MPV trend for 3#sigma layer outlier: cellid=%d, layer=%d chip=%d ch=%d;Run;MPV",
                    cellid, metadata->layer, metadata->chip, metadata->channel),
               runs.size(), 0.5, runs.size() + 0.5);
    frame.SetMinimum(frameYMin);
    frame.SetMaximum(frameYMax);
    frame.SetStats(false);
    for (size_t r = 0; r < runs.size(); ++r) {
      frame.GetXaxis()->SetBinLabel(r + 1, runs[r].label.c_str());
    }
    frame.GetXaxis()->LabelsOption("v");
    frame.GetYaxis()->SetTitleOffset(1.15);
    frame.Draw();

    TGraphErrors validGraph(x.size(), &x[0], &y[0], &ex[0], &ey[0]);
    validGraph.SetName(("g_layer_mpv_three_sigma_outlier_trend_" + channelKey).c_str());
    validGraph.SetMarkerStyle(20);
    validGraph.SetMarkerSize(1.0);
    validGraph.SetMarkerColor(kBlue + 1);
    validGraph.SetLineColor(kBlue + 1);
    validGraph.SetLineWidth(2);
    validGraph.Draw("PL SAME");

    TGraph missingGraph;
    missingGraph.SetName(("g_missing_layer_mpv_three_sigma_outlier_trend_" + channelKey).c_str());
    missingGraph.SetMarkerStyle(5);
    missingGraph.SetMarkerSize(1.4);
    missingGraph.SetMarkerColor(kRed + 1);
    missingGraph.SetLineColor(kRed + 1);
    for (size_t i = 0; i < missingX.size(); ++i) {
      missingGraph.SetPoint(i, missingX[i], missingYPlaceholder[i]);
    }
    if (!missingX.empty()) missingGraph.Draw("P SAME");

    TLine missingLine(0.5, missingY, runs.size() + 0.5, missingY);
    missingLine.SetLineColor(kRed + 1);
    missingLine.SetLineStyle(2);
    if (!missingX.empty()) missingLine.Draw();

    TLegend legend(0.58, 0.70, 0.94, 0.89);
    legend.SetBorderSize(0);
    legend.AddEntry(&validGraph, Form("valid MPV (entries >= %d)", minEntries), "pl");
    if (!missingX.empty()) {
      legend.AddEntry(&missingGraph, Form("missing / invalid / entries < %d", minEntries), "p");
    }
    legend.Draw();

    const std::string baseName = plotDir + "/" + channelKey + "_mpv_trend";
    canvas.SaveAs((baseName + ".pdf").c_str());
    canvas.SaveAs((baseName + ".png").c_str());
    rootFile.cd();
    canvas.Write();
    frame.Write();
    validGraph.Write();
    if (!missingX.empty()) missingGraph.Write();
  }
  rootFile.Close();
}

struct MpvRunInstabilityChannel {
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;
  size_t validRuns = 0;
  size_t missingRuns = 0;
  double meanMpv = 0.0;
  double rmsMpv = 0.0;
  double meanMpvError = 0.0;
  double rmsOverMeanError = 0.0;
  double minMpv = 0.0;
  double maxMpv = 0.0;
};

std::vector<MpvRunInstabilityChannel> writeMpvRunInstabilitySummary(
    const std::vector<RunData> &runs,
    const std::string &outDir,
    double minRmsOverMeanError = 3.0,
    int minEntries = 100,
    size_t minValidRuns = 3) {
  std::set<int> cellids;
  for (const RunData &run : runs) {
    for (const Record &record : run.records) cellids.insert(record.cellid);
  }

  std::vector<MpvRunInstabilityChannel> flagged;
  for (int cellid : cellids) {
    const Record *metadata = nullptr;
    std::vector<double> mpvs;
    std::vector<double> mpvErrors;
    size_t missingRuns = 0;
    for (const RunData &run : runs) {
      const auto index = run.indexByCell.find(cellid);
      if (index == run.indexByCell.end()) {
        ++missingRuns;
        continue;
      }
      const Record &record = run.records[index->second];
      if (!metadata) metadata = &record;
      if (!validValue(record.mpv) || record.entries < minEntries ||
          !std::isfinite(record.mpvError) || record.mpvError <= 0.0) {
        ++missingRuns;
        continue;
      }
      mpvs.push_back(record.mpv);
      mpvErrors.push_back(record.mpvError);
    }
    if (!metadata || mpvs.size() < minValidRuns) continue;
    const Summary mpvSummary = summarize(mpvs);
    const Summary errorSummary = summarize(mpvErrors);
    if (errorSummary.mean <= 1e-12) continue;
    const double ratio = mpvSummary.rms / errorSummary.mean;
    if (ratio < minRmsOverMeanError) continue;

    MpvRunInstabilityChannel item;
    item.cellid = cellid;
    item.layer = metadata->layer;
    item.chip = metadata->chip;
    item.channel = metadata->channel;
    item.validRuns = mpvs.size();
    item.missingRuns = missingRuns;
    item.meanMpv = mpvSummary.mean;
    item.rmsMpv = mpvSummary.rms;
    item.meanMpvError = errorSummary.mean;
    item.rmsOverMeanError = ratio;
    item.minMpv = mpvSummary.minimum;
    item.maxMpv = mpvSummary.maximum;
    flagged.push_back(item);
  }

  std::sort(flagged.begin(), flagged.end(),
            [](const MpvRunInstabilityChannel &a, const MpvRunInstabilityChannel &b) {
              return a.rmsOverMeanError > b.rmsOverMeanError;
            });

  std::ofstream csv((outDir + "/MPVRunInstabilitySummary.csv").c_str());
  csv << "cellid,layer,chip,channel,valid_runs,missing_or_invalid_runs,mean_mpv,"
      << "rms_mpv,mean_mpv_error,rms_over_mean_error,min_mpv,max_mpv\n";
  for (const MpvRunInstabilityChannel &item : flagged) {
    csv << item.cellid << "," << item.layer << "," << item.chip << ","
        << item.channel << "," << item.validRuns << "," << item.missingRuns << ","
        << std::setprecision(12) << item.meanMpv << "," << item.rmsMpv << ","
        << item.meanMpvError << "," << item.rmsOverMeanError << ","
        << item.minMpv << "," << item.maxMpv << "\n";
  }
  return flagged;
}

void drawMpvRunInstabilityTrends(const std::vector<RunData> &runs,
                                 const std::string &outDir,
                                 const std::vector<MpvRunInstabilityChannel> &channels,
                                 int minEntries = 100) {
  if (channels.empty()) return;
  const std::string plotDir = outDir + "/MPVRunInstabilityTrends";
  gSystem->mkdir(plotDir.c_str(), true);
  TFile rootFile((plotDir + "/mpv_run_instability_trends.root").c_str(), "RECREATE");

  for (const MpvRunInstabilityChannel &channel : channels) {
    std::vector<double> x, y, ex, ey;
    std::vector<double> missingX, missingYPlaceholder;
    for (size_t r = 0; r < runs.size(); ++r) {
      const auto index = runs[r].indexByCell.find(channel.cellid);
      if (index == runs[r].indexByCell.end()) {
        missingX.push_back(r + 1.0);
        continue;
      }
      const Record &record = runs[r].records[index->second];
      if (!validValue(record.mpv) || record.entries < minEntries ||
          !std::isfinite(record.mpvError) || record.mpvError <= 0.0) {
        missingX.push_back(r + 1.0);
        continue;
      }
      x.push_back(r + 1.0);
      y.push_back(record.mpv);
      ex.push_back(0.0);
      ey.push_back(record.mpvError);
    }
    if (y.empty()) continue;

    std::vector<double> yRangeValues = y;
    for (size_t i = 0; i < y.size(); ++i) {
      yRangeValues.push_back(y[i] - ey[i]);
      yRangeValues.push_back(y[i] + ey[i]);
    }
    std::pair<double, double> range = histogramRange(yRangeValues);
    const double span = std::max(1.0, range.second - range.first);
    const double missingY = range.first - 0.08 * span;
    const double frameYMin = range.first - 0.18 * span;
    const double frameYMax = range.second + 0.12 * span;
    missingYPlaceholder.assign(missingX.size(), missingY);

    const std::string channelKey =
        Form("cell%d_layer%02d_chip%02d_channel%02d",
             channel.cellid, channel.layer, channel.chip, channel.channel);
    TCanvas canvas(("c_mpv_run_instability_trend_" + channelKey).c_str(), "", 1200, 750);
    canvas.SetBottomMargin(0.22);
    canvas.SetLeftMargin(0.10);
    canvas.SetGrid();

    TH1D frame(("frame_mpv_run_instability_trend_" + channelKey).c_str(),
               Form("MPV run instability: cellid=%d, layer=%d chip=%d ch=%d, RMS/<err>=%.2f;Run;MPV",
                    channel.cellid, channel.layer, channel.chip, channel.channel,
                    channel.rmsOverMeanError),
               runs.size(), 0.5, runs.size() + 0.5);
    frame.SetMinimum(frameYMin);
    frame.SetMaximum(frameYMax);
    frame.SetStats(false);
    for (size_t r = 0; r < runs.size(); ++r) {
      frame.GetXaxis()->SetBinLabel(r + 1, runs[r].label.c_str());
    }
    frame.GetXaxis()->LabelsOption("v");
    frame.GetYaxis()->SetTitleOffset(1.15);
    frame.Draw();

    TGraphErrors validGraph(x.size(), &x[0], &y[0], &ex[0], &ey[0]);
    validGraph.SetName(("g_mpv_run_instability_trend_" + channelKey).c_str());
    validGraph.SetMarkerStyle(20);
    validGraph.SetMarkerSize(1.0);
    validGraph.SetMarkerColor(kBlue + 1);
    validGraph.SetLineColor(kBlue + 1);
    validGraph.SetLineWidth(2);
    validGraph.Draw("PL SAME");

    TGraph missingGraph;
    missingGraph.SetName(("g_missing_mpv_run_instability_trend_" + channelKey).c_str());
    missingGraph.SetMarkerStyle(5);
    missingGraph.SetMarkerSize(1.4);
    missingGraph.SetMarkerColor(kRed + 1);
    missingGraph.SetLineColor(kRed + 1);
    for (size_t i = 0; i < missingX.size(); ++i) {
      missingGraph.SetPoint(i, missingX[i], missingYPlaceholder[i]);
    }
    if (!missingX.empty()) missingGraph.Draw("P SAME");

    TLine missingLine(0.5, missingY, runs.size() + 0.5, missingY);
    missingLine.SetLineColor(kRed + 1);
    missingLine.SetLineStyle(2);
    if (!missingX.empty()) missingLine.Draw();

    TLegend legend(0.55, 0.68, 0.94, 0.89);
    legend.SetBorderSize(0);
    legend.AddEntry(&validGraph,
                    Form("MPV #pm error, RMS=%.3g, <err>=%.3g",
                         channel.rmsMpv, channel.meanMpvError),
                    "pl");
    if (!missingX.empty()) {
      legend.AddEntry(&missingGraph,
                      Form("missing / invalid / entries < %d / no fit error", minEntries),
                      "p");
    }
    legend.Draw();

    const std::string baseName = plotDir + "/" + channelKey + "_mpv_run_instability";
    canvas.SaveAs((baseName + ".pdf").c_str());
    canvas.SaveAs((baseName + ".png").c_str());
    rootFile.cd();
    canvas.Write();
    frame.Write();
    validGraph.Write();
    if (!missingX.empty()) missingGraph.Write();
  }
  rootFile.Close();
}

struct MpvDataLossChannel {
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;
  size_t validRuns = 0;
  size_t missingRuns = 0;
  std::string missingRunLabels;
};

std::vector<MpvDataLossChannel> writeMpvDataLossSummary(
    const std::vector<RunData> &runs,
    const std::string &outDir,
    int minEntries = 100,
    bool excludeLastRun = true) {
  const size_t runCount =
      excludeLastRun && runs.size() > 1 ? runs.size() - 1 : runs.size();
  std::set<int> cellids;
  for (size_t r = 0; r < runCount; ++r) {
    for (const Record &record : runs[r].records) cellids.insert(record.cellid);
  }

  std::vector<MpvDataLossChannel> flagged;
  for (int cellid : cellids) {
    const Record *metadata = nullptr;
    size_t validRuns = 0;
    std::vector<std::string> missingLabels;
    for (size_t r = 0; r < runCount; ++r) {
      const auto index = runs[r].indexByCell.find(cellid);
      if (index == runs[r].indexByCell.end()) {
        missingLabels.push_back(runs[r].label);
        continue;
      }
      const Record &record = runs[r].records[index->second];
      if (!metadata) metadata = &record;
      if (!validValue(record.mpv) || record.entries < minEntries) {
        missingLabels.push_back(runs[r].label);
        continue;
      }
      ++validRuns;
    }
    if (!metadata || missingLabels.empty()) continue;

    std::stringstream missingStream;
    for (size_t i = 0; i < missingLabels.size(); ++i) {
      if (i) missingStream << "|";
      missingStream << missingLabels[i];
    }

    MpvDataLossChannel item;
    item.cellid = cellid;
    item.layer = metadata->layer;
    item.chip = metadata->chip;
    item.channel = metadata->channel;
    item.validRuns = validRuns;
    item.missingRuns = missingLabels.size();
    item.missingRunLabels = missingStream.str();
    flagged.push_back(item);
  }

  std::sort(flagged.begin(), flagged.end(),
            [](const MpvDataLossChannel &a, const MpvDataLossChannel &b) {
              if (a.missingRuns != b.missingRuns) return a.missingRuns > b.missingRuns;
              if (a.layer != b.layer) return a.layer < b.layer;
              if (a.chip != b.chip) return a.chip < b.chip;
              return a.channel < b.channel;
            });

  std::ofstream csv((outDir + "/MPVDataLossChannels.csv").c_str());
  csv << "cellid,layer,chip,channel,valid_runs,missing_or_invalid_runs,"
      << "considered_runs,excluded_last_run,excluded_last_run_label,missing_run_labels\n";
  for (const MpvDataLossChannel &item : flagged) {
    csv << item.cellid << "," << item.layer << "," << item.chip << ","
        << item.channel << "," << item.validRuns << "," << item.missingRuns << ","
        << runCount << "," << excludeLastRun << ",";
    if (excludeLastRun && runs.size() > runCount) csv << csvQuote(runs.back().label);
    csv << "," << csvQuote(item.missingRunLabels) << "\n";
  }
  return flagged;
}

void drawMpvDataLossTrends(const std::vector<RunData> &runs,
                           const std::string &outDir,
                           const std::vector<MpvDataLossChannel> &channels,
                           int minEntries = 100,
                           bool excludeLastRun = true) {
  if (channels.empty()) return;
  const size_t runCount =
      excludeLastRun && runs.size() > 1 ? runs.size() - 1 : runs.size();
  if (!runCount) return;

  const std::string plotDir = outDir + "/MPVDataLossTrends";
  gSystem->mkdir(plotDir.c_str(), true);
  TFile rootFile((plotDir + "/mpv_data_loss_trends.root").c_str(), "RECREATE");

  for (const MpvDataLossChannel &channel : channels) {
    std::vector<double> x, y, ex, ey;
    std::vector<double> missingX, missingYPlaceholder;
    for (size_t r = 0; r < runCount; ++r) {
      const auto index = runs[r].indexByCell.find(channel.cellid);
      if (index == runs[r].indexByCell.end()) {
        missingX.push_back(r + 1.0);
        continue;
      }
      const Record &record = runs[r].records[index->second];
      if (!validValue(record.mpv) || record.entries < minEntries) {
        missingX.push_back(r + 1.0);
        continue;
      }
      x.push_back(r + 1.0);
      y.push_back(record.mpv);
      ex.push_back(0.0);
      ey.push_back(validValue(record.mpvError) ? record.mpvError : 0.0);
    }
    if (y.empty()) continue;

    std::vector<double> yRangeValues = y;
    for (size_t i = 0; i < y.size(); ++i) {
      yRangeValues.push_back(y[i] - ey[i]);
      yRangeValues.push_back(y[i] + ey[i]);
    }
    std::pair<double, double> range = histogramRange(yRangeValues);
    const double span = std::max(1.0, range.second - range.first);
    const double missingY = range.first - 0.08 * span;
    const double frameYMin = range.first - 0.18 * span;
    const double frameYMax = range.second + 0.12 * span;
    missingYPlaceholder.assign(missingX.size(), missingY);

    const std::string channelKey =
        Form("cell%d_layer%02d_chip%02d_channel%02d",
             channel.cellid, channel.layer, channel.chip, channel.channel);
    TCanvas canvas(("c_mpv_data_loss_trend_" + channelKey).c_str(), "", 1200, 750);
    canvas.SetBottomMargin(0.22);
    canvas.SetLeftMargin(0.10);
    canvas.SetGrid();

    TH1D frame(("frame_mpv_data_loss_trend_" + channelKey).c_str(),
               Form("MPV data-loss trend: cellid=%d, layer=%d chip=%d ch=%d;Run;MPV",
                    channel.cellid, channel.layer, channel.chip, channel.channel),
               runCount, 0.5, runCount + 0.5);
    frame.SetMinimum(frameYMin);
    frame.SetMaximum(frameYMax);
    frame.SetStats(false);
    for (size_t r = 0; r < runCount; ++r) {
      frame.GetXaxis()->SetBinLabel(r + 1, runs[r].label.c_str());
    }
    frame.GetXaxis()->LabelsOption("v");
    frame.GetYaxis()->SetTitleOffset(1.15);
    frame.Draw();

    TGraphErrors validGraph(x.size(), &x[0], &y[0], &ex[0], &ey[0]);
    validGraph.SetName(("g_mpv_data_loss_trend_" + channelKey).c_str());
    validGraph.SetMarkerStyle(20);
    validGraph.SetMarkerSize(1.0);
    validGraph.SetMarkerColor(kBlue + 1);
    validGraph.SetLineColor(kBlue + 1);
    validGraph.SetLineWidth(2);
    validGraph.Draw("PL SAME");

    TGraph missingGraph;
    missingGraph.SetName(("g_missing_mpv_data_loss_trend_" + channelKey).c_str());
    missingGraph.SetMarkerStyle(5);
    missingGraph.SetMarkerSize(1.4);
    missingGraph.SetMarkerColor(kRed + 1);
    missingGraph.SetLineColor(kRed + 1);
    for (size_t i = 0; i < missingX.size(); ++i) {
      missingGraph.SetPoint(i, missingX[i], missingYPlaceholder[i]);
    }
    missingGraph.Draw("P SAME");

    TLine missingLine(0.5, missingY, runCount + 0.5, missingY);
    missingLine.SetLineColor(kRed + 1);
    missingLine.SetLineStyle(2);
    missingLine.Draw();

    TLegend legend(0.55, 0.68, 0.94, 0.89);
    legend.SetBorderSize(0);
    legend.AddEntry(&validGraph, Form("valid MPV (entries >= %d)", minEntries), "pl");
    legend.AddEntry(&missingGraph,
                    Form("data loss: missing / invalid / entries < %d", minEntries),
                    "p");
    if (excludeLastRun && runs.size() > runCount) {
      legend.SetHeader(("Excluded last run: " + runs.back().label).c_str());
    }
    legend.Draw();

    const std::string baseName = plotDir + "/" + channelKey + "_mpv_data_loss";
    canvas.SaveAs((baseName + ".pdf").c_str());
    canvas.SaveAs((baseName + ".png").c_str());
    rootFile.cd();
    canvas.Write();
    frame.Write();
    validGraph.Write();
    missingGraph.Write();
  }
  rootFile.Close();
}

TH1 *loadAdcHistogram(const std::string &filePath, const Record &record,
                      const char *cloneName) {
  TFile file(filePath.c_str(), "READ");
  if (file.IsZombie()) {
    return nullptr;
  }

  const int layer = record.layer >= 0 ? record.layer : record.cellid / 100000;
  const int chip = record.chip >= 0 ? record.chip : (record.cellid / 10000) % 10;
  std::vector<std::string> candidates;
  candidates.push_back(Form("MIP/Layer%02d/Chip%d/hMIP_%d", layer, chip, record.cellid));
  candidates.push_back(Form("MIP/Layer%d/Chip%d/hMIP_%d", layer, chip, record.cellid));
  candidates.push_back(Form("MIP/Layer%02d/Chip%02d/hMIP_%d", layer, chip, record.cellid));
  candidates.push_back(Form("hMIP_%d", record.cellid));

  TH1 *hist = nullptr;
  for (const std::string &name : candidates) {
    hist = dynamic_cast<TH1 *>(file.Get(name.c_str()));
    if (hist) break;
  }

  TH1 *clone = nullptr;
  if (hist) {
    clone = dynamic_cast<TH1 *>(hist->Clone(cloneName));
    if (clone) {
      clone->SetDirectory(nullptr);
      clone->Rebin(4);
    }
  }
  file.Close();
  return clone;
}

struct LargeValueShift {
  int cellid = -1;
  size_t runIndex = 0;
  double referenceValue = -1.0;
  double currentValue = -1.0;
  double delta = 0.0;
};

void drawLargeValueShiftAdcHistograms(const std::vector<RunData> &runs,
                                      const std::string &outDir,
                                      const std::string &histBaseDir,
                                      const std::string &histFileName,
                                      double Record::*valueMember,
                                      const std::string &valueName,
                                      const std::string &valueTitle,
                                      double absDeltaThreshold,
                                      const std::string &plotSubdir,
                                      const std::string &csvFileName,
                                      const std::string &rootFileName) {
  if (runs.size() < 2 || absDeltaThreshold <= 0.0) return;

  std::map<int, std::vector<LargeValueShift> > shiftsByCell;
  for (size_t r = 1; r < runs.size(); ++r) {
    for (const Record &record : runs[r].records) {
      if (!validValue(record.*valueMember)) continue;
      const auto refIndex = runs[0].indexByCell.find(record.cellid);
      if (refIndex == runs[0].indexByCell.end()) continue;
      const Record &reference = runs[0].records[refIndex->second];
      if (!validValue(reference.*valueMember)) continue;
      const double delta = record.*valueMember - reference.*valueMember;
      if (std::abs(delta) <= absDeltaThreshold) continue;
      LargeValueShift shift;
      shift.cellid = record.cellid;
      shift.runIndex = r;
      shift.referenceValue = reference.*valueMember;
      shift.currentValue = record.*valueMember;
      shift.delta = delta;
      shiftsByCell[record.cellid].push_back(shift);
    }
  }
  if (shiftsByCell.empty()) return;

  const std::string plotDir = outDir + "/DeltaFromFirstRun/" + plotSubdir;
  gSystem->mkdir(plotDir.c_str(), true);
  TFile rootFile((plotDir + "/" + rootFileName).c_str(), "RECREATE");
  std::ofstream csv((outDir + "/" + csvFileName).c_str());
  csv << "reference_run,run,cellid,layer,chip,channel,reference_" << valueName
      << ",current_" << valueName << ",delta_" << valueName << ",hist_file\n";
  std::set<std::string> openedChannelPdfs;

  for (const auto &cellItem : shiftsByCell) {
    const int cellid = cellItem.first;
    const Record *reference = nullptr;
    const auto refIndex = runs[0].indexByCell.find(cellid);
    if (refIndex != runs[0].indexByCell.end()) reference = &runs[0].records[refIndex->second];
    if (!reference) continue;

    std::vector<std::pair<TH1 *, size_t> > histograms;

    for (size_t r = 0; r < runs.size(); ++r) {
      const auto index = runs[r].indexByCell.find(cellid);
      if (index == runs[r].indexByCell.end()) continue;
      const Record &record = runs[r].records[index->second];
      const std::string histPath = joinPath(joinPath(histBaseDir, runs[r].label), histFileName);
      if (!fileExists(histPath)) continue;
      TH1 *hist = loadAdcHistogram(histPath, record,
                                   Form("h_adc_cell%d_run%zu", cellid, r));
      if (!hist) continue;
      hist->SetStats(false);
      hist->SetLineColor(runColor(r));
      hist->SetMarkerColor(runColor(r));
      hist->SetLineWidth(2);
      hist->SetTitle(Form("ADC histogram for |#Delta %s| > %.1f;ADC;Entries",
                          valueTitle.c_str(), absDeltaThreshold));
      hist->Rebin(4);
      hist->GetXaxis()->SetRangeUser(-48, 1000.0);
      if (hist->Integral() > 0.0) hist->Scale(1.0 / hist->Integral());
      histograms.push_back(std::make_pair(hist, r));
      for (const LargeValueShift &shift : cellItem.second) {
        if (shift.runIndex != r) continue;
        csv << runs[0].label << "," << runs[r].label << "," << cellid << ","
            << record.layer << "," << record.chip << "," << record.channel << ","
            << std::setprecision(12) << shift.referenceValue << "," << shift.currentValue
            << "," << shift.delta << "," << histPath << "\n";
      }
    }
    if (histograms.empty()) continue;

    const std::string channelKey =
        Form("layer_%02d_chip_%02d_channel_%02d",
             reference->layer, reference->chip, reference->channel);
    const std::string channelPdf = plotDir + "/" + channelKey + "_" + valueName + "_adc_by_run.pdf";
    const bool firstChannelPage = openedChannelPdfs.insert(channelKey).second;
    if (firstChannelPage) {
      TCanvas opener(("c_open_adc_large_mpv_shift_" + channelKey).c_str(), "", 1000, 800);
      opener.Print((channelPdf + "[").c_str());
    }

    rootFile.cd();
    double maximum = 0.0;
    for (const auto &histItem : histograms) maximum = std::max(maximum, histItem.first->GetMaximum());
    TCanvas overlayCanvas(("c_adc_overlay_" + channelKey + "_" + valueName).c_str(), "", 1000, 800);
    TLegend overlayLegend(0.55, 0.56, 0.92, 0.89);
    overlayLegend.SetBorderSize(0);
    overlayLegend.SetHeader(Form("cellid=%d, layer=%d chip=%d ch=%d",
                                 cellid, reference->layer, reference->chip, reference->channel));
    for (size_t i = 0; i < histograms.size(); ++i) {
      TH1 *hist = histograms[i].first;
      const size_t r = histograms[i].second;
      const auto index = runs[r].indexByCell.find(cellid);
      if (index == runs[r].indexByCell.end()) continue;
      const Record &record = runs[r].records[index->second];
      const double delta = validValue(record.*valueMember) && validValue(reference->*valueMember)
                               ? record.*valueMember - reference->*valueMember
                               : 0.0;
      overlayLegend.AddEntry(hist,
                             Form("%s %s=%.1f #Delta=%.1f N=%.0f",
                                  runs[r].label.c_str(), valueTitle.c_str(),
                                  record.*valueMember, delta, hist->GetEntries()),
                             "l");
      if (i == 0) {
        hist->SetTitle(Form("ADC histograms for |#Delta %s| > %.1f;ADC;Entries",
                            valueTitle.c_str(), absDeltaThreshold));
        hist->SetMaximum(maximum * 1.2);
        hist->SetMinimum(0);
        hist->Draw("HIST");
      } else {
        hist->Draw("HIST SAME");
      }
    }
    overlayLegend.Draw();
    const std::string overlayBaseName =
        plotDir + "/" + channelKey + "_" + valueName + "_adc_overlay";
    overlayCanvas.SaveAs((overlayBaseName + ".pdf").c_str());
    overlayCanvas.SaveAs((overlayBaseName + ".png").c_str());
    overlayCanvas.Write();

    for (const auto &histItem : histograms) {
      TH1 *hist = histItem.first;
      const size_t r = histItem.second;
      const auto index = runs[r].indexByCell.find(cellid);
      if (index == runs[r].indexByCell.end()) continue;
      const Record &record = runs[r].records[index->second];
      const double delta = validValue(record.*valueMember) && validValue(reference->*valueMember)
                               ? record.*valueMember - reference->*valueMember
                               : 0.0;
      TCanvas canvas(Form("c_adc_large_mpv_shift_cell%d_run%zu", cellid, r), "", 1000, 800);
      TLegend legend(0.55, 0.70, 0.92, 0.89);
      legend.SetBorderSize(0);
      legend.SetHeader(Form("cellid=%d, layer=%d chip=%d ch=%d",
                            cellid, reference->layer, reference->chip, reference->channel));
      legend.AddEntry(hist,
                      Form("%s %s=%.1f #Delta=%.1f N=%.0f",
                           runs[r].label.c_str(), valueTitle.c_str(),
                           record.*valueMember, delta, hist->GetEntries()),
                      "l");
      hist->SetMaximum(hist->GetMaximum() * 1.2);
      hist->SetMinimum(0);
      hist->Draw("HIST");
      legend.Draw();
      canvas.Print(channelPdf.c_str());
      canvas.Write();
      hist->Write();
      delete hist;
    }
  }
  for (const std::string &channelKey : openedChannelPdfs) {
    TCanvas closer(("c_close_adc_large_mpv_shift_" + channelKey).c_str(), "", 1000, 800);
    const std::string channelPdf = plotDir + "/" + channelKey + "_" + valueName + "_adc_by_run.pdf";
    closer.Print((channelPdf + "]").c_str());
  }
  rootFile.Close();
}

} // namespace final_mip_monitor

void monitor_final_mip_calibration(
    const char *inputFile = "FinalMIPCalibration.root",
    const char *outDir = "final_mip_monitoring",
    const char *treeName = "mip_calibration",
    double relativeShiftThreshold = 0.10,
    double robustZThreshold = 5.0,
    double absDeltaMpvThreshold = 20.0,
    const char *histBaseDir = "../",
    const char *histFileName = "mip_neighborcheck_nofit.root") {
  using namespace final_mip_monitor;
  gSystem->mkdir(outDir, true);
  gStyle->SetOptStat(0);
  std::vector<RunData> runs;
  if (!loadRuns(inputFile, treeName, runs)) return;
  const std::string resolvedHistBaseDir = resolveHistBaseDir(
      inputFile, histBaseDir, runs,
      histFileName ? histFileName : "mip_neighborcheck_nofit.root");
  std::cout << "Loaded " << runs.size() << " runs from " << inputFile << std::endl;
  writeSummary(runs, outDir);
  writeStateSummary(runs, outDir);
  writeDecisionLayerSummary(runs, outDir);
  drawTrends(runs, outDir);
  drawOverlays(runs, outDir);
  writeAndDrawDeltas(runs, outDir, relativeShiftThreshold, robustZThreshold);
  const std::set<int> layerMpvThreeSigmaOutlierCells =
      writeLayerMpvThreeSigmaOutliers(runs, outDir);
  drawLayerMpvThreeSigmaOutlierTrends(runs, outDir, layerMpvThreeSigmaOutlierCells);
  const std::vector<MpvRunInstabilityChannel> mpvRunInstabilityChannels =
      writeMpvRunInstabilitySummary(runs, outDir);
  drawMpvRunInstabilityTrends(runs, outDir, mpvRunInstabilityChannels);
  const bool excludeLastRunForMpvDataLoss = true;
  const std::vector<MpvDataLossChannel> mpvDataLossChannels =
      writeMpvDataLossSummary(runs, outDir, 100, excludeLastRunForMpvDataLoss);
  drawMpvDataLossTrends(runs, outDir, mpvDataLossChannels, 100,
                        excludeLastRunForMpvDataLoss);
  // drawLargeValueShiftAdcHistograms(runs, outDir, resolvedHistBaseDir,
  //                                  histFileName ? histFileName : "mip_neighborcheck_nofit.root",
  //                                  &Record::mpv, "mpv", "MPV",
  //                                  absDeltaMpvThreshold,
  //                                  "ADCForLargeMPVShift",
  //                                  "LargeMPVShiftChannels_absDelta.csv",
  //                                  "large_mpv_shift_adc_histograms.root");
  // drawLargeValueShiftAdcHistograms(runs, outDir, resolvedHistBaseDir,
  //                                  histFileName ? histFileName : "mip_neighborcheck_nofit.root",
  //                                  &Record::directThreshold, "direct_threshold", "direct threshold",
  //                                  10.0,
  //                                  "ADCForLargeDirectThresholdShift",
  //                                  "LargeDirectThresholdShiftChannels_absDelta.csv",
  //                                  "large_direct_threshold_shift_adc_histograms.root");
  std::cout << "Monitoring outputs written to " << outDir << std::endl;
  std::cout << "  FinalValueSummary_byRun.csv\n"
            << "  FinalValueSummary_byRunState.csv\n"
            << "  FinalDecisionFraction_byRunLayer.csv\n"
            << "  FinalDecisionFraction_byRunLayer.root\n"
            << "  DeltaSummary_fromFirstRun.csv\n"
            << "  LargeShiftChannels.csv\n"
            << "  LayerMPVThreeSigmaOutliers.csv\n"
            << "  LayerMPVThreeSigmaOutlierTrends/\n"
            << "  MPVRunInstabilitySummary.csv\n"
            << "  MPVRunInstabilityTrends/\n"
            << "  MPVDataLossChannels.csv\n"
            << "  MPVDataLossTrends/\n"
            << "  LargeMPVShiftChannels_absDelta.csv\n"
            << "  LargeDirectThresholdShiftChannels_absDelta.csv\n"
            << "  RunLayerSuccessRate/, DecisionLayerSummaryStacks/, Trends/, OverlayDistributions/, LayerOverlayDistributions/, DeltaFromFirstRun/" << std::endl;
}

#ifndef __CLING__
int main(int argc, char **argv) {
  const char *input = argc > 1 ? argv[1] : "FinalMIPCalibration.root";
  const char *output = argc > 2 ? argv[2] : "final_mip_monitoring";
  const double relativeThreshold = argc > 3 ? std::atof(argv[3]) : 0.10;
  const double robustZThreshold = argc > 4 ? std::atof(argv[4]) : 5.0;
  const double absDeltaMpvThreshold = argc > 5 ? std::atof(argv[5]) : 20.0;
  const char *histBaseDir = argc > 6 ? argv[6] : "";
  const char *histFileName = argc > 7 ? argv[7] : "mip_neighborcheck_nofit.root";
  monitor_final_mip_calibration(input, output, "mip_calibration",
                                relativeThreshold, robustZThreshold,
                                absDeltaMpvThreshold, histBaseDir, histFileName);
  return 0;
}
#endif
