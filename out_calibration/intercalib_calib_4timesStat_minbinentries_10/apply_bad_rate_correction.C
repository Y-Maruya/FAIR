#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

// Use nlohmann/json for JSON output
#include <nlohmann/json.hpp>
using json = nlohmann::json;

// ============================================================================
// Data Structures
// ============================================================================

struct InterCalibEntry {
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;
  double intercept = -999.0;
  double intercept_err = -999.0;
  double slope = -999.0;
  double slope_err = -999.0;
  int n_points_used = 0;
  int n_fit_points_over500 = -1;
  double chi2_ndf = -1.0;
  uint32_t quality_flag = 0;
  double distance_rms = -1.0;
  double distance_sigma = -1.0;
  double x_mm = -999.0;
  double y_mm = -999.0;
  double hg_adc_saturation = -999.0; // Added to store HG saturation point if available
};

struct RunData {
  std::string dirName = "";
  int runNumber = -1;
  int runNumberEnd = -1; // For directories that represent a range of runs (e.g., "22286-22296")
  std::vector<InterCalibEntry> entries;
  std::map<int, InterCalibEntry> byCellId;
};

struct ChannelCorrectionData {
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;
  double bad_rate = 0.0;
  bool should_mask = false;
  double corrected_slope = -999.0;
  double corrected_intercept = -999.0;
  long long corrected_n_points = 0;
  double corrected_distance_rms = -1.0;
  double corrected_distance_sigma = -1.0;
};

const uint32_t HGLG_LOW_STAT = 1 << 0;        // Bit 0
const uint32_t HGLG_FIT_FAILED = 1 << 1;      // Bit 1
const uint32_t HGLG_BAD_RMSPERSIGMA = 1 << 2; // Bit 2
const uint32_t HGLG_BAD_LOWGAIN = 1 << 3;     // Bit 3

// ============================================================================
// Utility Functions
// ============================================================================

bool loadRunFromFile(const std::string &filePath, const std::string &file2Path, RunData &out) {
  out = RunData();
  out.dirName = filePath;

  TFile *file = TFile::Open(filePath.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Failed to open file: " << filePath << std::endl;
    return false;
  }

  TTree *tree = dynamic_cast<TTree *>(file->Get("HGLGCalibQuality"));
  if (!tree) {
    std::cerr << "Missing 'HGLGCalibQuality' tree in: " << filePath << std::endl;
    file->Close();
    return false;
  }

  TFile *file2 = TFile::Open(file2Path.c_str(), "READ");
  if (!file2 || file2->IsZombie()) {
    std::cerr << "Failed to open file: " << file2Path << std::endl;
    file->Close();
    return false;
  }

  TTree *tree2 = dynamic_cast<TTree *>(file2->Get("intercalib"));
  if (!tree2) {
    std::cerr << "Missing 'intercalib' tree in: " << file2Path << std::endl;
    file->Close();
    file2->Close();
    return false;
  }

  InterCalibEntry e;
  tree->SetBranchAddress("cellid", &e.cellid);
  tree->SetBranchAddress("layer", &e.layer);
  tree->SetBranchAddress("chip", &e.chip);
  tree->SetBranchAddress("channel", &e.channel);
  tree->SetBranchAddress("intercept", &e.intercept);
  tree->SetBranchAddress("intercept_error", &e.intercept_err);
  tree->SetBranchAddress("slope", &e.slope);
  tree->SetBranchAddress("slope_error", &e.slope_err);
  tree->SetBranchAddress("n_points_used", &e.n_points_used);
  tree->SetBranchAddress("chi2_ndf", &e.chi2_ndf);
  tree->SetBranchAddress("quality_flag", &e.quality_flag);
  tree2->SetBranchAddress("hg_adc_saturation", &e.hg_adc_saturation);
  bool hasNFitPointsOver500 = tree->GetBranch("n_fit_points_over500") != nullptr;
  if (hasNFitPointsOver500) {
    tree->SetBranchAddress("n_fit_points_over500", &e.n_fit_points_over500);
  } else {
    std::cout << "Warning: n_fit_points_over500 branch not found in "
              << filePath << std::endl;
  }
  bool hasDistanceRMS = tree->GetBranch("distance_rms") != nullptr;
  bool hasDistanceSigma = tree->GetBranch("distance_sigma") != nullptr;
  if (hasDistanceRMS) tree->SetBranchAddress("distance_rms", &e.distance_rms);
  if (hasDistanceSigma) tree->SetBranchAddress("distance_sigma", &e.distance_sigma);

  bool hasCoordinates = tree->GetBranch("x_mm") != nullptr && tree->GetBranch("y_mm") != nullptr;
  if (hasCoordinates) {
    tree->SetBranchAddress("x_mm", &e.x_mm);
    tree->SetBranchAddress("y_mm", &e.y_mm);
  }

  const Long64_t nEntries = tree->GetEntries();
  out.entries.reserve(nEntries);

  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    tree2->GetEntry(i); // To get hg_adc_saturation from the second tree
    out.entries.push_back(e);
    out.byCellId[e.cellid] = e;
  }

  file->Close();
  return true;
}

std::map<int, double> calculateRunThresholds(const std::vector<RunData> &allRuns) {
  std::map<int, double> run_thresholds;  // runNumber -> threshold
  
  for (size_t run_idx = 0; run_idx < allRuns.size(); ++run_idx) {
    const auto &run = allRuns[run_idx];
    std::vector<double> run_ratios;
    
    // Collect all distance_rms / distance_sigma ratios
    for (const auto &entry : run.entries) {
      if (entry.layer == 21 && entry.chip == 1) {
        // Skip known bad channel L21C1
        continue;
      }
      if (entry.layer == 13 && entry.chip == 5) {
        // Skip known bad channel L13C5
        continue;
      }
      if (entry.distance_sigma > 0 && entry.distance_rms >= 0) {
        double ratio = entry.distance_rms / entry.distance_sigma;
        run_ratios.push_back(ratio);
      }
    }
    
    // Calculate mean and rms of ratios
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
    
    double threshold = mean_ratio + rms_ratio * 8.0;
    run_thresholds[run.runNumber] = threshold;
    std::cout << "[DEBUG] Run " << run.runNumber << ": ratio mean=" << std::fixed << std::setprecision(4)
              << mean_ratio << " rms=" << rms_ratio << " threshold=" << threshold << std::endl;
  }

  return run_thresholds;
}

uint32_t getEffectiveQualityFlags(const InterCalibEntry &entry, double threshold) {
  uint32_t effective_flags =
      entry.quality_flag & (HGLG_LOW_STAT | HGLG_FIT_FAILED | HGLG_BAD_LOWGAIN);

  if (entry.distance_sigma > 0 && entry.distance_rms >= 0) {
    const double ratio = entry.distance_rms / entry.distance_sigma;
    if (ratio > threshold) {
      effective_flags |= HGLG_BAD_RMSPERSIGMA;
    }
  }

  if (entry.n_fit_points_over500 >= 0 && entry.n_fit_points_over500 <= 50) {
    effective_flags |= HGLG_BAD_LOWGAIN;
  }

  return effective_flags;
}

bool isBadEffectiveFlag(uint32_t effective_flags) {
  return (effective_flags & (HGLG_LOW_STAT | HGLG_FIT_FAILED |
                             HGLG_BAD_RMSPERSIGMA | HGLG_BAD_LOWGAIN)) != 0;
}

std::string describeEffectiveFlags(uint32_t effective_flags) {
  std::vector<std::string> names;
  if (effective_flags & HGLG_LOW_STAT) names.push_back("LOW_STAT");
  if (effective_flags & HGLG_FIT_FAILED) names.push_back("FIT_FAILED");
  if (effective_flags & HGLG_BAD_RMSPERSIGMA) names.push_back("BAD_RMSPERSIGMA");
  if (effective_flags & HGLG_BAD_LOWGAIN) names.push_back("BAD_LOWGAIN");
  if (names.empty()) return "GOOD";

  std::ostringstream os;
  for (size_t i = 0; i < names.size(); ++i) {
    if (i > 0) os << "|";
    os << names[i];
  }
  return os.str();
}

// Calculate bad rate for each channel, with BAD_RMSPERSIGMA recomputed from distance metrics
std::map<int, double> calculateBadRates(const std::vector<RunData> &allRuns,
                                        const std::map<int, double> &run_thresholds) {
  std::map<int, double> badRates;
  
  // Collect all unique cellids
  std::set<int> all_cellids;
  for (const auto &run : allRuns) {
    for (const auto &[cid, entry] : run.byCellId) {
      all_cellids.insert(cid);
    }
  }

  // Second pass: for each cellid, calculate bad rate using recomputed quality flags
  for (int cid : all_cellids) {
    int bad_count = 0;
    int total_runs = allRuns.size();

    for (size_t run_idx = 0; run_idx < allRuns.size(); ++run_idx) {
      const auto &run = allRuns[run_idx];
      auto it = run.byCellId.find(cid);
      if (it != run.byCellId.end()) {
        const auto &entry = it->second;
        const double threshold = run_thresholds.at(run.runNumber);
        const uint32_t effective_flags = getEffectiveQualityFlags(entry, threshold);
        
        // Count as bad if any of bits 0, 1, 2, or 3 are set
        if (isBadEffectiveFlag(effective_flags)) {
          bad_count++;
        }
      }
    }

    double bad_rate = 100.0 * bad_count / total_runs;
    if (bad_count > 0) {
      badRates[cid] = bad_rate;
    }
  }

  return badRates;
}

bool runContainsRunNumber(const RunData &run, int runNumber) {
  if (run.runNumberEnd >= run.runNumber) {
    return run.runNumber <= runNumber && runNumber <= run.runNumberEnd;
  }
  return run.runNumber == runNumber;
}

double getHgSaturationForChannelInRun(const RunData &run,
                                      int layer,
                                      int chip,
                                      int channel) {
  for (const auto &entry : run.entries) {
    if (entry.layer == layer && entry.chip == chip && entry.channel == channel &&
        entry.hg_adc_saturation > 0) {
      return entry.hg_adc_saturation;
    }
  }
  return -999.0;
}

double getMaxHgSaturationForChannel(const std::vector<RunData> &allRuns,
                                    int layer,
                                    int chip,
                                    int channel) {
  double maxValue = -999.0;
  for (const auto &run : allRuns) {
    const double value = getHgSaturationForChannelInRun(run, layer, chip, channel);
    if (value > maxValue) maxValue = value;
  }
  return maxValue;
}

double getPost22499AverageHgSaturationForChannel(const std::vector<RunData> &allRuns,
                                                 int layer,
                                                 int chip,
                                                 int channel) {
  std::vector<double> values;
  for (const auto &run : allRuns) {
    if (run.runNumber <= 22705) continue;
    const double value = getHgSaturationForChannelInRun(run, layer, chip, channel);
    if (value > 0) values.push_back(value);
  }
  if (values.empty()) return -999.0;
  return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

double getRunByRunHgSaturationForChannel(const std::vector<RunData> &allRuns,
                                         int sameRun,
                                         int layer,
                                         int chip,
                                         int channel) {
  for (const auto &run : allRuns) {
    if (!runContainsRunNumber(run, sameRun)) continue;
    const double value = getHgSaturationForChannelInRun(run, layer, chip, channel);
    if (value > 0) return value;
  }
  return -999.0;
}

std::vector<double> buildHgSaturationArrayForJson(const std::vector<RunData> &allRuns,
                                                  int sameRun,
                                                  int layer,
                                                  int nChannelsPerLayer) {
  std::vector<double> hgSaturation(nChannelsPerLayer, -999.0);

  for (int chip = 0; chip < 9; ++chip) {
    for (int channel = 0; channel < 36; ++channel) {
      const int idx = chip * 36 + channel;
      hgSaturation[idx] = getMaxHgSaturationForChannel(allRuns, layer, chip, channel);
    }
  }

  if (layer == 13) {
    for (int channel = 30; channel <= 35; ++channel) {
      const int idx = 5 * 36 + channel;
      double value = -999.0;
      if (sameRun >= 22499 && sameRun <= 22705) {
        value = getPost22499AverageHgSaturationForChannel(allRuns, 13, 5, channel);
      } else {
        value = getRunByRunHgSaturationForChannel(allRuns, sameRun, 13, 5, channel);
      }
      if (value > 0) hgSaturation[idx] = value;
    }
  }

  return hgSaturation;
}

// Apply correction: mask if bad_rate >= 50%, average otherwise
void applyCorrectionAndWriteJson(const std::vector<RunData> &allRuns,
                                 const std::map<int, double> &badRates,
                                 const std::map<int, double> &run_thresholds,
                                 const std::string &outDir) {
  gSystem->Exec(Form("mkdir -p %s", outDir.c_str()));

  const int n_channels_per_layer = 36 * 9;  // 9 chips, 36 channels each

  // For each layer, write corrected JSON
  for (int layer = 0; layer < 41; ++layer) {
    std::vector<double> p0_arr(n_channels_per_layer, -999.0);
    std::vector<double> p0_err_arr(n_channels_per_layer, -999.0);
    std::vector<double> p1_arr(n_channels_per_layer, -999.0);
    std::vector<double> p1_err_arr(n_channels_per_layer, -999.0);
    std::vector<long long> n_points_arr(n_channels_per_layer, 0);
    std::vector<double> distance_rms_arr(n_channels_per_layer, -1.0);
    std::vector<double> distance_sigma_arr(n_channels_per_layer, -1.0);
    std::vector<int> quality_flag_arr(n_channels_per_layer, 1);  // Default to masked (1), will set to 0 or 2 later
    std::vector<double> hg_adc_saturation_arr(n_channels_per_layer, -999.0);

    // For each channel in this layer
    for (int chip = 0; chip < 9; ++chip) {
      for (int channel = 0; channel < 36; ++channel) {
        const int idx = chip * 36 + channel;

        // Collect values from all runs for this channel
        std::vector<double> p0_vals, p1_vals;
        std::vector<double> p0_err_vals, p1_err_vals;
        std::vector<long long> n_points_vals;
        std::vector<double> rms_vals, sigma_vals;
        std::vector<int> cellids_found;
        std::vector<double> hg_adc_saturation_vals;

        for (const auto &run : allRuns) {
          for (const auto &entry : run.entries) {
            if (entry.layer == layer && entry.chip == chip && entry.channel == channel) {
              cellids_found.push_back(entry.cellid);
              
              const double threshold = run_thresholds.at(run.runNumber);
              const uint32_t effective_flags = getEffectiveQualityFlags(entry, threshold);
              
              // Only use entries with good effective flags
              if (!isBadEffectiveFlag(effective_flags)) {
                p0_vals.push_back(entry.intercept);
                p1_vals.push_back(entry.slope);
                p0_err_vals.push_back(entry.intercept_err);
                p1_err_vals.push_back(entry.slope_err);
                n_points_vals.push_back(entry.n_points_used);
                if (entry.distance_rms >= 0) rms_vals.push_back(entry.distance_rms);
                if (entry.distance_sigma > 0) sigma_vals.push_back(entry.distance_sigma);
                if (entry.hg_adc_saturation > 0) hg_adc_saturation_vals.push_back(entry.hg_adc_saturation);
              }
            }
          }
        }

        // Check if this channel is in bad rate map
        int cellid = -1;
        if (!cellids_found.empty()) {
          cellid = cellids_found[0];  // Use first found cellid
        }

        double bad_rate = 0.0;
        bool should_mask = false;

        if (cellid > 0) {
          auto it = badRates.find(cellid);
          if (it != badRates.end()) {
            bad_rate = it->second;
            should_mask = (bad_rate >= 50.0);
          }
        }

        // DEBUG: Print first few channels of layer 0
        if (layer == 8 && chip == 7 && channel == 13) {
          std::cout << "\n[DEBUG] Layer 8, Chip 7, Channel " << channel << ":" << std::endl;
          std::cout << "  cellids_found: " << cellids_found.size() << std::endl;
          std::cout << "  cellid: " << cellid << std::endl;
          std::cout << "  p0_vals (good): " << p0_vals.size() << std::endl;
          std::cout << "  bad_rate: " << std::fixed << std::setprecision(1) << bad_rate << "%" << std::endl;
          std::cout << "  should_mask: " << should_mask << std::endl;
        }

        if (should_mask) {
          // Mask: keep -999.0
          quality_flag_arr[idx] = 1;  // Flag as masked due to high bad rate
        } else if (!p0_vals.empty()) {
          // Average from good runs
          p0_arr[idx] = std::accumulate(p0_vals.begin(), p0_vals.end(), 0.0) / p0_vals.size();
          p1_arr[idx] = std::accumulate(p1_vals.begin(), p1_vals.end(), 0.0) / p1_vals.size();
          if (!p0_err_vals.empty()) {
            p0_err_arr[idx] = std::accumulate(p0_err_vals.begin(), p0_err_vals.end(), 0.0) / p0_err_vals.size();
          }
          if (!p1_err_vals.empty()) {
            p1_err_arr[idx] = std::accumulate(p1_err_vals.begin(), p1_err_vals.end(), 0.0) / p1_err_vals.size();
          }
          n_points_arr[idx] = std::accumulate(n_points_vals.begin(), n_points_vals.end(), 0LL) / p0_vals.size();
          
          if (!rms_vals.empty()) {
            distance_rms_arr[idx] = std::accumulate(rms_vals.begin(), rms_vals.end(), 0.0) / rms_vals.size();
          }

          if (!sigma_vals.empty()) {
            distance_sigma_arr[idx] = std::accumulate(sigma_vals.begin(), sigma_vals.end(), 0.0) / sigma_vals.size();
          }

          quality_flag_arr[idx] = (bad_rate > 0.0) ? 2 : 0;  // 2 = averaged, 0 = good
          if (!hg_adc_saturation_vals.empty()) {
            hg_adc_saturation_arr[idx] = std::accumulate(hg_adc_saturation_vals.begin(), hg_adc_saturation_vals.end(), 0.0) / hg_adc_saturation_vals.size();
          }
        }
      }
    }

    // Write JSON for this layer using nlohmann/json
    for (const auto &run : allRuns) {
      std::ifstream ref_json_file(Form("./%s/json_muon/run%d_intercalib_Layer%d.json", (std::to_string(run.runNumber)+std::string("-")+std::to_string(run.runNumberEnd)).c_str(), run.runNumber, layer));
      json ref_json;
      if (ref_json_file.is_open()) {
        try {
          ref_json = json::parse(ref_json_file);
        } catch (const std::exception &e) {
          std::cerr << "Failed to parse reference JSON for run " << run.runNumber << " layer " << layer << ": " << e.what() << std::endl;
        }
        ref_json_file.close();
      } else {
        std::cerr << "Reference JSON not found for run " << run.runNumber << " layer " << layer << ", proceeding with limited metadata" << std::endl;
        std::cerr << "Expected path: " << Form("./%s/json_muon/run%d_intercalib_Layer%d.json", (std::to_string(run.runNumber)+std::string("-")+std::to_string(run.runNumberEnd)).c_str(), run.runNumber, layer) << std::endl;
      }
      std::vector<int> same_data_runs = ref_json.contains("Summary") && ref_json["Summary"].contains("SameDataRuns") ? ref_json["Summary"]["SameDataRuns"].get<std::vector<int>>() : std::vector<int>{};
      for (int same_run : same_data_runs) {
        std::ifstream same_run_json_file(Form("./%s/json_muon/run%d_intercalib_Layer%d.json", (std::to_string(run.runNumber)+std::string("-")+std::to_string(run.runNumberEnd)).c_str(), same_run, layer));
        json same_run_json;
        if (same_run_json_file.is_open()) {
          try {
            same_run_json = json::parse(same_run_json_file);
          } catch (const std::exception &e) {
            std::cerr << "Failed to parse JSON for same run " << same_run << " layer " << layer << ": " << e.what() << std::endl;
          }
          same_run_json_file.close();
        } else {
          std::cerr << "JSON for same run " << same_run << " layer " << layer << " not found, proceeding with limited metadata" << std::endl;
          std::cerr << "Expected path: " << Form("./%s/json_muon/run%d_intercalib_Layer%d.json", (std::to_string(run.runNumber)+std::string("-")+std::to_string(run.runNumberEnd)).c_str(), same_run, layer) << std::endl;  
        }
        json j;
        j["RunNumber"] = same_run;
        j["Status"] = 0;  // Assuming 0 means "OK"
        j["TimeStamp"] = same_run_json.contains("TimeStamp") ? same_run_json["TimeStamp"].get<std::string>() : "";
        j["Layer"] = layer;
        j["CalibrationType"] = "CorrectedIntercalib";
        // j["Description"] = "Intercalib with bad rate correction: mask if bad_rate>=50%, average if <50%";
        j["Summary"]["Entries"] = ref_json.contains("Summary") && ref_json["Summary"].contains("Entries") ? ref_json["Summary"]["Entries"].get<int>() : 0;
        j["Summary"]["BadRateMaskedChannels"] = std::count_if(quality_flag_arr.begin(), quality_flag_arr.end(), [](int flag) { return flag == 1; });
        j["Summary"]["AveragedChannels"] = std::count_if(quality_flag_arr.begin(), quality_flag_arr.end(), [](int flag) { return flag == 2; });
        j["Summary"]["NumUsedRun"] = ref_json.contains("Summary") && ref_json["Summary"].contains("NumUsedRun") ? ref_json["Summary"]["NumUsedRun"].get<int>() : 0;
        j["Summary"]["SameDataRuns"] = ref_json.contains("Summary") && ref_json["Summary"].contains("SameDataRuns") ? ref_json["Summary"]["SameDataRuns"].get<std::vector<int>>() : std::vector<int>{};
        j["PerChannel"]["Intercept"] = p0_arr;
        j["PerChannel"]["Slope"] = p1_arr;
        j["PerChannel"]["InterceptError"] = p0_err_arr;
        j["PerChannel"]["SlopeError"] = p1_err_arr;
        j["PerChannel"]["NPoints"] = n_points_arr;
        j["PerChannel"]["DistanceRMS"] = distance_rms_arr;
        j["PerChannel"]["DistanceSigma"] = distance_sigma_arr;
        j["PerChannel"]["QualityFlag"] = quality_flag_arr;
        j["PerChannel"]["HG_SaturationPoint"] =
            buildHgSaturationArrayForJson(allRuns, same_run, layer, n_channels_per_layer);
    

        std::ostringstream filename;
        filename << outDir << "/run" << same_run << "_intercalib_corrected_Layer" << layer << ".json";

        std::ofstream jout(filename.str());
        if (!jout.is_open()) {
          std::cerr << "Cannot write JSON: " << filename.str() << std::endl;
          continue;
        }

        jout << j.dump(2) << std::endl;
        jout.close();

        std::cout << "Wrote corrected JSON: " << filename.str() << std::endl;
      }
    }
  }
}

void printMaskedChannelDetails(const std::vector<RunData> &allRuns,
                               const std::map<int, double> &badRates,
                               const std::map<int, double> &run_thresholds) {
  std::cout << "  Masked channel details:" << std::endl;

  bool found_masked = false;
  for (const auto &[cid, rate] : badRates) {
    if (rate < 50.0) continue;

    found_masked = true;
    const InterCalibEntry *representative = nullptr;
    std::vector<std::string> bad_run_details;
    int observed_runs = 0;
    int bad_runs = 0;

    for (const auto &run : allRuns) {
      auto it = run.byCellId.find(cid);
      if (it == run.byCellId.end()) continue;

      const auto &entry = it->second;
      if (!representative) representative = &entry;
      observed_runs++;

      const double threshold = run_thresholds.at(run.runNumber);
      const uint32_t effective_flags = getEffectiveQualityFlags(entry, threshold);
      if (isBadEffectiveFlag(effective_flags)) {
        bad_runs++;
        std::ostringstream detail;
        detail << run.runNumber;
        if (run.runNumberEnd >= run.runNumber && run.runNumberEnd != run.runNumber) {
          detail << "-" << run.runNumberEnd;
        }
        detail << "(" << describeEffectiveFlags(effective_flags) << ")";
        bad_run_details.push_back(detail.str());
      }
    }

    std::cout << "    CellID " << cid;
    if (representative) {
      std::cout << " layer=" << representative->layer
                << " chip=" << representative->chip
                << " channel=" << representative->channel;
    }
    std::cout << " bad_rate=" << std::fixed << std::setprecision(1) << rate << "%"
              << " bad_runs=" << bad_runs << "/" << allRuns.size();
    if (observed_runs != static_cast<int>(allRuns.size())) {
      std::cout << " observed_runs=" << observed_runs << "/" << allRuns.size();
    }
    std::cout << std::endl;

    if (!bad_run_details.empty()) {
      std::cout << "      bad run flags: ";
      for (size_t i = 0; i < bad_run_details.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << bad_run_details[i];
      }
      std::cout << std::endl;
    }
  }

  if (!found_masked) {
    std::cout << "    none" << std::endl;
  }
}

// ============================================================================
// Main Function
// ============================================================================

void apply_bad_rate_correction(const char *baseDir = ".",
                               const char *runList = "all",
                               const char *outDir = "intercalib_corrected",
                               bool writeJson = false) {
  // Parse directory ranges
  std::vector<std::string> dirNames;
  if (std::string(runList) == std::string("all")) {
    gSystem->Exec("ls -d */ | sed 's#/##' > dir_list_temp.txt 2>/dev/null");
    std::ifstream ifs("dir_list_temp.txt");
    std::string line;
    while (std::getline(ifs, line)) {
      line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
      if (line.empty()) continue;
      dirNames.push_back(line);
    }
    ifs.close();
    gSystem->Exec("rm -f dir_list_temp.txt");
  } else {
    std::stringstream ss(runList);
    std::string token;
    while (std::getline(ss, token, ',')) {
      token.erase(std::remove_if(token.begin(), token.end(), ::isspace), token.end());
      if (token.empty()) continue;
      dirNames.push_back(token);
    }
  }

  if (dirNames.empty()) {
    std::cerr << "No directories to process" << std::endl;
    return;
  }

  std::cout << "Loading data from " << dirNames.size() << " runs..." << std::endl;

  // Load all runs
  std::vector<RunData> allRuns;
  for (const auto &dirName : dirNames) {
    const std::string path = std::string(baseDir) + "/" + dirName + "/hglg_calib_quality_muon.root";
    const std::string path2 = std::string(baseDir) + "/" + dirName + "/intercalib_adj_muon.root";
    RunData data;
    if (loadRunFromFile(path, path2, data)) {
      // Extract run number from directory name (e.g., "22286-22296" -> use the first number)
      size_t dash_pos = dirName.find('-');
      if (dash_pos != std::string::npos) {
        data.runNumber = std::stoi(dirName.substr(0, dash_pos));
        data.runNumberEnd = std::stoi(dirName.substr(dash_pos + 1));
      } else {
        data.runNumber = std::stoi(dirName);
      }
      std::cout << "  " << dirName << " (Run " << data.runNumber << "): " << data.entries.size() << " entries" << std::endl;
      allRuns.push_back(data);
    } else {
      std::cerr << "  Failed to load: " << dirName << std::endl;
    }
  }

  if (allRuns.size() < 2) {
    std::cerr << "Need at least 2 runs for comparison" << std::endl;
    return;
  }

  std::cout << "Successfully loaded " << allRuns.size() << " runs" << std::endl;

  // DEBUG: Print sample entries from first run
  if (!allRuns.empty()) {
    std::cout << "\n[DEBUG] First run entries sample (first 5):" << std::endl;
    for (size_t i = 0; i < std::min(size_t(5), allRuns[0].entries.size()); ++i) {
      const auto &e = allRuns[0].entries[i];
      std::cout << "  Entry " << i << ": cellid=" << e.cellid << " layer=" << e.layer 
                << " chip=" << e.chip << " channel=" << e.channel 
                << " intercept=" << e.intercept << " slope=" << e.slope 
                << " quality_flag=" << std::hex << e.quality_flag << std::dec 
                << " n_fit_points_over500=" << e.n_fit_points_over500
                << " is_bad=" << ((e.quality_flag & 0xf) != 0 ? "yes" : "no") << std::endl;
    }
  }

  // Calculate bad rate information
  std::cout << "\nCalculating bad rates..." << std::endl;
  auto run_thresholds = calculateRunThresholds(allRuns);
  auto badRates = calculateBadRates(allRuns, run_thresholds);

  // Print bad rate statistics
  int mask_count = 0, average_count = 0;
  for (const auto &[cid, rate] : badRates) {
    if (rate >= 50.0) {
      mask_count++;
    } else {
      average_count++;
    }
  }

  std::cout << "  Total problematic channels: " << badRates.size() << std::endl;
  std::cout << "  Channels to MASK (BadRate >= 50%): " << mask_count << std::endl;
  printMaskedChannelDetails(allRuns, badRates, run_thresholds);
  std::cout << "  Channels to AVERAGE (BadRate < 50%): " << average_count << std::endl;
  
  // DEBUG: Print sample bad rates
  if (!badRates.empty()) {
    std::cout << "\n[DEBUG] Sample bad rates (first 5):" << std::endl;
    int count = 0;
    for (const auto &[cid, rate] : badRates) {
      std::cout << "  CellID " << cid << ": " << std::fixed << std::setprecision(1) << rate << "%" << std::endl;
      if (++count >= 5) break;
    }
  }

  if (writeJson) {
    // Apply correction and write corrected JSON
    std::cout << "\nWriting corrected JSON files..." << std::endl;
    applyCorrectionAndWriteJson(allRuns, badRates, run_thresholds, outDir);
    std::cout << "\nDone! Corrected JSON files written to: " << outDir << std::endl;
  } else {
    std::cout << "\nDone! JSON output skipped. Add --write-json when running the executable, "
              << "or pass true as the 4th macro argument, to write JSON files." << std::endl;
  }
}

int main(int argc, char **argv) {
  const char *baseDir = ".";
  const char *runList = "all";
  const char *outDir = "intercalib_corrected";
  bool writeJson = false;

  int positional = 0;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--write-json" || arg == "-j") {
      writeJson = true;
      continue;
    }

    if (positional == 0) {
      baseDir = argv[i];
    } else if (positional == 1) {
      runList = argv[i];
    } else if (positional == 2) {
      outDir = argv[i];
    } else {
      std::cerr << "Unknown argument: " << arg << std::endl;
      std::cerr << "Usage: " << argv[0]
                << " [baseDir] [runList] [outDir] [--write-json|-j]" << std::endl;
      return 1;
    }
    positional++;
  }

  apply_bad_rate_correction(baseDir, runList, outDir, writeJson);
  return 0;
}
