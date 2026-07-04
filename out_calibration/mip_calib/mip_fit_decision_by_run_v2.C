#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TLine.h"
#include "TGraphErrors.h"

// ============================================================================
// Directory / branch utilities, following the style of compare.C
// ============================================================================

std::string trimCopy(std::string s) {
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
  return s;
}

std::vector<std::string> parseDirectoryRanges(const char *csv) {
  std::vector<std::string> dirs;
  if (!csv) return dirs;

  if (std::string(csv) == "all") {
    gSystem->Exec("ls -d */ | sed 's#/##' > dir_list.txt");
    std::ifstream ifs("dir_list.txt");
    std::string line;
    while (std::getline(ifs, line)) {
      line = trimCopy(line);
      if (!line.empty()) dirs.push_back(line);
    }
    return dirs;
  }

  std::stringstream ss(csv);
  std::string token;
  while (std::getline(ss, token, ',')) {
    token = trimCopy(token);
    if (!token.empty()) dirs.push_back(token);
  }
  return dirs;
}

bool setBranchIfExists(TTree *tree, const char *name, void *addr) {
  if (!tree || !tree->GetBranch(name)) {
    std::cerr << "  Warning: branch '" << name << "' not found in tree '" << tree->GetName() << "'" << std::endl;
    return false;
  }
  tree->SetBranchAddress(name, addr);
  return true;
}

bool setFirstExistingBranch(TTree *tree,
                            const std::vector<std::string> &names,
                            void *addr,
                            std::string &usedName) {
  for (const auto &name : names) {
    if (setBranchIfExists(tree, name.c_str(), addr)) {
      usedName = name;
      return true;
    }
  }
  usedName = "";
  return false;
}

std::set<int> parseCellIdList(const char *csv) {
  std::set<int> cellids;
  if (!csv) return cellids;

  std::stringstream ss(csv);
  std::string token;
  while (std::getline(ss, token, ',')) {
    token = trimCopy(token);
    if (token.empty()) continue;
    cellids.insert(std::atoi(token.c_str()));
  }
  return cellids;
}

std::string fileSafeString(const std::string &text) {
  std::string safe;
  for (char c : text) {
    safe += std::isalnum(static_cast<unsigned char>(c)) ? c : '_';
  }
  while (safe.find("__") != std::string::npos) {
    safe.replace(safe.find("__"), 2, "_");
  }
  if (safe.size() > 80) safe.resize(80);
  return safe.empty() ? "condition" : safe;
}

// ============================================================================
// Decision definitions
// ============================================================================

struct MIPFitEntry {
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;

  int entries = 0;

  double MPV = -999.0;
  double MPV_error = -999.0;
  double width = -999.0;
  double width_error = -999.0;
  double TotalArea = -999.0;
  double TotalArea_error = -999.0;
  double chi2 = -1.0;
  int ndf = -1;
  int fit_status = 999;
  int FitClass = 0;
  double gaus_sigma = -999.0;
  double gaus_sigma_error = -999.0;

  double refit_mpv = -999.0;
  double refit_mpv_error = -999.0;
  double refit_width = -999.0;
  double refit_width_error = -999.0;
  double refit_total_area = -999.0;
  double refit_total_area_error = -999.0;
  double refit_gaus_sigma = -999.0;
  double refit_gaus_sigma_error = -999.0;
  double refit_chi2_ndf = -999.0;
  double refit_rebin4_chi2_ndf = -999.0;
  int refit_ndf = -1;
  int refit_fit_status = 999;
  int refit_Class = 0;

  double direct_threshold = -999.0;
  double direct_threshold_error = -999.0;
  double direct_width = -999.0;
  double direct_width_error = -999.0;
  double direct_chi2 = -999.0;
  int direct_ndf = -1;
  double direct_chi2_ndf = -999.0;
  double direct_chi2_ndf_bin4 = -999.0;
  int direct_status = 999;
  bool direct_ok = false;
  bool direct_parameter_at_limit = false;
  double direct_efficiency = -999.0;
  double direct_efficiency_range = -999.0;
  double direct_lg_scale = -999.0;
  double direct_target_area = -999.0;
  double direct_actual_area = -999.0;
  double direct_x_min = -999.0;
  double direct_x_max = -999.0;
  int direct_nbins = -1;

  // Optional. If absent, the final AC branch cannot be fully resolved.
  double efficiency_data_mc = -999.0;
  double efficiency_data = -999.0;
  double efficiency_mc = -999.0;
  double percentile_15 = -999.0;

  double ratio_under_threshold = -999.0;
};

struct BranchStatus {
  bool has_cellid = false;
  bool has_layer = false;
  bool has_chip = false;
  bool has_channel = false;

  bool has_entries = false;
  bool has_MPV = false;
  bool has_MPV_error = false;
  bool has_width = false;
  bool has_width_error = false;
  bool has_TotalArea = false;
  bool has_TotalArea_error = false;
  bool has_chi2 = false;
  bool has_ndf = false;
  bool has_fit_status = false;
  bool has_FitClass = false;
  bool has_gaus_sigma = false;
  bool has_gaus_sigma_error = false;

  bool has_efficiency_data = false;
  bool has_efficiency_mc = false;
  bool has_efficiency_data_mc = false;

  bool has_percentile_15 = false;

  bool has_refit_mpv = false;
  bool has_refit_mpv_error = false;
  bool has_refit_width = false;
  bool has_refit_width_error = false;
  bool has_refit_total_area = false;
  bool has_refit_total_area_error = false;
  bool has_refit_gaus_sigma = false;
  bool has_refit_gaus_sigma_error = false;
  bool has_refit_chi2_ndf = false;
  bool has_refit_rebin4_chi2_ndf = false;
  bool has_refit_ndf = false;
  bool has_refit_fit_status = false;
  bool has_refit_Class = false;

  bool has_direct_threshold = false;
  bool has_direct_threshold_error = false;
  bool has_direct_width = false;
  bool has_direct_width_error = false;
  bool has_direct_chi2 = false;
  bool has_direct_ndf = false;
  bool has_direct_chi2_ndf = false;
  bool has_direct_chi2_ndf_bin4 = false;
  bool has_direct_status = false;
  bool has_direct_ok = false;
  bool has_direct_parameter_at_limit = false;
  bool has_direct_efficiency = false;
  bool has_direct_efficiency_range = false;
  bool has_direct_lg_scale = false;
  bool has_direct_target_area = false;
  bool has_direct_actual_area = false;
  bool has_direct_x_min = false;
  bool has_direct_x_max = false;
  bool has_direct_nbins = false;


  bool has_ratio_under_threshold = false;

  std::string efficiency_branch_name;
  std::string ref_mpv_branch_name;
};

enum DecisionClass {
  kInsufficientStatistics = 0,
  kFirstFitSuccess,
  kThresholdFitSuccess,
  kGoodFitWithoutThreshold,
  kSecondandFirstFitFailed,
  kSecondFitFailedFirstFitSuccess,
  kThresholdFitFailedButNeeded,
  kNDecisionClasses
};

enum StateContainer {
  kInsufficientStatisticsState = 1<<0,
  kFirstFitFailed_BigChi2Ndf = 1<<1,
  kFirstFitFailed_ParameterAtLimit = 1<<2,
  kFirstFitFailed_FitStatusNotZero = 1<<3,
  kFirstFitFailed_Others = 1<<4,
  kFirstFitSuccess_LowGausSigma = 1<<5,
  kFirstFitSuccess_KnownEfficiencyDegraded = 1<<6,
  kFirstFitSuccessState = 1<<7,
  kSecondFitFailed_RefitStatusNotZero = 1<<8,
  kSecondFitFailed_BigChi2Ndf = 1<<9,
  kSecondFitFailed_ParameterAtLimit = 1<<10,
  kSecondFitFailed_Others = 1<<11,
  kSecondFitSuccess_LowMPV = 1<<12,
  kSecondFitSuccess_LowGausSigma = 1<<13,
  kSecondFitSuccess_LowWidth = 1<<14,
  kSecondFitSuccess_HighMPVError = 1<<15,
  kSecondFitSuccessState = 1<<16,
  kThresholdFitFailed_BigChi2Ndf = 1<<29,
  kThresholdFitFailed_BigWidth = 1<<17,
  kThresholdFitFailed_BigWidthRatio = 1<<28,
  kThresholdFitFailed_LowLeftMean = 1<<18,
  kThresholdFitFailed_LowRightMean = 1<<19,
  kThresholdFitFailed_ParameterAtLimit = 1<<20,
  kThresholdFitFailed_BigUncertainty = 1<<30,
  kThresholdFitFailed_Others = 1<<21,
  kThresholdFitFailed_and_GoodEffMPVStability = 1<<22,
  kThresholdFitFailed_and_BadEffMPVStability = 1<<23,
  kThresholdFitSuccessState = 1<<24,
  kEfficiencyRatioOver1p0 = 1<<25,
  kEfficiencyRatioOver0p98 = 1<<26,
  kEfficiencyRatioUnder0p98 = 1<<27,
  kNStateContainers = 0
};

struct StateInfo {
  int bit;
  const char *group;
  const char *key;
  const char *label;
};

const std::vector<StateInfo> &stateInfos() {
  static const std::vector<StateInfo> states = {
      {kInsufficientStatisticsState, "Statistics", "insufficient_statistics", "Insufficient statistics"},
      {kFirstFitFailed_BigChi2Ndf, "First fit failed", "first_fit_big_chi2_ndf", "Big chi2/ndf"},
      {kFirstFitFailed_ParameterAtLimit, "First fit failed", "first_fit_parameter_at_limit", "Parameter at limit"},
      {kFirstFitFailed_FitStatusNotZero, "First fit failed", "first_fit_status_not_zero", "Fit status != 0"},
      {kFirstFitFailed_Others, "First fit failed", "first_fit_others", "Others"},
      {kFirstFitSuccess_LowGausSigma, "First fit quality", "first_fit_low_gaus_sigma", "Low Gaussian sigma"},
      {kFirstFitSuccess_KnownEfficiencyDegraded, "First fit quality", "known_efficiency_degraded", "Known efficiency degraded"},
      {kFirstFitSuccessState, "First fit quality", "first_fit_success", "First fit success"},
      {kSecondFitFailed_RefitStatusNotZero, "Second fit failed", "second_fit_status_not_zero", "Refit status != 0"},
      {kSecondFitFailed_BigChi2Ndf, "Second fit failed", "second_fit_big_chi2_ndf", "Big chi2/ndf"},
      {kSecondFitFailed_ParameterAtLimit, "Second fit failed", "second_fit_parameter_at_limit", "Parameter at limit"},
      {kSecondFitFailed_Others, "Second fit failed", "second_fit_others", "Others"},
      {kSecondFitSuccess_LowMPV, "Second fit quality", "second_fit_low_mpv", "Low MPV"},
      {kSecondFitSuccess_LowGausSigma, "Second fit quality", "second_fit_low_gaus_sigma", "Low Gaussian sigma"},
      {kSecondFitSuccess_LowWidth, "Second fit quality", "second_fit_low_width", "Low width"},
      {kSecondFitSuccess_HighMPVError, "Second fit quality", "second_fit_high_mpv_error", "High MPV error"},
      {kSecondFitSuccessState, "Second fit quality", "second_fit_success", "Second fit success"},
      {kThresholdFitFailed_BigChi2Ndf, "Threshold fit", "threshold_fit_big_chi2_ndf", "Big chi2/ndf"},
      {kThresholdFitFailed_BigWidth, "Threshold fit", "threshold_fit_big_width", "Big width"},
      {kThresholdFitFailed_BigWidthRatio, "Threshold fit", "threshold_fit_big_width_ratio", "Big width ratio"},
      {kThresholdFitFailed_BigUncertainty, "Threshold fit", "threshold_fit_big_uncertainty", "Big uncertainty"},
      {kThresholdFitFailed_LowLeftMean, "Threshold fit", "threshold_fit_low_left_mean", "Low left mean"},
      {kThresholdFitFailed_LowRightMean, "Threshold fit", "threshold_fit_low_right_mean", "Low right mean"},
      {kThresholdFitFailed_ParameterAtLimit, "Threshold fit", "threshold_fit_parameter_at_limit", "Parameter at limit"},
      {kThresholdFitFailed_Others, "Threshold fit", "threshold_fit_status_not_zero", "Fit status != 0"},
      {kThresholdFitFailed_and_GoodEffMPVStability, "Threshold fit", "threshold_fit_failed_but_good_eff_mpv_stability", "Failed threshold fit but good efficiency/MPV stability"},
      {kThresholdFitFailed_and_BadEffMPVStability, "Threshold fit", "threshold_fit_failed_and_bad_eff_mpv_stability", "Failed threshold fit and bad efficiency/MPV stability"},
      {kThresholdFitSuccessState, "Threshold fit", "threshold_fit_success", "Threshold fit success"},
      {kEfficiencyRatioOver1p0, "EfficiencyRatio", "efficiency_ratio_over_1p0", "Efficiency Ratio > 1.0"},
      {kEfficiencyRatioOver0p98, "EfficiencyRatio", "efficiency_ratio_0p98_to_1p0", "0.98 < Efficiency Ratio <= 1.0"},
      {kEfficiencyRatioUnder0p98, "EfficiencyRatio", "efficiency_ratio_under_0p98", "Efficiency Ratio <= 0.98"}};
  return states;
}

int countStateBits(int mask) {
  int count = 0;
  for (const auto &state : stateInfos()) {
    if (mask & state.bit) ++count;
  }
  return count;
}

std::string stateList(int mask, const char *separator = " | ") {
  std::ostringstream ss;
  bool first = true;
  for (const auto &state : stateInfos()) {
    if (!(mask & state.bit)) continue;
    if (!first) ss << separator;
    ss << state.key;
    first = false;
  }
  return first ? "(none)" : ss.str();
}

const char *decisionName(int c) {
  switch (c) {
    case kInsufficientStatistics: return "Insufficient statistics";
    case kFirstFitSuccess: return "First Fit Success";
    // case kSecondFitFailed: return "Second Fit Failed";
    case kSecondandFirstFitFailed: return "Second and First Fit Failed";
    case kSecondFitFailedFirstFitSuccess: return "Second Fit Failed, but First Fit Success";
    case kThresholdFitSuccess: return "Threshold fit Success";
    case kGoodFitWithoutThreshold: return "Good fit without threshold";
    case kThresholdFitFailedButNeeded: return "Threshold fit failed, but threshold needed";
    default: return "Unknown";
  }
}

const char *decisionKey(int c) {
  switch (c) {
    case kInsufficientStatistics: return "insufficient_statistics";
    case kFirstFitSuccess: return "first_fit_success";
    // case kSecondFitFailed: return "second_fit_failed";
    case kSecondandFirstFitFailed: return "second_and_first_fit_failed";
    case kSecondFitFailedFirstFitSuccess: return "second_fit_failed_first_fit_success";
    case kThresholdFitSuccess: return "threshold_fit_success";
    case kGoodFitWithoutThreshold: return "good_fit_without_threshold";
    case kThresholdFitFailedButNeeded: return "threshold_fit_failed_but_needed";
    default: return "unknown";
  }
}

int decisionColor(int c) {
  switch (c) {
    case kInsufficientStatistics: return kGray + 1;
    case kFirstFitSuccess: return kGreen + 2;
    case kSecondandFirstFitFailed: return kRed + 1;
    case kSecondFitFailedFirstFitSuccess: return kRed + 4;
    case kThresholdFitSuccess: return kAzure + 1;
    case kGoodFitWithoutThreshold: return kSpring + 5;
    case kThresholdFitFailedButNeeded: return kOrange + 7;
    default: return kBlack;
  }
}

struct DecisionResult {
  int decision = kInsufficientStatistics;
  std::string reason;
  int stateContainer = kNStateContainers; // Optional detailed state container for further analysis, not used in main decision logic
};

// Bit definitions copied conceptually from fit_histograms_noplot.C.
// For first/refit judgement this macro recomputes the conditions from the
// stored branch values as much as possible, instead of relying on FitClass.
bool isFirstParameterAtLimit(const MIPFitEntry &e, const BranchStatus &bs) {
  // Same logic as calculateFitClassification(): parameter-at-limit is inferred
  // from fitted parameter values, not from FitClass, whenever the branches exist.
  bool atLimit = false;
  if (bs.has_MPV && e.MPV <= 51.0) atLimit = true;
  if (bs.has_width && (e.width <= 11.0 || e.width > 89.0)) atLimit = true;
  if (bs.has_TotalArea && e.TotalArea <= 101.0) atLimit = true;
  if (bs.has_gaus_sigma && (e.gaus_sigma <= 1.0 || e.gaus_sigma > 149.0)) atLimit = true;
  return atLimit;
}

bool isSecondParameterAtLimit(const MIPFitEntry &e, const BranchStatus &bs) {
  // Same logic as calculateReFitClassification(): again recomputed from the
  // stored refit branches rather than using refit_Class.
  bool atLimit = false;
  if (bs.has_refit_mpv && e.refit_mpv < 2.0) atLimit = true;
  if (bs.has_refit_width && e.refit_width < 2.0) atLimit = true;
  if (bs.has_refit_gaus_sigma && e.refit_gaus_sigma < 2.0) atLimit = true;
  if (bs.has_refit_total_area && e.refit_total_area < 10.0) atLimit = true;
  if (bs.has_refit_width && e.refit_width > 89.0) atLimit = true;
  if (bs.has_refit_gaus_sigma && e.refit_gaus_sigma > 149.0) atLimit = true;
  return atLimit;
}

bool isThresholdParameterAtLimit(const MIPFitEntry &e, const BranchStatus &bs) {
  if (bs.has_direct_parameter_at_limit && e.direct_parameter_at_limit) return true;
  if (bs.has_direct_threshold && bs.has_direct_x_min &&
      e.direct_threshold < e.direct_x_min + 1.0) return true;
  if (bs.has_direct_threshold && bs.has_direct_x_max &&
      e.direct_threshold > e.direct_x_max - 1.0) return true;
  if (bs.has_direct_width && (e.direct_width < 4.0 || e.direct_width > 99.0)) return true;
  return false;
}

bool isKnownEfficiencyDegradedChannel(const MIPFitEntry &e,
                                      const std::set<int> &knownCellIds,
                                      const BranchStatus &bs) {
  if (knownCellIds.empty()) return false;
  if (!bs.has_efficiency_data_mc) return false;
  return knownCellIds.count(e.cellid) && e.efficiency_data_mc < 1.0;
}

DecisionResult classifyEntry(const MIPFitEntry &e,
                             const BranchStatus &bs,
                             const std::set<int> &knownCellIds) {
  DecisionResult out;
  
  if (e.efficiency_data_mc > 1.0) {
    out.stateContainer |= kEfficiencyRatioOver1p0;
  } else if (e.efficiency_data_mc > 0.98) {
    out.stateContainer |= kEfficiencyRatioOver0p98;
  } else if (e.efficiency_data_mc >= 0.0) { // Assuming efficiency_data_mc < 0 means not available
    out.stateContainer |= kEfficiencyRatioUnder0p98;
  }

  // B: Statistics > 300?
  if (e.entries <= 300) {
    out.decision = kInsufficientStatistics;
    out.reason = "entries <= 300";
    out.stateContainer |= kInsufficientStatisticsState;
    return out;
  }

  // D: Fit without Total Area constraint
  const double chi2_ndf = (e.ndf > 0) ? e.chi2 / double(e.ndf) : 1e30;
  const bool firstFitFailed_BigChi2Ndf = bs.has_chi2 && bs.has_ndf && chi2_ndf > 3.0;
  const bool firstFitFailed_ParameterAtLimit = isFirstParameterAtLimit(e, bs);
  const bool firstFitFailed_Others = bs.has_fit_status && e.fit_status != 0;
  const bool firstFitFailed = firstFitFailed_BigChi2Ndf || firstFitFailed_ParameterAtLimit || firstFitFailed_Others;
  const bool firstLowGausSigma = bs.has_gaus_sigma && e.gaus_sigma < 20.0;
  const bool knownEfficiencyDegraded = isKnownEfficiencyDegradedChannel(e, knownCellIds, bs);
  out.stateContainer |= (firstFitFailed_BigChi2Ndf ? kFirstFitFailed_BigChi2Ndf : 0);
  out.stateContainer |= (firstFitFailed_ParameterAtLimit ? kFirstFitFailed_ParameterAtLimit : 0);
  out.stateContainer |= (firstFitFailed_Others ? kFirstFitFailed_FitStatusNotZero : 0);
  out.stateContainer |= (firstLowGausSigma ? kFirstFitSuccess_LowGausSigma : 0);
  out.stateContainer |= (knownEfficiencyDegraded ? kFirstFitSuccess_KnownEfficiencyDegraded : 0);
  // H: First Fit Success is terminal, as in the provided flowchart.
  if (!firstFitFailed && !firstLowGausSigma && !knownEfficiencyDegraded) {
    out.decision = kFirstFitSuccess;
    out.stateContainer |= kFirstFitSuccessState;
    out.reason = "first fit: Others";
    return out;
  }

  // I: Fit with Total Area constraint
  const bool secondFitFailed_BigChi2Ndf = bs.has_refit_rebin4_chi2_ndf && e.refit_rebin4_chi2_ndf > 9.0;
  const bool secondFitFailed_RefitStatusNotZero = bs.has_refit_fit_status && e.refit_fit_status != 0;
  const bool secondFitFailed_ParameterAtLimit = isSecondParameterAtLimit(e, bs);
  const bool secondFitFailed = secondFitFailed_RefitStatusNotZero || secondFitFailed_ParameterAtLimit;
  out.stateContainer |= (secondFitFailed_RefitStatusNotZero ? kSecondFitFailed_RefitStatusNotZero : 0);
  out.stateContainer |= (secondFitFailed_BigChi2Ndf ? kSecondFitFailed_BigChi2Ndf : 0);
  out.stateContainer |= (secondFitFailed_ParameterAtLimit ? kSecondFitFailed_ParameterAtLimit : 0);
  if (secondFitFailed && firstFitFailed) {
    out.decision = kSecondandFirstFitFailed;
    out.reason = "refit_status != 0 or rebin4_chi2_ndf > 9 or refit parameter at limit and first fit failed";
    return out;
  }else if (secondFitFailed && !firstFitFailed) {
    out.decision = kSecondFitFailedFirstFitSuccess;
    out.reason = "refit_status != 0 or rebin4_chi2_ndf > 9 or refit parameter at limit, but first fit success";
    return out;
  }

  const bool lowMPV = bs.has_refit_mpv && e.refit_mpv < 50.0;
  const bool lowGausSigma = bs.has_refit_gaus_sigma && e.refit_gaus_sigma < 10.0;
  const bool lowWidth = bs.has_refit_width && e.refit_width < 12.0;
  const bool highMPVError = bs.has_refit_mpv && bs.has_refit_mpv_error &&
                            std::abs(e.refit_mpv) > 1e-12 &&
                            e.refit_mpv_error / e.refit_mpv > 0.3;

  out.stateContainer |= (lowMPV ? kSecondFitSuccess_LowMPV : 0);
  out.stateContainer |= (lowGausSigma ? kSecondFitSuccess_LowGausSigma : 0);
  out.stateContainer |= (lowWidth ? kSecondFitSuccess_LowWidth : 0);
  out.stateContainer |= (highMPVError ? kSecondFitSuccess_HighMPVError : 0);

  // if (!secondFitFailed && !lowMPV && !lowGausSigma && !lowWidth && !highMPVError) {
  if (!secondFitFailed && !highMPVError) {
    // out.decision = kSecondFitSuccess;
    out.stateContainer |= kSecondFitSuccessState;
    // out.reason = "second fit: Others";
    // return out;
  }
  // AB: direct threshold fit
  const bool thresholdFitFailed_BigChi2Ndf = bs.has_direct_chi2_ndf && e.direct_chi2_ndf > 9.0;
  const bool thresholdFitFailed_ParameterAtLimit = isThresholdParameterAtLimit(e, bs);
  const bool thresholdFitFailed_FitStatusNotZero =
      (bs.has_direct_status && e.direct_status != 0) ||
      (bs.has_direct_ok && !e.direct_ok);

  const bool thresholdFitFailed_BigUncertainty = (bs.has_direct_threshold_error && bs.has_direct_threshold &&
                                              std::abs(e.direct_threshold) > 1e-12 &&
                                              e.direct_threshold_error / std::abs(e.direct_threshold) > 0.5 )||
                                               ( bs.has_direct_width_error && bs.has_direct_width && 
                                                  std::abs(e.direct_width) > 1e-12 &&
                                                  e.direct_width_error / std::abs(e.direct_width) > 0.5);


  out.stateContainer |= (thresholdFitFailed_BigChi2Ndf ? kThresholdFitFailed_BigChi2Ndf : 0);
  out.stateContainer |= (thresholdFitFailed_ParameterAtLimit ? kThresholdFitFailed_ParameterAtLimit : 0);
  out.stateContainer |= (thresholdFitFailed_FitStatusNotZero ? kThresholdFitFailed_Others : 0);
  out.stateContainer |= (thresholdFitFailed_BigUncertainty ? kThresholdFitFailed_BigUncertainty : 0);
  
  const bool bigWidthRatio = bs.has_direct_width &&
                        bs.has_direct_threshold &&
                        std::abs(e.direct_threshold) > 1e-12 &&
                        e.direct_width / e.direct_threshold > 0.5;

  out.stateContainer |= (bigWidthRatio ? kThresholdFitFailed_BigWidthRatio : 0);

  const bool bigWidth = bs.has_direct_width && e.direct_width > 35.0;

  out.stateContainer |= (bigWidth ? kThresholdFitFailed_BigWidth : 0);

  const bool thresholdFitFailed = thresholdFitFailed_ParameterAtLimit || thresholdFitFailed_FitStatusNotZero || bigWidth || bigWidthRatio || thresholdFitFailed_BigChi2Ndf || thresholdFitFailed_BigUncertainty;
  if (!thresholdFitFailed) {
    out.decision = kThresholdFitSuccess;
    out.reason = "threshold fit: Others";
    out.stateContainer |= kThresholdFitSuccessState;
    return out;
  }

  // AC: Threshold Fit Failed, then final judgement.
  const double refMPVForCheck = e.refit_mpv > 0 ? e.refit_mpv : e.MPV; // Use refit_mpv if available, else MPV

  if (bs.has_efficiency_data_mc && bs.has_MPV &&
      e.efficiency_data_mc > 0.98 && std::abs(refMPVForCheck - e.MPV) < 20.0) {
    out.decision = kGoodFitWithoutThreshold;
    out.stateContainer |= kThresholdFitFailed_and_GoodEffMPVStability;
    out.reason = bigWidth ? "Big Width -> AC, but efficiency high and MPV stable":
                 "AC, but efficiency high and MPV stable";
    return out;
  }

  if (bs.has_efficiency_data_mc && bs.has_MPV &&
      (e.efficiency_data_mc <= 0.98 || std::abs(refMPVForCheck - e.MPV) >= 20.0)) {
    out.decision = kThresholdFitFailedButNeeded;
    out.stateContainer |= kThresholdFitFailed_and_BadEffMPVStability;
    out.reason = bigWidth ? "Big Width -> AC and efficiency data/MC < 0.98"
                          : "AC and efficiency data/MC < 0.98 or MPV unstable";
    return out;
  }
  return out;
}

// ============================================================================
// Run handling and output
// ============================================================================

struct FinalCalibrationEntry {
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;
  int entries = 0;
  int decision = kInsufficientStatistics;
  int state = 0;
  double mpv = -1.0;
  double width = -1.0;
  double gaus_sigma = -1.0;
  double direct_threshold = -1.0;
  double direct_width = -1.0;
  double mpv_error = -1.0;
  double width_error = -1.0;
  double direct_threshold_error = -1.0;
  double direct_width_error = -1.0;
  double direct_chi2_ndf_bin4 = -1.0;
  double TotalArea_error = -1.0;
  double gaus_sigma_error = -1.0;
  double chi2 = -1.0;
  int ndf = -1;
};

struct RunDecisionData {
  std::string dirName;
  std::string label;
  std::string path;
  long long totalEntries = 0;
  std::vector<int> counts = std::vector<int>(kNDecisionClasses, 0);
  std::map<std::string, int> reasonCounts;
  std::map<int, int> stateCounts;
  std::map<int, int> stateCombinationCounts;
  std::map<int, int> stateMultiplicityCounts;
  std::map<std::pair<int, int>, int> statePairCounts;
  std::map<int, int> layerChannelCounts;
  std::map<int, std::map<int, int> > layerStateCounts;
  std::map<std::pair<int, int>, int> decisionByDetectorPosition;
  std::map<std::string, std::vector<std::map<std::string, double> > > directParameterGroups;
  std::map<int, std::vector<std::pair<double, double> > > efficiencyDeficitByState;
  // Fit values keyed by cellid.  Keeping this per run makes it possible to
  // compare the same channel between non-identical input trees.
  std::map<int, std::map<std::string, double> > fitResultsByCell;
  std::map<int, int> decisionByCell;
  std::vector<FinalCalibrationEntry> finalCalibrationEntries;
};

FinalCalibrationEntry makeFinalCalibrationEntry(const MIPFitEntry &e,
                                                const BranchStatus &bs,
                                                const DecisionResult &result) {
  FinalCalibrationEntry finalEntry;
  finalEntry.cellid = e.cellid;
  finalEntry.layer = e.layer;
  finalEntry.chip = e.chip;
  finalEntry.channel = e.channel;
  finalEntry.entries = e.entries;
  finalEntry.decision = result.decision;
  finalEntry.state = result.stateContainer;

  if (result.decision == kFirstFitSuccess) {
    if (bs.has_MPV && std::isfinite(e.MPV) && e.MPV > -900.0) finalEntry.mpv = e.MPV;
    if (bs.has_width && std::isfinite(e.width) && e.width > -900.0) finalEntry.width = e.width;
    if (bs.has_gaus_sigma && std::isfinite(e.gaus_sigma) && e.gaus_sigma > -900.0)
      finalEntry.gaus_sigma = e.gaus_sigma;
    if (bs.has_MPV_error && std::isfinite(e.MPV_error) && e.MPV_error > -900.0)
      finalEntry.mpv_error = e.MPV_error;
    if (bs.has_width_error && std::isfinite(e.width_error) && e.width_error > -900.0)
      finalEntry.width_error = e.width_error;
    if (bs.has_TotalArea_error && std::isfinite(e.TotalArea_error) && e.TotalArea_error > -900.0)
      finalEntry.TotalArea_error = e.TotalArea_error;
    if (bs.has_gaus_sigma_error && std::isfinite(e.gaus_sigma_error) && e.gaus_sigma_error > -900.0)
      finalEntry.gaus_sigma_error = e.gaus_sigma_error;
    if (bs.has_chi2 && std::isfinite(e.chi2) && e.chi2 > -900.0)
      finalEntry.chi2 = e.chi2;
    if (bs.has_ndf && std::isfinite(e.ndf) && e.ndf > -900.0)
      finalEntry.ndf = e.ndf;
  } else if (result.stateContainer & kSecondFitSuccessState) {
    if (bs.has_refit_mpv && std::isfinite(e.refit_mpv) && e.refit_mpv > -900.0)
      finalEntry.mpv = e.refit_mpv;
    if (bs.has_refit_width && std::isfinite(e.refit_width) && e.refit_width > -900.0)
      finalEntry.width = e.refit_width;
    if (bs.has_refit_gaus_sigma && std::isfinite(e.refit_gaus_sigma) &&
        e.refit_gaus_sigma > -900.0)
      finalEntry.gaus_sigma = e.refit_gaus_sigma;
    if (bs.has_refit_mpv_error && std::isfinite(e.refit_mpv_error) &&
        e.refit_mpv_error > -900.0)
      finalEntry.mpv_error = e.refit_mpv_error;
    if (bs.has_refit_width_error && std::isfinite(e.refit_width_error) &&
        e.refit_width_error > -900.0)
      finalEntry.width_error = e.refit_width_error;
    if (bs.has_refit_total_area_error && std::isfinite(e.refit_total_area_error) && e.refit_total_area_error > -900.0)
      finalEntry.TotalArea_error = e.refit_total_area_error;
    if (bs.has_refit_gaus_sigma_error && std::isfinite(e.refit_gaus_sigma_error) &&
        e.refit_gaus_sigma_error > -900.0)
      finalEntry.gaus_sigma_error = e.refit_gaus_sigma_error;
    if (bs.has_refit_chi2_ndf && std::isfinite(e.refit_chi2_ndf) && e.refit_chi2_ndf > -900.0)
      finalEntry.chi2 = e.refit_chi2_ndf * static_cast<double>(e.refit_ndf);
    if (bs.has_refit_ndf && std::isfinite(e.refit_ndf) && e.refit_ndf > -900.0)
      finalEntry.ndf = e.refit_ndf;

  }

  if (result.decision == kThresholdFitSuccess) {
    if (bs.has_direct_threshold && std::isfinite(e.direct_threshold) &&
        e.direct_threshold > -900.0)
      finalEntry.direct_threshold = e.direct_threshold;
    if (bs.has_direct_width && std::isfinite(e.direct_width) && e.direct_width > -900.0)
      finalEntry.direct_width = e.direct_width;
    if (bs.has_direct_threshold_error && std::isfinite(e.direct_threshold_error) &&
        e.direct_threshold_error > -900.0)
      finalEntry.direct_threshold_error = e.direct_threshold_error;
    if (bs.has_direct_width_error && std::isfinite(e.direct_width_error) &&
        e.direct_width_error > -900.0)
      finalEntry.direct_width_error = e.direct_width_error;
    if (bs.has_direct_chi2_ndf_bin4 && std::isfinite(e.direct_chi2_ndf_bin4) &&
        e.direct_chi2_ndf_bin4 > -900.0)
      finalEntry.direct_chi2_ndf_bin4 = e.direct_chi2_ndf_bin4;
  }
  return finalEntry;
}

std::map<std::string, double> fitResultValues(const MIPFitEntry &e,
                                              const BranchStatus &bs) {
  std::map<std::string, double> values;
#define STORE_FIT_VALUE(branch, member) \
  if (bs.branch && std::isfinite(static_cast<double>(e.member))) \
    values[#member] = static_cast<double>(e.member)
  STORE_FIT_VALUE(has_MPV, MPV);
  STORE_FIT_VALUE(has_MPV_error, MPV_error);
  STORE_FIT_VALUE(has_width, width);
  STORE_FIT_VALUE(has_width_error, width_error);
  STORE_FIT_VALUE(has_gaus_sigma, gaus_sigma);
  STORE_FIT_VALUE(has_gaus_sigma_error, gaus_sigma_error);
  STORE_FIT_VALUE(has_TotalArea, TotalArea);
  STORE_FIT_VALUE(has_TotalArea_error, TotalArea_error);
  STORE_FIT_VALUE(has_chi2, chi2);
  STORE_FIT_VALUE(has_ndf, ndf);
  STORE_FIT_VALUE(has_fit_status, fit_status);
  STORE_FIT_VALUE(has_refit_mpv, refit_mpv);
  STORE_FIT_VALUE(has_refit_mpv_error, refit_mpv_error);
  STORE_FIT_VALUE(has_refit_width, refit_width);
  STORE_FIT_VALUE(has_refit_width_error, refit_width_error);
  STORE_FIT_VALUE(has_refit_gaus_sigma, refit_gaus_sigma);
  STORE_FIT_VALUE(has_refit_gaus_sigma_error, refit_gaus_sigma_error);
  STORE_FIT_VALUE(has_refit_chi2_ndf, refit_chi2_ndf);
  STORE_FIT_VALUE(has_refit_rebin4_chi2_ndf, refit_rebin4_chi2_ndf);
  STORE_FIT_VALUE(has_refit_ndf, refit_ndf);
  STORE_FIT_VALUE(has_refit_fit_status, refit_fit_status);
  STORE_FIT_VALUE(has_direct_threshold, direct_threshold);
  STORE_FIT_VALUE(has_direct_threshold_error, direct_threshold_error);
  STORE_FIT_VALUE(has_direct_width, direct_width);
  STORE_FIT_VALUE(has_direct_width_error, direct_width_error);
  STORE_FIT_VALUE(has_direct_chi2_ndf_bin4, direct_chi2_ndf_bin4);
  STORE_FIT_VALUE(has_direct_status, direct_status);
  STORE_FIT_VALUE(has_efficiency_data_mc, efficiency_data_mc);
#undef STORE_FIT_VALUE
  if (bs.has_chi2 && bs.has_ndf && e.ndf > 0 && std::isfinite(e.chi2))
    values["chi2_ndf"] = e.chi2 / static_cast<double>(e.ndf);
  if (bs.has_direct_chi2 && bs.has_direct_ndf && e.direct_ndf > 0 &&
      std::isfinite(e.direct_chi2))
    values["direct_chi2_ndf"] = e.direct_chi2 / static_cast<double>(e.direct_ndf);
  else if (bs.has_direct_chi2_ndf && std::isfinite(e.direct_chi2_ndf))
    values["direct_chi2_ndf"] = e.direct_chi2_ndf;
  return values;
}

std::map<std::string, double> directParameterValues(const MIPFitEntry &e,
                                                    const BranchStatus &bs) {
  std::map<std::string, double> values;
#define STORE_DIRECT_VALUE(name) \
  if (bs.has_##name && std::isfinite(static_cast<double>(e.name))) \
    values[#name] = static_cast<double>(e.name)
  STORE_DIRECT_VALUE(direct_threshold);
  STORE_DIRECT_VALUE(direct_threshold_error);
  STORE_DIRECT_VALUE(direct_width);
  STORE_DIRECT_VALUE(direct_width_error);
  STORE_DIRECT_VALUE(direct_chi2_ndf_bin4);
  STORE_DIRECT_VALUE(direct_status);
  STORE_DIRECT_VALUE(direct_ok);
  STORE_DIRECT_VALUE(direct_parameter_at_limit);
  STORE_DIRECT_VALUE(direct_efficiency);
  STORE_DIRECT_VALUE(direct_efficiency_range);
  STORE_DIRECT_VALUE(direct_lg_scale);
  STORE_DIRECT_VALUE(direct_target_area);
  STORE_DIRECT_VALUE(direct_actual_area);
  STORE_DIRECT_VALUE(direct_x_min);
  STORE_DIRECT_VALUE(direct_x_max);
  STORE_DIRECT_VALUE(direct_nbins);
#undef STORE_DIRECT_VALUE
  if (bs.has_direct_chi2 && bs.has_direct_ndf && e.direct_ndf > 0 &&
      std::isfinite(e.direct_chi2)) {
    values["direct_chi2_ndf"] = e.direct_chi2 / static_cast<double>(e.direct_ndf);
  } else if (bs.has_direct_chi2_ndf && std::isfinite(e.direct_chi2_ndf)) {
    values["direct_chi2_ndf"] = e.direct_chi2_ndf;
  }
  if (bs.has_efficiency_data_mc && std::isfinite(e.efficiency_data_mc)) {
    values["efficiency_ratio"] = e.efficiency_data_mc;
  }
  if (bs.has_refit_chi2_ndf && std::isfinite(e.refit_chi2_ndf)) {
    values["refit_chi2_ndf"] = e.refit_chi2_ndf;
  }
  return values;
}

void storeDirectParameterGroups(const MIPFitEntry &e,
                                const BranchStatus &bs,
                                const DecisionResult &result,
                                RunDecisionData &out) {
  const std::map<std::string, double> values = directParameterValues(e, bs);
  if (values.empty()) return;

  const std::pair<int, const char *> groups[] = {
      {kThresholdFitSuccessState, "success"},
      {kThresholdFitFailed_ParameterAtLimit, "fail_parameter_at_limit"},
      {kThresholdFitFailed_Others, "fail_fit_status_not_zero"},
      {kThresholdFitFailed_BigWidth, "fail_big_width"},
      {kThresholdFitFailed_BigWidthRatio, "fail_big_width_ratio"},
      {kThresholdFitFailed_BigChi2Ndf, "fail_big_chi2_ndf"},
      {kThresholdFitFailed_BigUncertainty, "fail_big_uncertainty"}};
  for (const auto &group : groups) {
    if (result.stateContainer & group.first) {
      out.directParameterGroups[group.second].push_back(values);
    }
  }
}

bool setupBranches(TTree *tree, MIPFitEntry &e, BranchStatus &bs) {
  bs.has_cellid = setBranchIfExists(tree, "cellid", &e.cellid);
  bs.has_layer = setBranchIfExists(tree, "layer", &e.layer);
  bs.has_chip = setBranchIfExists(tree, "chip", &e.chip);
  bs.has_channel = setBranchIfExists(tree, "channel", &e.channel);

  bs.has_entries = setBranchIfExists(tree, "entries", &e.entries);
  bs.has_MPV = setBranchIfExists(tree, "MPV", &e.MPV);
  bs.has_MPV_error = setBranchIfExists(tree, "MPV_error", &e.MPV_error);
  bs.has_width = setBranchIfExists(tree, "width", &e.width);
  bs.has_width_error = setBranchIfExists(tree, "width_error", &e.width_error);
  bs.has_TotalArea = setBranchIfExists(tree, "TotalArea", &e.TotalArea);
  bs.has_TotalArea_error = setBranchIfExists(tree, "TotalArea_error", &e.TotalArea_error);
  bs.has_chi2 = setBranchIfExists(tree, "chi2", &e.chi2);
  bs.has_ndf = setBranchIfExists(tree, "ndf", &e.ndf);
  bs.has_fit_status = setBranchIfExists(tree, "fit_status", &e.fit_status);
  bs.has_FitClass = setBranchIfExists(tree, "FitClass", &e.FitClass);
  bs.has_gaus_sigma = setBranchIfExists(tree, "gaus_sigma", &e.gaus_sigma);
  bs.has_gaus_sigma_error = setBranchIfExists(tree, "gaus_sigma_error", &e.gaus_sigma_error);

  bs.has_refit_mpv = setBranchIfExists(tree, "refit_mpv", &e.refit_mpv);
  bs.has_refit_mpv_error = setBranchIfExists(tree, "refit_mpv_error", &e.refit_mpv_error);
  bs.has_refit_width = setBranchIfExists(tree, "refit_width", &e.refit_width);
  bs.has_refit_width_error = setBranchIfExists(tree, "refit_width_error", &e.refit_width_error);
  bs.has_refit_total_area = setBranchIfExists(tree, "refit_total_area", &e.refit_total_area);
  bs.has_refit_total_area_error = setBranchIfExists(tree, "refit_total_area_error", &e.refit_total_area_error);
  bs.has_refit_gaus_sigma_error = setBranchIfExists(tree, "refit_gaus_sigma_error", &e.refit_gaus_sigma_error);
  bs.has_refit_gaus_sigma = setBranchIfExists(tree, "refit_gaus_sigma", &e.refit_gaus_sigma);
  bs.has_refit_chi2_ndf = setBranchIfExists(tree, "refit_chi2_ndf", &e.refit_chi2_ndf);
  bs.has_refit_rebin4_chi2_ndf = setBranchIfExists(tree, "refit_rebin4_chi2_ndf", &e.refit_rebin4_chi2_ndf);
  bs.has_refit_ndf = setBranchIfExists(tree, "refit_ndf", &e.refit_ndf);
  bs.has_refit_fit_status = setBranchIfExists(tree, "refit_fit_status", &e.refit_fit_status);
  bs.has_refit_Class = setBranchIfExists(tree, "refit_Class", &e.refit_Class);

#define SET_DIRECT_BRANCH(name) \
  bs.has_##name = setBranchIfExists(tree, #name, &e.name)
  SET_DIRECT_BRANCH(direct_threshold);
  SET_DIRECT_BRANCH(direct_threshold_error);
  SET_DIRECT_BRANCH(direct_width);
  SET_DIRECT_BRANCH(direct_width_error);
  SET_DIRECT_BRANCH(direct_chi2);
  SET_DIRECT_BRANCH(direct_ndf);
  SET_DIRECT_BRANCH(direct_chi2_ndf);
  SET_DIRECT_BRANCH(direct_chi2_ndf_bin4);
  SET_DIRECT_BRANCH(direct_status);
  SET_DIRECT_BRANCH(direct_ok);
  SET_DIRECT_BRANCH(direct_parameter_at_limit);
  SET_DIRECT_BRANCH(direct_efficiency);
  SET_DIRECT_BRANCH(direct_efficiency_range);
  SET_DIRECT_BRANCH(direct_lg_scale);
  SET_DIRECT_BRANCH(direct_target_area);
  SET_DIRECT_BRANCH(direct_actual_area);
  SET_DIRECT_BRANCH(direct_x_min);
  SET_DIRECT_BRANCH(direct_x_max);
  SET_DIRECT_BRANCH(direct_nbins);
#undef SET_DIRECT_BRANCH

  bs.has_efficiency_data = setBranchIfExists(tree, "efficiency", &e.efficiency_data);
  bs.has_efficiency_mc = setBranchIfExists(tree, "MC_efficiency", &e.efficiency_mc);
  bs.has_efficiency_data_mc = setFirstExistingBranch(
      tree,
      {"efficiency_data_mc", "efficiency_dataMC", "eff_data_mc", "efficiency_ratio", "hit_efficiency_data_mc"},
      &e.efficiency_data_mc,
      bs.efficiency_branch_name);
  bs.has_percentile_15 = setBranchIfExists(tree, "percentile_15", &e.percentile_15);
  bs.has_ratio_under_threshold = setBranchIfExists(tree, "refit_ratio_under_threshold", &e.ratio_under_threshold);
  if (!bs.has_entries || !bs.has_cellid) {
    std::cerr << "Missing required branch: entries and/or cellid" << std::endl;
    return false;
  }
  return true;
}

bool processRun(const std::string &baseDir,
                const std::string &dirName,
                const std::string &fileName,
                const std::string &treeName,
                const std::set<int> &knownCellIds,
                const std::string &outDir,
                bool writePerEntryCsv,
                RunDecisionData &out) {
  out = RunDecisionData();
  out.dirName = dirName;
  out.label = dirName;
  out.path = baseDir + "/" + dirName + "/" + fileName;

  TFile *file = TFile::Open(out.path.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Failed to open file: " << out.path << std::endl;
    return false;
  }

  TTree *tree = dynamic_cast<TTree *>(file->Get(treeName.c_str()));
  if (!tree) {
    std::cerr << "Missing tree '" << treeName << "' in: " << out.path << std::endl;
    file->Close();
    return false;
  }

  MIPFitEntry e;
  BranchStatus bs;
  if (!setupBranches(tree, e, bs)) {
    file->Close();
    return false;
  }

  if (!bs.has_efficiency_data_mc) {
    std::cout << "  Warning " << dirName << ": efficiency data/MC branch not found. "
              << "AC final judgement may become unresolved." << std::endl;
  } else {
    std::cout << "  " << dirName << ": using efficiency branch '" << bs.efficiency_branch_name << "'" << std::endl;
  }
  if (!bs.has_width || !bs.has_TotalArea) {
    std::cout << "  Warning " << dirName << ": first-fit parameter-at-limit check is partial "
              << "because width and/or TotalArea branch is missing." << std::endl;
  }
  if (!bs.has_refit_total_area) {
    std::cout << "  Warning " << dirName << ": refit_total_area branch not found. "
              << "Refit parameter-at-limit check is partial." << std::endl;
  }

  std::ofstream perEntry;
  if (writePerEntryCsv) {
    perEntry.open((outDir + "/DecisionEntries_" + dirName + ".csv").c_str());
    perEntry << "run,cellid,layer,chip,channel,entries,decision,reason,state_mask,state_count,states\n";
  }

  const Long64_t nEntries = tree->GetEntries();
  out.totalEntries = nEntries;
  for (Long64_t i = 0; i < nEntries; ++i) {
    // Reset optional values to avoid stale values when branches are missing.
    e = MIPFitEntry();
    tree->GetEntry(i);

    DecisionResult result = classifyEntry(e, bs, knownCellIds);
    out.finalCalibrationEntries.push_back(makeFinalCalibrationEntry(e, bs, result));
    out.fitResultsByCell[e.cellid] = fitResultValues(e, bs);
    out.decisionByCell[e.cellid] = result.decision;
    out.counts[result.decision]++;
    out.reasonCounts[result.reason]++;
    storeDirectParameterGroups(e, bs, result, out);
    out.layerChannelCounts[e.layer]++;
    if (bs.has_layer && bs.has_chip && bs.has_channel &&
        e.layer >= 0 && e.chip >= 0 && e.channel >= 0) {
      out.decisionByDetectorPosition[std::make_pair(e.layer * 9 + e.chip, e.channel)] =
          result.decision;
    }
    for (const auto &state : stateInfos()) {
      if (result.stateContainer & state.bit) {
        out.stateCounts[state.bit]++;
        out.layerStateCounts[e.layer][state.bit]++;
        if (bs.has_efficiency_data_mc && bs.has_direct_efficiency &&
            std::isfinite(e.efficiency_data_mc) && std::isfinite(e.direct_efficiency) &&
            e.efficiency_data_mc >= 0.0 && e.direct_efficiency >= 0.0) {
          out.efficiencyDeficitByState[state.bit].push_back(
              std::make_pair(1.0 - e.efficiency_data_mc, 1.0 - e.direct_efficiency));
        }
      }
    }
    for (size_t a = 0; a < stateInfos().size(); ++a) {
      if (!(result.stateContainer & stateInfos()[a].bit)) continue;
      for (size_t b = a; b < stateInfos().size(); ++b) {
        if (!(result.stateContainer & stateInfos()[b].bit)) continue;
        out.statePairCounts[std::make_pair(stateInfos()[a].bit, stateInfos()[b].bit)]++;
      }
    }
    out.stateCombinationCounts[result.stateContainer]++;
    out.stateMultiplicityCounts[countStateBits(result.stateContainer)]++;

    if (writePerEntryCsv && perEntry) {
      perEntry << dirName << "," << e.cellid << "," << e.layer << "," << e.chip << ","
               << e.channel << "," << e.entries << ","
               << decisionKey(result.decision) << ",\"" << result.reason << "\","
               << result.stateContainer << "," << countStateBits(result.stateContainer)
               << ",\"" << stateList(result.stateContainer) << "\"\n";
    }
  }

  file->Close();
  return true;
}

void exportChannelsByStateAndCondition(
    const char *baseDir,
    const char *runList,
    const char *fileName,
    const char *treeName,
    const char *outDir,
    int requiredStateMask,
    const char *condition,
    const char *knownEfficiencyDegradedCellIds = "") {
  if (requiredStateMask == 0) {
    std::cerr << "exportChannelsByStateAndCondition: requiredStateMask must be non-zero" << std::endl;
    return;
  }
  if (!condition || std::string(condition).empty()) {
    std::cerr << "exportChannelsByStateAndCondition: condition must not be empty" << std::endl;
    return;
  }

  gSystem->mkdir(outDir, true);
  const std::string stem = "channels_state_mask_" + std::to_string(requiredStateMask) +
                           "_" + fileSafeString(condition);
  const std::string detailPath = std::string(outDir) + "/" + stem + "_details.txt";
  const std::string cellidPath = std::string(outDir) + "/" + stem + "_cellids.txt";
  std::ofstream details(detailPath.c_str());
  std::set<int> selectedCellIds;
  long long selectedEntries = 0;

  details << "# required_state_mask=" << requiredStateMask << "\n"
          << "# required_states=" << stateList(requiredStateMask) << "\n"
          << "# condition=" << condition << "\n"
          << "run,entry,cellid,layer,chip,channel,entries,decision,state_mask,state_count,states,"
          << "condition_value,efficiency_ratio,MPV,width,chi2_ndf,fit_status,"
          << "refit_mpv,refit_width,refit_chi2_ndf,refit_rebin4_chi2_ndf,refit_fit_status,"
          << "direct_threshold,direct_width,direct_chi2_ndf,direct_status,direct_ok,"
          << "direct_parameter_at_limit,direct_efficiency\n";

  const std::set<int> knownCellIds = parseCellIdList(knownEfficiencyDegradedCellIds);
  for (const auto &dirName : parseDirectoryRanges(runList)) {
    const std::string path = std::string(baseDir) + "/" + dirName + "/" + fileName;
    TFile *file = TFile::Open(path.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "  Warning: cannot open " << path << std::endl;
      delete file;
      continue;
    }
    TTree *tree = dynamic_cast<TTree *>(file->Get(treeName));
    if (!tree) {
      std::cerr << "  Warning: missing tree '" << treeName << "' in " << path << std::endl;
      file->Close();
      delete file;
      continue;
    }

    MIPFitEntry e;
    BranchStatus bs;
    if (!setupBranches(tree, e, bs)) {
      file->Close();
      delete file;
      continue;
    }
    const std::string formulaName = "selection_" + fileSafeString(dirName);
    TTreeFormula formula(formulaName.c_str(), condition, tree);
    if (formula.GetNdim() <= 0) {
      std::cerr << "  Invalid condition for " << dirName << ": " << condition << std::endl;
      file->Close();
      delete file;
      continue;
    }

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
      e = MIPFitEntry();
      tree->GetEntry(i);
      const DecisionResult result = classifyEntry(e, bs, knownCellIds);
      if ((result.stateContainer & requiredStateMask) != requiredStateMask) continue;
      const double conditionValue = formula.EvalInstance();
      if (conditionValue == 0.0 || !std::isfinite(conditionValue)) continue;

      const double chi2Ndf = e.ndf > 0 ? e.chi2 / static_cast<double>(e.ndf) : -999.0;
      const double directChi2Ndf =
          e.direct_ndf > 0 ? e.direct_chi2 / static_cast<double>(e.direct_ndf)
                           : e.direct_chi2_ndf;
      details << dirName << "," << i << "," << e.cellid << "," << e.layer << ","
              << e.chip << "," << e.channel << "," << e.entries << ","
              << decisionKey(result.decision) << "," << result.stateContainer << ","
              << countStateBits(result.stateContainer) << ",\"" << stateList(result.stateContainer)
              << "\"," << conditionValue << "," << e.efficiency_data_mc << ","
              << e.MPV << "," << e.width << "," << chi2Ndf << "," << e.fit_status << ","
              << e.refit_mpv << "," << e.refit_width << "," << e.refit_chi2_ndf << ","
              << e.refit_rebin4_chi2_ndf << "," << e.refit_fit_status << ","
              << e.direct_threshold << "," << e.direct_width << "," << directChi2Ndf << ","
              << e.direct_status << "," << e.direct_ok << "," << e.direct_parameter_at_limit
              << "," << e.direct_efficiency << "\n";
      selectedCellIds.insert(e.cellid);
      ++selectedEntries;
    }
    file->Close();
    delete file;
  }

  std::ofstream cellids(cellidPath.c_str());
  cellids << "cellid={";
  bool first = true;
  for (int cellid : selectedCellIds) {
    if (!first) cellids << ",";
    cellids << cellid;
    first = false;
  }
  cellids << "};\n";

  std::cout << "Selected " << selectedEntries << " entries and " << selectedCellIds.size()
            << " unique cellids\n"
            << "  details: " << detailPath << "\n"
            << "  cellids: " << cellidPath << std::endl;
}

void writeStateSummaryTables(const std::vector<RunDecisionData> &runs, const std::string &outDir) {
  std::ofstream stateCsv((outDir + "/StateSummary_byRun.csv").c_str());
  stateCsv << "run,total_channels,state_group,state_key,state_label,channels_with_state,fraction_of_channels\n";
  for (const auto &run : runs) {
    for (const auto &state : stateInfos()) {
      const int count = run.stateCounts.count(state.bit) ? run.stateCounts.at(state.bit) : 0;
      stateCsv << run.label << "," << run.totalEntries << "," << state.group << ","
               << state.key << ",\"" << state.label << "\"," << count << ","
               << (run.totalEntries > 0 ? double(count) / run.totalEntries : 0.0) << "\n";
    }
  }

  std::ofstream multiplicityCsv((outDir + "/StateOverlapMultiplicity_byRun.csv").c_str());
  multiplicityCsv << "run,total_channels,state_count_per_channel,count,fraction_of_channels\n";
  for (const auto &run : runs) {
    for (const auto &kv : run.stateMultiplicityCounts) {
      multiplicityCsv << run.label << "," << run.totalEntries << "," << kv.first << ","
                      << kv.second << ","
                      << (run.totalEntries > 0 ? double(kv.second) / run.totalEntries : 0.0) << "\n";
    }
  }

  std::ofstream combinationCsv((outDir + "/StateOverlapCombinations_byRun.csv").c_str());
  combinationCsv << "run,total_channels,state_mask,state_count,combination_count,fraction_of_channels,states\n";
  for (const auto &run : runs) {
    std::vector<std::pair<int, int> > combinations(run.stateCombinationCounts.begin(),
                                                   run.stateCombinationCounts.end());
    std::sort(combinations.begin(), combinations.end(),
              [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
                return a.second > b.second;
              });
    for (const auto &kv : combinations) {
      combinationCsv << run.label << "," << run.totalEntries << "," << kv.first << ","
                     << countStateBits(kv.first) << "," << kv.second << ","
                     << (run.totalEntries > 0 ? double(kv.second) / run.totalEntries : 0.0)
                     << ",\"" << stateList(kv.first) << "\"\n";
    }
  }

  std::ofstream layerCsv((outDir + "/StateSummary_byRunLayer.csv").c_str());
  layerCsv << "run,layer,total_channels,state_group,state_key,state_label,channels_with_state,"
              "fraction_of_layer_channels\n";
  for (const auto &run : runs) {
    for (const auto &layerCount : run.layerChannelCounts) {
      const int layer = layerCount.first;
      const int total = layerCount.second;
      for (const auto &state : stateInfos()) {
        int count = 0;
        const auto layerIt = run.layerStateCounts.find(layer);
        if (layerIt != run.layerStateCounts.end()) {
          const auto stateIt = layerIt->second.find(state.bit);
          if (stateIt != layerIt->second.end()) count = stateIt->second;
        }
        layerCsv << run.label << "," << layer << "," << total << "," << state.group << ","
                 << state.key << ",\"" << state.label << "\"," << count << ","
                 << (total > 0 ? double(count) / total : 0.0) << "\n";
      }
    }
  }

  std::ofstream txt((outDir + "/StateSummary_byRun.txt").c_str());
  txt << "================================================================================\n";
  txt << "STATE CONTAINER SUMMARY BY RUN\n";
  txt << "================================================================================\n";
  txt << "Each state fraction uses all channels as its denominator. States can overlap,\n";
  txt << "so the state fractions below do NOT sum to 100%. The multiplicity and exact\n";
  txt << "combination sections are exclusive and each sums to 100%.\n";
  txt << "States are evaluated only along the classification path reached by a channel.\n\n";

  for (const auto &run : runs) {
    long long totalTags = 0;
    for (const auto &kv : run.stateCounts) totalTags += kv.second;
    txt << "Run: " << run.label << "\n";
    txt << "Total channels: " << run.totalEntries << "\n";
    txt << "Total state tags (overlap included): " << totalTags << "\n";
    txt << "Average state tags/channel: " << std::fixed << std::setprecision(3)
        << (run.totalEntries > 0 ? double(totalTags) / run.totalEntries : 0.0) << "\n\n";

    txt << "Individual states (overlap included):\n";
    txt << std::left << std::setw(24) << "Group" << std::setw(38) << "State"
        << std::right << std::setw(20) << "Channels / Total" << std::setw(12) << "Fraction" << "\n";
    txt << std::string(94, '-') << "\n";
    for (const auto &state : stateInfos()) {
      const int count = run.stateCounts.count(state.bit) ? run.stateCounts.at(state.bit) : 0;
      txt << std::left << std::setw(24) << state.group << std::setw(38) << state.label
          << std::right << std::setw(9) << count << " / " << std::left << std::setw(8)
          << run.totalEntries << std::right << std::setw(11) << std::setprecision(2)
          << (run.totalEntries > 0 ? 100.0 * count / run.totalEntries : 0.0) << "%\n";
    }

    txt << "\nNumber of states per channel (exclusive):\n";
    for (const auto &kv : run.stateMultiplicityCounts) {
      txt << "  " << std::setw(2) << kv.first << " states: " << std::setw(8) << kv.second
          << "  " << std::setw(7) << std::setprecision(2)
          << (run.totalEntries > 0 ? 100.0 * kv.second / run.totalEntries : 0.0) << "%\n";
    }

    std::vector<std::pair<int, int> > combinations(run.stateCombinationCounts.begin(),
                                                   run.stateCombinationCounts.end());
    std::sort(combinations.begin(), combinations.end(),
              [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
                return a.second > b.second;
              });
    txt << "\nTop exact state combinations (exclusive, top 15):\n";
    const size_t nShown = std::min<size_t>(15, combinations.size());
    for (size_t i = 0; i < nShown; ++i) {
      txt << "  " << std::setw(8) << combinations[i].second << "  " << std::setw(7)
          << std::setprecision(2)
          << (run.totalEntries > 0 ? 100.0 * combinations[i].second / run.totalEntries : 0.0)
          << "%  " << stateList(combinations[i].first) << "\n";
    }
    txt << "\n";
  }
}

void writeSummaryTables(const std::vector<RunDecisionData> &runs, const std::string &outDir) {
  std::ofstream csv((outDir + "/DecisionSummary_byRun.csv").c_str());
  csv << "run,total";
  for (int c = 0; c < kNDecisionClasses; ++c) csv << "," << decisionKey(c);
  for (int c = 0; c < kNDecisionClasses; ++c) csv << ",frac_" << decisionKey(c);
  csv << "\n";

  for (const auto &run : runs) {
    int total = 0;
    for (int c = 0; c < kNDecisionClasses; ++c) total += run.counts[c];
    csv << run.label << "," << total;
    for (int c = 0; c < kNDecisionClasses; ++c) csv << "," << run.counts[c];
    for (int c = 0; c < kNDecisionClasses; ++c) {
      csv << "," << (total > 0 ? double(run.counts[c]) / total : 0.0);
    }
    csv << "\n";
  }

  std::ofstream txt((outDir + "/DecisionSummary_byRun.txt").c_str());
  txt << "================================================================================\n";
  txt << "MIP FIT DECISION SUMMARY BY RUN\n";
  txt << "================================================================================\n\n";

  for (const auto &run : runs) {
    int total = 0;
    for (int c = 0; c < kNDecisionClasses; ++c) total += run.counts[c];

    txt << "Run: " << run.label << "\n";
    txt << "Path: " << run.path << "\n";
    txt << "Total: " << total << "\n";
    txt << std::left << std::setw(48) << "Decision" << std::right
        << std::setw(10) << "Count" << std::setw(12) << "Fraction" << "\n";
    txt << std::string(70, '-') << "\n";
    for (int c = 0; c < kNDecisionClasses; ++c) {
      txt << std::left << std::setw(48) << decisionName(c) << std::right
          << std::setw(10) << run.counts[c]
          << std::setw(11) << std::fixed << std::setprecision(2)
          << (total > 0 ? 100.0 * run.counts[c] / total : 0.0) << "%\n";
    }

    txt << "\nReason counts:\n";
    for (const auto &kv : run.reasonCounts) {
      txt << "  " << std::left << std::setw(65) << kv.first << " " << kv.second << "\n";
    }
    txt << "\n";
  }
}

void drawDecisionStacks(const std::vector<RunDecisionData> &runs, const std::string &outDir) {
  if (runs.empty()) return;

  gStyle->SetOptStat(0);

  const int nRuns = static_cast<int>(runs.size());

  auto makeFrame = [&](const char *name, const char *title, const char *ytitle) {
    TH1D *h = new TH1D(name, title, nRuns, 0.5, nRuns + 0.5);
    h->GetYaxis()->SetTitle(ytitle);
    h->GetXaxis()->SetTitle("Run");
    for (int i = 0; i < nRuns; ++i) {
      h->GetXaxis()->SetBinLabel(i + 1, runs[i].label.c_str());
    }
    h->SetStats(0);
    return h;
  };

  // Count stack
  TCanvas *c1 = new TCanvas("c_decision_stack_counts", "Decision counts by run", 1200, 700);
  c1->SetBottomMargin(0.18);
  THStack *stackCounts = new THStack("stackCounts", "MIP fit decision counts by run;Run;Number of channels");
  TLegend *leg1 = new TLegend(0.58, 0.55, 0.92, 0.90);
  leg1->SetBorderSize(1);
  leg1->SetFillStyle(0);

  for (int c = 0; c < kNDecisionClasses; ++c) {
    TH1D *h = makeFrame(Form("h_count_%d", c), "", "Number of channels");
    for (int r = 0; r < nRuns; ++r) h->SetBinContent(r + 1, runs[r].counts[c]);
    h->SetFillColor(decisionColor(c));
    h->SetLineColor(kBlack);
    stackCounts->Add(h);
    leg1->AddEntry(h, decisionName(c), "f");
  }

  stackCounts->Draw("hist");
  stackCounts->GetXaxis()->SetTitle("Run");
  stackCounts->GetYaxis()->SetTitle("Number of channels");
  stackCounts->GetXaxis()->LabelsOption("v");
  leg1->Draw();
  c1->SaveAs((outDir + "/DecisionStack_counts_byRun.pdf").c_str());
  c1->SaveAs((outDir + "/DecisionStack_counts_byRun.png").c_str());

  // Fraction stack
  TCanvas *c2 = new TCanvas("c_decision_stack_fraction", "Decision fractions by run", 1200, 700);
  c2->SetBottomMargin(0.18);
  THStack *stackFrac = new THStack("stackFrac", "MIP fit decision fractions by run;Run;Fraction");
  TLegend *leg2 = new TLegend(0.58, 0.55, 0.92, 0.90);
  leg2->SetBorderSize(1);
  leg2->SetFillStyle(0);

  for (int c = 0; c < kNDecisionClasses; ++c) {
    TH1D *h = makeFrame(Form("h_frac_%d", c), "", "Fraction");
    for (int r = 0; r < nRuns; ++r) {
      int total = 0;
      for (int cc = 0; cc < kNDecisionClasses; ++cc) total += runs[r].counts[cc];
      h->SetBinContent(r + 1, total > 0 ? double(runs[r].counts[c]) / total : 0.0);
    }
    h->SetFillColor(decisionColor(c));
    h->SetLineColor(kBlack);
    stackFrac->Add(h);
    leg2->AddEntry(h, decisionName(c), "f");
  }

  stackFrac->Draw("hist");
  stackFrac->SetMinimum(0.0);
  stackFrac->SetMaximum(1.0);
  stackFrac->GetXaxis()->SetTitle("Run");
  stackFrac->GetYaxis()->SetTitle("Fraction");
  stackFrac->GetXaxis()->LabelsOption("v");
  leg2->Draw();
  c2->SaveAs((outDir + "/DecisionStack_fraction_byRun.pdf").c_str());
  c2->SaveAs((outDir + "/DecisionStack_fraction_byRun.png").c_str());
}

void drawDecisionDetectorMaps(const std::vector<RunDecisionData> &runs,
                              const std::string &outDir) {
  if (runs.empty()) return;

  gStyle->SetOptStat(0);
  const std::string plotDir = outDir + "/DecisionDetectorMap";
  gSystem->mkdir(plotDir.c_str(), true);
  TFile outputFile((plotDir + "/decision_detector_maps.root").c_str(), "RECREATE");

  int palette[kNDecisionClasses];
  double contours[kNDecisionClasses];
  for (int decision = 0; decision < kNDecisionClasses; ++decision) {
    palette[decision] = decisionColor(decision);
    contours[decision] = decision + 0.5;
  }
  gStyle->SetPalette(kNDecisionClasses, palette);

  for (const auto &run : runs) {
    if (run.decisionByDetectorPosition.empty()) continue;

    int minX = run.decisionByDetectorPosition.begin()->first.first;
    int maxX = minX;
    int minChannel = run.decisionByDetectorPosition.begin()->first.second;
    int maxChannel = minChannel;
    for (const auto &entry : run.decisionByDetectorPosition) {
      minX = std::min(minX, entry.first.first);
      maxX = std::max(maxX, entry.first.first);
      minChannel = std::min(minChannel, entry.first.second);
      maxChannel = std::max(maxChannel, entry.first.second);
    }

    const std::string histName = "h_decision_detector_map_" + run.label;
    const std::string title = "MIP fit decision detector map: " + run.label +
                              ";layer*9+chip;channel;Decision";
    TH2D *hist = new TH2D(histName.c_str(), title.c_str(),
                          maxX - minX + 1, minX - 0.5, maxX + 0.5,
                          maxChannel - minChannel + 1, minChannel - 0.5, maxChannel + 0.5);
    hist->SetMinimum(0.5);
    hist->SetMaximum(kNDecisionClasses + 0.5);
    hist->SetContour(kNDecisionClasses, contours);
    for (const auto &entry : run.decisionByDetectorPosition) {
      hist->SetBinContent(hist->FindBin(entry.first.first, entry.first.second), entry.second + 1);
    }

    const std::string canvasName = "c_decision_detector_map_" + run.label;
    TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 3000, 850);
    canvas->SetRightMargin(0.25);
    canvas->SetBottomMargin(0.12);
    hist->Draw("COLZ");

    TLegend *legend = new TLegend(0.78, 0.16, 0.99, 0.88);
    legend->SetHeader("Decision", "C");
    legend->SetBorderSize(1);
    std::vector<TH1D *> legendEntries;
    for (int decision = 0; decision < kNDecisionClasses; ++decision) {
      TH1D *entry = new TH1D(Form("h_decision_legend_%s_%d", run.label.c_str(), decision),
                            "", 1, 0.0, 1.0);
      entry->SetFillColor(decisionColor(decision));
      entry->SetLineColor(kBlack);
      legendEntries.push_back(entry);
      legend->AddEntry(entry, decisionName(decision), "f");
    }
    legend->Draw();

    const std::string outputStem = plotDir + "/DecisionDetectorMap_" + run.label;
    canvas->SaveAs((outputStem + ".pdf").c_str());
    canvas->SaveAs((outputStem + ".png").c_str());
    outputFile.cd();
    hist->Write();
    delete legend;
    for (TH1D *entry : legendEntries) delete entry;
    delete canvas;
    delete hist;
  }
  outputFile.Close();
  gStyle->SetPalette(kBird);
}

void drawStateSummaries(const std::vector<RunDecisionData> &runs, const std::string &outDir) {
  if (runs.empty()) return;

  gStyle->SetOptStat(0);
  const int nRuns = static_cast<int>(runs.size());
  const int nStates = static_cast<int>(stateInfos().size());

  TCanvas *c1 = new TCanvas("c_state_channel_count_heatmap", "Channels with each state by run",
                            std::max(1200, std::min(3000, 90 * nRuns + 500)),
                            std::max(850, 32 * nStates + 250));
  c1->SetLeftMargin(0.32);
  c1->SetRightMargin(0.14);
  c1->SetBottomMargin(0.24);
  TH2D *heatmap = new TH2D("h_state_channel_count_heatmap",
                           "Channels with each state (overlap included);Run (total channels);State",
                           nRuns, 0.5, nRuns + 0.5, nStates, 0.5, nStates + 0.5);
  for (int r = 0; r < nRuns; ++r) {
    const std::string runLabel = runs[r].label + " (N=" + std::to_string(runs[r].totalEntries) + ")";
    heatmap->GetXaxis()->SetBinLabel(r + 1, runLabel.c_str());
  }
  for (int s = 0; s < nStates; ++s) {
    const StateInfo &state = stateInfos()[s];
    const std::string axisLabel = std::string(state.group) + ": " + state.label;
    heatmap->GetYaxis()->SetBinLabel(nStates - s, axisLabel.c_str());
    for (int r = 0; r < nRuns; ++r) {
      const int count = runs[r].stateCounts.count(state.bit) ? runs[r].stateCounts.at(state.bit) : 0;
      heatmap->SetBinContent(r + 1, nStates - s, count);
    }
  }
  heatmap->GetZaxis()->SetTitle("Channels with state");
  heatmap->GetXaxis()->LabelsOption("v");
  heatmap->GetYaxis()->SetLabelSize(0.025);
  gStyle->SetPaintTextFormat(".0f");
  heatmap->Draw(nRuns <= 20 ? "COLZ TEXT" : "COLZ");
  c1->SaveAs((outDir + "/StateChannelCount_heatmap_byRun.pdf").c_str());
  c1->SaveAs((outDir + "/StateChannelCount_heatmap_byRun.png").c_str());

  int maxMultiplicity = 0;
  for (const auto &run : runs) {
    if (!run.stateMultiplicityCounts.empty()) {
      maxMultiplicity = std::max(maxMultiplicity, run.stateMultiplicityCounts.rbegin()->first);
    }
  }

  TCanvas *c2 = new TCanvas("c_state_overlap_multiplicity", "State overlap multiplicity by run", 1200, 700);
  c2->SetBottomMargin(0.18);
  c2->SetRightMargin(0.22);
  THStack *stack = new THStack(
      "stateMultiplicityStack",
      "Number of StateContainer flags per channel (exclusive);Run;Fraction of channels");
  TLegend *legend = new TLegend(0.80, 0.55, 0.97, 0.90);
  legend->SetBorderSize(1);
  legend->SetFillStyle(0);
  const int colors[] = {kGray, kAzure + 1, kGreen + 2, kOrange + 7, kRed + 1,
                        kViolet + 1, kCyan + 2, kPink + 1, kYellow + 2};
  const int nColors = sizeof(colors) / sizeof(colors[0]);

  for (int multiplicity = 0; multiplicity <= maxMultiplicity; ++multiplicity) {
    TH1D *h = new TH1D(Form("h_state_multiplicity_%d", multiplicity), "", nRuns, 0.5, nRuns + 0.5);
    for (int r = 0; r < nRuns; ++r) {
      h->GetXaxis()->SetBinLabel(r + 1, runs[r].label.c_str());
      const int count = runs[r].stateMultiplicityCounts.count(multiplicity)
                            ? runs[r].stateMultiplicityCounts.at(multiplicity)
                            : 0;
      h->SetBinContent(r + 1, runs[r].totalEntries > 0 ? double(count) / runs[r].totalEntries : 0.0);
    }
    h->SetFillColor(colors[multiplicity % nColors]);
    h->SetLineColor(kBlack);
    stack->Add(h);
    legend->AddEntry(h, Form("%d state%s", multiplicity, multiplicity == 1 ? "" : "s"), "f");
  }

  stack->Draw("hist");
  stack->SetMinimum(0.0);
  stack->SetMaximum(1.0);
  stack->GetXaxis()->LabelsOption("v");
  legend->Draw();
  c2->SaveAs((outDir + "/StateOverlapMultiplicity_fraction_byRun.pdf").c_str());
  c2->SaveAs((outDir + "/StateOverlapMultiplicity_fraction_byRun.png").c_str());
}

void drawStateLayerHeatmaps(const std::vector<RunDecisionData> &runs, const std::string &outDir) {
  if (runs.empty()) return;

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".0f");
  const int nStates = static_cast<int>(stateInfos().size());

  for (const auto &run : runs) {
    if (run.layerChannelCounts.empty()) continue;

    std::vector<int> layers;
    for (const auto &kv : run.layerChannelCounts) layers.push_back(kv.first);
    const int nLayers = static_cast<int>(layers.size());
    const std::string canvasName = "c_state_layer_heatmap_" + run.label;
    TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(),
                                  std::max(1400, std::min(3000, 75 * nLayers + 500)),
                                  std::max(850, 32 * nStates + 250));
    canvas->SetLeftMargin(0.32);
    canvas->SetRightMargin(0.14);
    canvas->SetBottomMargin(nLayers > 20 ? 0.30 : 0.24);

    const std::string histName = "h_state_layer_heatmap_" + run.label;
    const std::string title = "Channels with each state by layer: " + run.label +
                              " (overlap included);Layer (total channels);State";
    TH2D *heatmap = new TH2D(histName.c_str(), title.c_str(),
                             nLayers, 0.5, nLayers + 0.5, nStates, 0.5, nStates + 0.5);

    for (int l = 0; l < nLayers; ++l) {
      const int layer = layers[l];
      const std::string layerLabel = std::to_string(layer) + " (N=" +
                                     std::to_string(run.layerChannelCounts.at(layer)) + ")";
      heatmap->GetXaxis()->SetBinLabel(l + 1, layerLabel.c_str());
    }
    for (int s = 0; s < nStates; ++s) {
      const StateInfo &state = stateInfos()[s];
      const std::string axisLabel = std::string(state.group) + ": " + state.label;
      heatmap->GetYaxis()->SetBinLabel(nStates - s, axisLabel.c_str());
      for (int l = 0; l < nLayers; ++l) {
        const int layer = layers[l];
        int count = 0;
        const auto layerIt = run.layerStateCounts.find(layer);
        if (layerIt != run.layerStateCounts.end()) {
          const auto stateIt = layerIt->second.find(state.bit);
          if (stateIt != layerIt->second.end()) count = stateIt->second;
        }
        heatmap->SetBinContent(l + 1, nStates - s, count);
      }
    }

    heatmap->GetZaxis()->SetTitle("Channels with state");
    heatmap->GetXaxis()->LabelsOption("v");
    heatmap->GetXaxis()->SetLabelSize(nLayers > 20 ? 0.018 : 0.025);
    heatmap->GetXaxis()->SetTitleOffset(nLayers > 20 ? 3.6 : 2.4);
    heatmap->GetYaxis()->SetLabelSize(0.025);
    heatmap->Draw(nLayers <= 20 ? "COLZ TEXT" : "COLZ");
    canvas->SaveAs((outDir + "/StateChannelCount_heatmap_byLayer_" + run.label + ".pdf").c_str());
    canvas->SaveAs((outDir + "/StateChannelCount_heatmap_byLayer_" + run.label + ".png").c_str());
  }
}

void drawStatePairHeatmaps(const std::vector<RunDecisionData> &runs, const std::string &outDir) {
  if (runs.empty()) return;

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".0f");
  const int nStates = static_cast<int>(stateInfos().size());
  TFile outputFile((outDir + "/StatePairChannelCount_heatmaps.root").c_str(), "RECREATE");

  for (const auto &run : runs) {
    const std::string canvasName = "c_state_pair_heatmap_" + run.label;
    TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1800, 1550);
    canvas->SetLeftMargin(0.27);
    canvas->SetRightMargin(0.16);
    canvas->SetBottomMargin(0.28);

    const std::string histName = "h_state_pair_heatmap_" + run.label;
    const std::string title = "State i AND State j channel counts: " + run.label +
                              ";State i;State j";
    TH2D *heatmap = new TH2D(histName.c_str(), title.c_str(),
                             nStates, 0.5, nStates + 0.5, nStates, 0.5, nStates + 0.5);

    for (int i = 0; i < nStates; ++i) {
      const std::string label = std::string(stateInfos()[i].group) + ": " + stateInfos()[i].label;
      heatmap->GetXaxis()->SetBinLabel(i + 1, label.c_str());
      heatmap->GetYaxis()->SetBinLabel(nStates - i, label.c_str());
      for (int j = 0; j < nStates; ++j) {
        const int firstBit = stateInfos()[std::min(i, j)].bit;
        const int secondBit = stateInfos()[std::max(i, j)].bit;
        const auto it = run.statePairCounts.find(std::make_pair(firstBit, secondBit));
        const int count = it != run.statePairCounts.end() ? it->second : 0;
        heatmap->SetBinContent(i + 1, nStates - j, count);
      }
    }
    heatmap->GetZaxis()->SetTitle("Channels with state i AND state j");
    heatmap->GetXaxis()->LabelsOption("v");
    heatmap->GetXaxis()->SetLabelSize(0.018);
    heatmap->GetYaxis()->SetLabelSize(0.018);
    heatmap->SetMarkerSize(0.5);
    heatmap->GetXaxis()->SetTitleSize(0.0);
    heatmap->Draw("COLZ TEXT");

    TLine*line_h1 = new TLine(1.5, -0.5, 1.5, nStates + 0.5);
    TLine*line_h2 = new TLine(8.5, -0.5, 8.5, nStates + 0.5);
    TLine*line_h3 = new TLine(17.5, -0.5, 17.5, nStates + 0.5);
    TLine*line_h4 = new TLine(nStates + 1 - 3.5,-0.5,nStates + 1 - 3.5,nStates+0.5);

    line_h1->SetLineColor(kBlack);
    line_h1->SetLineWidth(2);
    line_h2->SetLineColor(kBlack);
    line_h2->SetLineWidth(2);
    line_h3->SetLineColor(kBlack);
    line_h3->SetLineWidth(2);
    line_h4->SetLineColor(kBlack);
    line_h4->SetLineWidth(2);
    line_h1->Draw("same");
    line_h2->Draw("same");
    line_h3->Draw("same");
    line_h4->Draw("same");
    
    TLine*line_v1 = new TLine(-0.5,3.5, nStates+0.5, 3.5);
    TLine*line_v2 = new TLine(-0.5,11.5, nStates+0.5, 11.5);
    TLine*line_v3 = new TLine(-0.5,20.5, nStates+0.5, 20.5);
    TLine*line_v4 = new TLine(-0.5, nStates + 1 - 1.5, nStates + 0.5, nStates + 1 - 1.5);
    line_v1->SetLineColor(kBlack);
    line_v1->SetLineWidth(2);
    line_v2->SetLineColor(kBlack);
    line_v2->SetLineWidth(2);
    line_v3->SetLineColor(kBlack);
    line_v3->SetLineWidth(2);
    line_v4->SetLineColor(kBlack);
    line_v4->SetLineWidth(2);
    line_v1->Draw("same");
    line_v2->Draw("same");
    line_v3->Draw("same");
    line_v4->Draw("same");
    canvas->SaveAs((outDir + "/StatePairChannelCount_heatmap_" + run.label + ".pdf").c_str());
    canvas->SaveAs((outDir + "/StatePairChannelCount_heatmap_" + run.label + ".png").c_str());
    outputFile.cd();
    heatmap->Write();
  }
  outputFile.Close();
}

struct DirectParameterInfo {
  const char *branch;
  const char *axisTitle;
};

const std::vector<DirectParameterInfo> &directParameterInfos() {
  static const std::vector<DirectParameterInfo> parameters = {
      {"direct_threshold", "Direct fit threshold"},
      {"direct_threshold_error", "Direct fit threshold error"},
      {"direct_width", "Direct fit width"},
      {"direct_width_error", "Direct fit width error"},
      {"direct_chi2_ndf", "Direct fit chi2/NDF"},
      {"direct_chi2_ndf_bin4", "Direct fit chi2/NDF in 4-bin fit range"},
      {"direct_status", "Direct fit status"},
      {"direct_ok", "Direct fit OK"},
      {"direct_parameter_at_limit", "Direct fit parameter at limit"},
      {"direct_efficiency", "Direct fit efficiency"},
      {"direct_efficiency_range", "Direct fit efficiency in fit range"},
      {"direct_lg_scale", "Direct fit LG scale"},
      {"direct_target_area", "Direct fit target area"},
      {"direct_actual_area", "Direct fit actual area"},
      {"direct_x_min", "Direct fit x min"},
      {"direct_x_max", "Direct fit x max"},
      {"direct_nbins", "Direct fit number of bins"}};
  return parameters;
}

const DirectParameterInfo &efficiencyRatioParameterInfo() {
  static const DirectParameterInfo parameter = {"efficiency_ratio", "Efficiency ratio (data/MC)"};
  return parameter;
}

const DirectParameterInfo &refitChi2NdfParameterInfo() {
  static const DirectParameterInfo parameter = {"refit_chi2_ndf", "Refit chi2/NDF"};
  return parameter;
}

std::vector<DirectParameterInfo> directParameter2DInfos() {
  std::vector<DirectParameterInfo> parameters = directParameterInfos();
  parameters.push_back(efficiencyRatioParameterInfo());
  parameters.push_back(refitChi2NdfParameterInfo());
  return parameters;
}

std::pair<double, double> paddedRange(const std::vector<double> &values) {
  const auto limits = std::minmax_element(values.begin(), values.end());
  double low = *limits.first;
  double high = *limits.second;
  if (low == high) {
    const double padding = std::max(1.0, std::abs(low) * 0.05);
    return std::make_pair(low - padding, high + padding);
  }
  const double padding = 0.05 * (high - low);
  return std::make_pair(low - padding, high + padding);
}

const std::pair<const char *, const char *> kDirectParameterGroups[] = {
    {"success", "Direct threshold fit success"},
    {"fail_parameter_at_limit", "Direct threshold fit fail: parameter at limit"},
    {"fail_fit_status_not_zero", "Direct threshold fit fail: status != 0 or not OK"},
    {"fail_big_width", "Direct threshold fit fail: big width"}};

std::vector<double> collectDirectParameterValues(const std::vector<RunDecisionData> &runs,
                                                 const std::string &groupKey,
                                                 const std::string &branch) {
  std::vector<double> collected;
  for (const auto &run : runs) {
    const auto groupIt = run.directParameterGroups.find(groupKey);
    if (groupIt == run.directParameterGroups.end()) continue;
    for (const auto &values : groupIt->second) {
      const auto valueIt = values.find(branch);
      if (valueIt != values.end()) collected.push_back(valueIt->second);
    }
  }
  return collected;
}

void drawDirectParameter1D(const std::vector<RunDecisionData> &runs,
                           const std::string &outDir) {
  if (runs.empty()) return;

  gStyle->SetOptStat(0);
  const std::string plotDir = outDir + "/DirectFitParameter1D";
  gSystem->mkdir(plotDir.c_str(), true);

  for (const auto &group : kDirectParameterGroups) {
    const std::string groupKey = group.first;
    const std::string groupTitle = group.second;
    const std::string groupDir = plotDir + "/" + groupKey;
    gSystem->mkdir(groupDir.c_str(), true);
    TFile outputFile((groupDir + "/direct_fit_parameter_1D.root").c_str(), "RECREATE");

    for (const auto &parameter : directParameterInfos()) {
      const std::vector<double> values =
          collectDirectParameterValues(runs, groupKey, parameter.branch);
      if (values.empty()) continue;

      const std::pair<double, double> range = paddedRange(values);
      const std::string histName = std::string("h_") + groupKey + "_" + parameter.branch;
      const std::string title = groupTitle + " (N=" + std::to_string(values.size()) + ");" +
                                parameter.axisTitle + ";Channels";
      TH1D *hist = new TH1D(histName.c_str(), title.c_str(), 100, range.first, range.second);
      for (double value : values) hist->Fill(value);

      const std::string canvasName = "c_" + histName;
      TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1000, 850);
      canvas->SetLeftMargin(0.12);
      canvas->SetBottomMargin(0.12);
      hist->Draw("HIST");

      const std::string outputStem = groupDir + "/" + parameter.branch;
      canvas->SaveAs((outputStem + ".pdf").c_str());
      canvas->SaveAs((outputStem + ".png").c_str());

      canvas->SetLogy(true);
      hist->SetMinimum(0.5);
      hist->SetMaximum(std::max(5.0, hist->GetMaximum() * 10.0));
      canvas->Modified();
      canvas->Update();
      canvas->SaveAs((outputStem + "_logy.pdf").c_str());
      canvas->SaveAs((outputStem + "_logy.png").c_str());

      outputFile.cd();
      hist->Write();
      delete canvas;
      delete hist;
    }
    outputFile.Close();
  }
}

void drawDirectParameter2D(const std::vector<RunDecisionData> &runs,
                           const std::string &outDir) {
  if (runs.empty()) return;

  gStyle->SetOptStat(0);
  const std::string plotDir = outDir + "/DirectFitParameter2D";
  gSystem->mkdir(plotDir.c_str(), true);

  const std::vector<DirectParameterInfo> parameters = directParameter2DInfos();
  for (const auto &group : kDirectParameterGroups) {
    const std::string groupKey = group.first;
    const std::string groupTitle = group.second;
    const std::string groupDir = plotDir + "/" + groupKey;
    gSystem->mkdir(groupDir.c_str(), true);
    TFile outputFile((groupDir + "/direct_fit_parameter_2D.root").c_str(), "RECREATE");

    for (size_t xIndex = 0; xIndex < parameters.size(); ++xIndex) {
      for (size_t yIndex = xIndex + 1; yIndex < parameters.size(); ++yIndex) {
        const DirectParameterInfo &xParameter = parameters[xIndex];
        const DirectParameterInfo &yParameter = parameters[yIndex];
        std::vector<double> xValues;
        std::vector<double> yValues;

        for (const auto &run : runs) {
          const auto groupIt = run.directParameterGroups.find(groupKey);
          if (groupIt == run.directParameterGroups.end()) continue;
          for (const auto &values : groupIt->second) {
            const auto xIt = values.find(xParameter.branch);
            const auto yIt = values.find(yParameter.branch);
            if (xIt == values.end() || yIt == values.end()) continue;
            xValues.push_back(xIt->second);
            yValues.push_back(yIt->second);
          }
        }
        if (xValues.empty()) continue;

        const std::pair<double, double> xRange = paddedRange(xValues);
        const std::pair<double, double> yRange = paddedRange(yValues);
        const std::string histName = std::string("h_") + groupKey + "_" +
                                     yParameter.branch + "_vs_" + xParameter.branch;
        const std::string title = groupTitle + " (N=" + std::to_string(xValues.size()) + ");" +
                                  xParameter.axisTitle + ";" + yParameter.axisTitle;
        TH2D *hist = new TH2D(histName.c_str(), title.c_str(),
                              100, xRange.first, xRange.second,
                              100, yRange.first, yRange.second);
        for (size_t i = 0; i < xValues.size(); ++i) hist->Fill(xValues[i], yValues[i]);

        const std::string canvasName = "c_" + histName;
        TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1000, 850);
        canvas->SetRightMargin(0.15);
        canvas->SetLeftMargin(0.12);
        canvas->SetBottomMargin(0.12);
        hist->GetZaxis()->SetTitle("Channels");
        hist->Draw("COLZ");

        const std::string outputStem = groupDir + "/" +
                                       yParameter.branch + "_vs_" + xParameter.branch;
        canvas->SaveAs((outputStem + ".pdf").c_str());
        canvas->SaveAs((outputStem + ".png").c_str());
        outputFile.cd();
        hist->Write();
        delete canvas;
        delete hist;
      }
    }
    outputFile.Close();
  }
}

void drawEfficiencyDeficitByState(const std::vector<RunDecisionData> &runs,
                                  const std::string &outDir) {
  if (runs.empty()) return;

  gStyle->SetOptStat(0);
  const std::string plotDir = outDir + "/EfficiencyDeficit2D_byState";
  gSystem->mkdir(plotDir.c_str(), true);
  TFile outputFile((plotDir + "/efficiency_deficit_2D_by_state.root").c_str(), "RECREATE");

  for (const auto &state : stateInfos()) {
    std::vector<double> xValues;
    std::vector<double> yValues;
    for (const auto &run : runs) {
      const auto stateIt = run.efficiencyDeficitByState.find(state.bit);
      if (stateIt == run.efficiencyDeficitByState.end()) continue;
      for (const auto &point : stateIt->second) {
        xValues.push_back(point.first);
        yValues.push_back(point.second);
      }
    }
    if (xValues.empty()) continue;

    const std::string stateDir = plotDir + "/" + state.key;
    gSystem->mkdir(stateDir.c_str(), true);
    const std::string titlePrefix = std::string(state.group) + ": " + state.label +
                                    " (N=" + std::to_string(xValues.size()) + ")";
    const std::string axisTitles = ";1 - efficiency ratio;1 - direct efficiency";

    const std::pair<double, double> xRange = paddedRange(xValues);
    const std::pair<double, double> yRange = paddedRange(yValues);
    TH2D *full = new TH2D(
        (std::string("h_efficiency_deficit_") + state.key).c_str(),
        (titlePrefix + axisTitles).c_str(),
        100, xRange.first, xRange.second, 100, yRange.first, yRange.second);
    TH2D *zoom = new TH2D(
        (std::string("h_efficiency_deficit_zoom_") + state.key).c_str(),
        (titlePrefix + " (0-0.1 zoom)" + axisTitles).c_str(),
        100, 0.0, 0.1, 100, 0.0, 0.1);

    std::vector<double> positiveX;
    std::vector<double> positiveY;
    for (size_t i = 0; i < xValues.size(); ++i) {
      full->Fill(xValues[i], yValues[i]);
      zoom->Fill(xValues[i], yValues[i]);
      if (xValues[i] > 0.0 && yValues[i] > 0.0) {
        positiveX.push_back(xValues[i]);
        positiveY.push_back(yValues[i]);
      }
    }

    const auto drawAndSave = [&](TH2D *hist, const std::string &stem, bool logAxes) {
      TCanvas *canvas = new TCanvas(
          (std::string("c_") + hist->GetName()).c_str(), hist->GetName(), 1000, 850);
      canvas->SetRightMargin(0.15);
      canvas->SetLeftMargin(0.12);
      canvas->SetBottomMargin(0.12);
      canvas->SetLogx(logAxes);
      canvas->SetLogy(logAxes);
      hist->GetZaxis()->SetTitle("Channels");
      hist->Draw("COLZ");
      canvas->SaveAs((stateDir + "/" + stem + ".pdf").c_str());
      canvas->SaveAs((stateDir + "/" + stem + ".png").c_str());
      delete canvas;
    };

    drawAndSave(full, "one_minus_efficiency_ratio_vs_one_minus_direct_efficiency", false);
    drawAndSave(zoom, "one_minus_efficiency_ratio_vs_one_minus_direct_efficiency_zoom_0p1", false);
    outputFile.cd();
    full->Write();
    zoom->Write();
    delete full;
    delete zoom;

    if (!positiveX.empty()) {
      std::pair<double, double> logXRange = paddedRange(positiveX);
      std::pair<double, double> logYRange = paddedRange(positiveY);
      logXRange.first = std::max(logXRange.first, *std::min_element(positiveX.begin(), positiveX.end()) * 0.5);
      logYRange.first = std::max(logYRange.first, *std::min_element(positiveY.begin(), positiveY.end()) * 0.5);
      TH2D *logHist = new TH2D(
          (std::string("h_efficiency_deficit_logxlogy_") + state.key).c_str(),
          (titlePrefix + " (positive entries, log X/Y)" + axisTitles).c_str(),
          100, logXRange.first, logXRange.second, 100, logYRange.first, logYRange.second);
      for (size_t i = 0; i < positiveX.size(); ++i) logHist->Fill(positiveX[i], positiveY[i]);
      drawAndSave(logHist, "one_minus_efficiency_ratio_vs_one_minus_direct_efficiency_logxlogy", true);
      outputFile.cd();
      logHist->Write();
      delete logHist;
    }
  }
  outputFile.Close();
}

bool isUsableFitResult(double value) {
  // Optional branches use -999 as their missing-value sentinel.
  return std::isfinite(value) && value > -900.0;
}

std::set<std::string> fitResultParameterNames(const std::vector<RunDecisionData> &runs) {
  std::set<std::string> names;
  for (const auto &run : runs)
    for (const auto &cell : run.fitResultsByCell)
      for (const auto &value : cell.second) names.insert(value.first);
  return names;
}

void writeFitResultChanges(const std::vector<RunDecisionData> &runs,
                           const std::string &outDir) {
  const std::set<std::string> parameters = fitResultParameterNames(runs);
  std::ofstream detail((outDir + "/FitResultChanges_byRunCell.csv").c_str());
  detail << "run,previous_valid_run,cellid,parameter,value,previous_value,delta,relative_delta\n";

  // cellid -> parameter -> (last valid run, value)
  std::map<int, std::map<std::string, std::pair<std::string, double> > > previous;
  for (const auto &run : runs) {
    for (const auto &cell : run.fitResultsByCell) {
      for (const auto &value : cell.second) {
        if (!isUsableFitResult(value.second)) continue;
        detail << run.label << ",";
        const auto previousCell = previous.find(cell.first);
        const bool hasPrevious = previousCell != previous.end() &&
            previousCell->second.count(value.first);
        if (hasPrevious) {
          const std::pair<std::string, double> &old = previousCell->second.at(value.first);
          const double delta = value.second - old.second;
          detail << old.first << "," << cell.first << "," << value.first << ","
                 << std::setprecision(12) << value.second << "," << old.second << "," << delta
                 << ",";
          if (std::abs(old.second) > 1e-12) detail << delta / old.second;
          detail << "\n";
        } else {
          detail << "," << cell.first << "," << value.first << ","
                 << std::setprecision(12) << value.second << ",,,\n";
        }
        previous[cell.first][value.first] = std::make_pair(run.label, value.second);
      }
    }
  }

  std::ofstream summary((outDir + "/FitResultSummary_byRun.csv").c_str());
  summary << "run,parameter,valid_channels,mean,rms,min,max\n";
  for (const auto &run : runs) {
    for (const auto &parameter : parameters) {
      std::vector<double> values;
      for (const auto &cell : run.fitResultsByCell) {
        const auto it = cell.second.find(parameter);
        if (it != cell.second.end() && isUsableFitResult(it->second)) values.push_back(it->second);
      }
      if (values.empty()) continue;
      double sum = 0.0;
      double sum2 = 0.0;
      for (double value : values) {
        sum += value;
        sum2 += value * value;
      }
      const double mean = sum / values.size();
      const double rms = std::sqrt(std::max(0.0, sum2 / values.size() - mean * mean));
      summary << run.label << "," << parameter << "," << values.size() << ","
              << std::setprecision(12) << mean << "," << rms << ","
              << *std::min_element(values.begin(), values.end()) << ","
              << *std::max_element(values.begin(), values.end()) << "\n";
    }
  }
}

void drawFitResultTrends(const std::vector<RunDecisionData> &runs,
                         const std::string &outDir) {
  if (runs.empty()) return;
  const std::string plotDir = outDir + "/FitResultTrends";
  gSystem->mkdir(plotDir.c_str(), true);
  TFile outputFile((plotDir + "/fit_result_trends.root").c_str(), "RECREATE");
  const std::set<std::string> parameters = fitResultParameterNames(runs);

  for (const auto &parameter : parameters) {
    if (parameter.find("status") != std::string::npos) continue;
    std::vector<double> x, means, xErrors, rmsValues;
    for (size_t r = 0; r < runs.size(); ++r) {
      std::vector<double> values;
      for (const auto &cell : runs[r].fitResultsByCell) {
        const auto it = cell.second.find(parameter);
        if (it != cell.second.end() && isUsableFitResult(it->second)) values.push_back(it->second);
      }
      if (values.empty()) continue;
      double sum = 0.0, sum2 = 0.0;
      for (double value : values) { sum += value; sum2 += value * value; }
      const double mean = sum / values.size();
      x.push_back(r + 1.0);
      means.push_back(mean);
      xErrors.push_back(0.0);
      rmsValues.push_back(std::sqrt(std::max(0.0, sum2 / values.size() - mean * mean)));
    }
    if (x.empty()) continue;

    const std::string safe = fileSafeString(parameter);
    TCanvas canvas(("c_fit_trend_" + safe).c_str(), parameter.c_str(), 1200, 750);
    TGraphErrors graph(x.size(), &x[0], &means[0], &xErrors[0], &rmsValues[0]);
    graph.SetName(("g_" + safe + "_mean_rms").c_str());
    graph.SetTitle((parameter + " by run (error bars: channel RMS);Run;" + parameter).c_str());
    graph.SetMarkerStyle(20);
    graph.SetMarkerSize(1.1);
    graph.SetLineWidth(2);
    double yMin = means[0] - rmsValues[0];
    double yMax = means[0] + rmsValues[0];
    for (size_t i = 1; i < means.size(); ++i) {
      yMin = std::min(yMin, means[i] - rmsValues[i]);
      yMax = std::max(yMax, means[i] + rmsValues[i]);
    }
    const double padding = yMax > yMin ? 0.08 * (yMax - yMin) : 1.0;
    TH1D frame(("frame_" + safe).c_str(),
               (parameter + " by run (error bars: channel RMS);Run;" + parameter).c_str(),
               runs.size(), 0.5, runs.size() + 0.5);
    for (size_t r = 0; r < runs.size(); ++r)
      frame.GetXaxis()->SetBinLabel(r + 1, runs[r].label.c_str());
    frame.GetXaxis()->LabelsOption("v");
    frame.SetMinimum(yMin - padding);
    frame.SetMaximum(yMax + padding);
    frame.Draw("AXIS");
    graph.Draw("P SAME");
    canvas.SetBottomMargin(0.22);
    canvas.SetGridx();
    canvas.SetGridy();
    canvas.Modified();
    canvas.Update();
    canvas.SaveAs((plotDir + "/" + safe + "_byRun.pdf").c_str());
    canvas.SaveAs((plotDir + "/" + safe + "_byRun.png").c_str());
    outputFile.cd();
    graph.Write();
  }
  outputFile.Close();
}

std::vector<double> collectFitResults(const RunDecisionData &run,
                                      const std::string &parameter,
                                      int decision) {
  std::vector<double> values;
  for (const auto &cell : run.fitResultsByCell) {
    const auto decisionIt = run.decisionByCell.find(cell.first);
    if (decisionIt == run.decisionByCell.end() || decisionIt->second != decision) continue;
    const auto valueIt = cell.second.find(parameter);
    if (valueIt != cell.second.end() && isUsableFitResult(valueIt->second))
      values.push_back(valueIt->second);
  }
  return values;
}

int fitTrendRunColor(size_t index) {
  static const int colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1,
                               kOrange + 7, kCyan + 2, kViolet + 1, kTeal + 3,
                               kPink + 7, kAzure + 7, kSpring + 5, kGray + 2};
  return colors[index % (sizeof(colors) / sizeof(colors[0]))];
}

void writeFitResultSummaryByDecision(const std::vector<RunDecisionData> &runs,
                                     const std::string &outDir) {
  std::ofstream csv((outDir + "/FitResultSummary_byRunDecision.csv").c_str());
  csv << "run,decision,parameter,valid_channels,mean,rms,min,max\n";
  const std::set<std::string> parameters = fitResultParameterNames(runs);
  for (const auto &run : runs) {
    for (int decision = 0; decision < kNDecisionClasses; ++decision) {
      for (const auto &parameter : parameters) {
        const std::vector<double> values = collectFitResults(run, parameter, decision);
        if (values.empty()) continue;
        double sum = 0.0, sum2 = 0.0;
        for (double value : values) { sum += value; sum2 += value * value; }
        const double mean = sum / values.size();
        const double rms = std::sqrt(std::max(0.0, sum2 / values.size() - mean * mean));
        csv << run.label << "," << decisionKey(decision) << "," << parameter << ","
            << values.size() << "," << std::setprecision(12) << mean << "," << rms << ","
            << *std::min_element(values.begin(), values.end()) << ","
            << *std::max_element(values.begin(), values.end()) << "\n";
      }
    }
  }
}

void drawFitResultTrendsByDecision(const std::vector<RunDecisionData> &runs,
                                   const std::string &outDir) {
  if (runs.empty()) return;
  const std::string basePlotDir = outDir + "/FitResultTrendsByDecision";
  gSystem->mkdir(basePlotDir.c_str(), true);
  TFile outputFile((basePlotDir + "/fit_result_trends_by_decision.root").c_str(), "RECREATE");
  const std::set<std::string> parameters = fitResultParameterNames(runs);

  for (int decision = 0; decision < kNDecisionClasses; ++decision) {
    const std::string decisionDir = basePlotDir + "/" + decisionKey(decision);
    gSystem->mkdir(decisionDir.c_str(), true);
    for (const auto &parameter : parameters) {
      if (parameter.find("status") != std::string::npos) continue;
      std::vector<double> x, means, xErrors, rmsValues;
      for (size_t r = 0; r < runs.size(); ++r) {
        const std::vector<double> values = collectFitResults(runs[r], parameter, decision);
        if (values.empty()) continue;
        double sum = 0.0, sum2 = 0.0;
        for (double value : values) { sum += value; sum2 += value * value; }
        const double mean = sum / values.size();
        x.push_back(r + 1.0);
        means.push_back(mean);
        xErrors.push_back(0.0);
        rmsValues.push_back(std::sqrt(std::max(0.0, sum2 / values.size() - mean * mean)));
      }
      if (x.empty()) continue;

      const std::string safe = fileSafeString(parameter);
      const std::string objectStem = std::string(decisionKey(decision)) + "_" + safe;
      TCanvas canvas(("c_decision_trend_" + objectStem).c_str(), parameter.c_str(), 1200, 750);
      TGraphErrors graph(x.size(), &x[0], &means[0], &xErrors[0], &rmsValues[0]);
      graph.SetName(("g_" + objectStem + "_mean_rms").c_str());
      graph.SetMarkerStyle(20);
      graph.SetMarkerSize(1.1);
      graph.SetLineWidth(2);
      graph.SetLineColor(decisionColor(decision));
      graph.SetMarkerColor(decisionColor(decision));
      double yMin = means[0] - rmsValues[0], yMax = means[0] + rmsValues[0];
      for (size_t i = 1; i < means.size(); ++i) {
        yMin = std::min(yMin, means[i] - rmsValues[i]);
        yMax = std::max(yMax, means[i] + rmsValues[i]);
      }
      const double padding = yMax > yMin ? 0.08 * (yMax - yMin) : 1.0;
      TH1D frame(("frame_" + objectStem).c_str(),
                 (std::string(decisionName(decision)) + ": " + parameter +
                  " by run (error bars: channel RMS);Run;" + parameter).c_str(),
                 runs.size(), 0.5, runs.size() + 0.5);
      for (size_t r = 0; r < runs.size(); ++r)
        frame.GetXaxis()->SetBinLabel(r + 1, runs[r].label.c_str());
      frame.GetXaxis()->LabelsOption("v");
      frame.SetMinimum(yMin - padding);
      frame.SetMaximum(yMax + padding);
      frame.Draw("AXIS");
      graph.Draw("P SAME");
      canvas.SetBottomMargin(0.22);
      canvas.SetGrid();
      canvas.SaveAs((decisionDir + "/" + safe + "_byRun.pdf").c_str());
      canvas.SaveAs((decisionDir + "/" + safe + "_byRun.png").c_str());
      outputFile.cd();
      graph.Write();
    }
  }
  outputFile.Close();
}

std::pair<double, double> paddedHistogramRange(const std::vector<double> &values) {
  const double minimum = *std::min_element(values.begin(), values.end());
  const double maximum = *std::max_element(values.begin(), values.end());
  const double padding = maximum > minimum ? 0.03 * (maximum - minimum)
                                          : std::max(1.0, 0.05 * std::abs(minimum));
  return std::make_pair(minimum - padding, maximum + padding);
}

void drawFitResultDistributionsByDecision(const std::vector<RunDecisionData> &runs,
                                          const std::string &outDir) {
  if (runs.empty()) return;
  const std::string basePlotDir = outDir + "/FitResultDistributionsByDecision";
  gSystem->mkdir(basePlotDir.c_str(), true);
  TFile outputFile((basePlotDir + "/fit_result_distributions_by_decision.root").c_str(),
                   "RECREATE");
  const std::set<std::string> parameters = fitResultParameterNames(runs);
  for (int decision = 0; decision < kNDecisionClasses; ++decision) {
    const std::string decisionDir = basePlotDir + "/" + decisionKey(decision);
    gSystem->mkdir(decisionDir.c_str(), true);
    for (const auto &parameter : parameters) {
      if (parameter.find("status") != std::string::npos) continue;
      std::vector<std::vector<double> > valuesByRun(runs.size());
      std::vector<double> allValues;
      for (size_t r = 0; r < runs.size(); ++r) {
        valuesByRun[r] = collectFitResults(runs[r], parameter, decision);
        allValues.insert(allValues.end(), valuesByRun[r].begin(), valuesByRun[r].end());
      }
      if (allValues.empty()) continue;
      const std::pair<double, double> range = paddedHistogramRange(allValues);
      const std::string safe = fileSafeString(parameter);
      const std::string objectStem = std::string(decisionKey(decision)) + "_" + safe;
      TCanvas canvas(("c_distribution_" + objectStem).c_str(), parameter.c_str(), 1000, 800);
      TLegend legend(0.68, 0.62, 0.91, 0.89);
      legend.SetBorderSize(0);
      if (runs.size() > 6) legend.SetNColumns(2);
      std::vector<TH1D *> histograms;
      double maximum = 0.0;
      for (size_t r = 0; r < runs.size(); ++r) {
        if (valuesByRun[r].empty()) continue;
        TH1D *hist = new TH1D(("h_distribution_" + objectStem + "_run_" +
                               std::to_string(r)).c_str(), "", 60, range.first, range.second);
        hist->SetDirectory(nullptr);
        hist->SetStats(false);
        for (double value : valuesByRun[r]) hist->Fill(value);
        if (hist->Integral() > 0.0) hist->Scale(1.0 / hist->Integral(), "width");
        hist->SetLineColor(fitTrendRunColor(r));
        hist->SetLineWidth(2);
        maximum = std::max(maximum, hist->GetMaximum());
        histograms.push_back(hist);
        legend.AddEntry(hist, (runs[r].label + " (N=" +
                               std::to_string(valuesByRun[r].size()) + ")").c_str(), "l");
      }
      if (histograms.empty()) continue;
      histograms[0]->SetTitle((std::string(decisionName(decision)) + ": " + parameter +
                               " distributions;" + parameter + ";Density").c_str());
      histograms[0]->SetMaximum(maximum * 1.18);
      histograms[0]->Draw("HIST");
      for (size_t i = 1; i < histograms.size(); ++i) histograms[i]->Draw("HIST SAME");
      legend.Draw();
      canvas.SaveAs((decisionDir + "/" + safe + "_overlay.pdf").c_str());
      canvas.SaveAs((decisionDir + "/" + safe + "_overlay.png").c_str());
      outputFile.cd();
      for (TH1D *hist : histograms) { hist->Write(); delete hist; }
    }
  }
  outputFile.Close();
}

void drawFitResultDeltaFromFirstRunByDecision(const std::vector<RunDecisionData> &runs,
                                              const std::string &outDir) {
  if (runs.size() < 2) return;
  const std::string basePlotDir = outDir + "/FitResultDeltaFromFirstRunByDecision";
  gSystem->mkdir(basePlotDir.c_str(), true);
  TFile outputFile((basePlotDir + "/fit_result_delta_from_first_run_by_decision.root").c_str(),
                   "RECREATE");
  const std::set<std::string> parameters = fitResultParameterNames(runs);
  for (int decision = 0; decision < kNDecisionClasses; ++decision) {
    const std::string decisionDir = basePlotDir + "/" + decisionKey(decision);
    gSystem->mkdir(decisionDir.c_str(), true);
    for (const auto &parameter : parameters) {
      if (parameter.find("status") != std::string::npos) continue;
      std::vector<std::vector<double> > deltasByRun(runs.size());
      std::vector<double> allDeltas;
      for (size_t r = 1; r < runs.size(); ++r) {
        for (const auto &cell : runs[r].fitResultsByCell) {
          const auto decisionIt = runs[r].decisionByCell.find(cell.first);
          if (decisionIt == runs[r].decisionByCell.end() || decisionIt->second != decision) continue;
          const auto currentIt = cell.second.find(parameter);
          const auto referenceCell = runs[0].fitResultsByCell.find(cell.first);
          if (currentIt == cell.second.end() || referenceCell == runs[0].fitResultsByCell.end())
            continue;
          const auto referenceDecisionIt = runs[0].decisionByCell.find(cell.first);
          if (referenceDecisionIt == runs[0].decisionByCell.end() ||
              referenceDecisionIt->second != decision)
            continue;
          const auto referenceIt = referenceCell->second.find(parameter);
          if (referenceIt == referenceCell->second.end() ||
              !isUsableFitResult(currentIt->second) || !isUsableFitResult(referenceIt->second))
            continue;
          deltasByRun[r].push_back(currentIt->second - referenceIt->second);
        }
        allDeltas.insert(allDeltas.end(), deltasByRun[r].begin(), deltasByRun[r].end());
      }
      if (allDeltas.empty()) continue;
      std::pair<double, double> range = paddedHistogramRange(allDeltas);
      range.first = std::min(range.first, 0.0);
      range.second = std::max(range.second, 0.0);
      const std::string safe = fileSafeString(parameter);
      const std::string objectStem = std::string(decisionKey(decision)) + "_" + safe;
      TCanvas canvas(("c_delta_" + objectStem).c_str(), parameter.c_str(), 1000, 800);
      TLegend legend(0.65, 0.62, 0.91, 0.89);
      legend.SetHeader(("Reference: " + runs[0].label + " (same decision)").c_str());
      legend.SetBorderSize(0);
      if (runs.size() > 7) legend.SetNColumns(2);
      std::vector<TH1D *> histograms;
      double maximum = 0.0;
      for (size_t r = 1; r < runs.size(); ++r) {
        if (deltasByRun[r].empty()) continue;
        TH1D *hist = new TH1D(("h_delta_" + objectStem + "_run_" +
                               std::to_string(r)).c_str(), "", 60, range.first, range.second);
        hist->SetDirectory(nullptr);
        hist->SetStats(false);
        for (double delta : deltasByRun[r]) hist->Fill(delta);
        if (hist->Integral() > 0.0) hist->Scale(1.0 / hist->Integral(), "width");
        hist->SetLineColor(fitTrendRunColor(r));
        hist->SetLineWidth(2);
        maximum = std::max(maximum, hist->GetMaximum());
        histograms.push_back(hist);
        legend.AddEntry(hist, (runs[r].label + " (N=" +
                               std::to_string(deltasByRun[r].size()) + ")").c_str(), "l");
      }
      if (histograms.empty()) continue;
      histograms[0]->SetTitle((std::string(decisionName(decision)) + ": #Delta " + parameter +
                               " from first run;#Delta " + parameter + ";Density").c_str());
      histograms[0]->SetMaximum(maximum * 1.18);
      histograms[0]->Draw("HIST");
      for (size_t i = 1; i < histograms.size(); ++i) histograms[i]->Draw("HIST SAME");
      TLine zero(0.0, 0.0, 0.0, maximum * 1.18);
      zero.SetLineStyle(2);
      zero.SetLineColor(kBlack);
      zero.Draw();
      legend.Draw();
      canvas.SaveAs((decisionDir + "/" + safe + "_delta_overlay.pdf").c_str());
      canvas.SaveAs((decisionDir + "/" + safe + "_delta_overlay.png").c_str());
      outputFile.cd();
      for (TH1D *hist : histograms) { hist->Write(); delete hist; }
    }
  }
  outputFile.Close();
}

std::pair<int, int> parseRunRange(const std::string &label) {
  std::vector<int> numbers;
  std::string digits;
  for (size_t i = 0; i <= label.size(); ++i) {
    const char c = i < label.size() ? label[i] : '\0';
    if (std::isdigit(static_cast<unsigned char>(c))) {
      digits += c;
    } else if (!digits.empty()) {
      numbers.push_back(std::atoi(digits.c_str()));
      digits.clear();
      if (numbers.size() == 2) break;
    }
  }
  if (numbers.empty()) return std::make_pair(-1, -1);
  if (numbers.size() == 1) return std::make_pair(numbers[0], numbers[0]);
  return std::make_pair(numbers[0], numbers[1]);
}

void writeFinalCalibrationTree(const std::vector<RunDecisionData> &runs,
                               const std::string &outDir) {
  const std::string outputPath = outDir + "/FinalMIPCalibration.root";
  TFile outputFile(outputPath.c_str(), "RECREATE");
  if (outputFile.IsZombie()) {
    std::cerr << "Failed to create final calibration file: " << outputPath << std::endl;
    return;
  }

  TTree tree("mip_calibration", "Final MIP calibration values selected by fit decision");
  int runStart = -1, runEnd = -1;
  int cellid = -1, layer = -1, chip = -1, channel = -1, entries = 0;
  int decision = kInsufficientStatistics, state = 0;
  double mpv = -1.0, width = -1.0, gausSigma = -1.0;
  double mpvError = -1.0, widthError = -1.0, totalAreaError = -1.0, gausSigmaError = -1.0;
  double chi2 = -1.0;
  int ndf = -1;
  double directThreshold = -1.0, directWidth = -1.0;
  double directThresholdError = -1.0, directWidthError = -1.0, directChi2NdfBin4 = -1.0;
  std::string runLabel, decisionLabel, stateLabels;

  tree.Branch("run_start", &runStart, "run_start/I");
  tree.Branch("run_end", &runEnd, "run_end/I");
  tree.Branch("run_label", &runLabel);
  tree.Branch("cellid", &cellid, "cellid/I");
  tree.Branch("layer", &layer, "layer/I");
  tree.Branch("chip", &chip, "chip/I");
  tree.Branch("channel", &channel, "channel/I");
  tree.Branch("entries", &entries, "entries/I");
  tree.Branch("decision", &decision, "decision/I");
  tree.Branch("decision_name", &decisionLabel);
  tree.Branch("state", &state, "state/I");
  tree.Branch("state_names", &stateLabels);
  tree.Branch("mpv", &mpv, "mpv/D");
  tree.Branch("width", &width, "width/D");
  tree.Branch("gaus_sigma", &gausSigma, "gaus_sigma/D");
  tree.Branch("direct_threshold", &directThreshold, "direct_threshold/D");
  tree.Branch("direct_width", &directWidth, "direct_width/D");
  tree.Branch("mpv_error", &mpvError, "mpv_error/D");
  tree.Branch("width_error", &widthError, "width_error/D");
  tree.Branch("TotalArea_error", &totalAreaError, "TotalArea_error/D");
  tree.Branch("gaus_sigma_error", &gausSigmaError, "gaus_sigma_error/D");
  tree.Branch("chi2", &chi2, "chi2/D");
  tree.Branch("ndf", &ndf, "ndf/I");
  tree.Branch("direct_threshold_error", &directThresholdError, "direct_threshold_error/D");
  tree.Branch("direct_width_error", &directWidthError, "direct_width_error/D");
  tree.Branch("direct_chi2_ndf_bin4", &directChi2NdfBin4, "direct_chi2_ndf_bin4/D");

  long long writtenEntries = 0;
  for (const auto &run : runs) {
    const std::pair<int, int> runRange = parseRunRange(run.label);
    runStart = runRange.first;
    runEnd = runRange.second;
    runLabel = run.label;
    if (runStart < 0 || runEnd < 0) {
      std::cerr << "  Warning: could not parse run_start/run_end from '" << run.label
                << "'; storing -1" << std::endl;
    }
    for (const auto &entry : run.finalCalibrationEntries) {
      cellid = entry.cellid;
      layer = entry.layer;
      chip = entry.chip;
      channel = entry.channel;
      entries = entry.entries;
      decision = entry.decision;
      decisionLabel = decisionKey(entry.decision);
      state = entry.state;
      stateLabels = stateList(entry.state);
      mpv = entry.mpv;
      width = entry.width;
      gausSigma = entry.gaus_sigma;
      directThreshold = entry.direct_threshold;
      directWidth = entry.direct_width;
      mpvError = entry.mpv_error;
      widthError = entry.width_error;
      totalAreaError = entry.TotalArea_error;
      gausSigmaError = entry.gaus_sigma_error;
      chi2 = entry.chi2;
      ndf = entry.ndf;
      directThresholdError = entry.direct_threshold_error;
      directWidthError = entry.direct_width_error;
      directChi2NdfBin4 = entry.direct_chi2_ndf_bin4;
      tree.Fill();
      ++writtenEntries;
    }
  }
  outputFile.cd();
  tree.Write();
  outputFile.Close();
  std::cout << "Wrote " << writtenEntries << " entries to " << outputPath
            << " (tree: mip_calibration)" << std::endl;
}

// ============================================================================
// Main macro
// ============================================================================

void mip_fit_decision_by_run(const char *baseDir = ".",
                             const char *runList = "all",
                             const char *fileName = "mip_neighborcheck_fitted_sher_direct_5.root",
                             const char *treeName = "mip_fit_results",
                             const char *outDir = "mip_fit_decision_by_run_out_direct",
                             const char *knownEfficiencyDegradedCellIds = "",
                             bool writePerEntryCsv = false) {
  gSystem->Exec(Form("mkdir -p %s", outDir));

  std::vector<std::string> dirNames = parseDirectoryRanges(runList);
  if (dirNames.empty()) {
    std::cerr << "No directories to process" << std::endl;
    return;
  }

  const std::set<int> knownCellIds = parseCellIdList(knownEfficiencyDegradedCellIds);

  std::cout << "Loading " << dirNames.size() << " directories" << std::endl;
  std::cout << "  baseDir  = " << baseDir << std::endl;
  std::cout << "  fileName = " << fileName << std::endl;
  std::cout << "  treeName = " << treeName << std::endl;
  std::cout << "  outDir   = " << outDir << std::endl;
  if (!knownCellIds.empty()) {
    std::cout << "  known efficiency-degraded cellids: ";
    for (int cid : knownCellIds) std::cout << cid << " ";
    std::cout << std::endl;
  }

  std::vector<RunDecisionData> runs;
  for (const auto &dirName : dirNames) {
    const std::string path = std::string(baseDir) + "/" + dirName + "/" + fileName;
    if (gSystem->AccessPathName(path.c_str()) != 0) {
      std::cerr << "  Warning: file not accessible: " << path << std::endl;
      continue;
    }

    RunDecisionData data;
    if (processRun(baseDir, dirName, fileName, treeName, knownCellIds, outDir, writePerEntryCsv, data)) {
      int total = 0;
      for (int c = 0; c < kNDecisionClasses; ++c) total += data.counts[c];
      std::cout << "  " << dirName << ": " << total << " entries classified" << std::endl;
      runs.push_back(data);
    } else {
      std::cerr << "  Failed to process directory: " << dirName << std::endl;
    }
  }

  if (runs.empty()) {
    std::cerr << "No valid runs were processed" << std::endl;
    return;
  }

  writeFinalCalibrationTree(runs, outDir);
  writeSummaryTables(runs, outDir);
  writeStateSummaryTables(runs, outDir);
  drawDecisionStacks(runs, outDir);
  drawDecisionDetectorMaps(runs, outDir);
  drawStateSummaries(runs, outDir);
  drawStateLayerHeatmaps(runs, outDir);
  drawStatePairHeatmaps(runs, outDir);
  drawDirectParameter1D(runs, outDir);
  drawDirectParameter2D(runs, outDir);
  drawEfficiencyDeficitByState(runs, outDir);
  writeFitResultChanges(runs, outDir);
  drawFitResultTrends(runs, outDir);
  writeFitResultSummaryByDecision(runs, outDir);
  drawFitResultTrendsByDecision(runs, outDir);
  drawFitResultDistributionsByDecision(runs, outDir);
  drawFitResultDeltaFromFirstRunByDecision(runs, outDir);

  std::cout << "\nDone. Outputs written to: " << outDir << std::endl;
  std::cout << "  - DecisionSummary_byRun.csv" << std::endl;
  std::cout << "  - FinalMIPCalibration.root (TTree: mip_calibration)" << std::endl;
  std::cout << "  - DecisionSummary_byRun.txt" << std::endl;
  std::cout << "  - DecisionStack_counts_byRun.pdf/png" << std::endl;
  std::cout << "  - DecisionStack_fraction_byRun.pdf/png" << std::endl;
  std::cout << "  - DecisionDetectorMap/DecisionDetectorMap_<run>.pdf/png" << std::endl;
  std::cout << "  - DecisionDetectorMap/decision_detector_maps.root" << std::endl;
  std::cout << "  - StateSummary_byRun.csv/txt" << std::endl;
  std::cout << "  - StateSummary_byRunLayer.csv" << std::endl;
  std::cout << "  - StateOverlapMultiplicity_byRun.csv" << std::endl;
  std::cout << "  - StateOverlapCombinations_byRun.csv" << std::endl;
  std::cout << "  - StateChannelCount_heatmap_byRun.pdf/png" << std::endl;
  std::cout << "  - StateChannelCount_heatmap_byLayer_<run>.pdf/png" << std::endl;
  std::cout << "  - StatePairChannelCount_heatmap_<run>.pdf/png" << std::endl;
  std::cout << "  - StatePairChannelCount_heatmaps.root" << std::endl;
  std::cout << "  - StateOverlapMultiplicity_fraction_byRun.pdf/png" << std::endl;
  std::cout << "  - DirectFitParameter1D/<success-or-fail-reason>/*.{pdf,png,root}" << std::endl;
  std::cout << "  - DirectFitParameter2D/<success-or-fail-reason>/*.pdf/png/root" << std::endl;
  std::cout << "  - EfficiencyDeficit2D_byState/<state>/*.pdf/png/root" << std::endl;
  std::cout << "  - FitResultSummary_byRun.csv" << std::endl;
  std::cout << "  - FitResultChanges_byRunCell.csv" << std::endl;
  std::cout << "  - FitResultTrends/<parameter>_byRun.pdf/png" << std::endl;
  std::cout << "  - FitResultSummary_byRunDecision.csv" << std::endl;
  std::cout << "  - FitResultTrendsByDecision/<decision>/<parameter>_byRun.pdf/png" << std::endl;
  std::cout << "  - FitResultDistributionsByDecision/<decision>/<parameter>_overlay.pdf/png" << std::endl;
  std::cout << "  - FitResultDeltaFromFirstRunByDecision/<decision>/<parameter>_delta_overlay.pdf/png"
            << std::endl;
}

int main(int argc, char *argv[]) {
  std::string baseDir = ".";
  std::string runList = "all";
  std::string fileName = "mip_neighborcheck_fitted_sher_direct_8_full.root";
  std::string treeName = "mip_fit_results";
  std::string outDir = "full_fit";
  std::string knownCellIds = "";
  std::string selectCondition = "";
  int selectStateMask = 0;
  bool exportOnly = false;
  bool writePerEntryCsv = true;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--baseDir" && i + 1 < argc) baseDir = argv[++i];
    else if (arg == "--runList" && i + 1 < argc) runList = argv[++i];
    else if (arg == "--fileName" && i + 1 < argc) fileName = argv[++i];
    else if (arg == "--treeName" && i + 1 < argc) treeName = argv[++i];
    else if (arg == "--outDir" && i + 1 < argc) outDir = argv[++i];
    else if (arg == "--knownCellIds" && i + 1 < argc) knownCellIds = argv[++i];
    else if (arg == "--selectStateMask" && i + 1 < argc)
      selectStateMask = static_cast<int>(std::strtol(argv[++i], nullptr, 0));
    else if (arg == "--selectCondition" && i + 1 < argc) selectCondition = argv[++i];
    else if (arg == "--exportOnly") exportOnly = true;
    else if (arg == "--writePerEntryCsv") writePerEntryCsv = true;
    else if (arg == "--help" || arg == "-h") {
      std::cout << "Usage: " << argv[0] << " [options]\n"
                << "  --baseDir <path>       Base directory containing run subdirectories [.]\n"
                << "  --runList <list>       Comma-separated run directories or all [all]\n"
                << "  --fileName <name>      ROOT file name in each run directory [mip_neighborcheck_fitted_8.root]\n"
                << "  --treeName <name>      TTree name [mip_fit_results]\n"
                << "  --outDir <path>        Output directory [mip_fit_decision_by_run_out_8]\n"
                << "  --knownCellIds <csv>   CellIDs for known efficiency degradation, e.g. 210001,210002\n"
                << "  --selectStateMask <n>  Export channels containing all state bits; decimal or 0x...\n"
                << "  --selectCondition <e>  TTreeFormula condition, e.g. direct_chi2_ndf>3\n"
                << "  --exportOnly           Skip summary plots and only export selected channels\n"
                << "  --writePerEntryCsv     Also write DecisionEntries_<run>.csv\n";
      return 0;
    }
  }

  if (!exportOnly) {
    mip_fit_decision_by_run(baseDir.c_str(), runList.c_str(), fileName.c_str(), treeName.c_str(),
                            outDir.c_str(), knownCellIds.c_str(), writePerEntryCsv);
  }
  if (selectStateMask != 0 && !selectCondition.empty()) {
    exportChannelsByStateAndCondition(
        baseDir.c_str(), runList.c_str(), fileName.c_str(), treeName.c_str(), outDir.c_str(),
        selectStateMask, selectCondition.c_str(), knownCellIds.c_str());
  } else if (selectStateMask != 0 || !selectCondition.empty()) {
    std::cerr << "Both --selectStateMask and --selectCondition are required for channel export"
              << std::endl;
    return 1;
  }
  return 0;
}
