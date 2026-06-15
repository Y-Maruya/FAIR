#pragma once
#include <string>
#include <yaml-cpp/yaml.h>
#include "common/config/YAMLUtil.hpp"
#include "common/Logger.hpp"
#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"

struct RunConfig {
  std::string input;                 // required
  std::string output = "out.root";
  std::string log_file = "app.log";
  std::string config_input_output;     // empty = timestamped input config beside log_file
  std::string config_output;           // empty = timestamped resolved config beside log_file
  int runNumber = 0;
  double luminosity = 0;              // in fb^-1
  int poolIndex = 0;          // for multiple pooling writers
  bool MC = false;                  // is MC data?
  long long nEvents = -1;            // -1 = until EOF (if Reader supports)
  std::string log_level = "info";    // spdlog level name
};

struct ConditionStore {
  // Placeholder for condition data
  std::vector<int> skipLayers = {0,2,14,28};
  int nTriggerLayers = 4;
  std::vector<int> triggerLayers = {9,19,29,38};
  double starttime = 0;
  double endtime = 0;
  std::string triggerLogic = "";
  std::vector<double> thresholds = {0.1, 0.1, 0.1, 0.1, -0.075, -0.075};
  std::vector<int> triggerStretch = {3,3,3,3,20,20}; // in units of 6.25ns
  std::vector<int> triggerDelay = {2,2,2,2,0,0}; // in units of 6.25ns
  int calibRate = 0; 
};

struct RunContext {
  RunConfig config;
  ConditionStore conditions;
  bool operator!() const {
    return config.input.empty();
  }
};

// void parse_run_config(const YAML::Node& n, RunConfig& cfg) {
//   if (!has_node(n, "runinfo")){
//     LOG_ERROR("RunContext::parse_run_config: missing 'runinfo' node");
//     throw std::runtime_error("Missing 'runinfo' node in config");
//   }
//   const YAML::Node& runinfo = require_node(n, "runinfo");
//   cfg.input = require_string(runinfo, "input");
//   cfg.output = get_or<std::string>(runinfo,"output", "out.root");
//   cfg.runNumber = get_or<int>(runinfo,"runNumber", 0);
//   cfg.poolIndex = get_or<int>(runinfo,"poolIndex", 0);
//   cfg.MC = get_or<bool>(runinfo,"MC", false);
//   cfg.nEvents = get_or<long long>(runinfo,"nEvents", -1);
//   cfg.log_level = get_or<std::string>(runinfo,"log_level", "info");
// }

inline std::vector<FieldDesc> describe(const RunConfig*) {
  return {
    field("input", &RunConfig::input),
    field("output", &RunConfig::output),
    field("log_file", &RunConfig::log_file),
    field("config_input_output", &RunConfig::config_input_output),
    field("config_output", &RunConfig::config_output),
    field("runNumber", &RunConfig::runNumber),
    field("luminosity", &RunConfig::luminosity),
    field("poolIndex", &RunConfig::poolIndex),
    field("isMC", &RunConfig::MC),
    field("nEvents", &RunConfig::nEvents),
    field("log_level", &RunConfig::log_level)
  };
}

inline std::vector<FieldDesc> describe(const ConditionStore*) {
  return {
    field("skipLayers", &ConditionStore::skipLayers),
    field("nTriggerLayers", &ConditionStore::nTriggerLayers),
    field("triggerLayers", &ConditionStore::triggerLayers),
    field("starttime", &ConditionStore::starttime),
    field("endtime", &ConditionStore::endtime),
    field("triggerLogic", &ConditionStore::triggerLogic),
    field("thresholds", &ConditionStore::thresholds),
    field("triggerStretch", &ConditionStore::triggerStretch),
    field("triggerDelay", &ConditionStore::triggerDelay),
    field("calibRate", &ConditionStore::calibRate)
  };
}
AHCAL_REGISTER_IO_STRUCT(RunConfig, "RunConfig");
AHCAL_REGISTER_IO_STRUCT(ConditionStore, "ConditionStore");
