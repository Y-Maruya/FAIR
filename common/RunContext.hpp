#pragma once
#include <string>
#include <yaml-cpp/yaml.h>
#include "common/config/YAMLUtil.hpp"
#include "common/Logger.hpp"
struct RunConfig {
  std::string input;                 // required
  std::string output = "out.root";
  std::string log_file = "app.log";
  int runNumber = 0;
  int poolIndex = 0;          // for multiple pooling writers
  bool MC = false;                  // is MC data?
  long long nEvents = -1;            // -1 = until EOF (if Reader supports)
  std::string log_level = "info";    // spdlog level name
};

struct ConditionStore {
  // Placeholder for condition data
  std::vector<int> skipLayers = {0,2,14};
};

struct RunContext {
  RunConfig config;
  ConditionStore conditions;
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