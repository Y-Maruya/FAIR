#pragma once
#include <yaml-cpp/yaml.h>
#include "common/RunContext.hpp"
#include "common/config/YAMLUtil.hpp"

inline RunConfig parse_run_config(const YAML::Node& root) {
    RunConfig rc;
    const auto& run = require_node(root, "run");
    rc.input     = require_string(run, "input");
    rc.log_file  = require_string(run, "log_file");
    if (has_node(run, "MC")) {
        rc.MC = run["MC"].as<bool>();
    }
    rc.output    = require_string(run, "output");
    if (has_node(run, "nEvents")){
        rc.nEvents = run["nEvents"].as<long long>();
    }
    rc.log_level = require_string(run, "log_level");
    if (has_node(run, "runNumber")) {
        rc.runNumber = run["runNumber"].as<int>();
    }
    if (has_node(run, "poolIndex")) {
        rc.poolIndex = run["poolIndex"].as<int>();
    }
    return rc;
}