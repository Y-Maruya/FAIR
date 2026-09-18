#pragma once

#include "common/RunContext.hpp"
#include "IO/Descriptor.hpp"

inline std::vector<FieldDesc> describe(const RunConfig*) {
    return {
        field("input", &RunConfig::input),
        field("output", &RunConfig::output),
        field("log_file", &RunConfig::log_file),
        field("runNumber", &RunConfig::runNumber),
        field("poolIndex", &RunConfig::poolIndex),
        field("MC", &RunConfig::MC),
        field("nEvents", &RunConfig::nEvents),
        field("log_level", &RunConfig::log_level),
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
        field("calibRate", &ConditionStore::calibRate),
    };
}
