#include "TriggersReader.hpp"
#include <stdexcept>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include <yaml-cpp/yaml.h>
#include "common/Logger.hpp"
#include "CalibDBIO/Query.hpp"
#include "common/config/ParseRunConfig.hpp"
#include "common/AHCALGeometry.hpp"
namespace CalibDBIO {
TriggersReader::TriggersReader(int runNumber) : m_runNumber(runNumber) {
    // Constructor  
}

TriggersReader::~TriggersReader() {
    // Destructor
}

Triggers TriggersReader::getTriggers() {
    readTriggers(m_runNumber);
    return triggers_;
}

void TriggersReader::readTriggers(int runNumber) {
    try {
        nlohmann::json response = QueryRun(runNumber, "Triggers", -1, false, false);
        if (response.empty()) {
            LOG_ERROR("No triggers data found for run {}", runNumber);
            throw std::runtime_error("No triggers data found for run " + std::to_string(runNumber));
        }
        if (response.is_array() && response.size() > 0) {
            response = response[0]; // Take the first object if it's an array
        }
        if (!response.is_object()) {
            LOG_ERROR("Invalid triggers data format for run {}: expected an object", runNumber);
            throw std::runtime_error("Invalid triggers data format for run " + std::to_string(runNumber) + ": expected an object");
        }
        auto summary = response["Summary"];
        triggers_.AnyVetoed = summary["AnyVetoed"].get<int>();
        triggers_.NonVetoed = summary["NonVetoed"].get<int>();
        triggers_.NonPhysical = summary["NonPhysical"].get<int>();
        triggers_.Physical = summary["Physical"].get<int>();
        triggers_.PreVeto = summary["PreVeto"].get<int>();
        triggers_.PostVeto = summary["PostVeto"].get<int>();
        triggers_.Recorded = summary["Recorded"].get<int>();
        LOG_DEBUG("Successfully read triggers for run {}: AnyVetoed={}, NonVetoed={}, NonPhysical={}, Physical={}, PreVeto={}, PostVeto={}, Recorded={}", 
            runNumber, triggers_.AnyVetoed, triggers_.NonVetoed, triggers_.NonPhysical, triggers_.Physical, triggers_.PreVeto, triggers_.PostVeto, triggers_.Recorded);
    } catch (const std::exception& e) {
        LOG_ERROR("Error reading triggers for run {}: {}", runNumber, e.what());
        throw; // rethrow the exception after logging
    }
    LOG_INFO("Finished reading triggers for run {}", runNumber);
}
} // namespace CalibDBIO