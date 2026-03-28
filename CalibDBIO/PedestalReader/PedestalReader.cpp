#include "PedestalReader.hpp"
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

PedestalReader::PedestalReader(int runNumber) {
    // Constructor  
    readPedestals(runNumber);
}

PedestalReader::~PedestalReader() {
    // Destructor
}

Pedestal PedestalReader::getPedestal(int cellID) {
    if (pedestalMap_.find(cellID) != pedestalMap_.end()) {
        return pedestalMap_[cellID];
    } else {
        LOG_ERROR("CellID not found in pedestal map: {}", cellID);
        throw std::runtime_error("CellID not found in pedestal map: " + std::to_string(cellID));
    }
}

void PedestalReader::readPedestals(int runNumber) {
    int pedestalRun = selectNearestPedestalRun(runNumber);
    LOG_INFO("Selected pedestal run {} for physics run {}", pedestalRun, runNumber);
    for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer) {
        LOG_DEBUG("Reading pedestal data for Layer {}", layer);
        // Query the calibration database for pedestal data of the selected run
        auto response = CalibDBIO::QueryRun(pedestalRun, "Pedestal", layer, true, false);
        if (response.empty()) {
            LOG_ERROR("No pedestal data found for run {}", pedestalRun);
            throw std::runtime_error("No pedestal data found for run " + std::to_string(pedestalRun));
        }

        // Parse the response and fill the pedestalMap_
        try {
            if (!response.contains("PerChannel") || !response["PerChannel"].is_array()) {
                LOG_ERROR("Invalid pedestal data format: 'PerChannel' field missing or not an array");
                throw std::runtime_error("Invalid pedestal data format: 'PerChannel' field missing or not an array");
            }
            if (response["Layer"] != layer) {
                LOG_ERROR("Pedestal data layer mismatch: expected {}, got {}", layer, response["Layer"].get<int>());
                throw std::runtime_error("Pedestal data layer mismatch: expected " + std::to_string(layer) + ", got " + std::to_string(response["Layer"].get<int>()));
            }
            for (const auto& entry : response["PerChannel"]) {
                Pedestal ped;
                for (int chch = 0; chch < AHCALGeometry::chip_No*AHCALGeometry::channel_No; ++chch) {
                    int cellid = AHCALGeometry::CellID(layer, chch / AHCALGeometry::channel_No, chch % AHCALGeometry::channel_No);
                    ped.HighGainPeak = entry["HighGainPeak"][chch].get<double>();
                    ped.HighGainSigma = entry["HighGainSigma"][chch].get<double>();
                    ped.HighGainStatus = entry["HighGainStatus"][chch].get<int>();
                    ped.LowGainPeak = entry["LowGainPeak"][chch].get<double>();
                    ped.LowGainSigma = entry["LowGainSigma"][chch].get<double>();
                    ped.LowGainStatus = entry["LowGainStatus"][chch].get<int>();
                    if ((ped.HighGainStatus != 0 && ped.HighGainStatus != 999) || (ped.LowGainStatus != 0 && ped.LowGainStatus != 999)) {
                        LOG_WARN("Pedestal fit status not OK for cellID {}: HG status {}, LG status {}", cellid, ped.HighGainStatus, ped.LowGainStatus);
                    }
                    pedestalMap_[cellid] = ped;
                }
            }
            LOG_INFO("Successfully read {} pedestal entries for run {}", pedestalMap_.size(), pedestalRun);
        } catch (const std::exception& e) {
            LOG_ERROR("Error parsing pedestal data: {}", e.what());
            throw std::runtime_error("Error parsing pedestal data: " + std::string(e.what()));
        }
    }
}

int PedestalReader::selectNearestPedestalRun(int runNumber) {
    // If the previous run is a pedestal run, return it
    // If not select the nearest pedestal run based on the timestamp of the runs
    // make the table of pedestal runs and their timestamps at first. continue from the last content of the table if it is already made.
    // Then find the nearest pedestal run to the given runNumber based on the timestamp
    std::vector<std::pair<int, std::time_t>> pedestalRuns; // pair of (runNumber, timestamp)
    std::ifstream infile(pedestal_runs_file_, std::ios::in);
    if (infile.is_open()) {
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            int run;
            std::time_t timestamp;
            if (!(iss >> run >> timestamp)) { break; }  
            pedestalRuns.emplace_back(run, timestamp);
        }
        infile.close();
    } else {
        LOG_ERROR("Could not open {} to read pedestal run information", pedestal_runs_file_);
        throw std::runtime_error("Could not open " + pedestal_runs_file_ + " to read pedestal run information");
    }
    if (pedestalRuns.empty() || pedestalRuns.back().first < runNumber) {
        //load pedestal run info from the runinfo
        for (int r = runNumber - 1; r >= 0; --r) {
            if (std::find_if(pedestalRuns.begin(), pedestalRuns.end(), [r](const auto& p){ return p.first == r; }) != pedestalRuns.end()) {
                break; // already have this run in the table
            }
            try {
                auto response = CalibDBIO::QueryRun(r, "Pedestal", 0, false, false);
                if (!response.empty()) {
                    std::time_t timestamp = response["TimeStamp"].get<std::time_t>();
                    pedestalRuns.emplace_back(r, timestamp);
                    LOG_INFO("Found pedestal run {} with timestamp {} for physics run {}", r, timestamp, runNumber);
                    break;
                }
            } catch (const std::exception& e) {
                LOG_ERROR("Error querying pedestal run {}: {}", r, e.what());
            }
        }
    }
    std::sort(pedestalRuns.begin(), pedestalRuns.end(), [](const auto& a, const auto& b){ return a.second < b.second; });
    // output the pedestal run table
    std::ofstream outfile(pedestal_runs_file_, std::ios::out);
    if (outfile.is_open()) {
        for (const auto& p : pedestalRuns) {
            outfile << p.first << " " << p.second << "\n";
        }
        outfile.close();
    } else {
        LOG_ERROR("Could not open {} to write pedestal run information", pedestal_runs_file_);
        throw std::runtime_error("Could not open " + pedestal_runs_file_ + " to write pedestal run information");
    }
    if (std::find_if(pedestalRuns.begin(), pedestalRuns.end(), [runNumber](const auto& p){ return p.first == runNumber - 1 ; }) != pedestalRuns.end()) {
        return runNumber - 1; // the previous run is a pedestal run
    }
    // find the nearest pedestal run based on the timestamp
    std::time_t targetTime = 0;
    YAML::Node runInfo = load_json_from_url(runinfo_url + std::to_string(runNumber));
    double startTime = 0.0;
    double endTime = 0.0;
    if (has_node(runInfo, "starttime")) {
        startTime = parse_unix_time_from_string(runInfo["starttime"].as<std::string>());
    } else {
        LOG_WARN("Start time not found in run info for run {}, using current time as fallback", runNumber);
    }
    if (has_node(runInfo, "stoptime")) {
        endTime = parse_unix_time_from_string(runInfo["stoptime"].as<std::string>());
    } else {
        LOG_WARN("End time not found in run info for run {}, using current time as fallback", runNumber);
    }
    if (startTime > 0.0 && endTime > 0.0) {
        targetTime = static_cast<std::time_t>(std::floor((startTime + endTime) / 2.0));
    } else {
        targetTime = std::time(nullptr);
    }
    auto it = std::min_element(pedestalRuns.begin(), pedestalRuns.end(), [targetTime](const auto& a, const auto& b){
        return std::abs(a.second - targetTime) < std::abs(b.second - targetTime);
    });
    if (it != pedestalRuns.end()) {
        LOG_INFO("Selected pedestal run {} with timestamp {} for physics run {}", it->first, it->second, runNumber);
        return it->first;
    } else {
        LOG_ERROR("No pedestal runs available to select for physics run {}", runNumber);
        throw std::runtime_error("No pedestal runs available to select for physics run " + std::to_string(runNumber));
    }
    return -1; // should never reach here
}