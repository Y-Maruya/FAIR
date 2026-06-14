#include "PedestalReader.hpp"
#include <stdexcept>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <cstdio>
#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include <yaml-cpp/yaml.h>
#include "common/Logger.hpp"
#include "CalibDBIO/Query.hpp"
#include "common/config/ParseRunConfig.hpp"
#include "common/AHCALGeometry.hpp"
namespace CalibDBIO {
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
            // Response can be either an object or an array with one object
            nlohmann::json data = response.is_array() && response.size() > 0 ? response[0] : response;
            
            if (!data.contains("PerChannel") || !data["PerChannel"].is_object()) {
                LOG_ERROR("Invalid pedestal data format: 'PerChannel' field missing or not an object");
                throw std::runtime_error("Invalid pedestal data format: 'PerChannel' field missing or not an object");
            }
            if (data["Layer"] != layer) {
                LOG_ERROR("Pedestal data layer mismatch: expected {}, got {}", layer, data["Layer"].get<int>());
                throw std::runtime_error("Pedestal data layer mismatch: expected " + std::to_string(layer) + ", got " + std::to_string(data["Layer"].get<int>()));
            }
            
            const auto& perChannel = data["PerChannel"];
            int totalChannels = perChannel["HighGainPeak"].size();
            
            for (int chch = 0; chch < totalChannels; ++chch) {
                Pedestal ped;
                int cellid = AHCALGeometry::CellID(layer, chch / AHCALGeometry::channel_No, chch % AHCALGeometry::channel_No);
                ped.HighGainPeak = perChannel["HighGainPeak"][chch].get<double>();
                ped.HighGainSigma = perChannel["HighGainSigma"][chch].get<double>();
                ped.HighGainStatus = perChannel["HighGainStatus"][chch].get<int>();
                ped.LowGainPeak = perChannel["LowGainPeak"][chch].get<double>();
                ped.LowGainSigma = perChannel["LowGainSigma"][chch].get<double>();
                ped.LowGainStatus = perChannel["LowGainStatus"][chch].get<int>();
                if ((ped.HighGainStatus != 0 && ped.HighGainStatus != 999) || (ped.LowGainStatus != 0 && ped.LowGainStatus != 999)) {
                    LOG_WARN("Pedestal fit status not OK for cellID {}: HG status {}, LG status {}", cellid, ped.HighGainStatus, ped.LowGainStatus);
                }
                pedestalMap_[cellid] = ped;
            }
        } catch (const std::exception& e) {
            LOG_ERROR("Error parsing pedestal data: {}", e.what());
            throw std::runtime_error("Error parsing pedestal data: " + std::string(e.what()));
        }
    }
    LOG_INFO("Successfully read {} pedestal entries for run {}", pedestalMap_.size(), pedestalRun);
}

int PedestalReader::selectNearestPedestalRun(int runNumber) {
    // If the previous run is a pedestal run, return it
    // If not select the nearest pedestal run based on the timestamp of the runs
    // make the table of pedestal runs and their timestamps at first. continue from the last content of the table if it is already made.
    // Then find the nearest pedestal run to the given runNumber based on the timestamp
    std::vector<std::pair<int, std::time_t>> pedestalRuns; // pair of (runNumber, timestamp)
    std::ifstream infile(pedestal_runs_file_, std::ios::in);
    int i = 0;
    int max_retries = 5;
    while (!infile.is_open() && i < max_retries) {
        i++;
        LOG_WARN("Could not open {} to read pedestal run information, retrying in 1 second...", pedestal_runs_file_);
        std::this_thread::sleep_for(std::chrono::seconds(1));
        infile.open(pedestal_runs_file_, std::ios::in);
    }
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
        for (int r = runNumber - 1; r >= 20000; --r) {
            if (std::find_if(pedestalRuns.begin(), pedestalRuns.end(), [r](const auto& p){ return p.first == r; }) != pedestalRuns.end()) {
                break; // already have this run in the table
            }
            try {
                auto response = CalibDBIO::QueryRun(r, "Pedestal", 0, true, false);
                LOG_DEBUG("Queried pedestal run {}: response size {}", r, response.size());
                LOG_DEBUG("Response content: {}", response.dump());
                if (!response.empty()) {
                    // Response can be either an object or an array with one object
                    nlohmann::json data = response.is_array() && response.size() > 0 ? response[0] : response;
                    if (data.is_object() && data.contains("TimeStamp")) {
                        std::string timestamp_string = data["TimeStamp"].get<std::string>();
                        struct std::tm tm = {};
                        strptime(timestamp_string.c_str(), "%Y-%m-%d %H:%M:%S", &tm);
                        std::time_t timestamp = mktime(&tm);
                        pedestalRuns.emplace_back(r, timestamp);
                        LOG_INFO("Found pedestal run {} with timestamp {} for physics run {}", r, timestamp, runNumber);
                        // break;
                    }
                }
            } catch (const std::exception& e) {
                LOG_ERROR("Error querying pedestal run {}: {}", r, e.what());
            }
        }
    }
    std::sort(pedestalRuns.begin(), pedestalRuns.end(), [](const auto& a, const auto& b){ return a.second < b.second; });
    // output the pedestal run table using atomic rename for thread-safety
    std::string rand = std::to_string(std::rand());
    std::string tmpFile = pedestal_runs_file_ + ".tmp"+"."+rand;
    std::ofstream outfile(tmpFile, std::ios::out);
    if (outfile.is_open()) {
        for (const auto& p : pedestalRuns) {
            outfile << p.first << " " << p.second << "\n";
        }
        outfile.close();
        // atomic rename to ensure readers always see complete file
        if (std::rename(tmpFile.c_str(), pedestal_runs_file_.c_str()) != 0) {
            // LOG_ERROR("Could not rename {} to {}", tmpFile, pedestal_runs_file_);
            LOG_ERROR("Could not rename {} to {}, but the pedestal run information is updated in {}", tmpFile, pedestal_runs_file_, tmpFile);
            // throw std::runtime_error("Could not rename " + tmpFile + " to " + pedestal_runs_file_);
        }
    } else {
        LOG_ERROR("Could not open {} to write pedestal run information", tmpFile);
        throw std::runtime_error("Could not open " + tmpFile + " to write pedestal run information");
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
        LOG_ERROR("Invalid start/end time for run {}, using current time {} as target time for pedestal run selection", runNumber, targetTime);
        throw std::runtime_error("Invalid start/end time for run " + std::to_string(runNumber) + ", cannot select pedestal run");
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
} // namespace CalibDBIO