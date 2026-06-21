#include "PedestalReader.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <stdexcept>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <cstdio>
#include <filesystem>
#include <map>
#include <thread>
#include <utility>
#include <vector>
#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include <yaml-cpp/yaml.h>
#include "common/Logger.hpp"
#include "CalibDBIO/Query.hpp"
#include "common/config/ParseRunConfig.hpp"
#include "common/AHCALGeometry.hpp"
#include <unistd.h>
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
    std::map<int, std::time_t> pedestalRunMap; // runNumber to timestamp
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
            if (!(iss >> run >> timestamp)) { continue; }  
            pedestalRunMap[run] = timestamp;
        }
        infile.close();
    } else {
        LOG_WARN("Could not open {} to read pedestal run information; rebuilding cache from DB queries", pedestal_runs_file_);
    }
    const std::filesystem::path cache_dir = pedestal_runs_file_ + ".d";
    std::error_code ec;
    std::filesystem::create_directories(cache_dir, ec);
    if (ec) {
        LOG_WARN("Could not create pedestal run cache directory {}: {}", cache_dir.string(), ec.message());
        ec.clear();
    }
    if (std::filesystem::exists(cache_dir, ec) && std::filesystem::is_directory(cache_dir, ec)) {
        for (const auto& entry : std::filesystem::directory_iterator(cache_dir, ec)) {
            if (ec) {
                LOG_WARN("Could not scan pedestal run cache directory {}: {}", cache_dir.string(), ec.message());
                break;
            }
            if (!entry.is_regular_file(ec)) {
                continue;
            }
            std::ifstream cache_file(entry.path());
            int run;
            std::time_t timestamp;
            if (cache_file >> run >> timestamp) {
                pedestalRunMap[run] = timestamp;
            }
        }
    }
    const int max_cached_run = pedestalRunMap.empty()
        ? 0
        : std::max_element(pedestalRunMap.begin(), pedestalRunMap.end(),
                           [](const auto& a, const auto& b){ return a.first < b.first; })->first;
    if (pedestalRunMap.empty() || max_cached_run < runNumber) {
        //load pedestal run info from the runinfo
        for (int r = runNumber - 1; r >= 20000; --r) {
            if (pedestalRunMap.find(r) != pedestalRunMap.end()) {
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
                        pedestalRunMap[r] = timestamp;
                        LOG_INFO("Found pedestal run {} with timestamp {} for physics run {}", r, timestamp, runNumber);
                        if (std::filesystem::exists(cache_dir, ec) && std::filesystem::is_directory(cache_dir, ec)) {
                            ec.clear();
                            const auto final_path = cache_dir / (std::to_string(r) + ".txt");
                            const auto tmp_path = cache_dir / (std::to_string(r) + ".txt.tmp." + std::to_string(::getpid()));
                            std::ofstream cache_out(tmp_path);
                            if (cache_out.is_open()) {
                                cache_out << r << " " << timestamp << "\n";
                                cache_out.close();
                                std::error_code rename_ec;
                                std::filesystem::rename(tmp_path, final_path, rename_ec);
                                if (rename_ec) {
                                    LOG_WARN("Could not update pedestal cache {}: {}", final_path.string(), rename_ec.message());
                                }
                            } else {
                                LOG_WARN("Could not write pedestal cache {}", tmp_path.string());
                            }
                        }
                    }
                }
            } catch (const std::exception& e) {
                LOG_ERROR("Error querying pedestal run {}: {}", r, e.what());
            }
        }
    }
    std::vector<std::pair<int, std::time_t>> pedestalRuns;
    pedestalRuns.reserve(pedestalRunMap.size());
    for (const auto& [run, timestamp] : pedestalRunMap) {
        pedestalRuns.emplace_back(run, timestamp);
    }
    std::sort(pedestalRuns.begin(), pedestalRuns.end(), [](const auto& a, const auto& b){ return a.second < b.second; });
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
