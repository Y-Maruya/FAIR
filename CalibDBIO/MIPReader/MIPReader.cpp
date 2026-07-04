#include "MIPReader.hpp"
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
#include <mutex>
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

namespace {
std::mutex mipCacheMutex;
std::unordered_map<int, std::weak_ptr<const MIPMap>> mipMapCache;
std::unordered_map<int, std::weak_ptr<const ThresholdMap>> thresholdMapCache;
std::unordered_map<int, std::weak_ptr<const MIPMPVMap>> mipMPVMapCache;
}

MIPReader::MIPReader(int runNumber) {
    // Constructor 
    readMIPs(runNumber);
}

MIPReader::~MIPReader() {
    // Destructor
}

MIP MIPReader::getMIP(int cellID) const {
    auto it = mipMap_->find(cellID);
    if (it != mipMap_->end()) {
        return it->second;
    } else {
        LOG_ERROR("MIP data not found for cellID {}", cellID);
        throw std::runtime_error("MIP data not found for cellID " + std::to_string(cellID));
    }
}

Threshold MIPReader::getThreshold(int cellID) const {
    auto it = thresholdMap_->find(cellID);
    if (it != thresholdMap_->end()) {
        return it->second;
    } else {
        LOG_ERROR("Threshold data not found for cellID {}", cellID);
        throw std::runtime_error("Threshold data not found for cellID " + std::to_string(cellID));
    }
}

void MIPReader::readMIPs(int runNumber) {
    {
        std::lock_guard<std::mutex> lock(mipCacheMutex);
        auto cachedMIP = mipMapCache[runNumber].lock();
        auto cachedThreshold = thresholdMapCache[runNumber].lock();
        auto cachedMPV = mipMPVMapCache[runNumber].lock();
        if (cachedMIP && cachedThreshold && cachedMPV) {
            mipMap_ = std::move(cachedMIP);
            thresholdMap_ = std::move(cachedThreshold);
            mipMPVMap_ = std::move(cachedMPV);
            LOG_INFO("Reusing {} cached MIP entries and {} cached threshold entries for run {}",
                     mipMap_->size(), thresholdMap_->size(), runNumber);
            return;
        }
    }

    auto mipMap = std::make_shared<MIPMap>();
    auto thresholdMap = std::make_shared<ThresholdMap>();
    auto mipMPVMap = std::make_shared<MIPMPVMap>();

    for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer) {
        LOG_DEBUG("Reading MIP data for Layer {}", layer);
        // Query the calibration database for MIP data of the selected run
        auto response = CalibDBIO::QueryRun(runNumber, "MIPImputed", layer, true, false);
        if (response.empty()) {
            LOG_ERROR("No MIP data found for run {}", runNumber);
            throw std::runtime_error("No MIP data found for run " + std::to_string(runNumber));
        }

        // Parse the response and fill the mipMap_
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
            int totalChannels = perChannel["Decision"].size();
            int lossyChannelsMPV = 0;
            int lossyChannelsThreshold = 0;
            int ImputeBitMPV = data["Summary"]["ImputedBitDefinitions"]["MPV"].get<int>();
            int FullFitBitMPV = data["Summary"]["FullFitBitDefinitions"]["MPV"].get<int>();
            int RefBitMPV = data["Summary"]["RefBitDefinitions"]["MPV"].get<int>();
            int ImputeBitThreshold = data["Summary"]["ImputedBitDefinitions"]["Threshold"].get<int>();
            int FullFitBitThreshold = data["Summary"]["FullFitBitDefinitions"]["Threshold"].get<int>();
            int RefBitThreshold = data["Summary"]["RefBitDefinitions"]["Threshold"].get<int>();
            for (int chch = 0; chch < totalChannels; ++chch) {
                MIP mip;
                int cellid = AHCALGeometry::CellID(layer, chch / AHCALGeometry::channel_No, chch % AHCALGeometry::channel_No);
                mip.mpv = perChannel["MPV"][chch].get<double>();
                mip.mpverror = perChannel["MPVError"][chch].get<double>();
                std::pair<double, double> mpvPair = {mip.mpv, mip.mpverror};
                mip.width = perChannel["Width"][chch].get<double>();
                mip.widtherror = perChannel["WidthError"][chch].get<double>();
                mip.gaussigma = perChannel["GausSigma"][chch].get<double>();
                mip.gaussigmaerror = perChannel["GausSigmaError"][chch].get<double>();
                mip.decision = perChannel["Decision"][chch].get<int>();
                mip.state = perChannel["State"][chch].get<int>();
                mip.imputed = perChannel["Imputed"][chch].get<int>() & ImputeBitMPV;
                mip.fullfit = perChannel["FullFit"][chch].get<int>() & FullFitBitMPV;
                mip.ref = perChannel["Ref"][chch].get<int>() & RefBitMPV;
                if ((mip.mpv == -999.0 || mip.mpverror == -999.0 )&& (mip.decision != kNoData)) {
                    ++lossyChannelsMPV;
                    continue; // Skip channels with invalid MIP values
                }
                (*mipMap)[cellid] = mip;
                (*mipMPVMap)[cellid] = mpvPair;
                Threshold threshold;
                threshold.threshold = perChannel["Threshold"][chch].get<double>();
                threshold.thresholderror = perChannel["ThresholdError"][chch].get<double>();
                threshold.thresholdwidth = perChannel["ThresholdWidth"][chch].get<double>();
                threshold.thresholdwidtherror = perChannel["ThresholdWidthError"][chch].get<double>();
                threshold.decision = perChannel["Decision"][chch].get<int>();
                threshold.state = perChannel["State"][chch].get<int>();
                threshold.imputed = perChannel["Imputed"][chch].get<int>() & ImputeBitThreshold;
                threshold.fullfit = perChannel["FullFit"][chch].get<int>() & FullFitBitThreshold;
                threshold.ref = perChannel["Ref"][chch].get<int>() & RefBitThreshold;
                if ((threshold.threshold == -999.0 || threshold.thresholderror == -999.0)) {
                    ++lossyChannelsThreshold;
                    continue; // Skip channels with invalid threshold values
                }
                (*thresholdMap)[cellid] = threshold;
            }
            if (lossyChannelsMPV > 0) {
                LOG_WARN("Found {} channels with non-zero MPV fit status in layer {}", lossyChannelsMPV, layer);
            }
            if (lossyChannelsThreshold > 0) {
                LOG_INFO("Found {} channels with non-zero threshold fit status in layer {}", lossyChannelsThreshold, layer);
            }
        } catch (const std::exception& e) {
            LOG_ERROR("Error parsing mpv data: {}", e.what());
            throw std::runtime_error("Error parsing mpv data: " + std::string(e.what()));
        }
    }
    mipMap_ = std::move(mipMap);
    thresholdMap_ = std::move(thresholdMap);
    mipMPVMap_ = std::move(mipMPVMap);
    {
        std::lock_guard<std::mutex> lock(mipCacheMutex);
        mipMapCache[runNumber] = mipMap_;
        thresholdMapCache[runNumber] = thresholdMap_;
        mipMPVMapCache[runNumber] = mipMPVMap_;
    }
    LOG_INFO("Successfully read {} mpv entries for run {}", mipMap_->size(), runNumber);
    LOG_INFO("Successfully read {} threshold entries for run {}", thresholdMap_->size(), runNumber);
}

} // namespace CalibDBIO
