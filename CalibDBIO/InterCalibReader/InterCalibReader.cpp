#include "InterCalibReader.hpp"
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
std::mutex interCalibCacheMutex;
std::unordered_map<int, std::weak_ptr<const HGLGRatioMap>> hglgMapCache;
std::unordered_map<int, std::weak_ptr<const HGSaturationMap>> hgSaturationMapCache;
}

InterCalibReader::InterCalibReader(int runNumber) {
    // Constructor 
    readHGLGRatios(runNumber);
}

InterCalibReader::~InterCalibReader() {
    // Destructor
}

HGLGRatio InterCalibReader::getHGLGRatio(int cellID) const {
    auto it = hglgMap_->find(cellID);
    if (it != hglgMap_->end()) {
        return it->second;
    } else {
        LOG_ERROR("HGLG ratio not found for cellID {}", cellID);
        throw std::runtime_error("HGLG ratio not found for cellID " + std::to_string(cellID));
    }
}

double InterCalibReader::getHGSaturationPoint(int cellID) const {
    auto it = hgSaturationMap_->find(cellID);
    if (it != hgSaturationMap_->end()) {
        return it->second;
    } else {
        LOG_ERROR("HG saturation point not found for cellID {}", cellID);
        throw std::runtime_error("HG saturation point not found for cellID " + std::to_string(cellID));
    }
}


void InterCalibReader::readHGLGRatios(int runNumber) {
    {
        std::lock_guard<std::mutex> lock(interCalibCacheMutex);
        auto cachedHGLG = hglgMapCache[runNumber].lock();
        auto cachedSaturation = hgSaturationMapCache[runNumber].lock();
        if (cachedHGLG && cachedSaturation) {
            hglgMap_ = std::move(cachedHGLG);
            hgSaturationMap_ = std::move(cachedSaturation);
            LOG_INFO("Reusing {} cached HGLG entries and {} cached HG saturation entries for run {}",
                     hglgMap_->size(), hgSaturationMap_->size(), runNumber);
            return;
        }
    }

    auto hglgMap = std::make_shared<HGLGRatioMap>();
    auto hgSaturationMap = std::make_shared<HGSaturationMap>();

    for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer) {
        LOG_DEBUG("Reading HGLG ratios for Layer {}", layer);
        // Query the calibration database for HGLG ratios of the selected run
        auto response = CalibDBIO::QueryRun(runNumber, "CorrectedIntercalib", layer, true, false);
        if (response.empty()) {
            LOG_ERROR("No HGLG ratio data found for run {}", runNumber);
            throw std::runtime_error("No HGLG ratio data found for run " + std::to_string(runNumber));
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
            int totalChannels = perChannel["Slope"].size();
            int MaskedChannels = 0;
            int correctedChannels = 0;
            for (int chch = 0; chch < totalChannels; ++chch) {
                HGLGRatio ratio;
                int cellid = AHCALGeometry::CellID(layer, chch / AHCALGeometry::channel_No, chch % AHCALGeometry::channel_No);
                ratio.slope = perChannel["Slope"][chch].get<double>();
                ratio.slopeerror = perChannel["SlopeError"][chch].get<double>();
                ratio.intercept = perChannel["Intercept"][chch].get<double>();
                ratio.intercepterror = perChannel["InterceptError"][chch].get<double>();
                ratio.qualityflag = perChannel["QualityFlag"][chch].get<int>();
                if (ratio.qualityflag == 1 ){
                    // Masked channel
                    ++MaskedChannels;
                }else if (ratio.qualityflag == 2){
                    // Corrected channel
                    ++correctedChannels;
                }
                (*hglgMap)[cellid] = ratio;
                double saturationPoint = perChannel["HG_SaturationPoint"][chch].get<double>();
                if (saturationPoint > 0) {
                    (*hgSaturationMap)[cellid] = saturationPoint;
                }
            }
            if (MaskedChannels > 0) {
                LOG_INFO("Layer {}: {} channels are masked", layer, MaskedChannels);
            }
            if (correctedChannels > 0) {
                LOG_INFO("Layer {}: {} channels are corrected", layer, correctedChannels);
            }
        } catch (const std::exception& e) {
            LOG_ERROR("Error parsing hglg data: {}", e.what());
            throw std::runtime_error("Error parsing hglg data: " + std::string(e.what()));
        }
    }
    hglgMap_ = std::move(hglgMap);
    hgSaturationMap_ = std::move(hgSaturationMap);
    {
        std::lock_guard<std::mutex> lock(interCalibCacheMutex);
        hglgMapCache[runNumber] = hglgMap_;
        hgSaturationMapCache[runNumber] = hgSaturationMap_;
    }
    LOG_INFO("Successfully read {} hglg entries for run {}", hglgMap_->size(), runNumber);
}

} // namespace CalibDBIO
