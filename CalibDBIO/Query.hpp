#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include "common/Logger.hpp"

using json = nlohmann::json;
namespace CalibDBIO {
    inline std::string db_url = "https://ahcalib-calibrationdb.app.cern.ch/";
    // inline std::string db_url = "http://localhost:5000/";

    inline std::string trim_url(std::string value) {
        const auto first = value.find_first_not_of(" \t\n\r\f\v");
        if (first == std::string::npos) {
            return "";
        }
        const auto last = value.find_last_not_of(" \t\n\r\f\v");
        return value.substr(first, last - first + 1);
    }

    inline std::string normalized_db_url() {
        std::string value = trim_url(db_url);
        if (value.empty()) {
            LOG_ERROR("Calibration DB URL is empty");
            return "";
        }
        if (value.rfind("http://", 0) != 0 && value.rfind("https://", 0) != 0) {
            LOG_ERROR("Calibration DB URL must start with http:// or https://: {}", value);
            return "";
        }
        if (value.back() != '/') {
            value.push_back('/');
        }
        return value;
    }
    
    static size_t WriteCallback(void* contents, size_t size, size_t nmemb, void* userp) {
        size_t totalSize = size * nmemb;
        std::string* str = static_cast<std::string*>(userp);
        str->append(static_cast<char*>(contents), totalSize);
        return totalSize;
    }

    inline json QueryTime(std::string StartTime, std::string EndTime, std::string CalibrationType, int Layer, bool PerChannel = true, bool PerChip = false) {

        CURL* curl = curl_easy_init();
        if (!curl) {
            throw std::runtime_error("Failed to initialize CURL");
        }

        std::string response;
        std::ostringstream url;
        // URL-encode the input URL (only query parameters)
        const std::string base_url = normalized_db_url();
        if (base_url.empty()) {
            curl_easy_cleanup(curl);
            return json();
        }

        char* escStart = curl_easy_escape(curl, StartTime.c_str(), 0);
        char* escEnd   = curl_easy_escape(curl, EndTime.c_str(), 0);

        url << base_url << "Query"
            << "?CalibrationType=" << CalibrationType
            << "&startTime=" << (escStart ? escStart : "")
            << "&endTime="   << (escEnd ? escEnd : "")
            << "&Layer="     << Layer
            << "&PerChannel=" << (PerChannel ? "True" : "False")
            << "&PerChip=" << (PerChip ? "True" : "False");

        if (escStart) curl_free(escStart);
        if (escEnd) curl_free(escEnd);

        curl_easy_setopt(curl, CURLOPT_URL, url.str().c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);

        // HTTPS verification
        curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 1L);
        curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 2L);

        CURLcode res = curl_easy_perform(curl);

        long http_code = 0;
        curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);

        curl_easy_cleanup(curl);

        if (res != CURLE_OK) {
            LOG_ERROR("CURL request failed: {}", curl_easy_strerror(res));
            return json();
        }

        if (http_code != 200) {
            LOG_ERROR("HTTP error: {}", http_code);
            return json();
        }

        try {
            return json::parse(response);
        } catch (const std::exception& e) {
            LOG_ERROR("JSON parse error: {}", e.what());
            return json();
        }
    }

    inline json QueryRun(int RunNumber, std::string CalibrationType, int Layer, bool PerChannel = true, bool PerChip = false) {
        CURL* curl = curl_easy_init();
        if (!curl) {
            throw std::runtime_error("Failed to initialize CURL");
        }

        std::string response;
        std::ostringstream url;
        // URL-encode the input URL (only query parameters)
        const std::string base_url = normalized_db_url();
        if (base_url.empty()) {
            curl_easy_cleanup(curl);
            return json();
        }

        url << base_url << "Query"
            << "?CalibrationType=" << CalibrationType
            << "&RunNumber=" << RunNumber
            << "&Layer="     << Layer
            << "&PerChannel=" << (PerChannel ? "True" : "False")
            << "&PerChip=" << (PerChip ? "True" : "False");


        curl_easy_setopt(curl, CURLOPT_URL, url.str().c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);

        // HTTPS verification
        curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 1L);
        curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 2L);

        CURLcode res = curl_easy_perform(curl);

        long http_code = 0;
        curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);

        curl_easy_cleanup(curl);

        if (res != CURLE_OK) {
            LOG_ERROR("CURL request failed: {}", curl_easy_strerror(res));
            return json();
        }

        if (http_code != 200) {
            LOG_ERROR("HTTP error: {}", http_code);
            return json();
        }

        try {
            return json::parse(response);
        } catch (const std::exception& e) {
            LOG_ERROR("JSON parse error: {}", e.what());
            return json();
        }
    }
}
