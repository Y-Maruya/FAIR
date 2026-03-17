#pragma once
#include <yaml-cpp/yaml.h>
#include "common/RunContext.hpp"
#include "common/config/YAMLUtil.hpp"
#include <curl/curl.h>
#include <cstdio>
#include <ctime>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
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
inline std::string make_runinfo_url(int runNumber) {
  return "https://faser-runinfo.app.cern.ch/cgibin/getRunInfo.py?runno=" + std::to_string(runNumber);
}
inline void ensure_curl_initialized() {
  static const bool initialized = []() {
    const CURLcode code = curl_global_init(CURL_GLOBAL_DEFAULT);
    if (code != CURLE_OK) {
      throw std::runtime_error(std::string("curl_global_init failed: ") + curl_easy_strerror(code));
    }
    return true;
  }();
  (void)initialized;
}
inline size_t curl_write_to_string(char* ptr, size_t size, size_t nmemb, void* userdata) {
  auto* out = static_cast<std::string*>(userdata);
  out->append(ptr, size * nmemb);
  return size * nmemb;
}

inline std::string fetch_text_from_url(const std::string& url) {
  ensure_curl_initialized();

  std::unique_ptr<CURL, decltype(&curl_easy_cleanup)> handle(curl_easy_init(), &curl_easy_cleanup);
  if (!handle) {
    throw std::runtime_error("curl_easy_init failed");
  }

  std::string response;
  char error_buffer[CURL_ERROR_SIZE] = {0};

  curl_easy_setopt(handle.get(), CURLOPT_URL, url.c_str());
  curl_easy_setopt(handle.get(), CURLOPT_FOLLOWLOCATION, 1L);
  curl_easy_setopt(handle.get(), CURLOPT_FAILONERROR, 1L);
  curl_easy_setopt(handle.get(), CURLOPT_CONNECTTIMEOUT, 10L);
  curl_easy_setopt(handle.get(), CURLOPT_TIMEOUT, 30L);
  curl_easy_setopt(handle.get(), CURLOPT_USERAGENT, "FAIR/0.1");
  curl_easy_setopt(handle.get(), CURLOPT_WRITEFUNCTION, &curl_write_to_string);
  curl_easy_setopt(handle.get(), CURLOPT_WRITEDATA, &response);
  curl_easy_setopt(handle.get(), CURLOPT_ERRORBUFFER, error_buffer);

  const CURLcode code = curl_easy_perform(handle.get());
  if (code != CURLE_OK) {
    const std::string detail = error_buffer[0] != '\0' ? error_buffer : curl_easy_strerror(code);
    throw std::runtime_error("Failed to fetch URL '" + url + "': " + detail);
  }

  return response;
}
inline bool has_suffix(const std::string& s, const std::string& suffix) {
  return s.size() >= suffix.size() &&
         s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

inline std::string trim_copy(const std::string& s) {
  const auto first = s.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    return "";
  }
  const auto last = s.find_last_not_of(" \t\r\n");
  return s.substr(first, last - first + 1);
}

inline double parse_unix_time_from_string(const std::string& input) {
  std::string ts = trim_copy(input);
  if (ts.empty()) {
    throw std::runtime_error("Empty timestamp string");
  }

  int tz_offset_seconds = 3600; // default: CET = UTC+1

  if (has_suffix(ts, " CEST")) {
    tz_offset_seconds = 2 * 3600;
    ts.erase(ts.size() - 5);
  } else if (has_suffix(ts, " CET")) {
    tz_offset_seconds = 1 * 3600;
    ts.erase(ts.size() - 4);
  } else if (has_suffix(ts, " UTC")) {
    tz_offset_seconds = 0;
    ts.erase(ts.size() - 4);
  } else if (!ts.empty() && ts.back() == 'Z') {
    tz_offset_seconds = 0;
    ts.pop_back();
  }

  ts = trim_copy(ts);

  std::string base = ts;
  std::int64_t micros = 0;

  const auto dot_pos = ts.find('.');
  if (dot_pos != std::string::npos) {
    base = ts.substr(0, dot_pos);
    std::string frac = ts.substr(dot_pos + 1);

    if (frac.empty() || frac.find_first_not_of("0123456789") != std::string::npos) {
      throw std::runtime_error("Invalid fractional seconds in timestamp: " + input);
    }

    if (frac.size() > 6) {
      frac = frac.substr(0, 6);
    }
    while (frac.size() < 6) {
      frac.push_back('0');
    }

    micros = std::stoll(frac);
  }

  std::tm tm = {};
  char* end = strptime(base.c_str(), "%Y-%m-%dT%H:%M:%S", &tm);
  if (end == nullptr || *end != '\0') {
    throw std::runtime_error("Invalid timestamp format: " + input);
  }

  const std::time_t utc_seconds = timegm(&tm) - tz_offset_seconds;
  return static_cast<double>(utc_seconds) + static_cast<double>(micros) / 1'000'000.0;
}
inline YAML::Node load_json_from_url(const std::string& url) {
  const std::string payload = fetch_text_from_url(url);
  try {
    return YAML::Load(payload);
  } catch (const std::exception& ex) {
    throw std::runtime_error("Failed to parse JSON from '" + url + "': " + ex.what());
  }
}
inline ConditionStore parse_condition_store(const int runNumber){
    ConditionStore cs;
    try {
        const std::string url = make_runinfo_url(runNumber);
        const YAML::Node json = load_json_from_url(url);
        if (!has_node(json, "runnumber") || json["runnumber"].as<int>() != runNumber) {
            LOG_ERROR("RunContext::parse_condition_store: run number mismatch in response from {}: expected {}, got {}", url, runNumber, json["runnumber"].as<int>());
            throw std::runtime_error("Run number mismatch in response");
        }
        if (has_node(json, "starttime")) {
            std::string starttime_str = json["starttime"].as<std::string>();
            cs.starttime = parse_unix_time_from_string(starttime_str);
        }else {
            LOG_WARN("RunContext::parse_condition_store: 'starttime' not found in response from {}, setting to 0", url);
            cs.starttime = 0;
        }
        if (has_node(json, "stoptime")) {
            std::string endtime_str = json["stoptime"].as<std::string>();
            cs.endtime = parse_unix_time_from_string(endtime_str);
        }else {
            LOG_WARN("RunContext::parse_condition_store: 'stoptime' not found in response from {}, setting to 0", url);
            cs.endtime = 1000;
        }
        if (has_node(json, "configuration")) {
            const auto& config = json["configuration"];
            if (has_node(config, "components")) {
                const auto& components = config["components"];
                for (const auto& comp : components) {
                    if (has_node(comp, "name") && comp["name"].as<std::string>() == "tlureceiver00"){
                        if (has_node(comp, "modules")) {
                            const auto& modules = comp["modules"];
                            for (const auto& mod : modules) {
                                if (has_node(mod, "name") && mod["name"].as<std::string>() == "tlureceiver00") {
                                    if (has_node(mod, "settings")) {
                                        const auto& settings = mod["settings"];
                                        if (has_node(settings, "trigger_logic")) {
                                            cs.triggerLogic = settings["trigger_logic"].as<std::string>();
                                        }
                                        if (has_node(settings, "voltage_threshold")){
                                            cs.thresholds = settings["voltage_threshold"].as<std::vector<double>>();
                                        } 
                                        if (has_node(settings, "trigger_stretch")) {
                                            cs.triggerStretch = settings["trigger_stretch"].as<std::vector<int>>();
                                        }
                                        if (has_node(settings, "trigger_delay")) {
                                            cs.triggerDelay = settings["trigger_delay"].as<std::vector<int>>();
                                        }
                                        if (has_node(settings, "FASER_CalibRate")) {
                                            cs.calibRate = settings["FASER_CalibRate"].as<int>();
                                        }
                                    }else {
                                        LOG_WARN("RunContext::parse_condition_store: 'settings' not found for module 'tlureceiver00' in configuration from {}, skipping trigger logic and settings", url);
                                    }
                                }else {
                                    LOG_DEBUG("RunContext::parse_condition_store: skipping module without name 'tlureceiver00' in configuration from {}", url);
                                }
                            }
                        }else {
                            LOG_WARN("RunContext::parse_condition_store: skipping component 'tlureceiver00' without 'modules' in configuration from {}", url);
                        }
                    }else {
                        LOG_DEBUG("RunContext::parse_condition_store: skipping component without name 'tlureceiver00' in configuration from {}", url);
                    }
                }
            }else {
                LOG_WARN("RunContext::parse_condition_store: 'components' not found in 'configuration' from {}, skipping trigger logic and settings", url);
            }
        }else {
            LOG_WARN("RunContext::parse_condition_store: 'configuration' not found in response from {}, skipping trigger logic and settings", url);
        }
        LOG_INFO("Successfully parsed condition store for run number {} from {}", runNumber, url);
        LOG_INFO("ConditionStore for run number {}: starttime={}, stoptime={}, triggerLogic='{}', thresholds=[{}], triggerStretch=[{}], triggerDelay=[{}], calibRate={}", 
            runNumber, cs.starttime, cs.endtime, cs.triggerLogic, fmt::join(cs.thresholds, ", "), fmt::join(cs.triggerStretch, ", "), fmt::join(cs.triggerDelay, ", "), cs.calibRate);
    } catch (const std::exception& ex) {
        throw std::runtime_error("Failed to parse condition store for run number " + std::to_string(runNumber) + ": " + ex.what());
    }
    return cs;
}
