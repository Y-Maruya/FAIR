#pragma once
#include <yaml-cpp/yaml.h>
#include "common/RunContext.hpp"
#include "common/config/YAMLUtil.hpp"
#include <chrono>
#include <curl/curl.h>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstdint>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
inline RunConfig parse_run_config(const YAML::Node& root) {
    RunConfig rc;
    const auto& run = require_node(root, "run");
    rc.input     = require_string(run, "input");
    rc.log_file  = require_string(run, "log_file");
    if (has_node(run, "config_input_output")) {
        rc.config_input_output = run["config_input_output"].as<std::string>();
    }
    if (has_node(run, "config_output")) {
        rc.config_output = run["config_output"].as<std::string>();
    }
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

inline bool is_transient_fetch_failure(CURLcode code, long http_code) {
  if (code == CURLE_HTTP_RETURNED_ERROR) {
    return http_code == 429 || http_code >= 500;
  }
  return code == CURLE_COULDNT_RESOLVE_HOST ||
         code == CURLE_COULDNT_CONNECT ||
         code == CURLE_OPERATION_TIMEDOUT ||
         code == CURLE_SEND_ERROR ||
         code == CURLE_RECV_ERROR ||
         code == CURLE_GOT_NOTHING ||
         code == CURLE_PARTIAL_FILE;
}

inline std::string fetch_text_from_url(const std::string& url) {
  ensure_curl_initialized();

  constexpr int max_attempts = 20;
  for (int attempt = 1; attempt <= max_attempts; ++attempt) {
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
    if (code == CURLE_OK) {
      return response;
    }

    long http_code = 0;
    curl_easy_getinfo(handle.get(), CURLINFO_RESPONSE_CODE, &http_code);
    const std::string detail = error_buffer[0] != '\0' ? error_buffer : curl_easy_strerror(code);
    if (!is_transient_fetch_failure(code, http_code) || attempt == max_attempts) {
      throw std::runtime_error(
          "Failed to fetch URL '" + url + "' after " + std::to_string(attempt) +
          (attempt == 1 ? " attempt: " : " attempts: ") + detail);
    }

    const int delay_seconds = attempt;
    std::fprintf(stderr,
                 "Transient failure fetching '%s' (attempt %d/%d): %s. Retrying in %d second(s).\n",
                 url.c_str(), attempt, max_attempts, detail.c_str(), delay_seconds);
    std::this_thread::sleep_for(std::chrono::seconds(delay_seconds));
  }

  throw std::runtime_error("Failed to fetch URL '" + url + "'");
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

  int tz_offset_seconds = 0; // default: CET = UTC+1

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

inline std::string runinfo_snapshot_path() {
  if (const char* path = std::getenv("FAIR_RUNINFO_SNAPSHOT"); path != nullptr && path[0] != '\0') {
    return path;
  }
  if (const char* fair_base = std::getenv("FAIR_BASE"); fair_base != nullptr && fair_base[0] != '\0') {
    return std::string(fair_base) + "/config/local/runinfo_snapshot.json";
  }
  return "";
}

struct RunInfoSnapshot {
  bool available = false;
  YAML::Node root;
};

inline RunInfoSnapshot load_runinfo_snapshot(const std::string& path) {
  RunInfoSnapshot snapshot;
  if (path.empty()) {
    return snapshot;
  }
  try {
    snapshot.root = YAML::LoadFile(path);
    if (!snapshot.root.IsMap() ||
        !has_node(snapshot.root, "version") ||
        snapshot.root["version"].as<int>() != 1 ||
        !has_node(snapshot.root, "run_types") ||
        !snapshot.root["run_types"].IsMap() ||
        !has_node(snapshot.root, "details") ||
        !snapshot.root["details"].IsMap()) {
      throw std::runtime_error("expected version 1 with map nodes 'run_types' and 'details'");
    }
    snapshot.available = true;
  } catch (const std::exception& ex) {
    std::fprintf(stderr,
                 "RunInfo snapshot '%s' is unavailable: %s. Falling back to HTTP as needed.\n",
                 path.c_str(), ex.what());
  }
  return snapshot;
}

inline RunInfoSnapshot get_runinfo_snapshot() {
  const std::string path = runinfo_snapshot_path();
  static std::mutex mutex;
  static std::unordered_map<std::string, RunInfoSnapshot> snapshots;
  std::lock_guard<std::mutex> lock(mutex);
  const auto found = snapshots.find(path);
  if (found != snapshots.end()) {
    return found->second;
  }
  return snapshots.emplace(path, load_runinfo_snapshot(path)).first->second;
}

inline YAML::Node find_runinfo_detail_in_snapshot(int runNumber) {
  const RunInfoSnapshot snapshot = get_runinfo_snapshot();
  if (!snapshot.available) {
    return YAML::Node(YAML::NodeType::Undefined);
  }
  const YAML::Node detail = snapshot.root["details"][std::to_string(runNumber)];
  if (!detail || !detail.IsMap() || !has_node(detail, "runnumber") ||
      detail["runnumber"].as<int>() != runNumber) {
    return YAML::Node(YAML::NodeType::Undefined);
  }
  return detail;
}

inline std::optional<std::string> find_run_type_in_snapshot(int runNumber) {
  const RunInfoSnapshot snapshot = get_runinfo_snapshot();
  if (!snapshot.available) {
    return std::nullopt;
  }
  const YAML::Node run_type = snapshot.root["run_types"][std::to_string(runNumber)];
  if (!run_type || !run_type.IsScalar()) {
    return std::nullopt;
  }
  return run_type.as<std::string>();
}

struct RunInfoRecord {
  YAML::Node json;
  std::string source;
};

inline RunInfoRecord load_runinfo(int runNumber) {
  if (const YAML::Node detail = find_runinfo_detail_in_snapshot(runNumber); detail) {
    return {detail, "RunInfo snapshot '" + runinfo_snapshot_path() + "'"};
  }

  static std::mutex mutex;
  static std::unordered_map<int, YAML::Node> http_cache;
  {
    std::lock_guard<std::mutex> lock(mutex);
    const auto found = http_cache.find(runNumber);
    if (found != http_cache.end()) {
      return {found->second, "in-process HTTP cache"};
    }
  }

  const std::string url = make_runinfo_url(runNumber);
  const YAML::Node json = load_json_from_url(url);
  if (!has_node(json, "runnumber") || json["runnumber"].as<int>() != runNumber) {
    throw std::runtime_error("Run number mismatch in response from " + url);
  }
  {
    std::lock_guard<std::mutex> lock(mutex);
    http_cache[runNumber] = json;
  }
  return {json, url};
}

inline std::string get_run_type(int runNumber) {
  if (const auto run_type = find_run_type_in_snapshot(runNumber); run_type.has_value()) {
    return *run_type;
  }

  const RunInfoRecord record = load_runinfo(runNumber);
  if (!has_node(record.json, "runnumber") || record.json["runnumber"].as<int>() != runNumber) {
    throw std::runtime_error("Run number mismatch in response from " + record.source);
  }
  if (!has_node(record.json, "type")) {
    throw std::runtime_error("Missing run type in response from " + record.source);
  }
  return record.json["type"].as<std::string>();
}

inline ConditionStore parse_condition_store_from_runinfo(
    const YAML::Node& json, int runNumber, const std::string& source) {
    ConditionStore cs;
    if (!has_node(json, "runnumber") || json["runnumber"].as<int>() != runNumber) {
        const std::string actual = has_node(json, "runnumber")
            ? std::to_string(json["runnumber"].as<int>())
            : "<missing>";
        LOG_ERROR("RunContext::parse_condition_store: run number mismatch in response from {}: expected {}, got {}", source, runNumber, actual);
        throw std::runtime_error("Run number mismatch in response");
    }
    if (has_node(json, "starttime")) {
        std::string starttime_str = json["starttime"].as<std::string>();
        cs.starttime = parse_unix_time_from_string(starttime_str);
    }else {
        LOG_WARN("RunContext::parse_condition_store: 'starttime' not found in response from {}, setting to 0", source);
        cs.starttime = 0;
    }
    if (has_node(json, "stoptime")) {
        std::string endtime_str = json["stoptime"].as<std::string>();
        cs.endtime = parse_unix_time_from_string(endtime_str);
    }else {
        LOG_WARN("RunContext::parse_condition_store: 'stoptime' not found in response from {}, setting to 0", source);
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
                                    LOG_WARN("RunContext::parse_condition_store: 'settings' not found for module 'tlureceiver00' in configuration from {}, skipping trigger logic and settings", source);
                                }
                            }else {
                                LOG_DEBUG("RunContext::parse_condition_store: skipping module without name 'tlureceiver00' in configuration from {}", source);
                            }
                        }
                    }else {
                        LOG_WARN("RunContext::parse_condition_store: skipping component 'tlureceiver00' without 'modules' in configuration from {}", source);
                    }
                }else {
                    LOG_DEBUG("RunContext::parse_condition_store: skipping component without name 'tlureceiver00' in configuration from {}", source);
                }
            }
        }else {
            LOG_WARN("RunContext::parse_condition_store: 'components' not found in 'configuration' from {}, skipping trigger logic and settings", source);
        }
    }else {
        LOG_WARN("RunContext::parse_condition_store: 'configuration' not found in response from {}, skipping trigger logic and settings", source);
    }
    LOG_INFO("Successfully parsed condition store for run number {} from {}", runNumber, source);
    std::string thresholds_str = std::to_string(cs.thresholds[0]);
    std::string triggerStretch_str = std::to_string(cs.triggerStretch[0]);
    std::string triggerDelay_str = std::to_string(cs.triggerDelay[0]);
    for (size_t i = 1; i < cs.thresholds.size(); ++i) {
        thresholds_str += ", " + std::to_string(cs.thresholds[i]);
        triggerStretch_str += ", " + std::to_string(cs.triggerStretch[i]);
        triggerDelay_str += ", " + std::to_string(cs.triggerDelay[i]);
    }
    LOG_INFO("ConditionStore for run number {}: starttime={}, stoptime={}, triggerLogic='{}', thresholds=[{}], triggerStretch=[{}], triggerDelay=[{}], calibRate={}",
        runNumber, cs.starttime, cs.endtime, cs.triggerLogic, thresholds_str, triggerStretch_str, triggerDelay_str, cs.calibRate);
    return cs;
}

inline ConditionStore parse_condition_store(const int runNumber){
    try {
        const RunInfoRecord record = load_runinfo(runNumber);
        return parse_condition_store_from_runinfo(record.json, runNumber, record.source);
    } catch (const std::exception& ex) {
        throw std::runtime_error("Failed to parse condition store for run number " + std::to_string(runNumber) + ": " + ex.what());
    }
}
