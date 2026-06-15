#pragma once

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

#include <yaml-cpp/yaml.h>

#include "common/Logger.hpp"
#include "common/RunContext.hpp"

inline const std::string& config_execution_timestamp() {
    static const std::string timestamp = [] {
        const auto now = std::chrono::system_clock::now();
        const auto milliseconds =
            std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
        const std::time_t time = std::chrono::system_clock::to_time_t(now);
        std::tm local_time = {};
        localtime_r(&time, &local_time);

        std::ostringstream output;
        output << std::put_time(&local_time, "%Y%m%d_%H%M%S")
               << '_' << std::setfill('0') << std::setw(3) << milliseconds.count();
        return output.str();
    }();
    return timestamp;
}

inline std::filesystem::path config_output_directory(const RunConfig& cfg) {
    const std::filesystem::path log_path(cfg.log_file);
    return log_path.has_parent_path() ? log_path.parent_path() : std::filesystem::path(".");
}

inline std::filesystem::path input_config_path(const RunConfig& cfg) {
    if (!cfg.config_input_output.empty()) {
        return cfg.config_input_output;
    }
    return config_output_directory(cfg) /
           ("config_input_run" + std::to_string(cfg.runNumber) + "_" +
            config_execution_timestamp() + ".yaml");
}

inline std::filesystem::path resolved_config_path(const RunConfig& cfg) {
    if (!cfg.config_output.empty()) {
        return cfg.config_output;
    }
    return config_output_directory(cfg) /
           ("config_run" + std::to_string(cfg.runNumber) + "_" +
            config_execution_timestamp() + ".yaml");
}

inline void update_resolved_run_node(YAML::Node& root, const RunConfig& cfg) {
    YAML::Node run = root["run"];
    if (!run || !run.IsMap()) {
        run = root["run"] = YAML::Node(YAML::NodeType::Map);
    }

    run["input"] = cfg.input;
    run["output"] = cfg.output;
    run["log_file"] = cfg.log_file;
    if (!cfg.config_input_output.empty()) {
        run["config_input_output"] = cfg.config_input_output;
    }
    run["runNumber"] = cfg.runNumber;
    run["luminosity"] = cfg.luminosity;
    run["poolIndex"] = cfg.poolIndex;
    run["MC"] = cfg.MC;
    run["nEvents"] = cfg.nEvents;
    run["log_level"] = cfg.log_level;
    if (!cfg.config_output.empty()) {
        run["config_output"] = cfg.config_output;
    }
}

inline std::filesystem::path write_config_snapshot(const YAML::Node& config,
                                                   const std::filesystem::path& output_path,
                                                   const char* description) {
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    const std::filesystem::path temporary_path = output_path.string() + ".tmp";
    {
        std::ofstream output(temporary_path);
        if (!output) {
            throw std::runtime_error("Could not open config snapshot output: " +
                                     temporary_path.string());
        }
        output << config;
        if (!output) {
            throw std::runtime_error("Failed to write config snapshot: " +
                                     temporary_path.string());
        }
    }

    std::filesystem::rename(temporary_path, output_path);
    LOG_INFO("{} configuration written to {}", description, output_path.string());
    return output_path;
}

inline std::filesystem::path write_input_config(const YAML::Node& input_config,
                                                const RunConfig& run_config) {
    return write_config_snapshot(input_config, input_config_path(run_config), "Input");
}

inline std::filesystem::path write_resolved_config(const YAML::Node& parsed_config,
                                                   const RunConfig& run_config) {
    YAML::Node output_config = YAML::Clone(parsed_config);
    update_resolved_run_node(output_config, run_config);
    return write_config_snapshot(output_config, resolved_config_path(run_config), "Resolved");
}
