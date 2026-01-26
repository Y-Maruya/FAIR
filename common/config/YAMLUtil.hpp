#pragma once
#include <yaml-cpp/yaml.h>
#include <string>
#include <stdexcept>
#include "common/Logger.hpp"

template <class T>
inline T get_or(const YAML::Node& n, const char* key, T def) {
    if (!n[key]){
        LOG_WARN("YAMLUtil::get_or: missing key '{}', using default", key);
        return def;
    }
    return n[key].as<T>();
}

inline std::string require_string(const YAML::Node& n, const char* key) {
    if (!n[key]){
        LOG_ERROR("YAMLUtil::require_string: missing key '{}'", key);
        throw std::runtime_error(std::string("Missing key: ") + key);
    }
    return n[key].as<std::string>();
}

inline const YAML::Node require_node(const YAML::Node& n, const char* key) {
    YAML::Node child = n[key]; 
    if (!child) {
        LOG_ERROR("YAMLUtil::require_node: missing key '{}'", key);
        throw std::runtime_error(std::string("Missing node: ") + key);
    }
    return child;
}

inline bool has_node(const YAML::Node& n, const char* key) {
    return static_cast<bool>(n[key]);
}
