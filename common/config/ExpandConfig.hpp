#pragma once

#include <yaml-cpp/yaml.h>

#include <cstdlib>
#include <regex>
#include <stdexcept>
#include <string>

inline std::string expand_env_string(const std::string& input) {
  static const std::regex env_pattern(R"(\$\{([A-Za-z_][A-Za-z0-9_]*)\})");

  std::string output;
  std::size_t last_pos = 0;

  for (std::sregex_iterator it(input.begin(), input.end(), env_pattern), end; it != end; ++it) {
    const std::smatch& match = *it;
    output.append(input, last_pos, static_cast<std::size_t>(match.position()) - last_pos);

    const std::string var_name = match[1].str();
    const char* value = std::getenv(var_name.c_str());
    if (value == nullptr) {
      throw std::runtime_error(
          "Environment variable '" + var_name + "' is used in YAML but is not set.");
    }

    output += value;
    last_pos = static_cast<std::size_t>(match.position() + match.length());
  }

  output.append(input, last_pos, std::string::npos);
  return output;
}

inline void expand_env_in_yaml(YAML::Node node) {
  if (!node) {
    return;
  }

  if (node.IsScalar()) {
    const std::string value = node.as<std::string>();
    if (value.find("${") != std::string::npos) {
      node = expand_env_string(value);
    }
    return;
  }

  if (node.IsSequence()) {
    for (YAML::Node child : node) {
      expand_env_in_yaml(child);
    }
    return;
  }

  if (node.IsMap()) {
    for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
      expand_env_in_yaml(it->second);
    }
  }
}
