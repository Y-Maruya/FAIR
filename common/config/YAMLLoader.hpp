#pragma once
#include <yaml-cpp/yaml.h>
#include <string>
#include "common/YAMLUtils.hpp"
#include <filesystem>
#include <stdexcept>

static YAML::Node resolve_node(const YAML::Node& n,
                               const std::filesystem::path& base_dir);

YAML::Node LoadYAMLWithInclude(const std::string& filename) {
    std::filesystem::path p(filename);
    if (!std::filesystem::exists(p)) {
        LOG_ERROR("YAMLLoader::LoadYAMLWithInclude: file not found '{}'", filename);
        throw std::runtime_error("YAML not found: " + filename);
    }
    YAML::Node root = YAML::LoadFile(filename);
    return resolve_node(root, p.parent_path());
}

static YAML::Node resolve_node(const YAML::Node& n,
                               const std::filesystem::path& base_dir) {
    // Our app-defined include tag
    if (n.Tag() == "!include") {
        if (!n.IsScalar()) {
            LOG_ERROR("YAMLLoader::resolve_node: !include must be a scalar path");
            throw std::runtime_error("!include must be a scalar path");
        }
        const std::string rel = n.as<std::string>();
        std::filesystem::path inc_path = base_dir / rel;

        if (!std::filesystem::exists(inc_path)) {
            throw std::runtime_error("Included YAML not found: " + inc_path.string());
        }

        YAML::Node inc = YAML::LoadFile(inc_path.string());
        // important: base_dir becomes included file's dir
        return resolve_node(inc, inc_path.parent_path());
    }

    if (n.IsMap()) {
        YAML::Node out(YAML::NodeType::Map);
        for (auto it : n) out[it.first] = resolve_node(it.second, base_dir);
        return out;
    }

    if (n.IsSequence()) {
        YAML::Node out(YAML::NodeType::Sequence);
        for (auto it : n) out.push_back(resolve_node(it, base_dir));
        return out;
    }

    return n; // scalar/null
}