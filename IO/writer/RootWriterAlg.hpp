// IO/RootWriterAlg.hpp
#pragma once
#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include "IO/writer/WriterRegistry.hpp"
#include "common/Logger.hpp"
#include <unordered_set>

class RootWriterAlg final : public IAlg {
public:
  RootWriterAlg(RunContext& ctx, std::string name, std::string filename, WriterRegistry reg)
    : IAlg(ctx, std::move(name)), m_out(std::move(filename)), m_reg(std::move(reg)) {}

  void execute(EventStore& evt) override {
    for (const auto& key : evt.keys()) {
      LOG_DEBUG("Processing key='{}'", key);
      const auto& a = evt.any(key);
      if (!m_reg.can_write(a)) {
        LOG_DEBUG("Skip key='{}' (type={})", key, a.type().name());
        continue;
      }
      LOG_DEBUG("Writing key='{}' (type={})", key, a.type().name());
      m_reg.write_any(key, a, m_out);
    }
    m_out.fill();
  }
  void parse_cfg(const YAML::Node& n) override {
    // Optional: keylist parameter can override/extend the whitelist
    // (keys are automatically extracted from outputlist [type, key] pairs)
    if (n["keylist"] && n["keylist"].IsSequence()) {
      std::unordered_set<std::string> keys;
      for (const auto& key_node : n["keylist"]) {
        keys.insert(key_node.as<std::string>());
      }
      m_reg.set_key_whitelist(keys);
      LOG_INFO("Configured RootWriterAlg with {} keys in whitelist (from keylist)", keys.size());
    }
    // If no keylist is specified, the whitelist set by outputlist [type, key] pairs is used
  }

private:
  RootOutput m_out;
  WriterRegistry m_reg;
};
