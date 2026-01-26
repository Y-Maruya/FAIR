// IO/RootWriterAlg.hpp
#pragma once
#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include "IO/writer/WriterRegistry.hpp"
#include "common/Logger.hpp"

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
    // No specific config for RootWriterAlg
    // suppress unused parameter warning
    (void)n;
  }

private:
  RootOutput m_out;
  WriterRegistry m_reg;
};
