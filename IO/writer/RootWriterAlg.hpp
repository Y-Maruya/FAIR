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
    m_out.write_event_counter(evt.event_counter());
    update_run_context(evt.event_counter());

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
    m_last_event = evt.event_counter();
    m_out.fill();
  }

  void finalize() override {
    if (!m_has_active_context) return;
    const long long end_event = (m_last_event < m_context_start_event) ? m_context_start_event : m_last_event;
    m_out.write_run_context(m_active_context, m_context_start_event, end_event);
  }

  void parse_cfg(const YAML::Node& n) override {
    // No specific config for RootWriterAlg
    // suppress unused parameter warning
    (void)n;
  }

private:
  static bool same_run_config(const RunConfig& a, const RunConfig& b) {
    return a.input == b.input &&
           a.output == b.output &&
           a.log_file == b.log_file &&
           a.runNumber == b.runNumber &&
           a.poolIndex == b.poolIndex &&
           a.MC == b.MC &&
           a.nEvents == b.nEvents &&
           a.log_level == b.log_level;
  }

  static bool same_condition_store(const ConditionStore& a, const ConditionStore& b) {
    return a.skipLayers == b.skipLayers &&
           a.nTriggerLayers == b.nTriggerLayers &&
           a.triggerLayers == b.triggerLayers &&
           a.starttime == b.starttime &&
           a.endtime == b.endtime &&
           a.triggerLogic == b.triggerLogic &&
           a.thresholds == b.thresholds &&
           a.triggerStretch == b.triggerStretch &&
           a.triggerDelay == b.triggerDelay &&
           a.calibRate == b.calibRate;
  }

  static bool same_run_context(const RunContext& a, const RunContext& b) {
    return same_run_config(a.config, b.config) &&
           same_condition_store(a.conditions, b.conditions);
  }

  void update_run_context(long long event_counter) {
    if (!m_has_active_context) {
      m_active_context = ctx();
      m_context_start_event = event_counter;
      m_has_active_context = true;
      return;
    }

    if (!same_run_context(m_active_context, ctx())) {
      const long long end_event =
        (m_last_event < m_context_start_event) ? m_context_start_event : m_last_event;
      m_out.write_run_context(m_active_context, m_context_start_event, end_event);
      m_active_context = ctx();
      m_context_start_event = event_counter;
    }
  }

  RootOutput m_out;
  WriterRegistry m_reg;
  RunContext m_active_context{};
  long long m_context_start_event = 0;
  long long m_last_event = -1;
  bool m_has_active_context = false;
};
