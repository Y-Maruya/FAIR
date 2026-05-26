#pragma once
#include <vector>
#include "RunContext.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"

namespace AHCALRecoAlg {
    struct ContextPutterConfig {
        std::string m_out_RunConfig_key = "RunConfig";
        std::string m_out_ConditionStore_key = "ConditionStore";
    };
    class ContextPutter final : public IAlg {
    public:
        ContextPutter(RunContext& ctx, std::string name)
            : IAlg(ctx, name){ }
        ~ContextPutter() override = default;

        void initialize() override {}
        void finalize() override {}
        void execute(EventStore& evt) override {
            RunConfig cfg = ctx().config;
            ConditionStore cs = ctx().conditions;
            evt.put(m_config.m_out_RunConfig_key, std::move(cfg));
            evt.put(m_config.m_out_ConditionStore_key, std::move(cs));
        }
        void parse_cfg(const YAML::Node& cfg) override {
            m_config.m_out_RunConfig_key = get_or<std::string>(cfg, "out_RunConfig_key", m_config.m_out_RunConfig_key);
            m_config.m_out_ConditionStore_key = get_or<std::string>(cfg, "out_ConditionStore_key", m_config.m_out_ConditionStore_key);
        }
    private:
        ContextPutterConfig m_config;
    };
}

AHCAL_REGISTER_ALG(AHCALRecoAlg::ContextPutter, "ContextPutter")
