#pragma once
#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include "common/edm/EDM.hpp"
#include <yaml-cpp/yaml.h>
#include <string>
#include <memory>
#include <vector>

namespace AHCALRecoAlg {
    struct SimMimicAlgCfg {
        std::string in_simhit_key;
        std::string out_recohit_key;
        std::string out_rawhit_key;
    };
    class SimMimicAlg : public IAlg {
    public:
        SimMimicAlg(RunContext& ctx, std::string name)
            : IAlg(ctx, std::move(name)) {}
        void initialize() override;
        void finalize() override;
        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& cfg) override;
    private:
        SimMimicAlgCfg cfg_;
    };
}