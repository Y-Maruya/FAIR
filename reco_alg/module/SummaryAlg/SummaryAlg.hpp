#pragma once
#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include "common/edm/EDM.hpp"
#include <yaml-cpp/yaml.h>
#include <string>
#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <fstream>

namespace AHCALRecoAlg {
    struct SummaryAlgCfg{
        std::string in_recohit_key = "RecoHits";
        std::string out_summary_key = "EventSummary";
        bool output_summary_txt = true;
        std::string summary_txt_path = "EventSummary.txt";
    };

    class SummaryAlg final : public IAlg { // final to prevent inheritance
    public:
        SummaryAlg(RunContext& ctx, std::string name)
            : IAlg(ctx, std::move(name)) {}

        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& n) override;
        void initialize() override;
        void finalize() override;
    private:
        SummaryAlgCfg m_cfg;
        std::ofstream m_summary_txt_file;
        std::pair<Long64_t, int> m_max_nhit_event{0, 0}; // (event number, nHits)
        std::pair<Long64_t, double> m_max_energy_event{0, 0.0}; // (event number, TotalEnergy)
    };

} // namespace AHCALRecoAlg
