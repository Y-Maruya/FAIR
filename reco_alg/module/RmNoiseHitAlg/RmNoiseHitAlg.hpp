#pragma once
#include <yaml-cpp/yaml.h>
#include "common/EventStore.hpp"
#include "common/edm/EDM.hpp"
#include "common/IAlg.hpp"
#include "common/RunContext.hpp"
#include "CalibDBIO/PedestalReader/PedestalReader.hpp"
#include <string>
#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace AHCALRecoAlg {
    class NoiseCut {
    public:
        explicit NoiseCut(const std::string& criteria) {
            std::string e = criteria;

            parser_.DefineVar("HG", &hg_);
            parser_.DefineVar("LG", &lg_);
            parser_.SetExpr(e);

            LOG_INFO("NoiseCut test: HG=500 LG=100 => {}", eval(500, 100));
            LOG_INFO("NoiseCut test: HG=1500 LG=100 => {}", eval(1500, 100));
            LOG_INFO("NoiseCut test: HG=500 LG=10 => {}", eval(500, 10));
            LOG_INFO("NoiseCut test: HG=500 LG=400 => {}", eval(500, 400));
        }

        NoiseCut(const NoiseCut&) = delete;
        NoiseCut& operator=(const NoiseCut&) = delete;

        bool eval(double hg, double lg) {
            hg_ = hg;
            lg_ = lg;

            try {
                return parser_.Eval() != 0.0;
            } catch (mu::Parser::exception_type& e) {
                LOG_ERROR("Error evaluating noise cut expression: {}", e.GetMsg());
                return false;
            }
        }

    private:
        mu::Parser parser_;
        double hg_{0.0};
        double lg_{0.0};
    };
    struct RmNoiseHitAlgCfg{
        std::string in_rawhit_key = "RawHits";
        std::string out_rawhit_key = "RmNoiseRawHits";
        std::string criteria = "";
        bool output_removed_hits = false;
        std::string out_removedrawhit_key = "NoiseRawHits";
        std::string out_summary_key = "RmNoiseSummary";
        bool hittag1_only = true;
    };

    class RmNoiseHitAlg final : public IAlg { // final to prevent inheritance
    public:
        RmNoiseHitAlg(RunContext& ctx, std::string name)
            : IAlg(ctx, std::move(name)) {
            }
        void parse_cfg(const YAML::Node& n) override;
        void execute(EventStore& evt) override; // Main execution method
        void initialize() override {
            LOG_INFO("Initializing RmNoiseHitAlg with criteria: {}", m_cfg.criteria);
            m_criteria = m_cfg.criteria;
            cut = std::make_unique<NoiseCut>(m_criteria);
            LOG_INFO("Noise cut criteria set to: {}", m_criteria);
            loadPedestals();
        }
        void init_by_run() override {
            LOG_INFO("RmNoiseHitAlg: init_by_run called, runNumber={}, starttime={}, endtime={}",
                     ctx().config.runNumber, ctx().conditions.starttime, ctx().conditions.endtime);
            loadPedestals();
        }
        void loadPedestals(){
            CalibDBIO::PedestalReader reader(ctx().config.runNumber);
            ped_map_ = reader.getPedestalMap();
            LOG_INFO("InterCalibAlg: loaded {} pedestal entries from DB", ped_map_.size());
        }
        void computeNoiseFlags(const std::vector<AHCALRawHit>& hits, std::vector<int>& isoFlag, bool hittag1_only);
    private:
        RmNoiseHitAlgCfg m_cfg;
        std::unique_ptr<NoiseCut> cut;
        std::string m_criteria;
        std::unordered_map<int, CalibDBIO::Pedestal> ped_map_;
    };
} // namespace AHCALRecoAlg