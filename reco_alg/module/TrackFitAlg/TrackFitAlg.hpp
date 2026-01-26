#ifndef AHCAL_RECO_ALG_TRACK_FIT_ALG_HPP
#define AHCAL_RECO_ALG_TRACK_FIT_ALG_HPP
#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include "common/edm/EDM.hpp"
#include <yaml-cpp/yaml.h>
#include <string>
#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace AHCALRecoAlg {
    struct TrackFitAlgCfg{
        std::string in_recohit_key = "RecoHits";
        std::string out_track_key = "FittedTrack";
        double threshold_xy = 40./2;
    };

    class TrackFitAlg final : public IAlg { // final to prevent inheritance
    public:
        TrackFitAlg(RunContext& ctx, std::string name)
            : IAlg(ctx, std::move(name)) {}

        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& n) override;
        void setThresholdXY(double threshold) {
            m_cfg.threshold_xy = threshold;
        }
    private:
        TrackFitAlgCfg m_cfg;
    };

} // namespace AHCALRecoAlg

#endif // AHCAL_RECO_ALG_TRACK_FIT_ALG_HPP