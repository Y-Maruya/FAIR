#pragma once
#include <yaml-cpp/yaml.h>
#include "common/EventStore.hpp"
#include "common/edm/EDM.hpp"
#include "common/IAlg.hpp"
#include "common/RunContext.hpp"
#include <string>
#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace AHCALRecoAlg {
    struct RmIsolatedHitAlgCfg{
        std::string in_recohit_key = "RecoHits";
        std::string out_recohit_key = "RmIsolatedHits";
        bool   useNmipWindow = true;
        double nmipMin = 0.2;
        double nmipMax = 3.0;
        double xy_cone = 42.0; // mm
        int maxMissingLayers = 0;
    };

    class RmIsolatedHitAlg final : public IAlg { // final to prevent inheritance
    public:
        RmIsolatedHitAlg(RunContext& ctx, std::string name)
            : IAlg(ctx, std::move(name)) {
                m_skiplayers = ctx.conditions.skipLayers;
                computeNeighborLayers();
            }
        void parse_cfg(const YAML::Node& n) override;
        void execute(EventStore& evt) override; // Main execution method
        void computeIsolatedFlags(const std::vector<AHCALRecoHit>& hits,
                                  std::vector<int>& isoFlag) const;
        bool isSkippedLayer(int layer) const {
            return std::find(m_skiplayers.begin(), m_skiplayers.end(), layer) != m_skiplayers.end();
        }
        void computeNeighborLayers();
    private:
        RmIsolatedHitAlgCfg m_cfg;
        std::vector<int> m_skiplayers;
        std::array<std::pair<int,int>, AHCALGeometry::Layer_No> before_layers; // layer -> (prev, next physical layer), -1 if none
        std::array<std::pair<int,int>, AHCALGeometry::Layer_No> after_layers;  // layer -> (prev, next physical layer), -1 if none
    };
} // namespace AHCALRecoAlg