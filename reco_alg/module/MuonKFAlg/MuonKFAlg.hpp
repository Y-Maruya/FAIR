#pragma once
#include <bitset>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include "common/edm/EDM.hpp"



namespace AHCALRecoAlg {

    struct MuonKFAlgCfg {
        std::string in_recohit_key = "RecoHits";
        std::string out_track_key = "MuonKFTrack";
        // window
        int lastNLayers = 40;          // use only last N physical layers (based on layer index)
        int minUsedLayers = 10;        // require this many updated layers
        int maxConsecutiveSkips = 3;   // terminate if too many layers skipped in a row

        // hit selection
        bool   useNmipWindow = true;
        double nmipMin = 0.2;
        double nmipMax = 3.0;

        // measurement resolution (mm). If <=0, defaults to pitch/sqrt(12)
        double measSigmaXY_mm = 0.0;

        // process noise (very simple): sigma_theta [rad] added per propagation step
        double sigmaTheta = 0.004;     // 4 mrad start; tune

        // gating threshold on d2 = r^T S^-1 r
        double gateD2 = 9.0;           // ~3sigma in 2D

        // seeding
        int seedLayerGap = 4;          // choose two seed layers separated by >= gap
        int maxSeedHitsPerLayer = 8;   // keep top-K candidate hits per layer for seeds

        // skip layer mask (欠損レイヤー等)
        std::bitset<40> skipLayer;
    };

    // Run KF muon tagging on a single event
    bool find_muon_track_kf(const std::vector<AHCALRecoHit>& recoHits,
                            Track& bestOut,
                            const MuonKFAlgCfg& cfg);


    class MuonKFAlg final : public IAlg {
    public:
        MuonKFAlg(RunContext& ctx, std::string name) 
            : IAlg(ctx, std::move(name)) {}
        void parse_cfg(const YAML::Node& n) override;
        MuonKFAlgCfg& config() { return m_cfg; }
        const MuonKFAlgCfg& config() const { return m_cfg; }

        // event-by-event entry point
        void execute(EventStore& evt) override;

        // (optional) expose result (or you can evt.put it in cpp)
        const Track& last_track() const { return m_last; }
        bool last_found() const { return m_last_found; }

    private:
        // std::string m_in_recohit_key;
        // std::string m_out_track_key;
        MuonKFAlgCfg m_cfg;

        // cache last result
        Track m_last;
        bool m_last_found = false;
    };

} // namespace AHCALRecoAlg
