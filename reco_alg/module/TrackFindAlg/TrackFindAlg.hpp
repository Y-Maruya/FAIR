#pragma once

#include "common/EventStore.hpp"
#include "common/edm/EDM.hpp"
#include <string>
#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace AHCALRecoAlg {
    class TrackFindAlg final : public IAlg { // final to prevent inheritance
    public:
        TrackFindAlg(std::string in_recohit_key,
            std::string out_track_key)
            : m_in_recohit_key(std::move(in_recohit_key)),
            m_out_track_key(std::move(out_track_key)) {}

        void execute(EventStore& evt) override; // Main execution method

        // Optional knobs (macro 相当)
        void SetMipThresholdFrac(double frac) { m_thr_mip_frac = frac; } // default 0.5
        void SetUseTightOnly(bool v) { m_use_tight_only = v; }

        void SetSkipLayers(const std::vector<int>& layers) {
            m_skiplayers = layers;
        }

    private:
        struct TrackSummary {
            int label = 0;
            int nHits = 0;
            int nLayers = 0;
            double tanTheta = 999.0;
            double chi2x = 999.0;
            double chi2y = 999.0;
            double meanQ = 999.0;
            double rmsQ  = 999.0;
            bool isLoose = false;
            bool isTight = false;
        };

        // muon-finding core (macro の「muon を探すところ」)
        void computeIsolatedFlags(const std::vector<AHCALRecoHit>& hits,
                                  std::vector<int>& isoFlag) const;
        void buildTracks(const std::vector<AHCALRecoHit>& hits,
                         const std::vector<int>& isoFlag,
                         std::vector<int>& trackLabel,
                         std::vector<int>& trackLength,
                         int& numTrack) const;
        void summarizeTracks(const std::vector<AHCALRecoHit>& hits,
                             const std::vector<int>& trackLabel,
                             int numTrack,
                             std::vector<TrackSummary>& out) const;

        std::string m_in_recohit_key;
        std::string m_out_track_key;

        std::vector<int> m_skiplayers = {};

        // macro-like parameters
        double m_thr_mip_frac = 0.5;   // Nmip threshold
        bool   m_use_tight_only = false;

        // (left as-is; not used by the muon finder below)
        int m_number_of_configurations = 1;
        std::vector<int> m_isolation_distance = {42};//mm per layer (20mm)
        std::vector<int> m_min_hits_in_track = {6};
        std::vector<int> m_max_missing_layers = {2};
        std::vector<int> m_max_distance_xy = {42};//mm
        std::vector<int> m_num_later_layers_to_fit = {40};

        bool find = false;
    };

} // namespace AHCALRecoAlg