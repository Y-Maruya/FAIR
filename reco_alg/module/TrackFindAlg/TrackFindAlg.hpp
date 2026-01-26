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
        // Methods for track finding
        void SelectNonIsolatedHits(std::vector<AHCALRecoHit>& recohits, int config_index = 0);
        void SearchConnectedHits(std::vector<AHCALRecoHit>& recohits, int config_index = 0);
        void SetSkipLayers(const std::vector<int>& layers) {
            m_skiplayers = layers;
        }

    private:
        std::string m_in_recohit_key;
        std::string m_out_track_key;
        std::vector<AHCALRecoHit> m_recohits;
        std::vector<AHCALRecoHit> m_nonisolated_hits;
        std::vector<int> m_nonisolated_hit_indices;
        std::vector<std::vector<AHCALRecoHit>> m_connected_hits;
        std::vector<std::vector<int>> m_connected_hit_indices;
        std::vector<int> m_skiplayers = {};
        int m_number_of_configurations = 1;
        std::vector<int> m_isolation_distance = {42};//mm per layer (20mm)
        std::vector<int> m_min_hits_in_track = {6};
        std::vector<int> m_max_missing_layers = {2};
        std::vector<int> m_max_distance_xy = {42};//mm
        std::vector<int> m_num_later_layers_to_fit = {40};
        bool find = false;
    };

} // namespace AHCALRecoAlg