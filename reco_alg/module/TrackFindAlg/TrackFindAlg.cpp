#include "TrackFindAlg.hpp"


namespace AHCALRecoAlg {

    void TrackFindAlg::SelectNonIsolatedHits( int config_index) {
        m_nonisolated_hits.clear();
        for (const auto& hit : m_recohits) {
            bool is_isolated = true;
            for (const auto& other_hit : m_recohits) {
                if (hit.layer == other_hit.layer) continue;
                double dx = hit.x - other_hit.x;
                double dy = hit.y - other_hit.y;
                if (std::abs(dx) <= m_isolation_distance[config_index] &&
                    std::abs(dy) <= m_isolation_distance[config_index]) {
                    is_isolated = false;
                    break;
                }
            }
            if (!is_isolated) {
                m_nonisolated_hits.push_back(hit);
                m_nonisolated_hit_indices.push_back(&hit - &m_recohits[0]);
            }
        }
    }

    void TrackFindAlg::SearchConnectedHits(int config_index) {
        m_connected_hits.clear();
        std::vector<bool> used(m_nonisolated_hits.size(), false);

        // sort by layer
        std::sort(m_nonisolated_hit_indices.begin(), m_nonisolated_hit_indices.end(),
            [this](int a, int b) {
                return m_recohits[a].layer < m_recohits[b].layer;
            });
        // backward search
        for (size_t i = 0; i < m_nonisolated_hits.size(); ++i) {
            AHCALRecoHit current_hit = m_recohits[m_nonisolated_hit_indices[i]];
            int layer = current_hit.layer;
            if (std::find(m_skiplayers.begin(), m_skiplayers.end(), layer) != m_skiplayers.end()) {
                continue; // skip this layer
            }
            if (layer < AHCALGeometry::Layer_No - m_num_later_layers_to_fit[config_index]) {
                continue; // skip layers to search the track which escape from the hadron shower.
            }
            if (used[i]) continue;
            std::vector<AHCALRecoHit> track_hits;
            track_hits.push_back(current_hit);
            used[i] = true;
            int missing_layers = 0;

            for (int next_layer = layer - 1; next_layer >= 0; --next_layer) {
                if (std::find(m_skiplayers.begin(), m_skiplayers.end(), next_layer) != m_skiplayers.end()) {
                    continue; // skip this layer
                }
                bool found_hit_in_layer = false;
                for (size_t j = 0; j < m_nonisolated_hits.size(); ++j) {
                    if (used[j]) continue;
                    AHCALRecoHit next_hit = m_nonisolated_hits[j];
                    if (next_hit.layer != next_layer) continue;
                    double dx = next_hit.x - current_hit.x;
                    double dy = next_hit.y - current_hit.y;
                    if (std::abs(dx) <= m_max_distance_xy[config_index] &&
                        std::abs(dy) <= m_max_distance_xy[config_index]) {
                        track_hits.push_back(next_hit);
                        used[j] = true;
                        current_hit = next_hit;
                        found_hit_in_layer = true;
                        break;
                    }
                }
                if (!found_hit_in_layer) {
                    missing_layers++;
                    if (missing_layers > m_max_missing_layers[config_index]) {
                        break;
                    }
                }
            }

    }

