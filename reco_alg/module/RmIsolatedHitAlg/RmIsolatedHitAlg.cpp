#include "RmIsolatedHitAlg.hpp"

#include "common/edm/EDM.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"
AHCAL_REGISTER_ALG(AHCALRecoAlg::RmIsolatedHitAlg, "RmIsolatedHitAlg")
namespace AHCALRecoAlg {

void RmIsolatedHitAlg::computeNeighborLayers() {
    for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer) {
        // Find previous physical layer
        int prev = -1;
        for (int l = layer - 1; l >= 0; --l) {
            if (!isSkippedLayer(l)) {
                prev = l;
                break;
            }
        }
        // Find next physical layer
        int next = -1;
        for (int l = layer + 1; l < AHCALGeometry::Layer_No; ++l) {
            if (!isSkippedLayer(l)) {
                next = l;
                break;
            }
        }
        before_layers[layer] = std::make_pair(layer, prev);
        after_layers[layer] = std::make_pair(layer, next);
    }
}
void RmIsolatedHitAlg::computeIsolatedFlags(const std::vector<AHCALRecoHit>& hits,
                                              std::vector<int>& isoFlag) const {
    const double cone = m_cfg.xy_cone;
    const int maxMissingLayers = m_cfg.maxMissingLayers;
    const bool useNmipWindow = m_cfg.useNmipWindow;
    const double nmipMin = m_cfg.nmipMin;
    const double nmipMax = m_cfg.nmipMax;
    for (size_t i = 0; i < hits.size(); ++i) {
        const auto& hit_i = hits[i];
        if (useNmipWindow && (hit_i.Nmip < nmipMin || hit_i.Nmip > nmipMax)) {
            isoFlag[i] = 1; 
            continue;
        }
        const int li = hit_i.layer();
        if (isSkippedLayer(li)) {
            isoFlag[i] = 1; 
            continue;
        }
        bool hasNeighbor = false;
        const double xi = hit_i.Xpos();
        const double yi = hit_i.Ypos();
        int prevLayer = before_layers[li].second;
        int nextLayer = after_layers[li].second;
        for (size_t j = 0; j < hits.size(); ++j) {
            if (i == j) continue;
            const auto& hit_j = hits[j];
            if (useNmipWindow && (hit_j.Nmip < nmipMin || hit_j.Nmip > nmipMax)) continue;
            if (isSkippedLayer(hit_j.layer())) continue;
            bool in = false;
            if (hit_j.layer() == li) {
                in = (std::fabs(hit_j.Xpos() - xi) < cone && std::fabs(hit_j.Ypos() - yi) < cone);
            }
            for (int dl = 0; dl <= maxMissingLayers; ++dl){
                if (hit_j.layer() == prevLayer-dl) {
                    in = (std::fabs(hit_j.Xpos() - xi) < cone && std::fabs(hit_j.Ypos() - yi) < cone);
                    break;
                }
                if (hit_j.layer() == nextLayer+dl) {
                    in = (std::fabs(hit_j.Xpos() - xi) < cone && std::fabs(hit_j.Ypos() - yi) < cone);
                    break;
                }
            }
            if (in) {
                hasNeighbor = true; 
                break;
            }
        }
        isoFlag[i] = hasNeighbor ? 0 : 1;
    }
}
        


void RmIsolatedHitAlg::execute(EventStore& evt) {
    auto& hits = evt.get<std::vector<AHCALRecoHit>>(m_cfg.in_recohit_key);
    std::vector<AHCALRecoHit> output_hits;
    output_hits.reserve(hits.size());

    std::vector<int> isoFlag(hits.size(), 0);
    computeIsolatedFlags(hits, isoFlag);
    for (size_t i = 0; i < hits.size(); ++i) {
        if (isoFlag[i] == 0) { // Not isolated
            output_hits.push_back(hits[i]);
        }
    }
    evt.put(m_cfg.out_recohit_key, std::move(output_hits));
}



void RmIsolatedHitAlg::parse_cfg(const YAML::Node& n) {
    m_cfg.in_recohit_key = n["in_recohit_key"].as<std::string>(m_cfg.in_recohit_key);
    m_cfg.out_recohit_key = n["out_recohit_key"].as<std::string>(m_cfg.out_recohit_key);
    m_cfg.useNmipWindow = n["useNmipWindow"].as<bool>(m_cfg.useNmipWindow);
    m_cfg.nmipMin = n["nmipMin"].as<double>(m_cfg.nmipMin);
    m_cfg.nmipMax = n["nmipMax"].as<double>(m_cfg.nmipMax);
    m_cfg.xy_cone = n["xy_cone"].as<double>(m_cfg.xy_cone);
    m_cfg.maxMissingLayers = n["maxMissingLayers"].as<int>(m_cfg.maxMissingLayers);
}
}