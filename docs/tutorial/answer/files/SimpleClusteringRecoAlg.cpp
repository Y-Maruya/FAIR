#include "SimpleClusteringRecoAlg.hpp"
#include "common/AlgRegistry.hpp"
#include <cmath>

AHCAL_REGISTER_ALG(AHCALRecoAlg::SimpleClusteringRecoAlg, "SimpleClusteringRecoAlg")

namespace AHCALRecoAlg {

void TutorialCluster::updateSummary() {
    nHitInCluster = static_cast<int>(hits.size());
    totalEdep = 0.0;
    for (const auto& h : hits) totalEdep += h.Edep;
    avgEdep = nHitInCluster > 0 ? totalEdep / nHitInCluster : 0.0;
}

static double dist3D(const AHCALRecoHit& a, const AHCALRecoHit& b) {
    const double dx = a.Xpos() - b.Xpos();
    const double dy = a.Ypos() - b.Ypos();
    const double dz = a.Zpos() - b.Zpos();
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void SimpleClusteringRecoAlg::initialize() {}
void SimpleClusteringRecoAlg::finalize() {}

void SimpleClusteringRecoAlg::execute(EventStore& evt) {
    auto& hits = evt.get<std::vector<AHCALRecoHit>>(cfg_.in_recohit_key);
    std::vector<TutorialCluster> clusters;

    for (const auto& h : hits) {
        bool attached = false;
        for (auto& c : clusters) {
            // Attach if close to ANY hit already in the cluster (single-linkage style)
            for (const auto& existing : c.hits) {
                if (dist3D(h, existing) <= cfg_.cluster_distance_mm) {
                    c.hits.push_back(h);
                    c.index.push_back(h.index);
                    attached = true;
                    break;
                }
            }
            if (attached) break;
        }
        if (!attached) {
            TutorialCluster c;
            c.hits.push_back(h);
            c.index.push_back(h.index);
            clusters.push_back(std::move(c));
        }
    }

    for (auto& c : clusters) c.updateSummary();
    evt.put(cfg_.out_cluster_key, std::move(clusters));
}

void SimpleClusteringRecoAlg::parse_cfg(const YAML::Node& n) {
    cfg_.in_recohit_key = n["in_recohit_key"].as<std::string>(cfg_.in_recohit_key);
    cfg_.out_cluster_key = n["out_cluster_key"].as<std::string>(cfg_.out_cluster_key);
    cfg_.cluster_distance_mm = n["cluster_distance_mm"].as<double>(cfg_.cluster_distance_mm);
}

} // namespace AHCALRecoAlg
