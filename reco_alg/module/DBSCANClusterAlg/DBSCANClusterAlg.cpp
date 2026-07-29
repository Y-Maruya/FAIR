#include "DBSCANClusterAlg.hpp"
#include "common/AlgRegistry.hpp"
#include <cmath>

AHCAL_REGISTER_ALG(AHCALRecoAlg::DBSCANClusterAlg, "DBSCANClusterAlg")

namespace AHCALRecoAlg {


static double dist2D(const AHCALRecoHit& a, const AHCALRecoHit& b) {
    const double dx = a.Xpos() - b.Xpos();
    const double dy = a.Ypos() - b.Ypos();
    return std::sqrt(dx * dx + dy * dy);
}

void DBSCANClusterAlg::initialize() {}
void DBSCANClusterAlg::finalize() {}

void DBSCANClusterAlg::execute(EventStore& evt) {
    auto& hits = evt.get<std::vector<AHCALRecoHit>>(cfg_.in_recohit_key);
    std::vector<Cluster> clusters;

    const int nHits = static_cast<int>(hits.size());
    if (nHits == 0) {
        evt.put(cfg_.out_cluster_key, std::move(clusters));
        return;
    }

    const int minPts = cfg_.minPts;

    // -1: unassigned/noise
    // >=0: cluster id
    std::vector<int> cluster_id(nHits, -1);
    std::vector<bool> visited(nHits, false);

    auto isKnownSkip = [&](const AHCALRecoHit& h) -> bool {
        const int layer = h.layer();
        const int chip  = h.chip();
        const int cellid = h.cellID;

        for (const auto& l : cfg_.known_skip_layer) {
            if (layer == l) return true;
        }

        for (const auto& lc : cfg_.known_skip_layer_chip) {
            if (layer == lc.first && chip == lc.second) return true;
        }
        for (const auto& c : cfg_.known_skip_cellid) {
            if (cellid == c) return true;
        }

        return false;
    };
    auto isSkippedPosition =
        [&](int layer, double x, double y) -> bool {

        // Entire skipped layer
        for (const int skippedLayer : cfg_.known_skip_layer) {
            if (layer == skippedLayer) {
                return true;
            }
        }

        int chip = -1;
        int channel = -1;
        AHCALGeometry::inverse(x, y, chip, channel);

        // Skipped chip
        for (const auto& [skippedLayer, skippedChip]
            : cfg_.known_skip_layer_chip) {

            if (layer == skippedLayer &&
                chip == skippedChip) {
                return true;
            }
        }

        // Skipped cell
        const int cellid =
            AHCALGeometry::CellID(layer, chip, channel);

        for (const int skippedCellID : cfg_.known_skip_cellid) {
            if (cellid == skippedCellID) {
                return true;
            }
        }

        return false;
    };
    auto layerDistanceAllowingSkips =
        [&](const AHCALRecoHit& a,
            const AHCALRecoHit& b) -> int {

        const int la = a.layer();
        const int lb = b.layer();

        if (la == lb) {
            return 0;
        }

        const double xa = a.Xpos();
        const double ya = a.Ypos();
        const double xb = b.Xpos();
        const double yb = b.Ypos();

        const int step = (lb > la) ? 1 : -1;

        int effective_gap = 0;

        for (int layer = la + step;
            layer != lb + step;
            layer += step) {

            const double fraction =
                static_cast<double>(layer - la) /
                static_cast<double>(lb - la);

            const double x =
                xa + fraction * (xb - xa);

            const double y =
                ya + fraction * (yb - ya);

            if (!isSkippedPosition(layer, x, y)) {
                ++effective_gap;
            }
        }

        return effective_gap;
    };

    auto regionQuery = [&](int i) -> std::vector<int> {
        std::vector<int> neighbors;

        if (isKnownSkip(hits[i])) return neighbors;

        for (int j = 0; j < nHits; ++j) {
            if (i == j) continue;
            if (isKnownSkip(hits[j])) continue;

            const double d_xy = dist2D(hits[i], hits[j]);
            const int d_layer = layerDistanceAllowingSkips(hits[i], hits[j]);

            if (d_xy <= cfg_.epsilon_xy_mm &&
                d_layer <= cfg_.epsilon_z_layer) {
                neighbors.push_back(j);
            }
        }

        return neighbors;
    };


    int current_cluster_id = 0;

    for (int i = 0; i < nHits; ++i) {
        if (visited[i]) continue;
        if (isKnownSkip(hits[i])) continue;

        visited[i] = true;

        auto neighbors = regionQuery(i);

        // Not a core point
        if (static_cast<int>(neighbors.size()) < minPts) {
            continue;
        }

        Cluster cluster;

        cluster_id[i] = current_cluster_id;
        cluster.hits.push_back(hits[i]);
        cluster.index.push_back(hits[i].index);
        cluster.cluster_id = current_cluster_id;

        std::vector<int> seeds = neighbors;

        for (size_t seed_pos = 0; seed_pos < seeds.size(); ++seed_pos) {
            const int j = seeds[seed_pos];

            if (isKnownSkip(hits[j])) continue;

            if (!visited[j]) {
                visited[j] = true;

                auto neighbors_j = regionQuery(j);

                if (static_cast<int>(neighbors_j.size()) >= minPts) {
                    for (const int k : neighbors_j) {
                        bool already_in_seeds = false;
                        for (const int s : seeds) {
                            if (s == k) {
                                already_in_seeds = true;
                                break;
                            }
                        }

                        if (!already_in_seeds) {
                            seeds.push_back(k);
                        }
                    }
                }
            }

            if (cluster_id[j] == -1) {
                cluster_id[j] = current_cluster_id;
                cluster.hits.push_back(hits[j]);
                cluster.index.push_back(hits[j].index);
                cluster.cluster_id = current_cluster_id;
            }
        }

        cluster.updateSummary();
        clusters.push_back(std::move(cluster));

        ++current_cluster_id;
    }

    evt.put(cfg_.out_cluster_key, std::move(clusters));
}

void DBSCANClusterAlg::parse_cfg(const YAML::Node& n) {
    cfg_.in_recohit_key = n["in_recohit_key"].as<std::string>(cfg_.in_recohit_key);
    cfg_.out_cluster_key = n["out_cluster_key"].as<std::string>(cfg_.out_cluster_key);
    cfg_.epsilon_xy_mm = n["epsilon_xy_mm"].as<double>(cfg_.epsilon_xy_mm);
    cfg_.epsilon_z_layer = n["epsilon_z_layer"].as<int>(cfg_.epsilon_z_layer);
    cfg_.minPts = n["minPts"].as<int>(cfg_.minPts);
    cfg_.known_skip_layer = n["known_skip_layer"].as<std::vector<int>>(cfg_.known_skip_layer);
    cfg_.known_skip_layer_chip = n["known_skip_layer_chip"].as<std::vector<std::pair<int,int>>>(cfg_.known_skip_layer_chip);
    cfg_.known_skip_cellid = n["known_skip_cellid"].as<std::vector<int>>(cfg_.known_skip_cellid);
}

} // namespace AHCALRecoAlg
