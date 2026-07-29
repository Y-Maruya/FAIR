#pragma once

#include <yaml-cpp/yaml.h>

#include "common/EventStore.hpp"
#include "common/edm/EDM.hpp"
#include "common/edm/Cluster.hpp"
#include "common/IAlg.hpp"
#include "common/RunContext.hpp"

#include <string>
#include <vector>

namespace AHCALRecoAlg {

struct PCAClusterAlgCfg {
    std::string in_cluster_key = "Clusters";
    std::string out_cluster_key = "PCAClusters";

    int min_hits = 3;

    double track_linearity_threshold = 0.85;
    double track_width_threshold_mm = 60.0;

    bool energy_weighted = false;
};

class PCAClusterAlg final : public IAlg {
public:
    PCAClusterAlg(RunContext& rc, std::string name)
        : IAlg(rc, std::move(name)) {}

    void initialize() override;
    void execute(EventStore& evt) override;
    void finalize() override;
    void parse_cfg(const YAML::Node& n) override;

private:
    PCAClusterAlgCfg cfg_;

    void computePCA(Cluster& c) const;
};

} // namespace AHCALRecoAlg