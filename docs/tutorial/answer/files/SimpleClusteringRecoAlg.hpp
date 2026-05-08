#pragma once
#include <yaml-cpp/yaml.h>
#include "common/EventStore.hpp"
#include "common/edm/EDM.hpp"
#include "common/IAlg.hpp"
#include "common/RunContext.hpp"
#include <string>
#include <vector>

namespace AHCALRecoAlg {

struct TutorialCluster {
    std::vector<AHCALRecoHit> hits;
    std::vector<int> index;
    int nHitInCluster = 0;
    double totalEdep = 0.0;
    double avgEdep = 0.0;
    void updateSummary();
};

struct SimpleClusteringRecoAlgCfg {
    std::string in_recohit_key = "RecoHits";
    std::string out_cluster_key = "Clusters";
    double cluster_distance_mm = 45.0;
};

class SimpleClusteringRecoAlg final : public IAlg {
public:
    SimpleClusteringRecoAlg(RunContext& rc, std::string name)
      : IAlg(rc, std::move(name)) {}
    void initialize() override;
    void execute(EventStore& evt) override;
    void finalize() override;
    void parse_cfg(const YAML::Node& n) override;
private:
    SimpleClusteringRecoAlgCfg cfg_;
};

} // namespace AHCALRecoAlg
