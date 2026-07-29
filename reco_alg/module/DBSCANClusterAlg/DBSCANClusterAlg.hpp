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


struct DBSCANClusterAlgCfg {
    std::string in_recohit_key = "RecoHits";
    std::string out_cluster_key = "Clusters";
    double epsilon_xy_mm = 45.0;
    int epsilon_z_layer = 1;
    int minPts = 1;
    std::vector<int> known_skip_layer = {0,1,2,3,4,9,14,28};
    std::vector<std::pair<int,int>> known_skip_layer_chip = {{21,1},{24,3},{24,4},{24,5},{24,6},{24,7}};
    std::vector<int> known_skip_cellid = {140013,870013,1100034,2770003,3200003};
};
class DBSCANClusterAlg final : public IAlg {
public:
    DBSCANClusterAlg(RunContext& rc, std::string name)
      : IAlg(rc, std::move(name)) {}
    void initialize() override;
    void execute(EventStore& evt) override;
    void finalize() override;
    void parse_cfg(const YAML::Node& n) override;
private:
    DBSCANClusterAlgCfg cfg_;
};

} // namespace AHCALRecoAlg
