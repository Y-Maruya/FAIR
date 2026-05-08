#pragma once
#include <yaml-cpp/yaml.h>
#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include "common/RunContext.hpp"
#include "SimpleClusteringRecoAlg.hpp"
#include <memory>
#include <string>
#include <vector>

namespace AHCALRecoAlg {

struct ClusterBasicAnaAlgCfg {
    std::string in_cluster_key = "Clusters";
    double edep_threshold = 5.0;
    double threshold_scan_max = 20.0;
    double threshold_scan_step = 1.0;
    std::string out_root_filename = "cluster_ana.root";
    bool write_to_png = true;
    std::string out_png_dir = "cluster_ana_png";

};

class ClusterBasicAnaAlg final : public IAlg {
public:
    ClusterBasicAnaAlg(RunContext& rc, std::string name)
      : IAlg(rc, std::move(name)) {}
    ~ClusterBasicAnaAlg() override;
    void initialize() override;
    void execute(EventStore& evt) override;
    void finalize() override;
    void parse_cfg(const YAML::Node& n) override;
private:
    struct Impl;
    struct ImplDeleter { void operator()(Impl* p) const; };
    ClusterBasicAnaAlgCfg cfg_;
    std::unique_ptr<Impl, ImplDeleter> impl_;
};

} // namespace AHCALRecoAlg
