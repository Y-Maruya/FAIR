#pragma once

#include "calibration/module/MIP/MIPAlg.hpp"
#include "common/IAlg.hpp"

#include <memory>
#include <string>
#include <unordered_map>

class TH1D;

namespace AHCALRecoAlg {

// Runs the real-data MIP selection/ADC fit on digitized hits, then records
// visible SimHit energy for precisely the raw hits filled by that selection.
class MIPSimAlg final : public IAlg {
public:
    MIPSimAlg(RunContext& rc, std::string name);
    ~MIPSimAlg() override;

    void parse_cfg(const YAML::Node& cfg) override;
    void init_by_run() override;
    void execute(EventStore& evt) override;

private:
    void write_truth();

    MIPAlg mip_;
    std::string in_simhit_key_ = "SimHits";
    std::string in_rawhit_key_ = "SimRawHits";
    std::string selected_key_;
    std::string out_simhit_filename_ = "mip_sim_truth.root";
    int edep_nbin_ = 200;
    double edep_max_ = 2.0;  // visible MeV
    std::unordered_map<int, std::unique_ptr<TH1D>> edep_hist_;
    long long n_selected_ = 0;
    long long n_missing_simhit_ = 0;
};

} // namespace AHCALRecoAlg
