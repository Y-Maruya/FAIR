#pragma once
#include <yaml-cpp/yaml.h>
#include "common/EventStore.hpp"
#include "common/edm/EDM.hpp"
#include "common/IAlg.hpp"
#include "common/RunContext.hpp"
#include <string>
#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace AHCALRecoAlg{
    struct MuonEffAlgCfg{
        std::string in_hit_key = "RawHits";
        std::string in_track_key = "Tracks";

        std::string string_track_struct = "SimpleFittedTrack"; // or "Track"
        std::string track_selection_string = "";
        std::string out_filename = "muon_eff.root";
        bool write_to_png = false;
        std::string out_png_dir = "muon_eff_png";
        bool hittag_or_MIP_cut = false; // if true, use hit tag, otherwise use MIP cut
        double MIP_cut = 0.5;
        double xy_size_threshold = 17.0; 
    };
    class MuonEffAlg final : public IAlg{
    public:
        MuonEffAlg(RunContext& rc, std::string name)
        : IAlg(rc, std::move(name)){}
        ~MuonEffAlg() override;
        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& cfg) override;
    
    private:
        MuonEffAlgCfg cfg_;
        struct Impl;
        struct ImplDeleter {
            void operator()(Impl* p) const;
        };
        std::unique_ptr<Impl, ImplDeleter> impl_;
    };    
}