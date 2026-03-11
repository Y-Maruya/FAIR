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
#include <fstream>

namespace AHCALRecoAlg{
    struct MuonEffAlgCfg{
        std::string in_hit_key = "RawHits";
        std::string in_rawdata_key = "TLUData";
        std::string in_track_key = "Tracks";

        std::string string_track_struct = "SimpleFittedTrack"; // or "Track"
        std::string track_selection_string = "";
        std::string out_filename = "muon_eff.root";
        bool write_to_png = false;
        std::string out_png_dir = "muon_eff_png";
        bool write_to_txt = false;
        std::string out_txt_filename = "muon_eff.txt";
        bool hittag_or_MIP_cut = false; // if true, use hit tag, otherwise use MIP cut
        double MIP_cut = 0.5;
        double xy_size_threshold = 17.0; 
        int ntrigger_layer = 4;
        std::vector<int> trigger_layer = {9, 19, 29, 38}; 
    };
    class MuonEffAlg final : public IAlg{
    public:
        MuonEffAlg(RunContext& rc, std::string name)
        : IAlg(rc, std::move(name)){}
        ~MuonEffAlg() override;
        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& cfg) override;
        void initialize() override;
    private:
        MuonEffAlgCfg cfg_;
        bool is_hittag1_exist(std::vector<AHCALRawHit>& rawHits, int layer);
        struct Impl;
        struct ImplDeleter {
            void operator()(Impl* p) const;
        };
        std::unique_ptr<std::ofstream> out_file_;
        std::unique_ptr<Impl, ImplDeleter> impl_;
    };    
}