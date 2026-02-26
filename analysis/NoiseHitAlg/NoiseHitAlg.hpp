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
    struct NoiseHitAlgCfg{
        std::string in_raw_hit_key = "RawHits";
        std::string in_track_hit_key = "FittedTracks";

        std::string string_track_struct = "SimpleFittedTrack"; // or "Track" //
        std::string track_selection_string = "";  // "" means not use outtrack hit.
        int nHits_threshold = 0; // if > 0, only consider events with more than this number of hits
        std::string out_filename = "NoiseHit_out.root";
        bool write_to_png = false;
        bool write_to_txt = false;
        std::string out_txt_name = "NoiseHit_out.txt";
        std::string out_png_dir = "NoiseHit_png";
        double xmin = 0;
        double xmax = 1000;
        int nbin = 100;
    };
    class NoiseHitAlg final : public IAlg{
    public:
        NoiseHitAlg(RunContext& rc, std::string name)
        : IAlg(rc, std::move(name)){}
        ~NoiseHitAlg() override;
        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& cfg) override;
        void initialize() override;

    private:
        NoiseHitAlgCfg cfg_;
        struct Impl;
        struct ImplDeleter {
            void operator()(Impl* p) const;
        };
        std::unique_ptr<std::ofstream> out_file_;
        std::unique_ptr<Impl, ImplDeleter> impl_;
    };    
}