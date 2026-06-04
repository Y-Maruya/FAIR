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
    struct MIPAlgCfg{
        std::string in_rawhit_key = "RawHits";
        std::string in_track_key = "Tracks";

        std::string string_track_struct = "SimpleFittedTrack"; // or "Track"
        std::string track_selection_string = "";
        bool mip_to_file = true;
        // bool mip_to_DB = false;  // Not implemented yet
        std::string out_mip_filename = "mip.root";
        bool output_to_png = false;
        std::string out_png_dir = "png/";
        bool fit = true;
        int nbin = 2096;
        double xmin = -96.0;
        double xmax = 2000.0;

        int min_entries = 200;

        bool mip_to_json = false;
        std::string out_json_dirname = ".";
        bool calculate_fwhm = true;

        bool substrate_pedestal = true;
        bool read_pedestal_from_ROOT = false;
        std::string in_pedestal_file = "pedestal.root";
        bool read_pedestal_from_DB = true; 

        double xy_size_threshold = 15.0; // mm

    };

    class MIPAlg final : public IAlg{
    public:
        MIPAlg(RunContext& rc, std::string name)
        : IAlg(rc, std::move(name)){}
        ~MIPAlg() override;

        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& cfg) override;
        void init_by_run() override;
    private:
        MIPAlgCfg cfg_;

        void ensure_impl();
        
        struct Impl;
        struct ImplDeleter {
            void operator()(Impl* p) const;
        };
        std::unique_ptr<Impl, ImplDeleter> impl_;
    };    
}