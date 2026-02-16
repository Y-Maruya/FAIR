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
    struct VetoEffAlgCfg{
        std::string in_data_key = "TLUData";
        std::string in_track_key = "Tracks";

        std::string string_track_struct = "SimpleFittedTrack"; // or "Track"
        std::string track_selection_string = "";
        std::string out_filename = "veto_eff.root";
        bool write_to_png = false;
        bool write_to_txt = false;
        std::string out_txt_name = "veto_eff.txt";
        std::string out_png_dir = "veto_eff_png";
        double xy_size_threshold = 300;
        double x_center = 0; //mm, for track extrapolation
        double y_center = 0; //mm, for track extrapolation
        double z_pos_input4 = -100; //mm, for track extrapolation
        double z_pos_input5 = -100; //mm, for track extrapolation
        double threshold_input4 = 75.0; // for plotting
        double threshold_input5 = 75.0;
    };
    class VetoEffAlg final : public IAlg{
    public:
        VetoEffAlg(RunContext& rc, std::string name)
        : IAlg(rc, std::move(name)){}
        ~VetoEffAlg() override;
        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& cfg) override;
        void initialize() override;
    private:
        VetoEffAlgCfg cfg_;
        struct Impl;
        struct ImplDeleter {
            void operator()(Impl* p) const;
        };
        std::unique_ptr<std::ofstream> out_file_;
        std::unique_ptr<Impl, ImplDeleter> impl_;
    };    
}