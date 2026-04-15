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
    struct VetoAnaAlgCfg{
        std::string in_data_key = "TLUData";
        // std::string string_track_struct = "SimpleFittedTrack"; // or "Track"
        // std::string track_selection_string = "";
        std::string out_filename = "Input_eff.root";
        bool write_to_png = false;
        bool write_to_txt = false;
        std::string out_txt_name = "Input_eff.txt";
        std::string out_png_dir = "Input_eff_png";
        double assumed_deadtime_ms = 1.3; // in milliseconds, for rate calculation
        double assumed_auto_triggered_ms = 156000*25e-6;
        std::vector<int> trigger_layer = {9,19,29,38};
        int ntrigger_layer = 4;
        double time_window_s = 2; // in seconds, for rate calculation
    };
    class VetoAnaAlg final : public IAlg{
    public:
        VetoAnaAlg(RunContext& rc, std::string name)
        : IAlg(rc, std::move(name)){}
        ~VetoAnaAlg() override;
        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& cfg) override;
        void initialize() override;

    private:
        VetoAnaAlgCfg cfg_;
        struct Impl;
        struct ImplDeleter {
            void operator()(Impl* p) const;
        };
        std::unique_ptr<std::ofstream> out_file_;
        std::unique_ptr<Impl, ImplDeleter> impl_;
        int last_cycle_id = -1;
        int last_TLU_timestamp = -1;
        int last_trigger_id = -1;
    };    
}