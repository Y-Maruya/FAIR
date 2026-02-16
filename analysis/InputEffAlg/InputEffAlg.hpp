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
    struct InputEffAlgCfg{
        std::string in_data_key = "TLUData";
        std::string in_reco_hit_key = "RecoHits";

        // std::string string_track_struct = "SimpleFittedTrack"; // or "Track"
        // std::string track_selection_string = "";
        std::string out_filename = "Input_eff.root";
        bool write_to_png = false;
        bool write_to_txt = false;
        std::string out_txt_name = "Input_eff.txt";
        std::string out_png_dir = "Input_eff_png";
        std::vector<int> trigger_layer = {9,19,29,38};
        int ntrigger_layer = 4;
    };
    class InputEffAlg final : public IAlg{
    public:
        InputEffAlg(RunContext& rc, std::string name)
        : IAlg(rc, std::move(name)){}
        ~InputEffAlg() override;
        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& cfg) override;
        void initialize() override;

    private:
        InputEffAlgCfg cfg_;
        struct Impl;
        struct ImplDeleter {
            void operator()(Impl* p) const;
        };
        std::unique_ptr<std::ofstream> out_file_;
        std::unique_ptr<Impl, ImplDeleter> impl_;
        bool is_hittag1_exist(std::vector<AHCALRecoHit>& recoHits, int layer);
    };    
}