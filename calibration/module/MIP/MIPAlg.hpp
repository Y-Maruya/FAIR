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
        bool mip_to_DB = false;  // Not implemented yet
        std::string out_mip_filename = "mip.root";

        int nbin = 4096;
        double xmin = 0.0;
        double xmax = 4096.0;

        int min_entries = 200;

    };

    class MIPAlg final : public IAlg{
    public:
        MIPAlg(RunContext& rc, std::string name)
        : IAlg(rc, std::move(name)){}
        ~MIPAlg() override;

        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& cfg) override;
    
    private:
        MIPAlgCfg cfg_;
        struct Impl;
        struct ImplDeleter {
            void operator()(Impl* p) const;
        };
        std::unique_ptr<Impl, ImplDeleter> impl_;
    };    
}