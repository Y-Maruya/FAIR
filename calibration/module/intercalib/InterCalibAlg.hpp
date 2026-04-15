#pragma once
#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include "common/edm/EDM.hpp"

#include <memory>
#include <string>
#include <vector>

namespace AHCALRecoAlg {
    struct InterCalibAlgCfg {
        std::string in_rawhit_key = "RawHits";
        bool read_pedestal_from_ROOT = false;
        std::string in_pedestal_file = "pedestal.root";
        bool read_pedestal_from_DB = true;
        bool intercalib_to_file = true;
        bool intercalib_to_json = true;
        std::string out_intercalib_filename = "intercalib.root";
        std::string out_json_dirname = ".";

        // Fit parameters
        double hg_fit_max = 1800.0;
        double lg_fit_min = -30.0;
        double hg_fit_min = -100.0;
        int    min_points = 200;

        // Example channel for visualization
        std::vector<int> example_layers = {13};
        std::vector<int> example_chips = {3};
        std::vector<int> example_chs = {13};

        // HitTag filtering
        bool require_hittag = true;
        bool require_yperx_over1 = true;
        bool output_bad_cells_json = true;
        std::string bad_cells_json_filename = "bad_cells.json";

        bool use_specific_fit_range_toL39C8 = true; // if true, use specific fit range for L39C8 (based on observed saturation behavior), otherwise use global fit range for all channels
        double hg_fit_max_L39C8 = 1500.0; // specific fit max for L39C8
    };

    class InterCalibAlg final : public IAlg {
    public:
        InterCalibAlg(RunContext& rc, std::string name)
        : IAlg(rc, std::move(name)) {}
        ~InterCalibAlg() override;
        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& cfg) override;
        void init_by_run() override;
    private:
        InterCalibAlgCfg cfg_;

        void ensure_impl();

        struct Impl;

        // Custom deleter avoids requiring complete type in headers/TUs.
        struct ImplDeleter {
            void operator()(Impl* p) const;
        };
        std::unique_ptr<Impl, ImplDeleter> impl_;
    };
}
