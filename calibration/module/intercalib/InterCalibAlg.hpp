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
        double hg_bin_width = 25.0;
        int    min_hg_bins_for_fit = 8;
        double outlier_sigma_threshold = 3.0;

        // Example cells for visualization (using cellid directly)
        std::vector<int> example_cellids = {
            // Layer 13, Chip 3, Channel 13
            AHCALGeometry::CellID(13, 3, 13)
        };

        // Save mode
        bool save_all_channels_root = false;  // if true, save ROOT histograms for all channels; if false, save only example channels
        bool save_residual_histograms = false;  // if true, save LG-residual histograms for all fit channels

        // HitTag filtering
        bool require_hittag = true;
        bool require_yperx_over1 = true;
        bool output_bad_cells_json = true;
        std::string bad_cells_json_filename = "bad_cells.json";

        bool use_specific_fit_range_toL39C8 = true; // if true, use specific fit range for L39C8 (based on observed saturation behavior), otherwise use global fit range for all channels
        double hg_fit_max_L39C8 = 1500.0; // specific fit max for L39C8

        bool select_muon_hits = false; // if true, only fill histograms with events that have at least one hit with HG>1000 and LG>100 (muon-like hits)
        std::string in_track_key = "FittedTrack";
        std::string string_track_struct = "SimpleFittedTrack"; // or "Track"
        std::string track_selection_string = ""; // optional selection string for tracks, e.g. "nHits >= 5 && chi2/ndf < 5"
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
