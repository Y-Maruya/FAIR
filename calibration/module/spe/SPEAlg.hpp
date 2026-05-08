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

namespace AHCALRecoAlg {

    struct SPEAlgCfg {
        // Input keys
        std::string in_rawhit_key = "RawHits";
        std::string in_track_key = "Tracks";

        // Track handling
        std::string string_track_struct = "SimpleFittedTrack"; // or "Track"
        std::string track_selection_string = "";

        // Output modes
        bool mip_to_file = true;
        bool mip_to_json = false;
        bool output_to_png = false;

        // File paths
        std::string out_mip_filename = "mip_spe.root";
        std::string out_png_dir = "png/";
        std::string out_json_dirname = ".";

        // MIP histogram parameters
        int nbin = 512;
        double xmin = 0.0;
        double xmax = 1000.0;
        int min_entries = 200;

        // Fitting control
        bool fit = true;  // Perform MIP fitting

        // FFT parameters
        int npad_req = 8192;
        double gain_min = 5.0;
        double gain_max = 30.0;
        int min_ent_fft = 500;
        double edge_gain_tol = 0.01;
        int kedge_skip = 4;

        // Output control
        bool save_per_channel_hists = true;
        int save_max_channel_hists = 600;
        bool do_event_scan = false;
        bool save_examples = true;

        // Pedestal handling
        bool substrate_pedestal = true;
        bool read_pedestal_from_ROOT = false;
        std::string in_pedestal_file = "pedestal.root";
        bool read_pedestal_from_DB = true;

        // Processing mode
        int processing_mode = 2; // 0=raw, 1=muon, 2=both

        // SNR threshold for efficiency calculation
        double snr_threshold = 3.0;

        // Calculate FWHM during fit
        bool calculate_fwhm = true;
    };

    class SPEAlg final : public IAlg {
    public:
        SPEAlg(RunContext& rc, std::string name)
            : IAlg(rc, std::move(name)) {}
        ~SPEAlg() override;

        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& cfg) override;
        void init_by_run() override;

    private:
        SPEAlgCfg cfg_;

        void ensure_impl();

        struct Impl;
        struct ImplDeleter {
            void operator()(Impl* p) const;
        };
        std::unique_ptr<Impl, ImplDeleter> impl_;
    };

} // namespace AHCALRecoAlg
