#pragma once
#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include "common/edm/EDM.hpp"

#include <memory>
#include <string>

namespace AHCALRecoAlg{
    struct PedestalAlgCfg{
        std::string in_rawhit_key = "RawHits";

        bool pedestal_to_file = true;
        bool pedestal_to_DB = false;
        bool pedestal_to_json = true;
        std::string out_pedestal_filename = "pedestal.root";
        std::string out_json_dirname = ".";

        // Histogram / fit config (macro defaults)
        int    nbin = 800;
        double xmin = 0.0;
        double xmax = 800.0;

        int    min_entries = 200;

        double nsigma_win1 = 2.0;
        double nsigma_win2 = 1.5;

        double sigma_min = 0.5;
        double sigma_max = 200.0;

        bool use_hittag = true;
        int  select_hittag = 0;
        
        bool hg_cut_failed = true; // whether to exclude failed fits from output
        bool hg_cut_rms_per_sigma = false; // whether to exclude fits with RMS/sigma > threshold
        double hg_max_rms_per_sigma = 4.0; // threshold for RMS/sigma cut 
        bool hg_cut_max_sigma_and_max_rms = false; // whether to exclude fits with sigma > threshold AND rms > threshold
        double hg_max_sigma = 20.0; // threshold for max sigma cut
        double hg_max_rms = 35.0; // threshold for max rms cut
        bool hg_cut_sigma_vs_chi2 = false; // whether to exclude fits with sigma vs chi2 outside envelope
        double hg_sigma_intercept = 35.0; // y-axis intercept of sigma vs chi2 envelope line
        double hg_chi2_ndf_intercept = 80.0; // x-axis intercept of sigma vs chi2 envelope line

        bool lg_cut_failed = true; // whether to exclude failed fits from output
        bool lg_cut_max_sigma = false; // whether to exclude fits with sigma > threshold
        double lg_max_sigma = 4.0; // threshold for max sigma cut
    };

    class PedestalAlg final : public IAlg{
    public:
        PedestalAlg(RunContext& rc, std::string name)
        : IAlg(rc, std::move(name)){}
        ~PedestalAlg() override;

        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& cfg) override;

    private:
        PedestalAlgCfg cfg_;

        struct Impl;

        // Custom deleter avoids requiring complete type in headers/TUs.
        struct ImplDeleter {
            void operator()(Impl* p) const;
        };

        std::unique_ptr<Impl, ImplDeleter> impl_;
    };
}