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
        bool pedestal_to_DB = false; // (未実装)
        std::string out_pedestal_filename = "pedestal.root";

        // Histogram / fit config (macro defaults)
        int    nbin = 800;
        double xmin = 0.0;
        double xmax = 2000.0;

        int    min_entries = 200;

        double nsigma_win1 = 2.0;
        double nsigma_win2 = 1.5;

        double sigma_min = 0.5;
        double sigma_max = 200.0;

        bool use_hittag = true;
        int  select_hittag = 0;
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