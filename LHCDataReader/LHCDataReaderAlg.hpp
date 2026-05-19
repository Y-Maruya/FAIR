#pragma once
#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include "common/edm/EDM.hpp"
#include <yaml-cpp/yaml.h>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>


namespace AHCALRecoAlg {
    struct LHCDataReaderAlgCfg {
        std::string in_rawdata_key = "TLURawData";
        std::string out_lhcdata_key = "LHCData";
        bool verbose = false;
    };

    class LHCDataReaderAlg final : public IAlg {
    public:
        LHCDataReaderAlg(RunContext& ctx, std::string name)
            : IAlg(ctx, std::move(name)) {}
        
        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& n) override;
        void initialize() override;
        void finalize() override;
    private:
        LHCDataReaderAlgCfg m_cfg;
        cool::IDatabasePtr m_db;
        cool::IFolderPtr m_folder;
        std::string m_normalizedDbId;
    };
} // namespace AHCALRecoAlg