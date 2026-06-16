#pragma once
#include "CoolKernel/IDatabase.h"
#include "CoolKernel/IFolder.h"

#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include "common/edm/EDM.hpp"
#include "common/edm/LumiRecord.hpp"
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
    struct LumiRecordReaderAlgCfg {
        std::string in_rawdata_key = "TLURawData";
        std::string out_lumirecord_key = "LumiRecord";
        std::string lumi_tag = "OflLumiAcct-Run3-007";
        bool verbose = false;
    };

    class LumiRecordReaderAlg final : public IAlg {
    public:
        LumiRecordReaderAlg(RunContext& ctx, std::string name)
            : IAlg(ctx, std::move(name)) {}
        
        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& n) override;
        void initialize() override;
        void finalize() override;
    private:
        LumiRecordReaderAlgCfg m_cfg;
        cool::IDatabasePtr m_db;
        cool::IFolderPtr m_folder;
        std::string m_normalizedDbId;
        
        // Caching for performance: avoid redundant COOL queries
        // Store last valid IOV range and corresponding data
        cool::ValidityKey m_cached_since = 0;
        cool::ValidityKey m_cached_until = 0;
        LumiRecord m_cached_data;
        bool m_cache_valid = false;
    };
} // namespace AHCALRecoAlg
