#include "CoolApplication/DatabaseSvcFactory.h"
#include "CoolKernel/ChannelSelection.h"
#include "CoolKernel/Exception.h"
#include "CoolKernel/IDatabase.h"
#include "CoolKernel/IDatabaseSvc.h"
#include "CoolKernel/IFolder.h"
#include "CoolKernel/IObject.h"
#include "CoolKernel/IObjectIterator.h"
#include "CoolKernel/IRecord.h"
#include "CoolKernel/ValidityKey.h"

#include "LHCDataReaderAlg.hpp"
#include "common/AlgRegistry.hpp"
#include "common/Logger.hpp"
#include "common/edm/LHCData.hpp"
#include "common/edm/RawData.hpp"
#include "common/config/YAMLUtil.hpp"


namespace AHCALRecoAlg {

constexpr const char* kDefaultDb =
  "sqlite_file:///cvmfs/faser.cern.ch/repo/sw/database/DBRelease/current/"
  "sqlite200/ALLP200.db/CONDBR3";
constexpr const char* kFolder = "/LHC/BeamData";

struct Options {
  unsigned int run = 0;
  cool::ValidityKey timestampNs = 0;
  std::string db = kDefaultDb;
  bool hasRun = false;
  bool hasTimestamp = false;
  bool verbose = false;
};

template <typename T>
T parseUnsigned(const std::string& text, const std::string& optionName) {
  if (text.empty() || text[0] == '-') {
    LOG_ERROR("Invalid value for {}: {}", optionName, text);
    throw std::runtime_error("Invalid value for " + optionName + ": " + text);
  }
  std::size_t parsed = 0;
  unsigned long long value = 0;
  try {
    value = std::stoull(text, &parsed, 0);
  } catch (const std::exception&) {
    throw std::runtime_error("Invalid value for " + optionName + ": " + text);
  }

  if (parsed != text.size() ||
      value > static_cast<unsigned long long>(std::numeric_limits<T>::max())) {
    throw std::runtime_error("Invalid value for " + optionName + ": " + text);
  }

  return static_cast<T>(value);
}

std::string normalizeDbId(const std::string& db) {
  const std::string prefix = "sqlite_file://";
  if (db.compare(0, prefix.size(), prefix) != 0) {
    return db;
  }

  const std::size_t dbMarker = db.rfind(".db/");
  if (dbMarker == std::string::npos) {
    throw std::runtime_error(
      "sqlite_file DB must have the form sqlite_file:///path/file.db/DBNAME: " +
      db);
  }

  const std::string schema =
    db.substr(prefix.size(), dbMarker + 3 - prefix.size());
  const std::string dbName = db.substr(dbMarker + 4);
  if (schema.empty() || dbName.empty()) {
    throw std::runtime_error(
      "sqlite_file DB must have a schema path and dbname: " + db);
  }

  return "sqlite://;schema=" + schema + ";dbname=" + dbName;
}

LHCData readBeamData(const cool::IObject& object) {
  const cool::IRecord& payload = object.payload();
  LHCData data;
  data.since = object.since();
  data.until = object.until();
  data.fillNumber = payload["FillNumber"].data<int>();
  data.stableBeams = payload["StableBeams"].data<bool>();
  data.beamMode = payload["BeamMode"].data<std::string>();
  data.machineMode = payload["MachineMode"].data<std::string>();
  data.beamType1 = payload["BeamType1"].data<int>();
  data.beamType2 = payload["BeamType2"].data<int>();
  data.beamEnergyGeV = payload["BeamEnergyGeV"].data<int>();
  data.betaStar = payload["BetaStar"].data<float>();
  data.crossingAngle = payload["CrossingAngle"].data<float>();
  data.injectionScheme = payload["InjectionScheme"].data<std::string>();
  return data;
}
Options parseArgsFromEvent(const AHCALTLURawData& rawdata, bool verbose) {
    Options options;
    options.run = rawdata.run_number;
    options.timestampNs = static_cast<cool::ValidityKey>(rawdata.Timestamp_cpu_mus * 1000); // Convert microseconds to nanoseconds
    options.hasRun = true;
    options.hasTimestamp = true;
    options.verbose = verbose;
    return options;
}

void LHCDataReaderAlg::execute(EventStore& evt) {
    const auto& rawdata = evt.get<AHCALTLURawData>(m_cfg.in_rawdata_key);
    LHCData lhcdata;
    try {
        const Options options = parseArgsFromEvent(rawdata, m_cfg.verbose);
        
        // Check cache: if timestamp falls in cached IOV range, reuse the result
        if (m_cache_valid && 
            options.timestampNs >= m_cached_since && 
            options.timestampNs < m_cached_until) {
            lhcdata = m_cached_data;
            if (m_cfg.verbose) {
                LOG_DEBUG("Using cached LHC data for timestamp {}", options.timestampNs);
            }
        } else {
            // Cache miss: query COOL database
            cool::IObjectIteratorPtr it = 
                m_folder->browseObjects(options.timestampNs, options.timestampNs, cool::ChannelSelection::all());
            std::vector<LHCData> lhcdataList;
            while (it->goToNext()) {
                const cool::IObject& object = it->currentRef();
                if (object.since() <= options.timestampNs && object.until() > options.timestampNs) {
                    lhcdataList.push_back(readBeamData(object));
                }
            }
            if (lhcdataList.empty()) {
                LOG_ERROR("No LHC data found for run {} at timestamp {} ns", options.run, options.timestampNs);
                throw std::runtime_error("No LHC data found for the given run and timestamp");
            } else if (lhcdataList.size() > 1) {
                LOG_WARN("Multiple LHC data entries found for run {} at timestamp {} ns. Using the first one.", options.run, options.timestampNs);
            }
            lhcdata = lhcdataList.front();
            
            // Update cache
            m_cached_data = lhcdata;
            m_cached_since = lhcdata.since;
            m_cached_until = lhcdata.until;
            m_cache_valid = true;
        }       
    }
    catch (const cool::Exception& e) {
        LOG_ERROR("COOL exception while reading LHC data: {}", e.what());
        throw std::runtime_error("COOL exception while reading LHC data: " + std::string(e.what()));
    } catch (const std::exception& e) {
        LOG_ERROR("Error while reading LHC data: {}", e.what());
        throw std::runtime_error("Error while reading LHC data: " + std::string(e.what()));
    }
    evt.put(m_cfg.out_lhcdata_key, std::move(lhcdata));
}

void LHCDataReaderAlg::parse_cfg(const YAML::Node& n) {
    m_cfg.in_rawdata_key = get_or<std::string>(n, "in_rawdata_key", m_cfg.in_rawdata_key);
    m_cfg.out_lhcdata_key = get_or<std::string>(n, "out_lhcdata_key", m_cfg.out_lhcdata_key);
    m_cfg.verbose = get_or<bool>(n, "verbose", m_cfg.verbose);
}

void LHCDataReaderAlg::initialize() {
    try {
        m_normalizedDbId = normalizeDbId(kDefaultDb);
        cool::IDatabaseSvc& dbSvc = cool::DatabaseSvcFactory::databaseService();
        m_db = dbSvc.openDatabase(m_normalizedDbId, true);
        m_folder = m_db->getFolder(kFolder);
        if (m_cfg.verbose) {
            LOG_INFO("LHCDataReaderAlg initialized with database: {}", m_normalizedDbId);
        }
    } catch (const cool::Exception& e) {
        LOG_ERROR("COOL exception during initialization: {}", e.what());
        throw std::runtime_error("COOL exception during initialization: " + std::string(e.what()));
    } catch (const std::exception& e) {
        LOG_ERROR("Error during initialization: {}", e.what());
        throw;
    }
}

void LHCDataReaderAlg::finalize() {
    m_folder = nullptr;
    m_db = nullptr;
    m_cache_valid = false;
    if (m_cfg.verbose) {
        LOG_INFO("LHCDataReaderAlg finalized");
    }
}

} // namespace AHCALRecoAlg

AHCAL_REGISTER_ALG(AHCALRecoAlg::LHCDataReaderAlg, "LHCDataReaderAlg")