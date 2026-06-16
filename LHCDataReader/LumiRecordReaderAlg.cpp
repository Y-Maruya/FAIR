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

#include "LumiRecordReaderAlg.hpp"
#include "common/AlgRegistry.hpp"
#include "common/Logger.hpp"
#include "common/edm/LumiRecord.hpp"
#include "common/edm/RawData.hpp"
#include "common/config/YAMLUtil.hpp"


namespace AHCALRecoAlg {
namespace {

constexpr const char* kDefaultDb =
  "COOLOFL_TRIGGER/CONDBR2";
constexpr const char* kFolder = "/TRIGGER/OFLLUMI/LumiAccounting";

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

LumiRecord readBeamData(const cool::IObject& object) {
    const cool::IRecord& payload = object.payload();
    LumiRecord data;
    data.since = object.since();
    data.until = object.until();
    data.instLumi = static_cast<double>(payload["InstLumi"].data<float>());
    data.AvrageEventsPerBunchCrossing = static_cast<double>(payload["AvEvtsPerBX"].data<float>());
    data.runNumber = static_cast<Long64_t>(payload["Run"].data<uint32_t>());
    data.LumiBlockNumber = static_cast<Long64_t>(payload["LumiBlock"].data<uint32_t>());
    data.runlb = (data.runNumber << 32) | data.LumiBlockNumber;
    data.status = static_cast<int>(payload["Status"].data<uint32_t>());
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

} // namespace

void LumiRecordReaderAlg::execute(EventStore& evt) {
    const auto& rawdata = evt.get<AHCALTLURawData>(m_cfg.in_rawdata_key);
    LumiRecord lumiRecord;
    try {
        const Options options = parseArgsFromEvent(rawdata, m_cfg.verbose);
        
        // Check cache: if timestamp falls in cached IOV range, reuse the result
        if (m_cache_valid && 
            options.timestampNs >= m_cached_since && 
            options.timestampNs < m_cached_until) {
            lumiRecord = m_cached_data;
            if (m_cfg.verbose) {
                LOG_DEBUG("Using cached LHC data for timestamp {}", options.timestampNs);
            }
        } else {
            // Cache miss: query COOL database
            const cool::ValidityKey queryUntil =
                options.timestampNs == cool::ValidityKeyMax
                    ? options.timestampNs
                    : options.timestampNs + 1;
            cool::IObjectIteratorPtr it =
                m_folder->browseObjects(
                    options.timestampNs,
                    queryUntil,
                    cool::ChannelSelection::all(),
                    m_cfg.lumi_tag);
            std::vector<LumiRecord> LumiRecordList;
            while (it->goToNext()) {
                const cool::IObject& object = it->currentRef();
                if (object.since() <= options.timestampNs && object.until() > options.timestampNs) {
                    LumiRecordList.push_back(readBeamData(object));
                }
            }
            if (LumiRecordList.empty()) {
                LOG_ERROR("No LHC data found for run {} at timestamp {} ns with lumi tag {}",
                          options.run, options.timestampNs, m_cfg.lumi_tag);
                throw std::runtime_error("No LHC data found for the given run and timestamp");
            } else if (LumiRecordList.size() > 1) {
                LOG_WARN("Multiple LHC data entries found for run {} at timestamp {} ns. Using the first one.", options.run, options.timestampNs);
            }
            lumiRecord = LumiRecordList.front();
            
            // Update cache
            m_cached_data = lumiRecord;
            m_cached_since = lumiRecord.since;
            m_cached_until = lumiRecord.until;
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
    evt.put(m_cfg.out_lumirecord_key, std::move(lumiRecord));
}

void LumiRecordReaderAlg::parse_cfg(const YAML::Node& n) {
    m_cfg.in_rawdata_key = get_or<std::string>(n, "in_rawdata_key", m_cfg.in_rawdata_key);
    m_cfg.out_lumirecord_key = get_or<std::string>(n, "out_lumirecord_key", m_cfg.out_lumirecord_key);
    m_cfg.lumi_tag = get_or<std::string>(n, "lumi_tag", m_cfg.lumi_tag);
    m_cfg.verbose = get_or<bool>(n, "verbose", m_cfg.verbose);
}

void LumiRecordReaderAlg::initialize() {
    try {
        m_normalizedDbId = normalizeDbId(kDefaultDb);
        cool::IDatabaseSvc& dbSvc = cool::DatabaseSvcFactory::databaseService();
        m_db = dbSvc.openDatabase(m_normalizedDbId, true);
        m_folder = m_db->getFolder(kFolder);
        if (m_cfg.verbose) {
            LOG_INFO("LumiRecordReaderAlg initialized with database: {}, tag: {}",
                     m_normalizedDbId, m_cfg.lumi_tag);
        }
    } catch (const cool::Exception& e) {
        LOG_ERROR("COOL exception during initialization: {}", e.what());
        throw std::runtime_error("COOL exception during initialization: " + std::string(e.what()));
    } catch (const std::exception& e) {
        LOG_ERROR("Error during initialization: {}", e.what());
        throw;
    }
}

void LumiRecordReaderAlg::finalize() {
    m_folder = nullptr;
    m_db = nullptr;
    m_cache_valid = false;
    if (m_cfg.verbose) {
        LOG_INFO("LumiRecordReaderAlg finalized");
    }
}

} // namespace AHCALRecoAlg

AHCAL_REGISTER_ALG(AHCALRecoAlg::LumiRecordReaderAlg, "LumiRecordReaderAlg")
