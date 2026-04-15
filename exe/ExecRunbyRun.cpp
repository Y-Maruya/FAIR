// #include "reco_alg/module/MuonKFAlg/MuonKFAlg.hpp"
#include "adc_to_energy/AdcToEnergyReadTTreeAlg.hpp"
#include "reco_alg/module/TrackFitAlg/TrackFitAlg.hpp"
#include "analysis/VetoAnaAlg/VetoAnaAlg.hpp"
#include "common/Logger.hpp"
#include "common/edm/EDM.hpp"
#include "common/IAlg.hpp"
#include "common/RunContext.hpp"
#include "common/AlgFactory.hpp"
#include "common/config/ParseRunConfig.hpp"
#include "IO/reader/RootRawHitReader.hpp"
#include "IO/reader/BinaryRawHitReader.hpp"
#include "IO/reader/SimHitReader.hpp"
#include "IO/writer/RootWriterAlg.hpp"
#include "IO/writer/WriterRegistry.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <filesystem>
#include <glob.h>


std::string customPad(const std::string& str, size_t totalLength) {
    std::ostringstream oss;
    if (str.length() < totalLength) {
        oss << std::string(totalLength - str.length(), '0') << str;
    } else {
        oss << str;
    }
    return oss.str();
}
std::vector<int> parseRunList(const char *csv) {
  std::vector<int> runs;
  if (!csv) return runs;
  std::stringstream ss(csv);
  std::string token;
  while (std::getline(ss, token, ',')) {
    token.erase(std::remove_if(token.begin(), token.end(), ::isspace), token.end());
    if (token.empty()) continue;
    runs.push_back(std::atoi(token.c_str()));
  }
  return runs;
}

void replace_config_run_number(const std::string& config_path, int run_number) {
    std::system(("cp " + config_path + " /tmp/temp_config.yaml").c_str());
    std::system(("python3 scripts/changepath.py /tmp/temp_config.yaml -r " + std::to_string(run_number)).c_str());
    if (!std::filesystem::exists("/tmp/temp_config.yaml")) {
        LOG_ERROR("Failed to create temporary config file with updated run number.");
    }
}

std::vector<std::string> glob_inputs(const std::string& pattern) {
    glob_t glob_result;
    glob(pattern.c_str(), GLOB_TILDE, nullptr, &glob_result);
    std::vector<std::string> files;
    for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
        files.push_back(std::string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return files;
}
using namespace AHCALRecoAlg;
int main(int argc, char* argv[]) {
    if (argc < 2) {
        LOG_ERROR("Usage: {} <config_yaml> -r <run_csv> [-m <max_pool_size>]", argv[0]);
        return 1;
    }
    std::string input_run_csv;
    int max_pool_size = -1;
    for (int i = 2; i < argc; ++i) {
        if (std::string(argv[i]) == "-r") {
            input_run_csv = argv[i + 1];
            std::cout << "Input run CSV file specified: " << input_run_csv << std::endl;
        }else if (std::string(argv[i]) == "-m") {
            max_pool_size = std::atoi(argv[i + 1]);
            std::cout << "Max pool size specified: " << max_pool_size << std::endl;
        }
    }
    std::cout << "AHCAL Application started." << std::endl;
    std::vector<int> runNumbers = parseRunList(input_run_csv.c_str());
    if (runNumbers.empty()) {
        LOG_ERROR("No run numbers found in input CSV: {}", input_run_csv);
        return 1;
    }
    for (size_t i = 0; i < runNumbers.size(); ++i) {
        int run = runNumbers[i];
        std::cout << "Run number: " << run << "," << i <<"/" << runNumbers.size()<<" is processing..." << std::endl;
        replace_config_run_number(argv[1], run);
        YAML::Node config = YAML::LoadFile("/tmp/temp_config.yaml");
        RunContext ctx;
        ctx.config = parse_run_config(config);
        if (ctx.config.log_level == "debug" || ctx.config.log_level == "DEBUG") {
            FAIR::init_logger("AHCALApp", ctx.config.log_file, spdlog::level::debug);
        } else if (ctx.config.log_level == "info" || ctx.config.log_level == "INFO") {
            FAIR::init_logger("AHCALApp", ctx.config.log_file, spdlog::level::info);
        } else if (ctx.config.log_level == "warn" || ctx.config.log_level == "WARN") {
            FAIR::init_logger("AHCALApp", ctx.config.log_file, spdlog::level::warn);
        } else if (ctx.config.log_level == "error" || ctx.config.log_level == "ERROR") {
            FAIR::init_logger("AHCALApp", ctx.config.log_file, spdlog::level::err);
        } else {
            FAIR::init_logger("AHCALApp", ctx.config.log_file, spdlog::level::info);
        } 
        LOG_INFO("RunConfig parsed successfully.");
        auto algs = build_pipeline(ctx, config);
        for (auto& alg : algs) {
            alg->initialize();
        }
        std::vector <std::string> raw_inputs_path = glob_inputs(ctx.config.input);
        if (raw_inputs_path.empty()) {
            LOG_ERROR("No input files found for run {} with path pattern: {}", ctx.config.runNumber, ctx.config.input);
            continue;
        }
        YAML::Node reader_config = require_node(config, "reader");
        const std::string type = require_string(reader_config, "type");
        const YAML::Node cfg = reader_config["cfg"] ? reader_config["cfg"] : YAML::Node(YAML::NodeType::Map);
        ctx.conditions = parse_condition_store(ctx.config.runNumber);
        if (type == "RootRawHitReader") {
            // Initialize RootRawHitReader
            std::string input_key_hits = require_string(cfg, "out_rawhits_key");
            std::string input_key_tlu = require_string(cfg, "out_tlu_key");
            int nEvent = 0;
            int nInputs = 0;
            for (std::string input_file : raw_inputs_path) {
                nInputs++;
                if (max_pool_size > 0 && nInputs > max_pool_size) {
                    break;
                }
                ctx.config.input = input_file;
                RootRawHitReader rawHitReader(ctx.config.input, "Raw_Hit");
                LOG_INFO("RootRawHitReader created successfully for file: {}", ctx.config.input);
                int nEvent_current = 0;
                Long64_t total_entries = rawHitReader.entries();
                LOG_INFO("Total entries in input file {}: {}", ctx.config.input, total_entries);
                while (true) {
                    std::vector<AHCALRawHit> rawHits;
                    AHCALTLURawData tluData;
                    if (!rawHitReader.next(rawHits, tluData)) {
                        break; // No more events
                    }
                    if (ctx.config.nEvents > 0 && nEvent >= ctx.config.nEvents) {
                        break; // Reached the maximum number of events to process
                    }
                    EventStore eventStore;
                    eventStore.set_event_counter(nEvent);
                    eventStore.put(input_key_hits, std::move(rawHits));
                    eventStore.put(input_key_tlu, std::move(tluData));
                    for (auto& alg : algs) {
                        alg->execute(eventStore);
                    }
                    nEvent++;
                    nEvent_current++;
                    if (nEvent_current % 10000 == 0) {
                        LOG_INFO("Processed {}/{} events in file {}.", nEvent_current, total_entries, ctx.config.input);
                    }
                    eventStore.clear();
                }
            }
        } else if (type == "BinaryRawHitReader") {
            // Initialize BinaryRawHitReader
            int nEvent = 0;
            int nInputs = 0;
            std::string input_key_hits = require_string(cfg, "out_rawhits_key");
            std::string input_key_tlu = require_string(cfg, "out_tlu_key");
            for (std::string input_file : raw_inputs_path) {
                nInputs++;
                if (max_pool_size > 0 && nInputs > max_pool_size) { 
                    break;
                }
                ctx.config.input = input_file;
                BinaryRawHitReader rawHitReader(ctx.config.input);
                LOG_INFO("BinaryRawHitReader created successfully for file: {}", ctx.config.input);
                int nEvent_current = 0;
                while (true) {
                    std::vector<AHCALRawHit> rawHits;
                    AHCALTLURawData tluData;
                    if (!rawHitReader.next(rawHits, tluData)) {
                        break; // No more events or error
                    }
                    EventStore eventStore;
                    eventStore.set_event_counter(nEvent);
                    eventStore.put(input_key_hits, std::move(rawHits));
                    eventStore.put(input_key_tlu, std::move(tluData));
                    for (auto& alg : algs) {
                        alg->execute(eventStore);
                    }
                    nEvent++;
                    nEvent_current++;
                    if (nEvent_current % 10000 == 0) {
                        LOG_INFO("Processed {} events in file {}.", nEvent_current, ctx.config.input);
                    }
                    eventStore.clear();
                }
            }
        } else if (type == "SimHitReader") {
            // Initialize SimHitReader
            int nEvent = 0;
            int nInputs = 0;
            std::string input_key_simhits = require_string(cfg, "out_simhits_key");
            std::string input_key_simdata = require_string(cfg, "out_simdata_key");
            for (std::string input_file : raw_inputs_path) {
                nInputs++;
                if (max_pool_size > 0 && nInputs > max_pool_size) {
                    break;
                }
                ctx.config.input = input_file;
                SimHitReader rawHitReader(ctx.config.input);
                Long64_t total_entries = rawHitReader.entries();
                LOG_INFO("SimHitReader created successfully for file: {}", ctx.config.input);
                int nEvent_current = 0;
                while (true) {
                    std::vector<AHCALSimHit> simHits;
                    SimData simData;
                    if (!rawHitReader.next(simHits, simData)) {
                        break; // No more events or error
                    }
                    EventStore eventStore;
                    eventStore.set_event_counter(nEvent);
                    eventStore.put(input_key_simhits, std::move(simHits));
                    eventStore.put(input_key_simdata, std::move(simData));
                    for (auto& alg : algs) {
                        alg->execute(eventStore);
                    }
                    nEvent++;
                    nEvent_current++;
                    if (nEvent_current % 10000 == 0) {
                        LOG_INFO("Processed {}/{} events in file {}.", nEvent_current, total_entries, ctx.config.input);
                    }
                    eventStore.clear();
                }
            }
        } else if (type == "RootInput") {
            int nEvent = 0;
            int nInputs = 0;
            for (std::string input_file : raw_inputs_path) {
                nInputs++;
                if (max_pool_size > 0 && nInputs > max_pool_size) {
                    break;
                }
                ctx.config.input = input_file;
                RootInput in(ctx.config.input, "events");
                LOG_INFO("RootInput reader created successfully.");
                ReaderRegistry rr = parse_reader_registry(cfg);
                Long64_t total_entries = in.entries();
                LOG_INFO("Total entries in input file: {}", total_entries);
                int nEvent_current = 0;
                while (true) {
                    if (!in.next()) {
                        break; // No more events
                    }
                    if (ctx.config.nEvents > 0 && nEvent >= ctx.config.nEvents) {
                        break; // Reached the maximum number of events to process
                    }
                    EventStore eventStore;
                    eventStore.set_event_counter(nEvent);
                    readandput(cfg, eventStore, rr, in);
                    for (auto& alg : algs) {
                        alg->execute(eventStore);
                    }
                    nEvent++;
                    nEvent_current++;
                    if (nEvent_current % 10000 == 0) {
                        LOG_INFO("Processed {}/{} events.", nEvent_current, total_entries);
                    }
                    eventStore.clear();
                }
                LOG_INFO("Finished processing input file: {} (RunNumber: {}, PoolIndex: {})", ctx.config.input, ctx.config.runNumber, ctx.config.poolIndex);
                LOG_INFO("Total events processed so far: {}", nEvent);
            }
        } else {
            LOG_ERROR("Unsupported reader type: {}", type);
            return 1;
        }
        for (auto& alg : algs) {
            alg->finalize();
        }
        algs.clear();
        LOG_INFO("Finished processing run {}.", ctx.config.runNumber);
    }
    LOG_INFO("AHCAL Application finished.");
    return 0;
}