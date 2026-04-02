// #include "reco_alg/module/MuonKFAlg/MuonKFAlg.hpp"
#include "adc_to_energy/AdcToEnergyReadTTreeAlg.hpp"
#include "reco_alg/module/TrackFitAlg/TrackFitAlg.hpp"

#include "common/Logger.hpp"
#include "common/edm/EDM.hpp"
#include "common/IAlg.hpp"
#include "common/RunContext.hpp"
#include "common/AlgFactory.hpp"
#include "common/config/ParseRunConfig.hpp"
#include "common/config/YAMLUtil.hpp"
#include "IO/reader/RootRawHitReader.hpp"
#include "IO/reader/BinaryRawHitReader.hpp"
#include "IO/reader/SimHitReader.hpp"
#if FAIR_HAS_GEANT4
#include "IO/reader/Geant4SimReader.hpp"
#endif
#include "IO/writer/RootWriterAlg.hpp"
#include "IO/writer/WriterRegistry.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

std::string customPad(const std::string& str, size_t totalLength) {
    std::ostringstream oss;
    if (str.length() < totalLength) {
        oss << std::string(totalLength - str.length(), '0') << str;
    } else {
        oss << str;
    }
    return oss.str();
}
using namespace AHCALRecoAlg;
int main(int argc, char* argv[]) {
    if (argc < 2) {
        LOG_ERROR("Usage: {} <config_yaml> [-i <input text file>]", argv[0]);
        return 1;
    }
    YAML::Node config = YAML::LoadFile(argv[1]);
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
    std::string input_text_file;
    bool use_input_text_file = false;
    for (int i = 2; i < argc; ++i) {
        if (std::string(argv[i]) == "-i" && i + 1 < argc) {
            input_text_file = argv[i + 1];
            use_input_text_file = true;
            LOG_INFO("Input text file specified: {}", input_text_file);
            break;
        }
    }
    std::cout << "AHCAL Application started." << std::endl;
    std::cout << "Logger initialized with level: " << ctx.config.log_level << " and log file: " << ctx.config.log_file << std::endl;
    std::string outputfile = ctx.config.output;
    std::cout << "Output file: " << outputfile << std::endl;
    int ninputs = 1;
    std::vector<std::string> input_files;
    std::vector<int> runNumbers;
    std::vector<int> poolIndexes;
    if (use_input_text_file) {
        std::ifstream infile(input_text_file, std::ios::in);
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            std::string filename;
            int runNumber = 0;
            int poolIndex = 0;
            if (!(iss >> filename >> runNumber >> poolIndex)) { break; }
            input_files.push_back(filename);
            runNumbers.push_back(runNumber);
            poolIndexes.push_back(poolIndex);
        }
        ninputs = input_files.size();
        LOG_INFO("Number of input files to process: {}", ninputs);
    }else {
        ninputs = 1;
        input_files.push_back(ctx.config.input);
        runNumbers.push_back(ctx.config.runNumber);
        poolIndexes.push_back(ctx.config.poolIndex);
    }
    LOG_INFO("AHCAL Application started.");
    LOG_INFO("RunConfig parsed successfully.");
    auto algs = build_pipeline(ctx, config);
    for (auto& alg : algs) {
        alg->initialize();
    }
    YAML::Node reader_config = require_node(config, "reader");
    const std::string type = require_string(reader_config, "type");
    const YAML::Node cfg = reader_config["cfg"] ? reader_config["cfg"] : YAML::Node(YAML::NodeType::Map);
    if (type == "RootRawHitReader") {
        // Initialize RootRawHitReader
        std::string input_key_hits = require_string(cfg, "out_rawhits_key");
        std::string input_key_tlu = require_string(cfg, "out_tlu_key");
        int nEvent = 0;
        for (int iinput = 0; iinput < ninputs; ++iinput) {
            ctx.config.input = input_files[iinput];
            ctx.config.runNumber = runNumbers[iinput];
            ctx.config.poolIndex = poolIndexes[iinput];
            ctx.config.output = outputfile;
            ctx.conditions = parse_condition_store(runNumbers[iinput]);
            LOG_INFO("Processing input file: {} (RunNumber: {}, PoolIndex: {})", ctx.config.input, ctx.config.runNumber, ctx.config.poolIndex);
            LOG_INFO("Inputs = {} / {}", iinput + 1, ninputs);
            RootRawHitReader rawHitReader(ctx.config.input, "Raw_Hit");
            LOG_INFO("RootRawHitReader created successfully.");
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
                    LOG_INFO("Processed {}/{} events.", nEvent, total_entries);
                }
                eventStore.clear();
            }
        }
    } else if (type == "BinaryRawHitReader") {
        // Initialize BinaryRawHitReader
        int nEvent = 0;
        std::string input_key_hits = require_string(cfg, "out_rawhits_key");
        std::string input_key_tlu = require_string(cfg, "out_tlu_key");
        // Long64_t total_entries = rawHitReader.entries();
        // LOG_INFO("Total entries in input file: {}", total_entries);
        for (int iinput = 0; iinput < ninputs; ++iinput) {
            ctx.config.input = input_files[iinput];
            ctx.config.runNumber = runNumbers[iinput];
            ctx.config.poolIndex = poolIndexes[iinput];
            ctx.config.output = outputfile;
            ctx.conditions = parse_condition_store(runNumbers[iinput]);
            LOG_INFO("Processing input file: {} (RunNumber: {}, PoolIndex: {})", ctx.config.input, ctx.config.runNumber, ctx.config.poolIndex);
            LOG_INFO("Inputs = {} / {}", iinput + 1, ninputs);
            BinaryRawHitReader rawHitReader(ctx.config.input);
            LOG_INFO("BinaryRawHitReader created successfully.");
            int nEvent_current = 0;
            while (true) {
                std::vector<AHCALRawHit> rawHits;
                AHCALTLURawData tluData;
                if (!rawHitReader.next(rawHits, tluData)) {
                    break; // No more events or error
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
                if (nEvent % 10000 == 0) {
                    LOG_INFO("Processed {}", nEvent);
                }
                eventStore.clear();
            }
            LOG_INFO("Finished processing input file: {} (RunNumber: {}, PoolIndex: {})", ctx.config.input, ctx.config.runNumber, ctx.config.poolIndex);
            LOG_INFO("Total events processed so far: {}", nEvent);
        }
    }else if (type == "SimHitReader") {
        // Initialize SimHitReader
        int nEvent = 0;
        std::string input_key_simhits = require_string(cfg, "out_simhits_key");
        std::string input_key_simdata = require_string(cfg, "out_simdata_key");
        for (int iinput = 0; iinput < ninputs; ++iinput) {
            ctx.config.input = input_files[iinput];
            ctx.config.runNumber = runNumbers[iinput];
            ctx.config.poolIndex = poolIndexes[iinput];
            ctx.config.output = outputfile;
            // ctx.conditions = parse_condition_store(runNumbers[iinput]);
            LOG_INFO("Processing input file: {} (RunNumber: {}, PoolIndex: {})", ctx.config.input, ctx.config.runNumber, ctx.config.poolIndex);
            LOG_INFO("Inputs = {} / {}", iinput + 1, ninputs);
            SimHitReader simHitReader(ctx.config.input);
            LOG_INFO("SimHitReader created successfully.");
            int nEvent_current = 0;
            while (true) {
                std::vector<AHCALSimHit> simHits;
                SimData simData;
                if (!simHitReader.next(simHits, simData)) {
                    break; // No more events or error
                }
                if (ctx.config.nEvents > 0 && nEvent >= ctx.config.nEvents) {
                    break; // Reached the maximum number of events to process
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
                if (nEvent % 10000 == 0) {
                    LOG_INFO("Processed {}", nEvent);
                }
                eventStore.clear();
            }
            LOG_INFO("Finished processing input file: {} (RunNumber: {}, PoolIndex: {})", ctx.config.input, ctx.config.runNumber, ctx.config.poolIndex);
            LOG_INFO("Total events processed so far: {}", nEvent);
        }
    }
#if FAIR_HAS_GEANT4
    else if (type == "Geant4SimReader" || type == "AHCALSimReader") {
        int nEvent = 0;
        std::string input_key_simhits = get_or<std::string>(cfg, "out_simhits_key", "SimHits");
        std::string input_key_simdata = get_or<std::string>(cfg, "out_simdata_key", "SimData");
        Geant4SimReader simReader(cfg);
        LOG_INFO("Geant4SimReader created successfully.");

        while (true) {
            std::vector<AHCALSimHit> simHits;
            SimData simData;
            if (!simReader.next(simHits, simData)) {
                break;
            }
            if (ctx.config.nEvents > 0 && nEvent >= ctx.config.nEvents) {
                break;
            }
            EventStore eventStore;
            eventStore.set_event_counter(nEvent);
            eventStore.put(input_key_simhits, std::move(simHits));
            eventStore.put(input_key_simdata, std::move(simData));
            for (auto& alg : algs) {
                alg->execute(eventStore);
            }
            nEvent++;
            if (nEvent % 10000 == 0) {
                LOG_INFO("Processed {}", nEvent);
            }
            eventStore.clear();
        }
        LOG_INFO("Finished Geant4 simulation events: {}", nEvent);
    }
#endif
    else if (type == "RootInput") {
        int nEvent = 0;
        for (int iinput = 0; iinput < ninputs; ++iinput) {
            ctx.config.input = input_files[iinput];
            ctx.config.runNumber = runNumbers[iinput];
            ctx.config.poolIndex = poolIndexes[iinput];
            ctx.config.output = outputfile;
            ctx.conditions = parse_condition_store(runNumbers[iinput]);
            LOG_INFO("Inputs = {} / {}", iinput + 1, ninputs);
            // Initialize RootInput reader
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
        LOG_ERROR("Unknown reader type specified in config.");
        return 1;
    }
    for (auto& alg : algs) {
        alg->finalize();
    }
    algs.clear();
    LOG_INFO("AHCAL Application finished.");
    std::cout << "AHCAL Application finished." << std::endl;
    return 0;
}
