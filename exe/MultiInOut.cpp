// #include "reco_alg/module/MuonKFAlg/MuonKFAlg.hpp"
#include "adc_to_energy/AdcToEnergyReadTTreeAlg.hpp"
#include "reco_alg/module/TrackFitAlg/TrackFitAlg.hpp"

#include "common/Logger.hpp"
#include "common/edm/EDM.hpp"
#include "common/IAlg.hpp"
#include "common/RunContext.hpp"
#include "common/AlgFactory.hpp"
#include "common/config/ParseRunConfig.hpp"
#include "IO/reader/RootRawHitReader.hpp"
#include "IO/reader/BinaryRawHitReader.hpp"
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
    std::string outputfile = ctx.config.output;
    int ninputs = 1;
    std::vector<std::string> input_files;
    std::vector<int> runNumbers;
    std::vector<int> poolIndexes;
    std::vector<std::string> output_files;
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
            std::string o = outputfile;
            const char* ext = ".root";
            constexpr std::size_t ext_len = 5;
            if (o.size() >= ext_len && o.compare(o.size() - ext_len, ext_len, ext) == 0) {
                o.erase(o.size() - ext_len);
            }

            output_files.push_back(
                o + "-" + customPad(std::to_string(runNumber), 6) + "-" +
                customPad(std::to_string(poolIndex), 5) + ".root"
            );
        }
        ninputs = input_files.size();
        LOG_INFO("Number of input files to process: {}", ninputs);
    }else {
        ninputs = 1;
        input_files.push_back(ctx.config.input);
        runNumbers.push_back(ctx.config.runNumber);
        poolIndexes.push_back(ctx.config.poolIndex);
        output_files.push_back(outputfile);
    }
    LOG_INFO("AHCAL Application started.");
    LOG_INFO("RunConfig parsed successfully.");
    for (int iinput = 0; iinput < ninputs; ++iinput) {
        ctx.config.input = input_files[iinput];
        ctx.config.runNumber = runNumbers[iinput];
        ctx.config.poolIndex = poolIndexes[iinput];
        ctx.config.output = output_files[iinput];
        LOG_INFO("Processing input file: {} (RunNumber: {}, PoolIndex: {})", ctx.config.input, ctx.config.runNumber, ctx.config.poolIndex);
        LOG_INFO("Inputs = {} / {}", iinput + 1, ninputs);
        auto algs = build_pipeline(ctx, config);
        for (auto& alg : algs) {
            alg->initialize();
        }
        YAML::Node reader_config = require_node(config, "reader");
        const std::string type = require_string(reader_config, "type");
        const YAML::Node cfg = reader_config["cfg"] ? reader_config["cfg"] : YAML::Node(YAML::NodeType::Map);

        if (type == "RootRawHitReader") {
            // Initialize RootRawHitReader
            RootRawHitReader rawHitReader(ctx.config.input, "Raw_Hit");
            LOG_INFO("RootRawHitReader created successfully.");
            int nEvent = 0;
            std::string input_key_hits = require_string(cfg, "out_rawhits_key");
            std::string input_key_tlu = require_string(cfg, "out_tlu_key");
            Long64_t total_entries = rawHitReader.entries();
            LOG_INFO("Total entries in input file: {}", total_entries);
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
                eventStore.put(input_key_hits, std::move(rawHits));
                eventStore.put(input_key_tlu, std::move(tluData));
                for (auto& alg : algs) {
                    alg->execute(eventStore);
                }
                nEvent++;
                if (nEvent % 10000 == 0) {
                    LOG_INFO("Processed {}/{} events.", nEvent, total_entries);
                }
                eventStore.clear();
            }
        } else if (type == "BinaryRawHitReader") {
            // Initialize BinaryRawHitReader
            BinaryRawHitReader rawHitReader(ctx.config.input);
            LOG_INFO("BinaryRawHitReader created successfully.");
            int nEvent = 0;
            std::string input_key_hits = require_string(cfg, "out_rawhits_key");
            std::string input_key_tlu = require_string(cfg, "out_tlu_key");
            // Long64_t total_entries = rawHitReader.entries();
            // LOG_INFO("Total entries in input file: {}", total_entries);
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
                eventStore.put(input_key_hits, std::move(rawHits));
                eventStore.put(input_key_tlu, std::move(tluData));
                for (auto& alg : algs) {
                    alg->execute(eventStore);
                }
                nEvent++;
                if (nEvent % 10000 == 0) {
                    LOG_INFO("Processed {}", nEvent);
                }
                eventStore.clear();
            }
        }else if (type == "RootInput") {
            // Initialize RootInput reader
            RootInput in(ctx.config.input, "events");
            LOG_INFO("RootInput reader created successfully.");
            ReaderRegistry rr = parse_reader_registry(cfg);
            Long64_t total_entries = in.entries();
            LOG_INFO("Total entries in input file: {}", total_entries);
            int nEvent = 0;
            while (true) {
                if (!in.next()) {
                    break; // No more events
                }
                if (ctx.config.nEvents > 0 && nEvent >= ctx.config.nEvents) {
                    break; // Reached the maximum number of events to process
                }
                EventStore eventStore;
                readandput(cfg, eventStore, rr, in);
                for (auto& alg : algs) {
                    alg->execute(eventStore);
                }
                nEvent++;
                if (nEvent % 10000 == 0) {
                    LOG_INFO("Processed {}/{} events.", nEvent, total_entries);
                }
                eventStore.clear();
            }
        } else {
            LOG_ERROR("Unknown reader type specified in config.");
            return 1;
        }
        for (auto& alg : algs) {
            alg->finalize();
        }
        algs.clear();
    }
    LOG_INFO("AHCAL Application finished.");
    std::cout << "AHCAL Application finished." << std::endl;
    return 0;
}