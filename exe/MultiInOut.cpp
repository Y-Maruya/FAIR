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
void replaceRunNumberPoolIndexInFile(const std::string& filename, int runnumber, int poolindex) {
    // read file
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return;
    }

    std::stringstream buffer;
    buffer << infile.rdbuf();
    std::string content = buffer.str();
    infile.close();

    // replacement string
    std::string placeholder = "###runnumber###";
    std::string runStr = std::to_string(runnumber);

    // replace all occurrences
    size_t pos = 0;
    while ((pos = content.find(placeholder, pos)) != std::string::npos) {
        content.replace(pos, placeholder.length(), runStr);
        pos += runStr.length();
    }

    // replacement string
    std::string placeholder_pool = "###poolindex###";
    std::string poolStr = std::to_string(poolindex);

    // replace all occurrences
    pos = 0;
    while ((pos = content.find(placeholder_pool, pos)) != std::string::npos) {
        content.replace(pos, placeholder_pool.length(), poolStr);
        pos += poolStr.length();
    }
    // write back
    std::ofstream outfile("/tmp/"+std::to_string(runnumber)+".yaml");
    outfile << content;
}
using namespace AHCALRecoAlg;
int main(int argc, char* argv[]) {
    if (argc < 2) {
        LOG_ERROR("Usage: {} <config_yaml> [-i <input text file>]", argv[0]);
        return 1;
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
    // std::cout << "Logger initialized with level: " << ctx.config.log_level << " and log file: " << ctx.config.log_file << std::endl;
    // std::string outputfile = ctx.config.output;
    // std::cout << "Output file: " << outputfile << std::endl;
    int ninputs = 1;
    std::vector<std::string> input_files;
    std::vector<int> runNumbers;
    std::vector<int> poolIndexes;
    // std::vector<std::string> output_files;
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
            // std::string o = outputfile;
            // const char* ext = ".root";
            // constexpr std::size_t ext_len = 5;
            // if (o.size() >= ext_len && o.compare(o.size() - ext_len, ext_len, ext) == 0) {
            //     o.erase(o.size() - ext_len);
            // }

            // output_files.push_back(
            //     o + "-" + customPad(std::to_string(runNumber), 6) + "-" +
            //     customPad(std::to_string(poolIndex), 5) + ".root"
            // );
        }
        ninputs = input_files.size();
        LOG_INFO("Number of input files to process: {}", ninputs);
    }else {
        ninputs = 1;
        // input_files.push_back(ctx.config.input);
        // runNumbers.push_back(ctx.config.runNumber);
        // poolIndexes.push_back(ctx.config.poolIndex);
        // output_files.push_back(outputfile);
    }
    LOG_INFO("AHCAL Application started.");
    LOG_INFO("RunConfig parsed successfully.");
    
    for (int iinput = 0; iinput < ninputs; ++iinput) {
        replaceRunNumberPoolIndexInFile(argv[1],runNumbers[iinput],poolIndexes[iinput]);
        YAML::Node config_i = YAML::LoadFile("/tmp/"+std::to_string(runNumbers[iinput])+".yaml");
        RunContext ctx_i;
        ctx_i.config = parse_run_config(config_i);
        if (ctx_i.config.log_level == "debug" || ctx_i.config.log_level == "DEBUG") {
            FAIR::init_logger("AHCALApp", ctx_i.config.log_file, spdlog::level::debug);
        } else if (ctx_i.config.log_level == "info" || ctx_i.config.log_level == "INFO") {
            FAIR::init_logger("AHCALApp", ctx_i.config.log_file, spdlog::level::info);
        } else if (ctx_i.config.log_level == "warn" || ctx_i.config.log_level == "WARN") {
            FAIR::init_logger("AHCALApp", ctx_i.config.log_file, spdlog::level::warn);
        } else if (ctx_i.config.log_level == "error" || ctx_i.config.log_level == "ERROR") {
            FAIR::init_logger("AHCALApp", ctx_i.config.log_file, spdlog::level::err);
        } else {
            FAIR::init_logger("AHCALApp", ctx_i.config.log_file, spdlog::level::info);
        } 
        if (ninputs!=1){
            ctx_i.config.input = input_files[iinput];
            ctx_i.config.runNumber = runNumbers[iinput];
            ctx_i.config.poolIndex = poolIndexes[iinput];
        }
        // ctx_i.config.output = output_files[iinput];
        ctx_i.conditions = parse_condition_store(runNumbers[iinput]);
        LOG_INFO("Processing input file: {} (RunNumber: {}, PoolIndex: {})", ctx_i.config.input, ctx_i.config.runNumber, ctx_i.config.poolIndex);
        LOG_INFO("Inputs = {} / {}", iinput + 1, ninputs);
        auto algs = build_pipeline(ctx_i, config_i);
        for (auto& alg : algs) {
            alg->initialize();
        }
        YAML::Node reader_config = require_node(config_i, "reader");
        const std::string type = require_string(reader_config, "type");
        const YAML::Node cfg = reader_config["cfg"] ? reader_config["cfg"] : YAML::Node(YAML::NodeType::Map);

        if (type == "RootRawHitReader") {
            // Initialize RootRawHitReader
            RootRawHitReader rawHitReader(ctx_i.config.input, "Raw_Hit");
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
                if (ctx_i.config.nEvents > 0 && nEvent >= ctx_i.config.nEvents) {
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
                if (nEvent % 10000 == 0) {
                    LOG_INFO("Processed {}/{} events.", nEvent, total_entries);
                }
                eventStore.clear();
            }
        } else if (type == "BinaryRawHitReader") {
            // Initialize BinaryRawHitReader
            BinaryRawHitReader rawHitReader(ctx_i.config.input);
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
                if (ctx_i.config.nEvents > 0 && nEvent >= ctx_i.config.nEvents) {
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
                if (nEvent % 10000 == 0) {
                    LOG_INFO("Processed {}", nEvent);
                }
                eventStore.clear();
            }
        }else if (type == "RootInput") {
            // Initialize RootInput reader
            RootInput in(ctx_i.config.input, "events");
            LOG_INFO("RootInput reader created successfully.");
            ReaderRegistry rr = parse_reader_registry(cfg);
            Long64_t total_entries = in.entries();
            LOG_INFO("Total entries in input file: {}", total_entries);
            int nEvent = 0;
            while (true) {
                if (!in.next()) {
                    break; // No more events
                }
                if (ctx_i.config.nEvents > 0 && nEvent >= ctx_i.config.nEvents) {
                    break; // Reached the maximum number of events to process
                }
                long long event_counter = nEvent;
                if (in.has_branch("EventMeta.event_counter")) {
                    event_counter = *in.get_or_make_address<long long>("EventMeta.event_counter");
                }
                if (in.has_run_context()) {
                    in.apply_run_context_for_event(event_counter, ctx_i);
                }
                EventStore eventStore;
                eventStore.set_event_counter(event_counter);
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
