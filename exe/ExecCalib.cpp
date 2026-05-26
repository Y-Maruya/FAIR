// #include "reco_alg/module/MuonKFAlg/MuonKFAlg.hpp"
#include "adc_to_energy/AdcToEnergyReadTTreeAlg.hpp"
#include "reco_alg/module/TrackFitAlg/TrackFitAlg.hpp"
#include "analysis/VetoAnaAlg/VetoAnaAlg.hpp"
#include "common/Logger.hpp"
#include "common/edm/EDM.hpp"
#include "common/IAlg.hpp"
#include "common/RunContext.hpp"
#include "common/ContextPutter.hpp"
#include "common/AlgFactory.hpp"
#include "common/config/ParseRunConfig.hpp"
#include "IO/reader/RootRawHitReader.hpp"
#include "IO/reader/BinaryRawHitReader.hpp"
#include "IO/reader/SimHitReader.hpp"
#include "IO/writer/RootWriterAlg.hpp"
#include "IO/writer/WriterRegistry.hpp"
#include "CalibDBIO/TriggersReader/TriggersReader.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <stdexcept>
#include <filesystem>
#include <glob.h>
#include <sys/resource.h>


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
    token.erase(std::remove_if(token.begin(), token.end(), [](unsigned char c){ return std::isspace(c); }), token.end());
    if (token.empty()) continue;
    runs.push_back(std::atoi(token.c_str()));
  }
  return runs;
}

std::string shell_quote(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('\'');
    for (char c : s) {
        if (c == '\'') {
            out += "'\\''";
        } else {
            out.push_back(c);
        }
    }
    out.push_back('\'');
    return out;
}

void replace_config_run_numbers(const std::string& config_path, const std::string& run_numbers, std::string random_suffix = "") {
    std::filesystem::copy_file(config_path, "/tmp/temp_config" + random_suffix + ".yaml", std::filesystem::copy_options::overwrite_existing);
    const std::string cmd = "python3 scripts/changepath_multirun.py /tmp/temp_config" + random_suffix + ".yaml -rs " + shell_quote(run_numbers);
    const int rc = std::system(cmd.c_str());
    if (rc != 0) {
        LOG_ERROR("Failed to update run range via changepath_multirun.py. command={}, return_code={}", cmd, rc);
        throw std::runtime_error("Failed to update run range in temporary config");
    }
    if (!std::filesystem::exists("/tmp/temp_config" + random_suffix + ".yaml")) {
        LOG_ERROR("Failed to create temporary config file with updated run numbers.");
        throw std::runtime_error("Temporary config file not found after run range update");
    }
}
void replace_config_run_number(const std::string& config_path, int run_number, std::string random_suffix = "") {
    std::filesystem::copy_file(config_path, "/tmp/temp_config" + random_suffix + "_2.yaml", std::filesystem::copy_options::overwrite_existing);
    const std::string cmd = "python3 scripts/changepath.py /tmp/temp_config" + random_suffix + "_2.yaml -r " + std::to_string(run_number);
    const int rc = std::system(cmd.c_str());
    if (rc != 0) {
        LOG_ERROR("Failed to update run number via changepath.py. command={}, return_code={}", cmd, rc);
        throw std::runtime_error("Failed to update run number in temporary config");
    }
    if (!std::filesystem::exists("/tmp/temp_config" + random_suffix + "_2.yaml")) {
        LOG_ERROR("Failed to create temporary config file with updated run numbers.");
        throw std::runtime_error("Temporary config file not found after run update");
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
class AHCALRuns{
    public:
        AHCALRuns(std::string excluded_runs_file){
            init_excluded_runs(excluded_runs_file);
        }
        bool is_AHCAL_run(int run_number){
            // check the excluded runs first
            if (AHCALRuns::excluded_runs(run_number)) {
                return false;
            }
            const YAML::Node runlist = load_json_from_url(make_runinfo_url(run_number));
            const bool has_runnumber = has_node(runlist, "runnumber");
            if (!has_runnumber || runlist["runnumber"].as<int>() != run_number) {
                if (has_runnumber) {
                    LOG_ERROR("is_AHCAL_run: run number mismatch in response from {}: expected {}, got {}", make_runinfo_url(run_number), run_number, runlist["runnumber"].as<int>());
                } else {
                    LOG_ERROR("is_AHCAL_run: missing runnumber in response from {}. expected {}", make_runinfo_url(run_number), run_number);
                }
                throw std::runtime_error("Run number mismatch in response");
            }
            if (has_node(runlist, "type") && runlist["type"].as<std::string>() == "AHCAL") {
                return true;
            } else {
                return false;
            }
        }
        bool excluded_runs(int run_number){
            if (excluded_run_numbers_.empty()) {
                LOG_WARN("Excluded run numbers list is empty. No runs will be excluded.");
                return false;
            }
            if (std::find(excluded_run_numbers_.begin(), excluded_run_numbers_.end(), run_number) != excluded_run_numbers_.end()) {
                LOG_INFO("Run {} is in the excluded runs list.", run_number);
                return true;
            } else {
                return false;
            }
        }
        void init_excluded_runs(std::string excluded_runs_file){
            std::ifstream infile(excluded_runs_file);
            if (!infile.is_open()) {
                LOG_ERROR("Could not open {} to read excluded run information", excluded_runs_file);
                throw std::runtime_error("Could not open excluded_runs.txt to read excluded run information");
            }
            std::string line;
            while (std::getline(infile, line)) {
                std::istringstream iss(line);
                int run;
                if (!(iss >> run)) { continue; }
                excluded_run_numbers_.push_back(run);
            }
        }
        std::vector<int> exceed_veto_trigger_run(int start_run, int num_veto_events) {
            int run = start_run;
            int veto_events_count = 0;
            std::vector<int> exceeded_runs;
            while (veto_events_count < num_veto_events) {
                if (AHCALRuns::is_AHCAL_run(run)) {
                    CalibDBIO::TriggersReader reader(run);
                    reader.getTriggers();
                    int current_veto_events = reader.getTriggers().AnyVetoed;
                    veto_events_count += current_veto_events;
                    LOG_INFO("Run {}: AnyVetoed={}, Cumulative Veto Events={}", run, current_veto_events, veto_events_count);
                    exceeded_runs.push_back(run);
                }
                run++;
                if (run > start_run + 1000) { // safety check to prevent infinite loop
                    LOG_ERROR("Exceeded 1000 runs while looking for veto events. Stopping search.");
                    break;
                }
            }
            return exceeded_runs; // return the list of runs that exceed the veto trigger count
        }
    private:
        std::vector<int> excluded_run_numbers_;
};

using namespace AHCALRecoAlg;
int main(int argc, char* argv[]) {
    if (argc < 4) {
        LOG_ERROR("Usage: {} <config_yaml> -r <start_run> [-n <num_veto_events> -e <excluded_runs_file>]", argv[0]);
        return 1;
    }
    std::string input_run;
    int num_veto_events = 7000*18*18; // default value
    std::string excluded_runs_file = "excluded_runs.txt"; // default value
    for (int i = 2; i < argc; ++i) {
        if (i + 1 >= argc) {
            LOG_ERROR("Missing value for option {}", argv[i]);
            return 1;
        }
        if (std::string(argv[i]) == "-r") {
            input_run = argv[i + 1];
            std::cout << "Input run specified: " << input_run << std::endl;
        }else if (std::string(argv[i]) == "-n") {
            num_veto_events = std::atoi(argv[i + 1]);
            std::cout << "Number of veto events specified: " << num_veto_events << std::endl;
        }else if (std::string(argv[i]) == "-e") {
            excluded_runs_file = argv[i + 1];
            std::cout << "Excluded runs file specified: " << excluded_runs_file << std::endl;
        } else {
            LOG_ERROR("Unknown option: {}", argv[i]);
            return 1;
        }
        ++i;
    }
    std::cout << "AHCAL Application started." << std::endl;
    const std::vector<int> parsed_runs = parseRunList(input_run.c_str());
    if (parsed_runs.empty()) {
        LOG_ERROR("No valid start run provided with -r option.");
        return 1;
    }
    int start_runNumber = parsed_runs[0];
    std::cout << "Parsed start run number: " << start_runNumber << std::endl;
    if (start_runNumber <= 0) {
        LOG_ERROR("Invalid start run number: {}", start_runNumber);
        return 1;
    }
    if (num_veto_events <= 0) {
        LOG_ERROR("Invalid number of veto events: {}", num_veto_events);
        return 1;
    }
    AHCALRuns ahcal_runs(excluded_runs_file);
    std::vector<int> runNumbers = ahcal_runs.exceed_veto_trigger_run(start_runNumber, num_veto_events);
    if (runNumbers.empty()) {
        LOG_ERROR("No runs were selected from start run {} with requested veto event count {}.", start_runNumber, num_veto_events);
        return 1;
    }
    std::cout << "Total runs to process: " << runNumbers.size() << std::endl;
    int next_run = runNumbers.back() + 1;
    std::cout << "Next run after processing: " << next_run << std::endl;
    std::string run_numbers_str;
    run_numbers_str += std::to_string(runNumbers[0]);
    run_numbers_str += "-";
    run_numbers_str += std::to_string(runNumbers.back());
    std::cout << "Run numbers to process: " << run_numbers_str << std::endl;
    std::string random_suffix = "_" + std::to_string(std::rand() % 10000);
    replace_config_run_numbers(argv[1], run_numbers_str, random_suffix);
    replace_config_run_number("/tmp/temp_config" + random_suffix + ".yaml", start_runNumber, random_suffix);
    YAML::Node config = YAML::LoadFile("/tmp/temp_config" + random_suffix + "_2.yaml");
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
    std::cout << "LOG is initialized. Log file: " << ctx.config.log_file << std::endl;
    LOG_INFO("RunConfig parsed successfully.");
    auto algs = build_pipeline(ctx, config);
    for (auto& alg : algs) {
        alg->initialize();
    }
    for (size_t i = 0; i < runNumbers.size(); ++i) {
        replace_config_run_number("/tmp/temp_config" + random_suffix + ".yaml", runNumbers[i], random_suffix);
        YAML::Node config = YAML::LoadFile("/tmp/temp_config" + random_suffix + "_2.yaml");
        ctx.config = parse_run_config(config);
        ctx.conditions = parse_condition_store(ctx.config.runNumber);
        LOG_INFO("Processing run number: {}", ctx.config.runNumber);
        for (auto& alg : algs) {
            alg->init_by_run();
        }
        std::vector<std::string> raw_inputs_path = glob_inputs(ctx.config.input);
        if (raw_inputs_path.empty()) {
            LOG_WARN("No input files found for run {} with input pattern: {}", ctx.config.runNumber, ctx.config.input);
            continue;
        }
        YAML::Node reader_config = require_node(config, "reader");
        const std::string type = require_string(reader_config, "type");
        const YAML::Node cfg = reader_config["cfg"] ? reader_config["cfg"] : YAML::Node(YAML::NodeType::Map);

        if (type == "RootRawHitReader") {
            // Initialize RootRawHitReader
            std::string input_key_hits = require_string(cfg, "out_rawhits_key");
            std::string input_key_tlu = require_string(cfg, "out_tlu_key");
            int nEvent = 0;
            int nInputs = 0;
            for (std::string input_file : raw_inputs_path) {
                nInputs++;

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
                        LOG_INFO("Processed {} events in file {}.", nEvent_current, ctx.config.input);
                    }
                    if (nEvent % 1000 == 0) {
                        struct rusage usage;
                        getrusage(RUSAGE_SELF, &usage);
                        long mem_mb = usage.ru_maxrss / 1024;
                        LOG_INFO("Event {}: Memory = {} MB", nEvent, mem_mb);
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
        LOG_INFO("Finished processing run {}.", ctx.config.runNumber);
    }
    for (auto& alg : algs) {
        alg->finalize();
    }
    algs.clear();
    LOG_INFO("AHCAL Application finished.");
    return 0;
}