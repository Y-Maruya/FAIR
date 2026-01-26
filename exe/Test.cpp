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
#include "IO/writer/RootWriterAlg.hpp"
#include "IO/writer/WriterRegistry.hpp"
#include <iostream>
using namespace AHCALRecoAlg;
int main(int argc, char* argv[]) {
    if (argc < 2) {
        LOG_ERROR("Usage: {} <config_yaml>", argv[0]);
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
    LOG_INFO("AHCAL Application started.");
    LOG_INFO("RunConfig parsed successfully.");
    auto algs = build_pipeline(ctx, config);
    if (require_string(config,"reader") == "RootRawHitReader") {
        // Initialize RootRawHitReader
        RootRawHitReader rawHitReader(ctx.config.input, "Raw_Hit");
        LOG_INFO("RootRawHitReader created successfully.");
        while (true) {
            std::vector<AHCALRawHit> rawHits;
            if (!rawHitReader.next(rawHits)) {
                break; // No more events
            }
            EventStore eventStore;
            eventStore.put("RawHits", std::move(rawHits));
            for (auto& alg : algs) {
                alg->execute(eventStore);
            }
            // Here you can retrieve and process the final results from eventStore
            // LOG_INFO("Processed one event.");
            // LOG_INFO("Fitted track has {} hits.", eventStore.get<SimpleFittedTrack>("FittedTrack").nTotalHits);
            
        }
    }else if (require_string(config,"reader") == "SomeOtherReader") {
        // Implement other readers if needed
    } else {
        LOG_ERROR("Unknown reader type specified in config.");
        return 1;
    }
    LOG_INFO("AHCAL Application finished.");
    return 0;
}