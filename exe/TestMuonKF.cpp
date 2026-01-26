#include "reco_alg/module/MuonKFAlg/MuonKFAlg.hpp"
#include "adc_to_energy/AdcToEnergyReadTTreeAlg.hpp"
#include "common/Logger.hpp"
#include "common/edm/EDM.hpp"
#include "common/IAlg.hpp"
#include "IO/reader/RootRawHitReader.hpp"
#include "IO/writer/RootWriterAlg.hpp"
#include "IO/writer/WriterRegistry.hpp"
#include <iostream>
using namespace AHCALRecoAlg;
int main(int argc, char* argv[]) {
    FAIR::init_logger("AHCALApp", "ahcal_app_MuonKF.log", spdlog::level::debug);
    LOG_INFO("AHCAL Application started.");
    std::vector<std::unique_ptr<IAlg>> algs;
    if (argc < 2) {
        LOG_ERROR("Usage: {} <cfg_file>", argv[0]);
        return 1;
    }
    // Initialize RootRawHitReader
    RootRawHitReader rawHitReader(argv[1], "Raw_Hit");
    LOG_INFO("RootRawHitReader created successfully.");
    
    // Initialize AdcToEnergyReadTTreeAlg
    std::unique_ptr<AdcToEnergyReadTTreeAlg> adcToEnergyAlg = std::make_unique<AdcToEnergyReadTTreeAlg>("RawHits", "RecoHits");
    adcToEnergyAlg->initialize_mip("calibration/reffiles/mip.root", "gaus_sigma>20", 0);
    adcToEnergyAlg->initialize_ped("calibration/reffiles/ped.root", "", 0);
    adcToEnergyAlg->initialize_dac("calibration/reffiles/dac.root", "", 0);
    algs.emplace_back(std::move(adcToEnergyAlg));
    LOG_INFO("AdcToEnergyReadTTreeAlg created successfully.");
    
    // Initialize MuonKFAlg
    std::unique_ptr<MuonKFAlg> muonKFAlg = 
    std::make_unique<MuonKFAlg>("RecoHits", "MuonTrack");
    algs.emplace_back(std::move(muonKFAlg));
    LOG_INFO("MuonKFAlg created successfully.");

    // Initialize WriterRegistry and register structs
    WriterRegistry registry;
    registry.register_vector_struct<AHCALRawHit>();
    registry.register_vector_struct<AHCALRecoHit>();
    registry.register_struct<Track>();
    RootWriterAlg rootWriter("output.root",  registry);
    algs.emplace_back(std::make_unique<RootWriterAlg>("output.root", registry));
    LOG_INFO("RootWriterAlg created successfully.");
    // event loop
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
        // LOG_INFO("--------------------------------------------------");
    }
    LOG_INFO("AHCAL Application finished.");
    return 0;
}