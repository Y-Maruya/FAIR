#include "SimMimicAlg.hpp"
#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"
#include "calibration/RefValues.hpp"

AHCAL_REGISTER_ALG(AHCALRecoAlg::SimMimicAlg, "SimMimicAlg")
namespace AHCALRecoAlg {
    void SimMimicAlg::execute(EventStore& evt) {
        auto simhits = evt.get<std::vector<AHCALSimHit>>(cfg_.in_simhit_key);
        std::vector<AHCALRecoHit> recohits;
        std::vector<AHCALRawHit> rawhits;
        for (const auto& simhit : simhits) {
            AHCALRecoHit recohit;
            AHCALRawHit rawhit;
            recohit.cellID = simhit.cellID;
            recohit.Edep = simhit.Edep;
            recohit.Nmip = simhit.Edep /AHCALRefValues::MIP_E; // example conversion
            recohit.index = simhit.index;

            rawhit.cellID = simhit.cellID;
            rawhit.hg_adc = static_cast<int>(simhit.Nmip*AHCALRefValues::ref_MIP)+400; // example conversion
            rawhit.hittag = 1;
            rawhit.bcid = 0;
            rawhit.index = simhit.index;

            recohits.push_back(recohit);
            rawhits.push_back(rawhit);
        }
        evt.put(cfg_.out_recohit_key, std::move(recohits));
        evt.put(cfg_.out_rawhit_key, std::move(rawhits));
    }
    void SimMimicAlg::initialize() {
        LOG_INFO("SimMimicAlg initialized with in_simhit_key={}, out_recohit_key={}, out_rawhit_key={}", 
            cfg_.in_simhit_key, cfg_.out_recohit_key, cfg_.out_rawhit_key);
    }
    void SimMimicAlg::finalize() {
        LOG_INFO("SimMimicAlg finalized.");
    }
    void SimMimicAlg::parse_cfg(const YAML::Node& cfg) {
        cfg_.in_simhit_key = get_or<std::string>(cfg, "in_simhit_key", cfg_.in_simhit_key);
        cfg_.out_recohit_key = get_or<std::string>(cfg, "out_recohit_key", cfg_.out_recohit_key);
        cfg_.out_rawhit_key = get_or<std::string>(cfg, "out_rawhit_key", cfg_.out_rawhit_key);
        LOG_INFO("SimMimicAlg configuration: in_simhit_key={}, out_recohit_key={}, out_rawhit_key={}", 
            cfg_.in_simhit_key, cfg_.out_recohit_key, cfg_.out_rawhit_key);
    }
}
