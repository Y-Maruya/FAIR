#pragma once
#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include "CalibDBIO/MIPReader/MIPReader.hpp"
#include "CalibDBIO/PedestalReader/PedestalReader.hpp"
#include "CalibDBIO/InterCalibReader/InterCalibReader.hpp"
#include <memory>
#include <string>
#include <utility>

namespace AHCALRecoAlg {
    struct AdcToEnergyAlgCfg{
        std::string in_rawhit_key = "RawHits";
        std::string out_recohit_key = "RecoHits";
        int mask_level = 2; // 0: no mask -> 1: mask pedestal lost channel -> 2: mask pedestal lost channel and intercalib lost channel
    };
    class AdcToEnergyAlg final : public IAlg {
    public:
        AdcToEnergyAlg(RunContext& ctx, std::string name)
            : IAlg(ctx, std::move(name)) {}
        ~AdcToEnergyAlg();
        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& n);
        void initialize() override;
        void init_by_run() override;
    private:
        void loadCalibrations();
        void initialize_mip(int runNumber);
        void initialize_pedestal(int runNumber);
        void initialize_intercalib(int runNumber);

        std::string m_in_rawhit_key;
        std::string m_out_recohit_key;
        std::shared_ptr<const CalibDBIO::PedestalMap> ped_map_ = std::make_shared<CalibDBIO::PedestalMap>();
        std::shared_ptr<const CalibDBIO::MIPMap> mip_map_ = std::make_shared<CalibDBIO::MIPMap>();
        std::shared_ptr<const CalibDBIO::HGLGRatioMap> hglg_map_ = std::make_shared<CalibDBIO::HGLGRatioMap>();
        std::shared_ptr<const CalibDBIO::HGSaturationMap> hg_saturation_map_ = std::make_shared<CalibDBIO::HGSaturationMap>();

        int loaded_run_ = -1;
        AdcToEnergyAlgCfg m_cfg;
    };
}
