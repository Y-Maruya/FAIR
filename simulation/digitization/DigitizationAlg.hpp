#pragma once
#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include "common/edm/EDM.hpp"
#include <yaml-cpp/yaml.h>
#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <stdexcept>
#include <TRandom3.h>
#include <TH1I.h>

namespace AHCALRecoAlg {
    struct DigitizationAlgCfg{
        std::string in_simhit_key = "SimHits";
        std::string out_rawhit_key = "SimRawHits";
        bool use_mip_threshold = true;
        double nmip_threshold = 0.5;
        bool add_noise = true;
        std::string noise_file = "noise.root";

        bool use_root_calibfile = true;
        std::string root_pedestal_file = "pedestal.root";
        // std::string pedestal_cut_string = "";
        int pedestal_cellid_version = 1;

        std::string root_dac_file = "dac_v2.root";
        // std::string dac_cut_string = "";
        int dac_cellid_version = 1;

        std::string root_mip_file = "mip.root";
        std::string mip_cut_string = "";
        int mip_cellid_version = 1;

        std::string root_spe_file = "spe.root";
        std::string spe_cut_string = "chi2/NDF < 20 && NDF >3 && spe > 15";
        int spe_cellid_version = 1;

        std::string sipm_model_file = "sipm_model_0.0xt.root";

        int seed = 42;
    };
    struct CalibrationParameters {
        double ped_high;
        double ped_low;
        double rms_high;
        double rms_low;
        double gain_ratio;
        double gain_plat;
        double lowgain_plat;
        double MIP;
        double SPE;
        double chi2_ndf;
        double NDF;
        double sigma;
    };
    struct ADCs {
        int hg_adc;
        int lg_adc;
    };
    class DigitizationAlg final : public IAlg{
    public:
        DigitizationAlg(RunContext& rc, std::string name)
        : IAlg(rc, std::move(name)){}
        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& cfg) override;
        void initialize() override;
    
    private:
        DigitizationAlgCfg cfg_;
        std::unordered_map<int, CalibrationParameters> calib_params_;
        void set_default_calibration_parameters();
        double MIP_E = 0.461; // MeV
        int pixel = 7284;
        double pde = 0.32;
        void load_calibration_parameters_fromRoot();
        void load_calibration_parameters_fromDB();
        AHCALRawHit simulate_adc(const AHCALSimHit& simhit, const CalibrationParameters& calib);
        CalibrationParameters default_calib_ = {
            374.0,  // ped_high
            371.0,  // ped_low
            3.5,    // rms_high
            2.3,    // rms_low
            26.0,   // gain_ratio
            2927.0, // gain_plat
            2000.0, // lowgain_plat
            344.3,  // MIP
            27.7,   // SPE
            0.0,    // chi2_ndf
            0.0,    // NDF
            7.0     // sigma
        };
        std::unique_ptr<TH1I> h_sipm[15000];
        std::unique_ptr<TRandom3> rand_gen = std::make_unique<TRandom3>(cfg_.seed);
    };
}