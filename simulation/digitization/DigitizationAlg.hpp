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
        bool use_mip_threshold = false;
        bool implement_hardware_threshold_effect = true;
        double nmip_threshold = 0.5;
        bool add_noise = false;
        std::string noise_file = "noise.root";

        bool use_root_calibfile = true;
        std::string root_pedestal_file = "/afs/cern.ch/work/y/ymaruya/private/FASERlink/AHCAL/FAIR/out_calibration/pedestal_calib/22338/pedestal.root";
        // std::string pedestal_cut_string = "(fitStatus_hg==0 || fitStatus_hg==1) && (fitStatus_lg==0 || fitStatus_lg==1)";
        int pedestal_cellid_version = 1;

        bool use_dac_file = false;
        std::string root_dac_file = "dac_v2.root";
        // std::string dac_cut_string = "";
        int dac_cellid_version = 1;

        std::string root_intercalib_file = "intercalib.root";


        std::string root_mip_file = "mip.root";
        std::string mip_cut_string = "";
        int mip_cellid_version = 1;

        std::string root_spe_file = "spe.root";
        std::string spe_cut_string = "chi2/NDF < 20 && NDF >3 && spe > 15";
        int spe_cellid_version = 1;

        std::string sipm_model_file = "sipm_model_0.0xt.root";

        int seed = 42;
    };
    enum SimStatus {
        SIM_STATUS_NORMAL = 0,
        SIM_STATUS_NOISE = 1 << 0,
        SIM_STATUS_HG_PEDESTAL_FAILURE = 1 << 1,
        SIM_STATUS_LG_PEDESTAL_FAILURE = 1 << 2,
        SIM_STATUS_GAIN_RATIO_FAILURE = 1 << 3,
        SIM_STATUS_GAIN_RATIO_CORRECTED = 1 << 4,
        SIM_STATUS_MIP_FAILURE = 1 << 5,
        SIM_STATUS_MIP_IMPUTED_USED = 1 << 6,
        SIM_STATUS_MIP_FULLFIT_USED = 1 << 7,
        SIM_STATUS_MIP_REFERENCE_USED = 1 << 8,
        SIM_STATUS_SPE_FAILURE = 1 << 9,
        SIM_STATUS_SPE_IMPUTED_USED = 1 << 10,
        SIM_STATUS_SPE_FULLFIT_USED = 1 << 11,
        SIM_STATUS_SPE_REFERENCE_USED = 1 << 12,
        SIM_STATUS_THRESHOLD_SHOULD_ON = 1 << 13,
        SIM_STATUS_THRESHOLD_FAILURE = 1 << 14,
        SIM_STATUS_THRESHOLD_IMPUTED_USED = 1 << 15,
        SIM_STATUS_THRESHOLD_FULLFIT_USED = 1 << 16,
        SIM_STATUS_THRESHOLD_REFERENCE_USED = 1 << 17, // not used yet
        SIM_STATUS_NOCALIBDATA = 1 << 18
    };
    struct CalibrationParameters {
        double ped_high = 377.9; // from avg. run22074 
        double ped_high_error = 26.22; //standard deviation of the pedestal distribution
        double ped_low = 373.9; // from avg. run22074
        double ped_low_error = 25.96;
        double rms_high = 3.624;
        double rms_high_error = 1.375;
        double rms_low = 2.27;
        double rms_low_error = 0.1835;
        double gain_ratio = 29.9; // from avg. run21987, 
        double gain_ratio_error = 0.8793;
        double intercept = 113.5;
        double intercept_error = 75.69;
        double gain_plat = 3302.0;
        double lowgain_plat = 2000; // 
        double MIP = 326.3; // from avg. run21987 decision==1
        double MIP_error = 62.83;
        double Threshold = 0; // Not Used as the default
        double Threshold_error = 0;
        double ThresholdWidth = 0;
        double ThresholdWidth_error = 0;
        double SPE = 26.54; // from avg. run21987 snr>5 and gain<38
        double SPE_error = 4.021;
        int decision = 0; // decision on the MIP quality
        int sim_status = SIM_STATUS_NORMAL;
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
        // std::unordered_map<int, SimStatus> sim_status_map_;
        void set_default_calibration_parameters();
        double MIP_E = 0.461; // MeV
        int pixel = 7284;
        double pde = 0.32;
        void load_calibration_parameters_fromRoot();
        void load_calibration_parameters_fromDB();
        AHCALRawHit simulate_adc(const AHCALSimHit& simhit, const CalibrationParameters& calib);
        CalibrationParameters default_calib_ = {
            374.0,  // ped_high
            26.22,  // ped_high_error
            371.0,  // ped_low
            25.96,  // ped_low_error
            3.5,    // rms_high
            1.375,  // rms_high_error
            2.3,    // rms_low
            0.1835, // rms_low_error
            26.0,   // gain_ratio
            0.8793, // gain_ratio_error
            0.0,    // intercept
            0.0,    // intercept_error
            2927.0, // gain_plat
            2000.0, // lowgain_plat
            344.3,  // MIP
            62.83,  // MIP_error
            0.0,    // Threshold
            0.0,    // Threshold_error
            0.0,    // ThresholdWidth
            0.0,    // ThresholdWidth_error
            27.7,   // SPE
            4.021,  // SPE_error
            0,      // decision
            SIM_STATUS_NORMAL
        };
        std::unique_ptr<TH1I> h_sipm[15000];
        std::unique_ptr<TRandom3> rand_gen = std::make_unique<TRandom3>(cfg_.seed);
        std::unordered_map<int, int> nodata_cellid_map_;
    };
}
