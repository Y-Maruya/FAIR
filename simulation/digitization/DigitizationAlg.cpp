#include "DigitizationAlg.hpp"
#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"
#include "calibration/RefValues.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <H5Cpp.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TRandom.h>
#include <TEventList.h>
#include <TBranch.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <string>
#include <TF1.h>
#include <TH1.h>
#include <TStyle.h>
#include "TMath.h"
#include "TRandom3.h"
#include <TSpectrum.h>
#include "TLegend.h"
#include "TLine.h"
#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGaxis.h"
#include <ROOT/RDataFrame.hxx>

AHCAL_REGISTER_ALG(AHCALRecoAlg::DigitizationAlg, "DigitizationAlg")


namespace AHCALRecoAlg{

    void DigitizationAlg::set_default_calibration_parameters(){
        for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer) {
            for (int chip = 0; chip < AHCALGeometry::chip_No; ++chip) {
                for (int channel = 0; channel < AHCALGeometry::channel_No; ++channel) {
                    int cellid = layer*100000 + chip*10000 + channel;
                    calib_params_.emplace(cellid, default_calib_);
                    nodata_cellid_map_.emplace(cellid, 0);
                }
            }
        }
    }
    void DigitizationAlg::initialize(){
        set_default_calibration_parameters();
        if (cfg_.use_root_calibfile) {
            load_calibration_parameters_fromRoot();
        } else {
            load_calibration_parameters_fromDB();
        }
    }
    void DigitizationAlg::execute(EventStore& evt){
        auto simHits = evt.get<std::vector<AHCALSimHit>>(cfg_.in_simhit_key);
        std::vector<AHCALRawHit> rawHits;
        for (const auto& simHit : simHits) {
            int cellid = simHit.cellID;
            auto it = calib_params_.find(cellid);
            if (it == calib_params_.end()) {
                LOG_WARN("DigitizationAlg: no calibration parameters found for cellID {}, using default values", cellid);
                continue;
            }
            const auto& calib = it->second;
            AHCALRawHit rawHit = simulate_adc(simHit, calib);
            if (rawHit.sim_status & SIM_STATUS_NOCALIBDATA) {
                LOG_DEBUG("DigitizationAlg: no calibration data for cellID {}, skipping hit", cellid);
                continue;
            }
            rawHits.push_back(rawHit);
        }
        evt.put(cfg_.out_rawhit_key, std::move(rawHits));
    }
    void DigitizationAlg::parse_cfg(const YAML::Node& cfg){
        cfg_.in_simhit_key = get_or<std::string>(cfg, "in_simhit_key", cfg_.in_simhit_key);
        cfg_.out_rawhit_key = get_or<std::string>(cfg, "out_rawhit_key", cfg_.out_rawhit_key);
        cfg_.use_mip_threshold = get_or<bool>(cfg, "use_mip_threshold", cfg_.use_mip_threshold);
        cfg_.nmip_threshold = get_or<double>(cfg, "nmip_threshold", cfg_.nmip_threshold);
        cfg_.add_noise = get_or<bool>(cfg, "add_noise", cfg_.add_noise);
        cfg_.noise_file = get_or<std::string>(cfg, "noise_file", cfg_.noise_file);
        cfg_.use_root_calibfile = get_or<bool>(cfg, "use_root_calibfile", cfg_.use_root_calibfile);
        cfg_.implement_hardware_threshold_effect = get_or<bool>(cfg, "implement_hardware_threshold_effect", cfg_.implement_hardware_threshold_effect);
        cfg_.root_pedestal_file = get_or<std::string>(cfg, "root_pedestal_file", cfg_.root_pedestal_file);
        if (cfg_.use_root_calibfile) {
            cfg_.pedestal_cellid_version = get_or<int>(cfg, "pedestal_cellid_version", cfg_.pedestal_cellid_version);
            cfg_.root_dac_file = get_or<std::string>(cfg, "root_dac_file", cfg_.root_dac_file);
            cfg_.use_dac_file = get_or<bool>(cfg, "use_dac_file", cfg_.use_dac_file);
            cfg_.root_intercalib_file = get_or<std::string>(cfg, "root_intercalib_file", cfg_.root_intercalib_file);
            cfg_.dac_cellid_version = get_or<int>(cfg, "dac_cellid_version", cfg_.dac_cellid_version);
            cfg_.root_mip_file = get_or<std::string>(cfg, "root_mip_file", cfg_.root_mip_file);
            cfg_.mip_cut_string = get_or<std::string>(cfg, "mip_cut_string", cfg_.mip_cut_string);
            cfg_.mip_cellid_version = get_or<int>(cfg, "mip_cellid_version", cfg_.mip_cellid_version);
            cfg_.root_spe_file = get_or<std::string>(cfg, "root_spe_file", cfg_.root_spe_file);
            cfg_.spe_cut_string = get_or<std::string>(cfg, "spe_cut_string", cfg_.spe_cut_string);
            cfg_.spe_cellid_version = get_or<int>(cfg, "spe_cellid_version", cfg_.spe_cellid_version);
        }
        cfg_.sipm_model_file = get_or<std::string>(cfg, "sipm_model_file", cfg_.sipm_model_file);
        cfg_.seed = get_or<int>(cfg, "seed", cfg_.seed);
    }
    void DigitizationAlg::load_calibration_parameters_fromRoot(){
        std::unique_ptr<TFile> fin(TFile::Open(cfg_.root_pedestal_file.c_str(), "READ"));
        if (!fin || fin->IsZombie()) {
            LOG_ERROR("DigitizationAlg: cannot open pedestal file: {}", cfg_.root_pedestal_file);
            throw std::runtime_error("Failed to open pedestal ROOT file");
            return;
        }
        TTree* tree_in = dynamic_cast<TTree*>(fin->Get("pedestal"));
        if (!tree_in) {
            LOG_ERROR("DigitizationAlg: cannot get pedestal tree from file: {}", cfg_.root_pedestal_file);
            throw std::runtime_error("Failed to load pedestal tree from ROOT file");
            return;
        }
        int cellid;
        int fitOk_hg = 1, fitOk_lg = 1;
        double pedestal_high, pedestal_low, sigma_high, sigma_low;
        double pedestal_high_error, pedestal_low_error, sigma_high_error, sigma_low_error;
        int fitStatus_hg = 0, fitStatus_lg = 0;
        if (tree_in->GetBranch("cellid")) tree_in->SetBranchAddress("cellid", &cellid);
        if (tree_in->GetBranch("CellID")) tree_in->SetBranchAddress("CellID", &cellid);
        tree_in->SetBranchAddress("highgain_peak", &pedestal_high);
        tree_in->SetBranchAddress("lowgain_peak", &pedestal_low);
        if (tree_in->GetBranch("highgain_peak_error")) tree_in->SetBranchAddress("highgain_peak_error", &pedestal_high_error);
        if (tree_in->GetBranch("lowgain_peak_error")) tree_in->SetBranchAddress("lowgain_peak_error", &pedestal_low_error);
        if (tree_in->GetBranch("highgain_sigma_error")) tree_in->SetBranchAddress("highgain_sigma_error", &sigma_high_error);
        if (tree_in->GetBranch("lowgain_sigma_error")) tree_in->SetBranchAddress("lowgain_sigma_error", &sigma_low_error);
        tree_in->SetBranchAddress("highgain_sigma", &sigma_high);
        tree_in->SetBranchAddress("lowgain_sigma", &sigma_low);
        if (tree_in->GetBranch("fitOk_hg")) tree_in->SetBranchAddress("fitOk_hg", &fitOk_hg);
        if (tree_in->GetBranch("fitOk_lg")) tree_in->SetBranchAddress("fitOk_lg", &fitOk_lg);
        if (tree_in->GetBranch("fitStatus_hg")) tree_in->SetBranchAddress("fitStatus_hg", &fitStatus_hg);
        if (tree_in->GetBranch("fitStatus_lg")) tree_in->SetBranchAddress("fitStatus_lg", &fitStatus_lg);
        for (int i = 0; i < tree_in->GetEntries(); ++i) {
            tree_in->GetEntry(i);
            cellid = AHCALGeometry::cellid_conversion(cellid, cfg_.pedestal_cellid_version);
            nodata_cellid_map_[cellid] = 1;
            if (fitOk_hg == 1 && AHCALRefValues::HGPedestalStatus_is_ok(fitStatus_hg)) {
                if (pedestal_high<0 || sigma_high<=0) {
                    LOG_DEBUG("Invalid pedestal fit result for HG: cellID={} pedestal={} sigma={}", cellid, pedestal_high, sigma_high);
                    continue;
                }
                calib_params_[cellid].ped_high = pedestal_high;
                calib_params_[cellid].rms_high = sigma_high;
                calib_params_[cellid].ped_high_error = pedestal_high_error;
                calib_params_[cellid].rms_high_error = sigma_high_error;
            }else{
                calib_params_[cellid].sim_status |= SIM_STATUS_HG_PEDESTAL_FAILURE;
            }
            if (fitOk_lg == 1 && AHCALRefValues::LGPedestalStatus_is_ok(fitStatus_lg)) {
                if ( pedestal_low < 0 || sigma_low <= 0) {
                    LOG_DEBUG("Invalid pedestal fit result for LG: cellID={} pedestal={} sigma={}", cellid, pedestal_low, sigma_low);
                    continue;
                }
                calib_params_[cellid].ped_low = pedestal_low;
                calib_params_[cellid].rms_low = sigma_low;
                calib_params_[cellid].ped_low_error = pedestal_low_error;
                calib_params_[cellid].rms_low_error = sigma_low_error;
            }else{
                calib_params_[cellid].sim_status |= SIM_STATUS_LG_PEDESTAL_FAILURE;
            }
            if (fitOk_hg != 1 && fitOk_lg != 1) {
                LOG_DEBUG("Both HG and LG pedestal fits failed for cellID={}", cellid);
                calib_params_[cellid].sim_status |= SIM_STATUS_NOCALIBDATA;
            }
        }
        for (const auto& [cellid, has_data] : nodata_cellid_map_) {
            if (has_data == 0) {
                LOG_DEBUG("No pedestal data for cellID={}", cellid);
                calib_params_[cellid].sim_status |= SIM_STATUS_NOCALIBDATA;
            }
        }
        tree_in = nullptr;
        fin->Close();
        nodata_cellid_map_.clear();
        // read dac
        if (cfg_.use_dac_file) {
            fin.reset(TFile::Open(cfg_.root_dac_file.c_str(), "READ"));
            if (!fin || fin->IsZombie()) {  
                LOG_ERROR("DigitizationAlg: cannot open dac file: {}", cfg_.root_dac_file);
                return;
            }
            tree_in = dynamic_cast<TTree*>(fin->Get("dac"));
            if (!tree_in) {
                LOG_ERROR("DigitizationAlg: cannot get dac tree from file: {}", cfg_.root_dac_file);
                return;
            }
            double gain_ratio, gain_plat, lowgain_plat;
            if (tree_in->GetBranch("cellid")) tree_in->SetBranchAddress("cellid", &cellid);
            if (tree_in->GetBranch("CellID")) tree_in->SetBranchAddress("CellID", &cellid);
            tree_in->SetBranchAddress("slope", &gain_ratio);
            tree_in->SetBranchAddress("plat", &gain_plat);
            tree_in->SetBranchAddress("lowgain_satu_point", &lowgain_plat);
            for (int i = 0; i < tree_in->GetEntries(); ++i) {
                tree_in->GetEntry(i);
                cellid = AHCALGeometry::cellid_conversion(cellid, cfg_.dac_cellid_version);
                nodata_cellid_map_[cellid] = 1;
                if (gain_plat < 0 || lowgain_plat < 0 || gain_ratio < 0) {
                    LOG_DEBUG("Invalid DAC fit result: cellID={} gain_ratio={} gain_plat={} lowgain_plat={}", cellid, gain_ratio, gain_plat, lowgain_plat);
                    continue;
                }
                calib_params_[cellid].gain_ratio = gain_ratio;  
                calib_params_[cellid].gain_plat = gain_plat;
                calib_params_[cellid].lowgain_plat = lowgain_plat;
                if (lowgain_plat < 2000) {
                    LOG_DEBUG("Low gain saturation point detected: cellID={} lowgain_plat={}", cellid, lowgain_plat);
                    calib_params_[cellid].lowgain_plat = 2000; 
                }
            }
            for (const auto& [cellid, has_data] : nodata_cellid_map_) {
                if (has_data == 0) {
                    LOG_DEBUG("No DAC data for cellID={}", cellid);
                    calib_params_[cellid].sim_status |= SIM_STATUS_NOCALIBDATA;
                }
            }
            tree_in = nullptr;
            fin->Close();
            nodata_cellid_map_.clear();
        }else{
            fin.reset(TFile::Open(cfg_.root_intercalib_file.c_str(), "READ"));
            if (!fin || fin->IsZombie()) {
                LOG_ERROR("DigitizationAlg: cannot open intercalib file: {}", cfg_.root_intercalib_file);
                return;
            }
            tree_in = dynamic_cast<TTree*>(fin->Get("intercalib"));
            if (!tree_in) {
                LOG_ERROR("DigitizationAlg: cannot get intercalib tree from file: {}", cfg_.root_intercalib_file);
                return;
            }
            double slope, intercept, gain_plat;
            double slope_error, intercept_error;
            int QualityFlag;

            if (tree_in->GetBranch("cellid")) tree_in->SetBranchAddress("cellid", &cellid);
            if (tree_in->GetBranch("CellID")) tree_in->SetBranchAddress("CellID", &cellid);
            tree_in->SetBranchAddress("Slope", &slope);
            tree_in->SetBranchAddress("Intercept", &intercept);
            tree_in->SetBranchAddress("SlopeError", &slope_error);
            tree_in->SetBranchAddress("InterceptError", &intercept_error);
            tree_in->SetBranchAddress("QualityFlag", &QualityFlag);
            tree_in->SetBranchAddress("HG_SaturationPoint", &gain_plat);

            for (int i = 0; i < tree_in->GetEntries(); ++i) {
                tree_in->GetEntry(i);
                // cellid = AHCALGeometry::cellid_conversion(cellid, cfg_.dac_cellid_version);
                if (slope < 0 || gain_plat < 0) {
                    LOG_DEBUG("Invalid intercalib fit result: cellID={} slope={} gain_plat={}", cellid, slope, gain_plat);
                    continue;
                }
                nodata_cellid_map_[cellid] = 1;
                if (QualityFlag == 0) {
                    calib_params_[cellid].gain_ratio = slope;  
                    calib_params_[cellid].gain_plat = gain_plat;
                    calib_params_[cellid].intercept = intercept;
                    calib_params_[cellid].gain_ratio_error = slope_error;
                    calib_params_[cellid].intercept_error = intercept_error;
                }else if (QualityFlag == 1) {
                    LOG_DEBUG("Intercalib fit result flagged as masked: cellID={} slope={} gain_plat={}", cellid, slope, gain_plat);
                    calib_params_[cellid].sim_status |= SIM_STATUS_GAIN_RATIO_FAILURE;
                }else if (QualityFlag == 2) {
                    LOG_DEBUG("Intercalib fit result flagged as corrected: cellID={} slope={} gain_plat={}", cellid, slope, gain_plat);
                    calib_params_[cellid].sim_status |= SIM_STATUS_GAIN_RATIO_CORRECTED;
                    calib_params_[cellid].gain_ratio = slope;  
                    calib_params_[cellid].gain_plat = gain_plat;
                    calib_params_[cellid].intercept = intercept;
                    calib_params_[cellid].gain_ratio_error = slope_error;
                    calib_params_[cellid].intercept_error = intercept_error;
                }
            }
            for (const auto& [cellid, has_data] : nodata_cellid_map_) {
                if (has_data == 0) {
                    LOG_DEBUG("No intercalib data for cellID={}", cellid);
                    calib_params_[cellid].sim_status |= SIM_STATUS_NOCALIBDATA;
                }
            }
            tree_in = nullptr;
            fin->Close();
            nodata_cellid_map_.clear();
        }


        // read MIP
        fin.reset(TFile::Open(cfg_.root_mip_file.c_str(), "READ"));
        if (!fin || fin->IsZombie()) {
            LOG_ERROR("DigitizationAlg: cannot open MIP file: {}", cfg_.root_mip_file);
            return;
        }
        tree_in = dynamic_cast<TTree*>(fin->Get("mip"));
        if (!tree_in) {
            LOG_ERROR("DigitizationAlg: cannot get MIP tree from file: {}", cfg_.root_mip_file);
            return;
        }
        ROOT::RDataFrame df(*tree_in);
        ROOT::RDF::RNode df_cut = df;
        if (!cfg_.mip_cut_string.empty()) df_cut = df.Filter(cfg_.mip_cut_string);

        if (!cfg_.mip_cut_string.empty()) {
            if (df_cut.Count().GetValue() == 0) {
            LOG_ERROR("No entries left after applying cut: {}", cfg_.mip_cut_string);
            return;
            }
            LOG_INFO("Applied cut string: {}", cfg_.mip_cut_string);
        }

        const auto cellids = df_cut.Take<int>("cellid").GetValue();
        const auto mpvs    = df_cut.Take<double>("MPV").GetValue();
        const auto mpv_errors = df_cut.Take<double>("MPVError").GetValue();
        const auto Thresholds = df_cut.Take<double>("Threshold").GetValue();
        const auto Threshold_errors = df_cut.Take<double>("ThresholdError").GetValue();
        const auto ThresholdWidths = df_cut.Take<double>("ThresholdWidth").GetValue();
        const auto ThresholdWidth_errors = df_cut.Take<double>("ThresholdWidthError").GetValue();
        const auto decisions = df_cut.Take<int>("Decision").GetValue();
        const auto states = df_cut.Take<int>("State").GetValue();
        const auto fullfits = df_cut.Take<int>("FullFit").GetValue();
        const auto imputeds = df_cut.Take<int>("Imputed").GetValue();
        const auto refs = df_cut.Take<int>("Ref").GetValue();
        if (cellids.size() != mpvs.size() || cellids.size() != mpv_errors.size() || cellids.size() != Thresholds.size() || cellids.size() != Threshold_errors.size()) {
            LOG_ERROR("MIP tree size mismatch: cellid={} MPV={} MPV_error={} Threshold={} Threshold_error={}", cellids.size(), mpvs.size(), mpv_errors.size(), Thresholds.size(), Threshold_errors.size());
            return;
        }
        int num_mip_passed = 0;
        for (size_t i = 0; i < cellids.size(); ++i) {
            const int conv_cellid = AHCALGeometry::cellid_conversion(cellids[i], cfg_.mip_cellid_version);
            const int mip_fullfit = fullfits[i] & 1;
            const int mip_imputed = imputeds[i] & 1;
            const int mip_ref = refs[i] & 1;
            const int threshold_fullfit = fullfits[i] & 8;
            const int threshold_imputed = imputeds[i] & 8;
            const int threshold_ref = refs[i] & 8;
            calib_params_[conv_cellid].decision = decisions[i];
            if (mip_fullfit == 1 || mip_imputed == 1 || mip_ref == 1) {
                calib_params_[conv_cellid].sim_status |= SIM_STATUS_MIP_FAILURE;
            }
            if (mip_fullfit == 1){
                calib_params_[conv_cellid].sim_status |= SIM_STATUS_MIP_FULLFIT_USED;
            }
            if (mip_imputed == 1) {
                calib_params_[conv_cellid].sim_status |= SIM_STATUS_MIP_IMPUTED_USED;
            }
            if (mip_ref == 1) {
                calib_params_[conv_cellid].sim_status |= SIM_STATUS_MIP_REFERENCE_USED;
            }
            if (threshold_fullfit == 8 || threshold_imputed == 8 || threshold_ref == 8) {
                calib_params_[conv_cellid].sim_status |= SIM_STATUS_THRESHOLD_FAILURE;
            }
            if (threshold_fullfit == 8) {
                calib_params_[conv_cellid].sim_status |= SIM_STATUS_THRESHOLD_FULLFIT_USED;
            }
            if (threshold_imputed == 8) {
                calib_params_[conv_cellid].sim_status |= SIM_STATUS_THRESHOLD_IMPUTED_USED;
            }
            if (threshold_ref == 8) {
                calib_params_[conv_cellid].sim_status |= SIM_STATUS_THRESHOLD_REFERENCE_USED;
            }
            if (mpvs[i] <= 0 || mpv_errors[i] <= 0 ) {
                LOG_DEBUG("Invalid MIP detected: cellID={} MPV={} MPV_error={} Threshold={} Threshold_error={}", conv_cellid, mpvs[i], mpv_errors[i], Thresholds[i], Threshold_errors[i]);
                // calib_params_[conv_cellid].sim_status |= SIM_STATUS_MIP_FAILURE;
                calib_params_[conv_cellid].sim_status |= SIM_STATUS_NOCALIBDATA;
                // calib_params_[conv_cellid].sim_status |= SIM_STATUS_THRESHOLD_FAILURE;
                continue;
            }
            nodata_cellid_map_[conv_cellid] = 1;
            calib_params_[conv_cellid].MIP = mpvs[i];
            calib_params_[conv_cellid].MIP_error = mpv_errors[i];
            num_mip_passed++;
            bool threshold_good_nonneed = (decisions[i] == 2 && (states[i] & (1<<25) || states[i] & (1<<26))); //if efficiency is higher than 0.98, the threshold is considered good, but not needed for simulation
            if (decisions[i] == 2 && !threshold_good_nonneed) {
                calib_params_[conv_cellid].sim_status |= SIM_STATUS_THRESHOLD_SHOULD_ON;
            }
            if (decisions[i] == 6 ){
                calib_params_[conv_cellid].sim_status |= SIM_STATUS_THRESHOLD_SHOULD_ON;
                if (threshold_fullfit == 8 || threshold_imputed == 8 || threshold_ref == 8) {
                    calib_params_[conv_cellid].sim_status |= SIM_STATUS_THRESHOLD_FAILURE;
                }else{
                    // calib_params_[conv_cellid].sim_status |= SIM_STATUS_THRESHOLD_SHOULD_ON;
                    LOG_WARN("Threshold should be on for cellID={} but no threshold data available", conv_cellid);
                }
            }
            if ((decisions[i] == 2 && !threshold_good_nonneed) || threshold_fullfit == 8 || threshold_imputed == 8 || threshold_ref == 8) {
                if (Thresholds[i] <= 0 || Threshold_errors[i] <= 0 || ThresholdWidths[i] <= 0 || ThresholdWidth_errors[i] <= 0) {
                    LOG_WARN("Invalid threshold detected: cellID={} Threshold={} Threshold_error={} ThresholdWidth={} ThresholdWidth_error={}", conv_cellid, Thresholds[i], Threshold_errors[i], ThresholdWidths[i], ThresholdWidth_errors[i]);
                    LOG_WARN("decision={} threshold_fullfit={} threshold_imputed={} threshold_ref={}", decisions[i], threshold_fullfit, threshold_imputed, threshold_ref);
                    calib_params_[conv_cellid].sim_status |= SIM_STATUS_THRESHOLD_FAILURE;
                    continue;
                }
                calib_params_[conv_cellid].Threshold = Thresholds[i];
                calib_params_[conv_cellid].Threshold_error = Threshold_errors[i];
                calib_params_[conv_cellid].ThresholdWidth = ThresholdWidths[i];
                calib_params_[conv_cellid].ThresholdWidth_error = ThresholdWidth_errors[i];
            }
        }
        for (const auto& [cellid, has_data] : nodata_cellid_map_) {
            if (has_data == 0) {
                LOG_DEBUG("No MIP data for cellID={}", cellid);
                calib_params_[cellid].sim_status |= SIM_STATUS_NOCALIBDATA;
            }
        }
        LOG_INFO("Loaded {} MIP entries from {}", num_mip_passed, cfg_.root_mip_file);
        tree_in = nullptr;
        fin->Close();
        nodata_cellid_map_.clear();
        // read spe
        fin.reset(TFile::Open(cfg_.root_spe_file.c_str(), "READ"));
        if (!fin || fin->IsZombie()) {
            LOG_ERROR("DigitizationAlg: cannot open SPE file: {}", cfg_.root_spe_file);
            return;
        }
        tree_in = dynamic_cast<TTree*>(fin->Get("spe"));
        if (!tree_in) {
            LOG_ERROR("DigitizationAlg: cannot get SPE tree from file: {}", cfg_.root_spe_file);
            return;
        }
        ROOT::RDataFrame df_spe(*tree_in);
        ROOT::RDF::RNode df_spe_cut = df_spe;
        if (!cfg_.spe_cut_string.empty()) df_spe_cut = df_spe.Filter(cfg_.spe_cut_string);
        if (!cfg_.spe_cut_string.empty()) {
            if (df_spe_cut.Count().GetValue() == 0) {
            LOG_ERROR("No entries left after applying cut: {}", cfg_.spe_cut_string);
            return;
            }
            LOG_INFO("Applied cut string: {}", cfg_.spe_cut_string);
        }
        const auto spe_cellids = df_spe_cut.Take<int>("cellid").GetValue();
        const auto spes = df_spe_cut.Take<double>("gain").GetValue();
        const auto spe_errors = df_spe_cut.Take<double>("gainErr").GetValue();
        const auto snrs = df_spe_cut.Take<double>("snr").GetValue();
        if (spe_cellids.size() != spes.size()) {
            LOG_ERROR("SPE tree size mismatch: cellid={} spe={}", spe_cellids.size(), spes.size());
            return;
        }
        int num_spe_passed = 0;
        for (size_t i = 0; i < spe_cellids.size(); ++i) {
            const int conv_cellid = AHCALGeometry::cellid_conversion(spe_cellids[i], cfg_.spe_cellid_version);
            if (spes[i] <= 0 || spe_errors[i] <= 0) {
                LOG_DEBUG("Invalid SPE detected: cellID={} SPE={} SPE_error={}", conv_cellid, spes[i], spe_errors[i]);
                calib_params_[conv_cellid].sim_status |= SIM_STATUS_NOCALIBDATA;
                continue;
            }
            nodata_cellid_map_[conv_cellid] = 1;
            if (snrs[i] < 5 || spes[i] > 38) {// temporary
                LOG_DEBUG("Low SNR detected: cellID={} SNR={}", conv_cellid, snrs[i]);
                calib_params_[conv_cellid].sim_status |= SIM_STATUS_SPE_FAILURE;
                continue;
            }
            calib_params_[conv_cellid].SPE = spes[i];
            calib_params_[conv_cellid].SPE_error = spe_errors[i];
            num_spe_passed++;
        }
        for (const auto& [cellid, has_data] : nodata_cellid_map_) {
            if (has_data == 0) {
                LOG_DEBUG("No SPE data for cellID={}", cellid);
                calib_params_[cellid].sim_status |= SIM_STATUS_NOCALIBDATA;
            }
        }
        tree_in = nullptr;
        fin->Close();
        LOG_INFO("Loaded {} SPE entries from {}", num_spe_passed, cfg_.root_spe_file);
        // read sipm model
        fin.reset(TFile::Open(cfg_.sipm_model_file.c_str(), "READ"));
        if (!fin || fin->IsZombie()) {
            LOG_ERROR("DigitizationAlg: cannot open SiPM model file: {}", cfg_.sipm_model_file);
            return;
        }
        for (int i = 0; i < 15000; i++)
        {
            TString name = "hist/incident_" + TString(std::to_string(i + 1).c_str()) + "_photons";
            auto* hist = dynamic_cast<TH1I*>(fin->Get(name));
            if (!hist) {
                LOG_ERROR("DigitizationAlg: cannot get histogram {} from SiPM model file: {}", name.Data(), cfg_.sipm_model_file);
                throw std::runtime_error("Failed to load SiPM model histogram");
                return;
            }
            h_sipm[i] = std::make_unique<TH1I>(*hist);
            h_sipm[i]->SetDirectory(nullptr); // Detach histogram from file to keep it in memory after file is closed
        }
        fin->Close();
    }
    void DigitizationAlg::load_calibration_parameters_fromDB(){
        //tobe implemented
    }
    AHCALRawHit DigitizationAlg::simulate_adc(const AHCALSimHit& simhit, const CalibrationParameters& calib){
        AHCALRawHit rawHit;
        rawHit.index = simhit.index;
        rawHit.cellID = simhit.cellID;
        rawHit.hg_adc = 0;
        rawHit.lg_adc = 0;
        rawHit.sim_status = calib.sim_status;
        if (calib.sim_status & SIM_STATUS_NOCALIBDATA) {
            return rawHit;
        }
        double energy = rand_gen->Gaus(simhit.Edep, simhit.Edep * 0.06);
        energy = rand_gen->Gaus(energy, energy * 0.015);
        energy = std::max(0.0, energy);
        if (calib.MIP <= 0.0 || calib.SPE <= 0.0 || calib.gain_ratio <= 0.0) {
            LOG_WARN("Simulated hit for cellID {} has invalid MIP, SPE, or gain_ratio value, using default values", simhit.cellID);
        }
        double mip = calib.MIP > 0.0 ? calib.MIP : default_calib_.MIP;
        double mip_error = calib.MIP_error > 0.0 ? calib.MIP_error : default_calib_.MIP_error;
        double spe = calib.SPE > 0.0 ? calib.SPE : default_calib_.SPE;
        double spe_error = calib.SPE_error > 0.0 ? calib.SPE_error : default_calib_.SPE_error;
        if (calib.sim_status & SIM_STATUS_MIP_FAILURE) {
            LOG_DEBUG("Simulated hit for cellID {} has MIP failure status: MIP={}", simhit.cellID, calib.MIP);
            // mip = default_calib_.MIP;
            // mip_error = default_calib_.MIP_error;
        }
        if (calib.sim_status & SIM_STATUS_SPE_FAILURE) {
            LOG_DEBUG("Simulated hit for cellID {} has SPE failure status: SPE={}", simhit.cellID, calib.SPE);
            spe = default_calib_.SPE;
            spe_error = default_calib_.SPE_error;
        }
        double gain_ratio = calib.gain_ratio > 0.0 ? calib.gain_ratio : default_calib_.gain_ratio;
        if (calib.sim_status & SIM_STATUS_GAIN_RATIO_FAILURE) {
            LOG_DEBUG("Simulated hit for cellID {} has gain ratio failure status: gain_ratio={}", simhit.cellID, calib.gain_ratio);
            gain_ratio = default_calib_.gain_ratio;
        }

        double n_photons = energy / MIP_E * mip / spe / pde;
        double n_photons_error = n_photons * std::sqrt(std::pow(mip_error/mip, 2) + std::pow(spe_error/spe, 2));

        n_photons = rand_gen->Poisson(n_photons);

        int n_photoelectrons = rand_gen->Binomial(static_cast<int>(n_photons), pde);
        if (n_photoelectrons >= 15000) {
            n_photoelectrons = 14999;
        }
        int n_fired;
        int n_fired_upper;
        int n_fired_lower;
        if (n_photoelectrons == 0) {
            n_fired = 0;
            n_fired_upper = 0;
            n_fired_lower = 0;
        } else if (!h_sipm[n_photoelectrons - 1]) {
            n_fired = n_photoelectrons;
            n_fired_upper = n_photoelectrons;
            n_fired_lower = n_photoelectrons;
            LOG_WARN("No SiPM model histogram for {} photoelectrons, using linear approximation", n_photoelectrons);
        } else {
            const int upper_photoelectrons =
                std::clamp(static_cast<int>(std::round(n_photoelectrons + n_photons_error)), 1, 14999);
            const int lower_photoelectrons =
                std::clamp(static_cast<int>(std::round(n_photoelectrons - n_photons_error)), 1, 14999);
            n_fired = int(h_sipm[n_photoelectrons - 1]->GetRandom());
            n_fired_upper = h_sipm[upper_photoelectrons - 1]
                ? int(h_sipm[upper_photoelectrons - 1]->GetRandom())
                : upper_photoelectrons;
            n_fired_lower = h_sipm[lower_photoelectrons - 1]
                ? int(h_sipm[lower_photoelectrons - 1]->GetRandom())
                : lower_photoelectrons;
        }
        double adc_value_upper = n_fired_upper * spe;
        double adc_value_lower = n_fired_lower * spe;

        double adc_value = n_fired * spe;
        double adc_value_error = std::max(std::abs(adc_value_upper - adc_value), std::abs(adc_value_lower - adc_value));
        LOG_DEBUG("Simulated hit for cellID {}: energy={} MeV, n_photons={} +/- {}, n_photoelectrons={} +/- {}, n_fired={} ({} - {}), adc_value={} +/- {}", simhit.cellID, energy, n_photons, n_photons_error, n_photoelectrons, std::sqrt(n_photoelectrons * pde * (1 - pde)), n_fired, n_fired_lower, n_fired_upper, adc_value, adc_value_error);
        double adc_noise = std::sqrt(n_fired) * 8.0; 

        adc_value = rand_gen->Gaus(adc_value, adc_noise);
        double ped_high = calib.ped_high;
        double rms_high = calib.rms_high;
        double ped_low = calib.ped_low;
        double rms_low = calib.rms_low;
        if (calib.ped_high < 0 || calib.rms_high <= 0 || calib.ped_low < 0 || calib.rms_low <= 0) {
            LOG_WARN("Simulated hit for cellID {} has invalid pedestal or RMS values, using default values", simhit.cellID);
            ped_high = default_calib_.ped_high;
            rms_high = default_calib_.rms_high;
            ped_low = default_calib_.ped_low;
            rms_low = default_calib_.rms_low;
        }
        // threshold effect need to implemented 
        double efficiency = 1.0;
        if (cfg_.implement_hardware_threshold_effect && (calib.sim_status & SIM_STATUS_THRESHOLD_SHOULD_ON) && (calib.Threshold > 0 && calib.ThresholdWidth > 0)) {
            double threshold = calib.Threshold;
            double threshold_width = calib.ThresholdWidth;
            if (threshold <= 0 || threshold_width <= 0) {
                LOG_WARN("Simulated hit for cellID {} has invalid threshold or threshold width values, skipping threshold effect", simhit.cellID);
            } else {
                efficiency = 0.5 * (1 + std::erf((adc_value - threshold) / (std::sqrt(2) * threshold_width)));
                LOG_DEBUG("Simulated hit for cellID {}: adc_value={} Threshold={} ThresholdWidth={} Efficiency={}", simhit.cellID, adc_value, threshold, threshold_width, efficiency);
            }
        }
        rawHit.hittag = 1;
        double HG = adc_value + rand_gen->Gaus(ped_high, 1.5 * rms_high);
        double LG = adc_value / gain_ratio + rand_gen->Gaus(ped_low, 5 * rms_low);
        if (HG > calib.gain_plat) {
            HG = calib.gain_plat;
        }
        if (LG > calib.lowgain_plat) {
            LG = calib.lowgain_plat;
        }
        rawHit.hg_adc = static_cast<int>(std::round(HG));
        rawHit.lg_adc = static_cast<int>(std::round(LG));
        if (rand_gen->Uniform() > efficiency) {
            LOG_DEBUG("Simulated hit for cellID {}: adc_value={} below threshold, hit rejected", simhit.cellID, adc_value);
            rawHit.hittag = 0;
            return rawHit;
        }
        return rawHit;
    }
}

        
