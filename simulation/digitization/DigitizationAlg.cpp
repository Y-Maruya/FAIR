#include "DigitizationAlg.hpp"
#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"

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
        cfg_.root_pedestal_file = get_or<std::string>(cfg, "root_pedestal_file", cfg_.root_pedestal_file);
        if (cfg_.use_root_calibfile) {
            cfg_.pedestal_cellid_version = get_or<int>(cfg, "pedestal_cellid_version", cfg_.pedestal_cellid_version);
            cfg_.root_dac_file = get_or<std::string>(cfg, "root_dac_file", cfg_.root_dac_file);
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
        std::unique_ptr<TFile> fin = std::make_unique<TFile>();
        fin->Open(cfg_.root_pedestal_file.c_str(), "READ");
        if (!fin || fin->IsZombie()) {
            LOG_ERROR("DigitizationAlg: cannot open pedestal file: {}", cfg_.root_pedestal_file);
            return;
        }
        std::unique_ptr<TTree> tree_in;
        tree_in = std::unique_ptr<TTree>((TTree*)fin->Get("pedestal"));
        if (!tree_in) {
            LOG_ERROR("DigitizationAlg: cannot get pedestal tree from file: {}", cfg_.root_pedestal_file);
            return;
        }
        int cellid;
        int fitOk_hg = 1, fitOk_lg = 1;
        double pedestal_high, pedestal_low, sigma_high, sigma_low;
        if (tree_in->GetBranch("cellid")) tree_in->SetBranchAddress("cellid", &cellid);
        if (tree_in->GetBranch("CellID")) tree_in->SetBranchAddress("CellID", &cellid);
        tree_in->SetBranchAddress("highgain_peak", &pedestal_high);
        tree_in->SetBranchAddress("lowgain_peak", &pedestal_low);
        tree_in->SetBranchAddress("highgain_sigma", &sigma_high);
        tree_in->SetBranchAddress("lowgain_sigma", &sigma_low);
        if (tree_in->GetBranch("fitOk_hg")) tree_in->SetBranchAddress("fitOk_hg", &fitOk_hg);
        if (tree_in->GetBranch("fitOk_lg")) tree_in->SetBranchAddress("fitOk_lg", &fitOk_lg);
        for (int i = 0; i < tree_in->GetEntries(); ++i) {
            tree_in->GetEntry(i);
            cellid = AHCALGeometry::cellid_conversion(cellid, cfg_.pedestal_cellid_version);
            if (fitOk_hg == 1){
                if (pedestal_high<0 || sigma_high<=0) {
                    LOG_DEBUG("Invalid pedestal fit result for HG: cellID={} pedestal={} sigma={}", cellid, pedestal_high, sigma_high);
                    continue;
                }
                calib_params_[cellid].ped_high = pedestal_high;
                calib_params_[cellid].rms_high = sigma_high;
            }
            if (fitOk_lg == 1){
                if ( pedestal_low < 0 || sigma_low <= 0) {
                    LOG_DEBUG("Invalid pedestal fit result for LG: cellID={} pedestal={} sigma={}", cellid, pedestal_low, sigma_low);
                    continue;
                }
                calib_params_[cellid].ped_low = pedestal_low;
                calib_params_[cellid].rms_low = sigma_low;
            }
        }
        tree_in->Delete();
        fin->Close();
        // read dac
        fin->Open(cfg_.root_dac_file.c_str(), "READ");
        if (!fin || fin->IsZombie()) {  
            LOG_ERROR("DigitizationAlg: cannot open dac file: {}", cfg_.root_dac_file);
            return;
        }
        tree_in = std::unique_ptr<TTree>((TTree*)fin->Get("dac"));
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
        tree_in->Delete();
        fin->Close();

        // read MIP
        fin->Open(cfg_.root_mip_file.c_str(), "READ");
        if (!fin || fin->IsZombie()) {
            LOG_ERROR("DigitizationAlg: cannot open MIP file: {}", cfg_.root_mip_file);
            return;
        }
        tree_in = std::unique_ptr<TTree>((TTree*)fin->Get("mip"));
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
        if (cellids.size() != mpvs.size()) {
            LOG_ERROR("MIP tree size mismatch: cellid={} MPV={}", cellids.size(),
                    mpvs.size());
            return;
        }
        int num_mip_passed = 0;
        for (size_t i = 0; i < cellids.size(); ++i) {
            const int conv_cellid = AHCALGeometry::cellid_conversion(cellids[i], cfg_.mip_cellid_version);
            calib_params_[conv_cellid].MIP = mpvs[i];
            if (mpvs[i] <= 100.0) {
                LOG_DEBUG("Low MIP value detected: cellID={} MIP={}", conv_cellid, mpvs[i]);
                calib_params_[conv_cellid].MIP = default_calib_.MIP;
                continue;
            }
            num_mip_passed++;
        }
        LOG_INFO("Loaded {} MIP entries from {}", num_mip_passed, cfg_.root_mip_file);
        tree_in->Delete();
        fin->Close();
        // read spe
        fin->Open(cfg_.root_spe_file.c_str(), "READ");
        if (!fin || fin->IsZombie()) {
            LOG_ERROR("DigitizationAlg: cannot open SPE file: {}", cfg_.root_spe_file);
            return;
        }
        tree_in = std::unique_ptr<TTree>(static_cast<TTree*>(fin->Get("spe")));
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
        const auto spes = df_spe_cut.Take<double>("spe").GetValue();
        if (spe_cellids.size() != spes.size()) {
            LOG_ERROR("SPE tree size mismatch: cellid={} spe={}", spe_cellids.size(), spes.size());
            return;
        }
        int num_spe_passed = 0;
        for (size_t i = 0; i < spe_cellids.size(); ++i) {
            const int conv_cellid = AHCALGeometry::cellid_conversion(spe_cellids[i], cfg_.spe_cellid_version);
            calib_params_[conv_cellid].SPE = spes[i];
            if (spes[i] <= 0) {
                LOG_DEBUG("Invalid SPE value detected: cellID={} SPE={}", conv_cellid, spes[i]);
                calib_params_[conv_cellid].SPE = default_calib_.SPE;
                continue;
            }
            num_spe_passed++;
        }
        tree_in->Delete();
        fin->Close();
        LOG_INFO("Loaded {} SPE entries from {}", num_spe_passed, cfg_.root_spe_file);
        // read sipm model
        fin->Open(cfg_.sipm_model_file.c_str(), "READ");
        if (!fin || fin->IsZombie()) {
            LOG_ERROR("DigitizationAlg: cannot open SiPM model file: {}", cfg_.sipm_model_file);
            return;
        }
        for (int i = 0; i < 15000; i++)
        {
            TString name = "hist/incident_" + TString(std::to_string(i + 1).c_str()) + "_photons";
            h_sipm[i] = std::make_unique<TH1I>(*dynamic_cast<TH1I*>(fin->Get(name)));
            if (!h_sipm[i]) {
                LOG_ERROR("DigitizationAlg: cannot get histogram {} from SiPM model file: {}", name.Data(), cfg_.sipm_model_file);
            }
            h_sipm[i]->SetDirectory(0); // Detach histogram from file to keep it in memory after file is closed
        }
        fin->Close();
    }
    void DigitizationAlg::load_calibration_parameters_fromDB(){
        //tobe implemented
    }
    AHCALRawHit DigitizationAlg::simulate_adc(const AHCALSimHit& simhit, const CalibrationParameters& calib){
        AHCALRawHit rawHit;
        rawHit.cellID = simhit.cellID;
        double energy = rand_gen->Gaus(simhit.Edep, simhit.Edep * 0.06);
        energy = rand_gen->Gaus(energy, energy * 0.015);
        energy = std::max(0.0, energy);
        double n_photons = energy / MIP_E * calib.MIP / calib.SPE / pde;
        n_photons = rand_gen->Poisson(n_photons);

        int n_photoelectrons = rand_gen->Binomial(static_cast<int>(n_photons), pde);
        if (n_photoelectrons >= 15000) {
            n_photoelectrons = 14999;
        }
        int n_fired;
        if (n_photoelectrons == 0) {
            n_fired = 0;
        } else {
            n_fired = int(h_sipm[n_photoelectrons - 1]->GetRandom());
        }
        double adc_value = n_fired * calib.SPE;
        double adc_noise = std::sqrt(n_fired) * 8.0; 

        adc_value = rand_gen->Gaus(adc_value, adc_noise);
        
        double HG = adc_value + rand_gen->Gaus(calib.ped_high, 1.5 * calib.rms_high);
        double LG = adc_value / calib.gain_ratio + rand_gen->Gaus(calib.ped_low, 5 * calib.rms_low);
        if (HG > calib.gain_plat) {
            HG = calib.gain_plat;
        }
        if (LG > calib.lowgain_plat) {
            LG = calib.lowgain_plat;
        }
        rawHit.hg_adc = static_cast<int>(std::round(HG));
        rawHit.lg_adc = static_cast<int>(std::round(LG));
        return rawHit;
    }
}

        

