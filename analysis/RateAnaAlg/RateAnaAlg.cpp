#include "RateAnaAlg.hpp"
#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"
#include "common/edm/SimpleFittedTrack.hpp"
#include "common/edm/Track.hpp"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLegend.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <fstream>

AHCAL_REGISTER_ALG(AHCALRecoAlg::RateAnaAlg, "RateAnaAlg")
const double XYMIN = -500;
const double XYMAX = +500;
constexpr int NBIN_XY = 100;
namespace AHCALRecoAlg{

    struct RateAnaAlg::Impl {
        explicit Impl(RateAnaAlgCfg cfg) : cfg_(std::move(cfg)) {}
    
        // void fill(){
        //     // implementation of the analysis logic goes here
        // }
        void write(int calibRate) {
            if (written_) return;
            written_ = true;

            auto fout = std::unique_ptr<TFile>(TFile::Open(cfg_.out_filename.c_str(), "RECREATE"));
            if (!fout || fout->IsZombie()) {
                LOG_ERROR("RateAnaAlg: cannot create output file: {}", cfg_.out_filename);
                return;
            }
            LOG_INFO("RateAnaAlg: writing output to {}", cfg_.out_filename);
            fout->cd();
            h_deltaCycleID->Write();
            h_deltaTimestamp->Write();
            h_deltaTriggerID->Write();
            h_full_triggerRate->Scale(1.0 / cfg_.time_window_s); // convert counts to rate
            h_full_triggerRate->Write();
            h_noveto_triggerRate->Scale(1.0 / cfg_.time_window_s);
            h_noveto_triggerRate->Write();
            h_veto_triggerRate->Scale(1.0 / cfg_.time_window_s);
            h_veto_triggerRate->Write();

            if (cfg_.write_to_png) {
                if (gSystem->AccessPathName(cfg_.out_png_dir.c_str())) {
                    if (gSystem->mkdir(cfg_.out_png_dir.c_str(), true) != 0) {
                        LOG_ERROR("RateAnaAlg: cannot create output PNG directory: {}", cfg_.out_png_dir);
                        return;
                    }
                }
                TCanvas c("c_trigger_rate", "Trigger Rate", 800, 600);
                h_full_triggerRate->Draw("hist");
                c.SaveAs(Form("%s/full_trigger_rate.png", cfg_.out_png_dir.c_str()));
                h_noveto_triggerRate->Draw("hist");
                c.SaveAs(Form("%s/noveto_trigger_rate.png", cfg_.out_png_dir.c_str()));
                h_veto_triggerRate->Draw("hist");
                c.SaveAs(Form("%s/veto_trigger_rate.png", cfg_.out_png_dir.c_str()));
            }
            TGraphAsymmErrors* g_full = new TGraphAsymmErrors();
            g_full->SetName("g_full");
            g_full->SetTitle("Full Trigger Rate;CalibRate;Rate (Hz)");
            g_full->SetPoint(0, calibRate, h_full_triggerRate->GetEntries()/cfg_.time_window_s/h_full_triggerRate->GetNbinsX());
            g_full->SetPointError(0, 0, 0, sqrt(h_full_triggerRate->GetEntries())/cfg_.time_window_s/h_full_triggerRate->GetNbinsX(), sqrt(h_full_triggerRate->GetEntries())/cfg_.time_window_s/h_full_triggerRate->GetNbinsX());
            g_full->Write();
            TGraphAsymmErrors* g_noveto = new TGraphAsymmErrors();
            g_noveto->SetName("g_noveto");
            g_noveto->SetTitle("No-veto Trigger Rate;CalibRate;Rate (Hz)");
            g_noveto->SetPoint(0, calibRate, h_noveto_triggerRate->GetEntries()/cfg_.time_window_s/h_noveto_triggerRate->GetNbinsX());
            g_noveto->SetPointError(0, 0, 0, sqrt(h_noveto_triggerRate->GetEntries())/cfg_.time_window_s/h_noveto_triggerRate->GetNbinsX(), sqrt(h_noveto_triggerRate->GetEntries())/cfg_.time_window_s/h_noveto_triggerRate->GetNbinsX());
            g_noveto->Write();
            TGraphAsymmErrors* g_veto = new TGraphAsymmErrors();
            g_veto->SetName("g_veto");
            g_veto->SetTitle("Veto Trigger Rate;CalibRate;Rate (Hz)");
            g_veto->SetPoint(0, calibRate, h_veto_triggerRate->GetEntries()/cfg_.time_window_s/h_veto_triggerRate->GetNbinsX());
            g_veto->SetPointError(0, 0, 0, sqrt(h_veto_triggerRate->GetEntries())/cfg_.time_window_s/h_veto_triggerRate->GetNbinsX(), sqrt(h_veto_triggerRate->GetEntries())/cfg_.time_window_s/h_veto_triggerRate->GetNbinsX());
            g_veto->Write();
            fout->Close();
            LOG_INFO("RateAnaAlg: wrote {}", cfg_.out_filename);
        }

        RateAnaAlgCfg cfg_;
        bool written_ = false;
        void initialize_histograms(int time_window_s, double start_time, double end_time) {
            h_deltaCycleID = std::make_unique<TH1D>("h_deltaCycleID", "Distribution of CycleID difference between last;#Delta CycleID;Events", 1000, 0, 1000);
            h_deltaTimestamp = std::make_unique<TH1D>("h_deltaTimestamp", "Distribution of Timestamp difference between last;#Delta Timestamp;Events", 1000, 0, 1e6);
            h_deltaTriggerID = std::make_unique<TH1D>("h_deltaTriggerID", "Distribution of TriggerID difference between last;#Delta TriggerID;Events", 5, -0.5, 4.5);
            h_deltaCycleID->SetDirectory(nullptr);
            h_deltaTimestamp->SetDirectory(nullptr);
            h_deltaTriggerID->SetDirectory(nullptr);
            // g_full_triggerRate->SetDirectory(nullptr);
            h_full_triggerRate = std::make_unique<TH1D>("h_full_triggerRate", "Trigger rate for all triggers;Time (s);Rate (Hz)", static_cast<int>((end_time - start_time) / time_window_s), 0, static_cast<int>((end_time - start_time) / time_window_s) * time_window_s);
            h_noveto_triggerRate = std::make_unique<TH1D>("h_noveto_triggerRate", "Trigger rate for no-veto triggers;Time (s);Rate (Hz)", static_cast<int>((end_time - start_time) / time_window_s), 0, static_cast<int>((end_time - start_time) / time_window_s) * time_window_s);
            h_veto_triggerRate = std::make_unique<TH1D>("h_veto_triggerRate", "Trigger rate for veto triggers;Time (s);Rate (Hz)", static_cast<int>((end_time - start_time) / time_window_s), 0, static_cast<int>((end_time - start_time) / time_window_s) * time_window_s);
        }
        std::unique_ptr<TH1D> h_deltaCycleID;
        std::unique_ptr<TH1D> h_deltaTimestamp;
        std::unique_ptr<TH1D> h_deltaTriggerID;
        std::unique_ptr<TH1D> h_full_triggerRate;
        std::unique_ptr<TH1D> h_noveto_triggerRate;
        std::unique_ptr<TH1D> h_veto_triggerRate;
    };
    void RateAnaAlg::ImplDeleter::operator()(Impl* p) const {
        delete p;
    }
    enum TriggerType : int {
        NoVeto = 1,
        Veto = 2
    };
    TriggerType getTriggerType(const std::vector<int>& Inputs){
        if (Inputs[4] == 1 || Inputs[5] == 1){
            return TriggerType::Veto;
        } else {
            return TriggerType::NoVeto;
        }
    }
    RateAnaAlg::~RateAnaAlg(){
        if (impl_) impl_->write(this->ctx().conditions.calibRate);
    }
    void RateAnaAlg::execute(EventStore& evt){
        if (!impl_){
            impl_.reset(new Impl(cfg_));
            impl_->initialize_histograms(cfg_.time_window_s, this->ctx().conditions.starttime, this->ctx().conditions.endtime);
        }
        auto rawdata = evt.get<AHCALTLURawData>(cfg_.in_data_key);
        TriggerType trigger_type = getTriggerType(rawdata.Inputs);
        int Event_Time = rawdata.Event_Time;
        impl_->h_full_triggerRate->Fill(Event_Time);
        if (trigger_type == TriggerType::NoVeto) {
            impl_->h_noveto_triggerRate->Fill(Event_Time);
        } else if (trigger_type == TriggerType::Veto) {
            impl_->h_veto_triggerRate->Fill(Event_Time);
        }
        int current_cycle_id = rawdata.CycleID;
        int current_TLU_timestamp = rawdata.Timestamp;
        int current_trigger_id = rawdata.TriggerID;
        if (last_cycle_id != -1 && last_TLU_timestamp != -1) {
            int delta_cycle_id = current_cycle_id - last_cycle_id;
            int delta_timestamp = current_TLU_timestamp - last_TLU_timestamp;
            int delta_trigger_id = current_trigger_id - last_trigger_id;
            impl_->h_deltaCycleID->Fill(delta_cycle_id);
            impl_->h_deltaTimestamp->Fill(delta_timestamp);
            impl_->h_deltaTriggerID->Fill(delta_trigger_id);
        }
        last_cycle_id = current_cycle_id;
        last_TLU_timestamp = current_TLU_timestamp;
        last_trigger_id = current_trigger_id;
    }
    void RateAnaAlg::parse_cfg(const YAML::Node& cfg){
        cfg_.in_data_key = get_or<std::string>(cfg, "in_data_key", cfg_.in_data_key);
        cfg_.out_filename = get_or<std::string>(cfg, "out_filename", cfg_.out_filename);
        cfg_.write_to_png = get_or<bool>(cfg, "write_to_png", cfg_.write_to_png);
        cfg_.out_png_dir = get_or<std::string>(cfg, "out_png_dir", cfg_.out_png_dir);
        cfg_.write_to_txt = get_or<bool>(cfg, "write_to_txt", cfg_.write_to_txt);
        cfg_.out_txt_name = get_or<std::string>(cfg, "out_txt_name", cfg_.out_txt_name);
        cfg_.ntrigger_layer = get_or<int>(cfg, "ntrigger_layer", cfg_.ntrigger_layer);
        cfg_.trigger_layer.clear(); //{1:9, 2:19, 3:29, 4:38}
        std::map<int, int> default_trigger_layer = {{0,9}, {1,19}, {2,29}, {3,38}};
        default_trigger_layer = get_or<std::map<int,int>>(cfg, "trigger_layer", default_trigger_layer);
        for (int i = 0; i < cfg_.ntrigger_layer; ++i) {
            if (default_trigger_layer.find(i) != default_trigger_layer.end()) {
                cfg_.trigger_layer.push_back(default_trigger_layer[i]);
            } else {
                LOG_ERROR("RateAnaAlg: trigger_layer config is missing entry for input index {}, skipping this input", i);
            }
        }
        cfg_.time_window_s = get_or<double>(cfg, "time_window_s", cfg_.time_window_s);
    }
    void RateAnaAlg::initialize() {
        if (cfg_.write_to_txt) {
            out_file_ = std::make_unique<std::ofstream>();
            out_file_->open(cfg_.out_txt_name, std::ios::out);
            if (!out_file_->is_open()) {
                LOG_ERROR("RateAnaAlg: cannot open output text file: {}", cfg_.out_txt_name);
                out_file_.reset();
            } else {
                LOG_INFO("RateAnaAlg: writing detailed event info to {}", cfg_.out_txt_name);
            }
        }
    }
} // namespace AHCALRecoAlg