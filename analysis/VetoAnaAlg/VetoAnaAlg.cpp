#include "VetoAnaAlg.hpp"
#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"
#include "common/edm/SimpleFittedTrack.hpp"
#include "common/edm/Track.hpp"
#include "CalibDBIO/TriggersReader/TriggersReader.hpp"
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

AHCAL_REGISTER_ALG(AHCALRecoAlg::VetoAnaAlg, "VetoAnaAlg")
const double XYMIN = -500;
const double XYMAX = +500;
constexpr int NBIN_XY = 100;
namespace AHCALRecoAlg{

    struct VetoAnaAlg::Impl {
        explicit Impl(VetoAnaAlgCfg cfg, RunContext& rc, std::ofstream* out_file)
            : cfg_(std::move(cfg)), rc_(rc), out_file_(out_file) {}

        // Solve R_reco = R_true(1-exp(-L*R_true)) / (1-exp(-L*R_true)+D*R_true)
        // for R_true with bisection. D and L are in seconds.
        static double solve_true_rate_from_reco(double r_reco, double d_s, double l_s) {
            if (!std::isfinite(r_reco) || r_reco <= 0.0) return 0.0;
            if (!std::isfinite(d_s) || d_s < 0.0 || !std::isfinite(l_s) || l_s <= 0.0) {
                return r_reco;
            }

            const auto reco_from_true = [d_s, l_s](double r_true) {
                const double one_minus_exp = 1.0 - std::exp(-l_s * r_true);
                const double denom = one_minus_exp + d_s * r_true;
                if (denom <= 0.0) return 0.0;
                return (r_true * one_minus_exp) / denom;
            };

            double lo = 0.0;
            double hi = std::max(1.0, r_reco * 2.0 + 1.0);
            for (int i = 0; i < 100 && reco_from_true(hi) < r_reco; ++i) {
                hi *= 2.0;
            }

            for (int i = 0; i < 120; ++i) {
                const double mid = 0.5 * (lo + hi);
                const double reco_mid = reco_from_true(mid);
                if (reco_mid < r_reco) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            return 0.5 * (lo + hi);
        }
    
        // void fill(){
        //     // implementation of the analysis logic goes here
        // }
        void write(){
            if (written_) return;
            written_ = true;
            hTriggerRates_->Scale(1.0 / cfg_.time_window_s); 
            hAnyVetoedRates_->Scale(1.0 / cfg_.time_window_s);
            hNonVetoedRates_->Scale(1.0 / cfg_.time_window_s);
            auto fout = std::unique_ptr<TFile>(TFile::Open(cfg_.out_filename.c_str(), "RECREATE"));
            if (!fout || fout->IsZombie()) {
                LOG_ERROR("VetoAnaAlg: cannot create output file: {}", cfg_.out_filename);
                return;
            }
            LOG_INFO("VetoAnaAlg: writing output to {}", cfg_.out_filename);
            fout->cd();
            hTriggerRates_->Write("TriggerRates");
            hAnyVetoedRates_->Write("AnyVetoedTriggerRates");
            hNonVetoedRates_->Write("NonVetoedTriggerRates");
            int calibRate = static_cast<int>(rc_.conditions.calibRate);
            CalibDBIO::TriggersReader reader(rc_.config.runNumber);
            auto triggers = reader.getTriggers();
            LOG_INFO("VetoAnaAlg: run {}, calibRate {}, triggers: AnyVetoed={}, NonVetoed={}, NonPhysical={}, Physical={}, PreVeto={}, PostVeto={}, Recorded={}", 
                rc_.config.runNumber, calibRate, triggers.AnyVetoed, triggers.NonVetoed, triggers.NonPhysical, triggers.Physical, triggers.PreVeto, triggers.PostVeto, triggers.Recorded);
            double correction_factor_from_NonPhysical = 1.0;
            if (triggers.NonPhysical != -1 && triggers.Physical != -1) {
                double total_time_s = (rc_.conditions.endtime - rc_.conditions.starttime);
                double deadtime_s = cfg_.assumed_deadtime_ms * (triggers.Physical + triggers.NonPhysical) / 1000.0; // convert ms to s
                double live_time_s = total_time_s - deadtime_s;
                if (live_time_s <= 0) {
                    LOG_ERROR("VetoAnaAlg: calculated live time is non-positive ({} s), check the assumed deadtime and trigger counts", live_time_s);
                    return;
                }
                double correction_factor = total_time_s / live_time_s;
                correction_factor_from_NonPhysical = correction_factor;
                LOG_INFO("VetoAnaAlg: total time = {} s, deadtime = {} s, live time = {} s, correction factor from NonPhysical triggers = {}", total_time_s, deadtime_s, live_time_s, correction_factor_from_NonPhysical);
                LOG_INFO("VetoAnaAlg: with assumed deadtime of {} ms", cfg_.assumed_deadtime_ms);              
            }
            TH1D* hCorrectedRates = new TH1D("hCorrectedRates", "Corrected Trigger Rates;Time [s];Rate [Hz]", static_cast<int>((rc_.conditions.endtime - rc_.conditions.starttime)/cfg_.time_window_s), 0, (rc_.conditions.endtime - rc_.conditions.starttime));
            TH1D* hCorrectedAnyVetoedRates = new TH1D("hCorrectedAnyVetoedRates", "Corrected Any Vetoed Trigger Rates;Time [s];Rate [Hz]", static_cast<int>((rc_.conditions.endtime - rc_.conditions.starttime)/cfg_.time_window_s), 0, (rc_.conditions.endtime - rc_.conditions.starttime));
            TH1D* hCorrectedNonVetoedRates = new TH1D("hCorrectedNonVetoedRates", "Corrected Non-Vetoed Trigger Rates;Time [s];Rate [Hz]", static_cast<int>((rc_.conditions.endtime - rc_.conditions.starttime)/cfg_.time_window_s), 0, (rc_.conditions.endtime - rc_.conditions.starttime));
            // double weighted_average_correction_factor = 0.0;
            for (int i = 1; i <= hCorrectedRates->GetNbinsX(); ++i) {
                // double time = hCorrectedRates->GetBinCenter(i);
                double r_reco = hTriggerRates_->GetBinContent(i);
                double r_true = solve_true_rate_from_reco(r_reco, cfg_.assumed_deadtime_ms / 1000.0, cfg_.assumed_auto_triggered_ms/1000.0);
                hCorrectedRates->SetBinContent(i, r_true);
                if (r_true <= 0) continue; // skip empty bins
                double r_any_vetoed_reco = hAnyVetoedRates_->GetBinContent(i);
                double r_any_vetoed_true = r_true/r_reco * r_any_vetoed_reco; // scale by the same factor as the total rate
                if (r_any_vetoed_true <= 0) r_any_vetoed_true = 0; // avoid negative rates due to statistical fluctuations
                hCorrectedAnyVetoedRates->SetBinContent(i, r_any_vetoed_true);
                double r_non_vetoed_reco = hNonVetoedRates_->GetBinContent(i);
                double r_non_vetoed_true = r_true/r_reco * r_non_vetoed_reco; // scale by the same factor as the total rate
                if (r_non_vetoed_true <= 0) r_non_vetoed_true = 0; // avoid negative rates due to statistical fluctuations
                hCorrectedNonVetoedRates->SetBinContent(i, r_non_vetoed_true);
                // average_correction_factor += (r_true / r_reco)*
            }
            // average_correction_factor /= hCorrectedRates->GetNbinsX();
            hCorrectedRates->Write("CorrectedTriggerRates");
            hCorrectedAnyVetoedRates->Write("CorrectedAnyVetoedTriggerRates");
            hCorrectedNonVetoedRates->Write("CorrectedNonVetoedTriggerRates");
            double correction_factor_from_Rates = hCorrectedRates->Integral() / hTriggerRates_->Integral();
            LOG_INFO("VetoAnaAlg: correction factor from rates = {}", correction_factor_from_Rates);  
            // LOG_INFO("VetoAnaAlg: total time = {} s, deadtime = {} s, live time = {} s, correction factor from NonPhysical = {}", total_time_s, deadtime_s, live_time_s, correction_factor_from_NonPhysical);
            // LOG_INFO("VetoAnaAlg: correction factor from NonPhysical triggers = {}", correction_factor_from_NonPhysical);
            double correction_factor = correction_factor_from_Rates; // use the correction factor derived from rates for the final correction
            double corrected_physical_triggers = hCorrectedRates->Integral() * cfg_.time_window_s; // convert back to counts
            double corrected_AnyVetoed_triggers = hCorrectedAnyVetoedRates->Integral() * cfg_.time_window_s;
            double corrected_NonVetoed_triggers = hCorrectedNonVetoedRates->Integral() * cfg_.time_window_s;
            // LOG_INFO("VetoAnaAlg: average correction factor = {}", average_correction_factor);
            if (calibRate > 0) {
                // need to correct the pre-scaling factor for the AnyVetoed, 
                // but the scaling seems to be complicated than 1/x
            }
            if (cfg_.write_to_txt && out_file_ && out_file_->is_open()) {
                (*out_file_) << "RunNumber: " << rc_.config.runNumber << "\n";
                (*out_file_) << "CalibRate: " << calibRate << "\n";
                // (*out_file_) << "TotalTime_s: " << total_time_s << "\n";
                // (*out_file_) << "DeadTime_s: " << deadtime_s << "\n";
                // (*out_file_) << "LiveTime_s: " << live_time_s << "\n";
                (*out_file_) << "CorrectionFactor: " << correction_factor << "\n";
                (*out_file_) << "CorrectionFactorFromNonPhysical: " << correction_factor_from_NonPhysical << "\n";
                (*out_file_) << "CorrectionFactorFromRates: " << correction_factor_from_Rates << "\n";
                (*out_file_) << "PhysicalTriggers: " << triggers.Physical << "\n";
                (*out_file_) << "CorrectedPhysicalTriggers: " << corrected_physical_triggers << "\n";
                (*out_file_) << "AnyVetoedTriggers: " << triggers.AnyVetoed << "\n";
                (*out_file_) << "CorrectedAnyVetoedTriggers: " << corrected_AnyVetoed_triggers << "\n";
                (*out_file_) << "NonVetoedTriggers: " << triggers.NonVetoed << "\n";
                (*out_file_) << "CorrectedNonVetoedTriggers: " << corrected_NonVetoed_triggers << "\n";
                LOG_INFO("VetoAnaAlg: written corrected trigger counts to {}", cfg_.out_txt_name);
            }
            if (cfg_.write_to_png) {
                TCanvas* cRates = new TCanvas("cRates", "Trigger Rates", 1200, 800);
                cRates->Divide(3,1);
                cRates->cd(1)->SetTopMargin(0.1);
                hTriggerRates_->Draw("HIST");
                cRates->cd(2)->SetTopMargin(0.1);
                hAnyVetoedRates_->Draw("HIST");
                cRates->cd(3)->SetTopMargin(0.1);
                hNonVetoedRates_->Draw("HIST");
                cRates->SaveAs((cfg_.out_png_dir + "/Trigger_Rates.png").c_str());
                TCanvas* cCorrectedRates = new TCanvas("cCorrectedRates", "Corrected Trigger Rates", 1200, 800);
                cCorrectedRates->Divide(3,1);
                cCorrectedRates->cd(1)->SetTopMargin(0.1);
                hCorrectedRates->Draw("HIST");
                cCorrectedRates->cd(2)->SetTopMargin(0.1);
                hCorrectedAnyVetoedRates->Draw("HIST");
                cCorrectedRates->cd(3)->SetTopMargin(0.1);
                hCorrectedNonVetoedRates->Draw("HIST");
                cCorrectedRates->SaveAs((cfg_.out_png_dir + "/Corrected_Trigger_Rates.png").c_str());
                TCanvas* cComparison = new TCanvas("cComparison", "Trigger Rates Comparison", 800, 600);
                hTriggerRates_->SetLineColor(kRed);
                hTriggerRates_->Draw("HIST");
                hCorrectedRates->SetLineColor(kBlue);
                hCorrectedRates->Draw("HIST SAME");
                TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
                legend->AddEntry(hTriggerRates_.get(), "Reconstructed Rate", "l");
                legend->AddEntry(hCorrectedRates, "Corrected Rate", "l");
                legend->Draw();
                cComparison->SaveAs((cfg_.out_png_dir + "/Trigger_Rates_Comparison.png").c_str());
            }
            h_deltaTiming_I1I2->Write("DeltaTiming_I1I2");
            h_deltaTiming_I2I3->Write("DeltaTiming_I2I3");
            h_deltaTiming_I3I4->Write("DeltaTiming_I3I4");
            h_deltaTiming_I1I3->Write("DeltaTiming_I1I3");
            h_deltaTiming_I1I4->Write("DeltaTiming_I1I4");
            h_deltaTiming_I2I4->Write("DeltaTiming_I2I4");

            fout->Close();
            LOG_INFO("VetoAnaAlg: wrote {}", cfg_.out_filename);
        }
        void initialize_histograms(int time) {
            hTriggerRates_ = std::make_unique<TH1D>("hTriggerRates", "Trigger Rates;Time [s];Rate [Hz]", static_cast<int>(time/cfg_.time_window_s), 0, time);
            hAnyVetoedRates_ = std::make_unique<TH1D>("hAnyVetoedRates", "Any Vetoed Trigger Rates;Time [s];Rate [Hz]", static_cast<int>(time/cfg_.time_window_s), 0, time);
            hNonVetoedRates_ = std::make_unique<TH1D>("hNonVetoedRates", "Non-Vetoed Trigger Rates;Time [s];Rate [Hz]", static_cast<int>(time/cfg_.time_window_s), 0, time);
            hTriggerRates_->SetDirectory(nullptr);
            hAnyVetoedRates_->SetDirectory(nullptr);
            hNonVetoedRates_->SetDirectory(nullptr);
            h_deltaTiming_I1I2 = std::make_unique<TH2D>("h_deltaTiming_I1I2", "Delta Timing I1 vs I2;#DeltaT(I1 - Veto);#DeltaT(I2 - Veto)", 150, -50, 100, 150, -50, 100);
            h_deltaTiming_I2I3 = std::make_unique<TH2D>("h_deltaTiming_I2I3", "Delta Timing I2 vs I3;#DeltaT(I2 - Veto);#DeltaT(I3 - Veto)", 150, -50, 100, 150, -50, 100);
            h_deltaTiming_I3I4 = std::make_unique<TH2D>("h_deltaTiming_I3I4", "Delta Timing I3 vs I4;#DeltaT(I3 - Veto);#DeltaT(I4 - Veto)", 150, -50, 100, 150, -50, 100);
            h_deltaTiming_I1I3 = std::make_unique<TH2D>("h_deltaTiming_I1I3", "Delta Timing I1 vs I3;#DeltaT(I1 - Veto);#DeltaT(I3 - Veto)", 150, -50, 100, 150, -50, 100);
            h_deltaTiming_I1I4 = std::make_unique<TH2D>("h_deltaTiming_I1I4", "Delta Timing I1 vs I4;#DeltaT(I1 - Veto);#DeltaT(I4 - Veto)", 150, -50, 100, 150, -50, 100);
            h_deltaTiming_I2I4 = std::make_unique<TH2D>("h_deltaTiming_I2I4", "Delta Timing I2 vs I4;#DeltaT(I2 - Veto);#DeltaT(I4 - Veto)", 150, -50, 100, 150, -50, 100);
            h_deltaTiming_I1I2->SetDirectory(nullptr);
            h_deltaTiming_I2I3->SetDirectory(nullptr);
            h_deltaTiming_I3I4->SetDirectory(nullptr);
            h_deltaTiming_I1I3->SetDirectory(nullptr);
            h_deltaTiming_I1I4->SetDirectory(nullptr);
            h_deltaTiming_I2I4->SetDirectory(nullptr);
        }
        VetoAnaAlgCfg cfg_;
        RunContext& rc_;
        std::ofstream* out_file_ = nullptr;
        bool written_ = false;
        std::unique_ptr<TH1D> hTriggerRates_ = nullptr;
        std::unique_ptr<TH1D> hAnyVetoedRates_ = nullptr;
        std::unique_ptr<TH1D> hNonVetoedRates_ = nullptr;
        std::unique_ptr<TH2D> h_deltaTiming_I1I2 = nullptr;
        std::unique_ptr<TH2D> h_deltaTiming_I2I3 = nullptr;
        std::unique_ptr<TH2D> h_deltaTiming_I3I4 = nullptr;
        std::unique_ptr<TH2D> h_deltaTiming_I1I3 = nullptr;
        std::unique_ptr<TH2D> h_deltaTiming_I1I4 = nullptr;
        std::unique_ptr<TH2D> h_deltaTiming_I2I4 = nullptr;

    };
    void VetoAnaAlg::ImplDeleter::operator()(Impl* p) const {
        delete p;
    }
    VetoAnaAlg::~VetoAnaAlg(){
        if (impl_) impl_->write();
    }    
    enum TriggerType : int {
        NoVeto = 1,
        Veto = 2
    };
    static TriggerType getTriggerType(const std::vector<int>& Inputs){
        if (Inputs[4] == 1 || Inputs[5] == 1){
            return TriggerType::Veto;
        } else {
            return TriggerType::NoVeto;
        }
    }
    enum InputType : int {
        I1I2 = 0,
        I2I3 = 1,
        I3I4 = 2,
        I1I3 = 3,
        I1I4 = 4,
        I2I4 = 5
    };
    static std::vector<InputType> getInputType(const std::vector<int>& Inputs){
        std::vector<InputType> result;
        if (Inputs[0] == 1 && Inputs[1] == 1){
            result.push_back(InputType::I1I2);
        } 
        if (Inputs[1] == 1 && Inputs[2] == 1){
            result.push_back(InputType::I2I3);
        }
        if (Inputs[2] == 1 && Inputs[3] == 1){
            result.push_back(InputType::I3I4);
        } 
        if (Inputs[0] == 1 && Inputs[2] == 1){
            result.push_back(InputType::I1I3);
        } 
        if (Inputs[0] == 1 && Inputs[3] == 1){
            result.push_back(InputType::I1I4);
        } 
        if (Inputs[1] == 1 && Inputs[3] == 1){
            result.push_back(InputType::I2I4);
        }
        return result;
    }
    void VetoAnaAlg::execute(EventStore& evt){
        if (!impl_){
            impl_.reset(new Impl(cfg_, this->ctx(), out_file_.get()));
            impl_->initialize_histograms(this->ctx().conditions.endtime - this->ctx().conditions.starttime);
        }
        auto rawdata = evt.get<AHCALTLURawData>(cfg_.in_data_key);
        TriggerType trigger_type = getTriggerType(rawdata.Inputs);
        int Event_Time = rawdata.Event_Time;
        impl_->hTriggerRates_->Fill(Event_Time);
        if (trigger_type == TriggerType::NoVeto) {
            impl_->hNonVetoedRates_->Fill(Event_Time);
        } else if (trigger_type == TriggerType::Veto) {
            impl_->hAnyVetoedRates_->Fill(Event_Time);
            auto input_types = getInputType(rawdata.Inputs);
            double fine_timestamp_veto = 0;
            if (rawdata.Inputs[4] == 1 || rawdata.Inputs[5] == 1) {
                fine_timestamp_veto = (rawdata.FineTimestamps[5]+rawdata.FineTimestamps[4])/2.0; // use the sum of the two veto fine timestamps as a proxy for the veto timing
            } else {
                fine_timestamp_veto = rawdata.Inputs[5]==1 ? rawdata.FineTimestamps[5] : rawdata.FineTimestamps[4]; // if only one of the veto inputs is active, use its fine timestamp
            }
            for (const auto& input_type : input_types) {
                switch (input_type) {
                    case InputType::I1I2:
                        impl_->h_deltaTiming_I1I2->Fill(rawdata.FineTimestamps[0] - fine_timestamp_veto, rawdata.FineTimestamps[1] - fine_timestamp_veto);
                        break;
                    case InputType::I2I3:
                        impl_->h_deltaTiming_I2I3->Fill(rawdata.FineTimestamps[1] - fine_timestamp_veto, rawdata.FineTimestamps[2] - fine_timestamp_veto);
                        break;
                    case InputType::I3I4:
                        impl_->h_deltaTiming_I3I4->Fill(rawdata.FineTimestamps[2] - fine_timestamp_veto, rawdata.FineTimestamps[3] - fine_timestamp_veto);
                        break;
                    case InputType::I1I3:
                        impl_->h_deltaTiming_I1I3->Fill(rawdata.FineTimestamps[0] - fine_timestamp_veto, rawdata.FineTimestamps[2] - fine_timestamp_veto);
                        break;
                    case InputType::I1I4:
                        impl_->h_deltaTiming_I1I4->Fill(rawdata.FineTimestamps[0] - fine_timestamp_veto, rawdata.FineTimestamps[3] - fine_timestamp_veto);
                        break;
                    case InputType::I2I4:
                        impl_->h_deltaTiming_I2I4->Fill(rawdata.FineTimestamps[1] - fine_timestamp_veto, rawdata.FineTimestamps[3] - fine_timestamp_veto);
                        break;
                    default:
                        break;
                }
            }
        }
    }
    void VetoAnaAlg::parse_cfg(const YAML::Node& cfg){
        cfg_.in_data_key = get_or<std::string>(cfg, "in_data_key", cfg_.in_data_key);
        cfg_.out_filename = get_or<std::string>(cfg, "out_filename", cfg_.out_filename);
        cfg_.write_to_png = get_or<bool>(cfg, "write_to_png", cfg_.write_to_png);
        cfg_.out_png_dir = get_or<std::string>(cfg, "out_png_dir", cfg_.out_png_dir);
        cfg_.write_to_txt = get_or<bool>(cfg, "write_to_txt", cfg_.write_to_txt);
        cfg_.out_txt_name = get_or<std::string>(cfg, "out_txt_name", cfg_.out_txt_name);
        cfg_.assumed_deadtime_ms = get_or<double>(cfg, "assumed_deadtime_ms", cfg_.assumed_deadtime_ms);
        cfg_.assumed_auto_triggered_ms = get_or<double>(cfg, "assumed_auto_triggered_ms", cfg_.assumed_auto_triggered_ms);
        cfg_.ntrigger_layer = get_or<int>(cfg, "ntrigger_layer", cfg_.ntrigger_layer);
        cfg_.trigger_layer.clear(); //{1:9, 2:19, 3:29, 4:38}
        std::map<int, int> default_trigger_layer = {{0,9}, {1,19}, {2,29}, {3,38}};
        default_trigger_layer = get_or<std::map<int,int>>(cfg, "trigger_layer", default_trigger_layer);
        for (int i = 0; i < cfg_.ntrigger_layer; ++i) {
            if (default_trigger_layer.find(i) != default_trigger_layer.end()) {
                cfg_.trigger_layer.push_back(default_trigger_layer[i]);
            } else {
                LOG_ERROR("VetoAnaAlg: trigger_layer config is missing entry for input index {}, skipping this input", i);
            }
        }
        cfg_.time_window_s = get_or<double>(cfg, "time_window_s", cfg_.time_window_s);

    }
    void VetoAnaAlg::initialize() {
        if (cfg_.write_to_png) {
            gStyle->SetOptStat(0);
            if (gSystem->AccessPathName(cfg_.out_png_dir.c_str())) {
                if (gSystem->mkdir(cfg_.out_png_dir.c_str(), true) != 0) {
                    LOG_ERROR("VetoAnaAlg: cannot create output PNG directory: {}", cfg_.out_png_dir);
                    cfg_.write_to_png = false;
                } else {
                    LOG_INFO("VetoAnaAlg: created output PNG directory: {}", cfg_.out_png_dir);
                }
            }
        }
        if (cfg_.write_to_txt) {
            out_file_ = std::make_unique<std::ofstream>();
            out_file_->open(cfg_.out_txt_name, std::ios::out);
            if (!out_file_->is_open()) {
                LOG_ERROR("VetoAnaAlg: cannot open output text file: {}", cfg_.out_txt_name);
                out_file_.reset();
            } else {
                LOG_INFO("VetoAnaAlg: writing detailed event info to {}", cfg_.out_txt_name);
            }
        }
    }
} // namespace AHCALRecoAlg