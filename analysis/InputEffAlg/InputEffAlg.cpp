#include "InputEffAlg.hpp"
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

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <fstream>

AHCAL_REGISTER_ALG(AHCALRecoAlg::InputEffAlg, "InputEffAlg")
const double XYMIN = -500;
const double XYMAX = +500;
constexpr int NBIN_XY = 100;
namespace AHCALRecoAlg{

    struct InputEffAlg::Impl {
        explicit Impl(InputEffAlgCfg cfg) : cfg_(std::move(cfg)) {}
    
        void fill_inputs(const std::vector<int> inputs) {
            h_input_sum->Fill(std::accumulate(inputs.begin(), inputs.end(), 0));
        }
        void fill_single(int input_index,bool is_has_hit) {
            h_full->Fill(input_index);
            if (is_has_hit) {
                h_passed->Fill(input_index);
            }
        }
        void fill_single_input1(int input_index, int nHits, unsigned long long event_number, const RunContext& ctx, bool is_has_hit, std::ofstream* out_file) {
            h_inputsumeqq1->Fill(input_index);
            h_nHits_sum_inputs_eq_1->Fill(nHits);
            if (cfg_.write_to_txt && out_file) {
                (*out_file) << ctx.config.runNumber << " " << ctx.config.poolIndex << " " << event_number << " " << input_index << " " << is_has_hit << " " << nHits << "\n";
            }
            if (nHits > 20){
                LOG_WARN("InputEffAlg::Impl::fill_single_input1: event {} has {} hits, which is unusually high. Check if this is expected.", event_number, nHits);
                std::cout << ctx.config.runNumber << " " << ctx.config.poolIndex << " " << event_number << " " << input_index << " " << nHits << "\n";
            }
        }
        void write() {
            if (written_) return;
            written_ = true;

            auto fout = std::unique_ptr<TFile>(TFile::Open(cfg_.out_filename.c_str(), "RECREATE"));
            if (!fout || fout->IsZombie()) {
                LOG_ERROR("InputEffAlg: cannot create output file: {}", cfg_.out_filename);
                return;
            }
            LOG_INFO("InputEffAlg: writing output to {}", cfg_.out_filename);
            fout->cd();
            h_full->Write();
            h_passed->Write();
            h_input_sum->Write();
            h_inputsumeqq1->Write();
            h_nHits_sum_inputs_eq_1->Write();
            std::unique_ptr<TEfficiency> eff = std::make_unique<TEfficiency>(*h_passed, *h_full);
            eff->SetName("input_efficiency");
            eff->SetTitle("HitTag Efficiency for Inputs;Input;Efficiency");
            eff->Draw("AP");
            eff->Write();
            gPad->Update();
            auto graph = eff->GetPaintedGraph();
            if (graph) {
                graph->GetXaxis()->SetBinLabel(1, Form("Input0 L%d", cfg_.trigger_layer[0]));
                graph->GetXaxis()->SetBinLabel(2, Form("Input1 L%d", cfg_.trigger_layer[1]));
                graph->GetXaxis()->SetBinLabel(3, Form("Input2 L%d", cfg_.trigger_layer[2]));
                graph->GetXaxis()->SetBinLabel(4, Form("Input3 L%d", cfg_.trigger_layer[3]));
                graph->GetXaxis()->SetBinLabel(5, "Input4 Veto1");
                graph->GetXaxis()->SetBinLabel(6, "Input5 Veto2");
            }
            if (cfg_.write_to_png) {
                if (gSystem->AccessPathName(cfg_.out_png_dir.c_str())) {
                    if (gSystem->mkdir(cfg_.out_png_dir.c_str(), true) != 0) {
                        LOG_ERROR("InputEffAlg: cannot create output PNG directory: {}", cfg_.out_png_dir);
                        return;
                    }
                }
                TCanvas c("c_input_efficiency", "HitTag Efficiency for Inputs", 800, 600);
                TH1D* h_dummy = new TH1D("h_dummy", eff->GetTitle(), 6, -0.5, 5.5);
                h_dummy->GetXaxis()->SetBinLabel(1, Form("Input0 L%d", cfg_.trigger_layer[0]));
                h_dummy->GetXaxis()->SetBinLabel(2, Form("Input1 L%d", cfg_.trigger_layer[1]));
                h_dummy->GetXaxis()->SetBinLabel(3, Form("Input2 L%d", cfg_.trigger_layer[2]));
                h_dummy->GetXaxis()->SetBinLabel(4, Form("Input3 L%d", cfg_.trigger_layer[3]));
                h_dummy->GetXaxis()->SetBinLabel(5, "Input4 Veto1");
                h_dummy->GetXaxis()->SetBinLabel(6, "Input5 Veto2");
                h_dummy->Draw();
                eff->Draw("Psame");
                c.SaveAs(Form("%s/input_efficiency.png", cfg_.out_png_dir.c_str()));
                h_full->Draw();
                h_full->GetXaxis()->SetBinLabel(1, Form("Input0 L%d", cfg_.trigger_layer[0]));
                h_full->GetXaxis()->SetBinLabel(2, Form("Input1 L%d", cfg_.trigger_layer[1]));
                h_full->GetXaxis()->SetBinLabel(3, Form("Input2 L%d", cfg_.trigger_layer[2]));
                h_full->GetXaxis()->SetBinLabel(4, Form("Input3 L%d", cfg_.trigger_layer[3]));
                h_full->GetXaxis()->SetBinLabel(5, "Input4 Veto1");
                h_full->GetXaxis()->SetBinLabel(6, "Input5 Veto2");
                h_full->Draw();
                c.SaveAs(Form("%s/input_full.png", cfg_.out_png_dir.c_str()));
                h_input_sum->Draw();
                c.SaveAs(Form("%s/input_sum.png", cfg_.out_png_dir.c_str()));
                h_inputsumeqq1->Draw();
                c.SaveAs(Form("%s/input_sum_eq_1.png", cfg_.out_png_dir.c_str()));
                h_nHits_sum_inputs_eq_1->Draw();
                c.SaveAs(Form("%s/nHits_sum_input_eq_1.png", cfg_.out_png_dir.c_str()));
            }
        
            fout->Close();
            LOG_INFO("InputEffAlg: wrote {}", cfg_.out_filename);
        }

        InputEffAlgCfg cfg_;
        bool written_ = false;
        void initialize_histograms() {
            h_full = std::make_unique<TH1D>("h_full", "Denominator for HitTag efficiency;Input;Events", 6, -0.5, 5.5);
            h_passed = std::make_unique<TH1D>("h_passed", "Numerator for HitTag efficiency;Input;Events", 6, -0.5, 5.5);
            h_input_sum = std::make_unique<TH1D>("h_input_sum", "Distribution of sum of Inputs;Sum of Inputs;Events", 6, -0.5, 5.5);
            h_inputsumeqq1 = std::make_unique<TH1D>("h_inputsumeqq1", "Distribution of HitTag for events with sum of Inputs == 1;Input;Events", 6, -0.5, 5.5);
            h_nHits_sum_inputs_eq_1 = std::make_unique<TH1D>("h_nHits_sum_inputs_eq_1", "Distribution of nHits for events with sum of Inputs == 1;nHits;Events", 100, 0, 100);
            h_full->SetDirectory(nullptr);
            h_passed->SetDirectory(nullptr);
            h_input_sum->SetDirectory(nullptr);
            h_inputsumeqq1->SetDirectory(nullptr);
            h_nHits_sum_inputs_eq_1->SetDirectory(nullptr);
        }

        std::unique_ptr<TH1D> h_full; // HitTag efficiency denominator 
        std::unique_ptr<TH1D> h_passed; // HitTag efficiency numerator 
        std::unique_ptr<TH1D> h_input_sum; // distribution of sum of inputs, for sanity check
        std::unique_ptr<TH1D> h_inputsumeqq1; // distribution of sum of inputs == 1, for sanity check
        std::unique_ptr<TH1D> h_nHits_sum_inputs_eq_1; // distribution of nHits for events with sum of inputs == 1, for sanity check
    };

    void InputEffAlg::ImplDeleter::operator()(Impl* p) const {
        delete p;
    }
    InputEffAlg::~InputEffAlg(){
        if (impl_) impl_->write();
    }
    bool InputEffAlg::is_hittag1_exist(std::vector<AHCALRecoHit>& recoHits, int layer) {
        for (const auto& hit : recoHits) {
            if (hit.layer() == layer) {
                return true;
            }
        }
        return false;
    } 
    void InputEffAlg::execute(EventStore& evt){
        if (!impl_){
            impl_.reset(new Impl(cfg_));
            impl_->initialize_histograms();
        }
        auto rawdata = evt.get<AHCALTLURawData>(cfg_.in_data_key);
        auto recoHits = evt.get<std::vector<AHCALRecoHit>>(cfg_.in_reco_hit_key);
        impl_->fill_inputs(rawdata.Inputs);
        // categorize the inputs 
        for (int l = 0; l < cfg_.ntrigger_layer; ++l){
            if (rawdata.Inputs[l]) {
                bool is_has_hit = is_hittag1_exist(recoHits, cfg_.trigger_layer[l]);
                impl_->fill_single(l, is_has_hit);
            }
        }
        // sum_inputs == 1 analysis
        if (std::accumulate(rawdata.Inputs.begin(), rawdata.Inputs.end(), 0) == 1) {
            for (int l = 0; l < cfg_.ntrigger_layer; ++l){
                if (rawdata.Inputs[l]) {
                    bool is_has_hit = is_hittag1_exist(recoHits, cfg_.trigger_layer[l]);
                    int nHits = recoHits.size();
                    impl_->fill_single_input1(l, nHits, evt.event_counter(), this->ctx(), is_has_hit, out_file_.get());
                }
            }
        }
    }
    void InputEffAlg::parse_cfg(const YAML::Node& cfg){
        cfg_.in_data_key = get_or<std::string>(cfg, "in_data_key", cfg_.in_data_key);
        cfg_.in_reco_hit_key = get_or<std::string>(cfg, "in_reco_hit_key", cfg_.in_reco_hit_key);
        // cfg_.string_track_struct = get_or<std::string>(cfg, "string_track_struct", cfg_.string_track_struct);
        // cfg_.track_selection_string = get_or<std::string>(cfg, "track_selection_string", cfg_.track_selection_string);
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
                LOG_ERROR("InputEffAlg: trigger_layer config is missing entry for input index {}, skipping this input", i);
            }
        }
    }
    void InputEffAlg::initialize() {
        if (cfg_.write_to_txt) {
            out_file_ = std::make_unique<std::ofstream>();
            out_file_->open(cfg_.out_txt_name, std::ios::out);
            if (!out_file_->is_open()) {
                LOG_ERROR("InputEffAlg: cannot open output text file: {}", cfg_.out_txt_name);
                out_file_.reset();
            } else {
                LOG_INFO("InputEffAlg: writing detailed event info to {}", cfg_.out_txt_name);
            }
        }
    }
} // namespace AHCALRecoAlg