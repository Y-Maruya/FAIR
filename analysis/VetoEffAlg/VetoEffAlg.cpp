#include "VetoEffAlg.hpp"
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

AHCAL_REGISTER_ALG(AHCALRecoAlg::VetoEffAlg, "VetoEffAlg")
const double XYMIN = -500;
const double XYMAX = +500;
constexpr int NBIN_XY = 100;
namespace AHCALRecoAlg{

    struct VetoEffAlg::Impl {
        explicit Impl(VetoEffAlgCfg cfg) : cfg_(std::move(cfg)) {}
    
        void fill(double x, double y, int veto, bool is_in_track = true) {
            if (veto != 0 && veto != 1) {
                LOG_ERROR("VetoEffAlg::Impl::fill: invalid veto value {}, expected 0 or 1", veto);
                return;
            }
            h2_fulls[veto]->Fill(x, y);
            if (abs(x-cfg_.x_center) <= cfg_.xy_size_threshold && abs(y-cfg_.y_center) <= cfg_.xy_size_threshold) {
                h_full->Fill(veto);
                if (is_in_track) h_passed->Fill(veto);
            }
            if (is_in_track) h2_passeds[veto]->Fill(x, y);
        }
        
        void write() {
            if (written_) return;
            written_ = true;

            auto fout = std::unique_ptr<TFile>(TFile::Open(cfg_.out_filename.c_str(), "RECREATE"));
            if (!fout || fout->IsZombie()) {
                LOG_ERROR("VetoEffAlg: cannot create output file: {}", cfg_.out_filename);
                return;
            }
            LOG_INFO("VetoEffAlg: writing output to {}", cfg_.out_filename);
            fout->cd();
            for (size_t i = 0; i < h2_fulls.size(); ++i) {
                h2_fulls[i]->Write();
                h2_passeds[i]->Write();
            }
            LOG_DEBUG("VetoEffAlg: h2_fulls and h2_passeds written to file");
            std::vector<std::unique_ptr<TEfficiency>> effs;
            effs.reserve(h2_fulls.size());
            for (size_t i = 0; i < h2_fulls.size(); ++i) {
                auto eff = std::make_unique<TEfficiency>(*h2_passeds[i], *h2_fulls[i]);
                // eff->SetDirectory(nullptr);
                effs.push_back(std::move(eff));
            }
            std::string title = "Veto Scintillator Efficiency";
            for (size_t i = 0; i < effs.size(); ++i) {
                title = "Veto Scintillator Efficiency: Input" + std::to_string(i+4) + " (threshold " + std::to_string(i==0 ? static_cast<int>(cfg_.threshold_input4) : static_cast<int>(cfg_.threshold_input5)) + ")";
                effs[i]->SetName(Form("eff_input%zu", i+4));
                effs[i]->SetTitle(Form("%s;X [mm];Y [mm]", title.c_str()));
            }
            for (size_t i = 0; i < effs.size(); ++i) {
                effs[i]->Write();
            }
            if (cfg_.write_to_png) {
                if (gSystem->AccessPathName(cfg_.out_png_dir.c_str())) {
                    if (gSystem->mkdir(cfg_.out_png_dir.c_str()) != 0) {
                        LOG_ERROR("VetoEffAlg: cannot create output PNG directory: {}", cfg_.out_png_dir);
                        return;
                    }
                }
                gStyle->SetOptStat(0);
                gStyle->SetPaintTextFormat("4.2f%%");
                for (size_t i = 0; i < effs.size(); ++i) {
                    TCanvas c(Form("c_eff_input%zu", i+4), Form("Efficiency for Input%zu", i+4), 800, 600);
                    effs[i]->Draw("COLZ");
                    // drawLayerLabel(static_cast<int>(i+4));
                    c.SaveAs(Form("%s/eff_input%zu.png", cfg_.out_png_dir.c_str(), i+4));
                }
            }
            std::unique_ptr<TEfficiency> overall_eff = std::make_unique<TEfficiency>(*h_passed, *h_full);
            // overall_eff->SetDirectory(nullptr);
            overall_eff->SetName("overall_efficiency");
            overall_eff->SetTitle(Form("Overall Veto Efficiency within |x|<%.1f mm, |y|<%.1f mm>;Veto;Efficiency", cfg_.xy_size_threshold, cfg_.xy_size_threshold));
            overall_eff->Draw("AP");
            overall_eff->Write();
            gPad->Update();
            auto graph = overall_eff->GetPaintedGraph();
            if (graph) {
                graph->GetXaxis()->SetBinLabel(1, "Input4");
                graph->GetXaxis()->SetBinLabel(2, "Input5");
            }
            if (cfg_.write_to_png) {
                TCanvas c("c_overall_efficiency", "Overall Veto Efficiency", 800, 600);
                TH1D* h_dummy = new TH1D("h_dummy", overall_eff->GetTitle(), 2, -0.5, 1.5);
                h_dummy->GetXaxis()->SetBinLabel(1, "Input4");
                h_dummy->GetXaxis()->SetBinLabel(2, "Input5");
                h_dummy->SetDirectory(nullptr);
                h_dummy->Draw();
                h_dummy->GetYaxis()->SetRangeUser(0.95, 1.0);
                graph->Draw("Psame");
                c.SaveAs(Form("%s/overall_efficiency.png", cfg_.out_png_dir.c_str()));
            }
            fout->Close();
            LOG_INFO("VetoEffAlg: wrote {}", cfg_.out_filename);
        }

        VetoEffAlgCfg cfg_;
        bool written_ = false;
        void initialize_histograms() {
            h2_fulls.clear();
            h2_passeds.clear();
            h2_fulls.reserve(2);
            h2_passeds.reserve(2);
            h2_fulls.push_back(std::make_unique<TH2D>("h2_full_input4", "Full input distribution;x[mm];y[mm]", NBIN_XY, XYMIN, XYMAX, NBIN_XY, XYMIN, XYMAX));
            h2_passeds.push_back(std::make_unique<TH2D>("h2_passed_input4", "Passed input distribution;x[mm];y[mm]", NBIN_XY, XYMIN, XYMAX, NBIN_XY, XYMIN, XYMAX));
            h2_fulls[0]->SetDirectory(nullptr);
            h2_passeds[0]->SetDirectory(nullptr);
            h2_fulls.push_back(std::make_unique<TH2D>("h2_full_input5", "Full input distribution;x[mm];y[mm]", NBIN_XY, XYMIN, XYMAX, NBIN_XY, XYMIN, XYMAX));
            h2_passeds.push_back(std::make_unique<TH2D>("h2_passed_input5", "Passed input distribution;x[mm];y[mm]", NBIN_XY, XYMIN, XYMAX, NBIN_XY, XYMIN, XYMAX));
            h2_fulls[1]->SetDirectory(nullptr);
            h2_passeds[1]->SetDirectory(nullptr);
            
            h_full = std::make_unique<TH1D>("h_full", "Full distribution;Veto;Entries", 2, -0.5, 1.5);
            h_passed = std::make_unique<TH1D>("h_passed", "Passed distribution;Veto;Entries", 2, -0.5, 1.5);
            h_full->SetDirectory(nullptr);
            h_passed->SetDirectory(nullptr);
        }
        std::vector<std::unique_ptr<TH2D>> h2_fulls;
        std::vector<std::unique_ptr<TH2D>> h2_passeds;
        std::unique_ptr<TH1D> h_full;
        std::unique_ptr<TH1D> h_passed;
    };

    void VetoEffAlg::ImplDeleter::operator()(Impl* p) const {
        delete p;
    }
    VetoEffAlg::~VetoEffAlg(){
        if (impl_) impl_->write();
    }

    void VetoEffAlg::execute(EventStore& evt){
        if (!impl_){
            impl_.reset(new Impl(cfg_));
            impl_->initialize_histograms();
        }
        if (cfg_.string_track_struct == "SimpleFittedTrack") {
            auto track = evt.get<SimpleFittedTrack>(cfg_.in_track_key);
            if (!cfg_.track_selection_string.empty()) {
                SimpleFittedTrackCut cut(cfg_.track_selection_string);
                if (!cut.eval(track)) {
                    LOG_DEBUG("VetoEffAlg: track did not pass selection cut '{}'", cfg_.track_selection_string);
                    return;
                }
            }
            auto rawdata = evt.get<AHCALTLURawData>(cfg_.in_data_key);
            for (int i = 0; i < 2; ++i){
                double x = track.init_pos_x + track.direction_x * (cfg_.z_pos_input4 + i*(cfg_.z_pos_input5-cfg_.z_pos_input4));
                double y = track.init_pos_y + track.direction_y * (cfg_.z_pos_input4 + i*(cfg_.z_pos_input5-cfg_.z_pos_input4));
                bool Input = (rawdata.Inputs[i+4] == 1);
                impl_->fill(x, y, i, Input);
            }
        } else if (cfg_.string_track_struct == "Track") {
            auto track = evt.get<Track>(cfg_.in_track_key);
            if (!cfg_.track_selection_string.empty()) {
                TrackCut cut(cfg_.track_selection_string);
                if (!cut.eval(track)) {
                    LOG_DEBUG("VetoEffAlg: track did not pass selection cut '{}'", cfg_.track_selection_string);
                    return;
                }
            }
            auto rawdata = evt.get<AHCALTLURawData>(cfg_.in_data_key);
            for (int i = 0; i < 2; ++i){
                double x = track.x + track.tx * (cfg_.z_pos_input4 + i*(cfg_.z_pos_input5-cfg_.z_pos_input4)-track.z);
                double y = track.y + track.ty * (cfg_.z_pos_input4 + i*(cfg_.z_pos_input5-cfg_.z_pos_input4)-track.z);
                bool Input = (rawdata.Inputs[i+4] == 1);
                impl_->fill(x, y, i, Input);
            }
        } else {
            LOG_ERROR("VetoEffAlg: unknown track struct string '{}'", cfg_.string_track_struct);
            return;
        }
    }
    void VetoEffAlg::parse_cfg(const YAML::Node& cfg){
        cfg_.in_data_key = get_or<std::string>(cfg, "in_data_key", cfg_.in_data_key);
        cfg_.in_track_key = get_or<std::string>(cfg, "in_track_key", cfg_.in_track_key);
        cfg_.string_track_struct = get_or<std::string>(cfg, "string_track_struct", cfg_.string_track_struct);
        cfg_.track_selection_string = get_or<std::string>(cfg, "track_selection_string", cfg_.track_selection_string);
        cfg_.out_filename = get_or<std::string>(cfg, "out_filename", cfg_.out_filename);
        cfg_.write_to_png = get_or<bool>(cfg, "write_to_png", cfg_.write_to_png);
        cfg_.out_png_dir = get_or<std::string>(cfg, "out_png_dir", cfg_.out_png_dir);
        cfg_.xy_size_threshold = get_or<double>(cfg, "xy_size_threshold", cfg_.xy_size_threshold);
        cfg_.x_center = get_or<double>(cfg, "x_center", cfg_.x_center);
        cfg_.y_center = get_or<double>(cfg, "y_center", cfg_.y_center);
        cfg_.z_pos_input4 = get_or<double>(cfg, "z_pos_input4", cfg_.z_pos_input4);
        cfg_.z_pos_input5 = get_or<double>(cfg, "z_pos_input5", cfg_.z_pos_input5);
        cfg_.threshold_input4 = get_or<double>(cfg, "threshold_input4", cfg_.threshold_input4);
        cfg_.threshold_input5 = get_or<double>(cfg, "threshold_input5", cfg_.threshold_input5);
    }
} // namespace AHCALRecoAlg