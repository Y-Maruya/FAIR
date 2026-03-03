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
#include <fstream>

AHCAL_REGISTER_ALG(AHCALRecoAlg::VetoEffAlg, "VetoEffAlg")
const double XYMIN = -500;
const double XYMAX = +500;
constexpr int NBIN_XY = 100;
namespace AHCALRecoAlg{

    struct VetoEffAlg::Impl {
        explicit Impl(VetoEffAlgCfg cfg) : cfg_(std::move(cfg)) {}
    
        void fill(double x, double y, int veto, int evtNum, RunContext& ctx, bool is_in_track, std::ofstream* out_file) {
            if (veto != 0 && veto != 1) {
                LOG_ERROR("VetoEffAlg::Impl::fill: invalid veto value {}, expected 0 or 1", veto);
                return;
            }
            h2_fulls[veto]->Fill(x, y);

            if (abs(x-cfg_.x_center) <= cfg_.xy_size_threshold && abs(y-cfg_.y_center) <= cfg_.xy_size_threshold) {
                h_full->Fill(veto);
                if (is_in_track) h_passed->Fill(veto);
                else if (cfg_.write_to_txt && out_file) {
                    (*out_file) << ctx.config.runNumber << " " << ctx.config.poolIndex << " " << evtNum << " " << x << " " << y << " " << veto << "\n";
                }
            }
            if (is_in_track) h2_passeds[veto]->Fill(x, y);
        }
        void fill_or(double x, double y, std::vector<bool> is_vetoed, int evtNum, RunContext& ctx, std::ofstream* out_file) {
            h2_full_or->Fill(x, y);
            if (is_vetoed[0] || is_vetoed[1]) h2_passed_or->Fill(x, y);
            if (abs(x-cfg_.x_center) <= cfg_.xy_size_threshold && abs(y-cfg_.y_center) <= cfg_.xy_size_threshold) {
                h_full->Fill(2); // OR full
                h_full->Fill(3); // AND full
                if (is_vetoed[0] || is_vetoed[1]) h_passed->Fill(2); // OR passed
                if (is_vetoed[0] && is_vetoed[1]) h_passed->Fill(3); // AND passed
                if (!(is_vetoed[0] || is_vetoed[1])) {
                    if (cfg_.write_to_txt && out_file) {
                        (*out_file) << ctx.config.runNumber << " " << ctx.config.poolIndex << " " << evtNum << " " << x << " " << y << " 2\n";
                    }
                }
            }
        }
        void fill_nonveto(double x, double y, int veto, int evtNum, RunContext& ctx, bool is_vetoed, std::ofstream* out_file) {
            if (veto != 0 && veto != 1) {
                LOG_ERROR("VetoEffAlg::Impl::fill_nonveto: invalid veto value {}, expected 0 or 1", veto);
                return;
            }
            if (abs(x-cfg_.x_center) > cfg_.xy_outofrange_threshold || abs(y-cfg_.y_center) > cfg_.xy_outofrange_threshold) {
                h2_nonveto_full[veto]->Fill(x, y);
                h_nonveto_full->Fill(veto);
                if (is_vetoed) {
                    h2_nonveto_vetoed[veto]->Fill(x, y);
                    h_nonveto_vetoed->Fill(veto);
                    if (cfg_.write_to_txt && out_file) {
                        (*out_file) << ctx.config.runNumber << " " << ctx.config.poolIndex << " " << evtNum << " " << x << " " << y << " " << veto << "\n";
                    }
                } 
            }
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
                h2_nonveto_full[i]->Write();
                h2_nonveto_vetoed[i]->Write();
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
            std::vector<std::unique_ptr<TEfficiency>> nonveto_effs;
            nonveto_effs.reserve(h2_nonveto_full.size());
            for (size_t i = 0; i < h2_nonveto_full.size(); ++i) {
                auto eff = std::make_unique<TEfficiency>(*h2_nonveto_vetoed[i], *h2_nonveto_full[i]);
                // eff->SetDirectory(nullptr);
                nonveto_effs.push_back(std::move(eff));
            }
            for (size_t i = 0; i < nonveto_effs.size(); ++i) {
                nonveto_effs[i]->SetName(Form("eff_nonveto_input%zu", i+4));
                nonveto_effs[i]->SetTitle(Form("Faluty veto Efficiency for Input%zu;X [mm];Y [mm]", i+4));
                nonveto_effs[i]->Write();
                if (cfg_.write_to_png) {
                    TCanvas c(Form("c_eff_nonveto_input%zu", i+4), Form("Non-veto Efficiency for Input%zu", i+4), 800, 600);
                    nonveto_effs[i]->Draw("COLZ");
                    c.SaveAs(Form("%s/eff_nonveto_input%zu.png", cfg_.out_png_dir.c_str(), i+4));
                }
            }

            std::unique_ptr<TEfficiency> overall_eff = std::make_unique<TEfficiency>(*h_passed, *h_full);
            // overall_eff->SetDirectory(nullptr);
            overall_eff->SetName("overall_efficiency");
            overall_eff->SetTitle(Form("Overall Veto Efficiency within |x|<%.1f mm, |y|<%.1f mm;Veto;Efficiency", cfg_.xy_size_threshold, cfg_.xy_size_threshold));
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
                TH1D* h_dummy = new TH1D("h_dummy", overall_eff->GetTitle(), 4, -0.5, 3.5);
                h_dummy->GetXaxis()->SetBinLabel(1, "Input4");
                h_dummy->GetXaxis()->SetBinLabel(2, "Input5");
                h_dummy->GetXaxis()->SetBinLabel(3, "OR");
                h_dummy->GetXaxis()->SetBinLabel(4, "AND");
                h_dummy->GetXaxis()->SetTitle("Veto");
                h_dummy->GetYaxis()->SetTitle("Efficiency");
                h_dummy->SetDirectory(nullptr);
                h_dummy->Draw();
                h_dummy->GetYaxis()->SetRangeUser(cfg_.draw_range_y_min, 1.0);
                graph->Draw("Psame");
                c.SaveAs(Form("%s/overall_efficiency.png", cfg_.out_png_dir.c_str()));
            }

            std::unique_ptr<TEfficiency> overall_nonveto_eff = std::make_unique<TEfficiency>(*h_nonveto_vetoed, *h_nonveto_full);
            // overall_nonveto_eff->SetDirectory(nullptr);
            overall_nonveto_eff->SetName("overall_nonveto_efficiency");
            overall_nonveto_eff->SetTitle(Form("Overall Faulty Veto Efficiency for |x|>%.1f mm or |y|>%.1f mm;Veto;Efficiency", cfg_.xy_outofrange_threshold, cfg_.xy_outofrange_threshold));
            overall_nonveto_eff->Draw("AP");
            overall_nonveto_eff->Write();
            gPad->Update();
            auto graph_nonveto = overall_nonveto_eff->GetPaintedGraph();
            if (graph_nonveto) {
                graph_nonveto->GetXaxis()->SetBinLabel(1, "Input4");
                graph_nonveto->GetXaxis()->SetBinLabel(2, "Input5");
            }
            if (cfg_.write_to_png) {
                TCanvas c("c_overall_nonveto_efficiency", "Overall Faulty Veto Efficiency", 800, 600);
                TH1D* h_dummy = new TH1D("h_dummy_nonveto", overall_nonveto_eff->GetTitle(), 2, -0.5, 1.5);
                h_dummy->GetXaxis()->SetBinLabel(1, "Input4");
                h_dummy->GetXaxis()->SetBinLabel(2, "Input5");
                h_dummy->GetXaxis()->SetTitle("Veto");
                h_dummy->GetYaxis()->SetTitle("Efficiency");
                h_dummy->SetDirectory(nullptr);
                h_dummy->Draw();
                h_dummy->GetYaxis()->SetRangeUser(0.0, 0.1);
                graph_nonveto->Draw("Psame");
                c.SaveAs(Form("%s/overall_nonveto_efficiency.png", cfg_.out_png_dir.c_str()));
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
            
            h2_nonveto_full.clear();
            h2_nonveto_vetoed.clear();
            h2_nonveto_full.reserve(2);
            h2_nonveto_vetoed.reserve(2);
            h2_nonveto_full.push_back(std::make_unique<TH2D>("h2_nonveto_full_input4", "Non-veto full input distribution;x[mm];y[mm]", NBIN_XY, XYMIN*2, XYMAX*2, NBIN_XY, XYMIN*2, XYMAX*2));
            h2_nonveto_vetoed.push_back(std::make_unique<TH2D>("h2_nonveto_vetoed_input4", "Non-veto vetoed input distribution;x[mm];y[mm]", NBIN_XY, XYMIN*2, XYMAX*2, NBIN_XY, XYMIN*2, XYMAX*2));
            h2_nonveto_full[0]->SetDirectory(nullptr);
            h2_nonveto_vetoed[0]->SetDirectory(nullptr);
            h2_nonveto_full.push_back(std::make_unique<TH2D>("h2_nonveto_full_input5", "Non-veto full input distribution;x[mm];y[mm]", NBIN_XY, XYMIN*2, XYMAX*2, NBIN_XY, XYMIN*2, XYMAX*2));
            h2_nonveto_vetoed.push_back(std::make_unique<TH2D>("h2_nonveto_vetoed_input5", "Non-veto vetoed input distribution;x[mm];y[mm]", NBIN_XY, XYMIN*2, XYMAX*2, NBIN_XY, XYMIN*2, XYMAX*2));
            h2_nonveto_full[1]->SetDirectory(nullptr);
            h2_nonveto_vetoed[1]->SetDirectory(nullptr);

            h_full = std::make_unique<TH1D>("h_full", "Full distribution;Veto;Entries", 4, -0.5, 3.5);
            h_passed = std::make_unique<TH1D>("h_passed", "Passed distribution;Veto;Entries", 4, -0.5, 3.5);
            h_full->SetDirectory(nullptr);
            h_passed->SetDirectory(nullptr);
            h2_full_or = std::make_unique<TH2D>("h2_full_or", "Full input distribution (OR);x[mm];y[mm]", NBIN_XY, XYMIN, XYMAX, NBIN_XY, XYMIN, XYMAX);
            h2_passed_or = std::make_unique<TH2D>("h2_passed_or", "Passed input distribution (OR);x[mm];y[mm]", NBIN_XY, XYMIN, XYMAX, NBIN_XY, XYMIN, XYMAX);
            h2_full_or->SetDirectory(nullptr);
            h2_passed_or->SetDirectory(nullptr);
            h_nonveto_full = std::make_unique<TH1D>("h_nonveto_full", "Non-veto full distribution;Veto;Entries", 2, -0.5, 1.5);
            h_nonveto_vetoed = std::make_unique<TH1D>("h_nonveto_vetoed", "Non-veto vetoed distribution;Veto;Entries", 2, -0.5, 1.5);
            h_nonveto_full->SetDirectory(nullptr);
            h_nonveto_vetoed->SetDirectory(nullptr);
        }
        std::vector<std::unique_ptr<TH2D>> h2_fulls;
        std::vector<std::unique_ptr<TH2D>> h2_passeds;
        std::unique_ptr<TH2D> h2_full_or;
        std::unique_ptr<TH2D> h2_passed_or;
        std::vector<std::unique_ptr<TH2D>> h2_nonveto_full;
        std::vector<std::unique_ptr<TH2D>> h2_nonveto_vetoed;
        std::unique_ptr<TH1D> h_full;
        std::unique_ptr<TH1D> h_passed;
        std::unique_ptr<TH1D> h_nonveto_full;
        std::unique_ptr<TH1D> h_nonveto_vetoed;
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
                }else {
                    LOG_DEBUG("VetoEffAlg: no track selection cut applied");
                    if (track.valid!=1) {
                        LOG_DEBUG("VetoEffAlg: track is not valid, skipping");
                        std::cerr << "VetoEffAlg: track is not valid, skipping" << std::endl;
                        return;
                    }
                }
            }
            auto rawdata = evt.get<AHCALTLURawData>(cfg_.in_data_key);
            std::vector<bool> is_vetoed(2, false);
            for (int i = 0; i < 2; ++i){
                double x = track.init_pos_x + track.direction_x * (cfg_.z_pos_input4 + i*(cfg_.z_pos_input5-cfg_.z_pos_input4));
                double y = track.init_pos_y + track.direction_y * (cfg_.z_pos_input4 + i*(cfg_.z_pos_input5-cfg_.z_pos_input4));
                bool Input = (rawdata.Inputs[i+4] == 1);
                impl_->fill(x, y, i, evt.event_counter(), this->ctx(), Input, out_file_.get());
                is_vetoed[i] = Input;
                impl_->fill_nonveto(x, y, i, evt.event_counter(), this->ctx(), Input, out_file_.get());
            }
            double x = track.init_pos_x + track.direction_x * (cfg_.z_pos_input4 + 0.5*(cfg_.z_pos_input5-cfg_.z_pos_input4));
            double y = track.init_pos_y + track.direction_y * (cfg_.z_pos_input4 + 0.5*(cfg_.z_pos_input5-cfg_.z_pos_input4));
            impl_->fill_or(x, y, is_vetoed, evt.event_counter(), this->ctx(), out_file_.get());
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
            std::vector<bool> is_vetoed(2, false);
            for (int i = 0; i < 2; ++i){
                double x = track.x + track.tx * (cfg_.z_pos_input4 + i*(cfg_.z_pos_input5-cfg_.z_pos_input4)-track.z);
                double y = track.y + track.ty * (cfg_.z_pos_input4 + i*(cfg_.z_pos_input5-cfg_.z_pos_input4)-track.z);
                bool Input = (rawdata.Inputs[i+4] == 1);
                impl_->fill(x, y, i, evt.event_counter(), this->ctx(), Input, out_file_.get());
                is_vetoed[i] = Input;
                impl_->fill_nonveto(x, y, i, evt.event_counter(), this->ctx(), Input, out_file_.get());
            }
            double x = track.x + track.tx * (cfg_.z_pos_input4 + 0.5*(cfg_.z_pos_input5-cfg_.z_pos_input4)-track.z);
            double y = track.y + track.ty * (cfg_.z_pos_input4 + 0.5*(cfg_.z_pos_input5-cfg_.z_pos_input4)-track.z);
            impl_->fill_or(x, y, is_vetoed, evt.event_counter(), this->ctx(), out_file_.get());
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
        cfg_.write_to_txt = get_or<bool>(cfg, "write_to_txt", cfg_.write_to_txt);
        cfg_.out_txt_name = get_or<std::string>(cfg, "out_txt_name", cfg_.out_txt_name);
        cfg_.xy_size_threshold = get_or<double>(cfg, "xy_size_threshold", cfg_.xy_size_threshold);
        cfg_.xy_outofrange_threshold = get_or<double>(cfg, "xy_outofrange_threshold", cfg_.xy_outofrange_threshold);
        cfg_.x_center = get_or<double>(cfg, "x_center", cfg_.x_center);
        cfg_.y_center = get_or<double>(cfg, "y_center", cfg_.y_center);
        cfg_.z_pos_input4 = get_or<double>(cfg, "z_pos_input4", cfg_.z_pos_input4);
        cfg_.z_pos_input5 = get_or<double>(cfg, "z_pos_input5", cfg_.z_pos_input5);
        cfg_.threshold_input4 = get_or<double>(cfg, "threshold_input4", cfg_.threshold_input4);
        cfg_.threshold_input5 = get_or<double>(cfg, "threshold_input5", cfg_.threshold_input5);
        cfg_.draw_range_y_min = get_or<double>(cfg, "draw_range_y_min", cfg_.draw_range_y_min);
    }
    void VetoEffAlg::initialize() {
        if (cfg_.write_to_txt) {
            out_file_ = std::make_unique<std::ofstream>();
            out_file_->open(cfg_.out_txt_name, std::ios::out);
            if (!out_file_->is_open()) {
                LOG_ERROR("VetoEffAlg: cannot open output text file: {}", cfg_.out_txt_name);
                out_file_.reset();
            } else {
                LOG_INFO("VetoEffAlg: writing detailed event info to {}", cfg_.out_txt_name);
            }
        }
    }
} // namespace AHCALRecoAlg