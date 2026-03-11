#include "MuonEffAlg.hpp"
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

AHCAL_REGISTER_ALG(AHCALRecoAlg::MuonEffAlg, "MuonEffAlg")
const double XYMIN = -AHCALGeometry::x_max;
const double XYMAX = +AHCALGeometry::x_max;
constexpr int NBIN_XY = 18;
namespace AHCALRecoAlg{

    static TDirectory* ensureDir(TDirectory* top, const char* name) {
        if (!top) return nullptr;
        auto* d = dynamic_cast<TDirectory*>(top->Get(name));
        if (!d) d = top->mkdir(name);
        return d;
    }
    static void drawLayerLabel(int layer, double x=0.10, double y=0.92) {
        TLatex l;
        l.SetNDC(true);
        l.SetTextSize(0.10);
        l.DrawLatex(x, y, Form("L%d", layer));
    }

    static inline void cellid_to_xy(int chip, int channel, double& x, double& y) {
        // Geometry helper expects channel_ID [0..35], chip_ID [0..8]
        x = AHCALGeometry::Pos_X(channel, chip);
        y = AHCALGeometry::Pos_Y(channel, chip);
    }
    struct MuonEffAlg::Impl {
        explicit Impl(MuonEffAlgCfg cfg) : cfg_(std::move(cfg)) {}
    
        void fill_channel(int cellid, bool is_in_track = true) {
            int layer = cellid/100000;
            int chip = cellid/10000 % 10;
            int channel = cellid % 10000;
            int channel_id = layer*AHCALGeometry::chip_No*AHCALGeometry::channel_No + chip*AHCALGeometry::channel_No + channel;
            h_full_channel->Fill(channel_id);
            if (is_in_track) {
                h_passed_channel->Fill(channel_id);
            }
        } 
        void fill_chip(int cellid, bool is_in_track = true) {
            int layer = cellid/100000;
            int chip = cellid/10000 % 10;
            int chip_id = layer*AHCALGeometry::chip_No + chip;
            h_full_chip->Fill(chip_id);
            if (is_in_track) {
                h_passed_chip->Fill(chip_id);
            }
        }
        
        void write() {
            if (written_) return;
            written_ = true;

            auto fout = std::unique_ptr<TFile>(TFile::Open(cfg_.out_filename.c_str(), "RECREATE"));
            if (!fout || fout->IsZombie()) {
                LOG_ERROR("MuonEffAlg: cannot create output file: {}", cfg_.out_filename);
                return;
            }
            LOG_INFO("MuonEffAlg: writing output to {}", cfg_.out_filename);
            fout->cd();
            h_full_channel->Write();
            h_passed_channel->Write();
            h_full_chip->Write();
            h_passed_chip->Write();
            h_full_trigger_layer->Write();
            h_passed_Input_trigger_layer->Write();
            h_passed_Input_and_Hittag_trigger_layer->Write();
            LOG_DEBUG("MuonEffAlg: h_full_channel entries = {}, h_passed_channel entries = {}, h_full_chip entries = {}, h_passed_chip entries = {}, h_full_trigger_layer entries = {}, h_passed_Input_trigger_layer entries = {}, h_passed_Input_and_Hittag_trigger_layer entries = {}", 
                h_full_channel->GetEntries(), h_passed_channel->GetEntries(), h_full_chip->GetEntries(), h_passed_chip->GetEntries(), h_full_trigger_layer->GetEntries(), h_passed_Input_trigger_layer->GetEntries(), h_passed_Input_and_Hittag_trigger_layer->GetEntries());


            TEfficiency* eff_channel = new TEfficiency(*h_passed_channel, *h_full_channel);
            eff_channel->SetName("eff_channel");
            std::string title = "Muon Efficiency";
            if (cfg_.hittag_or_MIP_cut) {
                title += ": HitTag==1";
            } else {
                title += ": Nmip>" + std::to_string(cfg_.MIP_cut);
            }
            std::string titles = title + ";40*9*Layer+chip*9+channel;Efficiency";
            eff_channel->SetTitle(titles.c_str());
            eff_channel->Write();
            TEfficiency* eff_chip = new TEfficiency(*h_passed_chip, *h_full_chip);
            eff_chip->SetName("eff_chip");
            titles = title + ";9*Layer+chip;Efficiency";
            eff_chip->SetTitle(titles.c_str());
            eff_chip->Write();
            LOG_DEBUG("MuonEffAlg: calculated efficiencies for channels and chips");

            TDirectory* dCan = ensureDir(fout.get(), "Canvas");
            gStyle->SetOptStat(0);
            if (dCan) dCan->cd();
            TCanvas c_eff_channel("c_eff_channel", "Muon Efficiency by Channel", 8000, 600);
            c_eff_channel.SetRightMargin(0.15);
            eff_channel->Draw();
            c_eff_channel.Update();
            auto graph = eff_channel->GetPaintedGraph();
            graph->GetXaxis()->SetRangeUser(0, AHCALGeometry::Layer_No*AHCALGeometry::chip_No*AHCALGeometry::channel_No);
            graph->Draw("AP");
            for (int i_layer = 0; i_layer < AHCALGeometry::Layer_No; ++i_layer){
                TLine *line = new TLine(i_layer*AHCALGeometry::chip_No*AHCALGeometry::channel_No, 0, i_layer*AHCALGeometry::chip_No*AHCALGeometry::channel_No, 1);
                line->SetLineColor(kRed);
                line->SetLineStyle(2);
                line->Draw("same");
                TLatex *latex_layer = new TLatex();
                latex_layer->SetTextSize(0.03);
                latex_layer->SetTextFont(42);
                latex_layer->SetTextColor(kBlack);
                latex_layer->DrawLatex(i_layer*AHCALGeometry::chip_No*AHCALGeometry::channel_No+ AHCALGeometry::chip_No*AHCALGeometry::channel_No*0.3, 1.05, Form("Layer %d", i_layer));
            }
            c_eff_channel.Update();
            c_eff_channel.Modified();
            c_eff_channel.Write();
            if (cfg_.write_to_png) {
                if (gSystem->AccessPathName(cfg_.out_png_dir.c_str())) {
                    gSystem->mkdir(cfg_.out_png_dir.c_str(), true);
                }
                c_eff_channel.SaveAs(Form("%s/eff_channel.png", cfg_.out_png_dir.c_str()));
            }
            LOG_DEBUG("MuonEffAlg: created efficiency graph for channels");

            TCanvas c_eff_chip("c_eff_chip", "Muon Efficiency by Chip", 4000, 600);
            c_eff_chip.SetRightMargin(0.15);
            eff_chip->Draw();
            c_eff_chip.Update();
            auto graph_chip = eff_chip->GetPaintedGraph();
            titles = title + ";9*Layer+chip;Efficiency";
            graph_chip->SetTitle(titles.c_str());
            graph_chip->GetXaxis()->SetRangeUser(0, AHCALGeometry::Layer_No*AHCALGeometry::chip_No);
            graph_chip->Draw("AP");
            for (int i_layer = 0; i_layer < AHCALGeometry::Layer_No; ++i_layer){
                TLine *line = new TLine(i_layer*AHCALGeometry::chip_No, 0, i_layer*AHCALGeometry::chip_No, 1);
                line->SetLineColor(kRed);
                line->SetLineStyle(2);
                line->Draw("same");
                TLatex *latex_layer = new TLatex();
                latex_layer->SetTextSize(0.03);
                latex_layer->SetTextFont(42);
                latex_layer->SetTextColor(kBlack);
                latex_layer->DrawLatex(i_layer*AHCALGeometry::chip_No+ AHCALGeometry::chip_No*0.15, 1.05, Form("Layer %d", i_layer));
            }
            c_eff_chip.Update();
            c_eff_chip.Modified();
            c_eff_chip.Write();
            if (cfg_.write_to_png) {
                c_eff_chip.SaveAs(Form("%s/eff_chip.png", cfg_.out_png_dir.c_str()));
            }
            LOG_DEBUG("MuonEffAlg: created efficiency graph for chips");

            std::vector<std::unique_ptr<TH2D>> h2_chip_channel_xy_vec;
            TDirectory* dLayer = ensureDir(fout.get(), "Layer");
            for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                TDirectory* dL = ensureDir(dLayer, Form("Layer%02d", L));
                dL->cd();
                TCanvas c_layer(Form("c_layer%02d", L), Form("Muon Efficiency for Layer %02d", L), 800, 600);
                graph->GetXaxis()->SetRangeUser(L*AHCALGeometry::chip_No*AHCALGeometry::channel_No, (L+1)*AHCALGeometry::chip_No*AHCALGeometry::channel_No);
                titles = title + Form(": Layer %d;Chip*Channel;Efficiency", L);
                graph->SetTitle(titles.c_str());
                graph->Draw("AP");
                for (int i_chip = 0; i_chip < AHCALGeometry::chip_No; ++i_chip){
                    TLine *line = new TLine(L*AHCALGeometry::chip_No*AHCALGeometry::channel_No + i_chip*AHCALGeometry::channel_No, 0, L*AHCALGeometry::chip_No*AHCALGeometry::channel_No + i_chip*AHCALGeometry::channel_No, 1);
                    line->SetLineColor(kRed);
                    line->SetLineStyle(2);
                    line->Draw("same");
                    TLatex *latex_chip = new TLatex();
                    latex_chip->SetTextSize(0.03);
                    latex_chip->SetTextFont(42);
                    latex_chip->SetTextColor(kBlack);
                    latex_chip->DrawLatex(L*AHCALGeometry::chip_No*AHCALGeometry::channel_No+ i_chip*AHCALGeometry::channel_No+ AHCALGeometry::channel_No*0.2, 1.05, Form("Chip %d", i_chip));
                }
                c_layer.Update();
                c_layer.Modified();
                c_layer.Write();
                if (cfg_.write_to_png) {
                    c_layer.SaveAs(Form("%s/eff_layer%02d.png", cfg_.out_png_dir.c_str(), L));
                }
                std::unique_ptr<TH2D> h2_chip_channel_xy = std::make_unique<TH2D>(Form("h2_chip_channel_xy_layer%02d", L), Form("Muon Efficiency for Layer %02d;x [mm];y [mm]", L), NBIN_XY, XYMIN, XYMAX, NBIN_XY, XYMIN, XYMAX);
                h2_chip_channel_xy->SetDirectory(nullptr);
                h2_chip_channel_xy->GetZaxis()->SetRangeUser(0, 1);
                for (int i_chip = 0; i_chip < AHCALGeometry::chip_No; ++i_chip) {
                    for (int i_channel = 0; i_channel < AHCALGeometry::channel_No; ++i_channel) {
                        int channel_id = L*AHCALGeometry::chip_No*AHCALGeometry::channel_No + i_chip*AHCALGeometry::channel_No + i_channel;
                        double eff = eff_channel->GetEfficiency(channel_id+1);
                        h2_chip_channel_xy->Fill(AHCALGeometry::Pos_X(i_channel, i_chip), AHCALGeometry::Pos_Y(i_channel, i_chip), eff);
                    }
                }
                TCanvas c_xy(Form("c_xy_layer%02d", L), Form("Muon Efficiency for Layer %02d in XY plane", L), 600, 600);
                h2_chip_channel_xy->SetStats(0);
                h2_chip_channel_xy->Draw("COLZ");
                c_xy.Update();
                c_xy.Modified();
                c_xy.Write();
                if (cfg_.write_to_png) {
                    c_xy.SaveAs(Form("%s/eff_layer%02d_xy.png", cfg_.out_png_dir.c_str(), L));
                }
                h2_chip_channel_xy_vec.push_back(std::move(h2_chip_channel_xy));
            }
            LOG_DEBUG("MuonEffAlg: created efficiency graphs for each layer and XY maps");

            if (dCan) dCan->cd();
            auto cAllchannel = std::make_unique<TCanvas>("cAllEff_channel", "Muon Efficiency maps (all layers)", 5600, 4200);
            cAllchannel->Divide(7, 6, 0.001, 0.001);
            for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                cAllchannel->cd(L + 1);
                gPad->SetMargin(0.08, 0.14, 0.08, 0.10);
                h2_chip_channel_xy_vec[L]->Draw("COLZ");
                drawLayerLabel(L);
            }        
            cAllchannel->Update();
            cAllchannel->Modified();
            cAllchannel->Write();
            if (cfg_.write_to_png) {
                cAllchannel->SaveAs(Form("%s/eff_channel_all_layers.png", cfg_.out_png_dir.c_str()));
            }
            fout->cd();
            TEfficiency* eff_trigger_layer = new TEfficiency(*h_passed_Input_trigger_layer, *h_full_trigger_layer);
            TEfficiency* eff_trigger_layer_with_hittag = new TEfficiency(*h_passed_Input_and_Hittag_trigger_layer, *h_full_trigger_layer);
            TEfficiency* eff_trigger_layer_hittag_input1 = new TEfficiency(*h_passed_Input_and_Hittag_trigger_layer, *h_passed_Input_trigger_layer);
            eff_trigger_layer->SetName("eff_trigger_layer");
            eff_trigger_layer_with_hittag->SetName("eff_trigger_layer_with_hittag");
            eff_trigger_layer_hittag_input1->SetName("eff_trigger_layer_hittag_input1");
            eff_trigger_layer->SetTitle("Muon Efficiency by Trigger Layer;Input;Efficiency");
            eff_trigger_layer_with_hittag->SetTitle("Muon Efficiency by Trigger Layer with HitTag Cut;Input;Efficiency");
            eff_trigger_layer_hittag_input1->SetTitle("Muon Efficiency of HitTag == 1 among Input=1 events;Input;Efficiency");
            eff_trigger_layer->Write();
            eff_trigger_layer_with_hittag->Write();
            eff_trigger_layer_hittag_input1->Write();
            if (cfg_.write_to_png) {
                TCanvas c_trigger_layer("c_trigger_layer", "Muon Efficiency by Trigger Layer", 800, 600);
                eff_trigger_layer->Draw();
                c_trigger_layer.Update();
                c_trigger_layer.Modified();
                c_trigger_layer.SaveAs(Form("%s/eff_trigger_layer.png", cfg_.out_png_dir.c_str()));

                TCanvas c_trigger_layer_hittag("c_trigger_layer_hittag", "Muon Efficiency by Trigger Layer with HitTag Cut", 800, 600);
                eff_trigger_layer_with_hittag->Draw();
                c_trigger_layer_hittag.Update();
                c_trigger_layer_hittag.Modified();
                c_trigger_layer_hittag.SaveAs(Form("%s/eff_trigger_layer_with_hittag.png", cfg_.out_png_dir.c_str()));

                TCanvas c_trigger_layer_hittag_input1("c_trigger_layer_hittag_input1", "Muon Efficiency of HitTag == 1 among Input=1 events", 800, 600);
                eff_trigger_layer_hittag_input1->Draw();
                c_trigger_layer_hittag_input1.Update();
                c_trigger_layer_hittag_input1.Modified();
                c_trigger_layer_hittag_input1.SaveAs(Form("%s/eff_trigger_layer_hittag_input1.png", cfg_.out_png_dir.c_str()));
            }
            fout->Close();
            LOG_INFO("MuonEffAlg: wrote {}", cfg_.out_filename);
            LOG_INFO("MuonEffAlg: total entries in channel histograms: {}", h_full_channel->GetEntries());
        }

        MuonEffAlgCfg cfg_;
        bool written_ = false;
        void initialize_histograms() {
            h_full_channel = std::make_unique<TH1D>("h_full_channel", "Entries;40*9*Layer+chip*9+channel;Entries", AHCALGeometry::Layer_No*AHCALGeometry::chip_No*AHCALGeometry::channel_No, 0, AHCALGeometry::Layer_No*AHCALGeometry::chip_No*AHCALGeometry::channel_No);
            h_passed_channel = std::make_unique<TH1D>("h_passed_channel", "Track Passed;40*9*Layer+chip*9+channel;Entries", AHCALGeometry::Layer_No*AHCALGeometry::chip_No*AHCALGeometry::channel_No, 0, AHCALGeometry::Layer_No*AHCALGeometry::chip_No*AHCALGeometry::channel_No);
            h_full_chip = std::make_unique<TH1D>("h_full_chip", "Entries;9*Layer+chip;Entries", AHCALGeometry::Layer_No*AHCALGeometry::chip_No, 0, AHCALGeometry::Layer_No*AHCALGeometry::chip_No);
            h_passed_chip = std::make_unique<TH1D>("h_passed_chip", "Track Passed;9*Layer+chip;Entries", AHCALGeometry::Layer_No*AHCALGeometry::chip_No, 0, AHCALGeometry::Layer_No*AHCALGeometry::chip_No);
            h_full_channel->SetDirectory(nullptr);
            h_passed_channel->SetDirectory(nullptr);
            h_full_chip->SetDirectory(nullptr);
            h_passed_chip->SetDirectory(nullptr);
            //
            h_full_trigger_layer = std::make_unique<TH1D>("h_full_trigger_layer", "Track Passed;Input;Entries", cfg_.ntrigger_layer, 0-0.5, cfg_.ntrigger_layer-0.5);
            h_passed_Input_trigger_layer = std::make_unique<TH1D>("h_passed_Input_trigger_layer", "Track Passed and Input=1;Input;Entries", cfg_.ntrigger_layer, 0-0.5, cfg_.ntrigger_layer-0.5);
            h_passed_Input_and_Hittag_trigger_layer = std::make_unique<TH1D>("h_passed_Input_and_Hittag_trigger_layer", "Track Passed and Input=1 and Hittag=1;Input;Entries", cfg_.ntrigger_layer, 0-0.5, cfg_.ntrigger_layer-0.5);
            h_full_trigger_layer->SetDirectory(nullptr);
            h_passed_Input_trigger_layer->SetDirectory(nullptr);
            h_passed_Input_and_Hittag_trigger_layer->SetDirectory(nullptr);
        }
        std::unique_ptr<TH1D> h_full_channel;
        std::unique_ptr<TH1D> h_passed_channel;
        std::unique_ptr<TH1D> h_full_chip;
        std::unique_ptr<TH1D> h_passed_chip;
        std::unique_ptr<TH1D> h_full_trigger_layer;
        std::unique_ptr<TH1D> h_passed_Input_trigger_layer;
        std::unique_ptr<TH1D> h_passed_Input_and_Hittag_trigger_layer;
    };
    bool MuonEffAlg::is_hittag1_exist(std::vector<AHCALRawHit>& rawHits, int layer) {
        for (const auto& hit : rawHits) {
            if (hit.layer() == layer && hit.hittag == 1) {
                return true;
            }
        }
        return false;
    } 
    void MuonEffAlg::ImplDeleter::operator()(Impl* p) const {
        delete p;
    }
    MuonEffAlg::~MuonEffAlg(){
        if (impl_) impl_->write();
        if (out_file_) {
            out_file_->close();
            LOG_INFO("MuonEffAlg: wrote {}", cfg_.out_txt_filename.c_str());
        }
    }

    void MuonEffAlg::execute(EventStore& evt){
        if (!impl_){
            impl_.reset(new Impl(cfg_));
            impl_->initialize_histograms();
        }
        auto rawData = evt.get<AHCALTLURawData>(cfg_.in_rawdata_key);
        if (cfg_.string_track_struct == "SimpleFittedTrack") {
            auto track = evt.get<SimpleFittedTrack>(cfg_.in_track_key);
            if (!cfg_.track_selection_string.empty()) {
                SimpleFittedTrackCut cut(cfg_.track_selection_string);
                if (!cut.eval(track)) {
                    LOG_DEBUG("MuonEffAlg: track did not pass selection cut '{}'", cfg_.track_selection_string);
                    return;
                }
            }

            // loop over layers
            if (cfg_.hittag_or_MIP_cut){
                // use hittag
                auto rawhits = evt.get<std::vector<AHCALRawHit>>(cfg_.in_hit_key);
                for (int i_layer = 0; i_layer < AHCALGeometry::Layer_No; ++i_layer) {
                    double x = track.init_pos_x + track.direction_x * (AHCALGeometry::Pos_Z(i_layer));
                    double y = track.init_pos_y + track.direction_y * (AHCALGeometry::Pos_Z(i_layer));
                    // double z = AHCALGeometry::Pos_Z(i_layer);
                    int chip = -1;
                    int channel = -1;
                    AHCALGeometry::inverse(x, y, chip, channel);
                    if (abs(x) > AHCALGeometry::x_max || abs(y) > AHCALGeometry::y_max) {
                        LOG_DEBUG("MuonEffAlg: track extrapolation out of bounds for layer {}: x={}, y={}", i_layer, x, y);
                        continue;
                    }
                    bool is_in_track = false;
                    bool is_in_chip = false;
                    for (const auto& rh : rawhits) {
                        if (rh.layer() != i_layer) continue;
                        if (rh.hittag == 1) {
                            double rh_x = AHCALGeometry::Pos_X(rh.channel(), rh.chip());
                            double rh_y = AHCALGeometry::Pos_Y(rh.channel(), rh.chip());
                            if (rh.chip() == chip) {
                                is_in_chip = true;
                            }
                            if (abs(x-rh_x) <= cfg_.xy_size_threshold && abs(y-rh_y) <= cfg_.xy_size_threshold) {
                                is_in_track = true;
                                chip = rh.chip();
                                channel = rh.channel();
                                break;
                            }
                        }
                    }
                    if (i_layer ==21 && chip == 1 && is_in_chip ==false){
                        if (cfg_.write_to_txt && out_file_) {
                            (*out_file_.get()) << this->ctx().config.runNumber << " " << this->ctx().config.poolIndex << " " << evt.event_counter() << " " << i_layer << " " << chip << " " << channel << " " << is_in_track << " " << is_in_chip << std::endl;
                        }
                    }
                    int cellid = i_layer*100000 + chip*10000 + channel;
                    impl_->fill_channel(cellid, is_in_track);
                    impl_->fill_chip(cellid, is_in_chip);
                    if (std::find(cfg_.trigger_layer.begin(), cfg_.trigger_layer.end(), i_layer) != cfg_.trigger_layer.end()) {
                        int i_input = std::distance(cfg_.trigger_layer.begin(), std::find(cfg_.trigger_layer.begin(), cfg_.trigger_layer.end(), i_layer));
                        impl_->h_full_trigger_layer->Fill(i_input);
                        if ( rawData.Inputs[i_input] == 1) {
                            impl_->h_passed_Input_trigger_layer->Fill(i_input);
                            if (is_hittag1_exist(rawhits, i_layer)) {
                                impl_->h_passed_Input_and_Hittag_trigger_layer->Fill(i_input);
                            }
                        }
                    }
                }
            } else {
                // use MIP cut
                auto recohits = evt.get<std::vector<AHCALRecoHit>>(cfg_.in_hit_key);
                for (int i_layer = 0; i_layer < AHCALGeometry::Layer_No; ++i_layer) {
                    double x = track.init_pos_x + track.direction_x * (AHCALGeometry::Pos_Z(i_layer));
                    double y = track.init_pos_y + track.direction_y * (AHCALGeometry::Pos_Z(i_layer));
                    // double z = AHCALGeometry::Pos_Z(i_layer);
                    int chip = -1;
                    int channel = -1;
                    AHCALGeometry::inverse(x, y, chip, channel);
                    if (abs(x) > AHCALGeometry::x_max || abs(y) > AHCALGeometry::y_max) {
                        LOG_DEBUG("MuonEffAlg: track extrapolation out of bounds for layer {}: x={}, y={}", i_layer, x, y);
                        continue;
                    }
                    bool is_in_track = false;
                    bool is_in_chip = false;
                    for (const auto& rh : recohits) {
                        if (rh.layer() != i_layer) continue;
                        if (rh.Nmip > cfg_.MIP_cut) {
                            double rh_x = AHCALGeometry::Pos_X(rh.channel(), rh.chip());
                            double rh_y = AHCALGeometry::Pos_Y(rh.channel(), rh.chip());
                            if (rh.chip() == chip) {
                                is_in_chip = true;
                            }
                            if (abs(x-rh_x) <= cfg_.xy_size_threshold && abs(y-rh_y) <= cfg_.xy_size_threshold) {
                                is_in_track = true;
                                chip = rh.chip();
                                channel = rh.channel();
                                break;
                            }
                        }
                    }
                    int cellid = i_layer*100000 + chip*10000 + channel;
                    impl_->fill_channel(cellid, is_in_track);
                    impl_->fill_chip(cellid, is_in_chip);
                }
            }
        } else if (cfg_.string_track_struct == "Track") {
            auto track = evt.get<Track>(cfg_.in_track_key);
            if (!cfg_.track_selection_string.empty()) {
                TrackCut cut(cfg_.track_selection_string);
                if (!cut.eval(track)) {
                    LOG_DEBUG("MuonEffAlg: track did not pass selection cut '{}'", cfg_.track_selection_string);
                    return;
                }
            }
            // similar loop as SimpleFittedTrack, but need to extrapolate track to each layer using track parameters
            // loop over layers
            if (cfg_.hittag_or_MIP_cut){
                // use hittag
                auto rawhits = evt.get<std::vector<AHCALRawHit>>(cfg_.in_hit_key);
                for (int i_layer = 0; i_layer < AHCALGeometry::Layer_No; ++i_layer) {
                    double x = track.x + track.tx * (AHCALGeometry::Pos_Z(i_layer)-track.z);
                    double y = track.y + track.ty * (AHCALGeometry::Pos_Z(i_layer)-track.z);
                    // double z = AHCALGeometry::Pos_Z(i_layer);
                    int chip = -1;
                    int channel = -1;
                    AHCALGeometry::inverse(x, y, chip, channel);
                    if (abs(x) > AHCALGeometry::x_max || abs(y) > AHCALGeometry::y_max) {
                        LOG_DEBUG("MuonEffAlg: track extrapolation out of bounds for layer {}: x={}, y={}", i_layer, x, y);
                        continue;
                    }
                    bool is_in_track = false;
                    bool is_in_chip = false;
                    for (const auto& rh : rawhits) {
                        if (rh.layer() != i_layer) continue;
                        if (rh.hittag == 1) {
                            double rh_x = AHCALGeometry::Pos_X(rh.channel(), rh.chip());
                            double rh_y = AHCALGeometry::Pos_Y(rh.channel(), rh.chip());
                            if (rh.chip() == chip) {
                                is_in_chip = true;
                            }
                            if (abs(x-rh_x) <= cfg_.xy_size_threshold && abs(y-rh_y) <= cfg_.xy_size_threshold) {
                                is_in_track = true;
                                chip = rh.chip();
                                channel = rh.channel();
                                break;
                            }
                        }
                    }
                    int cellid = i_layer*100000 + chip*10000 + channel;
                    impl_->fill_channel(cellid, is_in_track);
                    impl_->fill_chip(cellid, is_in_chip);
                    if (std::find(cfg_.trigger_layer.begin(), cfg_.trigger_layer.end(), i_layer) != cfg_.trigger_layer.end()) {
                        int i_input = std::distance(cfg_.trigger_layer.begin(), std::find(cfg_.trigger_layer.begin(), cfg_.trigger_layer.end(), i_layer));
                        impl_->h_full_trigger_layer->Fill(i_input);
                        if ( rawData.Inputs[i_input] == 1) {
                            impl_->h_passed_Input_trigger_layer->Fill(i_input);
                            if (is_hittag1_exist(rawhits, i_layer)) {
                                impl_->h_passed_Input_and_Hittag_trigger_layer->Fill(i_input);
                            }
                        }
                    }
                }
            } else {
                // use MIP cut
                auto recohits = evt.get<std::vector<AHCALRecoHit>>(cfg_.in_hit_key);
                for (int i_layer = 0; i_layer < AHCALGeometry::Layer_No; ++i_layer) {
                    double x = track.x + track.tx * (AHCALGeometry::Pos_Z(i_layer)-track.z);
                    double y = track.y + track.ty * (AHCALGeometry::Pos_Z(i_layer)-track.z);
                    // double z = AHCALGeometry::Pos_Z(i_layer);
                    int chip = -1;
                    int channel = -1;
                    AHCALGeometry::inverse(x, y, chip, channel);
                    if (abs(x) > AHCALGeometry::x_max || abs(y) > AHCALGeometry::y_max) {
                        LOG_DEBUG("MuonEffAlg: track extrapolation out of bounds for layer {}: x={}, y={}", i_layer, x, y);
                        continue;
                    }
                    bool is_in_track = false;
                    bool is_in_chip = false;
                    for (const auto& rh : recohits) {
                        if (rh.layer() != i_layer) continue;
                        if (rh.Nmip > cfg_.MIP_cut) {
                            double rh_x = AHCALGeometry::Pos_X(rh.channel(), rh.chip());
                            double rh_y = AHCALGeometry::Pos_Y(rh.channel(), rh.chip());
                            if (rh.chip() == chip) {
                                is_in_chip = true;
                            }
                            if (abs(x-rh_x) <= cfg_.xy_size_threshold && abs(y-rh_y) <= cfg_.xy_size_threshold) {
                                is_in_track = true;
                                chip = rh.chip();
                                channel = rh.channel();
                                break;
                            }
                        }
                    }
                    int cellid = i_layer*100000 + chip*10000 + channel;
                    impl_->fill_channel(cellid, is_in_track);
                    impl_->fill_chip(cellid, is_in_chip);
                }
            }
        } else {
            LOG_ERROR("MuonEffAlg: unknown track struct string '{}'", cfg_.string_track_struct);
            return;
        }
    }
    void MuonEffAlg::parse_cfg(const YAML::Node& cfg){
        cfg_.in_hit_key = get_or<std::string>(cfg, "in_hit_key", cfg_.in_hit_key);
        cfg_.in_track_key = get_or<std::string>(cfg, "in_track_key", cfg_.in_track_key);
        cfg_.in_rawdata_key = get_or<std::string>(cfg, "in_rawdata_key", cfg_.in_rawdata_key);
        cfg_.string_track_struct = get_or<std::string>(cfg, "string_track_struct", cfg_.string_track_struct);
        cfg_.track_selection_string = get_or<std::string>(cfg, "track_selection_string", cfg_.track_selection_string);
        cfg_.out_filename = get_or<std::string>(cfg, "out_filename", cfg_.out_filename);
        cfg_.write_to_png = get_or<bool>(cfg, "write_to_png", cfg_.write_to_png);
        cfg_.out_png_dir = get_or<std::string>(cfg, "out_png_dir", cfg_.out_png_dir);
        cfg_.write_to_txt = get_or<bool>(cfg, "write_to_txt", cfg_.write_to_txt);
        cfg_.out_txt_filename = get_or<std::string>(cfg, "out_txt_filename", cfg_.out_txt_filename);
        cfg_.hittag_or_MIP_cut = get_or<bool>(cfg, "hittag_or_MIP_cut", cfg_.hittag_or_MIP_cut);
        cfg_.MIP_cut = get_or<double>(cfg, "MIP_cut", cfg_.MIP_cut);
        cfg_.xy_size_threshold = get_or<double>(cfg, "xy_size_threshold", cfg_.xy_size_threshold);
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
    void MuonEffAlg::initialize(){
        if (cfg_.write_to_txt) {
            out_file_ = std::make_unique<std::ofstream>();
            out_file_->open(cfg_.out_txt_filename, std::ios::out);
            if (!out_file_->is_open()) {
                LOG_ERROR("VetoEffAlg: cannot open output text file: {}", cfg_.out_txt_filename);
                out_file_.reset();
            } else {
                LOG_INFO("VetoEffAlg: writing detailed event info to {}", cfg_.out_txt_filename);
            }
        }
    }
} // namespace AHCALRecoAlg