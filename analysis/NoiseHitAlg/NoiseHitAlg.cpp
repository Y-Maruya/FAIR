#include "NoiseHitAlg.hpp"
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

AHCAL_REGISTER_ALG(AHCALRecoAlg::NoiseHitAlg, "NoiseHitAlg")
const double XYMIN = -AHCALGeometry::x_max;
const double XYMAX = +AHCALGeometry::x_max;
constexpr int NBIN_XY = 18;
namespace AHCALRecoAlg{

    struct NoiseHitAlg::Impl {
        explicit Impl(NoiseHitAlgCfg cfg) : cfg_(std::move(cfg)) {}
        static TDirectory* ensureDir(TDirectory* top, const char* name) {
            if (!top) return nullptr;
            auto* d = dynamic_cast<TDirectory*>(top->Get(name));
            if (!d) d = top->mkdir(name);
            return d;
        }
        void fill(const AHCALRawHit& h) {
            const int cellid = h.cellID;

            int hg = h.hg_adc;
            if (hg >= cfg_.xmin && hg <= cfg_.xmax) {
            auto& ptr = hg_hist_[cellid];
            if (!ptr) ptr = make_hist(cellid, /*isHG=*/true);
            ptr->Fill(hg);
            }

            int lg = h.lg_adc;
            if (lg >= cfg_.xmin && lg <= cfg_.xmax) {
            auto& ptr = lg_hist_[cellid];
            if (!ptr) ptr = make_hist(cellid, /*isHG=*/false);
            ptr->Fill(lg);
            }
        }
        void write() {
            if (written_) return;
            written_ = true;

            auto fout = std::unique_ptr<TFile>(TFile::Open(cfg_.out_filename.c_str(), "RECREATE"));
            if (!fout || fout->IsZombie()) {
                LOG_ERROR("NoiseHitAlg: cannot create output file: {}", cfg_.out_filename);
                return;
            }
            LOG_INFO("NoiseHitAlg: writing output to {}", cfg_.out_filename);
            fout->cd();
            TDirectory* dHistograms = ensureDir(fout.get(), "NoiseHistograms");
            for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer){
                TDirectory* dLayer = ensureDir(dHistograms, Form("Layer_%02d", layer));
                for (int chip = 0; chip < AHCALGeometry::chip_No; ++chip){
                    TDirectory* dChip = ensureDir(dLayer, Form("Chip_%02d", chip));
                    TDirectory* dHG = ensureDir(dChip, "HG");
                    TDirectory* dLG = ensureDir(dChip, "LG");
                    for (int channel = 0; channel < AHCALGeometry::channel_No; ++channel){
                        int cellid = layer*100000 + chip*10000 + channel;
                        auto it_hg = hg_hist_.find(cellid);
                        if (it_hg != hg_hist_.end()) {
                            if (dHG) dHG->cd();
                            it_hg->second->Write();
                        }
                        auto it_lg = lg_hist_.find(cellid);
                        if (it_lg != lg_hist_.end()) {
                            if (dLG) dLG->cd();
                            it_lg->second->Write();
                        }
                    }
                }
            }
            
            auto makeMap = [&](const char* base, int L, const char* title) {
                auto h = std::make_unique<TH2D>(
                    Form("%s_L%02d", base, L),
                    Form("%s L%d;X [mm];Y [mm]", title, L),
                    NBIN_XY, XYMIN, XYMAX,
                    NBIN_XY, XYMIN, XYMAX
                );
                h->SetDirectory(nullptr);
                return h;
            };

            fout->cd();
            TDirectory* dMaps = ensureDir(fout.get(), "NoiseMap");
             for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                auto h2 = makeMap("h2NoiseMap", L, "Noise Map: Enties");
                if (dMaps) {
                    dMaps->cd();
                    // h2->Write();
                }
            }
            h2NoiseRatio_->Write("h2NoiseRatio");
            h2NoiseRatio_Layer->Write("hNoiseRatio_Layer");
            h2NoiseRatio_->Scale(1.0/nEvents_);
            h2NoiseRatio_Layer->Scale(1.0/nEvents_);
            h2NoiseRatio_->Write("h2NoiseRatio_Normalized");
            h2NoiseRatio_Layer->Write("hNoiseRatio_Layer_Normalized");
            if (cfg_.write_to_png) {
                gStyle->SetOptStat(0);
                TCanvas c("c", "c", 800, 600);
                for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                    auto h2 = makeMap("h2NoiseMap", L, "Noise Map: Enties");
                    for (int chip = 0; chip < AHCALGeometry::chip_No; ++chip) {
                        for (int channel = 0; channel < AHCALGeometry::channel_No; ++channel) {
                            h2->Fill(
                                AHCALGeometry::Pos_X(channel, chip),
                                AHCALGeometry::Pos_Y(channel, chip),
                                h2NoiseRatio_->GetBinContent(L*9 + chip + 1, channel + 1)
                            );
                        }
                    }
                    dMaps->cd();
                    h2->Write(Form("h2NoiseMap_L%02d", L));
                    if (cfg_.write_to_png) {
                        c.cd();
                        h2->Draw("COLZ");
                        if(gSystem->AccessPathName(cfg_.out_png_dir.c_str())) {
                            gSystem->mkdir(cfg_.out_png_dir.c_str(), true);
                        }
                        c.SaveAs(Form("%s/NoiseMap_L%02d.png", cfg_.out_png_dir.c_str(), L));
                    }
                }
                if (cfg_.write_to_png) {
                    h2NoiseRatio_Layer->Draw();
                    c.SaveAs(Form("%s/NoiseRatio_Layer.png", cfg_.out_png_dir.c_str()));
                }
            }

             
            fout->Close();
            LOG_INFO("NoiseHitAlg: wrote {}", cfg_.out_filename);
        }

        std::unique_ptr<TH1D> make_hist(int cellid, bool isHG) {
            int layer = cellid/100000;
            int chip = cellid/10000 % 10;
            const std::string name = (isHG ? "hNoiseHG_" : "hNoiseLG_") + std::to_string(cellid);
            const std::string title = "Layer " + std::to_string(layer) + " chip "+ std::to_string(chip) + (isHG ? " Noise HG;ADC;counts" : " Noise LG;ADC;counts");
            auto h = std::make_unique<TH1D>(name.c_str(), title.c_str(), cfg_.nbin, cfg_.xmin, cfg_.xmax);
            h->SetDirectory(nullptr);
            return h;
        }
        NoiseHitAlgCfg cfg_;
        bool written_ = false;
        void initialize_histograms() {
            h2NoiseRatio_ = std::make_unique<TH2D>("h2NoiseRatio", "Noise Ratio;9*layer+chip;channel", AHCALGeometry::Layer_No*9, -0.5, AHCALGeometry::Layer_No*9-0.5, AHCALGeometry::channel_No, -0.5, AHCALGeometry::channel_No-0.5);
            h2NoiseRatio_Layer = std::make_unique<TH1D>("hNoiseRatio_Layer", "Noise Ratio per Layer;Layer;Noise Ratio", AHCALGeometry::Layer_No, -0.5, AHCALGeometry::Layer_No-0.5);
            h2NoiseRatio_->SetDirectory(nullptr);
            h2NoiseRatio_Layer->SetDirectory(nullptr);
        }
        std::unordered_map<int, std::unique_ptr<TH1D>> hg_hist_;
        std::unordered_map<int, std::unique_ptr<TH1D>> lg_hist_;
        std::unique_ptr<TH2D> h2NoiseRatio_;
        std::unique_ptr<TH1D> h2NoiseRatio_Layer;
        int nEvents_ = 0;
    };

    void NoiseHitAlg::ImplDeleter::operator()(Impl* p) const {
        delete p;
    }
    NoiseHitAlg::~NoiseHitAlg(){
        if (impl_) impl_->write();
    }

    void NoiseHitAlg::execute(EventStore& evt){
        if (!impl_){
            impl_.reset(new Impl(cfg_));
            impl_->initialize_histograms();
        }
        auto rawHits = evt.get<std::vector<AHCALRawHit>>(cfg_.in_raw_hit_key);
        // track exist
        if (cfg_.string_track_struct == "SimpleFittedTrack") {
            auto track = evt.get<SimpleFittedTrack>(cfg_.in_track_hit_key);
            // track not exist
            if (!track.valid) {
                
                int nHits = 0;
                int layer_count[AHCALGeometry::Layer_No] = {0};
                for (const auto& rh : rawHits) {
                    if (rh.hittag == 1) {
                        nHits++;
                        impl_->h2NoiseRatio_->Fill(rh.layer()*9 + rh.chip(), rh.channel());
                        if (layer_count[rh.layer()] == 0) impl_->h2NoiseRatio_Layer->Fill(rh.layer(), 1);
                        layer_count[rh.layer()]++;
                    }
                }
                if (cfg_.nHits_threshold > 0 && nHits <= cfg_.nHits_threshold){
                    for (const auto& rh : rawHits){
                        if (rh.hittag == 1) {
                            impl_->fill(rh);
                        }
                    }
                }
                impl_->nEvents_++;
                return;
            }
            if (cfg_.track_selection_string.empty()) {
                // LOG_ERROR("NoiseHitAlg: track_selection_string is empty, cannot select hits. Please provide a valid selection string in the configuration.");
                return;
            }
            SimpleFittedTrackCut hit_cut(cfg_.track_selection_string);
            if (!hit_cut.eval(track)) {
                LOG_DEBUG("NoiseHitAlg: track did not pass hit selection cut '{}'", cfg_.track_selection_string);
                return;
            }
            for (const auto& rh : rawHits) {
                if (rh.hittag == 1) {
                    if (std::find(track.outTrackHitsIndices.begin(), track.outTrackHitsIndices.end(), rh.index) != track.outTrackHitsIndices.end()) {
                        impl_->fill(rh);
                    }
                }
            }
            return;
        } else if (cfg_.string_track_struct == "Track") {
            auto track = evt.get<Track>(cfg_.in_track_hit_key);
            // track not exist
            if (!track.valid) {
                int nHits = 0;
                int layer_count[AHCALGeometry::Layer_No] = {0};
                for (const auto& rh : rawHits) {
                    if (rh.hittag == 1) {
                        nHits++;
                        impl_->h2NoiseRatio_->Fill(rh.layer()*9 + rh.chip(), rh.channel());
                        if (layer_count[rh.layer()] == 0) impl_->h2NoiseRatio_Layer->Fill(rh.layer(), 1);
                        layer_count[rh.layer()]++;
                    }
                }
                if (cfg_.nHits_threshold > 0 && nHits <= cfg_.nHits_threshold){
                    for (const auto& rh : rawHits) {
                        if (rh.hittag == 1) {
                            impl_->fill(rh);
                        }
                    }
                }
                impl_->nEvents_++;
                return;
            }
            if (cfg_.track_selection_string.empty()) {
                // LOG_ERROR("NoiseHitAlg: track_selection_string is empty, cannot select hits. Please provide a valid selection string in the configuration.");
                return;
            }
            TrackCut hit_cut(cfg_.track_selection_string);
            if (!hit_cut.eval(track)) {
                LOG_DEBUG("NoiseHitAlg: track did not pass hit selection cut '{}'", cfg_.track_selection_string);
                return;
            }
            for (const auto& rh : rawHits) {
                if (rh.hittag == 1) {
                    if (std::find(track.outTrackHitsIndices.begin(), track.outTrackHitsIndices.end(), rh.index) != track.outTrackHitsIndices.end()) {
                        impl_->fill(rh);
                    }
                }
            }
            return;
        } else {
            LOG_ERROR("NoiseHitAlg: unknown track struct string '{}'", cfg_.string_track_struct);
            return;
        }
    }
    void NoiseHitAlg::parse_cfg(const YAML::Node& cfg){
        cfg_.in_raw_hit_key = get_or<std::string>(cfg, "in_raw_hit_key", cfg_.in_raw_hit_key);
        cfg_.in_track_hit_key = get_or<std::string>(cfg, "in_track_hit_key", cfg_.in_track_hit_key);
        cfg_.string_track_struct = get_or<std::string>(cfg, "string_track_struct", cfg_.string_track_struct);
        cfg_.track_selection_string = get_or<std::string>(cfg, "track_selection_string", cfg_.track_selection_string);
        cfg_.nHits_threshold = get_or<int>(cfg, "nHits_threshold", cfg_.nHits_threshold);
        cfg_.out_filename = get_or<std::string>(cfg, "out_filename", cfg_.out_filename);
        cfg_.write_to_png = get_or<bool>(cfg, "write_to_png", cfg_.write_to_png);
        cfg_.write_to_txt = get_or<bool>(cfg, "write_to_txt", cfg_.write_to_txt);
        cfg_.out_txt_name = get_or<std::string>(cfg, "out_txt_name", cfg_.out_txt_name);
        cfg_.out_png_dir = get_or<std::string>(cfg, "out_png_dir", cfg_.out_png_dir);
        cfg_.xmin = get_or<double>(cfg, "xmin", cfg_.xmin);
        cfg_.xmax = get_or<double>(cfg, "xmax", cfg_.xmax);
        cfg_.nbin = get_or<int>(cfg, "nbin", cfg_.nbin);
    }
    void NoiseHitAlg::initialize() {
        if (cfg_.write_to_txt) {
            out_file_ = std::make_unique<std::ofstream>();
            out_file_->open(cfg_.out_txt_name, std::ios::out);
            if (!out_file_->is_open()) {
                LOG_ERROR("NoiseHitAlg: cannot open output text file: {}", cfg_.out_txt_name);
                out_file_.reset();
            } else {
                LOG_INFO("NoiseHitAlg: writing detailed event info to {}", cfg_.out_txt_name);
            }
        }
    }
} // namespace AHCALRecoAlg