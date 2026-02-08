#include "MIPAlg.hpp"
#include "langaus.h"
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
#include <TF1.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

AHCAL_REGISTER_ALG(AHCALRecoAlg::MIPAlg, "MIPAlg")
const double XYMIN = -AHCALGeometry::x_max;
const double XYMAX = +AHCALGeometry::x_max;
constexpr int NBIN_XY = 18;
namespace AHCALRecoAlg{

    struct FitOut {
        double mpv = -1.0;
        double width = -1.0;
        double total_area = -1.0;
        double gaus_sigma = -1.0;
        double max_x = -1.0;
        double FWHM = -1.0;
        int entries = 0;
        double chi2 = -1.0;
        int ndf = 0;
        int fit_status = 999;
        bool fit_ok = false;
    };
    static FitOut fitLandauGaus(TH1D* h, 
                                int minEntries)
    {
        FitOut r;
        if (!h) return r;
        if (h->GetEntries() < minEntries){
            r.fit_status = -2;
            r.fit_ok = false;
            return r;
        }

        double fr[2];
        double sv[4], pllo[4], plhi[4], fps[4], fpe[4];
        double chisqr;
        int ndf;
        fr[0] = 0;fr[1] = 1500;
        pllo[0] = 0;  pllo[1] = 0; pllo[2] = 100; pllo[3] = 0;
        plhi[0] = 200;plhi[1] = 700;plhi[2] = 100000;plhi[3] = 300;
        sv[0] = 50; sv[1] = 300; sv[2] = 30000;   sv[3] = 60;

        TF1* fLandauGaus = langaufit(h, fr, sv, pllo, plhi, fps, fpe, &chisqr, &ndf);
        if (!fLandauGaus) return r;
        double maxx=0,fwhm;
        langaupro(fps,maxx,fwhm);
        r.mpv = fps[1];
        r.width = fps[0];
        r.total_area = fps[2];
        r.gaus_sigma = fps[3];
        r.entries = h->GetEntries();
        r.chi2 = chisqr;
        r.ndf = ndf;
        double x = h->FindBin(maxx);
        double chi2_ndf = (ndf > 0) ? (chisqr / double(ndf)) : -1.0;
        if (chi2_ndf > 30 || x > fr[1] - 40 || r.gaus_sigma < 30 || r.gaus_sigma > 160 || r.width < 15 || r.width > 80)
        {
            r.fit_status = -3;
            r.fit_ok = false;
        } else {
            r.fit_status = 1;
            r.fit_ok = true;
        }
        return r;
    }

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
    struct MIPAlg::Impl {
        explicit Impl(MIPAlgCfg cfg) : cfg_(std::move(cfg)) {}
    
        void fill(const AHCALRawHit& h) {
            const int cellid = h.cellID;

            int hg = h.hg_adc;
            if (hg >= cfg_.xmin && hg <= cfg_.xmax) {
                auto& ptr = hg_hist_[cellid];
                if (!ptr) ptr = make_hist(cellid);
                ptr->Fill(hg);
            }
        } 
        
        void write() {
            if (!cfg_.mip_to_file) return;
            if (written_) return;
            written_ = true;

            auto fout = std::unique_ptr<TFile>(TFile::Open(cfg_.out_mip_filename.c_str(), "RECREATE"));
            if (!fout || fout->IsZombie()) {
                LOG_ERROR("MIPAlg: cannot create output file: {}", cfg_.out_mip_filename);
                return;
            }

            // directories
            TDirectory* dHist = ensureDir(fout.get(), "MIP");
            for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                for (int C = 0; C < 10; ++C) {
                    if (dHist) dHist->mkdir(Form("MIP/Layer%02d/Chip%d", L, C));
                }
            }
            TDirectory* dMap = ensureDir(fout.get(), "Map");
            TDirectory* dCan = ensureDir(fout.get(), "Canvases");

            // union of keys 
            std::unordered_set<int> keys;
            keys.reserve(hg_hist_.size());
            for (const auto& [k, _] : hg_hist_) keys.insert(k);

            std::vector<std::unique_ptr<TH2D>> hMPV(AHCALGeometry::Layer_No), hEntries(AHCALGeometry::Layer_No), hWidth(AHCALGeometry::Layer_No), hTotal(AHCALGeometry::Layer_No), hGausSigma(AHCALGeometry::Layer_No);
            
            auto makeMap = [&](const char* base, int L, const char* title) {
                auto h = std::make_unique<TH2D>(
                    Form("h%s_L%02d", base, L), 
                    Form("%s L%d;X[mm];Y[mm]", title, L),
                    NBIN_XY, XYMIN, XYMAX,
                    NBIN_XY, XYMIN, XYMAX
                );
                h->SetDirectory(nullptr);
                return h;
            };

            for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                hMPV[L] = makeMap("MPV", L, "MIP MPV");
                hWidth[L] = makeMap("Width", L, "Landau Width");
                hEntries[L] = makeMap("Entries", L, "MIP Entries");
                hTotal[L] = makeMap("TotalArea", L, "MIP Total Area");
                hGausSigma[L] = makeMap("GausSigma", L, "MIP Gaussian Sigma");
            }

            TTree tp("mip", "MIP fit results");

            int cellid=-1;
            double mpv=-1.0, width=-1.0, total_area=-1.0, gaus_sigma=-1.0;
            int entries=0;
            double chi2=-1.0;
            int ndf=0;
            double chi2perndf=-1.0;
            int fit_status=999;
            int fit_ok = 0;
            double max_x=-1.0, FWHM=-1.0;
            double x_mm=-999, y_mm=-999;
            tp.Branch("cellid", &cellid);
            tp.Branch("MPV", &mpv);
            tp.Branch("width", &width);
            tp.Branch("TotalArea", &total_area);
            tp.Branch("entries", &entries);
            tp.Branch("gaus_sigma", &gaus_sigma);
            tp.Branch("max_x", &max_x);
            tp.Branch("FWHM", &FWHM);
            tp.Branch("chi2", &chi2);
            tp.Branch("ndf", &ndf);
            tp.Branch("chi2perndf", &chi2perndf);
            tp.Branch("fit_status", &fit_status);
            tp.Branch("fit_ok", &fit_ok);
            tp.Branch("x_mm", &x_mm);
            tp.Branch("y_mm", &y_mm);

            int nFitAll = 0, nFitOK = 0;
            LOG_INFO("MIPAlg: start fitting {} histograms", keys.size());
            int i = 0;
            for (int cid : keys) {
                cellid = cid;
                i++;
                if (i % 1000 == 0) {
                    LOG_INFO("MIPAlg: fitting histogram {}/{}", i, keys.size());
                }
                AHCALRawHit tmp;
                tmp.cellID = cellid;
                const int L  = tmp.layer();
                const int C  = tmp.chip();
                const int ch = tmp.channel();

                cellid_to_xy(C, ch, x_mm, y_mm);

                TH1D* hhg = nullptr;
                if (auto it = hg_hist_.find(cellid); it != hg_hist_.end()) hhg = it->second.get();

                FitOut fr;
                entries = hhg ? static_cast<int>(hhg->GetEntries()) : 0;
                if (hhg) {
                    nFitAll++;
                    fr = fitLandauGaus(hhg, cfg_.min_entries);
                    if (fr.fit_ok) nFitOK++;
                }
                mpv = fr.mpv;
                width = fr.width;
                total_area = fr.total_area;
                gaus_sigma = fr.gaus_sigma;
                chi2 = fr.chi2;
                ndf = fr.ndf;
                max_x = fr.max_x;
                FWHM = fr.FWHM;
                if (ndf > 0) chi2perndf = chi2 / double(ndf);
                fit_status = fr.fit_status;
                fit_ok = fr.fit_ok ? 1 : 0;

                if (L >= 0 && L < AHCALGeometry::Layer_No){
                    const int binx = hMPV[L]->GetXaxis()->FindBin(x_mm);
                    const int biny = hMPV[L]->GetYaxis()->FindBin(y_mm);
                    if (binx > 0 && binx <= NBIN_XY && biny > 0 && biny <= NBIN_XY){
                        if (fit_ok) {
                            hMPV[L]->SetBinContent(binx, biny, mpv);
                            hWidth[L]->SetBinContent(binx, biny, width);
                            hEntries[L]->SetBinContent(binx, biny, entries);
                            hTotal[L]->SetBinContent(binx, biny, total_area);
                            hGausSigma[L]->SetBinContent(binx, biny, gaus_sigma);
                        }
                    }
                }

                tp.Fill();
            }

            auto cAllMPV = std::make_unique<TCanvas>("cAllMPV_7x6", "MIP MPV maps (all layers)", 5600, 4200);
            cAllMPV->Divide(7, 6, 0.001, 0.001);

            for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                cAllMPV->cd(L + 1);
                gPad->SetMargin(0.08, 0.14, 0.08, 0.10);
                hMPV[L]->Draw("COLZ");
                drawLayerLabel(L);
            }
            if (dHist){
                for (auto& [cellid, h] : hg_hist_) {
                    fout->cd(Form("MIP/Layer%02d/Chip%d", cellid/100000, (cellid/10000)%10));
                    if (h) h->Write();
                }
            }
            if (dMap) dMap->cd();
            for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                hMPV[L]->Write();
                hWidth[L]->Write();
                hEntries[L]->Write();
                hTotal[L]->Write();
                hGausSigma[L]->Write();
            }
            if (dCan) dCan->cd();
            cAllMPV->Write();

            fout->cd();
            tp.Write();
            fout->Close();

            LOG_INFO("MIPAlg: wrote {}", cfg_.out_mip_filename);
            LOG_INFO("MIPAlg: fit {}/{} histograms successfully", nFitOK, nFitAll);
        }

        std::unique_ptr<TH1D> make_hist(int cellid){
            int layer = cellid/100000;
            int chip = cellid/10000 % 10;
            const std::string name = "hMIP_" + std::to_string(cellid);
            const std::string title = "Layer " + std::to_string(layer) + " chip "+ std::to_string(chip) + " MIP;ADC;counts";
            auto h = std::make_unique<TH1D>(name.c_str(), title.c_str(), cfg_.nbin, cfg_.xmin, cfg_.xmax);
            h->SetDirectory(nullptr);
            return h;
        }
        MIPAlgCfg cfg_;
        bool written_ = false;
        std::unordered_map<int, std::unique_ptr<TH1D>> hg_hist_;
    };

    void MIPAlg::ImplDeleter::operator()(Impl* p) const {
        delete p;
    }
    MIPAlg::~MIPAlg(){
        if (impl_) impl_->write();
    }

    void MIPAlg::execute(EventStore& evt){
        if (!impl_) impl_.reset(new Impl(cfg_));
        auto rawhits = evt.get<std::vector<AHCALRawHit>>(cfg_.in_rawhit_key);
        if (cfg_.string_track_struct == "SimpleFittedTrack") {
            auto track = evt.get<SimpleFittedTrack>(cfg_.in_track_key);
            if (!cfg_.track_selection_string.empty()) {
                SimpleFittedTrackCut cut(cfg_.track_selection_string);
                if (!cut.eval(track)) {
                    LOG_DEBUG("MIPAlg: track did not pass selection cut '{}'", cfg_.track_selection_string);
                    return;
                }
            }
            if (!track.inTrackHitsIndices.size()) {
                LOG_WARN("MIPAlg: track has no associated hits.");
                return;
            }
            for (const int index : track.inTrackHitsIndices) {
                if (index < 0 || index >= static_cast<int>(rawhits.size())) {
                    LOG_WARN("MIPAlg: track hit index {} out of range [0, {})", index, rawhits.size());
                    continue;
                }
                const auto& rh = rawhits.at(index);
                if (rh.index != index) {
                    LOG_WARN("MIPAlg: rawhit index {} does not match track hit index {}", rh.index, index);
                    continue;
                }
                impl_->fill(rh);
            }    
        } else if (cfg_.string_track_struct == "Track") {
            auto track = evt.get<Track>(cfg_.in_track_key);
            if (!cfg_.track_selection_string.empty()) {
                TrackCut cut(cfg_.track_selection_string);
                if (!cut.eval(track)) {
                    LOG_DEBUG("MIPAlg: track did not pass selection cut '{}'", cfg_.track_selection_string);
                    return;
                }
            }
            if (!track.inTrackHitsIndices.size()) {
                LOG_WARN("MIPAlg: no tracks found in the event.");
                return;
            }
            for (const int index : track.inTrackHitsIndices) {
                if (index < 0 || index >= static_cast<int>(rawhits.size())) {
                    LOG_WARN("MIPAlg: track hit index {} out of range [0, {})", index, rawhits.size());
                    continue;
                }
                const auto& rh = rawhits.at(index);
                if (rh.index != index) {
                    LOG_WARN("MIPAlg: rawhit index {} does not match track hit index {}", rh.index, index);
                    continue;
                }
                impl_->fill(rh);
            }
        } else {
            LOG_ERROR("MIPAlg: unknown track struct string '{}'", cfg_.string_track_struct);
            return;
        }
    }
    void MIPAlg::parse_cfg(const YAML::Node& cfg){
        cfg_.in_rawhit_key = get_or<std::string>(cfg, "in_rawhit_key", cfg_.in_rawhit_key);
        cfg_.in_track_key = get_or<std::string>(cfg, "in_track_key", cfg_.in_track_key);
        cfg_.string_track_struct = get_or<std::string>(cfg, "string_track_struct", cfg_.string_track_struct);
        cfg_.track_selection_string = get_or<std::string>(cfg, "track_selection_string", cfg_.track_selection_string);
        cfg_.mip_to_file = get_or<bool>(cfg, "mip_to_file", cfg_.mip_to_file);
        cfg_.mip_to_DB = get_or<bool>(cfg, "mip_to_DB", cfg_.mip_to_DB);
        cfg_.out_mip_filename = get_or<std::string>(cfg, "out_mip_filename", cfg_.out_mip_filename);
    }
} // namespace AHCALRecoAlg