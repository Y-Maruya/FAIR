#include "SPEAlg.hpp"
#include "fft_helper.h"
#include "langaus.h"
#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"
#include "CalibDBIO/PedestalReader/PedestalReader.hpp"
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
#include <TGraphErrors.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

AHCAL_REGISTER_ALG(AHCALRecoAlg::SPEAlg, "SPEAlg")

namespace AHCALRecoAlg {

    // ============ Fit output struct ============
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

    // ============ Landau-Gauss fitting ============
    static FitOut fitLandauGaus(TH1D* h, int minEntries, bool calculate_fwhm = true) {
        FitOut r;
        if (!h) return r;
        if (h->GetEntries() < minEntries) {
            r.fit_status = -2;
            r.fit_ok = false;
            return r;
        }

        double fr[2];
        double sv[4], pllo[4], plhi[4], fps[4], fpe[4];
        double chisqr;
        int ndf;
        fr[0] = 0;
        fr[1] = 1500;
        pllo[0] = 0;
        pllo[1] = 0;
        pllo[2] = 100;
        pllo[3] = 0;
        plhi[0] = 300;
        plhi[1] = 1200;
        plhi[2] = 100000;
        plhi[3] = 300;
        sv[0] = 50;
        sv[1] = 300;
        sv[2] = 30000;
        sv[3] = 60;

        TF1* fLandauGaus = langaufit(h, fr, sv, pllo, plhi, fps, fpe, &chisqr, &ndf);
        if (!fLandauGaus) return r;

        double maxx = 0, fwhm = -1.0;
        if (calculate_fwhm) {
            langaupro(fps, maxx, fwhm);
        }

        r.mpv = fps[1];
        r.width = fps[0];
        r.total_area = fps[2];
        r.gaus_sigma = fps[3];
        r.max_x = calculate_fwhm ? maxx : -1.0;
        r.FWHM = fwhm;
        r.entries = h->GetEntries();
        r.chi2 = chisqr;
        r.ndf = ndf;

        double x = calculate_fwhm ? h->FindBin(maxx) : -1.0;
        double chi2_ndf = (ndf > 0) ? (chisqr / double(ndf)) : -1.0;
        if (chi2_ndf > 30 || (x > 0 && x > fr[1] - 40) || r.gaus_sigma < 30 || 
            r.gaus_sigma > 160 || r.width < 15 || r.width > 80) {
            r.fit_status = -3;
            r.fit_ok = false;
        } else {
            r.fit_status = 1;
            r.fit_ok = true;
        }
        return r;
    }

    // ============ Directory utilities ============
    static TDirectory* ensureDir(TDirectory* top, const char* name) {
        if (!top) return nullptr;
        auto* d = dynamic_cast<TDirectory*>(top->Get(name));
        if (!d) d = top->mkdir(name);
        return d;
    }

    static void drawLayerLabel(int layer, double x = 0.10, double y = 0.92) {
        TLatex l;
        l.SetNDC(true);
        l.SetTextSize(0.10);
        l.DrawLatex(x, y, Form("L%d", layer));
    }

    static inline void cellid_to_xy(int chip, int channel, double& x, double& y) {
        x = AHCALGeometry::Pos_X(channel, chip);
        y = AHCALGeometry::Pos_Y(channel, chip);
    }

    static std::string formatTimestampFromRunContexts(const std::vector<RunContext>& run_contexts) {
        double start_time = 0.0;
        double end_time = 0.0;
        if (!run_contexts.empty()) {
            start_time = run_contexts.front().conditions.starttime;
            end_time = run_contexts.back().conditions.endtime;
        }

        double ref_time = (start_time > 0.0 && end_time > 0.0)
            ? (start_time + end_time) / 2.0
            : std::time(nullptr);

        std::time_t utc_time = static_cast<std::time_t>(ref_time);
        std::tm tm = {};
        gmtime_r(&utc_time, &tm);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
        return oss.str();
    }

    // ============ Implementation ============
    struct SPEAlg::Impl {
        struct FitResult {
            int cellid = -1;
            int layer = -1;
            int chip = -1;
            int channel = -1;
            int channel_index = -1;
            double x_mm = -999.0;
            double y_mm = -999.0;
            int entries = 0;
            FitOut fit_out;
            
            // FFT analysis results
            int npad = 0;
            double gain = -1.0;
            double gainErr = -1.0;
            double snr = -1.0;
            int kp = -1;
            double k0 = -1.0;
            double peakAmp = -1.0;
            double noiseMed = -1.0;
        };

        explicit Impl(SPEAlgCfg cfg, const RunContext& ctx)
            : cfg_(std::move(cfg)), ctx_(ctx) {
            loadPedestals();
            run_contexts_.push_back(ctx);
        }

        void loadPedestals() {
            ped_map_.clear();
            if (cfg_.read_pedestal_from_ROOT) {
                loadPedestals_fromROOT();
            } else if (cfg_.read_pedestal_from_DB) {
                loadPedestals_fromDB();
            } else {
                LOG_WARN("SPEAlg: no pedestal source specified");
            }
        }

        void loadPedestals_fromROOT() {
            std::unique_ptr<TFile> fped(TFile::Open(cfg_.in_pedestal_file.c_str(), "READ"));
            if (!fped || fped->IsZombie()) {
                LOG_ERROR("SPEAlg: cannot open pedestal file: {}", cfg_.in_pedestal_file);
                return;
            }

            TTree* tp = static_cast<TTree*>(fped->Get("pedestal"));
            if (!tp) {
                LOG_ERROR("SPEAlg: TTree 'pedestal' not found in {}", cfg_.in_pedestal_file);
                return;
            }

            int p_cellid = 0;
            double p_hgped = 0.0;
            double p_lgped = 0.0;
            tp->SetBranchAddress("cellid", &p_cellid);
            tp->SetBranchAddress("highgain_peak", &p_hgped);

            bool has_lg_ped = false;
            if (tp->GetBranch("lowgain_peak")) {
                tp->SetBranchAddress("lowgain_peak", &p_lgped);
                has_lg_ped = true;
            }

            ped_map_.reserve(tp->GetEntries());
            for (Long64_t i = 0; i < tp->GetEntries(); ++i) {
                tp->GetEntry(i);
                CalibDBIO::Pedestal ped{};
                ped.HighGainPeak = p_hgped;
                ped.LowGainPeak = has_lg_ped ? p_lgped : 0.0;
                ped_map_[p_cellid] = ped;
            }

            LOG_INFO("SPEAlg: loaded {} pedestal entries from ROOT", ped_map_.size());
        }

        void loadPedestals_fromDB() {
            CalibDBIO::PedestalReader reader(ctx_.config.runNumber);
            ped_map_ = reader.getPedestalMap();
            LOG_INFO("SPEAlg: loaded {} pedestal entries from DB", ped_map_.size());
        }

        void change_run(const RunContext& new_ctx) {
            ctx_ = new_ctx;
            run_contexts_.push_back(new_ctx);
        }

        void fill(const AHCALRawHit& h) {
            const int cellid = h.cellID;

            double hg_value = static_cast<double>(h.hg_adc);
            if (cfg_.substrate_pedestal) {
                auto itp = ped_map_.find(cellid);
                if (itp == ped_map_.end()) {
                    n_missing_ped_++;
                    return;
                }
                hg_value -= itp->second.HighGainPeak;
            }

            if (hg_value >= cfg_.xmin && hg_value <= cfg_.xmax) {
                auto& ptr = hg_hist_[cellid];
                if (!ptr) {
                    ptr = make_hist(cellid);
                }
                ptr->Fill(hg_value);
            }
        }

        void buildFitCache() {
            if (cache_built_) return;
            cache_built_ = true;

            fit_cache_.clear();
            fit_cache_.reserve(hg_hist_.size());

            nFitAll_ = 0;
            nFitOK_ = 0;

            for (const auto& [cid, hist_ptr] : hg_hist_) {
                FitResult res;
                res.cellid = cid;
                res.layer = cid / 100000;
                res.chip = (cid / 10000) % 10;
                res.channel = cid % 10000;
                res.channel_index = res.layer * (AHCALGeometry::chip_No * AHCALGeometry::channel_No) +
                                    res.chip * AHCALGeometry::channel_No + res.channel;
                cellid_to_xy(res.chip, res.channel, res.x_mm, res.y_mm);

                res.entries = hist_ptr ? static_cast<int>(hist_ptr->GetEntries()) : 0;
                if (hist_ptr && cfg_.fit) {
                    nFitAll_++;
                    res.fit_out = fitLandauGaus(hist_ptr.get(), cfg_.min_entries, cfg_.calculate_fwhm);
                    if (res.fit_out.fit_ok) {
                        nFitOK_++;
                    }
                }

                fit_cache_.emplace(cid, std::move(res));
            }

            LOG_INFO("SPEAlg: built fit cache with {}/{} successful fits", nFitOK_, nFitAll_);
        }

        void doFFTAnalysis() {
            if (fft_done_) return;
            fft_done_ = true;

            int count = 0;
            for (auto& [cid, res] : fit_cache_) {
                auto it = hg_hist_.find(cid);
                if (it == hg_hist_.end()) continue;
                TH1D* h = it->second.get();
                if (!h || h->GetEntries() < cfg_.min_ent_fft) continue;

                // Get histogram data
                int bmax = h->GetNbinsX();
                std::vector<double> v(bmax, 0.0);
                for (int b = 1; b <= bmax; ++b) {
                    v[b - 1] = h->GetBinContent(b);
                }

                // Mean-center
                double mu = 0;
                for (double x : v) mu += x;
                if (!v.empty()) mu /= (double)v.size();
                for (double& x : v) x -= mu;

                // Compute FFT
                int Npad = chooseNpad((int)v.size(), cfg_.npad_req);
                std::unique_ptr<TH1D> hFFT(doFFT(v, Npad, Form("hFFT_L%d_C%d_ch%d_k", 
                                                                 res.layer, res.chip, res.channel)));
                if (!hFFT) continue;

                // Pick peak and extract gain/SNR
                KPick kp = pickK0(hFFT.get(), Npad, cfg_.gain_min, cfg_.gain_max,
                                 cfg_.edge_gain_tol, cfg_.kedge_skip);

                res.npad = Npad;
                res.gain = kp.gain;
                res.gainErr = kp.gainErr;
                res.snr = kp.snr;
                res.kp = kp.kp;
                res.k0 = kp.k0;
                res.peakAmp = kp.peakAmp;
                res.noiseMed = kp.noiseMed;

                count++;
            }
            LOG_INFO("SPEAlg: FFT analysis completed for {} channels", count);
        }

        void fillMapsFromCache(std::vector<std::unique_ptr<TH2D>>& hMPV,
                                std::vector<std::unique_ptr<TH2D>>& hWidth,
                                std::vector<std::unique_ptr<TH2D>>& hEntries,
                                std::vector<std::unique_ptr<TH2D>>& hGain,
                                std::vector<std::unique_ptr<TH2D>>& hSNR) {
            for (const auto& [cid, res] : fit_cache_) {
                if (res.fit_out.fit_ok && res.layer >= 0 && res.layer < AHCALGeometry::Layer_No) {
                    const int binx = hMPV[res.layer]->GetXaxis()->FindBin(res.x_mm);
                    const int biny = hMPV[res.layer]->GetYaxis()->FindBin(res.y_mm);
                    if (binx > 0 && binx <= 18 && biny > 0 && biny <= 18) {
                        hMPV[res.layer]->SetBinContent(binx, biny, res.fit_out.mpv);
                        hWidth[res.layer]->SetBinContent(binx, biny, res.fit_out.width);
                        hEntries[res.layer]->SetBinContent(binx, biny, res.entries);
                    }
                }
                if (res.snr >= cfg_.snr_threshold && res.layer >= 0 && res.layer < AHCALGeometry::Layer_No) {
                    const int binx = hGain[res.layer]->GetXaxis()->FindBin(res.x_mm);
                    const int biny = hGain[res.layer]->GetYaxis()->FindBin(res.y_mm);
                    if (binx > 0 && binx <= 18 && biny > 0 && biny <= 18) {
                        hGain[res.layer]->SetBinContent(binx, biny, res.gain);
                        hSNR[res.layer]->SetBinContent(binx, biny, res.snr);
                    }
                }
            }
        }

        void writeRootOutput() {
            if (!cfg_.mip_to_file) return;

            auto fout = std::unique_ptr<TFile>(TFile::Open(cfg_.out_mip_filename.c_str(), "RECREATE"));
            if (!fout || fout->IsZombie()) {
                LOG_ERROR("SPEAlg: cannot create output file: {}", cfg_.out_mip_filename);
                return;
            }

            // Directories
            TDirectory* dMap = ensureDir(fout.get(), "Map");
            TDirectory* dCan = ensureDir(fout.get(), "Canvases");
            TDirectory* dHist = nullptr;
            if (cfg_.save_per_channel_hists) {
                dHist = ensureDir(fout.get(), "SPE");
                for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                    if (dHist) dHist->mkdir(Form("Layer%02d", L));
                    for (int C = 0; C < AHCALGeometry::chip_No; ++C) {
                        if (dHist) dHist->mkdir(Form("Layer%02d/Chip%d", L, C));
                    }
                }
            }

            // Build fit cache if needed
            if (!cache_built_) buildFitCache();
            if (!fft_done_) doFFTAnalysis();

            // Create output maps
            const double XYMIN = -AHCALGeometry::x_max;
            const double XYMAX = +AHCALGeometry::x_max;
            const int NBIN_XY = 18;

            std::vector<std::unique_ptr<TH2D>> hMPV(AHCALGeometry::Layer_No);
            std::vector<std::unique_ptr<TH2D>> hWidth(AHCALGeometry::Layer_No);
            std::vector<std::unique_ptr<TH2D>> hEntries(AHCALGeometry::Layer_No);
            std::vector<std::unique_ptr<TH2D>> hGain(AHCALGeometry::Layer_No);
            std::vector<std::unique_ptr<TH2D>> hSNR(AHCALGeometry::Layer_No);

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
                hGain[L] = makeMap("Gain", L, "FFT Gain");
                hSNR[L] = makeMap("SNR", L, "FFT SNR");
            }

            // Fill maps from cache
            fillMapsFromCache(hMPV, hWidth, hEntries, hGain, hSNR);

            // Build output tree
            TTree tp("spe", "SPE (Signal-to-noise ratio, Pedestal subtraction, Efficiency) analysis");
            int cellid = -1;
            int layer = -1, chip = -1, channel = -1;
            double mpv = -1.0, width = -1.0, total_area = -1.0, gaus_sigma = -1.0;
            int entries = 0;
            double chi2 = -1.0;
            int ndf = 0;
            double chi2perndf = -1.0;
            int fit_status = 999;
            int fit_ok = 0;
            double max_x = -1.0, FWHM = -1.0;
            double x_mm = -999, y_mm = -999;
            int npad = 0;
            double gain = -1.0, gainErr = -1.0, snr = -1.0;
            int kp = -1;
            double k0 = -1.0;
            double peakAmp = -1.0, noiseMed = -1.0;

            tp.Branch("cellid", &cellid);
            tp.Branch("layer", &layer);
            tp.Branch("chip", &chip);
            tp.Branch("channel", &channel);
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
            tp.Branch("npad", &npad);
            tp.Branch("gain", &gain);
            tp.Branch("gainErr", &gainErr);
            tp.Branch("snr", &snr);
            tp.Branch("kp", &kp);
            tp.Branch("k0", &k0);
            tp.Branch("peakAmp", &peakAmp);
            tp.Branch("noiseMed", &noiseMed);

            // Optional MuonSelect-compatible calibration tree
            TTree tMipCal("MIP_Calibration", "MIP calibration (MuonSelect compatible)");
            int calCellID = 0;
            double calMPV = 0;
            double calChi2 = 0;
            tMipCal.Branch("CellID", &calCellID, "CellID/I");
            tMipCal.Branch("MPV", &calMPV, "MPV/D");
            tMipCal.Branch("chi2perndf", &calChi2, "chi2perndf/D");

            // Layer efficiency stats
            int denLayer[AHCALGeometry::Layer_No] = {0};
            int numLayer[AHCALGeometry::Layer_No] = {0};
            double bestSNR = -1;
            int bestCellID = -1;

            int saved_hists = 0;

            for (const auto& [cid, res] : fit_cache_) {
                cellid = res.cellid;
                layer = res.layer;
                chip = res.chip;
                channel = res.channel;
                mpv = res.fit_out.mpv;
                width = res.fit_out.width;
                total_area = res.fit_out.total_area;
                gaus_sigma = res.fit_out.gaus_sigma;
                chi2 = res.fit_out.chi2;
                ndf = res.fit_out.ndf;
                max_x = res.fit_out.max_x;
                FWHM = res.fit_out.FWHM;
                if (ndf > 0) chi2perndf = res.fit_out.chi2 / double(ndf);
                fit_status = res.fit_out.fit_status;
                fit_ok = res.fit_out.fit_ok ? 1 : 0;
                entries = res.entries;
                x_mm = res.x_mm;
                y_mm = res.y_mm;
                npad = res.npad;
                gain = res.gain;
                gainErr = res.gainErr;
                snr = res.snr;
                kp = res.kp;
                k0 = res.k0;
                peakAmp = res.peakAmp;
                noiseMed = res.noiseMed;

                tp.Fill();

                // MIP calibration tree
                if (res.fit_out.fit_ok) {
                    calCellID = res.cellid;
                    calMPV = res.fit_out.mpv;
                    calChi2 = (ndf > 0) ? res.fit_out.chi2 / double(ndf) : -1.0;
                    tMipCal.Fill();
                }

                // Layer efficiency counts
                if (layer >= 0 && layer < AHCALGeometry::Layer_No) {
                    denLayer[layer]++;
                    if (res.snr >= cfg_.snr_threshold) {
                        numLayer[layer]++;
                    }
                }

                // Track best SNR channel for examples
                if (res.snr > bestSNR) {
                    bestSNR = res.snr;
                    bestCellID = res.cellid;
                }

                // Save per-channel histograms with FFT data
                if (cfg_.save_per_channel_hists && dHist && saved_hists < cfg_.save_max_channel_hists) {
                    auto it = hg_hist_.find(cid);
                    if (it != hg_hist_.end() && it->second) {
                        // Create nested directory structure
                        TDirectory* dNested = dynamic_cast<TDirectory*>(dHist->Get(Form("Layer%02d/Chip%d", layer, chip)));
                        if (!dNested) dNested = dHist;
                        
                        TH1D* hMIP = it->second.get();
                        if (dNested) dNested->cd();
                        hMIP->Write(Form("hMIP_L%d_C%d_ch%d", layer, chip, channel));
                        
                        // Compute and save FFT histograms if SNR is good
                        if (res.snr >= cfg_.snr_threshold && res.snr > 0) {
                            int bmax = hMIP->GetNbinsX();
                            if (res.fit_out.mpv > 0) {
                                int b = hMIP->FindBin(res.fit_out.mpv);
                                if (b >= 1) bmax = b;
                            }
                            
                            std::vector<double> v(bmax, 0.0);
                            for (int b = 1; b <= bmax; ++b) {
                                v[b - 1] = hMIP->GetBinContent(b);
                            }
                            
                            // Remove DC
                            double mu = 0;
                            for (double x : v) mu += x;
                            if (!v.empty()) mu /= (double)v.size();
                            for (double& x : v) x -= mu;
                            
                            // FFT analysis
                            std::unique_ptr<TH1D> hFFT_k(doFFT(v, res.npad, Form("hFFT_k_L%d_C%d_ch%d", layer, chip, channel)));
                            if (hFFT_k) {
                                hFFT_k->Write();
                                
                                // Convert to gain space
                                std::unique_ptr<TH1D> hFFT_gain(new TH1D(
                                    Form("hFFT_gain_L%d_C%d_ch%d", layer, chip, channel),
                                    "FFT spectrum vs gain;gain [ADC/p.e.];amplitude",
                                    250, cfg_.gain_min, cfg_.gain_max));
                                hFFT_gain->SetDirectory(nullptr);
                                
                                for (int k = 1; k <= hFFT_k->GetNbinsX(); ++k) {
                                    double a = hFFT_k->GetBinContent(k);
                                    if (k <= 0) continue;
                                    double g = (double)res.npad / (double)k;
                                    if (g < cfg_.gain_min || g > cfg_.gain_max) continue;
                                    int bg = hFFT_gain->FindBin(g);
                                    hFFT_gain->SetBinContent(bg, hFFT_gain->GetBinContent(bg) + a);
                                }
                                hFFT_gain->Write();
                            }
                        }
                        
                        saved_hists++;
                    }
                }
            }

            // Layer efficiency graph
            std::vector<double> xL, yEff, eEff;
            for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                if (denLayer[L] <= 0) continue;
                double eff = (double)numLayer[L] / (double)denLayer[L];
                double err = sqrt(eff * (1.0 - eff) / denLayer[L]);
                xL.push_back(L);
                yEff.push_back(eff);
                eEff.push_back(err);
            }

            TGraphErrors* gLayerEff = new TGraphErrors((int)xL.size());
            gLayerEff->SetName("gLayerEff");
            gLayerEff->SetTitle(Form("Layer efficiency (SNR>%.1f);Layer;Efficiency", cfg_.snr_threshold));
            for (int i = 0; i < (int)xL.size(); ++i) {
                gLayerEff->SetPoint(i, xL[i], yEff[i]);
                gLayerEff->SetPointError(i, 0.0, eEff[i]);
            }

            // Save best-SNR channel examples with full FFT analysis
            if (cfg_.save_examples && bestCellID >= 0) {
                auto it_mip = hg_hist_.find(bestCellID);
                auto it_fit = fit_cache_.find(bestCellID);
                if (it_mip != hg_hist_.end() && it_mip->second && it_fit != fit_cache_.end()) {
                    TH1D* hMIP = it_mip->second.get();
                    const FitResult& bestRes = it_fit->second;
                    
                    // Create examples subdirectory
                    auto dExamples = ensureDir(fout.get(), "examples/best_channel");
                    
                    dExamples->cd();
                    hMIP->Write("hMIP_best");
                    
                    // Compute and save FFT histograms for best channel
                    int bmax = hMIP->GetNbinsX();
                    if (bestRes.fit_out.mpv > 0) {
                        int b = hMIP->FindBin(bestRes.fit_out.mpv);
                        if (b >= 1) bmax = b;
                    }
                    
                    std::vector<double> v(bmax, 0.0);
                    for (int b = 1; b <= bmax; ++b) {
                        v[b - 1] = hMIP->GetBinContent(b);
                    }
                    
                    // Remove DC
                    double mu = 0;
                    for (double x : v) mu += x;
                    if (!v.empty()) mu /= (double)v.size();
                    for (double& x : v) x -= mu;
                    
                    // FFT analysis
                    std::unique_ptr<TH1D> hFFT_k(doFFT(v, bestRes.npad, "hFFT_best_k"));
                    if (hFFT_k) {
                        hFFT_k->Write();
                        
                        // Convert to gain space
                        std::unique_ptr<TH1D> hFFT_gain(new TH1D(
                            "hFFT_best_gain",
                            "FFT spectrum vs gain (best channel);gain [ADC/p.e.];amplitude",
                            250, cfg_.gain_min, cfg_.gain_max));
                        hFFT_gain->SetDirectory(nullptr);
                        
                        for (int k = 1; k <= hFFT_k->GetNbinsX(); ++k) {
                            double a = hFFT_k->GetBinContent(k);
                            if (k <= 0) continue;
                            double g = (double)bestRes.npad / (double)k;
                            if (g < cfg_.gain_min || g > cfg_.gain_max) continue;
                            int bg = hFFT_gain->FindBin(g);
                            hFFT_gain->SetBinContent(bg, hFFT_gain->GetBinContent(bg) + a);
                        }
                        hFFT_gain->Write();
                    }
                    
                    // Save PNG if enabled
                    if (cfg_.output_to_png) {
                        std::filesystem::path png_dir(cfg_.out_png_dir);
                        std::error_code ec;
                        std::filesystem::create_directories(png_dir, ec);
                        
                        auto c_mip = std::make_unique<TCanvas>("c_mip", "c_mip", 900, 650);
                        hMIP->Draw("E");
                        c_mip->SaveAs(Form("%s/h_MIP_best.png", cfg_.out_png_dir.c_str()));
                        
                        if (hFFT_k) {
                            auto c_fft = std::make_unique<TCanvas>("c_fft", "c_fft", 900, 650);
                            hFFT_k->Draw();
                            c_fft->SaveAs(Form("%s/h_FFT_best_k.png", cfg_.out_png_dir.c_str()));
                        }
                    }
                    
                    // Write canvas to ROOT file
                    auto can = std::make_unique<TCanvas>("cBestExample", "Best SNR channel", 800, 600);
                    hMIP->Draw("E");
                    if (dCan) dCan->cd();
                    can->Write();
                }
            }

            // Write canvases
            auto cAllMPV = std::make_unique<TCanvas>("cAllMPV_7x6", "MIP MPV maps (all layers)", 5600, 4200);
            cAllMPV->Divide(7, 6, 0.001, 0.001);
            for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                cAllMPV->cd(L + 1);
                gPad->SetMargin(0.08, 0.14, 0.08, 0.10);
                hMPV[L]->Draw("COLZ");
                drawLayerLabel(L);
            }

            if (cfg_.output_to_png) {
                std::filesystem::path png_dir(cfg_.out_png_dir);
                std::error_code ec;
                std::filesystem::create_directories(png_dir, ec);
                TString png_path = Form("%s/LayerEff_SNRgt%.1f.png", cfg_.out_png_dir.c_str(), cfg_.snr_threshold);
                auto can_eff = std::make_unique<TCanvas>("cLayerEff", "Layer Efficiency", 1000, 600);
                gLayerEff->Draw("AP");
                can_eff->SaveAs(png_path);
            }

            // Write to file
            if (dMap) dMap->cd();
            for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                hMPV[L]->Write();
                hWidth[L]->Write();
                hEntries[L]->Write();
                hGain[L]->Write();
                hSNR[L]->Write();
            }

            if (dCan) dCan->cd();
            cAllMPV->Write();
            gLayerEff->Write();

            fout->cd();
            tp.Write();
            tMipCal.Write();
            fout->Close();

            LOG_INFO("SPEAlg: wrote {} ({} per-channel hists)", cfg_.out_mip_filename, saved_hists);
        }

        void writeJsonOutput() {
            if (!cfg_.mip_to_json) return;
            if (!cache_built_) buildFitCache();

            const int n_channels_per_layer = AHCALGeometry::chip_No * AHCALGeometry::channel_No;
            std::filesystem::path out_dir(cfg_.out_json_dirname.empty() ? "." : cfg_.out_json_dirname);
            std::error_code ec;
            std::filesystem::create_directories(out_dir, ec);
            if (ec) {
                LOG_ERROR("SPEAlg: cannot create output directory {}: {}", out_dir.string(), ec.message());
                return;
            }

            const std::string timestamp = formatTimestampFromRunContexts(run_contexts_);
            std::vector<int> same_data_runs;
            for (const auto& rc : run_contexts_) {
                same_data_runs.push_back(rc.config.runNumber);
            }

            // Per run, per layer JSON (MIPAlg format)
            for (const auto& rc : run_contexts_) {
                for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer) {
                    // Initialize per-channel arrays
                    std::vector<double> mpv_arr(n_channels_per_layer, -1.0);
                    std::vector<double> width_arr(n_channels_per_layer, -1.0);
                    std::vector<double> total_area_arr(n_channels_per_layer, -1.0);
                    std::vector<double> gaus_sigma_arr(n_channels_per_layer, -1.0);
                    std::vector<double> max_x_arr(n_channels_per_layer, -1.0);
                    std::vector<double> fwhm_arr(n_channels_per_layer, -1.0);
                    std::vector<int> entries_arr(n_channels_per_layer, 0);
                    std::vector<double> chi2_arr(n_channels_per_layer, -1.0);
                    std::vector<int> ndf_arr(n_channels_per_layer, 0);
                    std::vector<int> fit_status_arr(n_channels_per_layer, 999);
                    std::vector<int> fit_ok_arr(n_channels_per_layer, 0);
                    
                    // SPE-specific arrays
                    std::vector<double> gain_arr(n_channels_per_layer, -1.0);
                    std::vector<double> gainErr_arr(n_channels_per_layer, -1.0);
                    std::vector<double> snr_arr(n_channels_per_layer, -1.0);
                    std::vector<double> peakAmp_arr(n_channels_per_layer, -1.0);
                    std::vector<double> noiseMed_arr(n_channels_per_layer, -1.0);

                    long long total_entries = 0;
                    int fit_failures = 0;

                    for (const auto& [_, res] : fit_cache_) {
                        if (res.layer != layer) continue;
                        const int idx = res.chip * AHCALGeometry::channel_No + res.channel;
                        if (idx < 0 || idx >= n_channels_per_layer) continue;

                        mpv_arr[idx] = res.fit_out.mpv;
                        width_arr[idx] = res.fit_out.width;
                        total_area_arr[idx] = res.fit_out.total_area;
                        gaus_sigma_arr[idx] = res.fit_out.gaus_sigma;
                        max_x_arr[idx] = res.fit_out.max_x;
                        fwhm_arr[idx] = res.fit_out.FWHM;
                        entries_arr[idx] = res.entries;
                        chi2_arr[idx] = res.fit_out.chi2;
                        ndf_arr[idx] = res.fit_out.ndf;
                        fit_status_arr[idx] = res.fit_out.fit_status;
                        fit_ok_arr[idx] = res.fit_out.fit_ok ? 1 : 0;
                        
                        // SPE-specific
                        gain_arr[idx] = res.gain;
                        gainErr_arr[idx] = res.gainErr;
                        snr_arr[idx] = res.snr;
                        peakAmp_arr[idx] = res.peakAmp;
                        noiseMed_arr[idx] = res.noiseMed;
                        
                        total_entries += res.entries;
                        if (!res.fit_out.fit_ok) {
                            fit_failures++;
                        }
                    }

                    json j;
                    j["RunNumber"] = rc.config.runNumber;
                    j["TimeStamp"] = timestamp;
                    j["Layer"] = layer;
                    j["CalibrationType"] = "SPE";
                    j["Summary"]["FitOK"] = nFitOK_;
                    j["Summary"]["FitAll"] = nFitAll_;
                    j["Summary"]["Entries"] = total_entries;
                    j["Summary"]["NumUsedRun"] = static_cast<int>(run_contexts_.size());
                    j["Summary"]["SameDataRuns"] = same_data_runs;
                    if (!cfg_.fit) {
                        fit_failures = 0;
                    }
                    j["Status"] = (fit_failures == 0) ? 0 : 1;

                    j["PerChannel"]["MPV"] = mpv_arr;
                    j["PerChannel"]["Width"] = width_arr;
                    j["PerChannel"]["TotalArea"] = total_area_arr;
                    j["PerChannel"]["GausSigma"] = gaus_sigma_arr;
                    j["PerChannel"]["MaxX"] = max_x_arr;
                    j["PerChannel"]["FWHM"] = fwhm_arr;
                    j["PerChannel"]["Entries"] = entries_arr;
                    j["PerChannel"]["Chi2"] = chi2_arr;
                    j["PerChannel"]["NDF"] = ndf_arr;
                    j["PerChannel"]["FitStatus"] = fit_status_arr;
                    j["PerChannel"]["FitOK"] = fit_ok_arr;
                    
                    // SPE-specific per-channel data
                    j["PerChannel"]["Gain"] = gain_arr;
                    j["PerChannel"]["GainErr"] = gainErr_arr;
                    j["PerChannel"]["SNR"] = snr_arr;
                    j["PerChannel"]["PeakAmp"] = peakAmp_arr;
                    j["PerChannel"]["NoiseMed"] = noiseMed_arr;

                    std::ostringstream filename;
                    if (!out_dir.empty()) {
                        filename << out_dir.string() << "/";
                    }
                    filename << "run" << rc.config.runNumber << "_spe_Layer" << layer << ".json";

                    const std::string out_filename = filename.str();
                    std::ofstream jout(out_filename);
                    if (!jout.is_open()) {
                        LOG_ERROR("SPEAlg: cannot write JSON {}", out_filename);
                        continue;
                    }
                    jout << j.dump(2) << std::endl;
                    jout.close();
                    LOG_INFO("SPEAlg: wrote JSON {}", out_filename);
                }
            }
        }

        void write() {
            if (!cfg_.mip_to_file && !cfg_.mip_to_json) return;
            if (written_) return;
            written_ = true;

            buildFitCache();
            doFFTAnalysis();

            if (cfg_.mip_to_file) {
                writeRootOutput();
            }
            if (cfg_.mip_to_json) {
                writeJsonOutput();
            }
        }

        std::unique_ptr<TH1D> make_hist(int cellid) {
            int layer = cellid / 100000;
            int chip = cellid / 10000 % 10;
            const std::string name = "hMIP_" + std::to_string(cellid);
            const std::string title = "Layer " + std::to_string(layer) + " chip " + 
                                     std::to_string(chip) + " MIP;ADC (ped-sub);counts";
            auto h = std::make_unique<TH1D>(name.c_str(), title.c_str(), cfg_.nbin, cfg_.xmin, cfg_.xmax);
            h->SetDirectory(nullptr);
            return h;
        }

        SPEAlgCfg cfg_;
        RunContext ctx_;
        bool written_ = false;
        bool cache_built_ = false;
        bool fft_done_ = false;
        int nFitAll_ = 0;
        int nFitOK_ = 0;
        long long n_missing_ped_ = 0;
        std::vector<RunContext> run_contexts_;
        std::unordered_map<int, CalibDBIO::Pedestal> ped_map_;
        std::unordered_map<int, std::unique_ptr<TH1D>> hg_hist_;
        std::unordered_map<int, FitResult> fit_cache_;
    };

    void SPEAlg::ImplDeleter::operator()(Impl* p) const {
        delete p;
    }

    void SPEAlg::ensure_impl() {
        if (!impl_) {
            impl_.reset(new Impl(cfg_, ctx()));
        }
    }

    SPEAlg::~SPEAlg() {
        if (impl_) impl_->write();
    }

    void SPEAlg::init_by_run() {
        const bool was_uninitialized = !impl_;
        ensure_impl();
        if (!was_uninitialized) {
            impl_->change_run(ctx());
            impl_->loadPedestals();
            LOG_INFO("SPEAlg: re-loaded pedestals for new run {}", ctx().config.runNumber);
        }
    }

    void SPEAlg::execute(EventStore& evt) {
        ensure_impl();
        auto rawhits = evt.get<std::vector<AHCALRawHit>>(cfg_.in_rawhit_key);

        if (cfg_.string_track_struct == "SimpleFittedTrack") {
            auto track = evt.get<SimpleFittedTrack>(cfg_.in_track_key);
            if (!cfg_.track_selection_string.empty()) {
                SimpleFittedTrackCut cut(cfg_.track_selection_string);
                if (!cut.eval(track)) {
                    return;
                }
            }
            if (!track.inTrackHitsIndices.size()) {
                LOG_WARN("SPEAlg: track has no associated hits.");
                return;
            }
            for (const int index : track.inTrackHitsIndices) {
                if (index < 0 || index >= static_cast<int>(rawhits.size())) {
                    continue;
                }
                const auto& rh = rawhits.at(index);
                if (rh.index != index) {
                    continue;
                }
                impl_->fill(rh);
            }
        } else if (cfg_.string_track_struct == "Track") {
            auto track = evt.get<Track>(cfg_.in_track_key);
            if (!cfg_.track_selection_string.empty()) {
                TrackCut cut(cfg_.track_selection_string);
                if (!cut.eval(track)) {
                    return;
                }
            }
            if (!track.inTrackHitsIndices.size()) {
                LOG_WARN("SPEAlg: track has no associated hits.");
                return;
            }
            for (const int index : track.inTrackHitsIndices) {
                if (index < 0 || index >= static_cast<int>(rawhits.size())) {
                    continue;
                }
                const auto& rh = rawhits.at(index);
                if (rh.index != index) {
                    continue;
                }
                impl_->fill(rh);
            }
        } else {
            LOG_ERROR("SPEAlg: unknown track struct string '{}'", cfg_.string_track_struct);
            return;
        }
    }

    void SPEAlg::parse_cfg(const YAML::Node& cfg) {
        cfg_.in_rawhit_key = get_or<std::string>(cfg, "in_rawhit_key", cfg_.in_rawhit_key);
        cfg_.in_track_key = get_or<std::string>(cfg, "in_track_key", cfg_.in_track_key);
        cfg_.string_track_struct = get_or<std::string>(cfg, "string_track_struct", cfg_.string_track_struct);
        cfg_.track_selection_string = get_or<std::string>(cfg, "track_selection_string", cfg_.track_selection_string);
        cfg_.mip_to_file = get_or<bool>(cfg, "mip_to_file", cfg_.mip_to_file);
        cfg_.mip_to_json = get_or<bool>(cfg, "mip_to_json", cfg_.mip_to_json);
        cfg_.output_to_png = get_or<bool>(cfg, "output_to_png", cfg_.output_to_png);
        cfg_.out_mip_filename = get_or<std::string>(cfg, "out_mip_filename", cfg_.out_mip_filename);
        cfg_.out_png_dir = get_or<std::string>(cfg, "out_png_dir", cfg_.out_png_dir);
        cfg_.out_json_dirname = get_or<std::string>(cfg, "out_json_dirname", cfg_.out_json_dirname);
        cfg_.nbin = get_or<int>(cfg, "nbin", cfg_.nbin);
        cfg_.xmin = get_or<double>(cfg, "xmin", cfg_.xmin);
        cfg_.xmax = get_or<double>(cfg, "xmax", cfg_.xmax);
        cfg_.min_entries = get_or<int>(cfg, "min_entries", cfg_.min_entries);
        cfg_.npad_req = get_or<int>(cfg, "npad_req", cfg_.npad_req);
        cfg_.gain_min = get_or<double>(cfg, "gain_min", cfg_.gain_min);
        cfg_.gain_max = get_or<double>(cfg, "gain_max", cfg_.gain_max);
        cfg_.min_ent_fft = get_or<int>(cfg, "min_ent_fft", cfg_.min_ent_fft);
        cfg_.edge_gain_tol = get_or<double>(cfg, "edge_gain_tol", cfg_.edge_gain_tol);
        cfg_.kedge_skip = get_or<int>(cfg, "kedge_skip", cfg_.kedge_skip);
        cfg_.save_per_channel_hists = get_or<bool>(cfg, "save_per_channel_hists", cfg_.save_per_channel_hists);
        cfg_.save_max_channel_hists = get_or<int>(cfg, "save_max_channel_hists", cfg_.save_max_channel_hists);
        cfg_.do_event_scan = get_or<bool>(cfg, "do_event_scan", cfg_.do_event_scan);
        cfg_.fit = get_or<bool>(cfg, "fit", cfg_.fit);
        cfg_.save_examples = get_or<bool>(cfg, "save_examples", cfg_.save_examples);
        cfg_.substrate_pedestal = get_or<bool>(cfg, "substrate_pedestal", cfg_.substrate_pedestal);
        cfg_.read_pedestal_from_ROOT = get_or<bool>(cfg, "read_pedestal_from_ROOT", cfg_.read_pedestal_from_ROOT);
        cfg_.in_pedestal_file = get_or<std::string>(cfg, "in_pedestal_file", cfg_.in_pedestal_file);
        cfg_.read_pedestal_from_DB = get_or<bool>(cfg, "read_pedestal_from_DB", cfg_.read_pedestal_from_DB);
        cfg_.processing_mode = get_or<int>(cfg, "processing_mode", cfg_.processing_mode);
        cfg_.snr_threshold = get_or<double>(cfg, "snr_threshold", cfg_.snr_threshold);
        cfg_.calculate_fwhm = get_or<bool>(cfg, "calculate_fwhm", cfg_.calculate_fwhm);
    }

} // namespace AHCALRecoAlg
