#include "SPEAlg.hpp"
#include "fft_helper.h"
#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"
#include "CalibDBIO/PedestalReader/PedestalReader.hpp"
#include "common/edm/SimpleFittedTrack.hpp"
#include "common/edm/Track.hpp"
#include "calibration/RefValues.hpp"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TCanvas.h>
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

    // ============ Directory utilities ============
    static TDirectory* ensureDir(TDirectory* top, const char* name) {
        if (!top) return nullptr;
        auto* d = dynamic_cast<TDirectory*>(top->Get(name));
        if (!d) d = top->mkdir(name);
        return d;
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
        struct ChannelResult {
            int cellid = -1;
            int layer = -1;
            int chip = -1;
            int channel = -1;
            int channel_index = -1;
            double x_mm = -999.0;
            double y_mm = -999.0;
            int entries = 0;
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
            ped_map_ = std::make_shared<CalibDBIO::PedestalMap>();
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

            auto ped_map = std::make_shared<CalibDBIO::PedestalMap>();
            ped_map->reserve(tp->GetEntries());
            for (Long64_t i = 0; i < tp->GetEntries(); ++i) {
                tp->GetEntry(i);
                CalibDBIO::Pedestal ped{};
                ped.HighGainPeak = p_hgped;
                ped.LowGainPeak = has_lg_ped ? p_lgped : 0.0;
                (*ped_map)[p_cellid] = ped;
            }
            ped_map_ = std::move(ped_map);

            LOG_INFO("SPEAlg: loaded {} pedestal entries from ROOT", ped_map_->size());
        }

        void loadPedestals_fromDB() {
            CalibDBIO::PedestalReader reader(ctx_.config.runNumber);
            ped_map_ = reader.getPedestalMapPtr();
            LOG_INFO("SPEAlg: loaded {} pedestal entries from DB", ped_map_->size());
        }

        void change_run(const RunContext& new_ctx) {
            ctx_ = new_ctx;
            run_contexts_.push_back(new_ctx);
        }

        void fill(const AHCALRawHit& h) {
            const int cellid = h.cellID;

            double hg_value = static_cast<double>(h.hg_adc);
            if (cfg_.substrate_pedestal) {
                auto itp = ped_map_->find(cellid);
                if (itp == ped_map_->end()) {
                    n_missing_ped_++;
                    return;
                }
                if (!AHCALRefValues::HGPedestalStatus_is_ok(itp->second.HighGainStatus)) {
                    n_missing_ped_++; // Treat non-ok pedestal as missing
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

        void buildResultCache() {
            if (cache_built_) return;
            cache_built_ = true;

            result_cache_.clear();
            result_cache_.reserve(hg_hist_.size());

            for (const auto& [cid, hist_ptr] : hg_hist_) {
                ChannelResult res;
                res.cellid = cid;
                res.layer = cid / 100000;
                res.chip = (cid / 10000) % 10;
                res.channel = cid % 10000;
                res.channel_index = res.layer * (AHCALGeometry::chip_No * AHCALGeometry::channel_No) +
                                    res.chip * AHCALGeometry::channel_No + res.channel;
                cellid_to_xy(res.chip, res.channel, res.x_mm, res.y_mm);

                res.entries = hist_ptr ? static_cast<int>(hist_ptr->GetEntries()) : 0;
                result_cache_.emplace(cid, std::move(res));
            }

            LOG_INFO("SPEAlg: built result cache for {} channels", result_cache_.size());
        }

        void doFFTAnalysis() {
            if (fft_done_) return;
            fft_done_ = true;

            int count = 0;
            for (auto& [cid, res] : result_cache_) {
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
                KPick kp = pickK0(hFFT.get(), Npad, h->GetBinWidth(1),
                                  cfg_.gain_min, cfg_.gain_max,
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

        void fillMapsFromCache(std::vector<std::unique_ptr<TH2D>>& hEntries,
                                std::vector<std::unique_ptr<TH2D>>& hGain,
                                std::vector<std::unique_ptr<TH2D>>& hSNR) {
            for (const auto& [_, res] : result_cache_) {
                if (res.layer >= 0 && res.layer < AHCALGeometry::Layer_No) {
                    const int binx = hGain[res.layer]->GetXaxis()->FindBin(res.x_mm);
                    const int biny = hGain[res.layer]->GetYaxis()->FindBin(res.y_mm);
                    if (binx > 0 && binx <= 18 && biny > 0 && biny <= 18) {
                        hEntries[res.layer]->SetBinContent(binx, biny, res.entries);
                        if (res.snr >= cfg_.snr_threshold) {
                            hGain[res.layer]->SetBinContent(binx, biny, res.gain);
                            hSNR[res.layer]->SetBinContent(binx, biny, res.snr);
                        }
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

            // Build channel result cache if needed
            if (!cache_built_) buildResultCache();
            if (!fft_done_) doFFTAnalysis();

            // Create output maps
            const double XYMIN = -AHCALGeometry::x_max;
            const double XYMAX = +AHCALGeometry::x_max;
            const int NBIN_XY = 18;

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
                hEntries[L] = makeMap("Entries", L, "SPE spectrum entries");
                hGain[L] = makeMap("Gain", L, "FFT Gain");
                hSNR[L] = makeMap("SNR", L, "FFT SNR");
            }

            // Fill maps from cache
            fillMapsFromCache(hEntries, hGain, hSNR);

            // Build output tree
            TTree tp("spe", "SPE (Signal-to-noise ratio, Pedestal subtraction, Efficiency) analysis");
            int cellid = -1;
            int layer = -1, chip = -1, channel = -1;
            int entries = 0;
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
            tp.Branch("entries", &entries);
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

            // Layer efficiency stats
            int denLayer[AHCALGeometry::Layer_No] = {0};
            int numLayer[AHCALGeometry::Layer_No] = {0};
            double bestSNR = -1;
            int bestCellID = -1;

            int saved_hists = 0;

            for (const auto& [cid, res] : result_cache_) {
                cellid = res.cellid;
                layer = res.layer;
                chip = res.chip;
                channel = res.channel;
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
                const bool save_all_channel_hists = (cfg_.save_max_channel_hists < 0);
                if (cfg_.save_per_channel_hists && dHist &&
                    (save_all_channel_hists || saved_hists < cfg_.save_max_channel_hists)) {
                    auto it = hg_hist_.find(cid);
                    if (it != hg_hist_.end() && it->second) {
                        // Create nested directory structure
                        TDirectory* dNested = dynamic_cast<TDirectory*>(dHist->Get(Form("Layer%02d/Chip%d", layer, chip)));
                        if (!dNested) dNested = dHist;
                        
                        TH1D* hSpectrum = it->second.get();
                        if (dNested) dNested->cd();
                        hSpectrum->Write(Form("hSPE_L%d_C%d_ch%d", layer, chip, channel));
                        
                        // Compute and save FFT histograms if SNR is good
                        if (res.snr >= cfg_.snr_threshold && res.snr > 0) {
                            int bmax = hSpectrum->GetNbinsX();
                            
                            std::vector<double> v(bmax, 0.0);
                            for (int b = 1; b <= bmax; ++b) {
                                v[b - 1] = hSpectrum->GetBinContent(b);
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
                                
                                std::unique_ptr<TH1D> hFFT_gain(makeFFTgainHist(
                                    hFFT_k.get(), res.npad, hSpectrum->GetBinWidth(1),
                                    cfg_.gain_min, cfg_.gain_max,
                                    Form("hFFT_gain_L%d_C%d_ch%d", layer, chip, channel)));
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
                auto it_spectrum = hg_hist_.find(bestCellID);
                auto it_result = result_cache_.find(bestCellID);
                if (it_spectrum != hg_hist_.end() && it_spectrum->second &&
                    it_result != result_cache_.end()) {
                    TH1D* hSpectrum = it_spectrum->second.get();
                    const ChannelResult& bestRes = it_result->second;
                    
                    // Create examples subdirectory
                    auto dExamples = ensureDir(fout.get(), "examples/best_channel");
                    
                    dExamples->cd();
                    hSpectrum->Write("hSPE_best");
                    
                    // Compute and save FFT histograms for best channel
                    int bmax = hSpectrum->GetNbinsX();
                    
                    std::vector<double> v(bmax, 0.0);
                    for (int b = 1; b <= bmax; ++b) {
                        v[b - 1] = hSpectrum->GetBinContent(b);
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
                        
                        std::unique_ptr<TH1D> hFFT_gain(makeFFTgainHist(
                            hFFT_k.get(), bestRes.npad, hSpectrum->GetBinWidth(1),
                            cfg_.gain_min, cfg_.gain_max, "hFFT_best_gain"));
                        hFFT_gain->Write();
                    }
                    
                    // Save PNG if enabled
                    if (cfg_.output_to_png) {
                        std::filesystem::path png_dir(cfg_.out_png_dir);
                        std::error_code ec;
                        std::filesystem::create_directories(png_dir, ec);
                        
                        auto c_spe = std::make_unique<TCanvas>("c_spe", "c_spe", 900, 650);
                        hSpectrum->Draw("E");
                        c_spe->SaveAs(Form("%s/h_SPE_best.png", cfg_.out_png_dir.c_str()));
                        
                        if (hFFT_k) {
                            auto c_fft = std::make_unique<TCanvas>("c_fft", "c_fft", 900, 650);
                            hFFT_k->Draw();
                            c_fft->SaveAs(Form("%s/h_FFT_best_k.png", cfg_.out_png_dir.c_str()));
                        }
                    }
                    
                    // Write canvas to ROOT file
                    auto can = std::make_unique<TCanvas>("cBestExample", "Best SNR channel", 800, 600);
                    hSpectrum->Draw("E");
                    if (dCan) dCan->cd();
                    can->Write();
                }
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
                hEntries[L]->Write();
                hGain[L]->Write();
                hSNR[L]->Write();
            }

            if (dCan) dCan->cd();
            gLayerEff->Write();

            fout->cd();
            tp.Write();
            fout->Close();

            LOG_INFO("SPEAlg: wrote {} ({} per-channel hists)", cfg_.out_mip_filename, saved_hists);
        }

        void writeJsonOutput() {
            if (!cfg_.mip_to_json) return;
            if (!cache_built_) buildResultCache();
            if (!fft_done_) doFFTAnalysis();

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

            // Per run, per layer SPE JSON
            for (const auto& rc : run_contexts_) {
                for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer) {
                    // Initialize per-channel arrays
                    std::vector<int> entries_arr(n_channels_per_layer, 0);
                    std::vector<double> gain_arr(n_channels_per_layer, -1.0);
                    std::vector<double> gainErr_arr(n_channels_per_layer, -1.0);
                    std::vector<double> snr_arr(n_channels_per_layer, -1.0);
                    std::vector<double> peakAmp_arr(n_channels_per_layer, -1.0);
                    std::vector<double> noiseMed_arr(n_channels_per_layer, -1.0);

                    long long total_entries = 0;
                    int successful_channels = 0;

                    for (const auto& [_, res] : result_cache_) {
                        if (res.layer != layer) continue;
                        const int idx = res.chip * AHCALGeometry::channel_No + res.channel;
                        if (idx < 0 || idx >= n_channels_per_layer) continue;

                        entries_arr[idx] = res.entries;
                        gain_arr[idx] = res.gain;
                        gainErr_arr[idx] = res.gainErr;
                        snr_arr[idx] = res.snr;
                        peakAmp_arr[idx] = res.peakAmp;
                        noiseMed_arr[idx] = res.noiseMed;
                        
                        total_entries += res.entries;
                        if (res.snr >= cfg_.snr_threshold) successful_channels++;
                    }

                    json j;
                    j["RunNumber"] = rc.config.runNumber;
                    j["TimeStamp"] = timestamp;
                    j["Layer"] = layer;
                    j["CalibrationType"] = "SPE";
                    j["Summary"]["SuccessfulChannels"] = successful_channels;
                    j["Summary"]["Entries"] = total_entries;
                    j["Summary"]["NumUsedRun"] = static_cast<int>(run_contexts_.size());
                    j["Summary"]["SameDataRuns"] = same_data_runs;
                    j["Status"] = 0;
                    j["PerChannel"]["Entries"] = entries_arr;
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

            buildResultCache();
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
            const std::string name = "hSPE_" + std::to_string(cellid);
            const std::string title = "Layer " + std::to_string(layer) + " chip " + 
                                     std::to_string(chip) + " SPE spectrum;ADC (ped-sub);counts";
            auto h = std::make_unique<TH1D>(name.c_str(), title.c_str(), cfg_.nbin, cfg_.xmin, cfg_.xmax);
            h->SetDirectory(nullptr);
            return h;
        }

        SPEAlgCfg cfg_;
        RunContext ctx_;
        bool written_ = false;
        bool cache_built_ = false;
        bool fft_done_ = false;
        long long n_missing_ped_ = 0;
        std::vector<RunContext> run_contexts_;
        std::shared_ptr<const CalibDBIO::PedestalMap> ped_map_ = std::make_shared<CalibDBIO::PedestalMap>();
        std::unordered_map<int, std::unique_ptr<TH1D>> hg_hist_;
        std::unordered_map<int, ChannelResult> result_cache_;
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
                // if (rh.index != index) {
                //     continue;
                // }
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
        cfg_.npad_req = get_or<int>(cfg, "npad_req", cfg_.npad_req);
        cfg_.gain_min = get_or<double>(cfg, "gain_min", cfg_.gain_min);
        cfg_.gain_max = get_or<double>(cfg, "gain_max", cfg_.gain_max);
        cfg_.min_ent_fft = get_or<int>(cfg, "min_ent_fft", cfg_.min_ent_fft);
        cfg_.edge_gain_tol = get_or<double>(cfg, "edge_gain_tol", cfg_.edge_gain_tol);
        cfg_.kedge_skip = get_or<int>(cfg, "kedge_skip", cfg_.kedge_skip);
        cfg_.save_per_channel_hists = get_or<bool>(cfg, "save_per_channel_hists", cfg_.save_per_channel_hists);
        cfg_.save_max_channel_hists = get_or<int>(cfg, "save_max_channel_hists", cfg_.save_max_channel_hists);
        cfg_.do_event_scan = get_or<bool>(cfg, "do_event_scan", cfg_.do_event_scan);
        cfg_.save_examples = get_or<bool>(cfg, "save_examples", cfg_.save_examples);
        cfg_.substrate_pedestal = get_or<bool>(cfg, "substrate_pedestal", cfg_.substrate_pedestal);
        cfg_.read_pedestal_from_ROOT = get_or<bool>(cfg, "read_pedestal_from_ROOT", cfg_.read_pedestal_from_ROOT);
        cfg_.in_pedestal_file = get_or<std::string>(cfg, "in_pedestal_file", cfg_.in_pedestal_file);
        cfg_.read_pedestal_from_DB = get_or<bool>(cfg, "read_pedestal_from_DB", cfg_.read_pedestal_from_DB);
        cfg_.processing_mode = get_or<int>(cfg, "processing_mode", cfg_.processing_mode);
        cfg_.snr_threshold = get_or<double>(cfg, "snr_threshold", cfg_.snr_threshold);
    }

} // namespace AHCALRecoAlg
