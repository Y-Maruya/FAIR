#include "InterCalibAlg.hpp"

#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"
#include "CalibDBIO/PedestalReader/PedestalReader.hpp"
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TString.h>
#include <TLine.h>

#include <algorithm>
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
AHCAL_REGISTER_ALG(AHCALRecoAlg::InterCalibAlg, "InterCalibAlg")

namespace AHCALRecoAlg {

namespace {


// Accumulation data for linear regression
struct Acc {
    // Store points for HG-binned profile fit:
    // fit LG = c*HG + d with sigma from RMS(LG) in each HG bin.
    std::vector<std::pair<double, double>> points_hg_lg; // (HG, LG)
    long long n  = 0;
};

// Linear fit result (HG = p0 + p1*LG)
struct FitResult {
    double p0 = -999.0;
    double p1 = -999.0;
    double chi2 = -1.0;
    int ndf = -1;
    double chi2_ndf = -1.0;
    int status = 999;
    bool ok = false;
};

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

static FitResult fitLinearFromHGBinnedProfile(const Acc& a, int minPoints, double hgBinWidth, int minBins) {
    FitResult r;
    if (a.n < minPoints || a.points_hg_lg.empty() || hgBinWidth <= 0.0) return r;

    struct BinStats { int n = 0; double sum_lg = 0.0; double sum_lg2 = 0.0; double sum_hg = 0.0; };
    std::unordered_map<int, BinStats> bins;
    bins.reserve(a.points_hg_lg.size() / 8 + 1);
    for (const auto& p : a.points_hg_lg) {
        const double hg = p.first;
        const double lg = p.second;
        const int ibin = static_cast<int>(std::floor(hg / hgBinWidth));
        auto& b = bins[ibin];
        b.n++;
        b.sum_hg += hg;
        b.sum_lg += lg;
        b.sum_lg2 += lg * lg;
    }
    if (static_cast<int>(bins.size()) < minBins) return r;

    // Weighted fit of LG = c*HG + d
    double S = 0.0, SX = 0.0, SY = 0.0, SXX = 0.0, SXY = 0.0;
    int used_bins = 0;
    for (const auto& [_, b] : bins) {
        if (b.n < 3) continue;
        const double mean_hg = b.sum_hg / b.n;
        const double mean_lg = b.sum_lg / b.n;
        const double var_lg = std::max(0.0, b.sum_lg2 / b.n - mean_lg * mean_lg);
        const double sigma_lg = std::sqrt(var_lg);
        const double sigma_eff = std::max(sigma_lg, 1.0); // avoid overweighting tiny-RMS bins
        const double w = 1.0 / (sigma_eff * sigma_eff);
        S += w;
        SX += w * mean_hg;
        SY += w * mean_lg;
        SXX += w * mean_hg * mean_hg;
        SXY += w * mean_hg * mean_lg;
        used_bins++;
    }
    if (used_bins < minBins) return r;
    const double denom = S * SXX - SX * SX;
    if (denom <= 0.0) {
        r.status = 1;
        return r;
    }
    const double c = (S * SXY - SX * SY) / denom;
    const double d = (SY * SXX - SX * SXY) / denom;
    if (std::abs(c) < 1e-9) {
        r.status = 2;
        return r;
    }
    // Convert LG=c*HG+d -> HG=(1/c)*LG + (-d/c)
    r.p1 = 1.0 / c;
    r.p0 = -d / c;
    r.ndf = used_bins - 2;
    if (r.ndf <= 0) return r;
    r.chi2 = 0.0;
    for (const auto& [_, b] : bins) {
        if (b.n < 3) continue;
        const double mean_hg = b.sum_hg / b.n;
        const double mean_lg = b.sum_lg / b.n;
        const double var_lg = std::max(0.0, b.sum_lg2 / b.n - mean_lg * mean_lg);
        const double sigma_lg = std::max(std::sqrt(var_lg), 1.0);
        const double resid = mean_lg - (c * mean_hg + d);
        r.chi2 += (resid * resid) / (sigma_lg * sigma_lg);
    }
    r.chi2_ndf = r.chi2 / r.ndf;
    r.status = 0;
    r.ok = true;
    return r;
}

} // namespace

struct InterCalibAlg::Impl {
    explicit Impl(InterCalibAlgCfg cfg, const RunContext& ctx)
        : cfg_(std::move(cfg)), ctx_(ctx) {
        loadPedestals();
        run_contexts_.push_back(ctx);
    }

    struct CellInterCalibResult {
        int cellid = -1;
        int layer = -1;
        int chip = -1;
        int channel = -1;
        int channel_index = -1;
        double x_mm = -999.0;
        double y_mm = -999.0;

        long long n_points = 0;
        double p0 = -999.0;
        double p1 = -999.0;
        double chi2 = -1.0;
        int ndf = -1;
        double chi2_ndf = -1.0;
        int fit_status = 999;
        bool fit_ok = false;
    };
    void loadPedestals() {
        if (cfg_.read_pedestal_from_ROOT) {
            loadPedestals_fromROOT();
        } else if (cfg_.read_pedestal_from_DB) {
            loadPedestals_fromDB();
        } else {
            LOG_ERROR("InterCalibAlg: no pedestal source specified");
        }
    }
    void loadPedestals_fromROOT() {
        std::unique_ptr<TFile> fped(TFile::Open(cfg_.in_pedestal_file.c_str(), "READ"));
        if (!fped || fped->IsZombie()) {
            LOG_ERROR("InterCalibAlg: cannot open pedestal file: {}", cfg_.in_pedestal_file);
            return;
        }

        TTree* tp = static_cast<TTree*>(fped->Get("pedestal"));
        if (!tp) {
            LOG_ERROR("InterCalibAlg: TTree 'pedestal' not found in {}", cfg_.in_pedestal_file);
            return;
        }

        int p_cellid = 0;
        double p_hgped = 0.0, p_lgped = 0.0;

        tp->SetBranchAddress("cellid", &p_cellid);
        tp->SetBranchAddress("highgain_peak", &p_hgped);

        bool has_lg_ped = false;
        if (tp->GetBranch("lowgain_peak")) {
            tp->SetBranchAddress("lowgain_peak", &p_lgped);
            has_lg_ped = true;
        }

        ped_map_.clear();
        ped_map_.reserve(tp->GetEntries());

        for (Long64_t i = 0; i < tp->GetEntries(); ++i) {
            tp->GetEntry(i);
            ped_map_[p_cellid].HighGainPeak = p_hgped;
            if (has_lg_ped) {
                ped_map_[p_cellid].LowGainPeak = p_lgped;
            }
        }

        LOG_INFO("InterCalibAlg: loaded {} pedestal entries", ped_map_.size());
    }
    void loadPedestals_fromDB() {
        CalibDBIO::PedestalReader reader(ctx_.config.runNumber);
        ped_map_ = reader.getPedestalMap();
        LOG_INFO("InterCalibAlg: loaded {} pedestal entries from DB", ped_map_.size());
    }
    void fill(const AHCALRawHit& h) {
        int cellid = h.cellID;

        // Find pedestal for this cell
        auto itp = ped_map_.find(cellid);
        if (itp == ped_map_.end()) {
            n_missing_ped_++;
            return;
        }

        int hg = h.hg_adc;
        int lg = h.lg_adc;

        // Pedestal subtraction
        double hg_sub = hg - itp->second.HighGainPeak;
        double lg_sub = lg - itp->second.LowGainPeak;

        // Keep macro-compatible fit domain.
        // This avoids bias from negative/zero pedestal-subtracted points.
        if (hg_sub <= cfg_.hg_fit_min) return;
        if (lg_sub <= cfg_.lg_fit_min) return;
        if (hg_sub > cfg_.hg_fit_max) return;
        if (cfg_.require_hittag && h.hittag != 1) return;
        if (cfg_.require_yperx_over1 && hg_sub <= lg_sub) return;
        if (hg <= 0 || lg <= 0 || hg >= 4095 || lg >= 4095) return;
        if (cfg_.use_specific_fit_range_toL39C8) {
            int tmp_layer, tmp_chip, tmp_channel;
            AHCALGeometry::CellIDToLC(cellid, tmp_layer, tmp_chip, tmp_channel);
            if (tmp_layer == 39 && tmp_chip == 8) {
                if (hg_sub > cfg_.hg_fit_max_L39C8) return;
            }
        }
        if (lg_sub >200 && hg_sub < cfg_.hg_fit_max){
            // LOG_WARN("InterCalibAlg: detected hit with large LG and small HG, likely saturation: cellid={}, hg={}, lg={}, hg_sub={}, lg_sub={}", cellid, hg, lg, hg_sub, lg_sub);
            int tmp_layer, tmp_chip, tmp_channel;
            AHCALGeometry::CellIDToLC(cellid, tmp_layer, tmp_chip, tmp_channel);
            if (cfg_.output_bad_cells_json) {
                if (std::find(bad_fit_cells_.begin(), bad_fit_cells_.end(), std::make_tuple(tmp_layer, tmp_chip, tmp_channel)) == bad_fit_cells_.end()) {
                    bad_fit_cells_.emplace_back(std::make_tuple(tmp_layer, tmp_chip, tmp_channel));
                }
            }
            // return;
        }
        // Fill example 2D histogram
        int tmp_layer, tmp_chip, tmp_channel;
        AHCALGeometry::CellIDToLC(cellid, tmp_layer, tmp_chip, tmp_channel);
        for (size_t i = 0; i < h_example_.size(); ++i) {
            if (tmp_layer == cfg_.example_layers[i] &&
                tmp_chip == cfg_.example_chips[i] &&
                tmp_channel == cfg_.example_chs[i]) {
                h_example_[i]->Fill(lg_sub, hg_sub);
                break;
            }
        }

        // Accumulate for linear regression
        Acc& a = acc_[cellid];
        a.points_hg_lg.emplace_back(hg_sub, lg_sub);
        a.n++;

        n_used_++;
    }

    void buildInterCalibCache() {
        if (cache_built_) return;
        cache_built_ = true;

        ic_cache_.clear();
        ic_cache_.reserve(acc_.size());

        n_fit_ok_ = 0;
        n_fit_all_ = 0;

        for (const auto& [cid, a] : acc_) {
            CellInterCalibResult res;
            res.cellid = cid;

            int tmp_layer, tmp_chip, tmp_channel;
            AHCALGeometry::CellIDToLC(cid, tmp_layer, tmp_chip, tmp_channel);
            res.layer = tmp_layer;
            res.chip = tmp_chip;
            res.channel = tmp_channel;
            res.channel_index = AHCALGeometry::channel_index_from_lcc(res.layer, res.chip, res.channel);
            res.x_mm = AHCALGeometry::Pos_X(tmp_channel, tmp_chip);
            res.y_mm = AHCALGeometry::Pos_Y(tmp_channel, tmp_chip);

            res.n_points = a.n;
            if (a.n >= cfg_.min_points) {
                n_fit_all_++;
                FitResult fr = fitLinearFromHGBinnedProfile(a, cfg_.min_points, cfg_.hg_bin_width, cfg_.min_hg_bins_for_fit);
                res.p0 = fr.p0;
                res.p1 = fr.p1;
                res.chi2 = fr.chi2;
                res.ndf = fr.ndf;
                res.chi2_ndf = fr.chi2_ndf;
                res.fit_status = fr.status;
                res.fit_ok = fr.ok;
                if (fr.ok) n_fit_ok_++;
            }

            ic_cache_.emplace(cid, std::move(res));
        }
    }

    void write() {
        if (!cfg_.intercalib_to_file && !cfg_.intercalib_to_json) return;
        if (written_) return;
        written_ = true;

        buildInterCalibCache();

        if (cfg_.intercalib_to_file) {
            writeRoot();
        }

        if (cfg_.intercalib_to_json) {
            writeJson();
        }
    }

    void writeRoot() {
        auto fout = std::unique_ptr<TFile>(TFile::Open(cfg_.out_intercalib_filename.c_str(), "RECREATE"));
        if (!fout || fout->IsZombie()) {
            LOG_ERROR("InterCalibAlg: cannot create output file: {}", cfg_.out_intercalib_filename);
            return;
        }

        // Directories
        // TDirectory* dIC = ensureDir(fout.get(), "InterCalib");
        TDirectory* dCan = ensureDir(fout.get(), "Canvases");

        // Per-layer 2D maps
        const int NX = AHCALGeometry::Layer_No * AHCALGeometry::chip_No;
        const int NCH = AHCALGeometry::channel_No;

        auto hMapP1 = std::make_unique<TH2D>("hIntercalibSlope_Map",
            "HG/LG slope p1; [layer*9+chip]; Channel; Slope p1",
            NX, 0, NX, NCH, 0, NCH);
        hMapP1->SetMinimum(0.0);
        hMapP1->SetMaximum(40.0);
        hMapP1->SetDirectory(nullptr);

        auto hMapP0 = std::make_unique<TH2D>("hIntercalibIntercept_Map",
            "HG/LG intercept p0; [layer*9+chip]; Channel; Intercept p0",
            NX, 0, NX, NCH, 0, NCH);
        hMapP0->SetDirectory(nullptr);

        // 1D distributions
        auto hSlope = std::make_unique<TH1D>("hIntercalibSlope_1D",
            "HG/LG slope p1; p1; Entries", 400, 0, 40);
        hSlope->SetDirectory(nullptr);

        auto hIntercept = std::make_unique<TH1D>("hIntercalibIntercept_1D",
            "HG/LG intercept p0; p0; Entries", 400, -500, 500);
        hIntercept->SetDirectory(nullptr);

        // Per-layer maps
        std::vector<std::unique_ptr<TH2D>> hLayerP1(AHCALGeometry::Layer_No);
        for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
            hLayerP1[L] = std::make_unique<TH2D>(
                Form("hIntercalibSlope_L%02d", L),
                Form("Slope p1 Layer %02d; Chip; Channel; p1", L),
                AHCALGeometry::chip_No, 0, AHCALGeometry::chip_No,
                NCH, 0, NCH);
            hLayerP1[L]->SetMinimum(0.0);
            hLayerP1[L]->SetMaximum(40.0);
            hLayerP1[L]->SetDirectory(nullptr);
        }
        std::vector<FitResult> fit_results_example(h_example_.size());
        // Output tree
        TTree tIC("intercalib", "HG/LG intercalibration (HG = p0 + p1*LG, pedestal-subtracted)");

        int o_cellid = -1, o_channel_index = -1;
        int o_layer = -1, o_chip = -1, o_ch = -1;
        long long o_n = 0;
        double o_p0 = -999.0, o_p1 = -999.0;
        double o_chi2 = -1.0;
        int o_ndf = -1;
        double o_chi2_ndf = -1.0;
        int o_fit_status = 999;
        int o_fit_ok = 0;
        double o_x_mm = -999.0, o_y_mm = -999.0;
        std::vector<int> o_same_data_runs;
        tIC.Branch("same_data_runs", &o_same_data_runs);
        tIC.Branch("cellid", &o_cellid, "cellid/I");
        tIC.Branch("channel_index", &o_channel_index, "channel_index/I");
        tIC.Branch("layer", &o_layer, "layer/I");
        tIC.Branch("chip", &o_chip, "chip/I");
        tIC.Branch("ch", &o_ch, "ch/I");
        tIC.Branch("n_points", &o_n, "n_points/L");
        tIC.Branch("intercept", &o_p0, "intercept/D");
        tIC.Branch("slope", &o_p1, "slope/D");
        tIC.Branch("chi2", &o_chi2, "chi2/D");
        tIC.Branch("ndf", &o_ndf, "ndf/I");
        tIC.Branch("chi2_ndf", &o_chi2_ndf, "chi2_ndf/D");
        tIC.Branch("fit_status", &o_fit_status, "fit_status/I");
        tIC.Branch("fit_ok", &o_fit_ok, "fit_ok/I");
        tIC.Branch("x_mm", &o_x_mm, "x_mm/D");
        tIC.Branch("y_mm", &o_y_mm, "y_mm/D");

        o_same_data_runs.clear();
        for (const auto& rc : run_contexts_) {
            o_same_data_runs.push_back(rc.config.runNumber);
        }
        for (const auto& [_, res] : ic_cache_) {
            o_cellid = res.cellid;
            o_channel_index = res.channel_index;
            o_layer = res.layer;
            o_chip = res.chip;
            o_ch = res.channel;
            o_n = res.n_points;
            o_p0 = res.p0;
            o_p1 = res.p1;
            o_chi2 = res.chi2;
            o_ndf = res.ndf;
            o_chi2_ndf = res.chi2_ndf;
            o_fit_status = res.fit_status;
            o_fit_ok = res.fit_ok ? 1 : 0;
            o_x_mm = res.x_mm;
            o_y_mm = res.y_mm;

            hMapP1->SetBinContent(res.layer * AHCALGeometry::chip_No + res.chip + 1, res.channel + 1, res.p1);
            hMapP0->SetBinContent(res.layer * AHCALGeometry::chip_No + res.chip + 1, res.channel + 1, res.p0);

            if (res.fit_ok) {
                hLayerP1[res.layer]->SetBinContent(res.chip + 1, res.channel + 1, res.p1);
                hSlope->Fill(res.p1);
                hIntercept->Fill(res.p0);
            }
            for (size_t i = 0; i < cfg_.example_layers.size(); ++i) {
                if (res.layer == cfg_.example_layers[i] &&
                    res.chip == cfg_.example_chips[i] &&
                    res.channel == cfg_.example_chs[i]) {
                    fit_results_example[i].p0 = res.p0;
                    fit_results_example[i].p1 = res.p1;
                    fit_results_example[i].status = res.fit_status;
                    fit_results_example[i].ok = res.fit_ok;
                    fit_results_example[i].chi2 = res.chi2;
                    fit_results_example[i].ndf = res.ndf;
                    fit_results_example[i].chi2_ndf = res.chi2_ndf;
                    break;
                }
            }
            if (cfg_.output_bad_cells_json && o_p1 < 10.0 && o_n > 0) {
                bad_fit_cells_.emplace_back(std::make_tuple(res.layer, res.chip, res.channel));
            }
            tIC.Fill();
        }

        // Canvas (7x6) for 40 layers
        auto cAllP1 = std::make_unique<TCanvas>("cIntercalibSlope_all_7x6",
            "HG/LG slope maps (all layers)", 5600, 4200);
        cAllP1->Divide(7, 6, 0.001, 0.001);

        for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
            cAllP1->cd(L + 1);
            gPad->SetRightMargin(0.12);
            gPad->SetLeftMargin(0.10);
            gPad->SetBottomMargin(0.10);
            gPad->SetTopMargin(0.10);

            hLayerP1[L]->Draw("COLZ");
            drawLayerLabel(L);
        }

        // Write to ROOT file
        fout->cd();
        hMapP1->Write();
        hMapP0->Write();
        hSlope->Write();
        hIntercept->Write();
        for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
            hLayerP1[L]->Write();
        }
        for (size_t i = 0; i < h_example_.size(); ++i) {
            if (h_example_[i]) {
                TCanvas c(Form("cExample_L%02d_C%02d_ch%02d", cfg_.example_layers[i], cfg_.example_chips[i], cfg_.example_chs[i]),
                    Form("Example HG vs LG L%d C%d ch%d", cfg_.example_layers[i], cfg_.example_chips[i], cfg_.example_chs[i]),
                    800, 600);
                c.SetRightMargin(0.12);
                c.SetLeftMargin(0.10);
                c.SetBottomMargin(0.10);
                c.SetTopMargin(0.10);
                h_example_[i]->Draw("COLZ");
                example_fit_lines_[i]->SetLineColor(kRed);
                example_fit_lines_[i]->SetLineWidth(2);
                if (h_example_[i]->GetEntries() > 0 && i < fit_results_example.size() && fit_results_example[i].ok) {
                    // Fit a line to the 2D histogram
                    example_fit_lines_[i]->SetX1(0);
                    example_fit_lines_[i]->SetY1(fit_results_example[i].p0);
                    example_fit_lines_[i]->SetX2(200); // extend to large LG
                    example_fit_lines_[i]->SetY2(fit_results_example[i].p0 + fit_results_example[i].p1 * 200);
                    example_fit_lines_[i]->Draw("SAME");
                }
                TLatex t;
                t.SetNDC(true);
                t.SetTextSize(0.04);
                // t.DrawLatex(0.15, 0.85, Form("Entries: %lld", h_example_[i]->GetEntries()));
                t.DrawLatex(0.15, 0.80, Form("Fit: p0=%.2f, p1=%.2f, chi2/ndf=%.2f", fit_results_example[i].p0, fit_results_example[i].p1, fit_results_example[i].chi2_ndf));
                c.Write();
                h_example_[i]->Write();
            }
        }
        tIC.Write();

        if (dCan) {
            dCan->cd();
            cAllP1->Write();
        }

        fout->Close();

        LOG_INFO("InterCalibAlg: wrote {}", cfg_.out_intercalib_filename);
    }

    void writeJson() {
        // Create output directory if it doesn't exist
        std::filesystem::path out_dir(cfg_.out_json_dirname);
        try {
            std::filesystem::create_directories(out_dir);
        } catch (const std::exception& e) {
            LOG_ERROR("InterCalibAlg: cannot create output directory: {} - {}",
                      cfg_.out_json_dirname, e.what());
            return;
        }

        // Generate ISO8601 timestamp
        
        double start_time = run_contexts_.front().conditions.starttime;
        double end_time = run_contexts_.back().conditions.endtime;
        double ref_time = (start_time > 0.0 && end_time > 0.0) ?
                          (start_time + end_time) / 2.0 : std::time(nullptr);
        time_t utc_time = static_cast<time_t>(ref_time);
        std::tm tm = {};
        gmtime_r(&utc_time, &tm);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
        std::string timestamp = oss.str();

        const int n_channels_per_layer = AHCALGeometry::chip_No * AHCALGeometry::channel_No;

        std::vector<int> same_data_runs;
        for (const auto& rc : run_contexts_) {
            same_data_runs.push_back(rc.config.runNumber);
        }
        
        // Loop over each layer and write separate JSON file
        for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer) {
            // Prepare per-layer arrays (indexed by channel_index = chip*36 + channel)
            std::vector<double> p0_arr(n_channels_per_layer, -999.0);
            std::vector<double> p1_arr(n_channels_per_layer, -999.0);
            std::vector<int> fit_status_arr(n_channels_per_layer, 999);
            std::vector<long long> n_points_arr(n_channels_per_layer, 0);

            // Fill arrays for this layer
            int fit_failures = 0;
            long long total_entries = 0;
            LOG_INFO("InterCalibAlg: preparing JSON with reading cache, cache size = {}", ic_cache_.size());
            for (const auto& [_, res] : ic_cache_) {
                if (res.layer != layer) continue;

                const int idx = res.chip * AHCALGeometry::channel_No + res.channel;
                if (idx < 0 || idx >= n_channels_per_layer) continue;

                p0_arr[idx] = res.p0;
                p1_arr[idx] = res.p1;
                fit_status_arr[idx] = res.fit_status;
                n_points_arr[idx] = res.n_points;
                total_entries += res.n_points;

                if (!res.fit_ok) fit_failures++;
            }
            // Build JSON object for this layer
            for (const auto& rc : run_contexts_) {
                json j;
                j["RunNumber"] = rc.config.runNumber;
                j["TimeStamp"] = timestamp;
                j["Layer"] = layer;
                j["CalibrationType"] = "HGLGIntercalib";
                j["Summary"]["FitFailures"] = fit_failures;
                j["Summary"]["Entries"] = total_entries;
                j["Summary"]["NumUsedRun"] = run_contexts_.size();
                j["Summary"]["SameDataRuns"] = same_data_runs;
                j["Status"] = (fit_failures == 0) ? 0 : 1;

                // Per-channel arrays
                j["PerChannel"]["Intercept"] = p0_arr;
                j["PerChannel"]["Slope"] = p1_arr;
                j["PerChannel"]["FitStatus"] = fit_status_arr;
                j["PerChannel"]["NPoints"] = n_points_arr;

                // Write to file
                std::ostringstream filename;
                filename << cfg_.out_json_dirname << "/run" << rc.config.runNumber << "_intercalib_Layer" << layer << ".json";
                std::string out_filename = filename.str();

                std::ofstream jout(out_filename);
                if (!jout.is_open()) {
                    LOG_ERROR("InterCalibAlg: cannot write JSON {}", out_filename);
                    continue;
                }
                jout << j.dump(2) << std::endl;
                jout.close();

                LOG_INFO("InterCalibAlg: wrote JSON {}", out_filename);
            }
        }
        if (cfg_.output_bad_cells_json) {
            json j_bad;
            std::vector<int> bad_layer, bad_chip, bad_channel;
            for (const auto& [layer, chip, channel] : bad_fit_cells_) {
                bad_layer.push_back(layer);
                bad_chip.push_back(chip);
                bad_channel.push_back(channel);
            }
            j_bad["example_layers"] = bad_layer;
            j_bad["example_chips"] = bad_chip;
            j_bad["example_channels"] = bad_channel;
            std::ofstream jout(cfg_.out_json_dirname + "/" + cfg_.bad_cells_json_filename);
            if (!jout.is_open()) {
                LOG_ERROR("InterCalibAlg: cannot write bad cells JSON {}", cfg_.bad_cells_json_filename);
            } else {
                jout << j_bad.dump(2) << std::endl;
                jout.close();
                LOG_INFO("InterCalibAlg: wrote bad cells JSON {}", cfg_.bad_cells_json_filename);
            }
        }
    }
    void initialize_histogram() {
        const size_t n_examples = std::min({cfg_.example_layers.size(), cfg_.example_chips.size(), cfg_.example_chs.size()});
        if (n_examples == 0) return;
        if (cfg_.example_layers.size() != cfg_.example_chips.size() || cfg_.example_layers.size() != cfg_.example_chs.size()) {
            LOG_WARN("InterCalibAlg: example array sizes differ (layers={}, chips={}, chs={}), using first {} entries",
                     cfg_.example_layers.size(), cfg_.example_chips.size(), cfg_.example_chs.size(), n_examples);
        }

        h_example_.clear();
        example_fit_lines_.clear();
        h_example_.reserve(n_examples);
        example_fit_lines_.reserve(n_examples);

        for (size_t i = 0; i < n_examples; ++i) {
            h_example_.push_back(std::make_unique<TH2D>(
                Form("h_example_hglg_L%02d_C%02d_ch%02d", cfg_.example_layers[i], cfg_.example_chips[i], cfg_.example_chs[i]),
                Form("L%d C%d ch%d; LG (ped-sub) [ADC]; HG (ped-sub) [ADC]",
                        cfg_.example_layers[i], cfg_.example_chips[i], cfg_.example_chs[i]),
                3000, -100, 2900, 400, -200, 3500));
            h_example_.back()->SetDirectory(nullptr);
            example_fit_lines_.push_back(std::make_unique<TLine>());
        }
    }
    void change_run(const RunContext& new_ctx) {
        ctx_ = new_ctx;
        run_contexts_.push_back(new_ctx);
    }
    InterCalibAlgCfg cfg_;
    RunContext ctx_;
    bool written_ = false;
    bool cache_built_ = false;
    long long n_used_ = 0;
    long long n_missing_ped_ = 0;
    int n_fit_ok_ = 0;
    int n_fit_all_ = 0;
    std::vector<RunContext> run_contexts_; // for storing contexts of all runs, for later use in writing JSON with run info

    std::unordered_map<int, CalibDBIO::Pedestal> ped_map_;
    std::unordered_map<int, Acc> acc_;
    std::unordered_map<int, CellInterCalibResult> ic_cache_;
    std::vector<std::unique_ptr<TH2D>> h_example_;
    std::vector<std::unique_ptr<TLine>> example_fit_lines_;
    std::vector<std::tuple<int, int, int>> bad_fit_cells_;
};

// Define deleter *after* Impl is a complete type in this TU.
void InterCalibAlg::ImplDeleter::operator()(InterCalibAlg::Impl* p) const {
    delete p;
}

void InterCalibAlg::ensure_impl() {
    if (!impl_) {
        impl_.reset(new Impl(cfg_, ctx()));
        impl_->initialize_histogram();
    }
}

InterCalibAlg::~InterCalibAlg() {
    if (impl_) impl_->write();
}


void InterCalibAlg::execute(EventStore& evt) {
    ensure_impl();

    auto raw_hits = evt.get<std::vector<AHCALRawHit>>(cfg_.in_rawhit_key);
    for (const auto& h : raw_hits) {
        impl_->fill(h);
    }
}

void InterCalibAlg::parse_cfg(const YAML::Node& n) {
    cfg_.in_rawhit_key = get_or<std::string>(n, "in_rawhit_key", cfg_.in_rawhit_key);
    cfg_.in_pedestal_file = get_or<std::string>(n, "in_pedestal_file", cfg_.in_pedestal_file);

    cfg_.intercalib_to_file = get_or<bool>(n, "intercalib_to_file", cfg_.intercalib_to_file);
    cfg_.intercalib_to_json = get_or<bool>(n, "intercalib_to_json", cfg_.intercalib_to_json);
    cfg_.out_intercalib_filename = get_or<std::string>(n, "out_intercalib_filename",
                                                        cfg_.out_intercalib_filename);
    cfg_.out_json_dirname = get_or<std::string>(n, "out_json_dirname", cfg_.out_json_dirname);

    cfg_.hg_fit_max = get_or<double>(n, "hg_fit_max", cfg_.hg_fit_max);
    // cfg_.lg_fit_max = get_or<double>(n, "lg_fit_max", cfg_.lg_fit_max);
    cfg_.hg_fit_min = get_or<double>(n, "hg_fit_min", cfg_.hg_fit_min);
    cfg_.lg_fit_min = get_or<double>(n, "lg_fit_min", cfg_.lg_fit_min);
    cfg_.min_points = get_or<int>(n, "min_points", cfg_.min_points);
    cfg_.hg_bin_width = get_or<double>(n, "hg_bin_width", cfg_.hg_bin_width);
    cfg_.min_hg_bins_for_fit = get_or<int>(n, "min_hg_bins_for_fit", cfg_.min_hg_bins_for_fit);

    cfg_.example_layers = get_or<std::vector<int> >(n, "example_layers", cfg_.example_layers);
    cfg_.example_chips = get_or<std::vector<int> >(n, "example_chips", cfg_.example_chips);
    cfg_.example_chs = get_or<std::vector<int> >(n, "example_chs", cfg_.example_chs);

    cfg_.require_hittag = get_or<bool>(n, "require_hittag", cfg_.require_hittag);
    cfg_.require_yperx_over1 = get_or<bool>(n, "require_yperx_over1", cfg_.require_yperx_over1);
    cfg_.output_bad_cells_json = get_or<bool>(n, "output_bad_cells_json", cfg_.output_bad_cells_json);
    cfg_.bad_cells_json_filename = get_or<std::string>(n, "bad_cells_json_filename", cfg_.bad_cells_json_filename);
    
    cfg_.use_specific_fit_range_toL39C8 = get_or<bool>(n, "use_specific_fit_range_toL39C8", cfg_.use_specific_fit_range_toL39C8);
    cfg_.hg_fit_max_L39C8 = get_or<double>(n, "hg_fit_max_L39C8", cfg_.hg_fit_max_L39C8);
}
void InterCalibAlg::init_by_run() {
    LOG_INFO("InterCalibAlg: init_by_run called, runNumber={}, starttime={}, endtime={}",
             ctx().config.runNumber, ctx().conditions.starttime, ctx().conditions.endtime);
    const bool was_uninitialized = !impl_;
    ensure_impl();
    if (!was_uninitialized) {
        impl_->change_run(ctx());
        impl_->loadPedestals();
        LOG_INFO("InterCalibAlg: re-loaded pedestals for new run");
    }
}
} // namespace AHCALRecoAlg
