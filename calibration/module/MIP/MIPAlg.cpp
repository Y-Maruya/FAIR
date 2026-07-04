#include "MIPAlg.hpp"
#include "langaus.h"
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
#include <TF1.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TString.h>
#include <TParameter.h>

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
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

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
        double fit_ms = 0.0;
        double pro_ms = 0.0;
    };
    static FitOut fitLandauGaus(TH1D* h, 
                                int minEntries,
                                bool calculate_fwhm = true)
    {
        using Clock = std::chrono::steady_clock;
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
        fr[0] = 50;fr[1] = 1000;
        pllo[0] = 0;  pllo[1] = 0; pllo[2] = 100; pllo[3] = 0;
        plhi[0] = 300;plhi[1] = 1200;plhi[2] = 100000;plhi[3] = 300;
        sv[0] = 50; sv[1] = 300; sv[2] = 30000;   sv[3] = 60;

        const auto fit_t0 = Clock::now();
        TF1* fLandauGaus = langaufit(h, fr, sv, pllo, plhi, fps, fpe, &chisqr, &ndf);
        const auto fit_t1 = Clock::now();
        r.fit_ms = std::chrono::duration<double, std::milli>(fit_t1 - fit_t0).count();
        if (!fLandauGaus) return r;

        double maxx=0,fwhm=-1.0;
        if (calculate_fwhm) {
            const auto pro_t0 = Clock::now();
            langaupro(fps,maxx,fwhm);
            const auto pro_t1 = Clock::now();
            r.pro_ms = std::chrono::duration<double, std::milli>(pro_t1 - pro_t0).count();
        } else {
            r.pro_ms = 0.0;
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
        if (chi2_ndf > 30 || (x > 0 && x > fr[1] - 40) || r.gaus_sigma < 30 || r.gaus_sigma > 160 || r.width < 15 || r.width > 80)
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

    struct MIPAlg::Impl {
        struct LayerCutflow {
            long long skipped_layer = 0;
            long long out_of_active_area = 0;
            long long skipped_chip = 0;
            long long outside_xy_threshold = 0;
            long long failed_neighbor_support = 0;
            long long neighbor_cell_hit = 0;
            long long selected_crossing = 0;
            long long target_hit_fill = 0;
            long long no_exact_target_hit = 0;
        };

        struct ProjectedCell {
            double x = 0.0;
            double y = 0.0;
            int chip = -1;
            int channel = -1;
        };

        struct FitResult {
            int cellid = -1;
            int layer = -1;
            int chip = -1;
            int channel = -1;
            int channel_index = -1;
            double x_mm = -999.0;
            double y_mm = -999.0;
            int entries = 0;
            int entries_adc_le50 = 0;
            double ratio_adc_le50 = -1.0;
            FitOut out;
        };

        explicit Impl(MIPAlgCfg cfg, const RunContext& ctx)
            : cfg_(std::move(cfg)), ctx_(ctx) {
            loadPedestals();
            run_contexts_.push_back(ctx);
        }

        bool isSkippedLayer(int layer) const {
            return std::find(cfg_.mip_skip_layers.begin(), cfg_.mip_skip_layers.end(), layer)
                != cfg_.mip_skip_layers.end();
        }

        bool isSkippedChip(int layer, int chip) const {
            return std::find(cfg_.mip_skip_chips.begin(), cfg_.mip_skip_chips.end(),
                             std::make_pair(layer, chip))
                != cfg_.mip_skip_chips.end();
        }

        bool projectTrackToCell(const SimpleFittedTrack& track, int layer,
                                ProjectedCell& projected) const {
            if (layer < 0 || layer >= AHCALGeometry::Layer_No || isSkippedLayer(layer)) {
                return false;
            }
            projected.x = track.init_pos_x + track.direction_x * AHCALGeometry::Pos_Z(layer);
            projected.y = track.init_pos_y + track.direction_y * AHCALGeometry::Pos_Z(layer);
            if (std::abs(projected.x) > AHCALGeometry::x_max ||
                std::abs(projected.y) > AHCALGeometry::y_max) {
                return false;
            }
            AHCALGeometry::inverse(projected.x, projected.y, projected.chip, projected.channel);
            return !isSkippedChip(layer, projected.chip);
        }

        bool hasExactHit(const std::vector<AHCALRawHit>& rawhits, int layer,
                         int chip, int channel) const {
            return std::any_of(rawhits.begin(), rawhits.end(), [&](const AHCALRawHit& rh) {
                return rh.layer() == layer && rh.chip() == chip &&
                       rh.channel() == channel && rh.hittag == 1;
            });
        }

        std::pair<int, int> requiredSupportLayers(const SimpleFittedTrack& track,
                                                  int target_layer) const {
            int required_up = std::max(0, cfg_.mip_neighbor_upstream_layers);
            int required_down = std::max(0, cfg_.mip_neighbor_downstream_layers);
            if (!cfg_.mip_edge_layer_policy) {
                return {required_up, required_down};
            }

            auto countValidLayers = [&](int first_layer, int step) {
                int count = 0;
                for (int layer = first_layer;
                     layer >= 0 && layer < AHCALGeometry::Layer_No;
                     layer += step) {
                    ProjectedCell projected;
                    if (projectTrackToCell(track, layer, projected)) {
                        ++count;
                    }
                }
                return count;
            };

            const int required_total = required_up + required_down;
            const int available_up = countValidLayers(target_layer - 1, -1);
            const int available_down = countValidLayers(target_layer + 1, +1);

            required_up = std::min(required_up, available_up);
            required_down = required_total - required_up;
            if (required_down > available_down) {
                required_down = available_down;
                required_up = required_total - required_down;
            }
            return {required_up, required_down};
        }

        bool hasRequiredSupport(const std::vector<AHCALRawHit>& rawhits,
                                const SimpleFittedTrack& track, int target_layer,
                                int required_up, int required_down) const {
            auto directionHasSupport = [&](int first_layer, int step, int required) {
                int valid_layers = 0;
                for (int layer = first_layer;
                     layer >= 0 && layer < AHCALGeometry::Layer_No && valid_layers < required;
                     layer += step) {
                    ProjectedCell projected;
                    if (!projectTrackToCell(track, layer, projected)) {
                        continue;
                    }
                    ++valid_layers;
                    if (!hasExactHit(rawhits, layer, projected.chip, projected.channel)) {
                        return false;
                    }
                }
                return valid_layers == required;
            };

            return directionHasSupport(target_layer - 1, -1, required_up) &&
                   directionHasSupport(target_layer + 1, +1, required_down);
        }

        bool hasNeighborCellHit(const std::vector<AHCALRawHit>& rawhits, int target_layer,
                                int target_chip, int target_channel) const {
            const double target_x = AHCALGeometry::Pos_X(target_channel, target_chip);
            const double target_y = AHCALGeometry::Pos_Y(target_channel, target_chip);
            const double window = std::max(0, cfg_.mip_neighbor_cell_search_radius) * 40.3;
            return std::any_of(rawhits.begin(), rawhits.end(), [&](const AHCALRawHit& rh) {
                if (rh.hittag != 1 || rh.layer() != target_layer ||
                    (rh.chip() == target_chip && rh.channel() == target_channel)) {
                    return false;
                }
                const double hit_x = AHCALGeometry::Pos_X(rh.channel(), rh.chip());
                const double hit_y = AHCALGeometry::Pos_Y(rh.channel(), rh.chip());
                return std::abs(hit_x - target_x) <= window &&
                       std::abs(hit_y - target_y) <= window;
            });
        }

        void loadPedestals() {
            ped_map_ = std::make_shared<CalibDBIO::PedestalMap>();
            if (cfg_.read_pedestal_from_ROOT) {
                loadPedestals_fromROOT();
            } else if (cfg_.read_pedestal_from_DB) {
                loadPedestals_fromDB();
            } else {
                LOG_WARN("MIPAlg: no pedestal source specified");
            }
        }

        void loadPedestals_fromROOT() {
            std::unique_ptr<TFile> fped(TFile::Open(cfg_.in_pedestal_file.c_str(), "READ"));
            if (!fped || fped->IsZombie()) {
                LOG_ERROR("MIPAlg: cannot open pedestal file: {}", cfg_.in_pedestal_file);
                return;
            }

            TTree* tp = static_cast<TTree*>(fped->Get("pedestal"));
            if (!tp) {
                LOG_ERROR("MIPAlg: TTree 'pedestal' not found in {}", cfg_.in_pedestal_file);
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

            LOG_INFO("MIPAlg: loaded {} pedestal entries from ROOT", ped_map_->size());
        }

        void loadPedestals_fromDB() {
            CalibDBIO::PedestalReader reader(ctx_.config.runNumber);
            ped_map_ = reader.getPedestalMapPtr();
            LOG_INFO("MIPAlg: loaded {} pedestal entries from DB", ped_map_->size());
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
                    n_missing_ped_++;
                    return;
                }
                hg_value -= itp->second.HighGainPeak;
            }

            // Count ADC <= 50 regardless of xmin/xmax range
            if (hg_value <= 50.0) {
                adc_le50_count_[cellid]++;
            }

            if (hg_value >= cfg_.xmin && hg_value <= cfg_.xmax) {
                auto& ptr = hg_hist_[cellid];
                if (!ptr) {
                    ptr = make_hist(cellid);
                    if (hg_hist_.size() == 1) {
                        LOG_INFO("MIPAlg: created first histogram for cellID={}", cellid);
                    } else if (hg_hist_.size() % 500 == 0) {
                        LOG_INFO("MIPAlg: created {} histograms so far (missing_ped={})",
                                 hg_hist_.size(), n_missing_ped_);
                    }
                }
                ptr->Fill(hg_value);
            }
            if (hg_hist_.size() % 1000 == 0) {
                LOG_INFO("MIPAlg: {} histograms created, estimated memory ≈ {} MB",
                        hg_hist_.size(), hg_hist_.size() * 45 / 1024);
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
                
                // Calculate ratio of ADC <= 50
                auto adc_le50_it = adc_le50_count_.find(cid);
                if (adc_le50_it != adc_le50_count_.end()) {
                    res.entries_adc_le50 = adc_le50_it->second;
                    if (res.entries > 0) {
                        res.ratio_adc_le50 = static_cast<double>(res.entries_adc_le50) / res.entries;
                    }
                }
                
                if (hist_ptr && cfg_.fit) {
                    nFitAll_++;
                    res.out = fitLandauGaus(hist_ptr.get(), cfg_.min_entries, cfg_.calculate_fwhm);
                    if (res.out.fit_ok) nFitOK_++;
                }

                fit_cache_.emplace(cid, std::move(res));
            }

            LOG_INFO("MIPAlg: built fit cache with {}/{} successful fits", nFitOK_, nFitAll_);
        }

        void fillMapsFromCache(std::vector<std::unique_ptr<TH2D>>& hMPV,
                                std::vector<std::unique_ptr<TH2D>>& hWidth,
                                std::vector<std::unique_ptr<TH2D>>& hEntries,
                                std::vector<std::unique_ptr<TH2D>>& hTotal,
                                std::vector<std::unique_ptr<TH2D>>& hGausSigma,
                                std::vector<std::unique_ptr<TH2D>>& hRatioADCLe50) {
            for (const auto& [cid, res] : fit_cache_) {
                if (res.layer >= 0 && res.layer < AHCALGeometry::Layer_No) {
                    const int binx = hMPV[res.layer]->GetXaxis()->FindBin(res.x_mm);
                    const int biny = hMPV[res.layer]->GetYaxis()->FindBin(res.y_mm);
                    if (binx > 0 && binx <= NBIN_XY && biny > 0 && biny <= NBIN_XY) {
                        // Fill ratio map for all cells
                        hRatioADCLe50[res.layer]->SetBinContent(binx, biny, res.ratio_adc_le50);
                        
                        if (res.out.fit_ok) {
                            hMPV[res.layer]->SetBinContent(binx, biny, res.out.mpv);
                            hWidth[res.layer]->SetBinContent(binx, biny, res.out.width);
                            hEntries[res.layer]->SetBinContent(binx, biny, res.entries);
                            hTotal[res.layer]->SetBinContent(binx, biny, res.out.total_area);
                            hGausSigma[res.layer]->SetBinContent(binx, biny, res.out.gaus_sigma);
                        }
                    }
                }
            }
        }

        void writeRootOutput() {
            if (!cfg_.mip_to_file) return;
            std::filesystem::path out_path(cfg_.out_mip_filename);
            if (out_path.has_parent_path()) {
                std::error_code ec;
                std::filesystem::create_directories(out_path.parent_path(), ec);
                if (ec) {
                    LOG_ERROR("MIPAlg: cannot create output directory {}: {}", out_path.parent_path().string(), ec.message());
                    return;
                }
            }
            auto fout = std::unique_ptr<TFile>(TFile::Open(cfg_.out_mip_filename.c_str(), "RECREATE"));
            if (!fout || fout->IsZombie()) {
                LOG_ERROR("MIPAlg: cannot create output file: {}", cfg_.out_mip_filename);
                std::cerr << "MIPAlg: cannot create output file: " << cfg_.out_mip_filename << std::endl;
                return;
            }

            // directories
            TDirectory* dHist = ensureDir(fout.get(), "MIP");
            for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                if (dHist) fout->mkdir(Form("MIP/Layer%02d", L));
                for (int C = 0; C < AHCALGeometry::chip_No; ++C) {
                    if (dHist) fout->mkdir(Form("MIP/Layer%02d/Chip%d", L, C));
                }
            }
            TDirectory* dMap = ensureDir(fout.get(), "Map");
            TDirectory* dCan = ensureDir(fout.get(), "Canvases");

            // Write histograms
            if (dHist) {
                TTree * tree_ntrack = new TTree("ntrack_pass_through_channel", "Number of tracks passing through each channel");
                int cellid_ = -1;
                int ntracks_ = 0;
                tree_ntrack->Branch("cellid", &cellid_, "cellid/I");
                tree_ntrack->Branch("ntracks", &ntracks_, "ntracks/I");
                for (const auto& [cellid, h] : hg_hist_) {
                    fout->cd(Form("MIP/Layer%02d/Chip%d", cellid/100000, (cellid/10000)%10));
                    if (h) h->Write();
                    TParameter<int> pNtracks_pass_through_channel(Form("Ntracks_pass_through_channel_%d", cellid), Ntrack_pass_through_channel_[cellid]);
                    pNtracks_pass_through_channel.Write();
                }
                for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer) {
                    for (int chip = 0; chip < AHCALGeometry::chip_No; ++chip) {
                        for (int channel = 0; channel < AHCALGeometry::channel_No; ++channel) {
                            cellid_ = AHCALGeometry::CellID(layer, chip, channel);
                            ntracks_ = Ntrack_pass_through_channel_[cellid_];
                            tree_ntrack->Fill();
                        }
                    }
                }
                fout->cd("MIP");
                tree_ntrack->Write();

                TTree tree_cutflow("mip_crossing_cutflow", "MIP crossing selection cutflow by layer");
                int layer_ = -1;
                long long skipped_layer_ = 0;
                long long out_of_active_area_ = 0;
                long long skipped_chip_ = 0;
                long long outside_xy_threshold_ = 0;
                long long failed_neighbor_support_ = 0;
                long long neighbor_cell_hit_ = 0;
                long long selected_crossing_ = 0;
                long long target_hit_fill_ = 0;
                long long no_exact_target_hit_ = 0;
                tree_cutflow.Branch("layer", &layer_);
                tree_cutflow.Branch("skipped_layer", &skipped_layer_);
                tree_cutflow.Branch("out_of_active_area", &out_of_active_area_);
                tree_cutflow.Branch("skipped_chip", &skipped_chip_);
                tree_cutflow.Branch("outside_xy_threshold", &outside_xy_threshold_);
                tree_cutflow.Branch("failed_neighbor_support", &failed_neighbor_support_);
                tree_cutflow.Branch("neighbor_cell_hit", &neighbor_cell_hit_);
                tree_cutflow.Branch("selected_crossing", &selected_crossing_);
                tree_cutflow.Branch("target_hit_fill", &target_hit_fill_);
                tree_cutflow.Branch("no_exact_target_hit", &no_exact_target_hit_);
                for (layer_ = 0; layer_ < AHCALGeometry::Layer_No; ++layer_) {
                    const auto& c = layer_cutflow_[layer_];
                    skipped_layer_ = c.skipped_layer;
                    out_of_active_area_ = c.out_of_active_area;
                    skipped_chip_ = c.skipped_chip;
                    outside_xy_threshold_ = c.outside_xy_threshold;
                    failed_neighbor_support_ = c.failed_neighbor_support;
                    neighbor_cell_hit_ = c.neighbor_cell_hit;
                    selected_crossing_ = c.selected_crossing;
                    target_hit_fill_ = c.target_hit_fill;
                    no_exact_target_hit_ = c.no_exact_target_hit;
                    tree_cutflow.Fill();
                }
                tree_cutflow.Write();

            }

            fout->cd();
            if (!cfg_.fit) {
                LOG_INFO("MIPAlg: wrote {} histograms without fitting", hg_hist_.size());
                fout->Close();
                return;
            }

            // Create output maps
            std::vector<std::unique_ptr<TH2D>> hMPV(AHCALGeometry::Layer_No);
            std::vector<std::unique_ptr<TH2D>> hWidth(AHCALGeometry::Layer_No);
            std::vector<std::unique_ptr<TH2D>> hEntries(AHCALGeometry::Layer_No);
            std::vector<std::unique_ptr<TH2D>> hTotal(AHCALGeometry::Layer_No);
            std::vector<std::unique_ptr<TH2D>> hGausSigma(AHCALGeometry::Layer_No);
            std::vector<std::unique_ptr<TH2D>> hRatioADCLe50(AHCALGeometry::Layer_No);
            
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
                hRatioADCLe50[L] = makeMap("RatioADCLe50", L, "Ratio ADC <= 50");
            }

            // Fill maps from cache
            fillMapsFromCache(hMPV, hWidth, hEntries, hTotal, hGausSigma, hRatioADCLe50);

            // Build and fill tree
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
            int entries_adc_le50 = 0;
            double ratio_adc_le50 = -1.0;
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
            tp.Branch("entries_adc_le50", &entries_adc_le50);
            tp.Branch("ratio_adc_le50", &ratio_adc_le50);

            for (const auto& [cid, res] : fit_cache_) {
                cellid = res.cellid;
                mpv = res.out.mpv;
                width = res.out.width;
                total_area = res.out.total_area;
                gaus_sigma = res.out.gaus_sigma;
                chi2 = res.out.chi2;
                ndf = res.out.ndf;
                max_x = res.out.max_x;
                FWHM = res.out.FWHM;
                if (ndf > 0) chi2perndf = res.out.chi2 / double(ndf);
                fit_status = res.out.fit_status;
                fit_ok = res.out.fit_ok ? 1 : 0;
                entries = res.entries;
                x_mm = res.x_mm;
                y_mm = res.y_mm;
                entries_adc_le50 = res.entries_adc_le50;
                ratio_adc_le50 = res.ratio_adc_le50;
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
            if (cfg_.output_to_png) {
                std::filesystem::path out_dir(cfg_.out_png_dir.empty() ? "." : cfg_.out_png_dir);
                std::error_code ec;
                std::filesystem::create_directories(out_dir, ec);
                if (ec) {
                    LOG_ERROR("MIPAlg: cannot create output directory {}: {}", out_dir.string(), ec.message());
                    return;
                }
                cAllMPV->SaveAs((cfg_.out_png_dir + "/MIP_MPV_AllLayers.png").c_str());
                for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                    auto c = std::make_unique<TCanvas>(Form("cMPV_L%02d", L), Form("MIP MPV Layer %d", L), 800, 700);
                    gPad->SetMargin(0.12, 0.14, 0.12, 0.10);
                    hMPV[L]->Draw("COLZ");
                    drawLayerLabel(L);
                    c->SaveAs((cfg_.out_png_dir + Form("/MIP_MPV_Layer%02d.png", L)).c_str());
                }
            }

            if (dMap) dMap->cd();
            for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
                hMPV[L]->Write();
                hWidth[L]->Write();
                hEntries[L]->Write();
                hTotal[L]->Write();
                hGausSigma[L]->Write();
                hRatioADCLe50[L]->Write();
            }
            if (dCan) dCan->cd();
            cAllMPV->Write();

            fout->cd();
            tp.Write();
            fout->Close();

            LOG_INFO("MIPAlg: wrote {}", cfg_.out_mip_filename);
        }

        void writeJsonOutput() {
            if (!cfg_.mip_to_json) return;

            const int n_channels_per_layer = AHCALGeometry::chip_No * AHCALGeometry::channel_No;
            std::filesystem::path out_dir(cfg_.out_json_dirname.empty() ? "." : cfg_.out_json_dirname);
            std::error_code ec;
            std::filesystem::create_directories(out_dir, ec);
            if (ec) {
                LOG_ERROR("MIPAlg: cannot create output directory {}: {}", out_dir.string(), ec.message());
                return;
            }

            const std::string timestamp = formatTimestampFromRunContexts(run_contexts_);
            std::vector<int> same_data_runs;
            for (const auto& rc : run_contexts_) {
                same_data_runs.push_back(rc.config.runNumber);
            }

            for (const auto& rc : run_contexts_) {
                for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer) {
                    std::vector<double> mpv_arr(n_channels_per_layer, -1.0);
                    std::vector<double> width_arr(n_channels_per_layer, -1.0);
                    std::vector<double> total_area_arr(n_channels_per_layer, -1.0);
                    std::vector<double> gaus_sigma_arr(n_channels_per_layer, -1.0);
                    std::vector<double> max_x_arr(n_channels_per_layer, -1.0);
                    std::vector<double> FWHM_arr(n_channels_per_layer, -1.0);
                    std::vector<int> entries_arr(n_channels_per_layer, 0);
                    std::vector<int> entries_adc_le50_arr(n_channels_per_layer, 0);
                    std::vector<double> ratio_adc_le50_arr(n_channels_per_layer, -1.0);
                    std::vector<double> chi2_arr(n_channels_per_layer, -1.0);
                    std::vector<int> ndf_arr(n_channels_per_layer, 0);
                    std::vector<int> fit_status_arr(n_channels_per_layer, 999);
                    std::vector<int> fit_ok_arr(n_channels_per_layer, 0);

                    long long total_entries = 0;
                    int fit_failures = 0;

                    for (const auto& [_, res] : fit_cache_) {
                        if (res.layer != layer) continue;
                        const int idx = res.chip * AHCALGeometry::channel_No + res.channel;
                        if (idx < 0 || idx >= n_channels_per_layer) continue;

                        mpv_arr[idx] = res.out.mpv;
                        width_arr[idx] = res.out.width;
                        total_area_arr[idx] = res.out.total_area;
                        gaus_sigma_arr[idx] = res.out.gaus_sigma;
                        max_x_arr[idx] = res.out.max_x;
                        FWHM_arr[idx] = res.out.FWHM;
                        entries_arr[idx] = res.entries;
                        entries_adc_le50_arr[idx] = res.entries_adc_le50;
                        ratio_adc_le50_arr[idx] = res.ratio_adc_le50;
                        chi2_arr[idx] = res.out.chi2;
                        ndf_arr[idx] = res.out.ndf;
                        fit_status_arr[idx] = res.out.fit_status;
                        fit_ok_arr[idx] = res.out.fit_ok ? 1 : 0;
                        total_entries += res.entries;
                        if (!res.out.fit_ok) {
                            fit_failures++;
                        }
                    }

                    json j;
                    j["RunNumber"] = rc.config.runNumber;
                    j["TimeStamp"] = timestamp;
                    j["Layer"] = layer;
                    j["CalibrationType"] = "MIP";
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
                    j["PerChannel"]["FWHM"] = FWHM_arr;
                    j["PerChannel"]["Entries"] = entries_arr;
                    j["PerChannel"]["EntriesADCLe50"] = entries_adc_le50_arr;
                    j["PerChannel"]["RatioADCLe50"] = ratio_adc_le50_arr;
                    j["PerChannel"]["Chi2"] = chi2_arr;
                    j["PerChannel"]["NDF"] = ndf_arr;
                    j["PerChannel"]["FitStatus"] = fit_status_arr;
                    j["PerChannel"]["FitOK"] = fit_ok_arr;

                    std::ostringstream filename;
                    if (!out_dir.empty()) {
                        filename << out_dir.string() << "/";
                    }
                    filename << "run" << rc.config.runNumber << "_mip_Layer" << layer << ".json";

                    const std::string out_filename = filename.str();
                    std::ofstream jout(out_filename);
                    if (!jout.is_open()) {
                        LOG_ERROR("MIPAlg: cannot write JSON {}", out_filename);
                        continue;
                    }
                    jout << j.dump(2) << std::endl;
                    jout.close();
                    LOG_INFO("MIPAlg: wrote JSON {}", out_filename);
                }
            }
        }

        void write() {
            if (!cfg_.mip_to_file && !cfg_.mip_to_json) return;
            if (written_) return;
            written_ = true;

            buildFitCache();

            if (cfg_.mip_to_file) {
                writeRootOutput();
            }
            if (cfg_.mip_to_json) {
                writeJsonOutput();
            }
        }

        std::unique_ptr<TH1D> make_hist(int cellid){
            int layer = cellid/100000;
            int chip = cellid/10000 % 10;
            const std::string name = "hMIP_" + std::to_string(cellid);
            const std::string title = "Layer " + std::to_string(layer) + " chip "+ std::to_string(chip) + " MIP;ADC (ped-sub);counts";
            auto h = std::make_unique<TH1D>(name.c_str(), title.c_str(), cfg_.nbin, cfg_.xmin, cfg_.xmax);
            h->SetDirectory(nullptr);
            return h;
        }

        void initialize(){
            Ntrack_pass_through_channel_.clear();
            for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer) {
                for (int chip = 0; chip < AHCALGeometry::chip_No; ++chip) {
                    for (int channel = 0; channel < AHCALGeometry::channel_No; ++channel) {
                        int cellid = layer * 100000 + chip * 10000 + channel;
                        Ntrack_pass_through_channel_[cellid] = 0;
                    }
                }
            }
        }
        MIPAlgCfg cfg_;
        RunContext ctx_;
        bool written_ = false;
        bool cache_built_ = false;
        int nFitAll_ = 0;
        int nFitOK_ = 0;
        long long n_missing_ped_ = 0;
        std::vector<RunContext> run_contexts_;
        std::shared_ptr<const CalibDBIO::PedestalMap> ped_map_ = std::make_shared<CalibDBIO::PedestalMap>();
        std::unordered_map<int, std::unique_ptr<TH1D>> hg_hist_;
        std::unordered_map<int, int> adc_le50_count_;
        std::unordered_map<int, FitResult> fit_cache_;
        std::unordered_map<int, int> Ntrack_pass_through_channel_;
        LayerCutflow layer_cutflow_[AHCALGeometry::Layer_No];
    };

    void MIPAlg::ImplDeleter::operator()(Impl* p) const {
        delete p;
    }

    void MIPAlg::ensure_impl() {
        if (!impl_) {
            impl_.reset(new Impl(cfg_, ctx()));
        }
    }

    MIPAlg::~MIPAlg(){
        if (impl_) impl_->write();
    }

    void MIPAlg::init_by_run() {
        const bool was_uninitialized = !impl_;
        ensure_impl();
        if (!was_uninitialized) {
            impl_->change_run(ctx());
            impl_->loadPedestals();
            LOG_INFO("MIPAlg: re-loaded pedestals for new run {}", ctx().config.runNumber);
        }
    }

    void MIPAlg::execute(EventStore& evt){
        ensure_impl();
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
            for (int i_layer = 0; i_layer < AHCALGeometry::Layer_No; ++i_layer) {
                auto& cutflow = impl_->layer_cutflow_[i_layer];
                if (impl_->isSkippedLayer(i_layer)) {
                    ++cutflow.skipped_layer;
                    LOG_DEBUG("MIPAlg: target layer {} skipped", i_layer);
                    continue;
                }

                const double x = track.init_pos_x + track.direction_x * AHCALGeometry::Pos_Z(i_layer);
                const double y = track.init_pos_y + track.direction_y * AHCALGeometry::Pos_Z(i_layer);
                if (std::abs(x) > AHCALGeometry::x_max || std::abs(y) > AHCALGeometry::y_max) {
                    ++cutflow.out_of_active_area;
                    LOG_DEBUG("MIPAlg: track extrapolation out of active area for layer {}: x={}, y={}",
                              i_layer, x, y);
                    continue;
                }

                int chip = -1;
                int channel = -1;
                AHCALGeometry::inverse(x, y, chip, channel);
                if (impl_->isSkippedChip(i_layer, chip)) {
                    ++cutflow.skipped_chip;
                    LOG_DEBUG("MIPAlg: predicted chip skipped for layer {} chip {}", i_layer, chip);
                    continue;
                }

                const double rh_x = AHCALGeometry::Pos_X(channel, chip);
                const double rh_y = AHCALGeometry::Pos_Y(channel, chip);
                if (std::abs(x - rh_x) > cfg_.xy_size_threshold ||
                    std::abs(y - rh_y) > cfg_.xy_size_threshold) {
                    ++cutflow.outside_xy_threshold;
                    LOG_DEBUG("MIPAlg: track extrapolation outside xy threshold for layer {}", i_layer);
                    continue;
                }

                if (cfg_.enable_neighbor_layer_crossing_selection) {
                    const auto [required_up, required_down] =
                        impl_->requiredSupportLayers(track, i_layer);
                    if (!impl_->hasRequiredSupport(rawhits, track, i_layer, required_up, required_down)) {
                        ++cutflow.failed_neighbor_support;
                        LOG_DEBUG("MIPAlg: failed neighboring-layer support for layer {}", i_layer);
                        continue;
                    }
                }

                if (cfg_.mip_reject_if_neighbor_cell_hit &&
                    impl_->hasNeighborCellHit(rawhits, i_layer, chip, channel)) {
                    ++cutflow.neighbor_cell_hit;
                    LOG_DEBUG("MIPAlg: rejected neighboring-cell hit for layer {}", i_layer);
                    continue;
                }

                const int cellid = AHCALGeometry::CellID(i_layer, chip, channel);
                ++impl_->Ntrack_pass_through_channel_[cellid];
                ++cutflow.selected_crossing;

                auto hit_it = std::find_if(rawhits.begin(), rawhits.end(), [&](const AHCALRawHit& rh) {
                    return rh.layer() == i_layer && rh.chip() == chip &&
                           rh.channel() == channel && rh.hittag == 1;
                });
                if (hit_it != rawhits.end()) {
                    impl_->fill(*hit_it);
                    ++cutflow.target_hit_fill;
                } else {
                    ++cutflow.no_exact_target_hit;
                    LOG_DEBUG("MIPAlg: no exact target hit for layer {} (x={}, y={})", i_layer, x, y);
                }
            }
            // if (!track.inTrackHitsIndices.size()) {
            //     LOG_WARN("MIPAlg: track has no associated hits.");
            //     return;
            // }
            // for (const int index : track.inTrackHitsIndices) {
            //     if (index < 0 || index >= static_cast<int>(rawhits.size())) {
            //         LOG_WARN("MIPAlg: track hit index {} out of range [0, {})", index, rawhits.size());
            //         continue;
            //     }
            //     const auto& rh = rawhits.at(index);
            //     if (rh.index != index) {
            //         LOG_WARN("MIPAlg: rawhit index {} does not match track hit index {}", rh.index, index);
            //         continue;
            //     }
            //     impl_->fill(rh);
            // }    
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
        // cfg_.mip_to_DB = get_or<bool>(cfg, "mip_to_DB", cfg_.mip_to_DB);
        cfg_.out_mip_filename = get_or<std::string>(cfg, "out_mip_filename", cfg_.out_mip_filename);
        cfg_.output_to_png = get_or<bool>(cfg, "output_to_png", cfg_.output_to_png);
        cfg_.out_png_dir = get_or<std::string>(cfg, "out_png_dir", cfg_.out_png_dir);
        cfg_.fit = get_or<bool>(cfg, "fit", cfg_.fit);
        cfg_.nbin = get_or<int>(cfg, "nbin", cfg_.nbin);
        cfg_.xmin = get_or<double>(cfg, "xmin", cfg_.xmin);
        cfg_.xmax = get_or<double>(cfg, "xmax", cfg_.xmax);
        cfg_.min_entries = get_or<int>(cfg, "min_entries", cfg_.min_entries);
        cfg_.mip_to_json = get_or<bool>(cfg, "mip_to_json", cfg_.mip_to_json);
        cfg_.out_json_dirname = get_or<std::string>(cfg, "out_json_dirname", cfg_.out_json_dirname);
        const std::string legacy_out_json_filename = get_or<std::string>(cfg, "out_json_filename", "");
        if (!legacy_out_json_filename.empty()) {
            std::filesystem::path legacy_path(legacy_out_json_filename);
            const std::string legacy_dir = legacy_path.parent_path().string();
            cfg_.out_json_dirname = legacy_dir.empty() ? "." : legacy_dir;
            LOG_WARN("MIPAlg: 'out_json_filename' is deprecated. Use 'out_json_dirname' instead.");
        }
        cfg_.calculate_fwhm = get_or<bool>(cfg, "calculate_fwhm", cfg_.calculate_fwhm);
        cfg_.substrate_pedestal = get_or<bool>(cfg, "substrate_pedestal", cfg_.substrate_pedestal);
        cfg_.read_pedestal_from_ROOT = get_or<bool>(cfg, "read_pedestal_from_ROOT", cfg_.read_pedestal_from_ROOT);
        cfg_.in_pedestal_file = get_or<std::string>(cfg, "in_pedestal_file", cfg_.in_pedestal_file);
        cfg_.read_pedestal_from_DB = get_or<bool>(cfg, "read_pedestal_from_DB", cfg_.read_pedestal_from_DB);
        cfg_.xy_size_threshold = get_or<double>(cfg, "xy_size_threshold", cfg_.xy_size_threshold);
        cfg_.enable_neighbor_layer_crossing_selection = get_or<bool>(
            cfg, "enable_neighbor_layer_crossing_selection", cfg_.enable_neighbor_layer_crossing_selection);
        cfg_.mip_neighbor_upstream_layers = get_or<int>(
            cfg, "mip_neighbor_upstream_layers", cfg_.mip_neighbor_upstream_layers);
        cfg_.mip_neighbor_downstream_layers = get_or<int>(
            cfg, "mip_neighbor_downstream_layers", cfg_.mip_neighbor_downstream_layers);
        cfg_.mip_edge_layer_policy = get_or<bool>(
            cfg, "mip_edge_layer_policy", cfg_.mip_edge_layer_policy);
        cfg_.mip_reject_if_neighbor_cell_hit = get_or<bool>(
            cfg, "mip_reject_if_neighbor_cell_hit", cfg_.mip_reject_if_neighbor_cell_hit);
        cfg_.mip_neighbor_cell_search_radius = get_or<int>(
            cfg, "mip_neighbor_cell_search_radius", cfg_.mip_neighbor_cell_search_radius);
        cfg_.mip_skip_layers = get_or<std::vector<int>>(cfg, "mip_skip_layers", cfg_.mip_skip_layers);
        cfg_.mip_skip_chips.clear();
        if (cfg["mip_skip_chips"]) {
            for (const auto& entry : cfg["mip_skip_chips"]) {
                if (!entry.IsSequence() || entry.size() != 2) {
                    throw std::runtime_error("MIPAlg: each mip_skip_chips entry must be [layer, chip]");
                }
                cfg_.mip_skip_chips.emplace_back(entry[0].as<int>(), entry[1].as<int>());
            }
        }
        if (cfg_.mip_neighbor_upstream_layers < 0 || cfg_.mip_neighbor_downstream_layers < 0 ||
            cfg_.mip_neighbor_cell_search_radius < 0) {
            throw std::runtime_error("MIPAlg: neighbor layer counts and cell search radius must be non-negative");
        }
    }
} // namespace AHCALRecoAlg
