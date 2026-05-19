#include "PedestalAlg.hpp"

#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"
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
AHCAL_REGISTER_ALG(AHCALRecoAlg::PedestalAlg, "PedestalAlg")
namespace AHCALRecoAlg {

namespace {

// Use geometry-defined full extent for maps
const double XYMIN = -AHCALGeometry::x_max;
const double XYMAX = +AHCALGeometry::x_max;
constexpr int NBIN_XY = 18;

struct FitOut {
  double mean  = -1.0;
  double sigma = -1.0;
  double rms = -1.0;
  double avg = -1.0;
  int    fitStatus = -999;
  int   status = -999;
  double chi2 = -1.0;
  int   ndf = -1;
  double tail_3sig_fraction = -1.0;
  double tail_5sig_fraction = -1.0;
  double deviation_squared_sum = -1.0;
  double deviation_squared_sum_fourth = -1.0;
  bool   ok = false;
};

static inline double peakFromMaxBin(TH1D* h) {
  const int b = h->GetMaximumBin();
  return h->GetBinCenter(b);
}

static FitOut fitPedestalGaussian(TH1D* h,
                                  int minEntries,
                                  double nsigma1,
                                  double nsigma2,
                                  double sigmaMin,
                                  double sigmaMax) {
  FitOut r;
  if (!h) return r;
  if (h->GetEntries() < minEntries) return r;

  double mu0 = peakFromMaxBin(h);
  double rms = h->GetRMS();
  r.rms = rms;
  r.avg = h->GetMean();
  if (!(rms > 0.0)) rms = 10.0;

  double sig0 = std::clamp(rms, sigmaMin, sigmaMax);

  double x1 = mu0 - nsigma1 * sig0;
  double x2 = mu0 + nsigma1 * sig0;
  x1 = std::max(x1, h->GetXaxis()->GetXmin());
  x2 = std::min(x2, h->GetXaxis()->GetXmax());
  if (x2 <= x1) return r;

  static long long ifit = 0;
  const std::string f1name = std::string("f1_") + h->GetName() + "_" + std::to_string(ifit);
  const std::string f2name = std::string("f2_") + h->GetName() + "_" + std::to_string(ifit);
  ++ifit;

  TF1 f1(f1name.c_str(), "gaus", x1, x2);
  f1.SetParameters(h->GetMaximum(), mu0, sig0);
  f1.SetParLimits(2, sigmaMin, sigmaMax);

  const int st1 = h->Fit(&f1, "RQ0");
  const double mu1 = (st1 == 0 ? f1.GetParameter(1) : mu0);
  double sg1 = (st1 == 0 ? f1.GetParameter(2) : sig0);
  sg1 = std::clamp(std::fabs(sg1), sigmaMin, sigmaMax);

  double y1 = mu1 - nsigma2 * sg1;
  double y2 = mu1 + nsigma2 * sg1;
  y1 = std::max(y1, h->GetXaxis()->GetXmin());
  y2 = std::min(y2, h->GetXaxis()->GetXmax());
  if (y2 <= y1) return r;

  TF1 f2(f2name.c_str(), "gaus", y1, y2);
  f2.SetParameters(h->GetMaximum(), mu1, sg1);
  f2.SetParLimits(2, sigmaMin, sigmaMax);

  const int st2 = h->Fit(&f2, "RQ0");
  if (st2 != 0) {
    r.mean = mu1;
    r.sigma = sg1;
    r.fitStatus = st2;
    r.ok = false;
    r.chi2 = f2.GetChisquare();
    r.ndf = f2.GetNDF();
    return r;
  }
  double sig3_upper = f2.GetParameter(1) + 3.0 * std::fabs(f2.GetParameter(2));
  double sig5_upper = f2.GetParameter(1) + 5.0 * std::fabs(f2.GetParameter(2));
  double sig3_lower = f2.GetParameter(1) - 3.0 * std::fabs(f2.GetParameter(2));
  double sig5_lower = f2.GetParameter(1) - 5.0 * std::fabs(f2.GetParameter(2));
  int entries = static_cast<int>(h->GetEntries());
  int tail_3sig = h->Integral(h->FindBin(sig3_upper), h->GetNbinsX()) + h->Integral(1, h->FindBin(sig3_lower));
  int tail_5sig = h->Integral(h->FindBin(sig5_upper), h->GetNbinsX()) + h->Integral(1, h->FindBin(sig5_lower));
  r.tail_3sig_fraction = (entries > 0) ? static_cast<double>(tail_3sig) / entries : -1.0;
  r.tail_5sig_fraction = (entries > 0) ? static_cast<double>(tail_5sig) / entries : -1.0;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    double x = h->GetBinCenter(b);
    double dev = (f2.GetParameter(2) > 0.0 ? (x - f2.GetParameter(1))/f2.GetParameter(2) : -1.0);
    dev *= h->GetBinContent(b);
    if (dev > 0) {
      r.deviation_squared_sum += dev * dev;
      r.deviation_squared_sum_fourth += dev * dev * dev * dev;
    }
  }
  r.chi2 = f2.GetChisquare();
  r.ndf = f2.GetNDF();
  r.deviation_squared_sum = (entries > 0) ? r.deviation_squared_sum / entries : -1.0;
  r.deviation_squared_sum_fourth = (entries > 0) ? r.deviation_squared_sum_fourth / entries : -1.0;
  r.mean = f2.GetParameter(1);
  r.sigma = std::fabs(f2.GetParameter(2));
  r.fitStatus = st2;
  r.ok = true;
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

} // namespace

struct PedestalAlg::Impl {
  explicit Impl(PedestalAlgCfg cfg, const RunContext& ctx) 
    : cfg_(std::move(cfg)), ctx_(ctx) {}

  struct CellPedestalResult {
    int cellid = -1;
    int layer = -1;
    int chip = -1;
    int channel = -1;
    int channel_index = -1;
    double x_mm = -999.0;
    double y_mm = -999.0;

    int entries_hg = 0;
    int entries_lg = 0;
    FitOut fit_hg;
    FitOut fit_lg;
  };

  void fill(const AHCALRawHit& h) {
    if (cfg_.use_hittag) {
      if (h.hittag != cfg_.select_hittag) return;
    }

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

  void buildPedestalCache() {
    if (cache_built_) return;
    cache_built_ = true;

    // union of keys (HG/LG)
    std::unordered_set<int> keys;
    keys.reserve(hg_hist_.size() + lg_hist_.size());
    for (const auto& [k, _] : hg_hist_) keys.insert(k);
    for (const auto& [k, _] : lg_hist_) keys.insert(k);

    ped_cache_.clear();
    ped_cache_.reserve(keys.size());

    nFitOK_HG_ = 0;
    nFitAll_HG_ = 0;
    nFitOK_LG_ = 0;
    nFitAll_LG_ = 0;

    for (int cid : keys) {
      CellPedestalResult res;
      res.cellid = cid;
      int tmp_layer, tmp_chip, tmp_channel;
      AHCALGeometry::CellIDToLC(cid, tmp_layer, tmp_chip, tmp_channel);
      res.layer = tmp_layer;
      res.chip = tmp_chip;
      res.channel = tmp_channel;
      res.channel_index = AHCALGeometry::channel_index_from_lcc(res.layer, res.chip, res.channel);
      cellid_to_xy(res.chip, res.channel, res.x_mm, res.y_mm);

      TH1D* hhg = nullptr;
      TH1D* hlg = nullptr;
      if (auto it = hg_hist_.find(cid); it != hg_hist_.end()) hhg = it->second.get();
      if (auto it = lg_hist_.find(cid); it != lg_hist_.end()) hlg = it->second.get();

      res.entries_hg = hhg ? static_cast<int>(hhg->GetEntries()) : 0;
      if (hhg) {
        nFitAll_HG_++;
        res.fit_hg = fitPedestalGaussian(hhg, cfg_.min_entries, cfg_.nsigma_win1, cfg_.nsigma_win2,
                                         cfg_.sigma_min, cfg_.sigma_max);
        int status_hg = res.fit_hg.fitStatus;
        if (!res.fit_hg.ok) {
          status_hg = res.fit_hg.fitStatus; // keep original fit status for failed fits
        }
        if (res.fit_hg.ok && cfg_.hg_cut_failed) {
          if (cfg_.hg_cut_rms_per_sigma) {
            double rms_per_sigma = (res.fit_hg.sigma > 0.0) ? res.fit_hg.rms / res.fit_hg.sigma : -1.0;
            if (rms_per_sigma > cfg_.hg_max_rms_per_sigma) {
              status_hg = AHCALRefValues::HGPedestal_cut_by_rms_per_sigma; // custom status for RMS/sigma cut
            }
          }
          if (cfg_.hg_cut_max_sigma_and_max_rms) {
            if (status_hg == AHCALRefValues::HGPedestal_cut_by_rms_per_sigma) {
              // If already failed RMS/sigma cut, we can directly set the combined cut status if it also fails max sigma AND max rms
              if (res.fit_hg.sigma > cfg_.hg_max_sigma && res.fit_hg.rms > cfg_.hg_max_rms) {
                status_hg = AHCALRefValues::HGPedestal_cut_by_both_rms_per_sigma_and_sigma_vs_chi2; // custom status for failing both cuts
              }
            } else {
              // If not already failed RMS/sigma cut, check max sigma AND max rms cut separately
              if (res.fit_hg.sigma > cfg_.hg_max_sigma && res.fit_hg.rms > cfg_.hg_max_rms) {
                status_hg = AHCALRefValues::HGPedestal_cut_by_max_sigma_and_max_rms; // custom status for max sigma AND max rms cut
              }
            }
          }
          if (cfg_.hg_cut_sigma_vs_chi2) {
            // Exclude fits with sigma vs chi2 outside envelope defined by line from (0, sigma_intercept) to (chi2_ndf_intercept, 0)
            double chi2_ndf = (res.fit_hg.ndf > 0) ? res.fit_hg.chi2 / res.fit_hg.ndf : -1.0;
            double sigma_cut = (cfg_.hg_chi2_ndf_intercept > 0.0) ?
              std::max(0.0, cfg_.hg_sigma_intercept * (1.0 - chi2_ndf / cfg_.hg_chi2_ndf_intercept)) : -1.0;
            if (res.fit_hg.sigma > sigma_cut) {
              if (status_hg == AHCALRefValues::HGPedestal_cut_by_rms_per_sigma){
                status_hg = AHCALRefValues::HGPedestal_cut_by_both_rms_per_sigma_and_sigma_vs_chi2; // custom status for failing both cuts
              } else {
                status_hg = AHCALRefValues::HGPedestal_cut_by_sigma_vs_chi2; // custom status for sigma vs chi2 cut
              }
            }
          }
        }
        res.fit_hg.status = status_hg; // update status with cut results
        if (status_hg!= AHCALRefValues::HGPedestal_OK) {
          LOG_WARN("CellID {}: HG pedestal fit failed or cut: fitStatus={}, status={}", cid, res.fit_hg.fitStatus, status_hg);
        }
        if (res.fit_hg.ok && AHCALRefValues::HGPedestalStatus_is_ok(status_hg)) nFitOK_HG_++;
      }

      res.entries_lg = hlg ? static_cast<int>(hlg->GetEntries()) : 0;
      if (hlg) {
        nFitAll_LG_++;
        res.fit_lg = fitPedestalGaussian(hlg, cfg_.min_entries, cfg_.nsigma_win1, cfg_.nsigma_win2,
                                         cfg_.sigma_min, cfg_.sigma_max);
        int status_lg = res.fit_lg.fitStatus;
        if (!res.fit_lg.ok) {
          status_lg = res.fit_lg.fitStatus; // keep original fit status for failed fits
        }
        if (res.fit_lg.ok && cfg_.lg_cut_failed) {
          if (cfg_.lg_cut_max_sigma) {
            if (res.fit_lg.sigma > cfg_.lg_max_sigma) {
              status_lg = AHCALRefValues::LGPedestal_cut_by_max_sigma; // custom status for max sigma cut
            }
          }
        }
        res.fit_lg.status = status_lg; // update status with cut results
        if (AHCALRefValues::LGPedestalStatus_is_ok(status_lg)) nFitOK_LG_++;
      }

      ped_cache_.emplace(cid, std::move(res));
    }
  }

  void write() {
    if (!cfg_.pedestal_to_file && !cfg_.pedestal_to_json) return;
    if (written_) return;
    written_ = true;

    buildPedestalCache();

    if (cfg_.pedestal_to_file) {
      auto fout = std::unique_ptr<TFile>(TFile::Open(cfg_.out_pedestal_filename.c_str(), "RECREATE"));
      if (!fout || fout->IsZombie()) {
        LOG_ERROR("PedestalAlg: cannot create output file: {}", cfg_.out_pedestal_filename);
      } else {
        // directories
        TDirectory* dHist = ensureDir(fout.get(), "Pedestal");
        TDirectory* dHG   = dHist ? ensureDir(dHist, "HG") : nullptr;
        TDirectory* dLG   = dHist ? ensureDir(dHist, "LG") : nullptr;

        TDirectory* dMap  = ensureDir(fout.get(), "PedMap2D");
        TDirectory* dCan  = ensureDir(fout.get(), "Canvases");

        // per-layer maps (HG/LG): fixed Layer_No from geometry
        std::vector<std::unique_ptr<TH2D>> hMeanHG(AHCALGeometry::Layer_No), hSigHG(AHCALGeometry::Layer_No), hEntHG(AHCALGeometry::Layer_No);
        std::vector<std::unique_ptr<TH2D>> hMeanLG(AHCALGeometry::Layer_No), hSigLG(AHCALGeometry::Layer_No), hEntLG(AHCALGeometry::Layer_No);
        std::unique_ptr<TH1D> hstatusHG = std::make_unique<TH1D>("hstatusHG", "HG pedestal fit status;Status;Entries", 7, -999.5, -992.5);
        hstatusHG->SetDirectory(nullptr);
        std::unique_ptr<TH1D> hstatusLG = std::make_unique<TH1D>("hstatusLG", "LG pedestal fit status;Status;Entries", 3, -999.5, -996.5);
        hstatusLG->SetDirectory(nullptr);
        hstatusHG->GetXaxis()->SetBinLabel( 1, "less entries");
        hstatusHG->GetXaxis()->SetBinLabel( 7, "fit failed");
        hstatusHG->GetXaxis()->SetBinLabel( 2, "high RMS/sigma");
        hstatusHG->GetXaxis()->SetBinLabel( 3, "sigma vs chi2");
        hstatusHG->GetXaxis()->SetBinLabel( 4, "both RMS/sigma and sigma vs chi2");
        hstatusHG->GetXaxis()->SetBinLabel( 5, "max sigma AND max RMS");
        hstatusHG->GetXaxis()->SetBinLabel( 6, "both max sigma AND max RMS AND RMS/sigma");

        hstatusLG->GetXaxis()->SetBinLabel( 1, "OK");
        hstatusLG->GetXaxis()->SetBinLabel( 3, "fit failed");
        hstatusLG->GetXaxis()->SetBinLabel( 2, "max sigma cut");
        std::unique_ptr<TH2D> h2_rms_sigma_HG = std::make_unique<TH2D>("h2_rms_sigma_HG", "RMS vs sigma for HG pedestal fits;Sigma;RMS", 100, 0.0, 50.0, 100, 0.0, 200.0);
        std::unique_ptr<TH2D> h2_rms_sigma_LG = std::make_unique<TH2D>("h2_rms_sigma_LG", "RMS vs sigma for LG pedestal fits;Sigma;RMS", 100, 0.0, 10.0, 100, 0.0, 40.0);
        h2_rms_sigma_HG->SetDirectory(nullptr);
        h2_rms_sigma_LG->SetDirectory(nullptr);
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

        for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
          hMeanHG[L] = makeMap("hPedMean2_HG", L, "Pedestal mean map (HG)");
          hSigHG[L]  = makeMap("hPedSigma2_HG", L, "Pedestal sigma map (HG)");
          hEntHG[L]  = makeMap("hPedEntries2_HG", L, "Pedestal entries map (HG)");

          hMeanLG[L] = makeMap("hPedMean2_LG", L, "Pedestal mean map (LG)");
          hSigLG[L]  = makeMap("hPedSigma2_LG", L, "Pedestal sigma map (LG)");
          hEntLG[L]  = makeMap("hPedEntries2_LG", L, "Pedestal entries map (LG)");
        }

        // tree
        TTree tp("pedestal", "Pedestal from RawHits (Gaussian fit)");

        int cellid=-1;
        int channel_index=-1;
        int entries_hg=0, entries_lg=0;
        int fitStatus_hg=999, fitStatus_lg=999;
        int fitOk_hg=0, fitOk_lg=0;
        int usable_hg=0, usable_lg=0;
        int status_hg=-999, status_lg=-999;

        // For AdcToEnergyReadTTreeAlg::initialize_ped()
        double highgain_peak=-1.0, lowgain_peak=-1.0;
        double highgain_rms=-1.0, lowgain_rms=-1.0;
        double highgain_avg=-1.0, lowgain_avg=-1.0;

        // Extra info
        double highgain_sigma=-1.0, lowgain_sigma=-1.0;
        double x_mm=-999.0, y_mm=-999.0;

        double highgain_chi2=-1.0, lowgain_chi2=-1.0;
        int highgain_ndf=-1, lowgain_ndf=-1;
        double highgain_chi2_ndf=-1.0, lowgain_chi2_ndf=-1.0;

        double highgain_tail_3sig_fraction=-1.0, lowgain_tail_3sig_fraction=-1.0;
        double highgain_tail_5sig_fraction=-1.0, lowgain_tail_5sig_fraction=-1.0;
        // double highgain_deviation_squared_sum=-1.0, lowgain_deviation_squared_sum=-1.0;
        // double highgain_deviation_quared_sum=-1.0, lowgain_deviation_quared_sum=-1.0;
        tp.Branch("cellid", &cellid, "cellid/I");
        tp.Branch("channel_index", &channel_index, "channel_index/I");
        tp.Branch("highgain_rms", &highgain_rms, "highgain_rms/D");
        tp.Branch("lowgain_rms", &lowgain_rms, "lowgain_rms/D");
        tp.Branch("highgain_avg", &highgain_avg, "highgain_avg/D");
        tp.Branch("lowgain_avg", &lowgain_avg, "lowgain_avg/D");
        tp.Branch("highgain_peak", &highgain_peak, "highgain_peak/D");
        tp.Branch("lowgain_peak", &lowgain_peak, "lowgain_peak/D");
        tp.Branch("highgain_sigma", &highgain_sigma, "highgain_sigma/D");
        tp.Branch("lowgain_sigma", &lowgain_sigma, "lowgain_sigma/D");
        tp.Branch("entries_hg", &entries_hg, "entries_hg/I");
        tp.Branch("entries_lg", &entries_lg, "entries_lg/I");
        tp.Branch("fitStatus_hg", &fitStatus_hg, "fitStatus_hg/I");
        tp.Branch("fitStatus_lg", &fitStatus_lg, "fitStatus_lg/I");
        tp.Branch("fitOk_hg", &fitOk_hg, "fitOk_hg/I");
        tp.Branch("fitOk_lg", &fitOk_lg, "fitOk_lg/I");
        tp.Branch("usable_hg", &usable_hg, "usable_hg/I");
        tp.Branch("usable_lg", &usable_lg, "usable_lg/I");
        tp.Branch("status_hg", &status_hg, "status_hg/I");
        tp.Branch("status_lg", &status_lg, "status_lg/I");
        tp.Branch("x_mm", &x_mm, "x_mm/D");
        tp.Branch("y_mm", &y_mm, "y_mm/D");
        tp.Branch("highgain_chi2", &highgain_chi2, "highgain_chi2/D");
        tp.Branch("lowgain_chi2", &lowgain_chi2, "lowgain_chi2/D");
        tp.Branch("highgain_ndf", &highgain_ndf, "highgain_ndf/I");
        tp.Branch("lowgain_ndf", &lowgain_ndf, "lowgain_ndf/I");
        tp.Branch("highgain_chi2_ndf", &highgain_chi2_ndf, "highgain_chi2_ndf/D");
        tp.Branch("lowgain_chi2_ndf", &lowgain_chi2_ndf, "lowgain_chi2_ndf/D");
        tp.Branch("highgain_tail_3sig_fraction", &highgain_tail_3sig_fraction, "highgain_tail_3sig_fraction/D");
        tp.Branch("lowgain_tail_3sig_fraction", &lowgain_tail_3sig_fraction, "lowgain_tail_3sig_fraction/D");
        tp.Branch("highgain_tail_5sig_fraction", &highgain_tail_5sig_fraction, "highgain_tail_5sig_fraction/D");
        tp.Branch("lowgain_tail_5sig_fraction", &lowgain_tail_5sig_fraction, "lowgain_tail_5sig_fraction/D");
        // tp.Branch("highgain_deviation_squared_sum", &highgain_deviation_squared_sum, "highgain_deviation_squared_sum/D");
        // tp.Branch("lowgain_deviation_squared_sum", &lowgain_deviation_squared_sum, "lowgain_deviation_squared_sum/D");
        // tp.Branch("highgain_deviation_quared_sum", &highgain_deviation_quared_sum, "highgain_deviation_quared_sum/D");
        // tp.Branch("lowgain_deviation_quared_sum", &lowgain_deviation_quared_sum, "lowgain_deviation_quared_sum/D");

        for (const auto& [_, res] : ped_cache_) {
          cellid = res.cellid;
          channel_index = res.channel_index;
          entries_hg = res.entries_hg;
          entries_lg = res.entries_lg;
          
          highgain_rms = res.fit_hg.rms;
          highgain_avg = res.fit_hg.avg;
          highgain_peak = res.fit_hg.mean;
          highgain_sigma = res.fit_hg.sigma;
          fitStatus_hg = res.fit_hg.status;
          fitOk_hg = res.fit_hg.ok ? 1 : 0;
          usable_hg = (fitOk_hg && AHCALRefValues::HGPedestalStatus_is_ok(res.fit_hg.status)) ? 1 : 0;
          status_hg = res.fit_hg.status;
          if (status_hg != AHCALRefValues::HGPedestal_OK) {
            if (status_hg > 1){
              hstatusHG->Fill(-993);
            } else {
              hstatusHG->Fill(status_hg);
            }
            // LOG_WARN("CellID {}: HG pedestal fit failed or cut: fitStatus={}, status={}", cellid, fitStatus_hg, status_hg);
          }
          h2_rms_sigma_HG->Fill(highgain_sigma, highgain_rms);
          highgain_chi2 = res.fit_hg.chi2;
          highgain_ndf = res.fit_hg.ndf;
          highgain_chi2_ndf = (res.fit_hg.ndf > 0) ? res.fit_hg.chi2 / res.fit_hg.ndf : -1.0;
          highgain_tail_3sig_fraction = res.fit_hg.tail_3sig_fraction;
          highgain_tail_5sig_fraction = res.fit_hg.tail_5sig_fraction;
          // highgain_deviation_squared_sum = res.fit_hg.deviation_squared_sum;
          // highgain_deviation_quared_sum = res.fit_hg.deviation_quared_sum;
          lowgain_rms = res.fit_lg.rms;
          lowgain_avg = res.fit_lg.avg;

          lowgain_peak = res.fit_lg.mean;
          lowgain_sigma = res.fit_lg.sigma;
          fitStatus_lg = res.fit_lg.status;
          fitOk_lg = res.fit_lg.ok ? 1 : 0;
          usable_lg = (fitOk_lg && AHCALRefValues::LGPedestalStatus_is_ok(res.fit_lg.status)) ? 1 : 0;
          status_lg = res.fit_lg.status;
          if (status_lg != AHCALRefValues::LGPedestal_OK) {
            // LOG_WARN("CellID {}: LG pedestal fit failed or cut: fitStatus={}, status={}", cellid, fitStatus_lg, status_lg);
            if (status_lg > 1){
              hstatusLG->Fill(-997); // using same histogram for LG fit status for simplicity, but with different underflow bin
            } else {
              hstatusLG->Fill(status_lg);
            }
          }
          h2_rms_sigma_LG->Fill(lowgain_sigma, lowgain_rms);
          lowgain_chi2 = res.fit_lg.chi2;
          lowgain_ndf = res.fit_lg.ndf;
          lowgain_chi2_ndf = (res.fit_lg.ndf > 0) ? res.fit_lg.chi2 / res.fit_lg.ndf : -1.0;
          lowgain_tail_3sig_fraction = res.fit_lg.tail_3sig_fraction;
          lowgain_tail_5sig_fraction = res.fit_lg.tail_5sig_fraction;
          // lowgain_deviation_squared_sum = res.fit_lg.deviation_squared_sum;
          // lowgain_deviation_quared_sum = res.fit_lg.deviation_quared_sum;

          x_mm = res.x_mm;
          y_mm = res.y_mm;

          // fill maps
          const int L = res.layer;
          if (L >= 0 && L < AHCALGeometry::Layer_No) {
            const int bx = hMeanHG[L]->GetXaxis()->FindBin(x_mm);
            const int by = hMeanHG[L]->GetYaxis()->FindBin(y_mm);

            if (bx >= 1 && bx <= hMeanHG[L]->GetNbinsX() && by >= 1 && by <= hMeanHG[L]->GetNbinsY()) {
              if (entries_hg > 0 && usable_hg) {
                hMeanHG[L]->SetBinContent(bx, by, highgain_peak);
                hSigHG[L]->SetBinContent(bx, by, highgain_sigma);
                hEntHG[L]->SetBinContent(bx, by, entries_hg);
              }
              if (entries_lg > 0 && usable_lg) {
                hMeanLG[L]->SetBinContent(bx, by, lowgain_peak);
                hSigLG[L]->SetBinContent(bx, by, lowgain_sigma);
                hEntLG[L]->SetBinContent(bx, by, entries_lg);
              }
            }
          }
          tp.Fill();
        }

        // canvases (7x6) for 40 layers
        auto cAllMeanHG = std::make_unique<TCanvas>("cPedMeanHG_all_7x6", "Pedestal mean maps HG (all layers)", 5600, 4200);
        cAllMeanHG->Divide(7, 6, 0.001, 0.001);

        auto cAllMeanLG = std::make_unique<TCanvas>("cPedMeanLG_all_7x6", "Pedestal mean maps LG (all layers)", 5600, 4200);
        cAllMeanLG->Divide(7, 6, 0.001, 0.001);

        for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
          cAllMeanHG->cd(L + 1);
          gPad->SetMargin(0.08, 0.14, 0.08, 0.10);
          hMeanHG[L]->Draw("COLZ");
          drawLayerLabel(L);

          cAllMeanLG->cd(L + 1);
          gPad->SetMargin(0.08, 0.14, 0.08, 0.10);
          hMeanLG[L]->Draw("COLZ");
          drawLayerLabel(L);
        }

        // write objects
        if (dHG) dHG->cd();
        for (auto& [_, h] : hg_hist_) if (h) h->Write();

        if (dLG) dLG->cd();
        for (auto& [_, h] : lg_hist_) if (h) h->Write();

        if (dMap) dMap->cd();
        for (int L = 0; L < AHCALGeometry::Layer_No; ++L) {
          hMeanHG[L]->Write(); hSigHG[L]->Write(); hEntHG[L]->Write();
          hMeanLG[L]->Write(); hSigLG[L]->Write(); hEntLG[L]->Write();
        }

        if (dCan) dCan->cd();
        cAllMeanHG->Write();
        cAllMeanLG->Write();
        fout->cd();
        auto canvasSummary = std::make_unique<TCanvas>("cPedestalSummary", "Pedestal summary", 1200, 800);
        canvasSummary->Divide(4, 1);
        canvasSummary->cd(1);
        gPad->SetMargin(0.12, 0.08, 0.12, 0.08);
        hstatusHG->Draw("HIST");
        canvasSummary->cd(2);
        gPad->SetMargin(0.12, 0.08, 0.12, 0.08);
        h2_rms_sigma_HG->Draw("COLZ");

        if (cfg_.hg_cut_failed) {
          if (cfg_.hg_cut_rms_per_sigma) {
            TLine line_rms_per_sigma_cut(0.0, 0.0, cfg_.hg_max_sigma, cfg_.hg_max_sigma * cfg_.hg_max_rms_per_sigma);
            line_rms_per_sigma_cut.SetLineColor(kRed);
            line_rms_per_sigma_cut.SetLineWidth(2);
            line_rms_per_sigma_cut.Draw("SAME");
          }
          if (cfg_.hg_cut_max_sigma_and_max_rms) {
            TLine line_max_sigma_and_max_rms_cut(0.0, cfg_.hg_max_rms, cfg_.hg_max_sigma, 0.0);
            line_max_sigma_and_max_rms_cut.SetLineColor(kRed);
            line_max_sigma_and_max_rms_cut.SetLineWidth(2);
            line_max_sigma_and_max_rms_cut.Draw("SAME");
            TLine line_max_sigma_and_max_rms_cut_2(cfg_.hg_max_sigma, 0.0, 0.0, cfg_.hg_max_rms);
            line_max_sigma_and_max_rms_cut_2.SetLineColor(kRed);
            line_max_sigma_and_max_rms_cut_2.SetLineWidth(2);
            line_max_sigma_and_max_rms_cut_2.Draw("SAME");
          }

          if (cfg_.hg_cut_sigma_vs_chi2) {
            TF1 f_cut_sigma_vs_chi2("f_cut_sigma_vs_chi2", Form("%f * (1 - x / %f)", cfg_.hg_sigma_intercept, cfg_.hg_chi2_ndf_intercept), 0.0, cfg_.hg_chi2_ndf_intercept);
            f_cut_sigma_vs_chi2.SetLineColor(kBlue);
            f_cut_sigma_vs_chi2.SetLineWidth(2);
            f_cut_sigma_vs_chi2.Draw("SAME");
          }
        }
        canvasSummary->cd(3);
        gPad->SetMargin(0.12, 0.08, 0.12, 0.08);
        hstatusLG->Draw("HIST");
        canvasSummary->cd(4);
        gPad->SetMargin(0.12, 0.08, 0.12, 0.08);
        h2_rms_sigma_LG->Draw("COLZ");
        canvasSummary->Write();
        h2_rms_sigma_LG->Write();
        h2_rms_sigma_HG->Write();
        hstatusLG->Write();
        hstatusHG->Write();
        // Skip writing hstatusHG to avoid ROOT file corruption issues with mixed status value types
        tp.Write();
        fout->Close();

        LOG_INFO("PedestalAlg: wrote {}", cfg_.out_pedestal_filename);
      }
    }

    LOG_INFO("PedestalAlg: HG OK/all = {}/{}", nFitOK_HG_, nFitAll_HG_);
    std::cout << "PedestalAlg: HG OK/all = " << nFitOK_HG_ << "/" << nFitAll_HG_ << std::endl;
    LOG_INFO("PedestalAlg: LG OK/all = {}/{}", nFitOK_LG_, nFitAll_LG_);
    std::cout << "PedestalAlg: LG OK/all = " << nFitOK_LG_ << "/" << nFitAll_LG_ << std::endl;

    if (cfg_.pedestal_to_json) writeJson(/*cache_ready=*/true);
  }

  std::unique_ptr<TH1D> make_hist(int cellid, bool isHG) {
    int layer = cellid/100000;
    int chip = cellid/10000 % 10;
    const std::string name = (isHG ? "hPedHG_" : "hPedLG_") + std::to_string(cellid);
    const std::string title = "Layer " + std::to_string(layer) + " chip "+ std::to_string(chip) + (isHG ? " Pedestal HG;ADC;counts" : " Pedestal LG;ADC;counts");
    auto h = std::make_unique<TH1D>(name.c_str(), title.c_str(), cfg_.nbin, cfg_.xmin, cfg_.xmax);
    h->SetDirectory(nullptr);
    return h;
  }

  void writeJson(bool cache_ready = false) {
    if (!cfg_.pedestal_to_json) return;

    if (!cache_ready) buildPedestalCache();

    // Create output directory if it doesn't exist
    std::filesystem::path out_dir(cfg_.out_json_dirname);
    try {
      std::filesystem::create_directories(out_dir);
    } catch (const std::exception& e) {
      LOG_ERROR("PedestalAlg: cannot create output directory: {} - {}", cfg_.out_json_dirname, e.what());
      return;
    }

    // Generate ISO8601 timestamp
    double start_time = ctx_.conditions.starttime;
    double end_time = ctx_.conditions.endtime;
    double ref_time = (start_time > 0.0 && end_time > 0.0) ? (start_time + end_time) / 2.0 : std::time(nullptr);
    time_t utc_time = static_cast<time_t>(ref_time);
    std::tm tm = {};
    gmtime_r(&utc_time, &tm);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    std::string timestamp = oss.str();

    const int n_channels_per_layer = AHCALGeometry::chip_No * AHCALGeometry::channel_No;

    // Loop over each layer and write separate JSON file
    for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer) {
      // Prepare per-layer arrays (indexed by channel_index = chip*36 + channel)
      std::vector<int> cellid_arr(n_channels_per_layer, -1);
      std::vector<double> hg_peak_arr(n_channels_per_layer, -1.0);
      std::vector<double> hg_sigma_arr(n_channels_per_layer, -1.0);
      std::vector<int>    hg_status_arr(n_channels_per_layer, -999);
      std::vector<int>    hg_usable_arr(n_channels_per_layer, 0);
      std::vector<double> lg_peak_arr(n_channels_per_layer, -1.0);
      std::vector<double> lg_sigma_arr(n_channels_per_layer, -1.0);
      std::vector<int>    lg_status_arr(n_channels_per_layer, -999);
      std::vector<int>    lg_usable_arr(n_channels_per_layer, 0);


      // Fill arrays for this layer
      int fit_failures_hg = 0;
      int fit_failures_lg = 0;
      int entries_sum = 0;
      for (const auto& [_, res] : ped_cache_) {
        if (res.layer != layer) continue;

        const int layer_idx = res.chip * AHCALGeometry::channel_No + res.channel;
        if (layer_idx < 0 || layer_idx >= n_channels_per_layer) continue;

        cellid_arr[layer_idx] = res.cellid;
        hg_peak_arr[layer_idx] = res.fit_hg.mean;
        hg_sigma_arr[layer_idx] = res.fit_hg.sigma;
        hg_status_arr[layer_idx] = res.fit_hg.status;
        hg_usable_arr[layer_idx] = (res.fit_hg.ok && AHCALRefValues::HGPedestalStatus_is_ok(res.fit_hg.status)) ? 1 : 0;
        lg_peak_arr[layer_idx] = res.fit_lg.mean;
        lg_sigma_arr[layer_idx] = res.fit_lg.sigma;
        lg_status_arr[layer_idx] = res.fit_lg.status;
        lg_usable_arr[layer_idx] = (res.fit_lg.ok && AHCALRefValues::LGPedestalStatus_is_ok(res.fit_lg.status)) ? 1 : 0;
        if (!res.fit_hg.ok) fit_failures_hg++;
        if (!res.fit_lg.ok) fit_failures_lg++;
        entries_sum = res.entries_hg;
      }

      // Fill missing cellids
      int dead_channels = 0;
      for (int idx = 0; idx < n_channels_per_layer; ++idx) {
        if (cellid_arr[idx] < 0) {
          dead_channels++;
          int chip = idx / AHCALGeometry::channel_No;
          int channel = idx % AHCALGeometry::channel_No;
          cellid_arr[idx] = AHCALGeometry::CellID(layer, chip, channel);
        }
      }

      // Build JSON object for this layer
      json j;

      j["RunNumber"] = ctx_.config.runNumber;
      j["TimeStamp"] = timestamp;
      j["Layer"] = layer;
      j["CalibrationType"] = "Pedestal";
      j["Summary"]["DeadChannels"] = dead_channels;
      j["Summary"]["HighGainFitFailures"] = fit_failures_hg;
      j["Summary"]["LowGainFitFailures"] = fit_failures_lg;
      j["Summary"]["Entries"] = entries_sum;
      j["Status"] = (fit_failures_hg == 0 && fit_failures_lg == 0) ? 0 : 1;

      // Per-channel arrays (following Pedestal.schema)
      // j["PerChannel"]["CellID"] = cellid_arr;
      j["PerChannel"]["HighGainPeak"] = hg_peak_arr;
      j["PerChannel"]["HighGainSigma"] = hg_sigma_arr;
      j["PerChannel"]["HighGainStatus"] = hg_status_arr;
      j["PerChannel"]["LowGainPeak"] = lg_peak_arr;
      j["PerChannel"]["LowGainSigma"] = lg_sigma_arr;
      j["PerChannel"]["LowGainStatus"] = lg_status_arr;
      j["PerChannel"]["HighGainUsable"] = hg_usable_arr;
      j["PerChannel"]["LowGainUsable"] = lg_usable_arr;

      // Write to file
      std::ostringstream filename;
      filename << cfg_.out_json_dirname << "/pedestal_Layer" << layer << ".json";
      std::string out_filename = filename.str();

      std::ofstream jout(out_filename);
      if (!jout.is_open()) {
        LOG_ERROR("PedestalAlg: cannot create JSON output file: {}", out_filename);
        continue;
      }
      jout << j.dump(2) << std::endl;
      jout.close();

      LOG_INFO("PedestalAlg: wrote JSON {}", out_filename);
    }
  }

  PedestalAlgCfg cfg_;
  const RunContext& ctx_;
  bool written_ = false;
  bool cache_built_ = false;
  int nFitOK_HG_ = 0;
  int nFitAll_HG_ = 0;
  int nFitOK_LG_ = 0;
  int nFitAll_LG_ = 0;
  std::unordered_map<int, std::unique_ptr<TH1D>> hg_hist_;
  std::unordered_map<int, std::unique_ptr<TH1D>> lg_hist_;
  std::unordered_map<int, CellPedestalResult> ped_cache_;
};

// Define deleter *after* Impl is a complete type in this TU.
void PedestalAlg::ImplDeleter::operator()(PedestalAlg::Impl* p) const {
  delete p;
}

PedestalAlg::~PedestalAlg() {
  if (impl_) impl_->write();
}

void PedestalAlg::execute(EventStore& evt) {
  if (!impl_) impl_.reset(new Impl(cfg_, ctx()));

  auto raw_hits = evt.get<std::vector<AHCALRawHit>>(cfg_.in_rawhit_key);
  for (const auto& h : raw_hits) {
    impl_->fill(h);
  }
}

void PedestalAlg::parse_cfg(const YAML::Node& n) {
  cfg_.in_rawhit_key = get_or<std::string>(n, "in_rawhit_key", cfg_.in_rawhit_key);

  cfg_.pedestal_to_file = get_or<bool>(n, "pedestal_to_file", cfg_.pedestal_to_file);
  cfg_.pedestal_to_json = get_or<bool>(n, "pedestal_to_json", cfg_.pedestal_to_json);
  cfg_.out_pedestal_filename = get_or<std::string>(n, "out_pedestal_filename", cfg_.out_pedestal_filename);
  cfg_.out_json_dirname = get_or<std::string>(n, "out_json_dirname", cfg_.out_json_dirname);

  cfg_.nbin = get_or<int>(n, "nbin", cfg_.nbin);
  cfg_.xmin = get_or<double>(n, "xmin", cfg_.xmin);
  cfg_.xmax = get_or<double>(n, "xmax", cfg_.xmax);

  cfg_.min_entries = get_or<int>(n, "min_entries", cfg_.min_entries);
  cfg_.nsigma_win1 = get_or<double>(n, "nsigma_win1", cfg_.nsigma_win1);
  cfg_.nsigma_win2 = get_or<double>(n, "nsigma_win2", cfg_.nsigma_win2);
  cfg_.sigma_min = get_or<double>(n, "sigma_min", cfg_.sigma_min);
  cfg_.sigma_max = get_or<double>(n, "sigma_max", cfg_.sigma_max);

  cfg_.use_hittag = get_or<bool>(n, "use_hittag", cfg_.use_hittag);
  cfg_.select_hittag = get_or<int>(n, "select_hittag", cfg_.select_hittag);

  cfg_.hg_cut_failed = get_or<bool>(n, "hg_cut_failed", cfg_.hg_cut_failed);
  cfg_.hg_cut_rms_per_sigma = get_or<bool>(n, "hg_cut_rms_per_sigma", cfg_.hg_cut_rms_per_sigma);
  cfg_.hg_max_rms_per_sigma = get_or<double>(n, "hg_max_rms_per_sigma", cfg_.hg_max_rms_per_sigma);
  cfg_.hg_cut_max_sigma_and_max_rms = get_or<bool>(n, "hg_cut_max_sigma_and_max_rms", cfg_.hg_cut_max_sigma_and_max_rms);
  cfg_.hg_max_sigma = get_or<double>(n, "hg_max_sigma", cfg_.hg_max_sigma);
  cfg_.hg_max_rms = get_or<double>(n, "hg_max_rms", cfg_.hg_max_rms);
  cfg_.hg_cut_sigma_vs_chi2 = get_or<bool>(n, "hg_cut_sigma_vs_chi2", cfg_.hg_cut_sigma_vs_chi2);
  cfg_.hg_sigma_intercept = get_or<double>(n, "hg_sigma_intercept", cfg_.hg_sigma_intercept);
  cfg_.hg_chi2_ndf_intercept = get_or<double>(n, "hg_chi2_ndf_intercept", cfg_.hg_chi2_ndf_intercept);

  
  cfg_.lg_cut_failed = get_or<bool>(n, "lg_cut_failed", cfg_.lg_cut_failed);
  cfg_.lg_cut_max_sigma = get_or<bool>(n, "lg_cut_max_sigma", cfg_.lg_cut_max_sigma);
  cfg_.lg_max_sigma = get_or<double>(n, "lg_max_sigma", cfg_.lg_max_sigma);
}

} // namespace AHCALRecoAlg