#include "PedestalAlg.hpp"

#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"

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
  int    status = 999;
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
    r.status = st2;
    r.ok = false;
    return r;
  }

  r.mean = f2.GetParameter(1);
  r.sigma = std::fabs(f2.GetParameter(2));
  r.status = st2;
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
  explicit Impl(PedestalAlgCfg cfg) : cfg_(std::move(cfg)) {}

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

  void write() {
    if (!cfg_.pedestal_to_file) return;
    if (written_) return;
    written_ = true;

    auto fout = std::unique_ptr<TFile>(TFile::Open(cfg_.out_pedestal_filename.c_str(), "RECREATE"));
    if (!fout || fout->IsZombie()) {
      LOG_ERROR("PedestalAlg: cannot create output file: {}", cfg_.out_pedestal_filename);
      return;
    }

    // directories
    TDirectory* dHist = ensureDir(fout.get(), "Pedestal");
    TDirectory* dHG   = dHist ? ensureDir(dHist, "HG") : nullptr;
    TDirectory* dLG   = dHist ? ensureDir(dHist, "LG") : nullptr;

    TDirectory* dMap  = ensureDir(fout.get(), "PedMap2D");
    TDirectory* dCan  = ensureDir(fout.get(), "Canvases");

    // union of keys (HG/LG)
    std::unordered_set<int> keys;
    keys.reserve(hg_hist_.size() + lg_hist_.size());
    for (const auto& [k, _] : hg_hist_) keys.insert(k);
    for (const auto& [k, _] : lg_hist_) keys.insert(k);

    // per-layer maps (HG/LG): fixed Layer_No from geometry
    std::vector<std::unique_ptr<TH2D>> hMeanHG(AHCALGeometry::Layer_No), hSigHG(AHCALGeometry::Layer_No), hEntHG(AHCALGeometry::Layer_No);
    std::vector<std::unique_ptr<TH2D>> hMeanLG(AHCALGeometry::Layer_No), hSigLG(AHCALGeometry::Layer_No), hEntLG(AHCALGeometry::Layer_No);

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
    int entries_hg=0, entries_lg=0;
    int fitStatus_hg=999, fitStatus_lg=999;
    int fitOk_hg=0, fitOk_lg=0;

    // For AdcToEnergyReadTTreeAlg::initialize_ped()
    double highgain_peak=-1.0, lowgain_peak=-1.0;

    // Extra info
    double highgain_sigma=-1.0, lowgain_sigma=-1.0;
    double x_mm=-999.0, y_mm=-999.0;

    tp.Branch("cellid", &cellid, "cellid/I");
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
    tp.Branch("x_mm", &x_mm, "x_mm/D");
    tp.Branch("y_mm", &y_mm, "y_mm/D");

    int nFitOK_HG = 0, nFitAll_HG = 0;
    int nFitOK_LG = 0, nFitAll_LG = 0;

    for (int cid : keys) {
      cellid = cid;

      // decode using EDM helpers (cellID format is defined in AHCALRawHit)
      AHCALRawHit tmp;
      tmp.cellID = cellid;
      const int L  = tmp.layer();
      const int C  = tmp.chip();
      const int ch = tmp.channel();

      cellid_to_xy(C, ch, x_mm, y_mm);

      TH1D* hhg = nullptr;
      TH1D* hlg = nullptr;
      if (auto it = hg_hist_.find(cellid); it != hg_hist_.end()) hhg = it->second.get();
      if (auto it = lg_hist_.find(cellid); it != lg_hist_.end()) hlg = it->second.get();

      FitOut fr_hg;
      entries_hg = hhg ? static_cast<int>(hhg->GetEntries()) : 0;
      if (hhg) {
        nFitAll_HG++;
        fr_hg = fitPedestalGaussian(hhg, cfg_.min_entries, cfg_.nsigma_win1, cfg_.nsigma_win2,
                                    cfg_.sigma_min, cfg_.sigma_max);
        if (fr_hg.ok) nFitOK_HG++;
      }
      highgain_peak  = fr_hg.mean;
      highgain_sigma = fr_hg.sigma;
      fitStatus_hg   = fr_hg.status;
      fitOk_hg       = fr_hg.ok ? 1 : 0;

      FitOut fr_lg;
      entries_lg = hlg ? static_cast<int>(hlg->GetEntries()) : 0;
      if (hlg) {
        nFitAll_LG++;
        fr_lg = fitPedestalGaussian(hlg, cfg_.min_entries, cfg_.nsigma_win1, cfg_.nsigma_win2,
                                    cfg_.sigma_min, cfg_.sigma_max);
        if (fr_lg.ok) nFitOK_LG++;
      }
      lowgain_peak  = fr_lg.mean;
      lowgain_sigma = fr_lg.sigma;
      fitStatus_lg  = fr_lg.status;
      fitOk_lg      = fr_lg.ok ? 1 : 0;

      // fill maps
      if (L >= 0 && L < AHCALGeometry::Layer_No) {
        const int bx = hMeanHG[L]->GetXaxis()->FindBin(x_mm);
        const int by = hMeanHG[L]->GetYaxis()->FindBin(y_mm);

        if (bx >= 1 && bx <= hMeanHG[L]->GetNbinsX() && by >= 1 && by <= hMeanHG[L]->GetNbinsY()) {
          if (entries_hg > 0) {
            hMeanHG[L]->SetBinContent(bx, by, highgain_peak);
            hSigHG[L]->SetBinContent(bx, by, highgain_sigma);
            hEntHG[L]->SetBinContent(bx, by, entries_hg);
          }
          if (entries_lg > 0) {
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
    tp.Write();
    fout->Close();

    LOG_INFO("PedestalAlg: wrote {}", cfg_.out_pedestal_filename);
    LOG_INFO("PedestalAlg: HG fit OK/all = {}/{}", nFitOK_HG, nFitAll_HG);
    LOG_INFO("PedestalAlg: LG fit OK/all = {}/{}", nFitOK_LG, nFitAll_LG);
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

  PedestalAlgCfg cfg_;
  bool written_ = false;
  std::unordered_map<int, std::unique_ptr<TH1D>> hg_hist_;
  std::unordered_map<int, std::unique_ptr<TH1D>> lg_hist_;
};

// Define deleter *after* Impl is a complete type in this TU.
void PedestalAlg::ImplDeleter::operator()(PedestalAlg::Impl* p) const {
  delete p;
}

PedestalAlg::~PedestalAlg() {
  if (impl_) impl_->write();
}

void PedestalAlg::execute(EventStore& evt) {
  if (!impl_) impl_.reset(new Impl(cfg_));

  auto raw_hits = evt.get<std::vector<AHCALRawHit>>(cfg_.in_rawhit_key);
  for (const auto& h : raw_hits) {
    impl_->fill(h);
  }
}

void PedestalAlg::parse_cfg(const YAML::Node& n) {
  cfg_.in_rawhit_key = get_or<std::string>(n, "in_rawhit_key", cfg_.in_rawhit_key);

  cfg_.pedestal_to_file = get_or<bool>(n, "pedestal_to_file", cfg_.pedestal_to_file);
  cfg_.out_pedestal_filename = get_or<std::string>(n, "out_pedestal_filename", cfg_.out_pedestal_filename);

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
}

} // namespace AHCALRecoAlg