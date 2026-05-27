#include "InterCalibAlg.hpp"

#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"
#include "CalibDBIO/PedestalReader/PedestalReader.hpp"
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
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
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <nlohmann/json.hpp>

using json = nlohmann::json;
AHCAL_REGISTER_ALG(AHCALRecoAlg::InterCalibAlg, "InterCalibAlg")

namespace AHCALRecoAlg {

// Quality flag bitmask for HG:LG calibration
enum HGLGQualityFlag : uint32_t {
    GOOD                = 0,
    LOW_STAT            = 1 << 0,
    FIT_FAILED          = 1 << 1,
    BAD_RMSPERSIGMA     = 1 << 2,
    PEDESTAL_MASKED     = 1 << 6,
    MANUAL_BAD          = 1 << 7,
    MANUAL_GOOD         = 1 << 8
};

// Per-channel HG:LG calibration result with quality metrics
struct HGLGCalibResult {
    int run = -1;
    int layer = -1;
    int chip = -1;
    int channel = -1;
    int cellid = -1;

    double hg_pedestal = NAN;
    double lg_pedestal = NAN;

    double slope = NAN;
    double intercept = NAN;
    double slope_error = NAN;
    double intercept_error = NAN;

    double chi2 = NAN;
    int ndf = -1;
    double chi2_ndf = NAN;
    int fit_status = -1;

    int n_points_total = 0;
    int n_points_used = 0;
    int n_points_outlier = 0;
    double outlier_fraction = -1.0;

    int n_fit_points = 0; // number of points used in the final fit after outlier removal
    int n_fit_points_over500 = 0; // number of fit points with HG>500

    double correlation = -999.0;

    double distance_avg = -1.0;
    double distance_rms = -1.0;
    double distance_sigma = -1.0;
    std::vector<double> distances;  // Distance of each raw point to the fitted line in HG-LG space
    
    int n_points_HGLG_cluster = 0; // lg>50 & HG<1000 & HG<2*LG
    
    double slope_median_layer = NAN;
    double slope_mad_layer = NAN;
    double slope_z_layer = -999.0;

    uint32_t quality_flag = 0;

    bool is_pedestal_masked = false;
    bool is_hglg_calib_bad = false;
    bool is_physics_mask_candidate = false;

    bool manual_bad = false;
    bool manual_good = false;
};

namespace {

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
  double asymmetry = -1.0;
};

static inline double peakFromMaxBin(TH1D* h) {
  const int b = h->GetMaximumBin();
  return h->GetBinCenter(b);
}

static FitOut fitSimpleGaussian(TH1D* h, int minEntries, double sigmaMin, double sigmaMax) {
    FitOut r;
    if (!h) return r;
    if (h->GetEntries() < minEntries) return r;

    double mu0 = peakFromMaxBin(h);
    double rms = h->GetRMS();
    r.rms = rms;
    r.avg = h->GetMean();
    if (!(rms > 0.0)) rms = 10.0;

    double sig0 = std::clamp(rms, sigmaMin, sigmaMax);

    const std::string fname = std::string("f_") + h->GetName();
    TF1 f(fname.c_str(), "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    f.SetParameters(h->GetMaximum(), mu0, sig0);
    f.SetParLimits(2, sigmaMin, sigmaMax);

    const int st = h->Fit(&f, "RQ0");
    if (st != 0) {
        r.mean = mu0;
        r.sigma = sig0;
        r.fitStatus = st;
        r.ok = false;
        r.chi2 = f.GetChisquare();
        r.ndf = f.GetNDF();
        return r;
    }
    r.mean = f.GetParameter(1);
    r.sigma = std::fabs(f.GetParameter(2));
    r.fitStatus = st;
    r.ok = true;
    r.chi2 = f.GetChisquare();
    r.ndf = f.GetNDF();
    double sig3_upper = f.GetParameter(1) + 3.0 * std::fabs(f.GetParameter(2));
    double sig5_upper = f.GetParameter(1) + 5.0 * std::fabs(f.GetParameter(2));
    double sig3_lower = f.GetParameter(1) - 3.0 * std::fabs(f.GetParameter(2));
    double sig5_lower = f.GetParameter(1) - 5.0 * std::fabs(f.GetParameter(2));
    int entries = static_cast<int>(h->GetEntries());
    int tail_3sig = h->Integral(h->FindBin(sig3_upper), h->GetNbinsX()) + h->Integral(1, h->FindBin(sig3_lower));
    int tail_5sig = h->Integral(h->FindBin(sig5_upper), h->GetNbinsX()) + h->Integral(1, h->FindBin(sig5_lower));
    r.tail_3sig_fraction = (entries > 0) ? static_cast<double>(tail_3sig) / entries : -1.0;
    r.tail_5sig_fraction = (entries > 0) ? static_cast<double>(tail_5sig) / entries : -1.0;
    r.asymmetry = (h->Integral(h->FindBin(f.GetParameter(1)), h->GetNbinsX()) - h->Integral(1, h->FindBin(f.GetParameter(1)))) / (entries > 0 ? static_cast<double>(entries) : 1.0);
    r.deviation_squared_sum =  0.0;
    r.deviation_squared_sum_fourth = 0.0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        double x = h->GetBinCenter(b);
        double dev = (f.GetParameter(2) > 0.0 ? (x - f.GetParameter(1)) : -1.0);
        dev *= h->GetBinContent(b);
        if (dev > 0) {
            r.deviation_squared_sum += dev * dev;
            r.deviation_squared_sum_fourth += dev * dev * dev * dev;
        }
    }
    r.deviation_squared_sum = (entries > 0) ? r.deviation_squared_sum / entries : -1.0;
    r.deviation_squared_sum_fourth = (entries > 0) ? r.deviation_squared_sum_fourth / entries : -1.0;
    return r;
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
  r.asymmetry = (h->Integral(h->FindBin(f2.GetParameter(1)), h->GetNbinsX()) - h->Integral(1, h->FindBin(f2.GetParameter(1)))) / (entries > 0 ? static_cast<double>(entries) : 1.0);
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

// Linear fit result (HG = p0 + p1*LG)
struct FitResult {
    double p0 = -999.0;
    double p1 = -999.0;
    double chi2 = -1.0;
    int ndf = -1;
    double chi2_ndf = -1.0;
    int status = 999;
    bool ok = false;
    int n_fit_points = 0;
    int n_fit_points_over500 = 0;
};

struct ProfilePoint {
    double hg = 0.0;
    double lg = 0.0;
    double lg_err = 0.0;
};

// Extract ProfilePoints from TH2D histogram
static std::vector<ProfilePoint> buildHGProfilePointsFromHist(const TH2D* hist, double hgBinWidth, int min_hg_bins_for_fit = 8) {
    std::vector<ProfilePoint> pts;
    if (!hist || hgBinWidth <= 0.0) return pts;
    
    struct BinStats { 
        int n = 0; 
        double sum_lg = 0.0; 
        double sum_lg2 = 0.0; 
        double sum_hg = 0.0; 
    };
    std::unordered_map<int, BinStats> bins;
    
    // Iterate through all histogram bins
    for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
        double hg_center = hist->GetYaxis()->GetBinCenter(iy);
        for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
            double content = hist->GetBinContent(ix, iy);
            if (content <= 0) continue;
            
            double lg_center = hist->GetXaxis()->GetBinCenter(ix);
            const int ibin = static_cast<int>(std::floor(hg_center / hgBinWidth));
            auto& b = bins[ibin];
            
            // Each bin count acts as weight
            for (int k = 0; k < static_cast<int>(content); ++k) {
                b.n++;
                b.sum_hg += hg_center;
                b.sum_lg += lg_center;
                b.sum_lg2 += lg_center * lg_center;
            }
        }
    }
    
    pts.reserve(bins.size());
    for (const auto& [_, b] : bins) {
        if (b.n < min_hg_bins_for_fit) continue;
        const double mean_hg = b.sum_hg / b.n;
        const double mean_lg = b.sum_lg / b.n;
        const double var_lg = std::max(0.0, b.sum_lg2 / b.n - mean_lg * mean_lg);
        pts.push_back(ProfilePoint{mean_hg, mean_lg, std::max(std::sqrt(var_lg), 1.0)});
    }
    std::sort(pts.begin(), pts.end(), [](const ProfilePoint& l, const ProfilePoint& r) { return l.hg < r.hg; });
    return pts;
}

static std::vector<ProfilePoint> buildHGProfilePoints(const std::vector<std::pair<double, double>>& points_hg_lg, double hgBinWidth, int min_hg_bins_for_fit = 8) {
    std::vector<ProfilePoint> pts;
    if (points_hg_lg.empty() || hgBinWidth <= 0.0) return pts;
    struct BinStats { int n = 0; double sum_lg = 0.0; double sum_lg2 = 0.0; double sum_hg = 0.0; };
    std::unordered_map<int, BinStats> bins;
    bins.reserve(points_hg_lg.size() / 8 + 1);
    for (const auto& p : points_hg_lg) {
        const int ibin = static_cast<int>(std::floor(p.first / hgBinWidth));
        auto& b = bins[ibin];
        b.n++;
        b.sum_hg += p.first;
        b.sum_lg += p.second;
        b.sum_lg2 += p.second * p.second;
    }
    pts.reserve(bins.size());
    for (const auto& [_, b] : bins) {
        if (b.n < min_hg_bins_for_fit) continue;
        const double mean_hg = b.sum_hg / b.n;
        const double mean_lg = b.sum_lg / b.n;
        const double var_lg = std::max(0.0, b.sum_lg2 / b.n - mean_lg * mean_lg);
        pts.push_back(ProfilePoint{mean_hg, mean_lg, std::max(std::sqrt(var_lg), 1.0)});
    }
    std::sort(pts.begin(), pts.end(), [](const ProfilePoint& l, const ProfilePoint& r) { return l.hg < r.hg; });
    return pts;
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

static FitResult fitLinearFromHGBinnedProfile(const std::vector<std::pair<double, double>>& points_hg_lg, long long n, int minPoints, double hgBinWidth, int minBins, int min_hg_bins_for_fit = 8) {
    FitResult r;
    r.n_fit_points = 0;
    r.n_fit_points_over500 = 0;
    for (const auto& p : points_hg_lg) {
        r.n_fit_points++;
        if (p.first > 500.0) r.n_fit_points_over500++;
    }
    if (n < minPoints || points_hg_lg.empty() || hgBinWidth <= 0.0) return r;

    const auto profile_pts = buildHGProfilePoints(points_hg_lg, hgBinWidth, min_hg_bins_for_fit);
    if (static_cast<int>(profile_pts.size()) < minBins) return r;

    // Weighted fit of LG = c*HG + d
    double S = 0.0, SX = 0.0, SY = 0.0, SXX = 0.0, SXY = 0.0;
    int used_bins = 0;
    for (const auto& p : profile_pts) {
        const double mean_hg = p.hg;
        const double mean_lg = p.lg;
        const double sigma_eff = p.lg_err;
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
    for (const auto& p : profile_pts) {
        const double mean_hg = p.hg;
        const double mean_lg = p.lg;
        const double sigma_lg = p.lg_err;
        const double resid = mean_lg - (c * mean_hg + d);
        r.chi2 += (resid * resid) / (sigma_lg * sigma_lg);
    }

    r.chi2_ndf = r.chi2 / r.ndf;
    r.status = 0;
    r.ok = true;
    return r;
}

static FitResult fitLinearFromHGBinnedProfile(const TH2D* hist, int minPoints, double hgBinWidth, int minBins, int min_hg_bins_for_fit = 8) {
    FitResult r;
    if (!hist || hgBinWidth <= 0.0) return r;
    
    long long total_entries = static_cast<long long>(hist->GetEntries());
    if (total_entries < minPoints) return r;
    
    const auto profile_pts = buildHGProfilePointsFromHist(hist, hgBinWidth, min_hg_bins_for_fit);
    if (static_cast<int>(profile_pts.size()) < minBins) return r;

    // Same weighted fit logic as before
    double S = 0.0, SX = 0.0, SY = 0.0, SXX = 0.0, SXY = 0.0;
    int used_bins = 0;
    for (const auto& p : profile_pts) {
        const double mean_hg = p.hg;
        const double mean_lg = p.lg;
        const double sigma_eff = p.lg_err;
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
    r.p1 = 1.0 / c;
    r.p0 = -d / c;
    r.ndf = used_bins - 2;
    if (r.ndf <= 0) return r;
    r.chi2 = 0.0;
    for (const auto& p : profile_pts) {
        const double mean_hg = p.hg;
        const double mean_lg = p.lg;
        const double sigma_lg = p.lg_err;
        const double resid = mean_lg - (c * mean_hg + d);
        r.chi2 += (resid * resid) / (sigma_lg * sigma_lg);
    }
    r.chi2_ndf = r.chi2 / r.ndf;
    r.status = 0;
    r.ok = true;
    return r;
}

static int countAndMarkOutliers(
    const std::vector<std::pair<double, double>>& raw_pts,
    double c, double d, double sigmaThr,
    const std::unordered_map<int, double>& sigma_by_bin, double hgBinWidth,
    std::vector<bool>& is_outlier) {
    is_outlier.assign(raw_pts.size(), false);
    int n_outlier = 0;
    if (sigmaThr <= 0.0) return 0;
    for (size_t i = 0; i < raw_pts.size(); ++i) {
        const double hg = raw_pts[i].first;
        const double lg = raw_pts[i].second;
        const int ibin = static_cast<int>(std::floor(hg / hgBinWidth));
        auto its = sigma_by_bin.find(ibin);
        const double sigma = (its != sigma_by_bin.end()) ? std::max(its->second, 1.0) : 5.0;
        const double pull = std::abs(lg - (c * hg + d)) / sigma;
        if (pull > sigmaThr) {
            is_outlier[i] = true;
            n_outlier++;
        }
    }
    return n_outlier;
}

static FitResult fitLinearFromHGBinnedProfileRejectOutliers(
    const std::vector<std::pair<double, double>>& points_hg_lg, long long n, int minPoints, double hgBinWidth, int minBins, double sigmaThr, int& n_outlier_points, int min_hg_bins_for_fit = 8) {
    n_outlier_points = 0;
    FitResult r;
    if (n < minPoints || points_hg_lg.empty() || hgBinWidth <= 0.0) return r;
    const auto profile_pts = buildHGProfilePoints(points_hg_lg, hgBinWidth, min_hg_bins_for_fit);
    if (static_cast<int>(profile_pts.size()) < minBins) return r;

    // First pass
    auto first = fitLinearFromHGBinnedProfile(points_hg_lg, n, minPoints, hgBinWidth, minBins, min_hg_bins_for_fit);
    if (!first.ok) return first;
    const double c1 = 1.0 / first.p1;
    const double d1 = -first.p0 / first.p1;
    std::unordered_map<int, double> sigma_by_bin;
    sigma_by_bin.reserve(profile_pts.size());
    for (const auto& p : profile_pts) {
        sigma_by_bin[static_cast<int>(std::floor(p.hg / hgBinWidth))] = p.lg_err;
    }
    std::vector<bool> is_outlier;
    n_outlier_points = countAndMarkOutliers(points_hg_lg, c1, d1, sigmaThr, sigma_by_bin, hgBinWidth, is_outlier);
    if (n_outlier_points == 0) return first;

    std::vector<std::pair<double, double>> inlier_points;
    inlier_points.reserve(points_hg_lg.size());
    for (size_t i = 0; i < points_hg_lg.size(); ++i) {
        if (is_outlier[i]) continue;
        inlier_points.push_back(points_hg_lg[i]);
    }
    if (static_cast<long long>(inlier_points.size()) < minPoints) return first;
    r.n_fit_points = static_cast<int>(inlier_points.size());
    for (const auto& p : inlier_points) {
        if (p.first > 500.0) r.n_fit_points_over500++;
    }
    return fitLinearFromHGBinnedProfile(inlier_points, static_cast<long long>(inlier_points.size()), minPoints, hgBinWidth, minBins, min_hg_bins_for_fit);
}

// Compute Pearson correlation coefficient between two vectors
static double computeCorrelation(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) return -999.0;
    
    double mean_x = 0.0, mean_y = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        mean_x += x[i];
        mean_y += y[i];
    }
    mean_x /= x.size();
    mean_y /= y.size();
    
    double cov = 0.0, var_x = 0.0, var_y = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        double dx = x[i] - mean_x;
        double dy = y[i] - mean_y;
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }
    
    double sigma_x = std::sqrt(var_x / x.size());
    double sigma_y = std::sqrt(var_y / y.size());
    
    if (sigma_x == 0.0 || sigma_y == 0.0) return -999.0;
    
    return cov / (x.size() * sigma_x * sigma_y);
}

// Compute median absolute deviation (MAD)
static std::pair<double, double> computeMedianAndMAD(std::vector<double> values) {
    if (values.empty()) return {NAN, NAN};
    
    std::sort(values.begin(), values.end());
    double median = values.size() % 2 == 0 
        ? (values[values.size() / 2 - 1] + values[values.size() / 2]) / 2.0
        : values[values.size() / 2];
    
    std::vector<double> abs_dev;
    abs_dev.reserve(values.size());
    for (double v : values) {
        abs_dev.push_back(std::abs(v - median));
    }
    std::sort(abs_dev.begin(), abs_dev.end());
    double mad = abs_dev.size() % 2 == 0
        ? (abs_dev[abs_dev.size() / 2 - 1] + abs_dev[abs_dev.size() / 2]) / 2.0
        : abs_dev[abs_dev.size() / 2];
    
    return {median, mad};
}

} // namespace

struct InterCalibAlg::Impl {
    explicit Impl(InterCalibAlgCfg cfg, const RunContext& ctx)
        : cfg_(std::move(cfg)), ctx_(ctx) {
        loadPedestals();
        run_contexts_.push_back(ctx);
        hg_saturation_map_.clear();
        lg_saturation_map_.clear();
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
        int n_outlier_points = 0;
        int n_hg_lt_2lg_aftercut = 0;
        int n_noise_points = 0;

        int n_fit_points = 0;
        int n_fit_points_over500 = 0;
        double hg_adc_saturation = 0.0;
        double lg_adc_saturation = 0.0;

        double distance_avg = -1.0;
        double distance_rms = -1.0;
        double distance_sigma = -1.0;
        std::vector<double> distances;  // Distance of each raw point to the fitted line in HG-LG space
        std::vector<double> residuals_lg;  // LG residuals for each raw point (measured - fitted)
        std::vector<double> residuals_hg;  // HG residuals for each raw point (measured - fitted)
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
    std::unique_ptr<TH2D> createAccumHistogram(int layer, int chip, int channel) {
        auto hist = std::make_unique<TH2D>(
            Form("h_accum_L%02d_C%02d_ch%02d", layer, chip, channel),
            Form("Accumulated L%d C%d ch%d; LG (ped-sub) [ADC]; HG (ped-sub) [ADC]",
                 layer, chip, channel),
            cfg_.accum_hist_x_nbins, cfg_.accum_hist_x_min, cfg_.accum_hist_x_max,
            cfg_.accum_hist_y_nbins, cfg_.accum_hist_y_min, cfg_.accum_hist_y_max);
        hist->SetDirectory(nullptr);
        return hist;
    }
    void fill(const AHCALRawHit& h) {
        int cellid = h.cellID;
        if (cfg_.require_hittag && h.hittag != 1) return;
        // Find pedestal for this cell
        auto itp = ped_map_.find(cellid);
        if (itp == ped_map_.end()) {
            n_missing_ped_++;
            return;
        }
        if (itp->second.HighGainStatus != 0 || itp->second.LowGainStatus != 0) {
            n_missing_ped_++;
            return;
        }

        int hg = h.hg_adc;
        int lg = h.lg_adc;
        if (hg_saturation_map_.find(cellid) == hg_saturation_map_.end()) {
            hg_saturation_map_[cellid] = hg;
        } else {
            hg_saturation_map_[cellid] = std::max(static_cast<int>(hg_saturation_map_[cellid]), static_cast<int>(hg));
        }
        if (lg_saturation_map_.find(cellid) == lg_saturation_map_.end()) {
            lg_saturation_map_[cellid] = lg;
        } else {            
            lg_saturation_map_[cellid] = std::max(static_cast<int>(lg_saturation_map_[cellid]), static_cast<int>(lg));
        }
        // Pedestal subtraction
        double hg_sub = hg - itp->second.HighGainPeak;
        double lg_sub = lg - itp->second.LowGainPeak;

        // Only check for valid ADC ranges and saturation, defer cfg-based cuts to fit time
        if (hg <= 0 || lg <= 0 || hg >= 4095 || lg >= 4095) return;

        // Detect saturation for later analysis but don't reject
        if (lg_sub > 200 && hg_sub < cfg_.hg_fit_max) {
            int tmp_layer, tmp_chip, tmp_channel;
            AHCALGeometry::CellIDToLC(cellid, tmp_layer, tmp_chip, tmp_channel);
            if (cfg_.output_bad_cells_json) {
                if (std::find(bad_fit_cells_.begin(), bad_fit_cells_.end(), std::make_tuple(tmp_layer, tmp_chip, tmp_channel)) == bad_fit_cells_.end()) {
                    bad_fit_cells_.emplace_back(std::make_tuple(tmp_layer, tmp_chip, tmp_channel));
                }
            }
        }
        
        // Fill example 2D histogram (all points collected, cuts applied at fit time)
        int tmp_layer, tmp_chip, tmp_channel;
        AHCALGeometry::CellIDToLC(cellid, tmp_layer, tmp_chip, tmp_channel);
        for (size_t i = 0; i < h_example_.size(); ++i) {
            if (cellid == cfg_.example_cellids[i]) {
                h_example_[i]->Fill(lg_sub, hg_sub);
                break;
            }
        }

        // Fill all channels histograms if in save_all_channels mode
        if (cfg_.save_all_channels_root) {
            auto it = h_all_channels_.find(cellid);
            if (it == h_all_channels_.end()) {
                // Create histogram for this channel on first fill
                auto hist = std::make_unique<TH2D>(
                    Form("h_hglg_L%02d_C%02d_ch%02d", tmp_layer, tmp_chip, tmp_channel),
                    Form("L%d C%d ch%d; LG (ped-sub) [ADC]; HG (ped-sub) [ADC]",
                         tmp_layer, tmp_chip, tmp_channel),
                    3000, -100, 2900, 400, -200, 3500);
                hist->SetDirectory(nullptr);
                h_all_channels_[cellid] = std::move(hist);
                all_channels_fit_lines_[cellid] = std::make_unique<TLine>();
                it = h_all_channels_.find(cellid);
            }
            it->second->Fill(lg_sub, hg_sub);
        }

        // Accumulate all points into TH2D histogram (cuts applied later during fit)
        auto it = h_accum_.find(cellid);
        if (it == h_accum_.end()) {
            // Create histogram on first fill
            int tmp_layer2, tmp_chip2, tmp_channel2;
            AHCALGeometry::CellIDToLC(cellid, tmp_layer2, tmp_chip2, tmp_channel2);
            auto hist = createAccumHistogram(tmp_layer2, tmp_chip2, tmp_channel2);
            h_accum_[cellid] = std::move(hist);
            it = h_accum_.find(cellid);
        }
        it->second->Fill(lg_sub, hg_sub);

        n_used_++;
    }
    bool applyFitCut(const AHCALRawHit* h_ptr, double hg_sub, double lg_sub, int tmp_layer, int tmp_chip, int tmp_channel, double hg_saturation) {
        // Apply cfg-based cuts that were deferred from fill()
        if (tmp_layer < 0 || tmp_chip < 0 || tmp_channel < 0) return false;
        if (hg_sub <= cfg_.hg_fit_min) return false;
        if (lg_sub <= cfg_.lg_fit_min) return false;
        if (cfg_.use_fit_max_saturation_minus_margin) {
            double local_hg_fit_max = hg_saturation - cfg_.saturation_margin;
            if (hg_sub > local_hg_fit_max) return false;
        } else {
            if (hg_sub > cfg_.hg_fit_max) return false;
        }
        if (cfg_.require_yperx_over1 && hg_sub <= lg_sub) return false;
        if (cfg_.remove_HGLG_outliers_before_fit) {
            if (hg_sub < 1000 && lg_sub > 50 && hg_sub < 5.0 * lg_sub) return false;
        }
        if (cfg_.use_specific_fit_range_toL39C8) {
            if (tmp_layer == 39 && tmp_chip == 8) {
                if (hg_sub > cfg_.hg_fit_max_L39C8) return false;
            }
        }
        return true;
    }

    void buildInterCalibCache() {
        if (cache_built_) return;
        cache_built_ = true;

        ic_cache_.clear();
        ic_cache_.reserve(h_accum_.size());

        n_fit_ok_ = 0;
        n_fit_all_ = 0;
        double min_p0 = 1e10, max_p0 = -1e10;
        double min_p1 = 1e10, max_p1 = -1e10;

        for (const auto& [cid, hist] : h_accum_) {
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

            // Recalculate saturation detection from histogram
            res.hg_adc_saturation = hg_saturation_map_[cid];
            res.lg_adc_saturation = lg_saturation_map_[cid];

            // Calculate HGLG cluster points regardless of fit status
            int n_hg_lt_2lg = 0;
            int n_noise_points_calc = 0;

            // Extract points from histogram with fit cuts applied
            std::vector<std::pair<double, double>> filtered_points;
            for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
                double hg_center = hist->GetYaxis()->GetBinCenter(iy);
                for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
                    double content = hist->GetBinContent(ix, iy);
                    if (content <= 0) continue;
                    double lg_center = hist->GetXaxis()->GetBinCenter(ix);
                    if (applyFitCut(nullptr, hg_center, lg_center, tmp_layer, tmp_chip, tmp_channel, res.hg_adc_saturation)) {
                        // Add this point 'content' times (content is the bin count)
                        for (int k = 0; k < static_cast<int>(content); ++k) {
                            filtered_points.push_back({hg_center, lg_center});
                            if (hg_center < 2.0 * lg_center) {
                                n_hg_lt_2lg++;
                            }
                        }
                    }
                    if (hg_center < 1000 && lg_center > 50 && hg_center < 5.0 * lg_center) {
                        n_noise_points_calc += static_cast<int>(content);
                    }
                }
            }

            // Mark bad cells for JSON output if saturation detected
            res.n_points = static_cast<long long>(filtered_points.size());
            
            res.n_hg_lt_2lg_aftercut = n_hg_lt_2lg;
            res.n_noise_points = n_noise_points_calc;
            
            if (res.n_points >= cfg_.min_points) {
                n_fit_all_++;
                int n_outlier_points = 0;

                FitResult fr = fitLinearFromHGBinnedProfileRejectOutliers(
                    filtered_points, res.n_points, cfg_.min_points, cfg_.hg_bin_width, cfg_.min_hg_bins_for_fit, cfg_.outlier_sigma_threshold, n_outlier_points, cfg_.min_hg_bins_for_fit);
                res.p0 = fr.p0;
                res.p1 = fr.p1;
                res.chi2 = fr.chi2;
                res.ndf = fr.ndf;
                res.chi2_ndf = fr.chi2_ndf;
                res.fit_status = fr.status;
                res.fit_ok = fr.ok;
                res.n_outlier_points = n_outlier_points;
                res.n_fit_points = fr.n_fit_points;
                res.n_fit_points_over500 = fr.n_fit_points_over500;
                if (fr.ok) {
                    n_fit_ok_++;
                    min_p0 = std::min(min_p0, fr.p0);
                    max_p0 = std::max(max_p0, fr.p0);
                    min_p1 = std::min(min_p1, fr.p1);
                    max_p1 = std::max(max_p1, fr.p1);
                }

                res.distance_avg = -1.0;
                res.distance_rms = -1.0;
                if (fr.ok) {
                    std::vector<double> distances;
                    distances.reserve(filtered_points.size());
                    for (const auto& p : filtered_points) {
                        const double hg_val = p.first;
                        const double lg_val = p.second;
                        // Distance from point to line in HG-LG space: |c*HG - LG + d| / sqrt(c^2 + 1), where c=1/p1 and d=-p0/p1
                        const double c = 1.0 / fr.p1;
                        const double d = -fr.p0 / fr.p1;
                        const double dist = (c * hg_val - lg_val + d) / std::sqrt(c * c + 1.0);
                        distances.push_back(dist);
                    }
                    res.distance_avg = std::accumulate(distances.begin(), distances.end(), 0.0) / distances.size();
                    double sq_sum = std::inner_product(distances.begin(), distances.end(), distances.begin(), 0.0);
                    res.distance_rms = std::sqrt(sq_sum / distances.size());
                    res.distances = std::move(distances);
                }
                // Calculate residuals if fit is ok and save_residual_histograms is enabled
                if (fr.ok && cfg_.save_residual_histograms) {
                    // Convert fit result: HG = p0 + p1*LG -> LG = (HG - p0) / p1
                    const double p1_inv = 1.0 / fr.p1;  // This is equivalent to c in the fit equation
                    const double p0_neg_p1_inv = -fr.p0 / fr.p1;  // This is equivalent to d
                    res.residuals_lg.reserve(filtered_points.size());
                    res.residuals_hg.reserve(filtered_points.size());
                    for (const auto& p : filtered_points) {
                        const double hg_val = p.first;
                        const double lg_val = p.second;
                        const double lg_fitted = p1_inv * hg_val + p0_neg_p1_inv;  // fitted LG = c*HG + d
                        const double residual = lg_val - lg_fitted;
                        res.residuals_lg.push_back(residual);
                        const double hg_fitted = fr.p0 + fr.p1 * lg_val;  // fitted HG = p0 + p1*LG
                        const double residual_hg = hg_val - hg_fitted;
                        res.residuals_hg.push_back(residual_hg);
                    }
                }
            }

            ic_cache_.emplace(cid, std::move(res));
        }
        
        // Adapt 1D histogram ranges based on successful fits
        slope_min_ = (min_p1 < 1e10) ? std::floor(min_p1 * 0.9) : 0.0;
        slope_max_ = (max_p1 > -1e10) ? std::ceil(max_p1 * 1.1) : 40.0;
        intercept_min_ = (min_p0 < 1e10) ? std::floor(min_p0 * 1.1) : -500.0;
        intercept_max_ = (max_p0 > -1e10) ? std::ceil(max_p0 * 1.1) : 500.0;
    }

    void buildQualityResults() {
        quality_results_.clear();
        quality_results_.reserve(ic_cache_.size());

        for (const auto& [cid, res] : ic_cache_) {
            HGLGCalibResult qr;
            int tmp_layer, tmp_chip, tmp_channel;
            AHCALGeometry::CellIDToLC(cid, tmp_layer, tmp_chip, tmp_channel);

            qr.run = ctx_.config.runNumber;
            qr.layer = tmp_layer;
            qr.chip = tmp_chip;
            qr.channel = tmp_channel;
            qr.cellid = cid;

            auto ped_it = ped_map_.find(cid);
            if (ped_it != ped_map_.end()) {
                qr.hg_pedestal = ped_it->second.HighGainPeak;
                qr.lg_pedestal = ped_it->second.LowGainPeak;
                qr.is_pedestal_masked = (ped_it->second.HighGainStatus != 0 || ped_it->second.LowGainStatus != 0);
            }

            qr.slope = res.p1;
            qr.intercept = res.p0;
            qr.chi2 = res.chi2;
            qr.ndf = res.ndf;
            qr.chi2_ndf = res.chi2_ndf;
            qr.fit_status = res.fit_status;

            qr.n_points_total = res.n_points;
            qr.n_points_used = res.n_points - res.n_outlier_points;
            qr.n_points_outlier = res.n_outlier_points;
            qr.n_points_HGLG_cluster = res.n_noise_points;
            if (res.n_points > 0) {
                qr.outlier_fraction = static_cast<double>(res.n_outlier_points) / res.n_points;
            }

            qr.n_fit_points = res.n_fit_points;
            qr.n_fit_points_over500 = res.n_fit_points_over500;

            qr.distance_avg = res.distance_avg;
            qr.distance_rms = res.distance_rms;
            qr.distances = res.distances;
            if (qr.distances.size() == 0) {
                LOG_WARN("InterCalibAlg: cellid {} has zero distances calculated, check fit and point extraction logic", cid);
            }
            // Compute correlation coefficient from accumulated histogram
            auto hist_it = h_accum_.find(cid);
            if (hist_it != h_accum_.end()) {
                std::vector<double> hg_vals, lg_vals;
                for (int iy = 1; iy <= hist_it->second->GetNbinsY(); ++iy) {
                    double hg_center = hist_it->second->GetYaxis()->GetBinCenter(iy);
                    for (int ix = 1; ix <= hist_it->second->GetNbinsX(); ++ix) {
                        double content = hist_it->second->GetBinContent(ix, iy);
                        if (content <= 0) continue;
                        double lg_center = hist_it->second->GetXaxis()->GetBinCenter(ix);
                        for (int k = 0; k < static_cast<int>(content); ++k) {
                            hg_vals.push_back(hg_center);
                            lg_vals.push_back(lg_center);
                        }
                    }
                }
                qr.correlation = computeCorrelation(hg_vals, lg_vals);
            }
            quality_results_.push_back(std::move(qr));
        }
    }

    void computeLayerStatistics() {
        if (quality_results_.empty()) return;

        std::unordered_map<int, std::vector<double>> slopes_by_layer;
        
        for (const auto& qr : quality_results_) {
            // Only include valid fitted channels in layer statistics
            if (qr.is_pedestal_masked) continue;
            if (qr.fit_status != 0) continue;
            if (qr.n_points_used < cfg_.quality_minimum_points) continue;
            if (!std::isfinite(qr.slope)) continue;
            if (qr.slope <= 0) continue;

            slopes_by_layer[qr.layer].push_back(qr.slope);
        }

        // Compute layer median and MAD, then assign z-scores
        std::unordered_map<int, std::pair<double, double>> layer_stats;
        for (auto& [layer, slopes] : slopes_by_layer) {
            auto [median, mad] = computeMedianAndMAD(slopes);
            layer_stats[layer] = {median, mad};
        }

        // Update quality results with layer statistics and z-scores
        for (auto& qr : quality_results_) {
            auto it = layer_stats.find(qr.layer);
            if (it != layer_stats.end()) {
                qr.slope_median_layer = it->second.first;
                qr.slope_mad_layer = it->second.second;

                if (it->second.second > 0.0 && std::isfinite(qr.slope)) {
                    qr.slope_z_layer = (qr.slope - it->second.first) / (1.4826 * it->second.second);
                }
            }
        }
    }

    void assignQualityFlags() {
        for (auto& qr : quality_results_) {
            qr.quality_flag = GOOD;

            if (qr.is_pedestal_masked) {
                qr.quality_flag |= PEDESTAL_MASKED;
            }

            if (qr.n_points_used < cfg_.quality_minimum_points) {
                qr.quality_flag |= LOW_STAT;
            }

            if (qr.fit_status != 0 || !std::isfinite(qr.slope)) {
                qr.quality_flag |= FIT_FAILED;
            }
            if (qr.distance_rms / qr.distance_sigma > cfg_.quality_bad_rms_persigma_threshold) {
                qr.quality_flag |= BAD_RMSPERSIGMA;
            }
            // Determine is_hglg_calib_bad based on flags
            if ((qr.quality_flag & LOW_STAT) ||
                (qr.quality_flag & FIT_FAILED) ||
                (qr.quality_flag & BAD_RMSPERSIGMA)) {
                qr.is_hglg_calib_bad = true;
            }

            // Determine is_physics_mask_candidate: only if pedestal-masked or manually marked bad
            qr.is_physics_mask_candidate = qr.is_pedestal_masked || qr.manual_bad;
        }
    }
    void assignQualityFlags_one(HGLGCalibResult& qr) {
        qr.quality_flag = GOOD;

        if (qr.is_pedestal_masked) {
            qr.quality_flag |= PEDESTAL_MASKED;
        }

        if (qr.n_points_used < cfg_.quality_minimum_points) {
            qr.quality_flag |= LOW_STAT;
        }

        if (qr.fit_status != 0 || !std::isfinite(qr.slope)) {
            qr.quality_flag |= FIT_FAILED;
        }
        if (qr.distance_rms / qr.distance_sigma > cfg_.quality_bad_rms_persigma_threshold) {
            qr.quality_flag |= BAD_RMSPERSIGMA;
        }
        // Determine is_hglg_calib_bad based on flags
        if ((qr.quality_flag & LOW_STAT) ||
            (qr.quality_flag & FIT_FAILED) ||
            (qr.quality_flag & BAD_RMSPERSIGMA)) {
            qr.is_hglg_calib_bad = true;
        }

        // Determine is_physics_mask_candidate: only if pedestal-masked or manually marked bad
        qr.is_physics_mask_candidate = qr.is_pedestal_masked || qr.manual_bad;
    }

    void writeQualityTree() {
        if (!cfg_.quality_output_ttree || quality_results_.empty()) return;

        std::unique_ptr<TFile> fout(TFile::Open(cfg_.quality_out_filename.c_str(), "RECREATE"));
        if (!fout || fout->IsZombie()) {
            LOG_ERROR("InterCalibAlg: cannot open quality output file: {}", cfg_.quality_out_filename);
            return;
        }
    
        TTree tree("HGLGCalibQuality", "HG:LG Calibration Quality Metrics");

        // Temporary variables for branches
        int run, layer, chip, channel, cellid;
        double hg_pedestal, lg_pedestal;
        double slope, intercept, slope_error, intercept_error;
        double chi2, chi2_ndf;
        int ndf, fit_status;
        int n_points_total, n_points_used, n_points_outlier;
        double outlier_fraction;
        double correlation;
        int n_points_HGLG_cluster = 0;
        double distance_avg, distance_rms, distance_sigma;
        double distance_asymmetry;
        double distance_deviation_squared_sum;
        double distance_deviation_squared_sum_fourth;
        double distance_3sig_tail_fraction;
        double distance_5sig_tail_fraction;
        double slope_median_layer, slope_mad_layer, slope_z_layer;
        uint32_t quality_flag;
        int is_pedestal_masked, is_hglg_calib_bad, is_physics_mask_candidate;
        int manual_bad, manual_good;

        int n_fit_points, n_fit_points_over500;


        tree.Branch("run", &run, "run/I");
        tree.Branch("layer", &layer, "layer/I");
        tree.Branch("chip", &chip, "chip/I");
        tree.Branch("channel", &channel, "channel/I");
        tree.Branch("cellid", &cellid, "cellid/I");

        tree.Branch("hg_pedestal", &hg_pedestal, "hg_pedestal/D");
        tree.Branch("lg_pedestal", &lg_pedestal, "lg_pedestal/D");

        tree.Branch("slope", &slope, "slope/D");
        tree.Branch("intercept", &intercept, "intercept/D");
        tree.Branch("slope_error", &slope_error, "slope_error/D");
        tree.Branch("intercept_error", &intercept_error, "intercept_error/D");

        tree.Branch("chi2", &chi2, "chi2/D");
        tree.Branch("ndf", &ndf, "ndf/I");
        tree.Branch("chi2_ndf", &chi2_ndf, "chi2_ndf/D");
        tree.Branch("fit_status", &fit_status, "fit_status/I");

        tree.Branch("n_points_total", &n_points_total, "n_points_total/I");
        tree.Branch("n_points_used", &n_points_used, "n_points_used/I");
        tree.Branch("n_points_outlier", &n_points_outlier, "n_points_outlier/I");
        tree.Branch("outlier_fraction", &outlier_fraction, "outlier_fraction/D");
        tree.Branch("n_points_HGLG_cluster", &n_points_HGLG_cluster, "n_points_HGLG_cluster/I");

        tree.Branch("n_fit_points", &n_fit_points, "n_fit_points/I");
        tree.Branch("n_fit_points_over500", &n_fit_points_over500, "n_fit_points_over500/I");
        
        tree.Branch("distance_avg", &distance_avg, "distance_avg/D");
        tree.Branch("distance_rms", &distance_rms, "distance_rms/D");

        tree.Branch("distance_sigma", &distance_sigma, "distance_sigma/D");
        tree.Branch("distance_asymmetry", &distance_asymmetry, "distance_asymmetry/D");
        tree.Branch("distance_deviation_squared_sum", &distance_deviation_squared_sum, "distance_deviation_squared_sum/D");
        tree.Branch("distance_deviation_squared_sum_fourth", &distance_deviation_squared_sum_fourth, "distance_deviation_squared_sum_fourth/D");
        tree.Branch("distance_3sig_tail_fraction", &distance_3sig_tail_fraction, "distance_3sig_tail_fraction/D");
        tree.Branch("distance_5sig_tail_fraction", &distance_5sig_tail_fraction, "distance_5sig_tail_fraction/D");

        tree.Branch("correlation", &correlation, "correlation/D");

        tree.Branch("slope_median_layer", &slope_median_layer, "slope_median_layer/D");
        tree.Branch("slope_mad_layer", &slope_mad_layer, "slope_mad_layer/D");
        tree.Branch("slope_z_layer", &slope_z_layer, "slope_z_layer/D");

        tree.Branch("quality_flag", &quality_flag, "quality_flag/i");

        tree.Branch("is_pedestal_masked", &is_pedestal_masked, "is_pedestal_masked/I");
        tree.Branch("is_hglg_calib_bad", &is_hglg_calib_bad, "is_hglg_calib_bad/I");
        tree.Branch("is_physics_mask_candidate", &is_physics_mask_candidate, "is_physics_mask_candidate/I");
        tree.Branch("manual_bad", &manual_bad, "manual_bad/I");
        tree.Branch("manual_good", &manual_good, "manual_good/I");
        if (cfg_.quality_output_distance_histograms) {
            fout->mkdir("DistanceHistograms");
            fout->cd("DistanceHistograms");
            for (int i= 0; i <AHCALGeometry::Layer_No; ++i) {
                for (int j = 0; j < AHCALGeometry::chip_No; ++j) {
                    fout->mkdir(Form("L%02d_C%02d", i, j));
                }
            }
        }

        for (const auto& qr : quality_results_) {
            run = qr.run;
            layer = qr.layer;
            chip = qr.chip;
            channel = qr.channel;
            cellid = qr.cellid;

            hg_pedestal = qr.hg_pedestal;
            lg_pedestal = qr.lg_pedestal;

            slope = qr.slope;
            intercept = qr.intercept;
            slope_error = qr.slope_error;
            intercept_error = qr.intercept_error;

            chi2 = qr.chi2;
            ndf = qr.ndf;
            chi2_ndf = qr.chi2_ndf;
            fit_status = qr.fit_status;

            n_points_total = qr.n_points_total;
            n_points_used = qr.n_points_used;
            n_points_outlier = qr.n_points_outlier;
            outlier_fraction = qr.outlier_fraction;

            n_points_HGLG_cluster = qr.n_points_HGLG_cluster;
            n_fit_points = qr.n_fit_points;
            n_fit_points_over500 = qr.n_fit_points_over500;

            distance_avg = qr.distance_avg;
            distance_rms = qr.distance_rms;
            // distance_sigma = qr.distance_sigma;
            if (qr.distances.empty()) {
                distance_sigma = -1.0;
                distance_asymmetry = -1.0;
                distance_deviation_squared_sum = -1.0;
                distance_deviation_squared_sum_fourth = -1.0;
                distance_3sig_tail_fraction = -1.0;
                distance_5sig_tail_fraction = -1.0;
            } 
            if (cfg_.quality_output_distance_histograms && !qr.distances.empty()) {
                std::unique_ptr<TH1D> h_distances = std::make_unique<TH1D>(
                    Form("h_distance_L%02d_C%02d_ch%02d", qr.layer, qr.chip, qr.channel),
                    Form("Distance of points to fit line for L%d C%d ch%d; Distance [ADC]; Entries",
                         qr.layer, qr.chip, qr.channel),
                    40, -10*distance_rms, 10 * distance_rms);
                for (double d : qr.distances) {
                    h_distances->Fill(d);
                }
                h_distances->SetDirectory(nullptr);
                fout->cd("DistanceHistograms");
                fout->cd(Form("L%02d_C%02d", qr.layer, qr.chip));
                h_distances->Write(Form("h_distance_L%02d_C%02d_ch%02d", qr.layer, qr.chip, qr.channel));
                FitOut fo = fitSimpleGaussian(h_distances.get(),200, 2, 5 * distance_rms);
                distance_sigma = fo.sigma;
                distance_asymmetry = fo.asymmetry;
                distance_deviation_squared_sum = fo.deviation_squared_sum;
                distance_deviation_squared_sum_fourth = fo.deviation_squared_sum_fourth;
                distance_3sig_tail_fraction = fo.tail_3sig_fraction;
                distance_5sig_tail_fraction = fo.tail_5sig_fraction;
            }
            

            correlation = qr.correlation;


            slope_median_layer = qr.slope_median_layer;
            slope_mad_layer = qr.slope_mad_layer;
            slope_z_layer = qr.slope_z_layer;
            assignQualityFlags_one(const_cast<HGLGCalibResult&>(qr));  // Update flags based on criteria
            quality_flag = qr.quality_flag;

            is_pedestal_masked = qr.is_pedestal_masked ? 1 : 0;
            is_hglg_calib_bad = qr.is_hglg_calib_bad ? 1 : 0;
            is_physics_mask_candidate = qr.is_physics_mask_candidate ? 1 : 0;
            manual_bad = qr.manual_bad ? 1 : 0;
            manual_good = qr.manual_good ? 1 : 0;

            tree.Fill();
        }

        fout->Write();
        LOG_INFO("InterCalibAlg: wrote quality metrics to {}", cfg_.quality_out_filename);
    }

    void write() {
        if (!cfg_.intercalib_to_file && !cfg_.intercalib_to_json && !cfg_.quality_output_ttree) return;
        if (written_) return;
        written_ = true;

        buildInterCalibCache();

        // Build and compute quality metrics
        if (cfg_.compute_quality_metrics) {
            buildQualityResults();
            computeLayerStatistics();
            assignQualityFlags();
            writeQualityTree();
        }

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

        // 1D distributions with adaptive ranges
        int n_slope_bins = std::max(100, static_cast<int>((slope_max_ - slope_min_) / 0.1));
        int n_intercept_bins = std::max(100, static_cast<int>((intercept_max_ - intercept_min_) / 2.5));
        
        auto hSlope = std::make_unique<TH1D>("hIntercalibSlope_1D",
            "HG/LG slope p1; p1; Entries", n_slope_bins, slope_min_, slope_max_);
        hSlope->SetDirectory(nullptr);

        auto hIntercept = std::make_unique<TH1D>("hIntercalibIntercept_1D",
            "HG/LG intercept p0; p0; Entries", n_intercept_bins, intercept_min_, intercept_max_);
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
        std::vector<std::unique_ptr<TGraphErrors>> g_example_profile(h_example_.size());
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
        int o_n_outlier_points = 0;
        int o_n_hg_lt_2lg_aftercut = 0;
        double o_x_mm = -999.0, o_y_mm = -999.0;
        double hg_rms = 0.0, lg_rms = 0.0;
        double hg_mean = 0.0, lg_mean = 0.0;
        double hg_sigma = 0.0, lg_sigma = 0.0;
        double hg_avg = 0.0, lg_avg = 0.0;
        double hg_chi2 = 0.0, lg_chi2 = 0.0;
        int hg_ndf = -1, lg_ndf = -1;
        double hg_chi2_ndf = 0.0, lg_chi2_ndf = 0.0;
        double hg_3sigma_fraction = 0.0, lg_3sigma_fraction = 0.0;
        double hg_5sigma_fraction = 0.0, lg_5sigma_fraction = 0.0;
        double hg_deviation_squared_sum = 0.0, lg_deviation_squared_sum = 0.0;
        double hg_asymmetry = 0.0, lg_asymmetry = 0.0;
        int n_noise_points = 0;
        double hg_saturation_points = 0;
        double lg_saturation_points = 0.0;
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
        tIC.Branch("n_outlier_points", &o_n_outlier_points, "n_outlier_points/I");
        tIC.Branch("n_hg_lt_2lg_aftercut", &o_n_hg_lt_2lg_aftercut, "n_hg_lt_2lg_aftercut/I");
        tIC.Branch("n_noise_points", &n_noise_points, "n_noise_points/I");
        tIC.Branch("x_mm", &o_x_mm, "x_mm/D");
        tIC.Branch("y_mm", &o_y_mm, "y_mm/D");
        tIC.Branch("hg_mean", &hg_mean, "hg_mean/D");
        tIC.Branch("hg_avg", &hg_avg, "hg_avg/D");
        tIC.Branch("hg_rms", &hg_rms, "hg_rms/D");
        tIC.Branch("hg_sigma", &hg_sigma, "hg_sigma/D");
        tIC.Branch("hg_chi2", &hg_chi2, "hg_chi2/D");
        tIC.Branch("hg_ndf", &hg_ndf, "hg_ndf/I");
        tIC.Branch("hg_chi2_ndf", &hg_chi2_ndf, "hg_chi2_ndf/D");
        tIC.Branch("hg_3sigma_fraction", &hg_3sigma_fraction, "hg_3sigma_fraction/D");
        tIC.Branch("hg_5sigma_fraction", &hg_5sigma_fraction, "hg_5sigma_fraction/D");
        tIC.Branch("hg_deviation_squared_sum", &hg_deviation_squared_sum, "hg_deviation_squared_sum/D");
        tIC.Branch("hg_asymmetry", &hg_asymmetry, "hg_asymmetry/D");
        tIC.Branch("lg_mean", &lg_mean, "lg_mean/D");
        tIC.Branch("lg_avg", &lg_avg, "lg_avg/D");
        tIC.Branch("lg_rms", &lg_rms, "lg_rms/D");
        tIC.Branch("lg_sigma", &lg_sigma, "lg_sigma/D");
        tIC.Branch("lg_chi2", &lg_chi2, "lg_chi2/D");
        tIC.Branch("lg_ndf", &lg_ndf, "lg_ndf/I");
        tIC.Branch("lg_chi2_ndf", &lg_chi2_ndf, "lg_chi2_ndf/D");
        tIC.Branch("lg_3sigma_fraction", &lg_3sigma_fraction, "lg_3sigma_fraction/D");
        tIC.Branch("lg_5sigma_fraction", &lg_5sigma_fraction, "lg_5sigma_fraction/D");
        tIC.Branch("lg_deviation_squared_sum", &lg_deviation_squared_sum, "lg_deviation_squared_sum/D");
        tIC.Branch("lg_asymmetry", &lg_asymmetry, "lg_asymmetry/D");
        tIC.Branch("hg_adc_saturation", &hg_saturation_points, "hg_adc_saturation/D");
        tIC.Branch("lg_adc_saturation", &lg_saturation_points, "lg_adc_saturation/D");
        o_same_data_runs.clear();
        for (const auto& rc : run_contexts_) {
            o_same_data_runs.push_back(rc.config.runNumber);
        }
        // Save residual histograms if enabled
        int n_residual_histograms = 0;
        auto canvases_dir = ensureDir(fout.get(), "ResidualCanvases");
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
            o_n_outlier_points = res.n_outlier_points;
            o_n_hg_lt_2lg_aftercut = res.n_hg_lt_2lg_aftercut;
            n_noise_points = res.n_noise_points;
            hg_saturation_points = res.hg_adc_saturation;
            lg_saturation_points = res.lg_adc_saturation;
            o_x_mm = res.x_mm;
            o_y_mm = res.y_mm;

            hMapP1->SetBinContent(res.layer * AHCALGeometry::chip_No + res.chip + 1, res.channel + 1, res.p1);
            hMapP0->SetBinContent(res.layer * AHCALGeometry::chip_No + res.chip + 1, res.channel + 1, res.p0);

            if (res.fit_ok) {
                hLayerP1[res.layer]->SetBinContent(res.chip + 1, res.channel + 1, res.p1);
                hSlope->Fill(res.p1);
                hIntercept->Fill(res.p0);
            }
            for (size_t i = 0; i < cfg_.example_cellids.size(); ++i) {
                if (res.cellid == cfg_.example_cellids[i]) {
                    fit_results_example[i].p0 = res.p0;
                    fit_results_example[i].p1 = res.p1;
                    fit_results_example[i].status = res.fit_status;
                    fit_results_example[i].ok = res.fit_ok;
                    fit_results_example[i].chi2 = res.chi2;
                    fit_results_example[i].ndf = res.ndf;
                    fit_results_example[i].chi2_ndf = res.chi2_ndf;
                    auto ithist = h_accum_.find(res.cellid);
                    if (ithist != h_accum_.end()) {
                        const auto prof_pts = buildHGProfilePointsFromHist(ithist->second.get(), cfg_.hg_bin_width);
                        if (!prof_pts.empty()) {
                            auto g = std::make_unique<TGraphErrors>(static_cast<int>(prof_pts.size()));
                            int layer_tmp, chip_tmp, ch_tmp;
                            AHCALGeometry::CellIDToLC(res.cellid, layer_tmp, chip_tmp, ch_tmp);
                            g->SetName(Form("gExample_profileFitPoints_L%02d_C%02d_ch%02d",
                                        layer_tmp, chip_tmp, ch_tmp));
                            g->SetTitle(Form("Profile points used for fit; LG (ped-sub) [ADC]; HG (ped-sub) [ADC]"));
                            for (size_t ip = 0; ip < prof_pts.size(); ++ip) {
                                g->SetPoint(ip, prof_pts[ip].lg, prof_pts[ip].hg);
                                g->SetPointError(ip, prof_pts[ip].lg_err, 0.0);
                            }
                            g_example_profile[i] = std::move(g);
                        }
                    }
                    break;
                }
            }
            if (cfg_.save_residual_histograms) {
                if (!res.fit_ok || res.residuals_lg.empty() || res.residuals_hg.empty()) continue;
                std::string dir_Layer_chip = Form("Layer_%02d/Chip_%02d", res.layer, res.chip);
                auto dLayerChip = ensureDir(canvases_dir, dir_Layer_chip.c_str());
                dLayerChip->cd();
                // Create TH1D histogram for residuals
                auto h_residual_lg = std::make_unique<TH1D>(
                    Form("hResidual_LG_L%02d_C%02d_ch%02d", res.layer, res.chip, res.channel),
                    Form("LG Residuals L%d C%d ch%d; Residual [ADC]; Entries",
                         res.layer, res.chip, res.channel),
                    600, -300.0, 300.0);
                h_residual_lg->SetDirectory(nullptr);
                auto h_residual_hg = std::make_unique<TH1D>(
                    Form("hResidual_HG_L%02d_C%02d_ch%02d", res.layer, res.chip, res.channel),
                    Form("HG Residuals L%d C%d ch%d; Residual [ADC]; Entries",
                         res.layer, res.chip, res.channel),
                    2000, -1000.0, 1000.0);
                h_residual_hg->SetDirectory(nullptr);

                // Fill histogram with residuals
                for (double residual : res.residuals_lg) {
                    h_residual_lg->Fill(residual);
                }

                for (double residual : res.residuals_hg) {
                    h_residual_hg->Fill(residual);
                }
                FitOut fit_lg = fitPedestalGaussian(h_residual_lg.get(),200, 2, 2, 0., 200.0);
                FitOut fit_hg = fitPedestalGaussian(h_residual_hg.get(),200, 2, 2, 0., 200.0);
                hg_chi2 = fit_hg.chi2;
                hg_ndf = fit_hg.ndf;
                hg_chi2_ndf = hg_ndf > 0 ? fit_hg.chi2 / fit_hg.ndf : -1.0;
                hg_3sigma_fraction = fit_hg.tail_3sig_fraction;
                hg_mean = fit_hg.mean;
                hg_rms = fit_hg.rms;
                hg_sigma = fit_hg.sigma;
                hg_avg = fit_hg.avg;
                hg_asymmetry = fit_hg.asymmetry;
                lg_chi2 = fit_lg.chi2;
                lg_ndf = fit_lg.ndf;
                lg_chi2_ndf = lg_ndf > 0 ? fit_lg.chi2 / fit_lg.ndf : -1.0;
                lg_3sigma_fraction = fit_lg.tail_3sig_fraction;
                lg_mean = fit_lg.mean;
                lg_rms = fit_lg.rms;
                lg_sigma = fit_lg.sigma;
                lg_avg = fit_lg.avg;
                lg_asymmetry = fit_lg.asymmetry;
                h_residual_lg->Write();
                h_residual_hg->Write();
                n_residual_histograms += 2;
            }
            if (cfg_.output_bad_cells_json && ((o_p1 < 10.0 && o_n > 0 )|| lg_rms > 5 || hg_avg > 40 || hg_avg < -50 || o_n_outlier_points>50 ||hg_chi2_ndf > 40 || lg_chi2_ndf > 40)) {
                bad_fit_cells_.emplace_back(std::make_tuple(res.layer, res.chip, res.channel));
            }
            tIC.Fill();
        }
        if (cfg_.output_bad_cells_json) {
            LOG_INFO("InterCalibAlg: identified {} bad fit cells for JSON output", bad_fit_cells_.size());
        }
        if (cfg_.save_residual_histograms) {
            LOG_INFO("InterCalibAlg: wrote {} residual histograms", n_residual_histograms);
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
                int layer, chip, channel;
                AHCALGeometry::CellIDToLC(cfg_.example_cellids[i], layer, chip, channel);
                TCanvas c(Form("cExample_L%02d_C%02d_ch%02d", layer, chip, channel),
                    Form("Example HG vs LG L%d C%d ch%d", layer, chip, channel),
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
                if (i < g_example_profile.size() && g_example_profile[i]) {
                    g_example_profile[i]->SetMarkerStyle(20);
                    g_example_profile[i]->SetMarkerColor(kBlack);
                    g_example_profile[i]->SetLineColor(kBlack);
                    g_example_profile[i]->Draw("P SAME");
                    g_example_profile[i]->Write();
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

        // Save all channels if enabled
        if (cfg_.save_all_channels_root) {
            for (const auto& [cellid, hist] : h_all_channels_) {
                if (!hist) continue;
                
                auto it_res = ic_cache_.find(cellid);
                if (it_res == ic_cache_.end()) continue;
                
                const auto& res = it_res->second;
                
                // Create canvas for this channel
                TCanvas c(Form("cAll_L%02d_C%02d_ch%02d", res.layer, res.chip, res.channel),
                    Form("All data HG vs LG L%d C%d ch%d", res.layer, res.chip, res.channel),
                    800, 600);
                c.SetRightMargin(0.12);
                c.SetLeftMargin(0.10);
                c.SetBottomMargin(0.10);
                c.SetTopMargin(0.10);
                hist->Draw("COLZ");
                
                // Draw fit line
                if (res.fit_ok) {
                    auto& fit_line = all_channels_fit_lines_[cellid];
                    fit_line->SetLineColor(kRed);
                    fit_line->SetLineWidth(2);
                    fit_line->SetX1(0);
                    fit_line->SetY1(res.p0);
                    fit_line->SetX2(200);
                    fit_line->SetY2(res.p0 + res.p1 * 200);
                    fit_line->Draw("SAME");
                }
                
                // Draw fit info
                TLatex t;
                t.SetNDC(true);
                t.SetTextSize(0.04);
                t.DrawLatex(0.15, 0.80, Form("Fit: p0=%.2f, p1=%.2f, chi2/ndf=%.2f", res.p0, res.p1, res.chi2_ndf));
                
                c.Write();
                hist->Write();
            }
        }
        
        // Save accumulated histograms if requested
        if (cfg_.save_accumulated_histograms && !h_accum_.empty()) {
            TDirectory* accum_dir = ensureDir(fout.get(), "AccumulatedHistograms");
            if (accum_dir) {
                accum_dir->cd();
                for (auto& [cellid, hist] : h_accum_) {
                    if (hist) hist->Write();
                }
                fout->cd();
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
            std::vector<int> n_outlier_points_arr(n_channels_per_layer, 0);

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
                n_outlier_points_arr[idx] = res.n_outlier_points;
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
                j["Summary"]["TotalOutlierPoints"] = std::accumulate(n_outlier_points_arr.begin(), n_outlier_points_arr.end(), 0LL);
                j["Summary"]["NumUsedRun"] = run_contexts_.size();
                j["Summary"]["SameDataRuns"] = same_data_runs;
                j["Status"] = (fit_failures == 0) ? 0 : 1;

                // Per-channel arrays
                j["PerChannel"]["Intercept"] = p0_arr;
                j["PerChannel"]["Slope"] = p1_arr;
                j["PerChannel"]["FitStatus"] = fit_status_arr;
                j["PerChannel"]["NPoints"] = n_points_arr;
                j["PerChannel"]["NOutlierPoints"] = n_outlier_points_arr;

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
            std::vector<int> bad_cellids;
            for (const auto& [layer, chip, channel] : bad_fit_cells_) {
                int cellid = AHCALGeometry::CellID(layer, chip, channel);
                bad_cellids.push_back(cellid);
            }
            j_bad["example_cellids"] = bad_cellids;
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
        const size_t n_examples = cfg_.example_cellids.size();
        if (n_examples == 0 && !cfg_.save_all_channels_root) return;
        
        if (cfg_.save_all_channels_root) {
            LOG_INFO("InterCalibAlg: save_all_channels_root enabled - will save ROOT histograms for all channels");
        }
        
        if (n_examples == 0) return;

        h_example_.clear();
        example_fit_lines_.clear();
        h_example_.reserve(n_examples);
        example_fit_lines_.reserve(n_examples);

        for (size_t i = 0; i < n_examples; ++i) {
            int layer, chip, channel;
            AHCALGeometry::CellIDToLC(cfg_.example_cellids[i], layer, chip, channel);
            h_example_.push_back(std::make_unique<TH2D>(
                Form("h_example_hglg_L%02d_C%02d_ch%02d", layer, chip, channel),
                Form("L%d C%d ch%d (cellid=%d); LG (ped-sub) [ADC]; HG (ped-sub) [ADC]",
                        layer, chip, channel, cfg_.example_cellids[i]),
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
    
    // Adaptive ranges for 1D histograms
    double slope_min_ = 0.0;
    double slope_max_ = 40.0;
    double intercept_min_ = -500.0;
    double intercept_max_ = 500.0;

    std::unordered_map<int, CalibDBIO::Pedestal> ped_map_;
    std::unordered_map<int, std::unique_ptr<TH2D>> h_accum_;  // cellid -> accumulated 2D histogram
    std::unordered_map<int, CellInterCalibResult> ic_cache_;
    std::unordered_map<int, int> hg_saturation_map_; // cellid -> saturation level in HG
    std::unordered_map<int, int> lg_saturation_map_; // cellid -> saturation level in LG 
    std::vector<std::unique_ptr<TH2D>> h_example_;
    std::vector<std::unique_ptr<TLine>> example_fit_lines_;
    std::vector<std::tuple<int, int, int>> bad_fit_cells_;

    // For save_all_channels_root mode
    std::unordered_map<int, std::unique_ptr<TH2D>> h_all_channels_;  // cellid -> 2D histogram
    std::unordered_map<int, std::unique_ptr<TLine>> all_channels_fit_lines_;  // cellid -> fit line
    
    // Quality metrics results
    std::vector<HGLGCalibResult> quality_results_;
    
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
    if (!cfg_.select_muon_hits) {
        for (const auto& h : raw_hits) {
            impl_->fill(h);
        }
    } else {
        if (cfg_.string_track_struct == "SimpleFittedTrack") {
            auto track = evt.get<SimpleFittedTrack>(cfg_.in_track_key);
            if (!cfg_.track_selection_string.empty()) {
                SimpleFittedTrackCut cut(cfg_.track_selection_string);
                if (!cut.eval(track)) {
                    LOG_DEBUG("track did not pass selection cut '{}'", cfg_.track_selection_string);
                    return;
                }
            }
            if (!track.inTrackHitsIndices.size()) {
                LOG_WARN("track has no associated hits.");
                return;
            }
            for (const int index : track.inTrackHitsIndices) {
                for (const auto& rh : raw_hits) {
                    if (rh.index == index) {
                        impl_->fill(rh);
                    }
                }
            }    
        } else if (cfg_.string_track_struct == "Track") {
            auto track = evt.get<Track>(cfg_.in_track_key);
            if (!cfg_.track_selection_string.empty()) {
                TrackCut cut(cfg_.track_selection_string);
                if (!cut.eval(track)) {
                    LOG_DEBUG("track did not pass selection cut '{}'", cfg_.track_selection_string);
                    return;
                }
            }
            if (!track.inTrackHitsIndices.size()) {
                LOG_WARN("no tracks found in the event.");
                return;
            }
            for (const int index : track.inTrackHitsIndices) {
                for (const auto& rh : raw_hits) {
                    if (rh.index == index) {
                        impl_->fill(rh);
                    }   
                }
            }
        } else {
            LOG_ERROR("unknown track struct string '{}'", cfg_.string_track_struct);
            return;
        }
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
    cfg_.use_fit_max_saturation_minus_margin = get_or<bool>(n, "use_fit_max_saturation_minus_margin", cfg_.use_fit_max_saturation_minus_margin);
    if (cfg_.use_fit_max_saturation_minus_margin) {
        LOG_INFO("InterCalibAlg: use_fit_max_saturation_minus_margin is enabled - will set hg_fit_max to (saturation level - margin) for each channel");
        cfg_.saturation_margin = get_or<double>(n, "saturation_margin", cfg_.saturation_margin);
    }
    // cfg_.lg_fit_max = get_or<double>(n, "lg_fit_max", cfg_.lg_fit_max);
    cfg_.hg_fit_min = get_or<double>(n, "hg_fit_min", cfg_.hg_fit_min);
    cfg_.lg_fit_min = get_or<double>(n, "lg_fit_min", cfg_.lg_fit_min);
    cfg_.min_points = get_or<int>(n, "min_points", cfg_.min_points);
    cfg_.hg_bin_width = get_or<double>(n, "hg_bin_width", cfg_.hg_bin_width);
    cfg_.min_hg_bins_for_fit = get_or<int>(n, "min_hg_bins_for_fit", cfg_.min_hg_bins_for_fit);
    cfg_.outlier_sigma_threshold = get_or<double>(n, "outlier_sigma_threshold", cfg_.outlier_sigma_threshold);

    cfg_.example_cellids = get_or<std::vector<int> >(n, "example_cellids", cfg_.example_cellids);
    if (cfg_.example_cellids.empty()) {
        LOG_WARN("InterCalibAlg: example_cellids is empty, no example channels will be visualized");
    }
    cfg_.save_all_channels_root = get_or<bool>(n, "save_all_channels_root", cfg_.save_all_channels_root);
    cfg_.save_residual_histograms = get_or<bool>(n, "save_residual_histograms", cfg_.save_residual_histograms);

    cfg_.require_hittag = get_or<bool>(n, "require_hittag", cfg_.require_hittag);
    cfg_.require_yperx_over1 = get_or<bool>(n, "require_yperx_over1", cfg_.require_yperx_over1);
    cfg_.output_bad_cells_json = get_or<bool>(n, "output_bad_cells_json", cfg_.output_bad_cells_json);
    cfg_.bad_cells_json_filename = get_or<std::string>(n, "bad_cells_json_filename", cfg_.bad_cells_json_filename);
    
    cfg_.use_specific_fit_range_toL39C8 = get_or<bool>(n, "use_specific_fit_range_toL39C8", cfg_.use_specific_fit_range_toL39C8);
    cfg_.hg_fit_max_L39C8 = get_or<double>(n, "hg_fit_max_L39C8", cfg_.hg_fit_max_L39C8);

    cfg_.select_muon_hits = get_or<bool>(n, "select_muon_hits", cfg_.select_muon_hits);
    cfg_.in_track_key = get_or<std::string>(n, "in_track_key", cfg_.in_track_key);
    cfg_.string_track_struct = get_or<std::string>(n, "string_track_struct", cfg_.string_track_struct);
    cfg_.track_selection_string = get_or<std::string>(n, "track_selection_string", cfg_.track_selection_string);

    // Quality metrics and TTree output
    cfg_.compute_quality_metrics = get_or<bool>(n, "compute_quality_metrics", cfg_.compute_quality_metrics);
    cfg_.quality_minimum_points = get_or<int>(n, "quality_minimum_points", cfg_.quality_minimum_points);
    // cfg_.quality_bad_correlation_threshold = get_or<double>(n, "quality_bad_correlation_threshold", cfg_.quality_bad_correlation_threshold);
    // cfg_.quality_outlier_fraction_threshold = get_or<double>(n, "quality_outlier_fraction_threshold", cfg_.quality_outlier_fraction_threshold);
    // cfg_.quality_slope_layer_z_threshold = get_or<double>(n, "quality_slope_layer_z_threshold", cfg_.quality_slope_layer_z_threshold);
    // cfg_.quality_chi2_ndf_threshold = get_or<double>(n, "quality_chi2_ndf_threshold", cfg_.quality_chi2_ndf_threshold);
    cfg_.quality_bad_rms_persigma_threshold = get_or<double>(n, "quality_bad_rms_persigma_threshold", cfg_.quality_bad_rms_persigma_threshold);
    cfg_.quality_output_ttree = get_or<bool>(n, "quality_output_ttree", cfg_.quality_output_ttree);
    cfg_.quality_out_filename = get_or<std::string>(n, "quality_out_filename", cfg_.quality_out_filename);
    cfg_.quality_save_summary_plots = get_or<bool>(n, "quality_save_summary_plots", cfg_.quality_save_summary_plots);
    cfg_.quality_output_distance_histograms = get_or<bool>(n, "quality_output_distance_histograms", cfg_.quality_output_distance_histograms);

    // Interactive mode
    cfg_.interactive_fit = get_or<bool>(n, "interactive_fit", cfg_.interactive_fit);
    cfg_.interactive_mode = get_or<std::string>(n, "interactive_mode", cfg_.interactive_mode);
    cfg_.interactive_layer = get_or<int>(n, "interactive_layer", cfg_.interactive_layer);
    cfg_.interactive_chip = get_or<int>(n, "interactive_chip", cfg_.interactive_chip);
    cfg_.interactive_channel = get_or<int>(n, "interactive_channel", cfg_.interactive_channel);
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
