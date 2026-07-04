/**
 * Standalone script to fit histograms from ROOT file (no draw)
 * Usage: root -b -q 'fit_histograms_noplot_direct_sharded.C("input.root", "output.root")'
 * 
 * Memory Management:
 * - TF1* fLandauGaus: Explicitly deleted inside the fitting helpers
 * - TH1D* histograms: Deleted in loop after processing
 * - TFile*, TTree*: Properly closed and deleted
 * - histograms vector: Cleared after processing all histograms
 * 
 * Persistence on Forced Kill:
 * - TTree is written to disk every 100 entries (every 100 histograms processed)
 * - This ensures partial results are preserved if process is terminated prematurely
 * - Each write uses TObject::kOverwrite flag to update existing entries
 */

#include "langaus.h"
#include <algorithm>
#include <cmath>
#include <TFile.h>
#include <TH1D.h>
#include <TDirectory.h>
#include <TIterator.h>
#include <TKey.h>
#include <TParameter.h>
#include <TTree.h>
#include <string>
#include <unordered_map>
#include <vector>

constexpr const char* kDefaultEfficiencyInput =
    "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out/SimMuon/mip_neighborcheck_nofit.root";

enum FitClass {
    FitClass_LessStatistics = -1,
    FitClass_GoodFit = 0b00000000,
    FitClass_fitFailed = 0b00000001,
    FitClass_BadChi2 =  0b00000010,
    FitClass_ParameterAtLimit = 0b00000100,
    FitClass_gaussigma_small_Threshold = 0b00001000,
    FitClass_gaussigma_small_Retry = 0b00010000
};
enum ReFitClass {
    ReFitClass_NotRefitted = -1,
    ReFitClass_GoodRefit = 0b00000000,
    ReFitClass_RefitFailed = 0b00000001,
    ReFitClass_BadChi2 = 0b00000010,
    ReFitClass_ParameterAtLimit = 0b00000100,
    ReFitClass_LowMPV = 0b00001000,
    ReFitClass_LowGausSigma = 0b00010000,
    ReFitClass_LowWidth = 0b00100000,
    ReFitClass_HighMPVFitError = 0b01000000
};
struct FitOut {
    double mpv = -1.0;
    double mpv_error = -1.0;
    double width = -1.0;
    double width_error = -1.0;
    double total_area = -1.0;
    double total_area_error = -1.0;
    double gaus_sigma = -1.0;
    double gaus_sigma_error = -1.0;
    std::vector<double> covarience_matrix; // Store covariance matrix elements for error analysis
    // |cov(width, width)   covar(width, mpv)   covar(width, total_area)   covar(width, gaus_sigma)|
    // |cov(mpv, width)    covar(mpv, mpv)    covar(mpv, total_area)    covar(mpv, gaus_sigma) |
    // |cov(total_area, width) covar(total_area, mpv) covar(total_area, total_area) covar(total_area, gaus_sigma)|
    // |cov(gaus_sigma, width) covar(gaus_sigma, mpv) covar(gaus_sigma, total_area) covar(gaus_sigma, gaus_sigma)|
    double max_x = -1.0;
    double FWHM = -1.0;
    int entries = 0;
    double chi2 = -1.0;
    int ndf = 0;
    int fit_status = 999;
    bool fit_ok = false;
    double percentile_1 = -1.0;
    double percentile_2 = -1.0;
    double percentile_3 = -1.0;
    double percentile_4 = -1.0;
    double percentile_5 = -1.0;
    double percentile_10 = -1.0;
    double percentile_15 = -1.0;
    double ratio_under_0ADC = -1.0;
    int entries_adc_le50 = 0;
    int fit_Class = -1;
};

struct FitResult {
    int cellid = -1;
    int layer = -1;
    int chip = -1;
    int channel = -1;
    int entries = 0;
    int entries_adc_le50 = 0;
    double ratio_adc_le50 = -1.0;
    double ratio_under_0ADC = -1.0;
    double mpv = -1.0;
    double mpv_error = -1.0;
    double width = -1.0;
    double width_error = -1.0;
    double total_area = -1.0;
    double total_area_error = -1.0;
    double gaus_sigma = -1.0;
    double gaus_sigma_error = -1.0;
    std::vector<double> covarience_matrix;
    double chi2 = -1.0;
    int ndf = 0;
    int fit_status = 999;
    int fit_ok = 0;
    double percentile_1 = -1.0;
    double percentile_2 = -1.0;
    double percentile_3 = -1.0;
    double percentile_4 = -1.0;
    double percentile_5 = -1.0;
    double percentile_10 = -1.0;
    double percentile_15 = -1.0;
    int fit_Class = -1;
    int ntrack_pass_through_channel = -1;
    // Refit results
    bool refit = false;
    double refit_mpv = -1.0;
    double refit_mpv_error = -1.0;
    double refit_width = -1.0;
    double refit_width_error = -1.0;
    double refit_total_area = -1.0;
    double refit_total_area_error = -1.0;
    double refit_gaus_sigma = -1.0;
    double refit_gaus_sigma_error = -1.0;
    std::vector<double> refit_covarience_matrix; // Store covariance matrix elements for error analysis
    double refit_chi2 = -1.0;
    int refit_ndf = 0;
    double refit_chi2_ndf = -1.0;
    double refit_rebin4_chi2 = -1.0;
    int refit_rebin4_ndf = 0;
    double refit_rebin4_chi2_ndf = -1.0;
    int refit_fit_status = 999;
    int refit_fit_ok = 0;
    double refit_percentile_1 = -1.0;
    double refit_percentile_2 = -1.0;
    double refit_percentile_3 = -1.0;
    double refit_percentile_4 = -1.0;
    double refit_percentile_5 = -1.0;
    double refit_percentile_10 = -1.0;
    double refit_percentile_15 = -1.0;
    int refit_Class = -1;
    double refit_ratio_under_0ADC = -1.0;
    double direct_threshold = -1.0;
    double direct_threshold_error = -1.0;
    double direct_width = -1.0;
    double direct_width_error = -1.0;
    double direct_chi2 = -1.0;
    double direct_rebin4_chi2 = -1.0;
    int direct_ndf = 0;
    int direct_rebin4_ndf = 0;
    double direct_chi2_ndf = -1.0;
    double direct_chi2_ndf_bin4 = -1.0;
    int direct_status = 999;
    bool direct_ok = false;
    bool direct_parameter_at_limit = false;
    double direct_efficiency = -1.0;
    double direct_efficiency_range = -1.0;
    double direct_lg_scale = -1.0;
    double direct_target_area = -1.0;
    double direct_actual_area = -1.0;
    double direct_x_min = -1.0;
    double direct_x_max = -1.0;
    int direct_nbins = 0;
};

struct RefitResult {
    double refit_mpv = -1.0;
    double refit_mpv_error = -1.0;
    double refit_width = -1.0;
    double refit_width_error = -1.0;
    double refit_total_area = -1.0;
    double refit_total_area_error = -1.0;
    double refit_gaus_sigma = -1.0;
    double refit_gaus_sigma_error = -1.0;
    std::vector<double> refit_covarience_matrix; // Store covariance matrix elements for error analysis
    double refit_chi2 = -1.0;
    double refit_rebin4_chi2 = -1.0;
    int refit_ndf = 0;
    int refit_rebin4_ndf = 0;
    double refit_chi2_ndf = -1.0;
    double refit_rebin4_chi2_ndf = -1.0;
    int refit_fit_status = 999;
    int refit_fit_ok = 0;
    double refit_percentile_1 = -1.0;
    double refit_percentile_2 = -1.0;
    double refit_percentile_3 = -1.0;
    double refit_percentile_4 = -1.0;
    double refit_percentile_5 = -1.0;
    double refit_percentile_10 = -1.0;
    double refit_percentile_15 = -1.0;
    int refit_Class = -1;
    double refit_ratio_under_0ADC = -1.0;
    double direct_threshold = -1.0;
    double direct_threshold_error = -1.0;
    double direct_width = -1.0;
    double direct_width_error = -1.0;
    double direct_chi2 = -1.0;
    double direct_rebin4_chi2 = -1.0;
    int direct_ndf = 0;
    int direct_rebin4_ndf = 0;
    double direct_chi2_ndf = -1.0;
    double direct_chi2_ndf_bin4 = -1.0;
    int direct_status = 999;
    bool direct_ok = false;
    bool direct_parameter_at_limit = false;
    double direct_efficiency = -1.0;
    double direct_efficiency_range = -1.0;
    double direct_lg_scale = -1.0;
    double direct_target_area = -1.0;
    double direct_actual_area = -1.0;
    double direct_x_min = -1.0;
    double direct_x_max = -1.0;
    int direct_nbins = 0;
};

void calculateFitClassification(FitOut& fitOut) {
    if (fitOut.entries < 200) {
        fitOut.fit_Class = FitClass_LessStatistics;
    } else {
        fitOut.fit_Class = FitClass_GoodFit;  // Default to good fit, will be updated if any issues are found
        if (fitOut.fit_status != 0) {
            fitOut.fit_Class |= FitClass_fitFailed;
        }
        if (fitOut.ndf > 0 && (fitOut.chi2 / double(fitOut.ndf) > 3 || fitOut.chi2 < 0)) {
            fitOut.fit_Class |= FitClass_BadChi2;
        }
        if (fitOut.mpv <= 51.0 || fitOut.width <= 11.0 || fitOut.total_area <= 101.0 || fitOut.gaus_sigma <= 1.0 || fitOut.gaus_sigma >149.0 || fitOut.width > 89.0 || fitOut.mpv < 51.0) {
            fitOut.fit_Class |= FitClass_ParameterAtLimit;
        }
        if (fitOut.gaus_sigma < 20 && fitOut.width <50) {
            fitOut.fit_Class |= FitClass_gaussigma_small_Threshold;
        }
        if (fitOut.gaus_sigma < 20 && fitOut.width >= 50) {
            fitOut.fit_Class |= FitClass_gaussigma_small_Retry;
        }
    }
}

void calculateReFitClassification(RefitResult& fitOut) {
    fitOut.refit_Class = ReFitClass_GoodRefit;
    if (fitOut.refit_fit_status != 0) {
        fitOut.refit_Class |= ReFitClass_RefitFailed;
    }
    if (fitOut.refit_ndf > 0 && (fitOut.refit_chi2 / double(fitOut.refit_ndf) > 9 || fitOut.refit_chi2 < 0)) {
        fitOut.refit_Class |= ReFitClass_BadChi2;
    }
    if (fitOut.refit_mpv < 2 || fitOut.refit_width < 2 || fitOut.refit_gaus_sigma < 2 || fitOut.refit_total_area < 10
        || fitOut.refit_width > 89 || fitOut.refit_gaus_sigma > 149)
    {
        fitOut.refit_Class |= ReFitClass_ParameterAtLimit;
    }
    if (fitOut.refit_mpv < 50) {
        fitOut.refit_Class |= ReFitClass_LowMPV;
    }
    if (fitOut.refit_gaus_sigma < 10) {
        fitOut.refit_Class |= ReFitClass_LowGausSigma;
    }
    if (fitOut.refit_width < 10) {
        fitOut.refit_Class |= ReFitClass_LowWidth;
    }
    if (fitOut.refit_mpv > 0 && fitOut.refit_mpv_error / fitOut.refit_mpv > 0.3) {
        fitOut.refit_Class |= ReFitClass_HighMPVFitError;
    }
}

double thresholdEfficiencyErf(double x, double threshold, double width) {
    if (width <= 0.0) return x >= threshold ? 1.0 : 0.0;
    return 0.5 * (1.0 + TMath::Erf(
        (x - threshold) / (std::sqrt(2.0) * width)));
}

double landauGausBinIntegralContent(TF1* fLandauGaus, TH1D* h, int bin) {
    if (!fLandauGaus || !h || bin < 1 || bin > h->GetNbinsX()) return 0.0;

    const double x_low = h->GetBinLowEdge(bin);
    const double bin_width = h->GetBinWidth(bin);
    if (bin_width <= 0.0) return 0.0;
    return fLandauGaus->Integral(x_low, x_low + bin_width) / bin_width;
}

double sumLandauGausIntegralInRange(
    TF1* fLandauGaus, TH1D* h, double x_min, double x_max, int* n_bins = nullptr
) {
    if (!fLandauGaus || !h) return -1.0;

    double sum = 0.0;
    int count = 0;
    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        const double x = h->GetBinCenter(bin);
        if (x <= x_min || x > x_max) continue;
        sum += landauGausBinIntegralContent(fLandauGaus, h, bin);
        ++count;
    }
    if (n_bins) *n_bins = count;
    return sum;
}

TF1* makeDirectThresholdModelIntegral(
    TF1* fLandauGaus,
    TH1D* h,
    double lg_scale,
    double x_min,
    double x_max,
    const char* name
) {
    if (!fLandauGaus || !h) return nullptr;

    std::vector<double> scaled_bin_integrals(h->GetNbinsX() + 2, 0.0);
    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        scaled_bin_integrals[bin] =
            lg_scale * landauGausBinIntegralContent(fLandauGaus, h, bin);
    }

    auto model = [h, scaled_bin_integrals](double* xx, double* par) {
        const int bin = h->FindBin(xx[0]);
        const double fit_bin =
            bin >= 0 && bin < static_cast<int>(scaled_bin_integrals.size())
                ? scaled_bin_integrals[bin]
                : 0.0;
        return fit_bin * thresholdEfficiencyErf(xx[0], par[0], par[1]);
    };

    TF1* f_model = new TF1(name, model, x_min, x_max, 2);
    f_model->SetParName(0, "threshold");
    f_model->SetParName(1, "width");
    return f_model;
}

double computeDirectThresholdEfficiencyIntegral(
    TF1* fLandauGaus,
    TH1D* h,
    double lg_scale,
    double threshold,
    double width,
    double x_min,
    double x_max
) {
    if (!fLandauGaus || !h || width <= 0.0) return -1.0;

    double denominator = 0.0;
    double numerator = 0.0;
    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        const double x = h->GetBinCenter(bin);
        if (x <= x_min || x > x_max) continue;
        const double fit =
            lg_scale * landauGausBinIntegralContent(fLandauGaus, h, bin);
        denominator += fit;
        numerator += fit * thresholdEfficiencyErf(x, threshold, width);
    }
    return denominator > 0.0 ? numerator / denominator : -1.0;
}


void calculatePercentiles(TH1D* h, FitOut& result, bool count_entries_below_threshold) {
    if (!h) return;

    double cumulative_entries = 0.0;
    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        const double bin_center = h->GetBinCenter(bin);
        const double bin_content = h->GetBinContent(bin);
        if (bin_center <= 50.0) {
            if (count_entries_below_threshold) {
                result.entries_adc_le50 += bin_content;
            }
            continue;
        }

        cumulative_entries += bin_content;
        if (result.percentile_1 < 0.0 && cumulative_entries >= result.entries * 0.01) {
            result.percentile_1 = bin_center;
        }
        if (result.percentile_2 < 0.0 && cumulative_entries >= result.entries * 0.02) {
            result.percentile_2 = bin_center;
        }
        if (result.percentile_3 < 0.0 && cumulative_entries >= result.entries * 0.03) {
            result.percentile_3 = bin_center;
        }
        if (result.percentile_4 < 0.0 && cumulative_entries >= result.entries * 0.04) {
            result.percentile_4 = bin_center;
        }
        if (result.percentile_5 < 0.0 && cumulative_entries >= result.entries * 0.05) {
            result.percentile_5 = bin_center;
        }
        if (result.percentile_10 < 0.0 && cumulative_entries >= result.entries * 0.10) {
            result.percentile_10 = bin_center;
        }
        if (cumulative_entries >= result.entries * 0.15) {
            result.percentile_15 = bin_center;
            break;
        }
    }
}

FitOut fitLandauGaus(TH1D* h, int minEntries = 200, bool calculate_fwhm = false, double min_gaus_sigma = 0.0, int fix_total_area = -1) {
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
    fr[0] = 50; fr[1] = 1000;
    if (h->GetBinCenter(h->GetMaximumBin()) > 350) {
        fr[0] = 50; fr[1] = 1500;
    }
    TH1D* h_copy_0 = (TH1D*)h->Clone("h_copy_0");
    h_copy_0->SetDirectory(nullptr);
    h_copy_0->Rebin(4);
    double x_max = h_copy_0->GetBinCenter(h_copy_0->GetMaximumBin());
    delete h_copy_0;

    pllo[0] = 10;    pllo[1] = 50;   pllo[2] = 100;  pllo[3] = min_gaus_sigma > 0.0 ? min_gaus_sigma : 0;
    plhi[0] = 90;  plhi[1] = 1200; plhi[2] = 100000; plhi[3] = 150;
    sv[0] = 50; sv[1] = x_max; sv[2] = 30000; sv[3] = 60;
    if (fix_total_area > 0) {
        sv[2] = fix_total_area;
        pllo[2] = fix_total_area;
        plhi[2] = fix_total_area;
    }
    int fitResult = 0;
    std::vector<double> covariance_matrix(16, 0.0); // Initialize covariance matrix with 16 elements (4x4)
    TF1* fLandauGaus = langaufit(h, fr, sv, pllo, plhi, fps, fpe, &chisqr, &ndf, &fitResult, &covariance_matrix[0]);
    if (!fLandauGaus) {
        std::cout << "Landau-Gaussian fit failed for histogram: " << h->GetName() << std::endl; 
        return r;
    }
    
    double maxx = 0, fwhm = -1.0;
    if (calculate_fwhm) {
        langaupro(fps, maxx, fwhm);
    }

    r.mpv = fps[1];
    r.mpv_error = fpe[1];
    r.width = fps[0];
    r.width_error = fpe[0];
    r.total_area = fps[2];
    r.total_area_error = fpe[2];
    r.gaus_sigma = fps[3];
    r.gaus_sigma_error = fpe[3];
    r.covarience_matrix = covariance_matrix; // Store covariance matrix for error analysis
    r.max_x = calculate_fwhm ? maxx : -1.0;
    r.FWHM = fwhm;
    r.entries = h->GetEntries();
    r.chi2 = chisqr;
    r.ndf = ndf;
    
    calculatePercentiles(h, r, true);

    double chi2_ndf = (ndf > 0) ? (chisqr / double(ndf)) : -1.0;
    if (chi2_ndf > 3 || chi2_ndf < 0) {
        r.fit_status = fitResult;
        r.fit_ok = false;
    } else {
        r.fit_status = fitResult;
        r.fit_ok = true;
    }
    calculateFitClassification(r);
    r.ratio_under_0ADC = (r.total_area > 0) ? double(fLandauGaus->Integral(-96, 0)) / double(r.total_area) : -1.0;
    delete fLandauGaus;
    fLandauGaus = nullptr;
    if (r.gaus_sigma < 1) {
        const FitOut retry = fitLandauGaus(h, minEntries, calculate_fwhm, 20);
        if (retry.fit_ok && retry.ndf > 0 && ndf > 0 &&
            retry.chi2 / retry.ndf < chisqr / ndf && retry.gaus_sigma > 20.1) {
            r = retry;
        }
    }
    return r;
}

RefitResult fitLandauGausWithThreshold(TH1D* h, int minEntries, double min_gaus_sigma, FitResult& res, double TotalArea ) {
    RefitResult rr;
    FitOut r;
    if (!h) return rr;
    if (h->GetEntries() < minEntries) {
        r.fit_status = -2;
        r.fit_ok = false;
        return rr;
    }

    double fr[2];
    double sv[4], pllo[4], plhi[4], fps[4], fpe[4];
    double chisqr;
    int ndf;
    TH1D* h_copy_0 = (TH1D*)h->Clone("h_copy_0");
    h_copy_0->SetDirectory(nullptr);
    h_copy_0->Rebin(4);
    double x_max = h_copy_0->GetBinCenter(h_copy_0->GetMaximumBin());
    delete h_copy_0;

    pllo[0] = 1;    pllo[1] = 1;   pllo[2] = 100;  pllo[3] = min_gaus_sigma > 0.0 ? min_gaus_sigma : 0;
    plhi[0] = 90;  plhi[1] = 1200; plhi[2] = 100000; plhi[3] = 150;
    sv[0] = 50; sv[1] = x_max; sv[2] = 30000; sv[3] = 60;
    sv[2] = TotalArea > 0 ? TotalArea*h->GetBinWidth(1)  : 30000;
    pllo[2] = TotalArea > 0 ? TotalArea*h->GetBinWidth(1)  : 100;
    plhi[2] = TotalArea > 0 ? TotalArea*h->GetBinWidth(1) : 100000;
    if (res.mpv <51){
        fr[0] = x_max > 0 ? x_max : 50;
        fr[1] = x_max > 0 ? x_max * 4 : 1000;
    }else{
        fr[0] = res.mpv > 0 ? res.mpv : 50;
        fr[1] = res.mpv > 0 ? res.mpv * 4 : 1000;
    }
    int fitResult = 0;
    std::vector<double> covariance_matrix(16, 0.0); // Initialize covariance matrix with 16 elements (4x4)
    TF1* fLandauGaus = langaufit(h, fr, sv, pllo, plhi, fps, fpe, &chisqr, &ndf, &fitResult, &covariance_matrix[0]);
    if (!fLandauGaus) {
        r.fit_status = -3;
        r.fit_ok = false;
        return rr;
    }
    r.mpv = fps[1];
    r.mpv_error = fpe[1];
    r.width = fps[0];
    r.width_error = fpe[0];
    r.total_area = fps[2];
    r.total_area_error = fpe[2];
    r.gaus_sigma = fps[3];
    r.gaus_sigma_error = fpe[3];
    r.covarience_matrix = covariance_matrix; // Store covariance matrix for error analysis
    r.entries = h->GetEntries();
    r.chi2 = chisqr;
    // recalculate chi2 with other rebin 4.
    TH1D* h_rebin = (TH1D*)h->Rebin(4, "h_rebin");
    double sum_chi2 = 0;
    int rebin_ndf = 0;
    for (int i = 1; i <= h_rebin->GetNbinsX(); ++i) {
        if (h_rebin->GetBinCenter(i) > fr[1]) continue;  // Skip bins outside fit range
        if (h_rebin->GetBinCenter(i) < fr[0]) continue;  // Skip bins outside fit range
        double expected = fLandauGaus->Eval(h_rebin->GetBinCenter(i)) * h_rebin->GetBinWidth(1)/h->GetBinWidth(1);  // Scale expected value by rebinning factor
        double observed = h_rebin->GetBinContent(i);
        if (expected > 0) {
            sum_chi2 += std::pow(observed - expected, 2) / expected;
        }
        rebin_ndf++;
    }
    rebin_ndf -= fLandauGaus->GetNpar();  // Subtract number of fit parameters from ndf
    rr.refit_rebin4_chi2 = sum_chi2;
    rr.refit_rebin4_ndf = rebin_ndf;
    rr.refit_rebin4_chi2_ndf = (rebin_ndf > 0) ? (sum_chi2 / double(rebin_ndf)) : -1.0;
    delete h_rebin;
    r.ndf = ndf;
    r.fit_status = fitResult;
    const double chi2_ndf = ndf > 0 ? chisqr / double(ndf) : -1.0;
    r.fit_ok = fitResult == 0;
    r.ratio_under_0ADC = (r.total_area > 0) ? double(fLandauGaus->Integral(-96, 0)) / double(r.total_area) : -1.0;
    calculatePercentiles(h, r, false);

    calculateFitClassification(r);
    // Fill RefitResult struct
    rr.refit_mpv = r.mpv;
    rr.refit_mpv_error = r.mpv_error;
    rr.refit_width = r.width;
    rr.refit_width_error = r.width_error;
    rr.refit_total_area = r.total_area;
    rr.refit_total_area_error = r.total_area_error;
    rr.refit_gaus_sigma = r.gaus_sigma;
    rr.refit_gaus_sigma_error = r.gaus_sigma_error;
    rr.refit_covarience_matrix = r.covarience_matrix;
    rr.refit_chi2 = r.chi2;
    rr.refit_ndf = r.ndf;
    rr.refit_chi2_ndf = chi2_ndf;
    rr.refit_fit_status = r.fit_status;
    rr.refit_fit_ok = r.fit_ok ? 1 : 0;
    rr.refit_percentile_1 = r.percentile_1;
    rr.refit_percentile_2 = r.percentile_2;
    rr.refit_percentile_3 = r.percentile_3;
    rr.refit_percentile_4 = r.percentile_4;
    rr.refit_percentile_5 = r.percentile_5;
    rr.refit_percentile_10 = r.percentile_10;
    rr.refit_percentile_15 = r.percentile_15;
    rr.refit_Class = r.fit_Class;
    rr.refit_ratio_under_0ADC = r.ratio_under_0ADC;
    if (rr.refit_Class & ReFitClass_RefitFailed) {
        delete fLandauGaus;
        fLandauGaus = nullptr;
        return rr;
    }

    rr.direct_x_min = 50.0;
    rr.direct_x_max = res.mpv > 0.0
        ? std::max(
            rr.direct_x_min + h->GetBinWidth(1),
            std::min(res.mpv * 4.0, 2000.0))
        : 2000.0;
    rr.direct_target_area = TotalArea;
    rr.direct_actual_area = sumLandauGausIntegralInRange(
        fLandauGaus, h, rr.direct_x_min, rr.direct_x_max,
        &rr.direct_nbins);
    if (rr.direct_target_area > 0.0 && rr.direct_actual_area > 0.0) {
        rr.direct_lg_scale =
            rr.direct_target_area / rr.direct_actual_area;
        TF1* f_direct = makeDirectThresholdModelIntegral(
            fLandauGaus,
            h,
            1.0,
            rr.direct_x_min,
            rr.direct_x_max,
            Form("f_direct_threshold_%d", res.cellid));
        if (f_direct) {
            double threshold_initial =
                res.percentile_1 > 0.0 ? res.percentile_1 : 0.6 * r.mpv;
            threshold_initial = std::max(
                rr.direct_x_min,
                std::min(threshold_initial, rr.direct_x_max));
            f_direct->SetParameters(threshold_initial, 20.0);
            f_direct->SetParLimits(
                0, rr.direct_x_min, rr.direct_x_max);
            f_direct->SetParLimits(1, 3.0, 100.0);

            rr.direct_status = h->Fit(f_direct, "RQN");
            rr.direct_threshold = f_direct->GetParameter(0);
            rr.direct_threshold_error = f_direct->GetParError(0);
            rr.direct_width = f_direct->GetParameter(1);
            rr.direct_width_error = f_direct->GetParError(1);
            rr.direct_chi2 = f_direct->GetChisquare();
            // rr.direct_chi2_bin4
            double chi2_bin4 = 0.0;
            int ndf_bin4 = 0;
            TH1D* h_rebin4 = (TH1D*)h->Rebin(4, "h_rebin4");
            for (int bin = 1; bin <= h_rebin4->GetNbinsX(); ++bin) {
                double x = h_rebin4->GetBinCenter(bin);
                if (x < rr.direct_x_min || x > rr.direct_x_max) continue;
                double expected = f_direct->Integral(x - h_rebin4->GetBinWidth(1)/2, x + h_rebin4->GetBinWidth(1)/2) / h_rebin4->GetBinWidth(1);  // Scale expected value by rebinning factor
                chi2_bin4 += (h_rebin4->GetBinContent(bin) - expected) * (h_rebin4->GetBinContent(bin) - expected) / (h_rebin4->GetBinError(bin) * h_rebin4->GetBinError(bin));
                ++ndf_bin4;
            }
            delete h_rebin4;
            rr.direct_rebin4_chi2 = chi2_bin4;
            rr.direct_rebin4_ndf = ndf_bin4 - f_direct->GetNpar();
            rr.direct_ndf = f_direct->GetNDF();
            rr.direct_chi2_ndf_bin4 = rr.direct_rebin4_chi2 / double(rr.direct_rebin4_ndf);
            rr.direct_chi2_ndf =
                rr.direct_ndf > 0
                    ? rr.direct_chi2 / double(rr.direct_ndf)
                    : -1.0;
            rr.direct_ok =
                rr.direct_status == 0 && rr.direct_ndf > 0;
            rr.direct_parameter_at_limit =
                rr.direct_threshold < rr.direct_x_min + 1.0 ||
                rr.direct_threshold > rr.direct_x_max - 1.0 ||
                rr.direct_width < 3.0 + 1.0 ||
                rr.direct_width > 100.0 - 1.0;
            rr.direct_efficiency_range =
                computeDirectThresholdEfficiencyIntegral(
                    fLandauGaus,
                    h,
                    rr.direct_lg_scale,
                    rr.direct_threshold,
                    rr.direct_width,
                    rr.direct_x_min,
                    rr.direct_x_max);

            double totalArea = fLandauGaus->Integral(-96, 4000);
            double totalAreaWithEff = f_direct->Integral(-96, 4000);
            rr.direct_efficiency =
                totalArea > 0.0 ? totalAreaWithEff / totalArea : -1.0;
            delete f_direct;
        }
    }

    calculateReFitClassification(rr);

    delete fLandauGaus;
    fLandauGaus = nullptr;
    
    return rr;
}

void collectHistograms(
    TDirectory* dir,
    std::vector<TH1D*>& histograms,
    std::unordered_map<int, int>& ntracks_by_cell,
    int shard_index,
    int shard_count
) {
    if (!dir) return;
    
    TIter next(dir->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        TObject* obj = dir->Get(key->GetName());
        if (!obj) continue;
        
        if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
            TDirectory* subdir = (TDirectory*)obj;
            collectHistograms(subdir, histograms, ntracks_by_cell, shard_index, shard_count);
        } else if (obj->IsA()->InheritsFrom(TH1D::Class())) {
            TH1D* h = (TH1D*)obj;
            int cellid = -1;
            if (sscanf(h->GetName(), "hMIP_%d", &cellid) != 1 ||
                cellid % shard_count != shard_index) {
                continue;
            }
            h->Rebin(4);
            h->SetDirectory(nullptr);
            histograms.push_back(h);
        } else if (obj->IsA()->InheritsFrom(TParameter<int>::Class())) {
            TParameter<int>* p = (TParameter<int>*)obj;
            std::string name = p->GetName();
            if (name.find("Ntracks_pass_through_channel_") == 0) {
                int cellid = std::stoi(name.substr(std::string("Ntracks_pass_through_channel_").length()));
                if (ntracks_by_cell.find(cellid) == ntracks_by_cell.end()) {
                    ntracks_by_cell[cellid] = p->GetVal();
                }
            }
        }
    }
}

bool isRefitAHCALCell(int cellID) {
    const int layer = cellID / 100000;
    const int chip = (cellID / 10000) % 10;
    // Bad full layers
    if (layer == 1 || layer == 3 || layer == 4 || layer == 5 || layer == 39) {
        return true;
    }

    // Bad specific layer-chip pairs
    if ((layer == 12 && chip == 1) ||
        (layer == 15 && chip == 0) ||
        (layer == 18 && chip == 4) ||
        (layer == 25 && chip == 0) ||
        (layer == 27 && (chip == 7 || chip == 8)) ||
        (layer == 32 && (chip == 1 || chip == 2)) ||
        (layer == 33 && (chip == 0 || chip == 1 || chip == 2))) {
        return true;
    }

    return false;
}

struct EfficiencyInput {
    std::unordered_map<int, Long64_t> entries;
    std::unordered_map<int, int> ntrack_pass_through;
};

void collectEfficiencyInput(TDirectory* dir, EfficiencyInput& input) {
    if (!dir) return;

    TIter next(dir->GetListOfKeys());
    TKey* key = nullptr;
    while ((key = static_cast<TKey*>(next()))) {
        TObject* obj = dir->Get(key->GetName());
        if (!obj) continue;

        if (obj->InheritsFrom(TDirectory::Class())) {
            collectEfficiencyInput(static_cast<TDirectory*>(obj), input);
        } else if (obj->InheritsFrom(TH1::Class())) {
            const std::string name = obj->GetName();
            const std::string prefix = "hMIP_";
            if (name.find(prefix) == 0) {
                const int cellid = std::stoi(name.substr(prefix.size()));
                input.entries[cellid] =
                    static_cast<TH1*>(obj)->GetEntries();
            }
        } else if (obj->InheritsFrom(TParameter<int>::Class())) {
            const std::string name = obj->GetName();
            const std::string prefix = "Ntracks_pass_through_channel_";
            if (name.find(prefix) == 0) {
                const int cellid = std::stoi(name.substr(prefix.size()));
                input.ntrack_pass_through[cellid] =
                    static_cast<TParameter<int>*>(obj)->GetVal();
            }
        }

    }
}

std::unordered_map<int, double> createTotalAreaMap(const char* input_file = kDefaultEfficiencyInput) {
    std::unordered_map<int, double> totalAreaMap;
    TFile* input_file_ptr = TFile::Open(input_file, "READ");
    if (!input_file_ptr || input_file_ptr->IsZombie()) {
        std::cerr << "Error: cannot open input file: " << input_file
                  << std::endl;
        delete input_file_ptr;
        return totalAreaMap;
    }

    EfficiencyInput input;
    TDirectory* mip_dir = input_file_ptr->GetDirectory("MIP");
    collectEfficiencyInput(mip_dir ? mip_dir : input_file_ptr, input);

    for (const auto& item : input.entries) {
        const auto ntracks_it = input.ntrack_pass_through.find(item.first);
        const int ntracks =
            ntracks_it != input.ntrack_pass_through.end() ? ntracks_it->second : -1;
        totalAreaMap[item.first] =
            ntracks > 0 ? static_cast<double>(item.second) / ntracks : -1.0;
    }

    input_file_ptr->Close();
    delete input_file_ptr;

    return totalAreaMap;
}

void copyFitResult(FitResult& result, const FitOut& fit) {
    result.mpv = fit.mpv;
    result.mpv_error = fit.mpv_error;
    result.width = fit.width;
    result.width_error = fit.width_error;
    result.total_area = fit.total_area;
    result.total_area_error = fit.total_area_error;
    result.gaus_sigma = fit.gaus_sigma;
    result.gaus_sigma_error = fit.gaus_sigma_error;
    result.covarience_matrix = fit.covarience_matrix;
    result.chi2 = fit.chi2;
    result.ndf = fit.ndf;
    result.fit_status = fit.fit_status;
    result.fit_ok = fit.fit_ok ? 1 : 0;
    result.entries_adc_le50 = fit.entries_adc_le50;
    result.percentile_1 = fit.percentile_1;
    result.percentile_2 = fit.percentile_2;
    result.percentile_3 = fit.percentile_3;
    result.percentile_4 = fit.percentile_4;
    result.percentile_5 = fit.percentile_5;
    result.percentile_10 = fit.percentile_10;
    result.percentile_15 = fit.percentile_15;
    result.ratio_under_0ADC = fit.ratio_under_0ADC;
    result.fit_Class = fit.fit_Class;
}

void copyRefitResult(FitResult& result, const RefitResult& refit) {
    result.refit = true;
    result.refit_mpv = refit.refit_mpv;
    result.refit_mpv_error = refit.refit_mpv_error;
    result.refit_width = refit.refit_width;
    result.refit_width_error = refit.refit_width_error;
    result.refit_total_area = refit.refit_total_area;
    result.refit_total_area_error = refit.refit_total_area_error;
    result.refit_gaus_sigma = refit.refit_gaus_sigma;
    result.refit_gaus_sigma_error = refit.refit_gaus_sigma_error;
    result.refit_covarience_matrix = refit.refit_covarience_matrix;
    result.refit_chi2 = refit.refit_chi2;
    result.refit_ndf = refit.refit_ndf;
    result.refit_chi2_ndf = refit.refit_chi2_ndf;
    result.refit_rebin4_chi2 = refit.refit_rebin4_chi2;
    result.refit_rebin4_ndf = refit.refit_rebin4_ndf;
    result.refit_rebin4_chi2_ndf = refit.refit_rebin4_chi2_ndf;
    result.refit_fit_status = refit.refit_fit_status;
    result.refit_fit_ok = refit.refit_fit_ok;
    result.refit_percentile_1 = refit.refit_percentile_1;
    result.refit_percentile_2 = refit.refit_percentile_2;
    result.refit_percentile_3 = refit.refit_percentile_3;
    result.refit_percentile_4 = refit.refit_percentile_4;
    result.refit_percentile_5 = refit.refit_percentile_5;
    result.refit_percentile_10 = refit.refit_percentile_10;
    result.refit_percentile_15 = refit.refit_percentile_15;
    result.refit_Class = refit.refit_Class;
    result.refit_ratio_under_0ADC = refit.refit_ratio_under_0ADC;
    result.direct_threshold = refit.direct_threshold;
    result.direct_threshold_error = refit.direct_threshold_error;
    result.direct_width = refit.direct_width;
    result.direct_width_error = refit.direct_width_error;
    result.direct_chi2 = refit.direct_chi2;
    result.direct_ndf = refit.direct_ndf;
    result.direct_chi2_ndf = refit.direct_chi2_ndf;
    result.direct_chi2_ndf_bin4 = refit.direct_chi2_ndf_bin4;
    result.direct_status = refit.direct_status;
    result.direct_ok = refit.direct_ok;
    result.direct_parameter_at_limit = refit.direct_parameter_at_limit;
    result.direct_efficiency = refit.direct_efficiency;
    result.direct_efficiency_range = refit.direct_efficiency_range;
    result.direct_lg_scale = refit.direct_lg_scale;
    result.direct_target_area = refit.direct_target_area;
    result.direct_actual_area = refit.direct_actual_area;
    result.direct_x_min = refit.direct_x_min;
    result.direct_x_max = refit.direct_x_max;
    result.direct_nbins = refit.direct_nbins;
}

bool shouldRefit(int fit_class, int cellid, double efficiency_ratio) {
    return fit_class > FitClass_GoodFit || ( isRefitAHCALCell(cellid) && efficiency_ratio < 1.0);
}

void fit_histograms_noplot_direct_sharded(
    const char* input_file = "mip.root",
    const char* output_file = "mip_fitted.root",
    const char* efficiency_input_file = kDefaultEfficiencyInput,
    int shard_index = 0,
    int shard_count = 1
) {
    if (shard_count <= 0 || shard_index < 0 || shard_index >= shard_count) {
        std::cerr << "Invalid shard selection: index=" << shard_index
                  << ", count=" << shard_count << std::endl;
        return;
    }

    std::cout << "Opening input file (no plot mode): " << input_file << std::endl;
    std::cout << "Processing shard " << shard_index << "/" << shard_count << std::endl;
    TFile* fin = TFile::Open(input_file, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: cannot open input file: " << input_file << std::endl;
        return;
    }

    std::vector<TH1D*> histograms;
    std::unordered_map<int, int> NtrackPassThroughChannel;
    std::cout << "Collecting histograms..." << std::endl;
    TDirectory* mipDir = (TDirectory*)fin->Get("MIP");
    if (mipDir) {
        collectHistograms(
            mipDir, histograms, NtrackPassThroughChannel, shard_index, shard_count);
    } else {
        collectHistograms(
            fin, histograms, NtrackPassThroughChannel, shard_index, shard_count);
    }
    std::cout << "Total histograms found: " << histograms.size() << std::endl;

    TFile* fout = TFile::Open(output_file, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: cannot create output file: " << output_file << std::endl;
        fin->Close();
        delete fin;
        return;
    }
    std::unordered_map<int, double> totalAreaMap;
    bool total_area_map_loaded = false;
    fout->cd();

    // Create result tree
    TTree* resultTree = new TTree("mip_fit_results", "MIP fit results");
    FitResult res;
    double MC_efficiency = -1.0;
    double efficiency = -1.0;
    double eff_ratio = -1.0;
    resultTree->Branch("cellid", &res.cellid);
    resultTree->Branch("layer", &res.layer);
    resultTree->Branch("chip", &res.chip);
    resultTree->Branch("channel", &res.channel);
    resultTree->Branch("entries", &res.entries);
    resultTree->Branch("entries_adc_le50", &res.entries_adc_le50);
    resultTree->Branch("ratio_adc_le50", &res.ratio_adc_le50);
    resultTree->Branch("MC_efficiency", &MC_efficiency);
    resultTree->Branch("efficiency", &efficiency);
    resultTree->Branch("efficiency_ratio", &eff_ratio);
    resultTree->Branch("MPV", &res.mpv);
    resultTree->Branch("MPV_error", &res.mpv_error);
    resultTree->Branch("width", &res.width);
    resultTree->Branch("width_error", &res.width_error);
    resultTree->Branch("TotalArea", &res.total_area);
    resultTree->Branch("TotalArea_error", &res.total_area_error);
    resultTree->Branch("gaus_sigma", &res.gaus_sigma);
    resultTree->Branch("gaus_sigma_error", &res.gaus_sigma_error);
    resultTree->Branch("covarience_matrix", &res.covarience_matrix);
    resultTree->Branch("chi2", &res.chi2);
    resultTree->Branch("ndf", &res.ndf);
    resultTree->Branch("fit_status", &res.fit_status);
    resultTree->Branch("fit_ok", &res.fit_ok);
    resultTree->Branch("percentile_1", &res.percentile_1);
    resultTree->Branch("percentile_2", &res.percentile_2);
    resultTree->Branch("percentile_3", &res.percentile_3);
    resultTree->Branch("percentile_4", &res.percentile_4);
    resultTree->Branch("percentile_5", &res.percentile_5);
    resultTree->Branch("percentile_10", &res.percentile_10);
    resultTree->Branch("percentile_15", &res.percentile_15);
    resultTree->Branch("ratio_under_0ADC", &res.ratio_under_0ADC);
    resultTree->Branch("ntrack_pass_through_channel", &res.ntrack_pass_through_channel);
    resultTree->Branch("FitClass", &res.fit_Class);
    // Refit branches
    resultTree->Branch("refit", &res.refit);
    resultTree->Branch("refit_mpv", &res.refit_mpv);
    resultTree->Branch("refit_mpv_error", &res.refit_mpv_error);
    resultTree->Branch("refit_width", &res.refit_width);
    resultTree->Branch("refit_width_error", &res.refit_width_error);
    resultTree->Branch("refit_total_area", &res.refit_total_area);
    resultTree->Branch("refit_total_area_error", &res.refit_total_area_error);
    resultTree->Branch("refit_gaus_sigma", &res.refit_gaus_sigma);
    resultTree->Branch("refit_gaus_sigma_error", &res.refit_gaus_sigma_error);
    resultTree->Branch("refit_covarience_matrix", &res.refit_covarience_matrix);
    resultTree->Branch("refit_chi2", &res.refit_chi2);
    resultTree->Branch("refit_ndf", &res.refit_ndf);
    resultTree->Branch("refit_chi2_ndf", &res.refit_chi2_ndf);
    resultTree->Branch("refit_rebin4_chi2", &res.refit_rebin4_chi2);
    resultTree->Branch("refit_rebin4_ndf", &res.refit_rebin4_ndf);
    resultTree->Branch("refit_rebin4_chi2_ndf", &res.refit_rebin4_chi2_ndf);
    resultTree->Branch("refit_fit_status", &res.refit_fit_status);
    resultTree->Branch("refit_fit_ok", &res.refit_fit_ok);
    resultTree->Branch("refit_percentile_1", &res.refit_percentile_1);
    resultTree->Branch("refit_percentile_2", &res.refit_percentile_2);
    resultTree->Branch("refit_percentile_3", &res.refit_percentile_3);
    resultTree->Branch("refit_percentile_4", &res.refit_percentile_4);
    resultTree->Branch("refit_percentile_5", &res.refit_percentile_5);
    resultTree->Branch("refit_percentile_10", &res.refit_percentile_10);
    resultTree->Branch("refit_Class", &res.refit_Class);
    resultTree->Branch("refit_ratio_under_0ADC", &res.refit_ratio_under_0ADC);
    resultTree->Branch("direct_threshold", &res.direct_threshold);
    resultTree->Branch("direct_threshold_error", &res.direct_threshold_error);
    resultTree->Branch("direct_width", &res.direct_width);
    resultTree->Branch("direct_width_error", &res.direct_width_error);
    resultTree->Branch("direct_chi2", &res.direct_chi2);
    resultTree->Branch("direct_ndf", &res.direct_ndf);
    resultTree->Branch("direct_chi2_ndf", &res.direct_chi2_ndf);
    resultTree->Branch("direct_rebin4_chi2", &res.direct_rebin4_chi2);
    resultTree->Branch("direct_rebin4_ndf", &res.direct_rebin4_ndf);
    resultTree->Branch("direct_chi2_ndf_bin4", &res.direct_chi2_ndf_bin4);
    resultTree->Branch("direct_status", &res.direct_status);
    resultTree->Branch("direct_ok", &res.direct_ok);
    resultTree->Branch("direct_parameter_at_limit", &res.direct_parameter_at_limit);
    resultTree->Branch("direct_efficiency", &res.direct_efficiency);
    resultTree->Branch("direct_efficiency_range", &res.direct_efficiency_range);
    resultTree->Branch("direct_lg_scale", &res.direct_lg_scale);
    resultTree->Branch("direct_target_area", &res.direct_target_area);
    resultTree->Branch("direct_actual_area", &res.direct_actual_area);
    resultTree->Branch("direct_x_min", &res.direct_x_min);
    resultTree->Branch("direct_x_max", &res.direct_x_max);
    resultTree->Branch("direct_nbins", &res.direct_nbins);

    int fitCount = 0;
    int fitOkCount = 0;
    int refitCount = 0;
    std::cout << "\n=== Fitting histograms (no plot) ===" << std::endl;
    
    for (TH1D* h : histograms) {
        res = FitResult{};
        MC_efficiency = -1.0;
        efficiency = -1.0;
        eff_ratio = -1.0;
        const char* histName = h->GetName();
        res.entries = h->GetEntries();

        sscanf(histName, "hMIP_%d", &res.cellid);
        res.layer = res.cellid / 100000;
        res.chip = (res.cellid / 10000) % 10;
        res.channel = res.cellid % 10000;
        const auto ntracks_it = NtrackPassThroughChannel.find(res.cellid);
        if (ntracks_it != NtrackPassThroughChannel.end()) {
            res.ntrack_pass_through_channel = ntracks_it->second;
        }
        FitOut fitOut = fitLandauGaus(h, 200, false);
        copyFitResult(res, fitOut);
        res.ratio_adc_le50 =
            res.entries > 0 ? double(res.entries_adc_le50) / res.entries : -1.0;
        if (total_area_map_loaded) {
            const auto efficiency_it = totalAreaMap.find(res.cellid);
            if (efficiency_it != totalAreaMap.end()) {
                MC_efficiency = efficiency_it->second;
            }
        }
        efficiency = (res.ntrack_pass_through_channel > 0 && res.entries > 0) ? double(res.entries) / res.ntrack_pass_through_channel : -1.0;
        eff_ratio = (MC_efficiency > 0 && efficiency > 0) ? efficiency / MC_efficiency : -1.0;
        if (shouldRefit(fitOut.fit_Class, res.cellid, eff_ratio)) {
            if (!total_area_map_loaded) {
                totalAreaMap = createTotalAreaMap(efficiency_input_file);
                total_area_map_loaded = true;
            }
            const auto efficiency_it = totalAreaMap.find(res.cellid);
            const int entries_above_threshold =
                res.entries - res.entries_adc_le50;
            double total_area_for_refit =
                entries_above_threshold > 0 && efficiency_it != totalAreaMap.end()
                    ? efficiency_it->second * entries_above_threshold / double(res.entries) * res.ntrack_pass_through_channel
                    : -1.0;
            if (eff_ratio >1.0){
                total_area_for_refit = entries_above_threshold;
            }
            const RefitResult refit =
                fitLandauGausWithThreshold(h, 200, 1, res, total_area_for_refit);
            copyRefitResult(res, refit);
            refitCount++;
        }
        
        resultTree->Fill();
        
        delete h;
        h = nullptr;
        
        fitCount++;
        if (fitOut.fit_ok) fitOkCount++;

        if (fitCount % 100 == 0) {
            std::cout << "Processed " << fitCount << " histograms, " << fitOkCount << " successful fits, " << refitCount << " refits" << std::endl;
            fout->cd();
            resultTree->Write(0, TObject::kOverwrite);
            fout->Flush();
        }
    }

    std::cout << "\nFitting complete:" << std::endl;
    std::cout << "  Total histograms: " << fitCount << std::endl;
    std::cout << "  Successful fits: " << fitOkCount << std::endl;
    std::cout << "  Refits performed: " << refitCount << std::endl;

    histograms.clear();

    fout->cd();
    resultTree->Write();
    delete resultTree;
    resultTree = nullptr;

    fout->Close();
    delete fout;
    fout = nullptr;
    
    fin->Close();
    delete fin;
    fin = nullptr;

    std::cout << "Output file created: " << output_file << std::endl;
}

void fit_histograms() {
    fit_histograms_noplot_direct_sharded(
        "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out_calibration/mip_calib/21987-22163/mip_neighborcheck_nofit.root",
        "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out_calibration/mip_calib/21987-22163/mip_neighborcheck_direct_fitted_part_0.root",
        kDefaultEfficiencyInput,
        0,
        5);
}

int main(int argc, char** argv) {
    if (argc == 1) {
        fit_histograms();
    } else if (argc == 6) {
        fit_histograms_noplot_direct_sharded(
            argv[1], argv[2], argv[3], std::stoi(argv[4]), std::stoi(argv[5]));
    } else {
        std::cout << "Usage: " << argv[0]
                  << " <input_file> <output_file> <efficiency_input_file>"
                  << " <shard_index> <shard_count>" << std::endl;
        std::cout << "Example: " << argv[0]
                  << " mip.root mip_fitted_part_0.root simulation_mip.root 0 5"
                  << std::endl;
    }
    return 0;
}
