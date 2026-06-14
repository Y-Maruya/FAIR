/**
 * Standalone script to fit histograms from ROOT file
 * Usage: root -b -q 'fit_histograms.C("input.root", "output.root")'
 * 
 * Memory Management:
 * - TF1* fLandauGaus: Explicitly deleted after use (both in draw and non-draw cases)
 * - TH1D* h_copy: Deleted after percentile calculation
 * - TH1D* histograms: Deleted in loop after processing
 * - TCanvas* cFit: Created once, cleared before each draw, deleted at the end
 * - TFile*, TTree*: Properly closed and deleted
 * - histograms vector: Cleared after processing all histograms
 * 
 * Persistence on Forced Kill:
 * - TTree is written to disk every 100 entries (every 100 histograms processed)
 * - This ensures partial results are preserved if process is terminated prematurely
 * - Each write uses TObject::kOverwrite flag to update existing entries
 */

#include "langaus.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TIterator.h>
#include <vector>
#include <unordered_map>
#include <set>
#include <fstream>
#include <sstream>
#include <string>
#include <TParameter.h>
#include <cmath>
#include <algorithm>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TFitResultPtr.h>
#include <TSystem.h>
#include <cctype>

 enum FitClass {
    FitClass_LessStatistics = -1,
    FitClass_GoodFit = 0,
    FitClass_fitFailed = 1,
    FitClass_BadChi2 = 2,
    FitClass_ParameterAtLimit = 3,
    FitClass_gaussigma_small_Threshold = 4,
    FitClass_gaussigma_small_Retry = 5
 };

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
    double percentile_1 = -1.0;
    double percentile_2 = -1.0;
    double percentile_3 = -1.0;
    double percentile_4 = -1.0;
    double percentile_5 = -1.0;
    double percentile_10 = -1.0;
    double percentile_15 = -1.0;
    int entries_adc_le50 = 0;
    TF1* fLandauGaus = nullptr;  // Add TF1 pointer to struct for delayed deletion
    int fit_Class = -1;
};

struct FitResult {
    int cellid = -1;
    int layer = -1;
    int chip = -1;
    int channel = -1;
    double x_mm = -999.0;
    double y_mm = -999.0;
    int entries = 0;
    int entries_adc_le50 = 0;
    double ratio_adc_le50 = -1.0;
    double mpv = -1.0;
    double width = -1.0;
    double total_area = -1.0;
    double gaus_sigma = -1.0;
    double max_x = -1.0;
    double FWHM = -1.0;
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
};

void calculateClassification(FitOut& fitOut) {
    if (fitOut.fit_status != 0) {
        fitOut.fit_Class = FitClass_fitFailed;
    } else if (fitOut.entries < 200) {
        fitOut.fit_Class = FitClass_LessStatistics;
    } else if (fitOut.ndf > 0 && fitOut.chi2 / double(fitOut.ndf) > 3 || fitOut.chi2 < 0) {
        fitOut.fit_Class = FitClass_BadChi2;
    } else if (fitOut.mpv <= 0 || fitOut.width <= 0 || fitOut.total_area <= 0 || fitOut.gaus_sigma <= 0 || (fitOut.gaus_sigma >149 && fitOut.width > 89)) {
        fitOut.fit_Class = FitClass_ParameterAtLimit;
    } else if (fitOut.gaus_sigma < 20 && fitOut.width <50) {
        fitOut.fit_Class = FitClass_gaussigma_small_Threshold;
    } else if (fitOut.gaus_sigma < 20 && fitOut.width >= 50) {
        fitOut.fit_Class = FitClass_gaussigma_small_Retry;
    } else {
        fitOut.fit_Class = FitClass_GoodFit;
    }
}

FitOut fitLandauGaus(TH1D* h, int minEntries = 200, bool calculate_fwhm = false, bool draw_fit = false, double min_gaus_sigma = 0.0, int fix_total_area = -1) {
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
    h_copy_0->Rebin(4);  // Rebin by factor of 4 to reduce fluctuations
    double x_max = h_copy_0->GetBinCenter(h_copy_0->GetMaximumBin());
    delete h_copy_0;  // Clean up after use

    pllo[0] = 10;    pllo[1] = 50;   pllo[2] = 100;  pllo[3] = min_gaus_sigma > 0.0 ? min_gaus_sigma : 0;  // Set minimum gaus_sigma if specified
    plhi[0] = 90;  plhi[1] = 1200; plhi[2] = 100000; plhi[3] = 150;
    sv[0] = 50; sv[1] = x_max; sv[2] = 30000; sv[3] = 60;
    if (fix_total_area > 0) {
        sv[2] = fix_total_area;
        pllo[2] = fix_total_area;
        plhi[2] = fix_total_area;
    }
    int fitResult = 0;
    TF1* fLandauGaus = langaufit(h, fr, sv, pllo, plhi, fps, fpe, &chisqr, &ndf, &fitResult);
    if (!fLandauGaus) {
        std::cout << "Landau-Gaussian fit failed for histogram: " << h->GetName() << std::endl; 
        return r;
    }
    
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
    
    // Calculate percentile_5 and entries_adc_le50
    TH1D* h_copy = (TH1D*)h->Clone("h_copy");
    h_copy->SetDirectory(nullptr);
    for (int i = 1; i <= h_copy->GetNbinsX(); ++i) {
        if (h_copy->GetBinCenter(i) <= 50) {
            r.entries_adc_le50 += h_copy->GetBinContent(i);
            h_copy->SetBinContent(i, 0); // Zero out bins <= 50 for ratio calculation
        }
    }
    int sum = 0;
    for (int i = 1; i <= h_copy->GetNbinsX(); ++i) {
        r.percentile_1 = -1.0;
        r.percentile_2 = -1.0;
        r.percentile_3 = -1.0;
        r.percentile_4 = -1.0;
        r.percentile_5 = -1.0;
        r.percentile_10 = -1.0;
        r.percentile_15 = -1.0;
        sum += h_copy->GetBinContent(i);
        if (sum >= r.entries * 0.01 && r.percentile_1 < 0) {
            r.percentile_1 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries * 0.02 && r.percentile_2 < 0) {
            r.percentile_2 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries * 0.03 && r.percentile_3 < 0) {
            r.percentile_3 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries * 0.04 && r.percentile_4 < 0) {
            r.percentile_4 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries*0.05 && r.percentile_5 < 0) {
            r.percentile_5 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries*0.10 && r.percentile_10 < 0) {
            r.percentile_10 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries*0.15 && r.percentile_15 < 0) {
            r.percentile_15 = h_copy->GetBinCenter(i);
             break; // No need to continue after 15%
        }
    }
    delete h_copy;  // Clean up after use
    h_copy = nullptr;

    double x = calculate_fwhm ? h->FindBin(maxx) : -1.0;
    double chi2_ndf = (ndf > 0) ? (chisqr / double(ndf)) : -1.0;
    if (chi2_ndf > 3 || chi2_ndf < 0) {
        r.fit_status = fitResult;
        r.fit_ok = false;
    } else {
        r.fit_status = fitResult;
        r.fit_ok = true;
    }
    calculateClassification(r);  // Determine fit classification based on results
    // Always store TF1 pointer - caller decides when to delete
    r.fLandauGaus = fLandauGaus;
    if (r.gaus_sigma < 1){
        FitOut FitO = fitLandauGaus(h, minEntries, calculate_fwhm, draw_fit, 20);  // Retry fit if gaus_sigma is unphysically small
        std::cout << "Initial fit had gaus_sigma = " << r.gaus_sigma << ", retrying with min_gaus_sigma = 20" << std::endl;
        if (FitO.fit_ok && FitO.chi2/FitO.ndf < chisqr/ndf && FitO.gaus_sigma > 20+0.1) {  // Accept retry fit if successful and has better chi2/ndf
            r = FitO;  // Use results from retry fit if successful
            std::cout << "Retry fit successful with gaus_sigma = " << FitO.gaus_sigma << std::endl;
        }
    }
    return r;
}

FitOut fitLandauGausWithThreshold(TH1D* h, int minEntries, double min_gaus_sigma, double& threshold_out, FitResult& res, TH1D*& h_ratio_out, double totalArea) {
    // Set the TotalArea is Ntracks pass through channel.
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
    TH1D* h_copy_0 = (TH1D*)h->Clone("h_copy_0");
    h_copy_0->SetDirectory(nullptr);
    h_copy_0->Rebin(4);  // Rebin by factor of 4 to reduce fluctuations
    double x_max = h_copy_0->GetBinCenter(h_copy_0->GetMaximumBin());
    delete h_copy_0;  // Clean up after use

    pllo[0] = 1;    pllo[1] = 1;   pllo[2] = 100;  pllo[3] = min_gaus_sigma > 0.0 ? min_gaus_sigma : 0;  // Set minimum gaus_sigma if specified
    plhi[0] = 90;  plhi[1] = 1200; plhi[2] = 100000; plhi[3] = 150;
    sv[0] = 50; sv[1] = x_max; sv[2] = 30000; sv[3] = 60;
    sv[2] = totalArea > 0 ? totalArea*h->GetBinWidth(1) : 30000;
    pllo[2] = totalArea > 0 ? totalArea*h->GetBinWidth(1)  : 100;
    plhi[2] = totalArea > 0 ? totalArea*h->GetBinWidth(1) : 100000;
    if (res.mpv <51){
        fr[0] = x_max > 0 ? x_max : 50;
        fr[1] = x_max > 0 ? x_max * 2.5 : 1000;
    }else{
        fr[0] = res.mpv > 0 ? res.mpv : 50;
        fr[1] = res.mpv > 0 ? res.mpv * 2.5 : 1000;
    }
    // fr[0] = res.mpv > 0 ? res.mpv : 50;  // Set fit range lower bound to MPV if available
    // fr[1] = res.mpv > 0 ? res.mpv * 2.5 : 1000;  // Set fit range upper bound to 2.5x MPV if available
    int fitResult = 0;
    TF1* fLandauGaus = langaufit(h, fr, sv, pllo, plhi, fps, fpe, &chisqr, &ndf, &fitResult);
    if (!fLandauGaus) {
        r.fit_status = -3;  // Mark as failed fit
        r.fit_ok = false;
        return r;
    }
    r.mpv = fps[1];
    r.width = fps[0];
    r.total_area = fps[2];
    r.gaus_sigma = fps[3];
    r.max_x = -1.0;
    r.FWHM = -1.0;  // FWHM calculation can be added if needed
    r.entries = h->GetEntries();
    r.chi2 = chisqr;
    // recalculate chi2 with other rebin 4.
    TH1D* h_rebin = (TH1D*)h->Rebin(4, "h_rebin");
    double sum_chi2 = 0;
    for (int i = 1; i <= h_rebin->GetNbinsX(); ++i) {
        if (h_rebin->GetBinCenter(i) > fr[1]) continue;  // Skip bins outside fit range
        if (h_rebin->GetBinCenter(i) < fr[0]) continue;  // Skip bins outside fit range
        double expected = fLandauGaus->Eval(h_rebin->GetBinCenter(i)) * h_rebin->GetBinWidth(1)/h->GetBinWidth(1);  // Scale expected value by rebinning factor
        double observed = h_rebin->GetBinContent(i);
        if (expected > 0) {
            sum_chi2 += std::pow(observed - expected, 2) / expected;
        }
    }
    // double chi2_rebin = fLandauGaus->Chi2Test(h_rebin, "WW CHI2/NDF"); // Get chi2 from rebinning to check fit quality
    r.ndf = (ndf + fLandauGaus->GetNpar())/4. - fLandauGaus->GetNpar();  // Adjust ndf for rebinning by factor of 4
    std::cout << "Refit results: MPV = " << r.mpv << ", Width = " << r.width << ", Total Area = " << r.total_area 
              << ", Gauss Sigma = " << r.gaus_sigma << ", Chi2 (rebin 4) = " << sum_chi2 << ", NDF (rebin 4) = " << r.ndf 
              << ", Fit Class = " << r.fit_Class << std::endl;
    delete h_rebin;
    r.ndf = ndf;
    r.fit_status = fitResult;
    r.fit_ok = (chisqr / double(ndf) <= 5 && chisqr >= 0) ? true : false;
    TH1D* h_ratio = (TH1D*)h->Clone("h_ratio");
    h_ratio->SetDirectory(nullptr);
    for (int i = 1; i <= h_ratio->GetNbinsX(); ++i) {
        if (fLandauGaus->Eval(h->GetBinCenter(i)) > 1 && h->GetBinCenter(i) > 50){
            h_ratio->SetBinContent(i, h->GetBinContent(i) / fLandauGaus->Eval(h->GetBinCenter(i)));
            h_ratio->SetBinError(i, h->GetBinError(i) / fLandauGaus->Eval(h->GetBinCenter(i)));
        }else if (h->GetBinContent(i) != 0){
            h_ratio->SetBinContent(i, -999);
            h_ratio->SetBinError(i, -999);
        }
    }
    // h_ratio->Rebin(4);  // Rebin ratio histogram to reduce fluctuations
    
    h_ratio_out = h_ratio;  // Return ratio histogram to caller
    
    // Calculate percentile_5 and entries_adc_le50
    TH1D* h_copy = (TH1D*)h->Clone("h_copy");
    h_copy->SetDirectory(nullptr);
    for (int i = 1; i <= h_copy->GetNbinsX(); ++i) {
        if (h_copy->GetBinCenter(i) <= 50) {
            r.entries_adc_le50 += h_copy->GetBinContent(i);
            h_copy->SetBinContent(i, 0); // Zero out bins <= 50 for ratio calculation
        }
    }
    int sum = 0;
    for (int i = 1; i <= h_copy->GetNbinsX(); ++i) {
        sum += h_copy->GetBinContent(i);
        if (sum >= r.entries * 0.01) {
            r.percentile_1 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries * 0.02) {
            r.percentile_2 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries * 0.03) {
            r.percentile_3 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries * 0.04) {
            r.percentile_4 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries*0.05) {
            r.percentile_5 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries*0.10) {
            r.percentile_10 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries*0.15) {
            r.percentile_15 = h_copy->GetBinCenter(i);
             break; // No need to continue after 15%
        }
    }
    delete h_copy;  // Clean up after use
    h_copy = nullptr;

    double chi2_ndf = (ndf > 0) ? (chisqr / double(ndf)) : -1.0;
    if (chi2_ndf > 5 || chi2_ndf < 0) {
        r.fit_status = fitResult;
        r.fit_ok = false;
    } else {
        r.fit_status = fitResult;
        r.fit_ok = true;
    }
    calculateClassification(r);  // Determine fit classification based on results
    // Always store TF1 pointer - caller decides when to delete
    r.fLandauGaus = fLandauGaus;
    
    return r;
}

double computeTotalEfficiencyBinned(TF1* fLandauGaus,
                                    int nbins,
                                    double T0,
                                    double sigma,
                                    double xmin,
                                    double xmax) {
    if (!fLandauGaus || sigma <= 0) return -1.0;

    double numer = 0.0;
    double denom = 0.0;

    for (int ibin = 1; ibin <= nbins; ++ibin) {
        double xlow = xmin + (ibin - 1) * (xmax - xmin) / nbins;
        double xhigh = xmin + ibin * (xmax - xmin) / nbins;
        double xcenter = 0.5 * (xlow + xhigh);

        if (xhigh < xmin) continue;
        if (xlow > xmax) break;

        double expected = fLandauGaus->Integral(xlow, xhigh);
        double eff = 0.5 * (1.0 + TMath::Erf(
            (xcenter - T0) / (std::sqrt(2.0) * sigma)
        ));

        denom += expected;
        numer += expected * eff;
    }

    if (denom <= 0) return -1.0;
    return numer / denom;
}


void collectHistograms(TDirectory* dir, const std::string& path, std::vector<TH1D*>& histograms, std::unordered_map<int, int>& NtrackPassThroughChannel) {
    if (!dir) return;
    
    TIter next(dir->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        std::string fullPath = path + "/" + key->GetName();
        TObject* obj = dir->Get(key->GetName());
        
        if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
            TDirectory* subdir = (TDirectory*)obj;
            collectHistograms(subdir, fullPath, histograms, NtrackPassThroughChannel);
        } else if (obj->IsA()->InheritsFrom(TH1D::Class())) {
            TH1D* h = (TH1D*)obj;
            h->Rebin(4);  // Rebin by factor of 4 to reduce fluctuations
            h->SetDirectory(nullptr);  // Prevent automatic deletion
            histograms.push_back(h);
            // std::cout << "Found histogram: " << fullPath << " (entries: " << h->GetEntries() << ")" << std::endl;
        } else if (obj->IsA()->InheritsFrom(TParameter<int>::Class())) {
            dir->cd();
            TParameter<int>* p = (TParameter<int>*)obj;
            std::string name = p->GetName();
            if (name.find("Ntracks_pass_through_channel_") == 0) {
                int cellid = std::stoi(name.substr(std::string("Ntracks_pass_through_channel_").length()));
                if (NtrackPassThroughChannel.find(cellid) != NtrackPassThroughChannel.end()) {
                    std::cerr << "Warning: duplicate Ntracks_pass_through_channel for cellid " << cellid << std::endl;
                } else {
                    NtrackPassThroughChannel[cellid] = p->GetVal();
                    // p->Print("v");
                    // std::cout << "Found Ntracks_pass_through_channel : "<<fullPath<<" for cellid " << cellid << ": " << p->GetVal() << std::endl;
                }
            }
        }
    }
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
                input.entries[cellid] = static_cast<TH1*>(obj)->GetEntries();
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


TGraphErrors* makeCumulativeDeficitGraph(TH1D* h, TF1* fFit, double x_min, double x_max, const char* name = "gDeficit") {
    if (!h || !fFit) return nullptr;

    TGraphErrors* g = new TGraphErrors();
    g->SetName(name);
    g->SetTitle("Cumulative Deficit D(T);ADC threshold candidate T;D(T) = #frac{#int_{0}^{T}(Fit - Data)dADC}{#int_{0}^{T}Fit dADC}");

    double cum_data = 0.0;
    double cum_fit = 0.0;
    int ip = 0;

    for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
        const double x = h->GetBinCenter(ibin);
        if (x <= 50) continue;  // skip pedestal/no-hold dominated region, same convention as percentile calculation
        if (x < x_min) continue;
        if (x > x_max) break;

        const double data = h->GetBinContent(ibin);
        const double fit = fFit->Integral(x - 0.5 * h->GetBinWidth(ibin), x + 0.5 * h->GetBinWidth(ibin))/h->GetBinWidth(ibin);  // Get fit value for this bin by integrating over the bin width and normalizing by bin width
        if (fit <= 0.0) continue;

        cum_data += data;
        cum_fit += fit;
        if (cum_fit <= 0.0) continue;

        const double D = (cum_fit - cum_data) / cum_fit;
        const double err = std::sqrt(std::max(cum_data, 1.0)) / cum_fit;

        g->SetPoint(ip, x, D);
        g->SetPointError(ip, 0.0, err);
        ++ip;
    }

    return g;
}

double thresholdEfficiencyErf(double x, double threshold, double width) {
    if (width <= 0.0) {
        return x >= threshold ? 1.0 : 0.0;
    }

    return 0.5 * (1.0 + TMath::Erf(
        (x - threshold) / (std::sqrt(2.0) * width)
    ));
}

TF1* makeCumulativeDeficitIntegralModel(
    TF1* fLandauGaus,
    double x0,
    double x_fit_min,
    double x_fit_max,
    const char* name = "f_cumdef_integral"
) {
    if (!fLandauGaus) return nullptr;

    /*
      Parameters:
        [0] threshold T0
        [1] threshold width sigma
        [2] high-ADC residual offset B

      Model:
        D(T) = B+(1-B)*Integral_x0^T f_LG(x) * (1 - eff(x;T0,sigma)) dx
              /Integral_x0^T f_LG(x) dx

      This is NOT a simple erf fit to cumulative deficit.
      T0 and sigma are differential hit-efficiency parameters.
    */

    auto model = [fLandauGaus, x0](double* xx, double* par) {
        const double T = xx[0];

        const double threshold = par[0];
        const double width     = par[1];
        // const double offset    = par[2];

        if (T <= x0) return 1.0;

        const double xmin = x0;
        const double xmax = T;

        TF1 f_loss(
            "f_loss_tmp",
            [fLandauGaus, threshold, width](double* x, double*) {
                const double adc = x[0];
                const double eff = thresholdEfficiencyErf(adc, threshold, width);
                return fLandauGaus->Eval(adc) * (1.0 - eff);
            },
            xmin,
            xmax,
            0
        );

        const double denom = fLandauGaus->Integral(xmin, xmax);
        if (denom <= 0.0 || !std::isfinite(denom)) return 1.0;

        const double numer = f_loss.Integral(xmin, xmax);
        double core = numer / denom;

        if (!std::isfinite(core)) core = 1.0;

        // return offset + (1.0 - offset) * core;
        return core;  // Use a fixed offset for now
    };

    TF1* f_model = new TF1(name, model, x_fit_min, x_fit_max, 3);

    f_model->SetParName(0, "threshold");
    f_model->SetParName(1, "width");
    f_model->SetParName(2, "offset");

    return f_model;
}

std::unordered_map<int, double> createTotalAreaMap(const char* input_file) {
    std::unordered_map<int, double> totalAreaMap;
    TFile* input_file_ptr = TFile::Open(input_file, "READ");
    if (!input_file_ptr || input_file_ptr->IsZombie()) {
        std::cerr << "Error: cannot open input file: " << input_file << std::endl;
        delete input_file_ptr;
        return totalAreaMap;
    }

    EfficiencyInput input;
    TDirectory* mip_dir = input_file_ptr->GetDirectory("MIP");
    collectEfficiencyInput(mip_dir ? mip_dir : input_file_ptr, input);
    std::cout << "Collected entries and ntrack_pass_through_channel values:" << std::endl;
    for (const auto& item : input.entries) {
        const int cellid = item.first;
        const Long64_t entries = item.second;
        const int ntrack_pass_through_channel =
            input.ntrack_pass_through.count(cellid) > 0
                ? input.ntrack_pass_through[cellid]
                : -1;
        totalAreaMap[cellid] = ntrack_pass_through_channel > 0
            ? static_cast<double>(entries) / ntrack_pass_through_channel
            : -1.0;
        std::cout << "CellID: " << cellid
                  << ", Entries: " << entries
                  << ", Ntracks_pass_through_channel: " << ntrack_pass_through_channel
                  << ", TotalArea: " << totalAreaMap[cellid] << std::endl;
    }

    input_file_ptr->Close();
    delete input_file_ptr;
    return totalAreaMap;
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

double computeThresholdEfficiencyFromLandauGaus(
    TF1* fLandauGaus,
    double threshold,
    double width,
    double xmin,
    double xmax
) {
    if (!fLandauGaus || width <= 0.0) return -1.0;

    TF1 f_weighted(
        "f_weighted_tmp",
        [fLandauGaus, threshold, width](double* x, double*) {
            const double adc = x[0];
            const double eff = thresholdEfficiencyErf(adc, threshold, width);
            return fLandauGaus->Eval(adc) * eff;
        },
        xmin,
        xmax,
        0
    );

    const double denom = fLandauGaus->Integral(xmin, xmax);
    const double numer = f_weighted.Integral(xmin, xmax);

    if (denom <= 0.0) return -1.0;
    return numer / denom;
}



// Return the expected number of entries in one histogram bin from the Landau-Gauss TF1.
// IMPORTANT: in this macro the langaus total-area parameter is set with
//   sv[2] = totalAreaCounts * h->GetBinWidth(1)
// Therefore TF1::Integral over a bin has the unit "counts * binWidth".
// Dividing by the bin width gives the expected bin content.
double landauGausBinIntegralContent(TF1* fLandauGaus, TH1D* h, int ibin) {
    if (!fLandauGaus || !h || ibin < 1 || ibin > h->GetNbinsX()) return 0.0;

    const double xlow = h->GetBinLowEdge(ibin);
    const double bw   = h->GetBinWidth(ibin);
    const double xhigh = xlow + bw;
    if (bw <= 0.0) return 0.0;

    return fLandauGaus->Integral(xlow, xhigh) / bw;
}

// Sum the Landau-Gauss expectation in the exact histogram-bin convention used
// by the direct threshold fit: bin integral divided by bin width.
double sumLandauGausIntegralInRange(TF1* fLandauGaus, TH1D* h, double xmin, double xmax) {
    if (!fLandauGaus || !h) return -1.0;

    double sum = 0.0;
    for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
        const double x = h->GetBinCenter(ibin);
        if (x <= xmin) continue;
        if (x > xmax) continue;
        sum += landauGausBinIntegralContent(fLandauGaus, h, ibin);
    }
    return sum;
}

double sumDataInRange(TH1D* h, double xmin, double xmax) {
    if (!h) return -1.0;

    double sum = 0.0;
    for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
        const double x = h->GetBinCenter(ibin);
        if (x <= xmin) continue;
        if (x > xmax) continue;
        sum += h->GetBinContent(ibin);
    }
    return sum;
}

TF1* makeDirectThresholdModelIntegral(
    TF1* fLandauGaus,
    TH1D* h,
    double lg_scale,
    double xmin,
    double xmax,
    const char* name = "f_direct_threshold"
) {
    if (!fLandauGaus || !h) return nullptr;

    auto model = [fLandauGaus, h, lg_scale](double* xx, double* par) {
        const double adc = xx[0];
        const double threshold = par[0];
        const double width = par[1];

        const int ibin = h->FindBin(adc);
        const double fit_bin = lg_scale * landauGausBinIntegralContent(fLandauGaus, h, ibin);
        return fit_bin * thresholdEfficiencyErf(adc, threshold, width);
    };

    TF1* f_model = new TF1(name, model, xmin, xmax, 2);
    f_model->SetParName(0, "threshold");
    f_model->SetParName(1, "width");
    return f_model;
}

TF1* makeScaledLandauGausModelIntegral(
    TF1* fLandauGaus,
    TH1D* h,
    double lg_scale,
    double xmin,
    double xmax,
    const char* name = "f_scaled_lg"
) {
    if (!fLandauGaus || !h) return nullptr;

    auto model = [fLandauGaus, h, lg_scale](double* xx, double*) {
        const int ibin = h->FindBin(xx[0]);
        return lg_scale * landauGausBinIntegralContent(fLandauGaus, h, ibin);
    };

    TF1* f_model = new TF1(name, model, xmin, xmax, 0);
    return f_model;
}

double computeDirectThresholdEfficiencyIntegral(
    TF1* fLandauGaus,
    TH1D* h,
    double lg_scale,
    double threshold,
    double width,
    double xmin,
    double xmax
) {
    if (!fLandauGaus || !h || width <= 0.0) return -1.0;

    double denom = 0.0;
    double numer = 0.0;

    for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
        const double x = h->GetBinCenter(ibin);
        if (x <= xmin) continue;
        if (x > xmax) continue;

        const double fit = lg_scale * landauGausBinIntegralContent(fLandauGaus, h, ibin);
        const double eff = thresholdEfficiencyErf(x, threshold, width);
        denom += fit;
        numer += fit * eff;
    }

    if (denom <= 0.0) return -1.0;
    return numer / denom;
}
// void collectNtrack(TDirectory* dir, std::unordered_map<int, int>& NtrackPassThroughChannel) {

//     TTree *tree = (TTree*)dir->Get("ntrack_pass_through_channel");
//     int cellid, ntracks;
//     if (tree) {
//         tree->SetBranchAddress("cellid", &cellid);
//         tree->SetBranchAddress("ntracks", &ntracks);
//         for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
//             tree->GetEntry(i);
//             NtrackPassThroughChannel[cellid] = ntracks;
//         }
//     }
// }

std::set<int> loadCellIds(const char* cellid_file) {
    std::set<int> cellids;
    std::ifstream input(cellid_file);
    if (!input) {
        std::cerr << "Error: cannot open cellid file: " << cellid_file << std::endl;
        return cellids;
    }
    std::stringstream buffer;
    buffer << input.rdbuf();
    const std::string text = buffer.str();
    std::string number;
    for (char c : text) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            number += c;
        } else if (!number.empty()) {
            cellids.insert(std::stoi(number));
            number.clear();
        }
    }
    if (!number.empty()) cellids.insert(std::stoi(number));
    return cellids;
}

void writePngGallery(const char* plot_dir, const std::set<int>& cellids) {
    std::ofstream html(std::string(plot_dir) + "/index.html");
    html << "<!doctype html><html><head><meta charset=\"utf-8\"><title>Direct integral fits</title>"
         << "<style>body{font-family:sans-serif}section{border-top:1px solid #ccc;padding:12px}"
         << "img{max-width:48%;height:auto;margin:4px}</style></head><body>"
         << "<h1>Direct integral fit plots</h1>\n";
    for (int cellid : cellids) {
        html << "<section><h2>cellid " << cellid << "</h2>"
             << "<a href=\"fit_example_cellid_" << cellid << ".png\"><img src=\"fit_example_cellid_"
             << cellid << ".png\"></a>"
             << "<a href=\"fit_retry_cellid_" << cellid << ".png\"><img src=\"fit_retry_cellid_"
             << cellid << ".png\"></a></section>\n";
    }
    html << "</body></html>\n";
}

void fit_histograms_A(const char* input_file = "mip.root",
                      const char* output_file = "mip_fitted.root",
                      const char* cellid_file = "cellids.txt",
                      const char* plot_dir = "direct_integral_plots",
                      const char* efficiency_input_file = "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out/SimMuon/mip_neighborcheck_nofit.root") {
    const std::set<int> targetCellids = loadCellIds(cellid_file);
    if (targetCellids.empty()) {
        std::cerr << "Error: no cellids loaded from " << cellid_file << std::endl;
        return;
    }
    gSystem->mkdir(plot_dir, true);
    std::cout << "Opening input file: " << input_file << std::endl;
    TFile* fin = TFile::Open(input_file, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: cannot open input file: " << input_file << std::endl;
        return;
    }
    std::vector<TH1D*> histograms;
        std::unordered_map<int, int> NtrackPassThroughChannel;  // Map to store ntrack_pass_through_channel values
    std::cout << "Collecting histograms from MIP directory..." << std::endl;
    TDirectory* mipDir = (TDirectory*)fin->Get("MIP");
    mipDir->cd();
    if (mipDir) {
        collectHistograms(mipDir, "MIP", histograms, NtrackPassThroughChannel);
    } else {
        std::cout << "Warning: MIP directory not found, searching all histograms..." << std::endl;
        collectHistograms(fin, "", histograms, NtrackPassThroughChannel);
    }

    std::cout << "Total histograms found: " << histograms.size() << std::endl;

    // Create output file
    std::cout << "Creating output file: " << output_file << std::endl;
    TFile* fout = TFile::Open(output_file, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: cannot create output file: " << output_file << std::endl;
        fin->Close();
        delete fin;
        return;
    }
    std::unordered_map<int, double> totalAreaMap = createTotalAreaMap(efficiency_input_file);
    fout->cd();

    // Create result tree
    TTree* resultTree = new TTree("mip_fit_results", "MIP fit results");
    FitResult res;
    int ntrack_pass_through_channel = -1;  // Default value if not found in map
    resultTree->Branch("cellid", &res.cellid);
    resultTree->Branch("layer", &res.layer);
    resultTree->Branch("chip", &res.chip);
    resultTree->Branch("channel", &res.channel);
    resultTree->Branch("x_mm", &res.x_mm);
    resultTree->Branch("y_mm", &res.y_mm);
    resultTree->Branch("entries", &res.entries);
    resultTree->Branch("entries_adc_le50", &res.entries_adc_le50);
    resultTree->Branch("ratio_adc_le50", &res.ratio_adc_le50);
    resultTree->Branch("MPV", &res.mpv);
    resultTree->Branch("width", &res.width);
    resultTree->Branch("TotalArea", &res.total_area);
    resultTree->Branch("gaus_sigma", &res.gaus_sigma);
    resultTree->Branch("max_x", &res.max_x);
    resultTree->Branch("FWHM", &res.FWHM);
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
    resultTree->Branch("ntrack_pass_through_channel", &ntrack_pass_through_channel);
    resultTree->Branch("FitClass", &res.fit_Class);
    FitResult res_refit;
    TTree* refit_resultTree = new TTree("mip_refit_results", "MIP refit results");
    // Fit histograms
    int fitCount = 0;
    int fitOkCount = 0;
    std::cout << "\n=== Fitting histograms ===" << std::endl;
    
    // Create canvas for drawing (with memory management)
    TCanvas* cFit = new TCanvas("cFit", "MIP Fit Example", 800, 600);
    // cFit->SetBatch(kTRUE);  // Batch mode to avoid displaying windows
    random_shuffle(histograms.begin(), histograms.end());  // Shuffle to get a random sample for testing

    for (TH1D* h : histograms) {
        const char* histName = h->GetName();
        res.entries = h->GetEntries();

        // Extract cellid from histogram name: hMIP_CELLID
        sscanf(histName, "hMIP_%d", &res.cellid);
        res.layer = res.cellid / 100000;
        res.chip = (res.cellid / 10000) % 10;
        res.channel = res.cellid % 10000;
        int ntrack_pass_through_channel = -1;  // Reset for each histogram
        if (NtrackPassThroughChannel.find(res.cellid) != NtrackPassThroughChannel.end()) {
            ntrack_pass_through_channel = NtrackPassThroughChannel[res.cellid];
        }
        if (!targetCellids.count(res.cellid)) continue;
#if 0
        // int cellids[] = {140009,150035,150021,150013,150020,150033,150015,150030,150000,150028,150008,150002,150001,150029,150009,150026,150018,150022,150034,150025,150032,150012,150016,150003,150005,150024,150019,150027,150017,150007,150031,150006,170023,170019,170024,170030,180016,180032,180026,180024,180003,180025,180004,180000,180033,180017,180035,180007,180002,180015,180034,300010,300016,300015,300022,300034,300021,300032,300035,300025,300024,300014,300029,300007,300027,300009,300033,300012,300031,300023,300017,300004,300001,300018,300000,300003,300006,300026,300005,300008,300028,300013,300002,300011,300020,310029,310033,310011,310005,310035,310015,310007,310018,310021,310032,310025,310022,310006,310024,310020,310031,310030,310034,310009,310008,310019,310028,310004,310016,310001,310013,310017,310027,310012,310000,310023,310014,310003,310002,320031,320032,320029,320020,320017,320034,320016,320013,320033,320035,320010,320008,320007,320018,320014,320009,320025,320030,320019,320021,320003,320001,320005,320015,320022,320026,320006,320027,320024,320028,320023,320011,320004,320002,320012,330029,330031,330030,330008,330011,330023,330034,330010,330032,330012,330009,330025,330026,330005,330001,330033,330014,330002,330024,330035,330027,330004,330013,330022,330021,330007,330000,330006,330028,330015,340034,340014,340026,340024,340025,340028,340007,340012,340019,340022,340023,340027,340032,340010,340031,340035,340008,340029,340016,340000,340018,340015,340021,340003,340001,340017,340030,340011,340004,340002,340005,340033,340009,340006,350009,350021,350020,350016,350019,350008,350013,350035,350022,350017,350023,350011,350030,350033,350012,350029,350015,350028,350002,350026,350010,350003,350004,350005,350024,350027,350007,350025,350031,350000,350034,350001,350006,360004,360033,360023,360025,360035,360031,360014,360013,360011,360018,360030,360019,360024,360034,360029,360010,360001,360026,360022,360006,360003,360002,360028,360005,360009,360032,360027,360016,360021,360008,360017,360000,370034,370017,370008,370025,370013,370023,370029,370011,370016,370007,370010,370031,370019,370027,370028,370035,370006,370005,370004,370003,370022,370024,370001,370015,370020,370026,370030,370012,370002,370009,370014,370000,370032,370033,380011,380020,380007,380028,380029,380031,380016,380015,380009,380014,380032,380019,380022,380026,380030,380013,380021,380012,380003,380005,380025,380002,380008,380001,380006,380010,380004,380000,380034,380023,380018,380033,380017,380035,380027,400024,400007,400035,400033,400012,400010,400031,400016,400017,400023,400019,400014,400032,400015,400001,400030,400018,400022,400000,400025,400008,400006,400003,400005,400004,400034,400028,400013,400009,400021,400011,400002,400026,400027,400020,400029,410012,410024,410020,410029,410013,410021,410017,410030,410019,410031,410025,410009,410008,410028,410034,410035,410004,410016,410033,410015,410032,410005,410027,410000,410010,410001,410002,410026,410011,410014,410022,410006,410007,410018,410023,410003,420031,420023,420033,420017,420020,420019,420014,420012,420009,420018,420027,420025,420003,420032,420030,420029,420034,420021,420035,420028,420004,420005,420022,420000,420015,420016,420010,420024,420001,420026,420006,420011,420002,420013,420007,420008,430012,430016,430023,430019,430008,430026,430014,430032,430002,430001,430024,430021,430027,430029,430035,430009,430003,430031,430034,430004,430013,430020,430033,430022,430018,430011,430007,430010,430030,430025,430000,430015,430005,430006,430028,430017,440023,440034,440012,440009,440010,440014,440035,440026,440027,440029,440022,440018,440020,440028,440002,440008,440021,440015,440025,440032,440003,440033,440017,440016,440011,440030,440005,440031,440000,440001,440024,440013,440007,440019,440004,440006,450035,450014,450020,450032,450029,450033,450015,450030,450028,450002,450001,450023,450009,450026,450013,450018,450022,450034,450012,450011,450003,450016,450010,450004,450005,450019,450027,450024,450017,450007,450025,450031,450000,450021,450008,450006,460015,460035,460013,460023,460012,460024,460034,460010,460030,460026,460002,460022,460006,460031,460029,460003,460025,460028,460004,460014,460020,460001,460005,460019,460009,460016,460018,460032,460027,460033,460000,460017,460021,460011,460008,460007,470033,470023,470035,470019,470029,470024,470005,470008,470012,470021,470011,470006,470017,470022,470018,470009,470003,470004,470015,470028,470020,470032,470026,470030,470001,470002,470014,470010,470025,470034,470027,470016,470000,470031,470013,470007,480016,480029,480032,480018,480014,480026,480013,480031,480017,480021,480003,480019,480005,480010,480008,480025,480024,480012,480028,480011,480002,480006,480001,480000,480033,480034,480022,480035,480007,480030,480004,480015,480009,480023,480027,480020,500007,500008,500035,500033,500012,500010,500031,500005,500001,500017,500023,500019,500000,500030,500014,500015,500003,500032,500018,500022,500011,500025,500006,500002,500004,500016,500034,500013,500009,500021,500024,500027,500020,500029,510012,510024,510020,510017,510008,510032,510019,510011,510025,510004,510028,510034,510031,510016,510027,510033,510013,510015,510005,510000,510010,510001,510026,510023,510022,510006,510007,510002,510018,510003,520017,520020,520002,520009,520014,520012,520018,520027,520025,520032,520019,520029,520034,520021,520003,520028,520035,520031,520004,520005,520000,520026,520030,520015,520016,520010,520024,520022,520008,520006,520011,520001,520013,520023,520007,520033,530019,530016,530032,530001,530023,530012,530029,530005,530035,530021,530020,530011,530015,530006,540035,540027,540005,540015,540013,540017,540007,540019,550021,550020,550029,550033,550015,550002,550023,550009,550013,550034,550011,550003,550001,550005,550024,550027,550017,550031,550000,550006,550008,560035,560015,560007,560031,560029,560025,560001,560019,560032,570019,570022,570011,570006,570017,570015,570025,570007,580016,580007,580014,580031,580021,580012,580025,580024,580028,580011,580002,580005,580034,580003,580023,580020,800012,800001,800004,800027,800020,810005,810029,810013,810020,810021,810016,810027,810033,810026,810002,810003,820027,820003,820005,820029,820033,830031,840020,840011,850020,860034,860028,860024,860023,870034,870002,870024,870006,870028,870033,880024,1500022,1500007,1500012,1500013,1500016,1500001,1500011,1500003,1500032,1500000,1500006,1500004,1500015,1500034,1500009,1500031,1500033,1500002,1500029,1500019,1500020,1500017,1520034,1520000,1520029,1520001,1520002,1530020,1530015,1540016,1550002,1550000,1550008,1560035,1560016,1560034,1560032,1560031,1560028,1560022,1560006,1560019,1560023,1560020,1560018,1570034,1570024,1570033,1570035,1570021,1570020,1570023,1570029,1570019,1580032,1580034,1580020,1700016,1710026,1720012,1750029,1750033,1770019,1780024,1960031,2910032,2940015,3500012,3550003,3550033,3800012,3880006,3880029};
        int cellids[] =
{
130022,
180035,
180026,
3210011,
3210013,
3220002
};
        int numCellids = sizeof(cellids) / sizeof(cellids[0]);
        bool foundCellid = false;
        for (int i = 0; i < numCellids; ++i) {
            if (res.cellid == cellids[i]) {
                std::cout << "Processing histogram: " << histName << " (entries: " << res.entries << ")" << std::endl;
                foundCellid = true;
                break;
            }
        }
        if (!foundCellid) {
            // std::cout << "Skipping histogram: " << histName << " (cellid not in target list)" << std::endl;
            continue;  // Skip histograms that are not in the target list
        }
#endif
    
        FitOut fitOut = fitLandauGaus(h, 300, false, false);
        std::cout << "Fit results for histogram " << histName << ": MPV = " << fitOut.mpv << ", Width = " << fitOut.width 
                  << ", Gaus Sigma = " << fitOut.gaus_sigma << ", Chi2/NDF = " << fitOut.chi2 << "/" << fitOut.ndf 
                  << ", Fit OK = " << (fitOut.fit_ok ? "Yes" : "No") << ", Fit Class = " 
                  << ((fitOut.fit_Class == FitClass_LessStatistics) ? "LessStatistics" :
                      (fitOut.fit_Class == FitClass_GoodFit) ? "GoodFit" :
                      (fitOut.fit_Class == FitClass_fitFailed) ? "fitFailed" :
                      (fitOut.fit_Class == FitClass_BadChi2) ? "BadChi2" :
                      (fitOut.fit_Class == FitClass_ParameterAtLimit) ? "ParameterAtLimit" :
                      (fitOut.fit_Class == FitClass_gaussigma_small_Threshold) ? "gaussigma_small_Threshold" :
                      (fitOut.fit_Class == FitClass_gaussigma_small_Retry) ? "gaussigma_small_Retry" : "Unknown") 
                  << std::endl;
        
        // Draw on canvas
        if (!cFit) {
            cFit = new TCanvas("cFit", "MIP Fit Example", 800, 600);
        }
        cFit->cd();
        cFit->Clear();  // Clear canvas before drawing
        std::string className = (fitOut.fit_Class == FitClass_LessStatistics) ? "LessStatistics" :
                                (fitOut.fit_Class == FitClass_GoodFit) ? "GoodFit" :
                                (fitOut.fit_Class == FitClass_fitFailed) ? "fitFailed" :
                                (fitOut.fit_Class == FitClass_BadChi2) ? "BadChi2" :
                                (fitOut.fit_Class == FitClass_ParameterAtLimit) ? "ParameterAtLimit" :
                                (fitOut.fit_Class == FitClass_gaussigma_small_Threshold) ? "gaussigma_small_Threshold" :
                                (fitOut.fit_Class == FitClass_gaussigma_small_Retry) ? "gaussigma_small_Retry" : "Unknown";
        std::cout << "Drawing histogram " << histName << " with fit class: " << className << std::endl;
        h->SetTitle(Form("Class = %s, Entries: %d", className.c_str(), res.entries));
        h->Draw();
        std::cout << "drawed histogram " << histName << std::endl;
        h->GetXaxis()->SetRangeUser(-96, 1000);
        std::cout << "Set x-axis range for histogram " << histName << std::endl;
        // if (fitOut.fit_ok) continue;  // Skip drawing fit if successful
        if (fitOut.fLandauGaus && h) {
            std::cout << "Drawing fit function for " << histName << std::endl;
            try {
                fitOut.fLandauGaus->SetLineColor(kRed);
                fitOut.fLandauGaus->SetLineWidth(2);
                fitOut.fLandauGaus->Draw("same");
            } catch (const std::exception& e) {
                std::cerr << "Error drawing fit function for histogram " << histName << ": " << e.what() << std::endl;
            }
            std::cout << "Drawed fit function for " << histName << std::endl;
            double max_val = h->GetMaximum();
            if (max_val > 0 && fitOut.mpv > 0) {
                TLine* line_tmp = new TLine(fitOut.mpv, 0, fitOut.mpv, max_val);
                line_tmp->SetLineColor(kBlue);
                line_tmp->SetLineStyle(kDashed);
                line_tmp->Draw("same");
            }
            std::cout << "Finished drawing fit function for " << histName << std::endl;
        }else {
            TLatex latex_title;
            latex_title.SetNDC();
            latex_title.SetTextSize(0.03);
            latex_title.DrawLatex(0.55, 0.85, "Fit failed");
            std::cout << "Fit failed for histogram: " << histName << std::endl;
        }
        
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.03);
        
        // Check for valid values before drawing
        double mpv_val = (fitOut.mpv > -999 && fitOut.mpv < 10000) ? fitOut.mpv : 0.0;
        double width_val = (fitOut.width > 0 && fitOut.width < 10000) ? fitOut.width : 0.0;
        double gaus_sigma_val = (fitOut.gaus_sigma > -999 && fitOut.gaus_sigma < 10000) ? fitOut.gaus_sigma : 0.0;
        double chi2_val = (fitOut.chi2 > -999 && fitOut.chi2 < 1e6) ? fitOut.chi2 : 0.0;
        double percentile_val = (fitOut.percentile_5 > -999 && fitOut.percentile_5 < 10000) ? fitOut.percentile_5 : 0.0;
        
        latex.DrawLatex(0.55, 0.85, Form("MPV = %.1f", mpv_val));
        latex.DrawLatex(0.55, 0.80, Form("Width = %.1f", width_val));
        latex.DrawLatex(0.55, 0.75, Form("Gaus Sigma = %.1f", gaus_sigma_val));
        latex.DrawLatex(0.55, 0.70, Form("Chi2/NDF = %.1f/%d", chi2_val, fitOut.ndf));
        latex.DrawLatex(0.55, 0.65, Form("5 Percentile = %.1f", percentile_val));
        latex.DrawLatex(0.55, 0.60, Form("fit status = %d", fitOut.fit_status));
        latex.DrawLatex(0.55, 0.55, Form("Ntracks pass through channel = %d (%.1f \%%)", ntrack_pass_through_channel, ntrack_pass_through_channel > 0 ? double(fitOut.entries)/double(ntrack_pass_through_channel)*100.0 : -1.0));
        std::cout << "Drawing histogram " << histName << " with fit class: " << className << std::endl;
        cFit->Update();
        cFit->SaveAs(Form("%s/fit_example_cellid_%d.png", plot_dir, res.cellid));
        // Delete TF1 after SaveAs
        if (fitOut.fLandauGaus) {
            delete fitOut.fLandauGaus;
            fitOut.fLandauGaus = nullptr;
        }
        
        res.mpv = fitOut.mpv;
        res.width = fitOut.width;
        res.total_area = fitOut.total_area;
        res.gaus_sigma = fitOut.gaus_sigma;
        res.max_x = fitOut.max_x;
        res.FWHM = fitOut.FWHM;
        res.chi2 = fitOut.chi2;
        res.ndf = fitOut.ndf;
        res.fit_status = fitOut.fit_status;
        res.fit_ok = fitOut.fit_ok ? 1 : 0;
        res.entries_adc_le50 = fitOut.entries_adc_le50;
        if (res.entries > 0) {
            res.ratio_adc_le50 = double(res.entries_adc_le50) / double(res.entries);
        } else {
            res.ratio_adc_le50 = -1.0;
        }
        res.percentile_1 = fitOut.percentile_1;
        res.percentile_2 = fitOut.percentile_2;
        res.percentile_3 = fitOut.percentile_3;
        res.percentile_4 = fitOut.percentile_4;
        res.percentile_5 = fitOut.percentile_5;
        res.percentile_10 = fitOut.percentile_10;
        res.percentile_15 = fitOut.percentile_15;
        res.fit_Class = fitOut.fit_Class;
        res.ntrack_pass_through_channel = ntrack_pass_through_channel;
        double old_mpv = res.mpv;
        // if (fitOut.fit_Class == FitClass_gaussigma_small_Threshold || fitOut.fit_Class == FitClass_gaussigma_small_Retry) {
        // if (fitOut.fit_Class == 0 || fitOut.fit_Class == FitClass_gaussigma_small_Threshold || fitOut.fit_Class == FitClass_gaussigma_small_Retry) {
        // if (fitOut.fit_Class == FitClass_gaussigma_small_Threshold || fitOut.fit_Class == FitClass_gaussigma_small_Retry || isRefitAHCALCell(res.cellid)) {
        if (true){ 
        // std::cout << "Histogram " << histName << " has small gaus_sigma = " << fitOut.gaus_sigma << ", retrying fit with min_gaus_sigma = 20" << std::endl;
            TH1D* h_ratio = nullptr;
            double threshold_out = -1.0;
            const double totalAreaForRefit =
                (res.entries - res.entries_adc_le50) > 0 &&
                totalAreaMap.find(res.cellid) != totalAreaMap.end()
                    ? totalAreaMap[res.cellid] *
                        (res.entries - res.entries_adc_le50)/double(res.entries) * res.ntrack_pass_through_channel 
                    : -1.0;
            if (totalAreaForRefit <= 0) {
                std::cout << "Total area for refit is invalid for histogram " << histName << ", skipping refit." << std::endl;
                throw std::runtime_error("Invalid total area for refit");
            }
            fitOut = fitLandauGausWithThreshold(h, 300, 1, threshold_out, res, h_ratio, totalAreaForRefit);
            TCanvas* cRetry = new TCanvas(Form("cRetry_%d", res.cellid), "MIP Fit Retry", 800, 1200);
            cRetry->Divide(1,2);
            cRetry->cd(1);
            h->SetTitle(Form("Retry Fit - Class = %s, Entries: %d", 
                (fitOut.fit_Class == FitClass_LessStatistics) ? "LessStatistics" :
                (fitOut.fit_Class == FitClass_GoodFit) ? "GoodFit" :
                (fitOut.fit_Class == FitClass_fitFailed) ? "fitFailed" :
                (fitOut.fit_Class == FitClass_BadChi2) ? "BadChi2" :
                (fitOut.fit_Class == FitClass_ParameterAtLimit) ? "ParameterAtLimit" :
                (fitOut.fit_Class == FitClass_gaussigma_small_Threshold) ? "gaussigma_small_Threshold" :
                (fitOut.fit_Class == FitClass_gaussigma_small_Retry) ? "gaussigma_small_Retry" : "Unknown", 
                int(h->GetEntries())));
            h->Draw();
            h->GetXaxis()->SetRangeUser(-96, 1000);
            if (fitOut.fLandauGaus) {
                fitOut.fLandauGaus->SetLineColor(kRed);
                fitOut.fLandauGaus->SetLineWidth(2);
                fitOut.fLandauGaus->Draw("same");
                TLine* line = new TLine(res.mpv, 0, res.mpv, h->GetMaximum());
                line->SetLineColor(kBlue);
                line->SetLineStyle(kDashed);
                line->SetLineWidth(2);
                line->Draw("same");
                TLine* line2 = new TLine(std::min(res.mpv*4, 1000.0), 0, std::min(res.mpv*4, 1000.0), h->GetMaximum());
                line2->SetLineColor(kBlue);
                line2->SetLineStyle(kDashed);
                line2->SetLineWidth(2);
                line2->Draw("same");
                
            }else {
                TLatex latex;
                latex.SetNDC();
                latex.SetTextSize(0.03);
                latex.DrawLatex(0.55, 0.85, "Fit failed");
                std::cout << "Retry fit failed for histogram: " << histName << std::endl;
            }
            TLatex latexRetry;
            latexRetry.SetNDC();
            latexRetry.SetTextSize(0.03);
            latexRetry.DrawLatex(0.55, 0.85, Form("MPV = %.1f", fitOut.mpv));
            latexRetry.DrawLatex(0.55, 0.80, Form("Width = %.1f", fitOut.width));
            latexRetry.DrawLatex(0.55, 0.75, Form("Gaus Sigma = %.1f", fitOut.gaus_sigma));
            latexRetry.DrawLatex(0.55, 0.70, Form("Chi2/NDF = %.1f/%d", fitOut.chi2, fitOut.ndf));
            latexRetry.DrawLatex(0.55, 0.65, Form("5 Percentile = %.1f", fitOut.percentile_5));
            latexRetry.DrawLatex(0.55, 0.60, Form("fit status = %d", fitOut.fit_status));
            latexRetry.DrawLatex(0.55, 0.55, Form("Ntracks pass through channel = %d (%.1f \%%)", ntrack_pass_through_channel, ntrack_pass_through_channel > 0 ? double(fitOut.entries)/double(ntrack_pass_through_channel)*100.0 : -1.0));
            latexRetry.DrawLatex(0.55, 0.50, Form("MC efficiency = %.1f \%%", totalAreaMap[res.cellid]*100.0));
            latexRetry.DrawLatex(0.55, 0.40, Form("Total Area (Ntracks) = %.1f", fitOut.total_area));
             if (threshold_out > 0) {
                latexRetry.DrawLatex(0.55, 0.35, Form("Threshold = %.1f", threshold_out));
            }
            cRetry->cd(2);
            if (fitOut.fLandauGaus) {
                const double x_fit_min = std::max(50.0, res.mpv*0.2);
                const double x_fit_max = old_mpv > 0 ? std::min(old_mpv*4, 2000.0) : 2000.0;

                // Direct fit replacing the cumulative-deficit fit.
                // The second Landau-Gauss refit above is intentionally unchanged.
                // Here we fit the observed MIP spectrum directly with:
                //   Data(ADC) = scale * LandauGaus(ADC) * thresholdEfficiencyErf(ADC; T0, width)
                // where scale is chosen so that the Landau-Gauss bin-sum reproduces
                // the target area A_def used for efficiency reproduction.
                const double x_direct_min = 50.0;
                const double x_direct_max = old_mpv > 0 ? std::min(old_mpv*4, 2000.0) : 2000.0;

                double N_all = h->GetEntries();
                double N_over50 = sumDataInRange(h, x_direct_min, x_direct_max);
                double mc_efficiency = totalAreaMap[res.cellid];
                double E = ntrack_pass_through_channel * mc_efficiency;
                double A_def = E * double(N_over50) / double(N_all);

                double A_actual = sumLandauGausIntegralInRange(
                    fitOut.fLandauGaus, h, x_direct_min, x_direct_max);

                if (A_def <= 0.0 || A_actual <= 0.0) {
                    TLatex latexDirect;
                    latexDirect.SetNDC();
                    latexDirect.SetTextSize(0.03);
                    latexDirect.DrawLatex(0.15, 0.70, "Direct threshold fit skipped: invalid normalization");
                    std::cout << "Direct threshold fit skipped: invalid normalization for "
                              << histName << std::endl;
                } else {
                    const double lg_scale = A_def / A_actual;

                    TF1* f_direct = makeDirectThresholdModelIntegral(
                        fitOut.fLandauGaus,
                        h,
                        1,
                        x_direct_min,
                        x_direct_max,
                        Form("f_direct_threshold_%d", res.cellid)
                    );

                    TF1* f_scaled_lg = makeScaledLandauGausModelIntegral(
                        fitOut.fLandauGaus,
                        h,
                        lg_scale,
                        x_direct_min,
                        x_direct_max,
                        Form("f_scaled_lg_%d", res.cellid)
                    );

                    if (!f_direct || !f_scaled_lg) {
                        TLatex latexDirect;
                        latexDirect.SetNDC();
                        latexDirect.SetTextSize(0.03);
                        latexDirect.DrawLatex(0.15, 0.70, "Direct threshold model creation failed");
                        std::cout << "Direct threshold model creation failed for "
                                  << histName << std::endl;
                    } else {
                        double threshold_init = res.percentile_1 > 0 ? res.percentile_1 : 0.6 * fitOut.mpv;
                        threshold_init = std::max(x_direct_min, std::min(threshold_init, x_direct_max));

                        f_direct->SetParameters(threshold_init, 20.0);
                        f_direct->SetParLimits(0, x_direct_min, x_direct_max);
                        f_direct->SetParLimits(1, 3.0, 100.0);
                        f_direct->SetLineColor(kBlue);
                        f_direct->SetLineWidth(2);

                        f_scaled_lg->SetLineColor(kRed);
                        f_scaled_lg->SetLineStyle(kDashed);
                        f_scaled_lg->SetLineWidth(2);

                        // Fit the data spectrum directly. Use the same x range as the model.
                        TFitResultPtr direct_fit_result = h->Fit(f_direct, "RQS");
                        threshold_out = f_direct->GetParameter(0);
                        const double threshold_width = f_direct->GetParameter(1);

                        h->Draw("E");
                        h->GetXaxis()->SetRangeUser(0, x_direct_max);
                        h->SetTitle("Direct threshold fit;ADC;Entries / bin");
                        // f_scaled_lg->Draw("same");
                        f_direct->Draw("same");

                        const double eff_threshold = computeDirectThresholdEfficiencyIntegral(
                            fitOut.fLandauGaus,
                            h,
                            lg_scale,
                            threshold_out,
                            threshold_width,
                            x_direct_min,
                            x_direct_max
                        );

                        double totalArea = fitOut.fLandauGaus->Integral(-96,4000);
                        double totalAreaWithEff = f_direct->Integral(-96,4000);
                        
                        double eff_threshold_integral = totalArea > 0.0 ? totalAreaWithEff / totalArea : -1.0;

                        // These are the quantities relevant for efficiency reproduction.
                        // eff_data_MC and eff_data_Fit_def are algebraically equal with the A_def definition.
                        const double eff_data_MC = N_all / E;
                        const double eff_data_Fit_def = N_over50 / A_def;

                        double fit_sum_scaled = 0.0;
                        double fit_eff_sum_scaled = 0.0;
                        double data_sum_direct = 0.0;
                        for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
                            const double x = h->GetBinCenter(ibin);
                            if (x <= x_direct_min) continue;
                            if (x > x_direct_max) continue;

                            const double data = h->GetBinContent(ibin);
                            const double fit = lg_scale * landauGausBinIntegralContent(fitOut.fLandauGaus, h, ibin);
                            const double eff_thr = thresholdEfficiencyErf(x, threshold_out, threshold_width);

                            data_sum_direct += data;
                            fit_sum_scaled += fit;
                            fit_eff_sum_scaled += fit * eff_thr;
                        }

                        const double eff_data_actual = fit_sum_scaled > 0.0 ? data_sum_direct / fit_sum_scaled : -1.0;
                        const double eff_threshold_actual = fit_sum_scaled > 0.0 ? fit_eff_sum_scaled / fit_sum_scaled : -1.0;

                        std::cout << "Direct threshold fit results" << std::endl;
                        std::cout << "A_def                 = " << A_def << std::endl;
                        std::cout << "A_actual_sum_integral     = " << A_actual << std::endl;
                        std::cout << "LG scale              = " << lg_scale << std::endl;
                        std::cout << "N_all                 = " << N_all << std::endl;
                        std::cout << "N_over50_in_range     = " << N_over50 << std::endl;
                        std::cout << "E=Ntrack*MCeff        = " << E << std::endl;
                        std::cout << "eff_data_MC           = " << eff_data_MC << std::endl;
                        std::cout << "eff_data_Fit_def      = " << eff_data_Fit_def << std::endl;
                        std::cout << "eff_data_actual_range = " << eff_data_actual << std::endl;
                        std::cout << "eff_threshold         = " << eff_threshold << std::endl;
                        std::cout << "eff_threshold_integral = " << eff_threshold_integral << std::endl;
                        std::cout << "eff_threshold_actual  = " << eff_threshold_actual << std::endl;
                        std::cout << "direct_fit_status     = " << int(direct_fit_result) << std::endl;
                        std::cout << "direct_chi2_ndf       = "
                                  << f_direct->GetChisquare() << "/" << f_direct->GetNDF()
                                  << std::endl;

                        TLatex latexDirect;
                        latexDirect.SetNDC();
                        latexDirect.SetTextSize(0.03);
                        latexDirect.DrawLatex(0.55, 0.85, Form("Threshold = %.1f", threshold_out));
                        latexDirect.DrawLatex(0.55, 0.80, Form("Width = %.1f", threshold_width));
                        latexDirect.DrawLatex(0.55, 0.75, Form("Eff data/MC = %.1f %%", eff_data_MC*100.0));
                        latexDirect.DrawLatex(0.55, 0.70, Form("Eff threshold = %.1f %%", eff_threshold*100.0));
                        latexDirect.DrawLatex(0.55, 0.65, Form("Eff threshold integral = %.1f %%", eff_threshold_integral*100.0));
                        latexDirect.DrawLatex(0.55, 0.60, Form("Eff data/fit range = %.1f %%", eff_data_actual*100.0));
                        latexDirect.DrawLatex(0.55, 0.55, Form("Scale = %.3f", lg_scale));
                        latexDirect.DrawLatex(0.55, 0.50, Form("Chi2/NDF = %.1f/%d", f_direct->GetChisquare(), f_direct->GetNDF()));
                    }
                }                
            }
            cRetry->Update();
            cRetry->SaveAs(Form("%s/fit_retry_cellid_%d.png", plot_dir, res.cellid));
            delete cRetry;
            cRetry = nullptr;
            if (h_ratio) {
                delete h_ratio;
                h_ratio = nullptr;
            }
            cFit->cd();
        }
        resultTree->Fill();
        
        // Clean up histogram after processing
        delete h;
        h = nullptr;
        
        fitCount++;
        if (fitOut.fit_ok) fitOkCount++;

        if (fitCount % 100 == 0) {
            std::cout << "Processed " << fitCount << " histograms, " << fitOkCount << " successful fits" << std::endl;
            // Write partial results to disk for persistence on forced kill
            fout->cd();
            resultTree->Write(0, TObject::kOverwrite);
            fout->Flush();
        }
    }

    std::cout << "\nFitting complete:" << std::endl;
    std::cout << "  Total histograms: " << fitCount << std::endl;
    std::cout << "  Successful fits: " << fitOkCount << std::endl;

    // Clean up canvas
    delete cFit;
    cFit = nullptr;
    
    // Clear histogram vector
    histograms.clear();

    // Write results
    fout->cd();
    resultTree->Write();
    delete resultTree;
    resultTree = nullptr;

    // Copy original histograms to output file (if needed)
    std::cout << "\nCopying directory structure to output file..." << std::endl;
    if (mipDir) {
        TIter next(fin->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)next())) {
            TObject* obj = fin->Get(key->GetName());
            if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
                fout->mkdir(key->GetName());
            }
        }
    }

    fout->Close();
    delete fout;
    fout = nullptr;
    
    fin->Close();
    delete fin;
    fin = nullptr;

    std::cout << "\nOutput file created: " << output_file << std::endl;
    writePngGallery(plot_dir, targetCellids);
    std::cout << "PNG gallery created: " << plot_dir << "/index.html" << std::endl;
}

// Wrapper for command line usage
void fit_histograms_directfit_integral(
    const char* input_file = "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out_calibration/mip_calib/21987-22163/mip_neighborcheck_nofit.root",
    const char* output_file = "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out_calibration/mip_calib/21987-22163/mip_neighborcheck_fitted_tmps.root",
    const char* cellid_file = "cellids.txt",
    const char* plot_dir = "direct_integral_plots",
    const char* efficiency_input_file = "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out/SimMuon/mip_neighborcheck_nofit.root") {
    fit_histograms_A(input_file, output_file, cellid_file, plot_dir, efficiency_input_file);
}

// Version without drawing (faster, lower memory footprint)
void fit_histograms_noplot(const char* input_file = "mip.root", const char* output_file = "mip_fitted.root") {
    std::cout << "Opening input file (no plot mode): " << input_file << std::endl;
    TFile* fin = TFile::Open(input_file, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: cannot open input file: " << input_file << std::endl;
        return;
    }

    std::vector<TH1D*> histograms;
    std::unordered_map<int, int> NtrackPassThroughChannel;  // Map to store ntrack_pass_through_channel values
    std::cout << "Collecting histograms..." << std::endl;
    TDirectory* mipDir = (TDirectory*)fin->Get("MIP");
    if (mipDir) {
        collectHistograms(mipDir, "MIP", histograms, NtrackPassThroughChannel);
    } else {
        collectHistograms(fin, "", histograms, NtrackPassThroughChannel);
    }
    std::cout << "Total histograms found: " << histograms.size() << std::endl;

    TFile* fout = TFile::Open(output_file, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: cannot create output file: " << output_file << std::endl;
        fin->Close();
        delete fin;
        return;
    }

    TTree* resultTree = new TTree("mip_fit_results", "MIP fit results");
    FitResult res;
    int ntrack_pass_through_channel = -1;  // Default value if not found in map
    resultTree->Branch("cellid", &res.cellid);
    resultTree->Branch("layer", &res.layer);
    resultTree->Branch("chip", &res.chip);
    resultTree->Branch("channel", &res.channel);
    resultTree->Branch("entries", &res.entries);
    resultTree->Branch("entries_adc_le50", &res.entries_adc_le50);
    resultTree->Branch("ratio_adc_le50", &res.ratio_adc_le50);
    resultTree->Branch("MPV", &res.mpv);
    resultTree->Branch("width", &res.width);
    resultTree->Branch("TotalArea", &res.total_area);
    resultTree->Branch("gaus_sigma", &res.gaus_sigma);
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
    resultTree->Branch("ntrack_pass_through_channel", &ntrack_pass_through_channel);
    resultTree->Branch("FitClass", &res.fit_Class);
    int fitCount = 0;
    int fitOkCount = 0;
    std::cout << "\n=== Fitting histograms (no plot) ===" << std::endl;
    
    for (TH1D* h : histograms) {
        const char* histName = h->GetName();
        res.entries = h->GetEntries();

        sscanf(histName, "hMIP_%d", &res.cellid);
        res.layer = res.cellid / 100000;
        res.chip = (res.cellid / 10000) % 10;
        res.channel = res.cellid % 10000;
        ntrack_pass_through_channel = -1;  // Reset for each histogram
        if (NtrackPassThroughChannel.find(res.cellid) != NtrackPassThroughChannel.end()) {
            ntrack_pass_through_channel = NtrackPassThroughChannel[res.cellid];
        }

        FitOut fitOut = fitLandauGaus(h, 300, false, false);  // draw_fit=false
        
        res.mpv = fitOut.mpv;
        res.width = fitOut.width;
        res.total_area = fitOut.total_area;
        res.gaus_sigma = fitOut.gaus_sigma;
        res.chi2 = fitOut.chi2;
        res.ndf = fitOut.ndf;
        res.fit_status = fitOut.fit_status;
        res.fit_ok = fitOut.fit_ok ? 1 : 0;
        res.entries_adc_le50 = fitOut.entries_adc_le50;
        res.percentile_1 = fitOut.percentile_1;
        res.percentile_2 = fitOut.percentile_2;
        res.percentile_3 = fitOut.percentile_3;
        res.percentile_4 = fitOut.percentile_4;
        res.percentile_5 = fitOut.percentile_5;
        res.percentile_10 = fitOut.percentile_10;
        res.percentile_15 = fitOut.percentile_15;

        if (res.entries > 0) {
            res.ratio_adc_le50 = double(res.entries_adc_le50) / double(res.entries);
        } else {
            res.ratio_adc_le50 = -1.0;
        }
        res.fit_Class = fitOut.fit_Class;
        resultTree->Fill();
        
        delete h;
        h = nullptr;
        
        fitCount++;
        if (fitOut.fit_ok) fitOkCount++;

        if (fitCount % 100 == 0) {
            std::cout << "Processed " << fitCount << " histograms, " << fitOkCount << " successful fits" << std::endl;
            // Write partial results to disk for persistence on forced kill
            fout->cd();
            resultTree->Write(0, TObject::kOverwrite);
            fout->Flush();
        }
    }

    std::cout << "\nFitting complete:" << std::endl;
    std::cout << "  Total histograms: " << fitCount << std::endl;
    std::cout << "  Successful fits: " << fitOkCount << std::endl;

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

int main(int argc, char** argv) {
    if (argc == 1) {
        fit_histograms_directfit_integral();  // Default mode with plotting
    } else if (argc == 3) {
         fit_histograms_noplot(argv[1], argv[2]);
    } else if (argc == 5 || argc == 6) {
        fit_histograms_directfit_integral(
            argv[1], argv[2], argv[3], argv[4],
            argc == 6 ? argv[5] : "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out/SimMuon/mip_neighborcheck_nofit.root");
    } else {
        std::cout << "Usage:\n"
                  << "  " << argv[0] << " <input_file> <output_file>\n"
                  << "  " << argv[0] << " <input_file> <output_file> <cellid_file> <plot_dir> [efficiency_input_file]\n";
    }
    return 0;
}
