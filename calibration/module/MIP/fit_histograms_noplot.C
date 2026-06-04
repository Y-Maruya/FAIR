/**
 * Standalone script to fit histograms from ROOT file (no draw)
 * Usage: root -b -q 'fit_histograms_noplot.C("input.root", "output.root")'
 * 
 * Memory Management:
 * - TF1* fLandauGaus: Explicitly deleted inside the fitting helpers
 * - TH1D* h_copy: Deleted after percentile calculation
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
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TIterator.h>
#include <vector>
#include <string>
#include <TParameter.h>
#include <cmath>
#include <TGraphErrors.h>
#include <unordered_map>

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
    double ratio_under_0ADC = -1.0;
    int entries_adc_le50 = 0;
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
    double ratio_under_0ADC = -1.0;
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
    // Refit results
    bool refit = false;
    double refit_mpv = -1.0;
    double refit_width = -1.0;
    double refit_total_area = -1.0;
    double refit_gaus_sigma = -1.0;
    double refit_chi2 = -1.0;
    int refit_ndf = 0;
    double refit_chi2_ndf = -1.0;
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
    double threshold = -1.0;
    double threshold_width = -1.0;
    double threshold_chi2 = -1.0;
    int threshold_ndf = 0;
    double threshold_chi2_ndf = -1.0;
    double refit_ratio_under_0ADC = -1.0;
    int n_refit_threshold_fit = 0;
    int num_points_in_transition_region = 0;
};

struct RefitResult {
    int cellid = -1;
    double refit_mpv = -1.0;
    double refit_width = -1.0;
    double refit_total_area = -1.0;
    double refit_gaus_sigma = -1.0;
    double refit_chi2 = -1.0;
    int refit_ndf = 0;
    double refit_chi2_ndf = -1.0;
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
    double threshold = -1.0;
    double threshold_width = -1.0;
    double threshold_chi2 = -1.0;
    int threshold_ndf = 0;
    double threshold_chi2_ndf = -1.0;
    double refit_ratio_under_0ADC = -1.0;
    int n_refit_threshold_fit = 0;
    int num_points_in_transition_region = 0;
};

void calculateClassification(FitOut& fitOut) {
    if (fitOut.fit_status != 0) {
        fitOut.fit_Class = FitClass_fitFailed;
    } else if (fitOut.entries < 200) {
        fitOut.fit_Class = FitClass_LessStatistics;
    } else if (fitOut.ndf > 0 && fitOut.chi2 / double(fitOut.ndf) > 3 || fitOut.chi2 < 0) {
        fitOut.fit_Class = FitClass_BadChi2;
    } else if (fitOut.mpv <= 0 || fitOut.width <= 0 || fitOut.total_area <= 0 || fitOut.gaus_sigma <= 0 || (fitOut.gaus_sigma > 149 && fitOut.width > 89)) {
        fitOut.fit_Class = FitClass_ParameterAtLimit;
    } else if (fitOut.gaus_sigma < 20 && fitOut.width < 50) {
        fitOut.fit_Class = FitClass_gaussigma_small_Threshold;
    } else if (fitOut.gaus_sigma < 20 && fitOut.width >= 50) {
        fitOut.fit_Class = FitClass_gaussigma_small_Retry;
    } else {
        fitOut.fit_Class = FitClass_GoodFit;
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
    
    // Calculate percentiles and entries_adc_le50
    TH1D* h_copy = (TH1D*)h->Clone("h_copy");
    h_copy->SetDirectory(nullptr);
    for (int i = 1; i <= h_copy->GetNbinsX(); ++i) {
        if (h_copy->GetBinCenter(i) <= 50) {
            r.entries_adc_le50 += h_copy->GetBinContent(i);
            h_copy->SetBinContent(i, 0);
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
        if (sum >= r.entries * 0.05) {
            r.percentile_5 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries * 0.10) {
            r.percentile_10 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries * 0.15) {
            r.percentile_15 = h_copy->GetBinCenter(i);
             break;
        }
    }
    delete h_copy;
    h_copy = nullptr;

    double chi2_ndf = (ndf > 0) ? (chisqr / double(ndf)) : -1.0;
    if (chi2_ndf > 3 || chi2_ndf < 0) {
        r.fit_status = fitResult;
        r.fit_ok = false;
    } else {
        r.fit_status = fitResult;
        r.fit_ok = true;
    }
    calculateClassification(r);
    // std::cout<< "before integral under 0 ADC: total_area = " << r.total_area << std::endl;
    r.ratio_under_0ADC = (r.total_area > 0) ? double(fLandauGaus->Integral(-96, 0)) / double(r.total_area) : -1.0;
    // std::cout << "Fit results: MPV = " << r.mpv << ", Width = " << r.width << ", Total Area = " << r.total_area 
    //           << ", Gauss Sigma = " << r.gaus_sigma << ", Chi2/NDF = " << chi2_ndf 
    //           << ", Fit Class = " << r.fit_Class << std::endl;
    delete fLandauGaus;
    fLandauGaus = nullptr;
    if (r.gaus_sigma < 1){
        FitOut FitO = fitLandauGaus(h, minEntries, calculate_fwhm, 20);
        // std::cout << "Initial fit had gaus_sigma = " << r.gaus_sigma << ", retrying with min_gaus_sigma = 20" << std::endl;
        if (FitO.fit_ok && FitO.chi2/FitO.ndf < chisqr/ndf && FitO.gaus_sigma > 20+0.1) {
            r = FitO;
            // std::cout << "Retry fit successful with gaus_sigma = " << FitO.gaus_sigma << std::endl;
        }
    }
    return r;
}

RefitResult fitLandauGausWithThreshold(TH1D* h, int minEntries, double min_gaus_sigma, FitResult& res) {
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
    sv[2] = res.ntrack_pass_through_channel > 0 ? res.ntrack_pass_through_channel * h->GetBinWidth(1) : 30000;
    pllo[2] = res.ntrack_pass_through_channel > 0 ? res.ntrack_pass_through_channel * h->GetBinWidth(1) : 100;
    plhi[2] = res.ntrack_pass_through_channel > 0 ? res.ntrack_pass_through_channel * h->GetBinWidth(1) : 100000;

    fr[0] = res.mpv > 0 ? res.mpv : 50;
    fr[1] = res.mpv > 0 ? res.mpv * 2.5 : 1000;
    int fitResult = 0;
    // std::cout << "Refitting histogram: " << h->GetName() << " with initial parameters: MPV = " << fr[0] << ", Max Range = " << fr[1] 
        //   << ", Total Area initial guess = " << sv[2] << std::endl;
    TF1* fLandauGaus = langaufit(h, fr, sv, pllo, plhi, fps, fpe, &chisqr, &ndf, &fitResult);
    if (!fLandauGaus) {
        r.fit_status = -3;
        r.fit_ok = false;
        return rr;
    }
    r.mpv = fps[1];
    r.width = fps[0];
    r.total_area = fps[2];
    r.gaus_sigma = fps[3];
    r.entries = h->GetEntries();
    r.chi2 = chisqr;
    r.ndf = ndf;
    r.fit_status = fitResult;
    r.fit_ok = (chisqr / double(ndf) <= 5 && chisqr >= 0) ? true : false;
    // std::cout << "Refit results: MPV = " << r.mpv << ", Width = " << r.width << ", Total Area = " << r.total_area 
            //   << ", Gauss Sigma = " << r.gaus_sigma << ", Chi2/NDF = " << ((ndf > 0) ? (chisqr / double(ndf)) : -1.0) 
            //   << ", Fit Class = " << r.fit_Class << std::endl;
    r.ratio_under_0ADC = (r.total_area > 0) ? double(fLandauGaus->Integral(-96, 0)) / double(r.total_area) : -1.0;
    // std::cout << "Refit ratio under 0 ADC: " << r.ratio_under_0ADC << std::endl;
    // Calculate percentiles
    TH1D* h_copy = (TH1D*)h->Clone("h_copy");
    h_copy->SetDirectory(nullptr);
    for (int i = 1; i <= h_copy->GetNbinsX(); ++i) {
        if (h_copy->GetBinCenter(i) <= 50) {
            h_copy->SetBinContent(i, 0);
        }
    }
    int sum = 0;
    for (int i = 1; i <= h_copy->GetNbinsX(); ++i) {
        sum += h_copy->GetBinContent(i);
        if (sum >= r.entries * 0.01) {
            r.percentile_1 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries * 0.05) {
            r.percentile_5 = h_copy->GetBinCenter(i);
        }
        if (sum >= r.entries * 0.10) {
            r.percentile_10 = h_copy->GetBinCenter(i);
             break;
        }
    }
    delete h_copy;
    h_copy = nullptr;

    // Create ratio histogram for threshold calculation
    TH1D* h_ratio = (TH1D*)h->Clone("h_ratio");
    h_ratio->SetDirectory(nullptr);
    for (int i = 1; i <= h_ratio->GetNbinsX(); ++i) {
        if (fLandauGaus->Eval(h->GetBinCenter(i)) > 1 && h->GetBinCenter(i) > 50){
            h_ratio->SetBinContent(i, h->GetBinContent(i) / fLandauGaus->Eval(h->GetBinCenter(i)));
            h_ratio->SetBinError(i, h->GetBinError(i) / fLandauGaus->Eval(h->GetBinCenter(i)));
        }else{
            h_ratio->SetBinContent(i, 0);
            h_ratio->SetBinError(i, 0);
        }
    }
    
    double chi2_ndf = (ndf > 0) ? (chisqr / double(ndf)) : -1.0;
    if (chi2_ndf > 5 || chi2_ndf < 0) {
        r.fit_status = fitResult;
        r.fit_ok = false;
    } else {
        r.fit_status = fitResult;
        r.fit_ok = true;
    }
    calculateClassification(r);
    // Fill RefitResult struct
    rr.cellid = res.cellid;
    rr.refit_mpv = r.mpv;
    rr.refit_width = r.width;
    rr.refit_total_area = r.total_area;
    rr.refit_gaus_sigma = r.gaus_sigma;
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
    
    // Fit error function to ratio histogram to find threshold
    TF1 *f_erf = new TF1("f_erf", "0.5*(1+TMath::Erf((x-[0])/(sqrt(2)*[1])))", 10, h->GetBinCenter(h->GetMaximumBin()+50));
    f_erf->SetParameters(r.percentile_1, 10);
    // f_erf->SetParLimits(1, 1, 50);
    
    // Fit using TGraphErrors to get better fit quality
    int rebin = 2;
    TH1D* h_ratio_rebinned = (TH1D*)h_ratio->Rebin(rebin, "h_ratio_rebinned");
    h_ratio_rebinned->SetDirectory(nullptr);
    TGraphErrors* gRatio = new TGraphErrors(h_ratio_rebinned);
    for (int i = 0; i < gRatio->GetN(); ++i) {
        double x, y;
        gRatio->GetPoint(i, x, y);
        gRatio->SetPoint(i, x, y/rebin); // Normalize by rebin factor
        gRatio->SetPointError(i, 0, h_ratio_rebinned->GetBinError(i+1)/rebin);
    }   
    gRatio->Fit(f_erf, "RQ");
    double chi2_ndf_threshold = f_erf->GetNDF() > 0 ? f_erf->GetChisquare() / f_erf->GetNDF() : -1.0;
    rr.n_refit_threshold_fit = 1;
    int num_points_in_transition = 0;
    if (chi2_ndf_threshold > 3){
        TH1D* h_ratio_rebinned_2 = (TH1D*)h_ratio_rebinned->Rebin(2, "h_ratio_rebinned_2");
        h_ratio_rebinned_2->SetDirectory(nullptr);
        TGraphErrors* gRatio2 = new TGraphErrors(h_ratio_rebinned_2);
        for (int i = 0; i < gRatio2->GetN(); ++i) {
            double x, y;
            gRatio2->GetPoint(i, x, y);
            gRatio2->SetPoint(i, x, y/4); // Normalize by total rebin factor
            gRatio2->SetPointError(i, 0, h_ratio_rebinned_2->GetBinError(i+1)/4);
        }
        TF1 *f_erf2 = new TF1("f_erf2", "0.5*(1+TMath::Erf((x-[0])/(sqrt(2)*[1])))", 10, h->GetBinCenter(h->GetMaximumBin()+50));
        f_erf2->SetParameters(r.percentile_1, 10);
        // f_erf2->SetParLimits(1, 1, 50);
        gRatio2->Fit(f_erf2, "RQ");
        // gRatio2->Fit(f_erf, "RQ");
        
        // if the number of points in the transition region is too small after rebinning,
        // the fit may fail or give unreliable results, so we can choose to skip the threshold calculation in that case
        // num_points_in_transition = 0;
        // for (int i = 0; i < gRatio2->GetN(); ++i) {
        //     double x, y;
        //     gRatio2->GetPoint(i, x, y);
        //     if (x > f_erf2->GetParameter(0) - 2*f_erf2->GetParameter(1) && x < f_erf2->GetParameter(0) + 2*f_erf2->GetParameter(1)) {
        //         num_points_in_transition++;
        //     }
        // }
        // if (num_points_in_transition > 3) { // Require at least 4 points in the transition region for reliable fit
        delete f_erf;
        f_erf = (TF1*)f_erf2->Clone("f_erf");
        rr.n_refit_threshold_fit = 2;
        delete gRatio2;
        gRatio2 = nullptr;
        delete h_ratio_rebinned_2;
        h_ratio_rebinned_2 = nullptr;
        delete f_erf2;
        f_erf2 = nullptr;
    }
    
    rr.threshold = f_erf->GetParameter(0);
    rr.threshold_width = f_erf->GetParameter(1);
    rr.threshold_chi2 = f_erf->GetChisquare();
    rr.threshold_ndf = f_erf->GetNDF();
    rr.threshold_chi2_ndf = (f_erf->GetNDF() > 0) ? (f_erf->GetChisquare() / f_erf->GetNDF()) : -1.0; 
    rr.num_points_in_transition_region = 0;
    for (int i = 0; i < gRatio->GetN(); ++i) {
        double x, y;
        gRatio->GetPoint(i, x, y);
        if (x > rr.threshold - 2*rr.threshold_width && x < rr.threshold + 2*rr.threshold_width) {
            rr.num_points_in_transition_region++;
        }
    }
    // std::cout << "Refit threshold: " << rr.threshold << ", width: " << rr.threshold_width 
    //           << ", chi2/ndf: " << rr.threshold_chi2_ndf << std::endl;
    
    delete gRatio;
    gRatio = nullptr;
    delete f_erf;
    f_erf = nullptr;
    delete h_ratio;
    h_ratio = nullptr;
    delete h_ratio_rebinned;
    h_ratio_rebinned = nullptr;

    delete fLandauGaus;
    fLandauGaus = nullptr;
    
    return rr;
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
            h->Rebin(4);
            h->SetDirectory(nullptr);
            histograms.push_back(h);
        } else if (obj->IsA()->InheritsFrom(TParameter<int>::Class())) {
            dir->cd();
            TParameter<int>* p = (TParameter<int>*)obj;
            std::string name = p->GetName();
            if (name.find("Ntracks_pass_through_channel_") == 0) {
                int cellid = std::stoi(name.substr(std::string("Ntracks_pass_through_channel_").length()));
                if (NtrackPassThroughChannel.find(cellid) == NtrackPassThroughChannel.end()) {
                    NtrackPassThroughChannel[cellid] = p->GetVal();
                }
            }
        }
    }
}

void fit_histograms_noplot(const char* input_file = "mip.root", const char* output_file = "mip_fitted.root") {
    std::cout << "Opening input file (no plot mode): " << input_file << std::endl;
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

    // Create result tree
    TTree* resultTree = new TTree("mip_fit_results", "MIP fit results");
    FitResult res;
    int ntrack_pass_through_channel = -1;
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
    resultTree->Branch("ratio_under_0ADC", &res.ratio_under_0ADC);
    resultTree->Branch("ntrack_pass_through_channel", &ntrack_pass_through_channel);
    resultTree->Branch("FitClass", &res.fit_Class);
    // Refit branches
    resultTree->Branch("refit", &res.refit);
    resultTree->Branch("refit_mpv", &res.refit_mpv);
    resultTree->Branch("refit_width", &res.refit_width);
    resultTree->Branch("refit_total_area", &res.refit_total_area);
    resultTree->Branch("refit_gaus_sigma", &res.refit_gaus_sigma);
    resultTree->Branch("refit_chi2", &res.refit_chi2);
    resultTree->Branch("refit_ndf", &res.refit_ndf);
    resultTree->Branch("refit_chi2_ndf", &res.refit_chi2_ndf);
    resultTree->Branch("refit_fit_status", &res.refit_fit_status);
    resultTree->Branch("refit_fit_ok", &res.refit_fit_ok);
    resultTree->Branch("refit_percentile_1", &res.refit_percentile_1);
    resultTree->Branch("refit_percentile_2", &res.refit_percentile_2);
    resultTree->Branch("refit_percentile_3", &res.refit_percentile_3);
    resultTree->Branch("refit_percentile_4", &res.refit_percentile_4);
    resultTree->Branch("refit_percentile_5", &res.refit_percentile_5);
    resultTree->Branch("refit_percentile_10", &res.refit_percentile_10);
    resultTree->Branch("refit_Class", &res.refit_Class);
    resultTree->Branch("threshold", &res.threshold);
    resultTree->Branch("threshold_width", &res.threshold_width);
    resultTree->Branch("threshold_chi2", &res.threshold_chi2);
    resultTree->Branch("threshold_ndf", &res.threshold_ndf);
    resultTree->Branch("threshold_chi2_ndf", &res.threshold_chi2_ndf);
    resultTree->Branch("refit_ratio_under_0ADC", &res.refit_ratio_under_0ADC);
    resultTree->Branch("refit_Class", &res.refit_Class);
    resultTree->Branch("n_refit_threshold_fit", &res.n_refit_threshold_fit);
    resultTree->Branch("num_points_in_transition_region", &res.num_points_in_transition_region);

    int fitCount = 0;
    int fitOkCount = 0;
    int refitCount = 0;
    std::cout << "\n=== Fitting histograms (no plot) ===" << std::endl;
    
    for (TH1D* h : histograms) {
        const char* histName = h->GetName();
        res.entries = h->GetEntries();

        sscanf(histName, "hMIP_%d", &res.cellid);
        res.layer = res.cellid / 100000;
        res.chip = (res.cellid / 10000) % 10;
        res.channel = res.cellid % 10000;
        ntrack_pass_through_channel = -1;
        if (NtrackPassThroughChannel.find(res.cellid) != NtrackPassThroughChannel.end()) {
            ntrack_pass_through_channel = NtrackPassThroughChannel[res.cellid];
        }

        FitOut fitOut = fitLandauGaus(h, 200, false);
        
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
        res.ratio_under_0ADC = fitOut.ratio_under_0ADC;
        // std::cout << "CellID: " << res.cellid << ", MPV: " << res.mpv << ", Width: " << res.width 
        //           << ", Total Area: " << res.total_area << ", Gauss Sigma: " << res.gaus_sigma 
        //           << ", Chi2: " << res.chi2 << ", NDF: " << res.ndf 
        //           << ", Fit OK: " << res.fit_ok 
        //           << ", Entries <= 50 ADC: " << res.entries_adc_le50 
        //           << ", Ratio under 0 ADC: " << res.ratio_under_0ADC 
        //           << std::endl;

        if (res.entries > 0) {
            res.ratio_adc_le50 = double(res.entries_adc_le50) / double(res.entries);
        } else {
            res.ratio_adc_le50 = -1.0;
        }
        res.fit_Class = fitOut.fit_Class;
        res.ntrack_pass_through_channel = ntrack_pass_through_channel;
        
        // Initialize refit fields
        res.refit = false;
        res.refit_mpv = -1.0;
        res.refit_width = -1.0;
        res.refit_total_area = -1.0;
        res.refit_gaus_sigma = -1.0;
        res.refit_chi2 = -1.0;
        res.refit_ndf = 0;
        res.refit_chi2_ndf = -1.0;
        res.refit_fit_status = 999;
        res.refit_fit_ok = 0;
        res.refit_percentile_1 = -1.0;
        res.refit_percentile_5 = -1.0;
        res.refit_percentile_10 = -1.0;
        res.threshold = -1.0;
        res.threshold_width = -1.0;
        res.threshold_chi2 = -1.0;
        res.threshold_ndf = 0;
        res.threshold_chi2_ndf = -1.0;
        res.refit_Class = -1;
        res.n_refit_threshold_fit = 0;
        res.num_points_in_transition_region = 0;
        // std::cout << "Initial fit class: " << res.fit_Class << std::endl;
        // Perform refit if small gaus_sigma
        if (fitOut.fit_Class == FitClass_gaussigma_small_Threshold || fitOut.fit_Class == FitClass_gaussigma_small_Retry) {

            std::cout << "Performing refit for cellid " << res.cellid << " with small gaus_sigma = " << fitOut.gaus_sigma << std::endl;
            RefitResult rr_result = fitLandauGausWithThreshold(h, 200, 1, res);
            res.refit = true;
            res.refit_mpv = rr_result.refit_mpv;
            res.refit_width = rr_result.refit_width;
            res.refit_total_area = rr_result.refit_total_area;
            res.refit_gaus_sigma = rr_result.refit_gaus_sigma;
            res.refit_chi2 = rr_result.refit_chi2;
            res.refit_ndf = rr_result.refit_ndf;
            res.refit_chi2_ndf = rr_result.refit_chi2_ndf;
            res.refit_fit_status = rr_result.refit_fit_status;
            res.refit_fit_ok = rr_result.refit_fit_ok;
            res.refit_percentile_1 = rr_result.refit_percentile_1;
            res.refit_percentile_2 = rr_result.refit_percentile_2;
            res.refit_percentile_3 = rr_result.refit_percentile_3;
            res.refit_percentile_4 = rr_result.refit_percentile_4;
            res.refit_percentile_5 = rr_result.refit_percentile_5;
            res.refit_percentile_10 = rr_result.refit_percentile_10;
            res.refit_percentile_15 = rr_result.refit_percentile_15;
            res.refit_ratio_under_0ADC = rr_result.refit_ratio_under_0ADC;
            res.threshold = rr_result.threshold;
            res.threshold_width = rr_result.threshold_width;
            res.threshold_chi2 = rr_result.threshold_chi2;
            res.threshold_ndf = rr_result.threshold_ndf;
            res.threshold_chi2_ndf = rr_result.threshold_chi2_ndf;
            res.refit_Class = rr_result.refit_Class;
            res.n_refit_threshold_fit = rr_result.n_refit_threshold_fit;
            res.num_points_in_transition_region = rr_result.num_points_in_transition_region;
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
    fit_histograms_noplot("/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out_calibration/mip_calib/21987-22133/mip_nofit.root", "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out_calibration/mip_calib/21987-22133/mip_fitted_noplot_7.root");
}

int main(int argc, char** argv) {
    if (argc == 1) {
        fit_histograms();
    } else if (argc == 3) {
        fit_histograms_noplot(argv[1], argv[2]);
    } else {
        std::cout << "Usage: " << argv[0] << " <input_file> <output_file>" << std::endl;
        std::cout << "Example: " << argv[0] << " mip.root mip_fitted.root" << std::endl;
    }
    return 0;
}
