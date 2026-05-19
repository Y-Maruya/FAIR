#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"

// ============================================================================
// Data Structures
// ============================================================================

struct PedestalData {
    int cellid = -1;
    double highgain_sigma = 0.0;
    double highgain_rms = 0.0;
    double lowgain_sigma = 0.0;
    double lowgain_rms = 0.0;
};

struct IntercalibData {
    int cellid = -1;
    double hg_sigma = 0.0;
    double hg_rms = 0.0;
    double lg_sigma = 0.0;
    double lg_rms = 0.0;
};

struct ComparisonStats {
    int n_matches = 0;
    double ped_mean = 0.0;
    double ped_rms = 0.0;
    double ped_min = 1e9;
    double ped_max = -1e9;
    double ic_mean = 0.0;
    double ic_rms = 0.0;
    double ic_min = 1e9;
    double ic_max = -1e9;
    double correlation = 0.0;
};

// ============================================================================
// Main Analysis Function
// ============================================================================

void compare_calibration(
    const char* pedestal_file = "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out_calibration/pedestal_calib/22140/pedestal.root",
    const char* intercalib_file = "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out_calibration/intercalib_calib/22137-22142/intercalib_adj.root",
    const char* output_dir = "./compare_calibration_results"
) {
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Calibration Comparison Analysis" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // ========================================================================
    // Phase 1: Load Pedestal Data
    // ========================================================================
    
    std::cout << "[1/4] Loading pedestal data from: " << pedestal_file << std::endl;
    
    TFile* ped_file = TFile::Open(pedestal_file, "READ");
    if (!ped_file || ped_file->IsZombie()) {
        std::cerr << "ERROR: Could not open pedestal file!" << std::endl;
        return;
    }
    
    TTree* ped_tree = (TTree*)ped_file->Get("pedestal");
    if (!ped_tree) {
        std::cerr << "ERROR: Could not find 'pedestal' TTree!" << std::endl;
        return;
    }
    
    std::map<int, PedestalData> pedestalMap;
    
    int ped_cellid;
    double ped_hg_sigma, ped_hg_rms, ped_lg_sigma, ped_lg_rms;
    
    ped_tree->SetBranchAddress("cellid", &ped_cellid);
    ped_tree->SetBranchAddress("highgain_sigma", &ped_hg_sigma);
    ped_tree->SetBranchAddress("highgain_rms", &ped_hg_rms);
    ped_tree->SetBranchAddress("lowgain_sigma", &ped_lg_sigma);
    ped_tree->SetBranchAddress("lowgain_rms", &ped_lg_rms);
    
    Long64_t n_ped_entries = ped_tree->GetEntries();
    std::cout << "  → Found " << n_ped_entries << " pedestal entries" << std::endl;
    
    for (Long64_t i = 0; i < n_ped_entries; i++) {
        ped_tree->GetEntry(i);
        PedestalData data;
        data.cellid = ped_cellid;
        data.highgain_sigma = ped_hg_sigma;
        data.highgain_rms = ped_hg_rms;
        data.lowgain_sigma = ped_lg_sigma;
        data.lowgain_rms = ped_lg_rms;
        pedestalMap[ped_cellid] = data;
    }
    
    std::cout << "  ✓ Loaded " << pedestalMap.size() << " unique pedestal channels" << std::endl;
    
    // ========================================================================
    // Phase 1b: Load Intercalib Data
    // ========================================================================
    
    std::cout << "\n[2/4] Loading intercalib data from: " << intercalib_file << std::endl;
    
    TFile* ic_file = TFile::Open(intercalib_file, "READ");
    if (!ic_file || ic_file->IsZombie()) {
        std::cerr << "ERROR: Could not open intercalib file!" << std::endl;
        return;
    }
    
    TTree* ic_tree = (TTree*)ic_file->Get("intercalib");
    if (!ic_tree) {
        std::cerr << "ERROR: Could not find 'intercalib' TTree!" << std::endl;
        return;
    }
    
    std::map<int, IntercalibData> intercalibMap;
    
    int ic_cellid;
    double ic_hg_sigma, ic_hg_rms, ic_lg_sigma, ic_lg_rms;
    
    ic_tree->SetBranchAddress("cellid", &ic_cellid);
    ic_tree->SetBranchAddress("hg_sigma", &ic_hg_sigma);
    ic_tree->SetBranchAddress("hg_rms", &ic_hg_rms);
    ic_tree->SetBranchAddress("lg_sigma", &ic_lg_sigma);
    ic_tree->SetBranchAddress("lg_rms", &ic_lg_rms);
    
    Long64_t n_ic_entries = ic_tree->GetEntries();
    std::cout << "  → Found " << n_ic_entries << " intercalib entries" << std::endl;
    
    for (Long64_t i = 0; i < n_ic_entries; i++) {
        ic_tree->GetEntry(i);
        IntercalibData data;
        data.cellid = ic_cellid;
        data.hg_sigma = ic_hg_sigma;
        data.hg_rms = ic_hg_rms;
        data.lg_sigma = ic_lg_sigma;
        data.lg_rms = ic_lg_rms;
        intercalibMap[ic_cellid] = data;
    }
    
    std::cout << "  ✓ Loaded " << intercalibMap.size() << " unique intercalib channels" << std::endl;
    
    // ========================================================================
    // Phase 2: Find Matching Channels and Collect Data
    // ========================================================================
    
    std::cout << "\n[3/4] Matching channels from both files..." << std::endl;
    
    // Storage for matching data
    std::vector<double> ped_hg_sigma_vec, ic_hg_sigma_vec;
    std::vector<double> ped_hg_rms_vec, ic_hg_rms_vec;
    std::vector<double> ped_lg_sigma_vec, ic_lg_sigma_vec;
    std::vector<double> ped_lg_rms_vec, ic_lg_rms_vec;
    
    int n_matches = 0;
    
    // Find matching channels
    for (const auto& ped_pair : pedestalMap) {
        int cellid = ped_pair.first;
        if (intercalibMap.find(cellid) != intercalibMap.end()) {
            n_matches++;
            
            const PedestalData& ped_data = ped_pair.second;
            const IntercalibData& ic_data = intercalibMap[cellid];
            
            ped_hg_sigma_vec.push_back(ped_data.highgain_sigma);
            ic_hg_sigma_vec.push_back(ic_data.hg_sigma);
            
            ped_hg_rms_vec.push_back(ped_data.highgain_rms);
            ic_hg_rms_vec.push_back(ic_data.hg_rms);
            
            ped_lg_sigma_vec.push_back(ped_data.lowgain_sigma);
            ic_lg_sigma_vec.push_back(ic_data.lg_sigma);
            
            ped_lg_rms_vec.push_back(ped_data.lowgain_rms);
            ic_lg_rms_vec.push_back(ic_data.lg_rms);
        }
    }
    
    std::cout << "  ✓ Found " << n_matches << " matching channels" << std::endl;
    
    if (n_matches == 0) {
        std::cerr << "ERROR: No matching channels found!" << std::endl;
        return;
    }
    
    // ========================================================================
    // Phase 3: Create Histograms
    // ========================================================================
    
    std::cout << "\n[4/4] Creating and filling comparison histograms..." << std::endl;
    
    gStyle->SetOptStat(111111);
    gStyle->SetOptFit(1111);
    
    // Determine bin ranges for each comparison
    auto get_range = [](const std::vector<double>& v) -> std::pair<double, double> {
        double vmin = *std::min_element(v.begin(), v.end());
        double vmax = *std::max_element(v.begin(), v.end());
        double margin = (vmax - vmin) * 0.1;
        return {vmin - margin, vmax + margin};
    };
    
    auto range_hg_sigma = get_range(ped_hg_sigma_vec);
    auto range_ic_hg_sigma = get_range(ic_hg_sigma_vec);
    auto range_hg_rms = get_range(ped_hg_rms_vec);
    auto range_ic_hg_rms = get_range(ic_hg_rms_vec);
    auto range_lg_sigma = get_range(ped_lg_sigma_vec);
    auto range_ic_lg_sigma = get_range(ic_lg_sigma_vec);
    auto range_lg_rms = get_range(ped_lg_rms_vec);
    auto range_ic_lg_rms = get_range(ic_lg_rms_vec);
    
    // Create 2D histograms
    TH2D* h_hg_sigma = new TH2D("h_hg_sigma",
        "High Gain Sigma: Pedestal vs Intercalib;Pedestal HG Sigma;Intercalib HG Sigma",
        100, range_hg_sigma.first, range_hg_sigma.second,
        100, range_ic_hg_sigma.first, range_ic_hg_sigma.second);
    
    TH2D* h_hg_rms = new TH2D("h_hg_rms",
        "High Gain RMS: Pedestal vs Intercalib;Pedestal HG RMS;Intercalib HG RMS",
        100, range_hg_rms.first, range_hg_rms.second,
        100, range_ic_hg_rms.first, range_ic_hg_rms.second);
    
    TH2D* h_lg_sigma = new TH2D("h_lg_sigma",
        "Low Gain Sigma: Pedestal vs Intercalib;Pedestal LG Sigma;Intercalib LG Sigma",
        100, range_lg_sigma.first, range_lg_sigma.second,
        100, range_ic_lg_sigma.first, range_ic_lg_sigma.second);
    
    TH2D* h_lg_rms = new TH2D("h_lg_rms",
        "Low Gain RMS: Pedestal vs Intercalib;Pedestal LG RMS;Intercalib LG RMS",
        100, range_lg_rms.first, range_lg_rms.second,
        100, range_ic_lg_rms.first, range_ic_lg_rms.second);
    
    // Fill histograms
    for (int i = 0; i < n_matches; i++) {
        h_hg_sigma->Fill(ped_hg_sigma_vec[i], ic_hg_sigma_vec[i]);
        h_hg_rms->Fill(ped_hg_rms_vec[i], ic_hg_rms_vec[i]);
        h_lg_sigma->Fill(ped_lg_sigma_vec[i], ic_lg_sigma_vec[i]);
        h_lg_rms->Fill(ped_lg_rms_vec[i], ic_lg_rms_vec[i]);
    }
    
    std::cout << "  ✓ Created and filled 4 histograms" << std::endl;
    
    // ========================================================================
    // Phase 4: Calculate Statistics
    // ========================================================================
    
    auto calc_stats = [](const std::vector<double>& v1, const std::vector<double>& v2) -> ComparisonStats {
        ComparisonStats stats;
        stats.n_matches = v1.size();
        
        // Calculate means
        double sum1 = 0, sum2 = 0;
        double min1 = 1e9, max1 = -1e9;
        double min2 = 1e9, max2 = -1e9;
        
        for (int i = 0; i < stats.n_matches; i++) {
            sum1 += v1[i];
            sum2 += v2[i];
            min1 = std::min(min1, v1[i]);
            max1 = std::max(max1, v1[i]);
            min2 = std::min(min2, v2[i]);
            max2 = std::max(max2, v2[i]);
        }
        
        stats.ped_mean = sum1 / stats.n_matches;
        stats.ic_mean = sum2 / stats.n_matches;
        stats.ped_min = min1;
        stats.ped_max = max1;
        stats.ic_min = min2;
        stats.ic_max = max2;
        
        // Calculate RMS
        double sum_sq1 = 0, sum_sq2 = 0;
        for (int i = 0; i < stats.n_matches; i++) {
            sum_sq1 += (v1[i] - stats.ped_mean) * (v1[i] - stats.ped_mean);
            sum_sq2 += (v2[i] - stats.ic_mean) * (v2[i] - stats.ic_mean);
        }
        
        stats.ped_rms = std::sqrt(sum_sq1 / stats.n_matches);
        stats.ic_rms = std::sqrt(sum_sq2 / stats.n_matches);
        
        // Calculate correlation
        double numerator = 0;
        for (int i = 0; i < stats.n_matches; i++) {
            numerator += (v1[i] - stats.ped_mean) * (v2[i] - stats.ic_mean);
        }
        numerator /= stats.n_matches;
        
        if (stats.ped_rms > 0 && stats.ic_rms > 0) {
            stats.correlation = numerator / (stats.ped_rms * stats.ic_rms);
        }
        
        return stats;
    };
    
    ComparisonStats stat_hg_sigma = calc_stats(ped_hg_sigma_vec, ic_hg_sigma_vec);
    ComparisonStats stat_hg_rms = calc_stats(ped_hg_rms_vec, ic_hg_rms_vec);
    ComparisonStats stat_lg_sigma = calc_stats(ped_lg_sigma_vec, ic_lg_sigma_vec);
    ComparisonStats stat_lg_rms = calc_stats(ped_lg_rms_vec, ic_lg_rms_vec);
    
    // ========================================================================
    // Phase 5: Save Results to ROOT File
    // ========================================================================
    
    std::cout << "\n[5/5] Saving results..." << std::endl;
    
    std::string output_root = std::string(output_dir) + "/calibration_comparison.root";
    TFile* output_file = TFile::Open(output_root.c_str(), "RECREATE");
    
    h_hg_sigma->Write();
    h_hg_rms->Write();
    h_lg_sigma->Write();
    h_lg_rms->Write();
    
    output_file->Close();
    std::cout << "  ✓ Saved histograms to: " << output_root << std::endl;
    
    // ========================================================================
    // Phase 6: Generate PNG Plots
    // ========================================================================
    
    std::cout << "\nGenerating PNG plots..." << std::endl;
    
    gStyle->SetPalette(55);  // Rainbow palette
    
    // HG Sigma plot
    TCanvas* c_hg_sigma = new TCanvas("c_hg_sigma", "HG Sigma Comparison", 800, 600);
    h_hg_sigma->Draw("COLZ");
    h_hg_sigma->SetStats(1);
    std::string png_hg_sigma = std::string(output_dir) + "/hg_sigma_comparison.png";
    c_hg_sigma->SaveAs(png_hg_sigma.c_str());
    std::cout << "  ✓ Saved: " << png_hg_sigma << std::endl;
    
    // HG RMS plot
    TCanvas* c_hg_rms = new TCanvas("c_hg_rms", "HG RMS Comparison", 800, 600);
    h_hg_rms->Draw("COLZ");
    h_hg_rms->SetStats(1);
    std::string png_hg_rms = std::string(output_dir) + "/hg_rms_comparison.png";
    c_hg_rms->SaveAs(png_hg_rms.c_str());
    std::cout << "  ✓ Saved: " << png_hg_rms << std::endl;
    
    // LG Sigma plot
    TCanvas* c_lg_sigma = new TCanvas("c_lg_sigma", "LG Sigma Comparison", 800, 600);
    h_lg_sigma->Draw("COLZ");
    h_lg_sigma->SetStats(1);
    std::string png_lg_sigma = std::string(output_dir) + "/lg_sigma_comparison.png";
    c_lg_sigma->SaveAs(png_lg_sigma.c_str());
    std::cout << "  ✓ Saved: " << png_lg_sigma << std::endl;
    
    // LG RMS plot
    TCanvas* c_lg_rms = new TCanvas("c_lg_rms", "LG RMS Comparison", 800, 600);
    h_lg_rms->Draw("COLZ");
    h_lg_rms->SetStats(1);
    std::string png_lg_rms = std::string(output_dir) + "/lg_rms_comparison.png";
    c_lg_rms->SaveAs(png_lg_rms.c_str());
    std::cout << "  ✓ Saved: " << png_lg_rms << std::endl;
    
    // ========================================================================
    // Phase 7: Print Summary Statistics
    // ========================================================================
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Summary Statistics (by variable)" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    auto print_stats = [](const char* name, const ComparisonStats& stats) {
        std::cout << name << ":" << std::endl;
        std::cout << "  Matching channels: " << stats.n_matches << std::endl;
        std::cout << "  Pedestal  → Mean: " << stats.ped_mean << " ± " << stats.ped_rms 
                  << " (Range: " << stats.ped_min << " - " << stats.ped_max << ")" << std::endl;
        std::cout << "  Intercalib → Mean: " << stats.ic_mean << " ± " << stats.ic_rms 
                  << " (Range: " << stats.ic_min << " - " << stats.ic_max << ")" << std::endl;
        std::cout << "  Correlation: " << stats.correlation << std::endl;
        std::cout << std::endl;
    };
    
    print_stats("High Gain Sigma", stat_hg_sigma);
    print_stats("High Gain RMS", stat_hg_rms);
    print_stats("Low Gain Sigma", stat_lg_sigma);
    print_stats("Low Gain RMS", stat_lg_rms);
    
    std::cout << "========================================" << std::endl;
    std::cout << "Analysis Complete!" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // Cleanup
    ped_file->Close();
    ic_file->Close();
    delete ped_file;
    delete ic_file;
}

// ============================================================================
// Helper function for command-line execution
// ============================================================================

#if !defined(__CINT__) && !defined(__CLING__)
int main(int argc, char** argv) {
    const char* ped_file = argc > 1 ? argv[1] : 
        "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/calibration/pedestal_calib/run21722/pedestal.root";
    const char* ic_file = argc > 2 ? argv[2] :
        "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out_calibration/intercalib_calib/run22133/intercalib.root";
    const char* out_dir = argc > 3 ? argv[3] : ".";
    
    compare_calibration(ped_file, ic_file, out_dir);
    return 0;
}
#endif
