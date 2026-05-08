#include "ClusterBasicAnaAlg.hpp"
#include "common/AlgRegistry.hpp"
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <filesystem>

AHCAL_REGISTER_ALG(AHCALRecoAlg::ClusterBasicAnaAlg, "ClusterBasicAnaAlg")

namespace AHCALRecoAlg {

struct ClusterBasicAnaAlg::Impl {
    std::unique_ptr<TFile> outFile;
    TH1D* hNcluster = nullptr;
    TH1D* hTotalEdep = nullptr;
    TH2D* hNclusterVsTotalEdep = nullptr;
    TH1D* hNclusterAboveThreshold = nullptr;
    TH2D* hThresholdVsNcluster = nullptr;
};

void ClusterBasicAnaAlg::ImplDeleter::operator()(Impl* p) const { delete p; }
ClusterBasicAnaAlg::~ClusterBasicAnaAlg() = default;

void ClusterBasicAnaAlg::initialize() {
    impl_.reset(new Impl());
    impl_->hNcluster = new TH1D("hNcluster", "Ncluster", 100, 0, 100);
    impl_->hTotalEdep = new TH1D("hTotalEdep", "TotalEdep", 200, 0, 2000);
    impl_->hNclusterVsTotalEdep = new TH2D("hNclusterVsTotalEdep", "Ncluster vs TotalEdep", 100, 0, 100, 200, 0, 2000);
    impl_->hNclusterAboveThreshold = new TH1D("hNclusterAboveThreshold", "Ncluster(Edep>threshold)", 100, 0, 100);
    impl_->hThresholdVsNcluster = new TH2D("hThresholdVsNcluster", "threshold vs Ncluster;threshold [MeV];Ncluster", 200, 0, 50, 100, 0, 100);

    const auto parent = std::filesystem::path(cfg_.out_root_filename).parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    impl_->outFile.reset(TFile::Open(cfg_.out_root_filename.c_str(), "RECREATE"));
}

void ClusterBasicAnaAlg::execute(EventStore& evt) {
    auto& clusters = evt.get<std::vector<TutorialCluster>>(cfg_.in_cluster_key);
    const int ncluster = static_cast<int>(clusters.size());

    double totalEdepEvent = 0.0;
    int nclusterAboveThreshold = 0;
    for (const auto& c : clusters) {
        totalEdepEvent += c.totalEdep;
        if (c.totalEdep >= cfg_.edep_threshold) ++nclusterAboveThreshold;
    }

    impl_->hNcluster->Fill(ncluster);
    impl_->hTotalEdep->Fill(totalEdepEvent);
    impl_->hNclusterVsTotalEdep->Fill(ncluster, totalEdepEvent);
    impl_->hNclusterAboveThreshold->Fill(nclusterAboveThreshold);

    for (double thr = 0.0; thr <= cfg_.threshold_scan_max; thr += cfg_.threshold_scan_step) {
        int nAbove = 0;
        for (const auto& c : clusters) {
            if (c.totalEdep >= thr) ++nAbove;
        }
        impl_->hThresholdVsNcluster->Fill(thr, nAbove);
    }
}

void ClusterBasicAnaAlg::finalize() {
    if (impl_ && impl_->outFile) {
        impl_->outFile->cd();
        impl_->hNcluster->Write();
        impl_->hTotalEdep->Write();
        impl_->hNclusterVsTotalEdep->Write();
        impl_->hNclusterAboveThreshold->Write();
        impl_->hThresholdVsNcluster->Write();
        impl_->outFile->Write();
        impl_->outFile->Close();
    }

    if (cfg_.write_to_png && impl_) {
        std::filesystem::create_directories(cfg_.out_png_dir);
        impl_->hNcluster->SaveAs((cfg_.out_png_dir + "/hNcluster.png").c_str());
        impl_->hTotalEdep->SaveAs((cfg_.out_png_dir + "/hTotalEdep.png").c_str());
        impl_->hNclusterVsTotalEdep->SaveAs((cfg_.out_png_dir + "/hNclusterVsTotalEdep.png").c_str());
        impl_->hNclusterAboveThreshold->SaveAs((cfg_.out_png_dir + "/hNclusterAboveThreshold.png").c_str());
        impl_->hThresholdVsNcluster->SaveAs((cfg_.out_png_dir + "/hThresholdVsNcluster.png").c_str());
    }
}

void ClusterBasicAnaAlg::parse_cfg(const YAML::Node& n) {
    cfg_.in_cluster_key = n["in_cluster_key"].as<std::string>(cfg_.in_cluster_key);
    cfg_.edep_threshold = n["edep_threshold"].as<double>(cfg_.edep_threshold);
    cfg_.threshold_scan_max = n["threshold_scan_max"].as<double>(cfg_.threshold_scan_max);
    cfg_.threshold_scan_step = n["threshold_scan_step"].as<double>(cfg_.threshold_scan_step);
    cfg_.out_root_filename = n["out_root_filename"].as<std::string>(cfg_.out_root_filename);
    cfg_.write_to_png = n["write_to_png"].as<bool>(cfg_.write_to_png);
    cfg_.out_png_dir = n["out_png_dir"].as<std::string>(cfg_.out_png_dir);
}

} // namespace AHCALRecoAlg
