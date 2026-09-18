#include "MIPSimAlg.hpp"

#include "common/AlgRegistry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/edm/RawHit.hpp"
#include "common/edm/SimHit.hpp"

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TParameter.h>

#include <filesystem>
#include <stdexcept>
#include <vector>

AHCAL_REGISTER_ALG(AHCALRecoAlg::MIPSimAlg, "MIPSimAlg")

namespace AHCALRecoAlg {

MIPSimAlg::MIPSimAlg(RunContext& rc, std::string name)
    : IAlg(rc, name), mip_(rc, name + "/MIPAlg"),
      selected_key_(name + "::selected_rawhit_indices") {}

MIPSimAlg::~MIPSimAlg() { write_truth(); }

void MIPSimAlg::parse_cfg(const YAML::Node& cfg) {
    in_simhit_key_ = get_or<std::string>(cfg, "in_simhit_key", in_simhit_key_);
    in_rawhit_key_ = get_or<std::string>(cfg, "in_rawhit_key", in_rawhit_key_);
    out_simhit_filename_ = get_or<std::string>(cfg, "out_simhit_filename", out_simhit_filename_);
    edep_nbin_ = get_or<int>(cfg, "edep_nbin", edep_nbin_);
    edep_max_ = get_or<double>(cfg, "edep_max", edep_max_);
    if (edep_nbin_ <= 0 || edep_max_ <= 0.0 ||
        out_simhit_filename_ == get_or<std::string>(cfg, "out_mip_filename", "mip_sim_adc.root")) {
        throw std::runtime_error("MIPSimAlg: invalid truth histogram configuration or output filename collision");
    }

    // Every MIPAlg selection and fitting parameter is read from this same cfg.
    YAML::Node mip_cfg = YAML::Clone(cfg);
    mip_cfg["in_rawhit_key"] = in_rawhit_key_;
    if (!mip_cfg["in_track_key"]) mip_cfg["in_track_key"] = "FittedTrack";
    if (!mip_cfg["out_mip_filename"]) mip_cfg["out_mip_filename"] = "mip_sim_adc.root";
    mip_cfg["out_selected_rawhit_indices_key"] = selected_key_;
    mip_.parse_cfg(mip_cfg);
}

void MIPSimAlg::init_by_run() { mip_.init_by_run(); }

void MIPSimAlg::execute(EventStore& evt) {
    const auto& simhits = evt.get<std::vector<AHCALSimHit>>(in_simhit_key_);
    // Check both inputs before MIPAlg publishes its selection into the event.
    const auto& rawhits = evt.get<std::vector<AHCALRawHit>>(in_rawhit_key_);
    mip_.execute(evt);
    const auto& selected = evt.get<std::vector<int>>(selected_key_);

    // SimHitReader normally gives one entry per cell; summing also handles
    // an input that contains multiple deposits in the same cell.
    std::unordered_map<int, double> edep_by_cell;
    edep_by_cell.reserve(simhits.size());
    for (const auto& hit : simhits) edep_by_cell[hit.cellID] += hit.Edep;

    for (int raw_index : selected) {
        const int cellid = rawhits.at(raw_index).cellID;
        ++n_selected_;
        const auto it = edep_by_cell.find(cellid);
        if (it == edep_by_cell.end()) {
            ++n_missing_simhit_; // possible for added electronic noise
            continue;
        }
        auto& hist = edep_hist_[cellid];
        if (!hist) {
            const auto hist_name = "hSimEdep_" + std::to_string(cellid);
            hist = std::make_unique<TH1D>(hist_name.c_str(),
                "Selected visible SimHit energy;E_{vis} [MeV];counts",
                edep_nbin_, 0.0, edep_max_);
            hist->SetDirectory(nullptr);
        }
        hist->Fill(it->second);
    }
}

void MIPSimAlg::write_truth() {
    const std::filesystem::path path(out_simhit_filename_);
    if (path.has_parent_path()) {
        std::error_code ec;
        std::filesystem::create_directories(path.parent_path(), ec);
        if (ec) {
            LOG_ERROR("MIPSimAlg: cannot create directory {}: {}", path.parent_path().string(), ec.message());
            return;
        }
    }
    std::unique_ptr<TFile> output(TFile::Open(out_simhit_filename_.c_str(), "RECREATE"));
    if (!output || output->IsZombie()) {
        LOG_ERROR("MIPSimAlg: cannot write {}", out_simhit_filename_);
        return;
    }
    TParameter<long long>("selected_rawhits", n_selected_).Write();
    TParameter<long long>("selected_without_simhit", n_missing_simhit_).Write();
    for (const auto& [cellid, hist] : edep_hist_) {
        const auto layer = "Layer" + std::to_string(cellid / 100000);
        auto* dir = output->GetDirectory(layer.c_str());
        if (!dir) dir = output->mkdir(layer.c_str());
        dir->cd();
        hist->Write();
        output->cd();
    }
    output->Close();
    LOG_INFO("MIPSimAlg: wrote {} truth spectra; {} of {} selected hits had no SimHit",
             edep_hist_.size(), n_missing_simhit_, n_selected_);
}

} // namespace AHCALRecoAlg
