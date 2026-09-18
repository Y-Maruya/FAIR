#include "calibration/module/MIPSim/MIPSimAlg.hpp"
#include "common/AHCALGeometry.hpp"
#include "common/EventStore.hpp"
#include "common/RunContext.hpp"
#include "common/edm/EDM.hpp"

#include <TFile.h>
#include <TH1D.h>
#include <TParameter.h>

#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void expect(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

void check_case(bool add_neighbor, bool reject_track) {
    const auto temp = std::filesystem::temp_directory_path();
    const auto suffix = add_neighbor ? "neighbor" : (reject_track ? "track" : "accepted");
    const auto adc_file = temp / (std::string("fair_mipsim_adc_") + suffix + ".root");
    const auto truth_file = temp / (std::string("fair_mipsim_truth_") + suffix + ".root");
    std::filesystem::remove(adc_file);
    std::filesystem::remove(truth_file);

    constexpr int chip = 4;
    constexpr int channel = 7;
    const int target = AHCALGeometry::CellID(10, chip, channel);
    RunContext ctx;
    ctx.config.runNumber = 1;
    {
        AHCALRecoAlg::MIPSimAlg alg(ctx, "MIPSimAlgTest");
        const auto cfg = YAML::Load(
            "in_simhit_key: SimHits\n"
            "in_rawhit_key: SimRawHits\n"
            "in_track_key: FittedTrack\n"
            "track_selection_string: valid>0\n"
            "mip_to_file: true\n"
            "out_mip_filename: " + adc_file.string() + "\n"
            "out_simhit_filename: " + truth_file.string() + "\n"
            "fit: false\n"
            "substrate_pedestal: false\n"
            "read_pedestal_from_DB: false\n"
            "xy_size_threshold: 1\n"
            "mip_neighbor_upstream_layers: 2\n"
            "mip_neighbor_downstream_layers: 2\n"
            "mip_reject_if_neighbor_cell_hit: true\n");
        alg.parse_cfg(cfg);
        EventStore event;
        std::vector<AHCALRawHit> raws;
        std::vector<AHCALSimHit> sims;
        for (int layer : {8, 9, 10, 11, 12}) {
            AHCALRawHit raw{};
            raw.cellID = AHCALGeometry::CellID(layer, chip, channel);
            raw.hg_adc = 250;
            raw.hittag = 1;
            raws.push_back(raw);
            AHCALSimHit sim{};
            sim.cellID = raw.cellID;
            sim.Edep = 0.461;
            sims.push_back(sim);
        }
        if (add_neighbor) {
            AHCALRawHit neighbor{};
            neighbor.cellID = AHCALGeometry::CellID(10, chip, channel + 1);
            neighbor.hittag = 1;
            raws.push_back(neighbor);
        }
        event.put("SimRawHits", std::move(raws));
        event.put("SimHits", std::move(sims));
        SimpleFittedTrack track{};
        track.init_pos_x = AHCALGeometry::Pos_X(channel, chip);
        track.init_pos_y = AHCALGeometry::Pos_Y(channel, chip);
        track.valid = !reject_track;
        event.put("FittedTrack", track);
        alg.execute(event);
    }
    TFile adc(adc_file.c_str(), "READ");
    TFile truth(truth_file.c_str(), "READ");
    expect(!adc.IsZombie() && !truth.IsZombie(), "MIPSim output missing");
    auto* adc_hist = adc.Get<TH1D>(("MIP/Layer10/Chip4/hMIP_" + std::to_string(target)).c_str());
    auto* truth_hist = truth.Get<TH1D>(("Layer10/hSimEdep_" + std::to_string(target)).c_str());
    const bool accepted = !add_neighbor && !reject_track;
    expect((adc_hist != nullptr) == accepted, "ADC selection mismatch");
    expect((truth_hist != nullptr) == accepted, "SimHit selection mismatch");
    if (accepted) {
        expect(adc_hist->GetEntries() == 1 && truth_hist->GetEntries() == 1,
               "selected hit did not fill both histograms");
        expect(truth_hist->GetMean() > 0.45 && truth_hist->GetMean() < 0.48,
               "visible energy not propagated");
    }
    auto* selected = truth.Get<TParameter<long long>>("selected_rawhits");
    expect(selected && selected->GetVal() >= (accepted ? 1 : 0), "selection summary missing");
    adc.Close();
    truth.Close();
    std::filesystem::remove(adc_file);
    std::filesystem::remove(truth_file);
}
} // namespace

int main() {
    check_case(false, false);
    check_case(true, false);
    check_case(false, true);
}
