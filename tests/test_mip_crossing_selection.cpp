#include "calibration/module/MIP/MIPAlg.hpp"
#include "common/AHCALGeometry.hpp"
#include "common/EventStore.hpp"
#include "common/RunContext.hpp"
#include "common/edm/EDM.hpp"

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Result {
    long long selected_crossing = 0;
    long long out_of_active_area = 0;
    long long outside_xy_threshold = 0;
    long long failed_neighbor_support = 0;
    long long neighbor_cell_hit = 0;
    long long target_hit_fill = 0;
    long long no_exact_target_hit = 0;
    int target_ntracks = 0;
    double target_entries = 0.0;
};

void expect(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

AHCALRawHit makeHit(int layer, int chip, int channel) {
    AHCALRawHit hit{};
    hit.cellID = AHCALGeometry::CellID(layer, chip, channel);
    hit.hg_adc = 250;
    hit.lg_adc = 0;
    hit.hittag = 1;
    return hit;
}

SimpleFittedTrack makeStraightTrack(int chip, int channel) {
    SimpleFittedTrack track{};
    track.init_pos_x = AHCALGeometry::Pos_X(channel, chip);
    track.init_pos_y = AHCALGeometry::Pos_Y(channel, chip);
    track.direction_x = 0.0;
    track.direction_y = 0.0;
    track.valid = true;
    return track;
}

Result runCase(const std::string& name, int target_layer,
               const std::vector<AHCALRawHit>& rawhits,
               const std::string& extra_cfg = "", double track_x_offset = 0.0) {
    constexpr int chip = 4;
    constexpr int channel = 7;
    const std::filesystem::path output =
        std::filesystem::temp_directory_path() / ("fair_" + name + ".root");
    std::filesystem::remove(output);

    RunContext ctx;
    ctx.config.runNumber = 1;
    {
        AHCALRecoAlg::MIPAlg alg(ctx, "MIPAlgTest");
        YAML::Node cfg = YAML::Load(
            "in_rawhit_key: RawHits\n"
            "in_track_key: FittedTrack\n"
            "string_track_struct: SimpleFittedTrack\n"
            "track_selection_string: ''\n"
            "mip_to_file: true\n"
            "out_mip_filename: " + output.string() + "\n"
            "fit: false\n"
            "mip_to_json: false\n"
            "substrate_pedestal: false\n"
            "read_pedestal_from_ROOT: false\n"
            "read_pedestal_from_DB: false\n"
            "xy_size_threshold: 1.0\n"
            "enable_neighbor_layer_crossing_selection: true\n"
            "mip_neighbor_upstream_layers: 2\n"
            "mip_neighbor_downstream_layers: 2\n"
            "mip_edge_layer_policy: true\n"
            "mip_reject_if_neighbor_cell_hit: true\n"
            "mip_neighbor_cell_search_radius: 1\n"
            "mip_skip_layers: []\n"
            "mip_skip_chips: []\n");
        if (!extra_cfg.empty()) {
            const YAML::Node overrides = YAML::Load(extra_cfg);
            for (const auto& entry : overrides) {
                cfg[entry.first.as<std::string>()] = entry.second;
            }
        }
        alg.parse_cfg(cfg);

        EventStore evt;
        evt.put("RawHits", rawhits);
        auto track = makeStraightTrack(chip, channel);
        track.init_pos_x += track_x_offset;
        evt.put("FittedTrack", track);
        alg.execute(evt);
    }

    TFile file(output.c_str(), "READ");
    expect(!file.IsZombie(), "failed to open test ROOT output");

    Result result;
    auto* cutflow = file.Get<TTree>("MIP/mip_crossing_cutflow");
    expect(cutflow != nullptr, "missing crossing cutflow tree");
    int layer = -1;
    cutflow->SetBranchAddress("layer", &layer);
    cutflow->SetBranchAddress("selected_crossing", &result.selected_crossing);
    cutflow->SetBranchAddress("out_of_active_area", &result.out_of_active_area);
    cutflow->SetBranchAddress("outside_xy_threshold", &result.outside_xy_threshold);
    cutflow->SetBranchAddress("failed_neighbor_support", &result.failed_neighbor_support);
    cutflow->SetBranchAddress("neighbor_cell_hit", &result.neighbor_cell_hit);
    cutflow->SetBranchAddress("target_hit_fill", &result.target_hit_fill);
    cutflow->SetBranchAddress("no_exact_target_hit", &result.no_exact_target_hit);
    for (Long64_t entry = 0; entry < cutflow->GetEntries(); ++entry) {
        cutflow->GetEntry(entry);
        if (layer == target_layer) break;
    }

    auto* ntrack = file.Get<TTree>("MIP/ntrack_pass_through_channel");
    expect(ntrack != nullptr, "missing all-cell ntrack tree");
    expect(ntrack->GetEntries() == AHCALGeometry::Layer_No * AHCALGeometry::chip_No *
                                      AHCALGeometry::channel_No,
           "ntrack tree does not contain every cell");
    int cellid = -1;
    int ntracks = 0;
    ntrack->SetBranchAddress("cellid", &cellid);
    ntrack->SetBranchAddress("ntracks", &ntracks);
    const int target_cellid = AHCALGeometry::CellID(target_layer, chip, channel);
    for (Long64_t entry = 0; entry < ntrack->GetEntries(); ++entry) {
        ntrack->GetEntry(entry);
        if (cellid == target_cellid) {
            result.target_ntracks = ntracks;
            break;
        }
    }

    auto* hist = file.Get<TH1D>(
        ("MIP/Layer" + (target_layer < 10 ? std::string("0") : std::string()) +
         std::to_string(target_layer) + "/Chip4/hMIP_" + std::to_string(target_cellid)).c_str());
    result.target_entries = hist ? hist->GetEntries() : 0.0;
    file.Close();
    std::filesystem::remove(output);
    return result;
}

std::vector<AHCALRawHit> hitsAt(const std::vector<int>& layers) {
    std::vector<AHCALRawHit> hits;
    for (int layer : layers) {
        hits.push_back(makeHit(layer, 4, 7));
    }
    return hits;
}

void testCentralAndMissingTargetHit() {
    const Result result = runCase("mip_central", 10, hitsAt({8, 9, 11, 12}));
    expect(result.selected_crossing == 1, "central crossing was not selected");
    expect(result.target_ntracks == 1, "missing target hit did not increment denominator");
    expect(result.no_exact_target_hit == 1, "missing target hit was not diagnosed");
    expect(result.target_entries == 0.0, "missing target hit filled the MIP histogram");
}

void testEdgePolicy() {
    const Result result = runCase("mip_edge", 0, hitsAt({1, 2, 3, 4}));
    expect(result.selected_crossing == 1, "edge layer 0 did not require four downstream hits");

    const Result skipped_upstream = runCase(
        "mip_edge_skipped_upstream", 1, hitsAt({3, 4, 5, 6}),
        "mip_skip_layers: [0, 2, 14, 28]\n");
    expect(skipped_upstream.selected_crossing == 1,
           "skipped layer 0 was not transferred to downstream support for layer 1");

    const Result missing_transferred_support = runCase(
        "mip_edge_missing_transferred_support", 1, hitsAt({3, 4, 5}),
        "mip_skip_layers: [0, 2, 14, 28]\n");
    expect(missing_transferred_support.failed_neighbor_support == 1,
           "layer 1 did not require four downstream hits after skipping layer 0");
}

void testSkippedLayersAndChips() {
    Result result = runCase("mip_skip_layer", 10, hitsAt({7, 8, 11, 12}),
                            "mip_skip_layers: [9]\n");
    expect(result.selected_crossing == 1, "skipped layer was not replaced by next valid layer");

    result = runCase("mip_skip_chip", 10, hitsAt({7, 8, 11, 12}),
                     "mip_skip_chips: [[9, 4]]\n");
    expect(result.selected_crossing == 1, "skipped chip was not replaced by next valid layer");
}

void testRejectionsAndTargetFill() {
    Result result = runCase("mip_support_fail", 10, hitsAt({8, 9, 11}));
    expect(result.failed_neighbor_support == 1, "missing support hit was not rejected");
    expect(result.target_ntracks == 0, "rejected crossing incremented denominator");

    auto neighbor_hits = hitsAt({8, 9, 11, 12});
    neighbor_hits.push_back(makeHit(10, 4, 8));
    result = runCase("mip_neighbor", 10, neighbor_hits);
    expect(result.neighbor_cell_hit == 1, "neighbor cell hit was not rejected");

    result = runCase("mip_target_fill", 10, hitsAt({8, 9, 10, 11, 12}));
    expect(result.target_hit_fill == 1, "exact target hit was not filled");
    expect(result.target_entries == 1.0, "target MIP histogram has wrong entries");
}

void testGeometryRejections() {
    Result result = runCase("mip_xy_threshold", 10, {}, "", 5.0);
    expect(result.outside_xy_threshold == 1, "xy threshold rejection was not diagnosed");

    const double out_of_area_offset =
        AHCALGeometry::x_max - AHCALGeometry::Pos_X(7, 4) + 1.0;
    result = runCase("mip_out_of_area", 10, {}, "", out_of_area_offset);
    expect(result.out_of_active_area == 1, "active-area rejection was not diagnosed");
}

void testLegacySelectionSwitches() {
    const Result result = runCase(
        "mip_legacy", 10, {},
        "enable_neighbor_layer_crossing_selection: false\n"
        "mip_reject_if_neighbor_cell_hit: false\n");
    expect(result.selected_crossing == 1, "legacy switches did not preserve geometric crossing");
    expect(result.target_ntracks == 1, "legacy crossing did not increment denominator");
}

}  // namespace

int main() {
    testCentralAndMissingTargetHit();
    testEdgePolicy();
    testSkippedLayersAndChips();
    testRejectionsAndTargetFill();
    testGeometryRejections();
    testLegacySelectionSwitches();
    return 0;
}
