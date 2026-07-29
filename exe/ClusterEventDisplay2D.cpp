#include <TCanvas.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TString.h>
#include <TPad.h>
#include <TLine.h>
#include <TMarker.h>
#include <TBox.h>
#include <TROOT.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "common/RunContext.hpp"
#include "common/config/ParseRunConfig.hpp"
#include "common/edm/Cluster.hpp"
#include "common/edm/EDM.hpp"
#include "IO/reader/ReaderRegistry.hpp"
#include "IO/reader/RootInput.hpp"

namespace {

const double xMin = -AHCALGeometry::x_max;
const double xMax = +AHCALGeometry::x_max;
const double yMin = -AHCALGeometry::y_max;
const double yMax = +AHCALGeometry::y_max;
const double zMin = -9.0;
const double zMax = 1191.0 + 15.0;
const double xyCellSize = AHCALGeometry::x_max * 2.0 / 18.0;
const double zCellSize = (1191.0 + 24.0) / 135.0;

struct Options {
    std::string input;
    Long64_t event_id = -1;
    std::string output;
    std::string cluster_key = "PCAClusters";
    std::string hit_key = "RecoHits";
    std::string tlu_key = "TLURawData";
    std::string tree_name = "events";
    bool by_tlu_eventid = false;
};

bool ends_with(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool is_yaml_file(const std::string& path) {
    return ends_with(path, ".yaml") || ends_with(path, ".yml");
}

void usage(const char* argv0) {
    std::cerr
        << "Usage: " << argv0 << " <config.yaml|input.root> <event_id> [output.png]\n"
        << "       [--cluster-key PCAClusters] [--hit-key RecoHits]\n"
        << "       [--tlu-key TLURawData] [--tree events] [--by-tlu-eventid]\n\n"
        << "  event_id is the ROOT TTree entry index by default.\n"
        << "  With --by-tlu-eventid, event_id is matched to TLURawData.EventID.\n";
}

Options parse_args(int argc, char** argv) {
    if (argc < 3) {
        usage(argv[0]);
        throw std::runtime_error("missing required arguments");
    }

    Options opt;
    opt.input = argv[1];
    opt.event_id = std::stoll(argv[2]);

    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const std::string& name) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error("missing value for " + name);
            }
            return argv[++i];
        };

        if (arg == "--cluster-key") {
            opt.cluster_key = require_value(arg);
        } else if (arg.rfind("--cluster-key=", 0) == 0) {
            opt.cluster_key = arg.substr(std::string("--cluster-key=").size());
        } else if (arg == "--hit-key") {
            opt.hit_key = require_value(arg);
        } else if (arg.rfind("--hit-key=", 0) == 0) {
            opt.hit_key = arg.substr(std::string("--hit-key=").size());
        } else if (arg == "--tlu-key") {
            opt.tlu_key = require_value(arg);
        } else if (arg.rfind("--tlu-key=", 0) == 0) {
            opt.tlu_key = arg.substr(std::string("--tlu-key=").size());
        } else if (arg == "--tree") {
            opt.tree_name = require_value(arg);
        } else if (arg.rfind("--tree=", 0) == 0) {
            opt.tree_name = arg.substr(std::string("--tree=").size());
        } else if (arg == "--by-tlu-eventid") {
            opt.by_tlu_eventid = true;
        } else if (arg.rfind("-", 0) == 0) {
            throw std::runtime_error("unknown option: " + arg);
        } else if (opt.output.empty()) {
            opt.output = arg;
        } else {
            throw std::runtime_error("unexpected argument: " + arg);
        }
    }

    if (opt.output.empty()) {
        opt.output = Form("eventdisplay_cluster/evt%lld_%s.png",
                          static_cast<long long>(opt.event_id),
                          opt.cluster_key.c_str());
    }

    return opt;
}

std::string resolve_input_root(const std::string& input) {
    if (!is_yaml_file(input)) return input;

    YAML::Node config = YAML::LoadFile(input);
    RunContext ctx;
    ctx.config = parse_run_config(config);
    return ctx.config.output;
}

int cluster_color(int i) {
    static const int colors[] = {
        kRed + 1, kAzure + 2, kGreen + 2, kMagenta + 1, kOrange + 7,
        kCyan + 2, kViolet + 7, kPink + 6, kSpring + 5, kBlue + 1,
        kTeal + 3, kOrange + 2
    };
    return colors[i % (sizeof(colors) / sizeof(colors[0]))];
}

bool has_pca_position(const AHCALRecoAlg::Cluster& c) {
    return c.nHitsInCluster >= 6 &&
           c.pca_lambda1 > 0.0 &&
           std::isfinite(c.pca_mean_x) &&
           std::isfinite(c.pca_mean_y) &&
           std::isfinite(c.pca_mean_z);
}

double pca_draw_length(const AHCALRecoAlg::Cluster& c) {
    return std::clamp(2.0 * std::sqrt(std::max(c.pca_lambda1, 0.0)), 35.0, 220.0);
}

using ClusterHitRefs = std::vector<std::vector<const AHCALRecoHit*>>;

ClusterHitRefs make_cluster_hit_refs(const std::vector<AHCALRecoHit>& hits,
                                     const std::vector<AHCALRecoAlg::Cluster>& clusters) {
    std::unordered_map<int, const AHCALRecoHit*> hit_by_index;
    hit_by_index.reserve(hits.size());
    for (const auto& hit : hits) {
        hit_by_index.emplace(hit.index, &hit);
    }

    ClusterHitRefs refs(clusters.size());
    for (std::size_t ic = 0; ic < clusters.size(); ++ic) {
        const auto& c = clusters[ic];
        auto& out = refs[ic];

        if (!c.hits.empty()) {
            out.reserve(c.hits.size());
            for (const auto& h : c.hits) out.push_back(&h);
            continue;
        }

        out.reserve(c.index.size());
        for (const int idx : c.index) {
            const auto it = hit_by_index.find(idx);
            if (it != hit_by_index.end()) {
                out.push_back(it->second);
            }
        }
    }

    return refs;
}

std::size_t count_cluster_hit_refs(const ClusterHitRefs& refs) {
    std::size_t n = 0;
    for (const auto& v : refs) n += v.size();
    return n;
}

void draw_box_xy(const AHCALRecoHit& hit, int color) {
    const double x = hit.Xpos();
    const double y = hit.Ypos();
    auto* box = new TBox(x - xyCellSize / 2.0, y - xyCellSize / 2.0,
                         x + xyCellSize / 2.0, y + xyCellSize / 2.0);
    box->SetLineColor(color);
    box->SetLineWidth(2);
    box->SetFillStyle(0);
    box->Draw("SAME");
}

void draw_box_xz(const AHCALRecoHit& hit, int color) {
    const double x = hit.Xpos();
    const double z = hit.Zpos();
    auto* box = new TBox(z - zCellSize / 2.0, x - xyCellSize / 2.0,
                         z + zCellSize / 2.0, x + xyCellSize / 2.0);
    box->SetLineColor(color);
    box->SetLineWidth(2);
    box->SetFillStyle(0);
    box->Draw("SAME");
}

void draw_box_yz(const AHCALRecoHit& hit, int color) {
    const double y = hit.Ypos();
    const double z = hit.Zpos();
    auto* box = new TBox(z - zCellSize / 2.0, y - xyCellSize / 2.0,
                         z + zCellSize / 2.0, y + xyCellSize / 2.0);
    box->SetLineColor(color);
    box->SetLineWidth(2);
    box->SetFillStyle(0);
    box->Draw("SAME");
}

void draw_cluster_hit_boxes_xy(const ClusterHitRefs& refs) {
    for (std::size_t ic = 0; ic < refs.size(); ++ic) {
        const int color = cluster_color(static_cast<int>(ic));
        for (const auto* hit : refs[ic]) {
            if (hit) draw_box_xy(*hit, color);
        }
    }
}

void draw_cluster_hit_boxes_xz(const ClusterHitRefs& refs) {
    for (std::size_t ic = 0; ic < refs.size(); ++ic) {
        const int color = cluster_color(static_cast<int>(ic));
        for (const auto* hit : refs[ic]) {
            if (hit) draw_box_xz(*hit, color);
        }
    }
}

void draw_cluster_hit_boxes_yz(const ClusterHitRefs& refs) {
    for (std::size_t ic = 0; ic < refs.size(); ++ic) {
        const int color = cluster_color(static_cast<int>(ic));
        for (const auto* hit : refs[ic]) {
            if (hit) draw_box_yz(*hit, color);
        }
    }
}

void draw_cluster_overlay_xy(const std::vector<AHCALRecoAlg::Cluster>& clusters) {
    TLatex label;
    label.SetTextSize(0.022);
    label.SetTextColor(kBlack);

    int drawn = 0;
    for (std::size_t i = 0; i < clusters.size(); ++i) {
        const auto& c = clusters[i];
        if (!has_pca_position(c)) continue;

        const int color = cluster_color(static_cast<int>(i));
        const double len = pca_draw_length(c);
        const double x1 = c.pca_mean_x - c.pca_axis_x * len;
        const double y1 = c.pca_mean_y - c.pca_axis_y * len;
        const double x2 = c.pca_mean_x + c.pca_axis_x * len;
        const double y2 = c.pca_mean_y + c.pca_axis_y * len;

        auto* line = new TLine(x1, y1, x2, y2);
        line->SetLineColor(color);
        line->SetLineWidth(3);
        line->Draw("SAME");

        auto* marker = new TMarker(c.pca_mean_x, c.pca_mean_y, 29);
        marker->SetMarkerColor(color);
        marker->SetMarkerSize(1.6);
        marker->Draw("SAME");

        if (drawn < 20) {
            label.SetTextColor(color);
            label.DrawLatex(c.pca_mean_x + 5.0, c.pca_mean_y + 5.0,
                            Form("#%d", c.cluster_id));
        }
        ++drawn;
    }
}

void draw_cluster_overlay_xz(const std::vector<AHCALRecoAlg::Cluster>& clusters) {
    int drawn = 0;
    TLatex label;
    label.SetTextSize(0.022);

    for (std::size_t i = 0; i < clusters.size(); ++i) {
        const auto& c = clusters[i];
        if (!has_pca_position(c)) continue;

        const int color = cluster_color(static_cast<int>(i));
        const double len = pca_draw_length(c);
        const double z1 = c.pca_mean_z - c.pca_axis_z * len;
        const double x1 = c.pca_mean_x - c.pca_axis_x * len;
        const double z2 = c.pca_mean_z + c.pca_axis_z * len;
        const double x2 = c.pca_mean_x + c.pca_axis_x * len;

        auto* line = new TLine(z1, x1, z2, x2);
        line->SetLineColor(color);
        line->SetLineWidth(3);
        line->Draw("SAME");

        auto* marker = new TMarker(c.pca_mean_z, c.pca_mean_x, 29);
        marker->SetMarkerColor(color);
        marker->SetMarkerSize(1.6);
        marker->Draw("SAME");

        if (drawn < 20) {
            label.SetTextColor(color);
            label.DrawLatex(c.pca_mean_z + 8.0, c.pca_mean_x + 5.0,
                            Form("#%d", c.cluster_id));
        }
        ++drawn;
    }
}

void draw_cluster_overlay_yz(const std::vector<AHCALRecoAlg::Cluster>& clusters) {
    int drawn = 0;
    TLatex label;
    label.SetTextSize(0.022);

    for (std::size_t i = 0; i < clusters.size(); ++i) {
        const auto& c = clusters[i];
        if (!has_pca_position(c)) continue;

        const int color = cluster_color(static_cast<int>(i));
        const double len = pca_draw_length(c);
        const double z1 = c.pca_mean_z - c.pca_axis_z * len;
        const double y1 = c.pca_mean_y - c.pca_axis_y * len;
        const double z2 = c.pca_mean_z + c.pca_axis_z * len;
        const double y2 = c.pca_mean_y + c.pca_axis_y * len;

        auto* line = new TLine(z1, y1, z2, y2);
        line->SetLineColor(color);
        line->SetLineWidth(3);
        line->Draw("SAME");

        auto* marker = new TMarker(c.pca_mean_z, c.pca_mean_y, 29);
        marker->SetMarkerColor(color);
        marker->SetMarkerSize(1.6);
        marker->Draw("SAME");

        if (drawn < 20) {
            label.SetTextColor(color);
            label.DrawLatex(c.pca_mean_z + 8.0, c.pca_mean_y + 5.0,
                            Form("#%d", c.cluster_id));
        }
        ++drawn;
    }
}

void draw_panel(Long64_t entry,
                Long64_t requested_event_id,
                bool by_tlu_eventid,
                const std::string& input_root,
                const std::string& cluster_key,
                const std::vector<AHCALRecoHit>& hits,
                const std::vector<AHCALRecoAlg::Cluster>& clusters,
                const ClusterHitRefs& cluster_hit_refs,
                const AHCALTLURawData* tlu) {
    TLatex text;
    text.SetNDC(true);
    text.SetTextSize(0.033);
    text.DrawLatex(0.06, 0.94, Form("Event entry: %lld", static_cast<long long>(entry)));
    text.DrawLatex(0.06, 0.89, Form("RecoHits: %zu", hits.size()));
    text.DrawLatex(0.06, 0.84, Form("%s: %zu", cluster_key.c_str(), clusters.size()));
    if (tlu) {
        text.DrawLatex(0.06, 0.79, Form("TLU EventID=%d  TriggerID=%d",
                                        tlu->EventID, tlu->TriggerID));
    } else if (by_tlu_eventid) {
        text.DrawLatex(0.06, 0.79, Form("Requested TLU EventID=%lld",
                                        static_cast<long long>(requested_event_id)));
    }

    text.SetTextSize(0.026);
    text.DrawLatex(0.06, 0.72, "Heatmap: all RecoHits, color = Edep [MeV]");
    text.DrawLatex(0.06, 0.68, "Colored boxes: member RecoHits in each cluster");
    text.DrawLatex(0.06, 0.64, "Star: cluster PCA mean");
    text.DrawLatex(0.06, 0.60, "Line: first PCA axis");
    text.DrawLatex(0.06, 0.56, Form("Matched boxed hits: %zu",
                                    count_cluster_hit_refs(cluster_hit_refs)));

    text.SetTextSize(0.024);
    text.DrawLatex(0.06, 0.49, "Top clusters by totalEdep:");
    text.DrawLatex(0.06, 0.45, "id  nHit   Edep   lin.   width   track");

    std::vector<std::size_t> order(clusters.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return clusters[a].totalEdep > clusters[b].totalEdep;
    });

    const std::size_t nshow = std::min<std::size_t>(order.size(), 12);
    double y = 0.405;
    for (std::size_t row = 0; row < nshow; ++row) {
        const std::size_t idx = order[row];
        const auto& c = clusters[idx];
        text.SetTextColor(cluster_color(static_cast<int>(idx)));
        text.DrawLatex(0.06, y,
                       Form("%2d  %4d  %7.2f  %5.2f  %6.1f  %s",
                            c.cluster_id,
                            c.nHitsInCluster,
                            c.totalEdep,
                            c.pca_linearity,
                            c.pca_width,
                            c.is_track_like ? "yes" : "no"));
        y -= 0.032;
    }

    text.SetTextColor(kBlack);
    text.SetTextSize(0.018);
    std::string shown = input_root;
    if (shown.size() > 72) shown = "..." + shown.substr(shown.size() - 69);
    text.DrawLatex(0.06, 0.03, shown.c_str());
}

} // namespace

int main(int argc, char** argv) {
    try {
        FAIR::init_logger("AHCALClusterEventDisplay", "", spdlog::level::warn);
        gROOT->SetBatch(kTRUE);

        const Options opt = parse_args(argc, argv);
        const std::string input_root = resolve_input_root(opt.input);

        ReaderRegistry rr;
        rr.register_vector_struct<AHCALRecoHit>("vector<AHCALRecoHit>");
        rr.register_vector_struct<AHCALRecoAlg::Cluster>("vector<Cluster>");
        rr.register_struct<AHCALTLURawData>("AHCALTLURawData");

        RootInput in(input_root, opt.tree_name);
        Long64_t entry = opt.event_id;
        AHCALTLURawData tlu{};
        bool have_tlu = false;

        if (opt.by_tlu_eventid) {
            bool found = false;
            for (Long64_t i = 0; i < in.entries(); ++i) {
                in.read_entry(i);
                const auto probe =
                    rr.read<AHCALTLURawData>("AHCALTLURawData", opt.tlu_key, in);
                if (probe.EventID == opt.event_id) {
                    entry = i;
                    tlu = probe;
                    have_tlu = true;
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cerr << "TLURawData.EventID not found: " << opt.event_id
                          << " (entries=" << in.entries() << ")\n";
                return 1;
            }
        } else if (!in.read_entry(entry)) {
            std::cerr << "Event entry out of range: " << entry
                      << " (entries=" << in.entries() << ")\n";
            return 1;
        }

        const auto hits =
            rr.read<std::vector<AHCALRecoHit>>("vector<AHCALRecoHit>", opt.hit_key, in);
        const auto clusters =
            rr.read<std::vector<AHCALRecoAlg::Cluster>>("vector<Cluster>", opt.cluster_key, in);
        const auto cluster_hit_refs = make_cluster_hit_refs(hits, clusters);

        AHCALTLURawData* tlu_ptr = nullptr;
        if (have_tlu) {
            tlu_ptr = &tlu;
        } else {
            try {
                tlu = rr.read<AHCALTLURawData>("AHCALTLURawData", opt.tlu_key, in);
                tlu_ptr = &tlu;
            } catch (const std::exception&) {
                tlu_ptr = nullptr;
            }
        }

        std::cout << "Event entry " << entry;
        if (opt.by_tlu_eventid) {
            std::cout << " (TLU EventID " << opt.event_id << ")";
        }
        std::cout << ": " << hits.size()
                  << " RecoHits, " << clusters.size() << " clusters from "
                  << opt.cluster_key << "\n";

        gStyle->SetOptStat(0);
        gStyle->SetNumberContours(100);
        gStyle->SetPalette(kViridis);

        TH2D hxy("cluster_hxy",
                 Form("RecoHits + %s XY;X [mm];Y [mm]", opt.cluster_key.c_str()),
                 18, xMin, xMax, 18, yMin, yMax);
        TH2D hxz("cluster_hxz",
                 Form("RecoHits + %s XZ;Z [mm];X [mm]", opt.cluster_key.c_str()),
                 135, zMin, zMax, 18, xMin, xMax);
        TH2D hyz("cluster_hyz",
                 Form("RecoHits + %s YZ;Z [mm];Y [mm]", opt.cluster_key.c_str()),
                 135, zMin, zMax, 18, yMin, yMax);

        for (const auto& hit : hits) {
            hxy.Fill(hit.Xpos(), hit.Ypos(), hit.Edep);
            hxz.Fill(hit.Zpos(), hit.Xpos(), hit.Edep);
            hyz.Fill(hit.Zpos(), hit.Ypos(), hit.Edep);
        }
        hxy.SetMinimum(0.0);
        hxz.SetMinimum(0.0);
        hyz.SetMinimum(0.0);

        TCanvas canvas("c_cluster", "AHCAL Cluster Event Display", 1800, 1100);
        canvas.Divide(2, 2);

        canvas.cd(1);
        gPad->SetRightMargin(0.14);
        hxy.Draw("COLZ");
        draw_cluster_hit_boxes_xy(cluster_hit_refs);
        draw_cluster_overlay_xy(clusters);

        canvas.cd(2);
        gPad->SetRightMargin(0.14);
        hxz.Draw("COLZ");
        draw_cluster_hit_boxes_xz(cluster_hit_refs);
        draw_cluster_overlay_xz(clusters);

        canvas.cd(3);
        gPad->SetRightMargin(0.14);
        hyz.Draw("COLZ");
        draw_cluster_hit_boxes_yz(cluster_hit_refs);
        draw_cluster_overlay_yz(clusters);

        canvas.cd(4);
        draw_panel(entry, opt.event_id, opt.by_tlu_eventid,
                   input_root, opt.cluster_key, hits, clusters,
                   cluster_hit_refs, tlu_ptr);

        const std::filesystem::path out_path(opt.output);
        if (!out_path.parent_path().empty()) {
            std::filesystem::create_directories(out_path.parent_path());
        }
        canvas.SaveAs(opt.output.c_str());
        std::cout << "Saved: " << opt.output << "\n";
    } catch (const std::exception& e) {
        std::cerr << "ClusterEventDisplay2D error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
