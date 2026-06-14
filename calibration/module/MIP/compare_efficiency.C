/**
 * Compare channel efficiencies from data and MC pre-fit MIP ROOT files.
 *
 * Usage:
 *   root -b -q 'compare_efficiency.C("data.root", "mc.root",
 *                                    "efficiency_comparison.root",
 *                                    "efficiency_plots",
 *                                    "data_fitted.root")'
 */

#include <TCanvas.h>
#include <TDirectory.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TKey.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace {

constexpr int kNLayers = 40;
constexpr int kNChips = 9;
constexpr int kNChannels = 36;
constexpr int kNChannelsPerLayer = kNChips * kNChannels;

struct Counts {
    Long64_t passed = 0;
    Long64_t total = 0;
    bool has_passed = false;
    bool has_total = false;
};

struct Sample {
    std::string name;
    std::string path;
    int color;
    int marker;
    std::map<int, Counts> counts;
    int invalid_bins = 0;
};

struct LayerObjects {
    std::unique_ptr<TH1D> passed;
    std::unique_ptr<TH1D> total;
    std::unique_ptr<TEfficiency> efficiency;
};

struct FitClassCategory {
    int key;
    const char* label;
    int color;
};

const std::vector<FitClassCategory> kFitClassCategories = {
    {0, "FitClass 0", kBlack},
    {4, "FitClass 4", kRed + 1},
    {5, "FitClass 5", kBlue + 1},
    {99, "Other FitClass", kGreen + 2},
};

int fitClassCategory(int fit_class) {
    return fit_class == 0 || fit_class == 4 || fit_class == 5 ? fit_class : 99;
}

bool loadFitClasses(const char* file_name, std::map<int, int>& fit_classes) {
    std::unique_ptr<TFile> file(TFile::Open(file_name, "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open Data fitted file: " << file_name
                  << std::endl;
        return false;
    }
    TTree* tree = nullptr;
    file->GetObject("mip_fit_results", tree);
    if (!tree || !tree->GetBranch("cellid") || !tree->GetBranch("FitClass")) {
        std::cerr << "Error: cellid/FitClass is missing in " << file_name
                  << std::endl;
        return false;
    }
    int cellid = -1;
    int fit_class = -999;
    tree->SetBranchAddress("cellid", &cellid);
    tree->SetBranchAddress("FitClass", &fit_class);
    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
        tree->GetEntry(entry);
        fit_classes[cellid] = fit_class;
    }
    return true;
}

void collectCounts(TDirectory* dir, std::map<int, Counts>& counts) {
    if (!dir) return;

    TIter next(dir->GetListOfKeys());
    TKey* key = nullptr;
    while ((key = static_cast<TKey*>(next()))) {
        TObject* obj = dir->Get(key->GetName());
        if (!obj) continue;

        if (obj->InheritsFrom(TDirectory::Class())) {
            collectCounts(static_cast<TDirectory*>(obj), counts);
        } else if (obj->InheritsFrom(TH1::Class())) {
            const std::string name = obj->GetName();
            const std::string prefix = "hMIP_";
            if (name.find(prefix) == 0) {
                const int cellid = std::stoi(name.substr(prefix.size()));
                counts[cellid].passed =
                    static_cast<Long64_t>(std::llround(
                        static_cast<TH1*>(obj)->GetEntries()));
                counts[cellid].has_passed = true;
            }
        } else if (obj->InheritsFrom(TParameter<int>::Class())) {
            const std::string name = obj->GetName();
            const std::string prefix = "Ntracks_pass_through_channel_";
            if (name.find(prefix) == 0) {
                const int cellid = std::stoi(name.substr(prefix.size()));
                counts[cellid].total =
                    static_cast<TParameter<int>*>(obj)->GetVal();
                counts[cellid].has_total = true;
            }
        }
    }
}

bool loadSample(Sample& sample) {
    std::unique_ptr<TFile> file(TFile::Open(sample.path.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open " << sample.name << " file: "
                  << sample.path << std::endl;
        return false;
    }

    TDirectory* mip_dir = file->GetDirectory("MIP");
    collectCounts(mip_dir ? mip_dir : file.get(), sample.counts);
    std::cout << sample.name << ": collected " << sample.counts.size()
              << " cell IDs" << std::endl;
    return true;
}

bool validCounts(const Counts& counts) {
    return counts.has_passed && counts.has_total && counts.total > 0 &&
           counts.passed >= 0 && counts.passed <= counts.total;
}

std::vector<LayerObjects> makeLayerObjects(Sample& sample) {
    std::vector<LayerObjects> layers;
    layers.reserve(kNLayers);

    for (int layer = 0; layer < kNLayers; ++layer) {
        LayerObjects objects;
        objects.passed = std::make_unique<TH1D>(
            Form("passed_%s_layer%02d", sample.name.c_str(), layer),
            Form("%s passed, layer %d;chip #times 36 + channel;Entries",
                 sample.name.c_str(), layer),
            kNChannelsPerLayer, -0.5, kNChannelsPerLayer - 0.5);
        objects.total = std::make_unique<TH1D>(
            Form("total_%s_layer%02d", sample.name.c_str(), layer),
            Form("%s total, layer %d;chip #times 36 + channel;Tracks",
                 sample.name.c_str(), layer),
            kNChannelsPerLayer, -0.5, kNChannelsPerLayer - 0.5);
        objects.passed->SetDirectory(nullptr);
        objects.total->SetDirectory(nullptr);

        for (int chip = 0; chip < kNChips; ++chip) {
            for (int channel = 0; channel < kNChannels; ++channel) {
                const int cellid =
                    layer * 100000 + chip * 10000 + channel;
                const int bin = chip * kNChannels + channel + 1;
                const auto it = sample.counts.find(cellid);
                if (it == sample.counts.end()) continue;
                if (!validCounts(it->second)) {
                    ++sample.invalid_bins;
                    continue;
                }
                objects.passed->SetBinContent(bin, it->second.passed);
                objects.total->SetBinContent(bin, it->second.total);
            }
        }

        objects.efficiency =
            std::make_unique<TEfficiency>(*objects.passed, *objects.total);
        objects.efficiency->SetName(
            Form("efficiency_%s_layer%02d", sample.name.c_str(), layer));
        objects.efficiency->SetTitle(
            Form("Layer %d efficiency;chip #times 36 + channel;Efficiency",
                 layer));
        objects.efficiency->SetLineColor(sample.color);
        objects.efficiency->SetMarkerColor(sample.color);
        objects.efficiency->SetMarkerStyle(sample.marker);
        objects.efficiency->SetMarkerSize(0.7);
        layers.push_back(std::move(objects));
    }
    return layers;
}

std::vector<LayerObjects> makeChipObjects(
    const Sample& sample, const std::vector<LayerObjects>& channel_layers) {
    std::vector<LayerObjects> layers;
    layers.reserve(kNLayers);

    for (int layer = 0; layer < kNLayers; ++layer) {
        LayerObjects objects;
        objects.passed = std::make_unique<TH1D>(
            Form("passed_%s_chip_layer%02d", sample.name.c_str(), layer),
            Form("%s passed by chip, layer %d;Chip;Entries",
                 sample.name.c_str(), layer),
            kNChips, -0.5, kNChips - 0.5);
        objects.total = std::make_unique<TH1D>(
            Form("total_%s_chip_layer%02d", sample.name.c_str(), layer),
            Form("%s total by chip, layer %d;Chip;Tracks",
                 sample.name.c_str(), layer),
            kNChips, -0.5, kNChips - 0.5);
        objects.passed->SetDirectory(nullptr);
        objects.total->SetDirectory(nullptr);

        for (int chip = 0; chip < kNChips; ++chip) {
            const int first_bin = chip * kNChannels + 1;
            const int last_bin = first_bin + kNChannels - 1;
            objects.passed->SetBinContent(
                chip + 1,
                channel_layers[layer].passed->Integral(first_bin, last_bin));
            objects.total->SetBinContent(
                chip + 1,
                channel_layers[layer].total->Integral(first_bin, last_bin));
        }

        objects.efficiency =
            std::make_unique<TEfficiency>(*objects.passed, *objects.total);
        objects.efficiency->SetName(
            Form("efficiency_%s_chip_layer%02d", sample.name.c_str(), layer));
        objects.efficiency->SetTitle(
            Form("Layer %d chip efficiency;Chip;Efficiency", layer));
        objects.efficiency->SetLineColor(sample.color);
        objects.efficiency->SetMarkerColor(sample.color);
        objects.efficiency->SetMarkerStyle(sample.marker);
        layers.push_back(std::move(objects));
    }
    return layers;
}

std::unique_ptr<TEfficiency> makeLayerSummary(
    const Sample& sample, const std::vector<LayerObjects>& layers) {
    TH1D passed(Form("passed_%s_by_layer", sample.name.c_str()),
                Form("%s passed by layer;Layer;Entries", sample.name.c_str()),
                kNLayers, -0.5, kNLayers - 0.5);
    TH1D total(Form("total_%s_by_layer", sample.name.c_str()),
               Form("%s total by layer;Layer;Tracks", sample.name.c_str()),
               kNLayers, -0.5, kNLayers - 0.5);

    for (int layer = 0; layer < kNLayers; ++layer) {
        passed.SetBinContent(layer + 1, layers[layer].passed->Integral());
        total.SetBinContent(layer + 1, layers[layer].total->Integral());
    }

    auto efficiency = std::make_unique<TEfficiency>(passed, total);
    efficiency->SetName(Form("efficiency_%s_by_layer", sample.name.c_str()));
    efficiency->SetTitle("Efficiency comparison by layer;Layer;Efficiency");
    efficiency->SetLineColor(sample.color);
    efficiency->SetMarkerColor(sample.color);
    efficiency->SetMarkerStyle(sample.marker);
    return efficiency;
}

LayerObjects makeCombinedObjects(const Sample& sample,
                                 const std::vector<LayerObjects>& layers,
                                 int bins_per_layer, const char* type,
                                 const char* x_title) {
    const int n_bins = kNLayers * bins_per_layer;
    LayerObjects objects;
    objects.passed = std::make_unique<TH1D>(
        Form("passed_%s_%s", sample.name.c_str(), type),
        Form("%s passed by %s;%s;Entries", sample.name.c_str(), type, x_title),
        n_bins, -0.5, n_bins - 0.5);
    objects.total = std::make_unique<TH1D>(
        Form("total_%s_%s", sample.name.c_str(), type),
        Form("%s total by %s;%s;Tracks", sample.name.c_str(), type, x_title),
        n_bins, -0.5, n_bins - 0.5);
    objects.passed->SetDirectory(nullptr);
    objects.total->SetDirectory(nullptr);

    for (int layer = 0; layer < kNLayers; ++layer) {
        for (int local_bin = 1; local_bin <= bins_per_layer; ++local_bin) {
            const int global_bin = layer * bins_per_layer + local_bin;
            objects.passed->SetBinContent(
                global_bin, layers[layer].passed->GetBinContent(local_bin));
            objects.total->SetBinContent(
                global_bin, layers[layer].total->GetBinContent(local_bin));
        }
    }

    objects.efficiency =
        std::make_unique<TEfficiency>(*objects.passed, *objects.total);
    objects.efficiency->SetName(
        Form("efficiency_%s_%s", sample.name.c_str(), type));
    objects.efficiency->SetTitle(
        Form("%s efficiency by %s;%s;Efficiency", sample.name.c_str(), type,
             x_title));
    objects.efficiency->SetLineColor(sample.color);
    objects.efficiency->SetMarkerColor(sample.color);
    objects.efficiency->SetMarkerStyle(sample.marker);
    objects.efficiency->SetMarkerSize(0.7);
    return objects;
}

double symmetricError(const TEfficiency& efficiency, int bin) {
    return std::max(efficiency.GetEfficiencyErrorLow(bin),
                    efficiency.GetEfficiencyErrorUp(bin));
}

std::unique_ptr<TGraphErrors> makeRatioGraph(
    const TEfficiency& data, const TEfficiency& mc, int n_bins,
    const char* name, const char* title, int color, int marker) {
    auto graph = std::make_unique<TGraphErrors>();
    graph->SetName(name);
    graph->SetTitle(title);
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(0.7);

    int point = 0;
    for (int bin = 1; bin <= n_bins; ++bin) {
        if (data.GetTotalHistogram()->GetBinContent(bin) <= 0 ||
            mc.GetTotalHistogram()->GetBinContent(bin) <= 0) {
            continue;
        }
        const double data_efficiency = data.GetEfficiency(bin);
        const double mc_efficiency = mc.GetEfficiency(bin);
        if (data_efficiency <= 0 || mc_efficiency <= 0) continue;

        const double ratio = data_efficiency / mc_efficiency;
        const double relative_data_error =
            symmetricError(data, bin) / data_efficiency;
        const double relative_mc_error =
            symmetricError(mc, bin) / mc_efficiency;
        const double error =
            ratio * std::hypot(relative_data_error, relative_mc_error);
        graph->SetPoint(point, bin - 1, ratio);
        graph->SetPointError(point, 0.0, error);
        ++point;
    }
    return graph;
}

double averageEfficiency(const TEfficiency& efficiency) {
    double sum = 0.0;
    int valid_bins = 0;
    const int n_bins = efficiency.GetTotalHistogram()->GetNbinsX();
    for (int bin = 1; bin <= n_bins; ++bin) {
        if (efficiency.GetTotalHistogram()->GetBinContent(bin) <= 0) continue;
        sum += efficiency.GetEfficiency(bin);
        ++valid_bins;
    }
    return valid_bins > 0 ? sum / valid_bins : -1.0;
}

void drawEfficiencyAverages(const TEfficiency& data, const TEfficiency& mc) {
    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.035);

    const double data_average = averageEfficiency(data);
    label.SetTextColor(data.GetLineColor());
    label.DrawLatex(0.14, 0.34,
                    data_average >= 0
                        ? Form("Data mean = %.4f", data_average)
                        : "Data mean = N/A");

    const double mc_average = averageEfficiency(mc);
    label.SetTextColor(mc.GetLineColor());
    label.DrawLatex(0.14, 0.28,
                    mc_average >= 0 ? Form("MC mean = %.4f", mc_average)
                                    : "MC mean = N/A");
}

using EfficiencyCategoryGraphs =
    std::map<int, std::unique_ptr<TGraphAsymmErrors>>;
using RatioCategoryGraphs = std::map<int, std::unique_ptr<TGraphErrors>>;

EfficiencyCategoryGraphs makeEfficiencyCategoryGraphs(
    const TEfficiency& efficiency, int layer,
    const std::map<int, int>& fit_classes) {
    EfficiencyCategoryGraphs graphs;
    for (const auto& category : kFitClassCategories) {
        auto graph = std::make_unique<TGraphAsymmErrors>();
        graph->SetName(
            Form("data_efficiency_layer%02d_fitClassCategory%d", layer,
                 category.key));
        graph->SetLineColor(category.color);
        graph->SetMarkerColor(category.color);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(0.7);
        graphs[category.key] = std::move(graph);
    }

    for (int index = 0; index < kNChannelsPerLayer; ++index) {
        const int bin = index + 1;
        if (efficiency.GetTotalHistogram()->GetBinContent(bin) <= 0) continue;
        const int chip = index / kNChannels;
        const int channel = index % kNChannels;
        const int cellid = layer * 100000 + chip * 10000 + channel;
        const auto fit_it = fit_classes.find(cellid);
        const int category =
            fitClassCategory(fit_it == fit_classes.end() ? -999 : fit_it->second);
        auto& graph = graphs[category];
        const int point = graph->GetN();
        graph->SetPoint(point, index, efficiency.GetEfficiency(bin));
        graph->SetPointError(point, 0.0, 0.0,
                             efficiency.GetEfficiencyErrorLow(bin),
                             efficiency.GetEfficiencyErrorUp(bin));
    }
    return graphs;
}

RatioCategoryGraphs makeRatioCategoryGraphs(
    const TEfficiency& data, const TEfficiency& mc, int layer,
    const std::map<int, int>& fit_classes) {
    RatioCategoryGraphs graphs;
    for (const auto& category : kFitClassCategories) {
        auto graph = std::make_unique<TGraphErrors>();
        graph->SetName(Form("data_over_mc_layer%02d_fitClassCategory%d", layer,
                            category.key));
        graph->SetLineColor(category.color);
        graph->SetMarkerColor(category.color);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(0.7);
        graphs[category.key] = std::move(graph);
    }

    for (int index = 0; index < kNChannelsPerLayer; ++index) {
        const int bin = index + 1;
        if (data.GetTotalHistogram()->GetBinContent(bin) <= 0 ||
            mc.GetTotalHistogram()->GetBinContent(bin) <= 0) {
            continue;
        }
        const double data_efficiency = data.GetEfficiency(bin);
        const double mc_efficiency = mc.GetEfficiency(bin);
        if (data_efficiency <= 0 || mc_efficiency <= 0) continue;
        const double ratio = data_efficiency / mc_efficiency;
        const double error =
            ratio * std::hypot(symmetricError(data, bin) / data_efficiency,
                               symmetricError(mc, bin) / mc_efficiency);
        const int chip = index / kNChannels;
        const int channel = index % kNChannels;
        const int cellid = layer * 100000 + chip * 10000 + channel;
        const auto fit_it = fit_classes.find(cellid);
        const int category =
            fitClassCategory(fit_it == fit_classes.end() ? -999 : fit_it->second);
        auto& graph = graphs[category];
        const int point = graph->GetN();
        graph->SetPoint(point, index, ratio);
        graph->SetPointError(point, 0.0, error);
    }
    return graphs;
}

EfficiencyCategoryGraphs makeCombinedEfficiencyCategoryGraphs(
    const TEfficiency& efficiency, const std::map<int, int>& fit_classes) {
    EfficiencyCategoryGraphs graphs;
    for (const auto& category : kFitClassCategories) {
        auto graph = std::make_unique<TGraphAsymmErrors>();
        graph->SetName(
            Form("data_efficiency_channel_fitClassCategory%d", category.key));
        graph->SetLineColor(category.color);
        graph->SetMarkerColor(category.color);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(0.7);
        graphs[category.key] = std::move(graph);
    }

    for (int index = 0; index < kNLayers * kNChannelsPerLayer; ++index) {
        const int bin = index + 1;
        if (efficiency.GetTotalHistogram()->GetBinContent(bin) <= 0) continue;
        const int layer = index / kNChannelsPerLayer;
        const int local_index = index % kNChannelsPerLayer;
        const int chip = local_index / kNChannels;
        const int channel = local_index % kNChannels;
        const int cellid = layer * 100000 + chip * 10000 + channel;
        const auto fit_it = fit_classes.find(cellid);
        const int category =
            fitClassCategory(fit_it == fit_classes.end() ? -999 : fit_it->second);
        auto& graph = graphs[category];
        const int point = graph->GetN();
        graph->SetPoint(point, index, efficiency.GetEfficiency(bin));
        graph->SetPointError(point, 0.0, 0.0,
                             efficiency.GetEfficiencyErrorLow(bin),
                             efficiency.GetEfficiencyErrorUp(bin));
    }
    return graphs;
}

RatioCategoryGraphs makeCombinedRatioCategoryGraphs(
    const TEfficiency& data, const TEfficiency& mc,
    const std::map<int, int>& fit_classes) {
    RatioCategoryGraphs graphs;
    for (const auto& category : kFitClassCategories) {
        auto graph = std::make_unique<TGraphErrors>();
        graph->SetName(
            Form("data_over_mc_channel_fitClassCategory%d", category.key));
        graph->SetLineColor(category.color);
        graph->SetMarkerColor(category.color);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(0.7);
        graphs[category.key] = std::move(graph);
    }

    for (int index = 0; index < kNLayers * kNChannelsPerLayer; ++index) {
        const int bin = index + 1;
        if (data.GetTotalHistogram()->GetBinContent(bin) <= 0 ||
            mc.GetTotalHistogram()->GetBinContent(bin) <= 0) {
            continue;
        }
        const double data_efficiency = data.GetEfficiency(bin);
        const double mc_efficiency = mc.GetEfficiency(bin);
        if (data_efficiency <= 0 || mc_efficiency <= 0) continue;
        const double ratio = data_efficiency / mc_efficiency;
        const double error =
            ratio * std::hypot(symmetricError(data, bin) / data_efficiency,
                               symmetricError(mc, bin) / mc_efficiency);
        const int layer = index / kNChannelsPerLayer;
        const int local_index = index % kNChannelsPerLayer;
        const int chip = local_index / kNChannels;
        const int channel = local_index % kNChannels;
        const int cellid = layer * 100000 + chip * 10000 + channel;
        const auto fit_it = fit_classes.find(cellid);
        const int category =
            fitClassCategory(fit_it == fit_classes.end() ? -999 : fit_it->second);
        auto& graph = graphs[category];
        const int point = graph->GetN();
        graph->SetPoint(point, index, ratio);
        graph->SetPointError(point, 0.0, error);
    }
    return graphs;
}

void styleEfficiencyGraph(TEfficiency& efficiency, double x_max,
                          double y_min = 0.0, double y_max = 1.05) {
    efficiency.Draw("AP");
    gPad->Update();
    TGraphAsymmErrors* graph = efficiency.GetPaintedGraph();
    if (!graph) return;
    graph->GetXaxis()->SetLimits(-0.5, x_max + 0.5);
    graph->GetYaxis()->SetRangeUser(y_min, y_max);
    graph->GetXaxis()->SetTitleOffset(1.1);
    graph->GetYaxis()->SetTitleOffset(1.1);
}

void drawLayerBoundaries(int bins_per_layer, double y_min, double y_max,
                         bool draw_labels) {
    for (int layer = 0; layer < kNLayers; ++layer) {
        TLine* line =
            new TLine(layer * bins_per_layer - 0.5, y_min,
                      layer * bins_per_layer - 0.5, y_max);
        line->SetLineColor(kRed);
        line->SetLineStyle(2);
        line->Draw("SAME");
        if (draw_labels) {
            TLatex* label = new TLatex(
                layer * bins_per_layer + bins_per_layer * 0.15,
                y_max - 0.04 * (y_max - y_min), Form("Layer %d", layer));
            label->SetTextSize(0.025);
            label->SetTextFont(42);
            label->SetTextColor(kBlack);
            label->Draw();
        }
    }
}

void drawCombinedChannelComparison(
    TEfficiency& mc, TGraphErrors& ratio,
    EfficiencyCategoryGraphs& data_categories,
    RatioCategoryGraphs& ratio_categories, TDirectory* canvas_dir,
    const std::string& output_dir) {
    constexpr int n_bins = kNLayers * kNChannelsPerLayer;
    TCanvas canvas("c_eff_channel", "Efficiency comparison by channel", 8000,
                   1200);
    canvas.Divide(1, 2);

    canvas.cd(1);
    gPad->SetGrid();
    gPad->SetBottomMargin(0.03);
    mc.SetTitle(
        "Efficiency comparison by channel;layer*9*36 + chip*36 + channel;"
        "Efficiency");
    styleEfficiencyGraph(mc, n_bins - 1, 0.0, 1.10);
    for (const auto& category : kFitClassCategories) {
        data_categories[category.key]->Draw("P SAME");
    }
    drawLayerBoundaries(kNChannelsPerLayer, 0.0, 1.10, true);
    TLegend legend(0.96, 0.15, 0.995, 0.42);
    legend.AddEntry(&mc, "MC", "lp");
    for (const auto& category : kFitClassCategories) {
        legend.AddEntry(data_categories[category.key].get(), category.label,
                        "lp");
    }
    legend.Draw();

    canvas.cd(2);
    gPad->SetGrid();
    gPad->SetTopMargin(0.03);
    ratio.SetLineColor(0);
    ratio.SetMarkerColor(0);
    ratio.Draw("AP");
    ratio.GetXaxis()->SetLimits(-0.5, n_bins - 0.5);
    ratio.GetXaxis()->SetTitle("layer*9*36 + chip*36 + channel");
    ratio.GetYaxis()->SetTitle("Data / MC efficiency");
    for (const auto& category : kFitClassCategories) {
        ratio_categories[category.key]->Draw("P SAME");
    }
    drawLayerBoundaries(kNChannelsPerLayer, gPad->GetUymin(), gPad->GetUymax(),
                        false);
    TLine unity(-0.5, 1.0, n_bins - 0.5, 1.0);
    unity.SetLineStyle(2);
    unity.Draw();

    canvas.Modified();
    canvas.Update();
    if (canvas_dir) {
        canvas_dir->cd();
        canvas.Write();
    }
    canvas.SaveAs(Form("%s/eff_channel.pdf", output_dir.c_str()));
    canvas.SaveAs(Form("%s/eff_channel.png", output_dir.c_str()));
}

void drawCombinedChipComparison(TEfficiency& data, TEfficiency& mc,
                                TGraphErrors& ratio,
                                TDirectory* canvas_dir,
                                const std::string& output_dir) {
    constexpr int n_bins = kNLayers * kNChips;
    TCanvas canvas("c_eff_chip", "Efficiency comparison by chip", 4000, 1200);
    canvas.Divide(1, 2);

    canvas.cd(1);
    gPad->SetGrid();
    gPad->SetBottomMargin(0.03);
    data.SetTitle("Efficiency comparison by chip;9*Layer + chip;Efficiency");
    styleEfficiencyGraph(data, n_bins - 1, 0.0, 1.10);
    mc.Draw("P SAME");
    drawLayerBoundaries(kNChips, 0.0, 1.10, true);
    TLegend legend(0.93, 0.15, 0.99, 0.30);
    legend.AddEntry(&data, "Data", "lp");
    legend.AddEntry(&mc, "MC", "lp");
    legend.Draw();

    canvas.cd(2);
    gPad->SetGrid();
    gPad->SetTopMargin(0.03);
    ratio.Draw("AP");
    ratio.GetXaxis()->SetLimits(-0.5, n_bins - 0.5);
    ratio.GetXaxis()->SetTitle("9*Layer + chip");
    ratio.GetYaxis()->SetTitle("Data / MC efficiency");
    drawLayerBoundaries(kNChips, gPad->GetUymin(), gPad->GetUymax(), false);
    TLine unity(-0.5, 1.0, n_bins - 0.5, 1.0);
    unity.SetLineStyle(2);
    unity.Draw();

    canvas.Modified();
    canvas.Update();
    if (canvas_dir) {
        canvas_dir->cd();
        canvas.Write();
    }
    canvas.SaveAs(Form("%s/eff_chip.pdf", output_dir.c_str()));
    canvas.SaveAs(Form("%s/eff_chip.png", output_dir.c_str()));
}

void drawLayerComparison(int layer, TEfficiency& data, TEfficiency& mc,
                         TGraphErrors& ratio,
                         EfficiencyCategoryGraphs& data_categories,
                         RatioCategoryGraphs& ratio_categories,
                         TDirectory* canvas_dir,
                         const std::string& output_dir) {
    TCanvas canvas(Form("comparison_layer%02d", layer),
                   Form("Efficiency comparison, layer %d", layer),
                   1600, 1000);
    canvas.Divide(1, 2);

    canvas.cd(1);
    gPad->SetGrid();
    gPad->SetBottomMargin(0.03);
    styleEfficiencyGraph(mc, kNChannelsPerLayer - 1);
    for (const auto& category : kFitClassCategories) {
        data_categories[category.key]->Draw("P SAME");
    }
    drawEfficiencyAverages(data, mc);
    TLegend legend(0.72, 0.15, 0.89, 0.35);
    legend.AddEntry(&mc, "MC", "lp");
    for (const auto& category : kFitClassCategories) {
        legend.AddEntry(data_categories[category.key].get(), category.label,
                        "lp");
    }
    legend.Draw();

    canvas.cd(2);
    gPad->SetGrid();
    gPad->SetTopMargin(0.03);
    ratio.SetLineColor(0);
    ratio.SetMarkerColor(0);
    ratio.Draw("AP");
    ratio.GetXaxis()->SetLimits(-0.5, kNChannelsPerLayer - 0.5);
    ratio.GetYaxis()->SetTitle("Data / MC efficiency");
    ratio.GetXaxis()->SetTitle("chip #times 36 + channel");
    for (const auto& category : kFitClassCategories) {
        ratio_categories[category.key]->Draw("P SAME");
    }
    TLine unity(-0.5, 1.0, kNChannelsPerLayer - 0.5, 1.0);
    unity.SetLineStyle(2);
    unity.Draw();
    TLegend ratio_legend(0.72, 0.15, 0.89, 0.32);
    for (const auto& category : kFitClassCategories) {
        ratio_legend.AddEntry(ratio_categories[category.key].get(),
                              category.label, "lp");
    }
    ratio_legend.Draw();

    canvas.Modified();
    canvas.Update();
    if (canvas_dir) {
        canvas_dir->cd();
        canvas.Write();
    }
    canvas.SaveAs(
        Form("%s/efficiency_layer%02d.pdf", output_dir.c_str(), layer));
    canvas.SaveAs(
        Form("%s/efficiency_layer%02d.png", output_dir.c_str(), layer));
}

void drawChipComparison(int layer, TEfficiency& data, TEfficiency& mc,
                        TGraphErrors& ratio, TDirectory* canvas_dir,
                        const std::string& output_dir) {
    TCanvas canvas(Form("chip_comparison_layer%02d", layer),
                   Form("Chip efficiency comparison, layer %d", layer),
                   1200, 900);
    canvas.Divide(1, 2);

    canvas.cd(1);
    gPad->SetGrid();
    styleEfficiencyGraph(data, kNChips - 1);
    mc.Draw("P SAME");
    drawEfficiencyAverages(data, mc);
    TLegend legend(0.72, 0.15, 0.89, 0.35);
    legend.AddEntry(&data, "Data", "lp");
    legend.AddEntry(&mc, "MC", "lp");
    legend.Draw();

    canvas.cd(2);
    gPad->SetGrid();
    ratio.Draw("AP");
    ratio.GetXaxis()->SetLimits(-0.5, kNChips - 0.5);
    ratio.GetXaxis()->SetTitle("Chip");
    ratio.GetYaxis()->SetTitle("Data / MC efficiency");
    TLine unity(-0.5, 1.0, kNChips - 0.5, 1.0);
    unity.SetLineStyle(2);
    unity.Draw();
    TLegend ratio_legend(0.72, 0.15, 0.89, 0.32);
    ratio_legend.AddEntry(&ratio, "Data / MC", "lp");
    ratio_legend.Draw();

    canvas.Modified();
    canvas.Update();
    if (canvas_dir) {
        canvas_dir->cd();
        canvas.Write();
    }
    canvas.SaveAs(
        Form("%s/efficiency_chip_layer%02d.pdf", output_dir.c_str(), layer));
    canvas.SaveAs(
        Form("%s/efficiency_chip_layer%02d.png", output_dir.c_str(), layer));
}

}  // namespace

void compare_efficiency(const char* data_file = "/eos/user/y/ymaruya/FASER/AHCAL/FAIR/out_calibration/mip_calib/21987-22133/mip_neighborcheck_nofit.root", 
                        const char* mc_file = "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out/SimMuon/mip_neighborcheck_nofit.root",
                        const char* output_root = "efficiency_comparison.root",
                        const char* output_dir = "efficiency_plots",
                        const char* data_fitted_file = "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out_calibration/mip_calib/21987-22133/mip_fitted_noplot_7.root") {
    gStyle->SetOptStat(0);
    gSystem->mkdir(output_dir, true);

    Sample data{"data", data_file, kBlack, 20};
    Sample mc{"mc", mc_file, kRed + 1, 24};
    if (!loadSample(data) || !loadSample(mc)) return;
    std::map<int, int> fit_classes;
    if (!loadFitClasses(data_fitted_file, fit_classes)) return;

    auto data_layers = makeLayerObjects(data);
    auto mc_layers = makeLayerObjects(mc);
    auto data_chips = makeChipObjects(data, data_layers);
    auto mc_chips = makeChipObjects(mc, mc_layers);
    auto data_channels_all =
        makeCombinedObjects(data, data_layers, kNChannelsPerLayer, "channel",
                            "layer*9*36 + chip*36 + channel");
    auto mc_channels_all =
        makeCombinedObjects(mc, mc_layers, kNChannelsPerLayer, "channel",
                            "layer*9*36 + chip*36 + channel");
    auto data_chips_all =
        makeCombinedObjects(data, data_chips, kNChips, "chip",
                            "9*Layer + chip");
    auto mc_chips_all =
        makeCombinedObjects(mc, mc_chips, kNChips, "chip",
                            "9*Layer + chip");
    auto data_summary = makeLayerSummary(data, data_layers);
    auto mc_summary = makeLayerSummary(mc, mc_layers);

    std::unique_ptr<TFile> output(TFile::Open(output_root, "RECREATE"));
    if (!output || output->IsZombie()) {
        std::cerr << "Error: cannot create output ROOT file: " << output_root
                  << std::endl;
        return;
    }

    TDirectory* data_dir = output->mkdir("data");
    TDirectory* mc_dir = output->mkdir("mc");
    TDirectory* comparison_dir = output->mkdir("comparison");
    TDirectory* canvas_dir = output->mkdir("canvas");

    int cellid = -1;
    int layer = -1;
    int chip = -1;
    int channel = -1;
    double data_efficiency = -1.0;
    double mc_efficiency = -1.0;
    double data_over_mc = -1.0;
    TTree comparison_tree("efficiency_comparison",
                          "Per-channel efficiency comparison");
    comparison_tree.SetDirectory(comparison_dir);
    comparison_tree.Branch("cellid", &cellid);
    comparison_tree.Branch("layer", &layer);
    comparison_tree.Branch("chip", &chip);
    comparison_tree.Branch("channel", &channel);
    comparison_tree.Branch("data_efficiency", &data_efficiency);
    comparison_tree.Branch("mc_efficiency", &mc_efficiency);
    comparison_tree.Branch("data_over_mc", &data_over_mc);

    int chip_layer = -1;
    int chip_number = -1;
    double chip_data_efficiency = -1.0;
    double chip_mc_efficiency = -1.0;
    double chip_data_over_mc = -1.0;
    TTree chip_comparison_tree("chip_efficiency_comparison",
                               "Per-chip efficiency comparison");
    chip_comparison_tree.SetDirectory(comparison_dir);
    chip_comparison_tree.Branch("layer", &chip_layer);
    chip_comparison_tree.Branch("chip", &chip_number);
    chip_comparison_tree.Branch("data_efficiency", &chip_data_efficiency);
    chip_comparison_tree.Branch("mc_efficiency", &chip_mc_efficiency);
    chip_comparison_tree.Branch("data_over_mc", &chip_data_over_mc);

    for (layer = 0; layer < kNLayers; ++layer) {
        data_dir->cd();
        data_layers[layer].passed->Write();
        data_layers[layer].total->Write();
        data_layers[layer].efficiency->Write();
        data_chips[layer].passed->Write();
        data_chips[layer].total->Write();
        data_chips[layer].efficiency->Write();
        mc_dir->cd();
        mc_layers[layer].passed->Write();
        mc_layers[layer].total->Write();
        mc_layers[layer].efficiency->Write();
        mc_chips[layer].passed->Write();
        mc_chips[layer].total->Write();
        mc_chips[layer].efficiency->Write();

        auto ratio = makeRatioGraph(
            *data_layers[layer].efficiency, *mc_layers[layer].efficiency,
            kNChannelsPerLayer, Form("data_over_mc_layer%02d", layer),
            Form("Layer %d Data / MC;chip #times 36 + channel;Ratio",
                 layer),
            mc.color, mc.marker);
        comparison_dir->cd();
        ratio->Write();
        auto data_categories = makeEfficiencyCategoryGraphs(
            *data_layers[layer].efficiency, layer, fit_classes);
        auto ratio_categories = makeRatioCategoryGraphs(
            *data_layers[layer].efficiency, *mc_layers[layer].efficiency,
            layer, fit_classes);
        for (const auto& category : kFitClassCategories) {
            data_categories[category.key]->Write();
            ratio_categories[category.key]->Write();
        }

        drawLayerComparison(layer, *data_layers[layer].efficiency,
                            *mc_layers[layer].efficiency, *ratio,
                            data_categories, ratio_categories,
                            canvas_dir, output_dir);

        auto chip_ratio = makeRatioGraph(
            *data_chips[layer].efficiency, *mc_chips[layer].efficiency,
            kNChips, Form("data_over_mc_chip_layer%02d", layer),
            Form("Layer %d chip Data / MC;Chip;Ratio", layer),
            mc.color, mc.marker);
        comparison_dir->cd();
        chip_ratio->Write();
        drawChipComparison(layer, *data_chips[layer].efficiency,
                           *mc_chips[layer].efficiency, *chip_ratio,
                           canvas_dir, output_dir);

        chip_layer = layer;
        for (chip_number = 0; chip_number < kNChips; ++chip_number) {
            const int bin = chip_number + 1;
            chip_data_efficiency =
                data_chips[layer].efficiency->GetTotalHistogram()
                            ->GetBinContent(bin) > 0
                    ? data_chips[layer].efficiency->GetEfficiency(bin)
                    : -1.0;
            chip_mc_efficiency =
                mc_chips[layer].efficiency->GetTotalHistogram()
                            ->GetBinContent(bin) > 0
                    ? mc_chips[layer].efficiency->GetEfficiency(bin)
                    : -1.0;
            chip_data_over_mc =
                chip_data_efficiency > 0 && chip_mc_efficiency > 0
                    ? chip_data_efficiency / chip_mc_efficiency
                    : -1.0;
            chip_comparison_tree.Fill();
        }
        
        for (chip = 0; chip < kNChips; ++chip) {
            for (channel = 0; channel < kNChannels; ++channel) {
                const int bin = chip * kNChannels + channel + 1;
                cellid = layer * 100000 + chip * 10000 + channel;
                data_efficiency =
                    data_layers[layer].efficiency->GetTotalHistogram()
                                ->GetBinContent(bin) > 0
                        ? data_layers[layer].efficiency->GetEfficiency(bin)
                        : -1.0;
                mc_efficiency =
                    mc_layers[layer].efficiency->GetTotalHistogram()
                                ->GetBinContent(bin) > 0
                        ? mc_layers[layer].efficiency->GetEfficiency(bin)
                        : -1.0;
                data_over_mc =
                    data_efficiency > 0 && mc_efficiency > 0
                        ? data_efficiency / mc_efficiency
                        : -1.0;
                comparison_tree.Fill();
            }
        }
    }

    data_dir->cd();
    data_channels_all.passed->Write();
    data_channels_all.total->Write();
    data_channels_all.efficiency->Write();
    data_chips_all.passed->Write();
    data_chips_all.total->Write();
    data_chips_all.efficiency->Write();
    data_summary->Write();
    mc_dir->cd();
    mc_channels_all.passed->Write();
    mc_channels_all.total->Write();
    mc_channels_all.efficiency->Write();
    mc_chips_all.passed->Write();
    mc_chips_all.total->Write();
    mc_chips_all.efficiency->Write();
    mc_summary->Write();

    auto channel_ratio_all = makeRatioGraph(
        *data_channels_all.efficiency, *mc_channels_all.efficiency,
        kNLayers * kNChannelsPerLayer, "data_over_mc_channel",
        "Data / MC by channel;layer*9*36 + chip*36 + channel;Ratio", mc.color,
        mc.marker);
    auto channel_data_categories = makeCombinedEfficiencyCategoryGraphs(
        *data_channels_all.efficiency, fit_classes);
    auto channel_ratio_categories = makeCombinedRatioCategoryGraphs(
        *data_channels_all.efficiency, *mc_channels_all.efficiency,
        fit_classes);
    comparison_dir->cd();
    channel_ratio_all->Write();
    for (const auto& category : kFitClassCategories) {
        channel_data_categories[category.key]->Write();
        channel_ratio_categories[category.key]->Write();
    }
    drawCombinedChannelComparison(
        *mc_channels_all.efficiency, *channel_ratio_all,
        channel_data_categories, channel_ratio_categories, canvas_dir,
        output_dir);

    auto chip_ratio_all = makeRatioGraph(
        *data_chips_all.efficiency, *mc_chips_all.efficiency,
        kNLayers * kNChips, "data_over_mc_chip",
        "Data / MC by chip;9*Layer + chip;Ratio", mc.color, mc.marker);
    comparison_dir->cd();
    chip_ratio_all->Write();
    drawCombinedChipComparison(*data_chips_all.efficiency,
                               *mc_chips_all.efficiency, *chip_ratio_all,
                               canvas_dir, output_dir);

    auto summary_ratio = makeRatioGraph(
        *data_summary, *mc_summary, kNLayers, "data_over_mc_by_layer",
        "Data / MC by layer;Layer;Ratio", mc.color, mc.marker);
    comparison_dir->cd();
    summary_ratio->Write();
    comparison_tree.Write();
    chip_comparison_tree.Write();

    TCanvas summary_canvas("comparison_by_layer",
                           "Efficiency comparison by layer", 1400, 900);
    summary_canvas.Divide(1, 2);
    summary_canvas.cd(1);
    gPad->SetGrid();
    styleEfficiencyGraph(*data_summary, kNLayers - 1);
    mc_summary->Draw("P SAME");
    TLegend summary_legend(0.72, 0.15, 0.89, 0.35);
    summary_legend.AddEntry(data_summary.get(), "Data", "lp");
    summary_legend.AddEntry(mc_summary.get(), "MC", "lp");
    summary_legend.Draw();
    summary_canvas.cd(2);
    gPad->SetGrid();
    summary_ratio->Draw("AP");
    summary_ratio->GetXaxis()->SetLimits(-0.5, kNLayers - 0.5);
    summary_ratio->GetYaxis()->SetTitle("Data / MC efficiency");
    TLine unity(-0.5, 1.0, kNLayers - 0.5, 1.0);
    unity.SetLineStyle(2);
    unity.Draw();
    summary_canvas.Update();
    canvas_dir->cd();
    summary_canvas.Write();
    summary_canvas.SaveAs(
        Form("%s/efficiency_comparison_by_layer.pdf", output_dir));
    summary_canvas.SaveAs(
        Form("%s/efficiency_comparison_by_layer.png", output_dir));

    output->Close();

    std::cout << "Output ROOT file: " << output_root << std::endl;
    std::cout << "Layer comparison plots: " << output_dir << std::endl;
    std::cout << "Invalid bins excluded from TEfficiency: data="
              << data.invalid_bins << ", mc=" << mc.invalid_bins << std::endl;
}
