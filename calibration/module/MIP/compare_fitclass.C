/**
 * Compare Data and MC efficiencies grouped by the Data FitClass.
 * MC fit information is not used.
 *
 * Usage:
 *   root -b -q 'compare_fitclass.C("data_fitted.root", "data_nofit.root",
 *                                  "mc_nofit.root", "fitclass_comparison.root",
 *                                  "fitclass_plots")'
 */

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TKey.h>
#include <TLegend.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {

const std::vector<int> kFitClasses = {-1, 0, 1, 2, 3, 4, 5};

const std::map<int, std::string> kFitClassNames = {
    {-1, "LessStatistics"},
    {0, "GoodFit"},
    {1, "FitFailed"},
    {2, "BadChi2"},
    {3, "ParameterAtLimit"},
    {4, "GausSigmaSmallThreshold"},
    {5, "GausSigmaSmallRetry"},
};

std::string classLabel(int fit_class) {
    const auto it = kFitClassNames.find(fit_class);
    return it != kFitClassNames.end() ? it->second
                                      : "Class" + std::to_string(fit_class);
}

std::unique_ptr<TH1D> normalizedClone(const TH1D& source,
                                      const std::string& name, int color,
                                      int style) {
    auto clone = std::unique_ptr<TH1D>(
        static_cast<TH1D*>(source.Clone(name.c_str())));
    clone->SetDirectory(nullptr);
    const double integral = clone->Integral(0, clone->GetNbinsX() + 1);
    if (integral > 0) clone->Scale(1.0 / integral);
    clone->SetLineColor(color);
    clone->SetMarkerColor(color);
    clone->SetLineStyle(style);
    clone->SetLineWidth(2);
    clone->SetStats(false);
    return clone;
}

}  // namespace

namespace fitclass_efficiency {

struct Counts {
    Long64_t entries = 0;
    Long64_t ntrack = -1;
    bool has_entries = false;
    bool has_ntrack = false;
};

struct DataChannel {
    int fit_class = -999;
};

bool valid(const Counts& counts) {
    return counts.has_entries && counts.has_ntrack && counts.ntrack > 0 &&
           counts.entries >= 0 && counts.entries <= counts.ntrack;
}

double getEfficiency(const Counts& counts) {
    return valid(counts) ? static_cast<double>(counts.entries) / counts.ntrack
                         : -1.0;
}

bool loadDataFitClasses(const char* file_name,
                        std::map<int, DataChannel>& channels) {
    std::unique_ptr<TFile> file(TFile::Open(file_name, "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open Data fitted file: " << file_name
                  << std::endl;
        return false;
    }
    TTree* tree = nullptr;
    file->GetObject("mip_fit_results", tree);
    if (!tree) {
        std::cerr << "Error: mip_fit_results is missing in " << file_name
                  << std::endl;
        return false;
    }
    const std::vector<std::string> branches = {"cellid", "FitClass"};
    for (const auto& branch : branches) {
        if (!tree->GetBranch(branch.c_str())) {
            std::cerr << "Error: Data branch " << branch << " is missing"
                      << std::endl;
            return false;
        }
    }

    int cellid = -1;
    int fit_class = -999;
    tree->SetBranchAddress("cellid", &cellid);
    tree->SetBranchAddress("FitClass", &fit_class);
    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
        tree->GetEntry(entry);
        channels[cellid].fit_class = fit_class;
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
                counts[cellid].entries =
                    static_cast<Long64_t>(static_cast<TH1*>(obj)->GetEntries());
                counts[cellid].has_entries = true;
            }
        } else if (obj->InheritsFrom(TParameter<int>::Class())) {
            const std::string name = obj->GetName();
            const std::string prefix = "Ntracks_pass_through_channel_";
            if (name.find(prefix) == 0) {
                const int cellid = std::stoi(name.substr(prefix.size()));
                counts[cellid].ntrack =
                    static_cast<TParameter<int>*>(obj)->GetVal();
                counts[cellid].has_ntrack = true;
            }
        }
    }
}

bool loadCounts(const char* file_name, const char* sample,
                std::map<int, Counts>& counts) {
    std::unique_ptr<TFile> file(TFile::Open(file_name, "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open " << sample
                  << " efficiency file: " << file_name
                  << std::endl;
        return false;
    }
    TDirectory* mip_dir = file->GetDirectory("MIP");
    collectCounts(mip_dir ? mip_dir : file.get(), counts);
    return true;
}

using EffHists = std::map<int, std::unique_ptr<TH1D>>;

struct Category {
    int key;
    const char* label;
    int color;
};

const std::vector<Category> kCategories = {
    {0, "FitClass 0", kBlack},
    {4, "FitClass 4", kRed + 1},
    {5, "FitClass 5", kBlue + 1},
    {99, "Other FitClass", kGreen + 2},
};

EffHists makeHists(const char* sample) {
    EffHists histograms;
    for (int fit_class : kFitClasses) {
        auto histogram = std::make_unique<TH1D>(
            Form("%s_efficiency_dataFitClass%d", sample, fit_class),
            Form("%s efficiency, Data FitClass %d;Efficiency;Channels", sample,
                 fit_class),
            105, 0.0, 1.05);
        histogram->SetDirectory(nullptr);
        histograms[fit_class] = std::move(histogram);
    }
    return histograms;
}

EffHists makeRatioHists() {
    EffHists histograms;
    for (int fit_class : kFitClasses) {
        auto histogram = std::make_unique<TH1D>(
            Form("data_over_mc_dataFitClass%d", fit_class),
            Form("Data / MC, Data FitClass %d;Data efficiency / MC efficiency;"
                 "Channels",
                 fit_class),
            150, 0.0, 1.5);
        histogram->SetDirectory(nullptr);
        histograms[fit_class] = std::move(histogram);
    }
    return histograms;
}

EffHists makeCategoryHists(const EffHists& source, const char* prefix) {
    EffHists categories;
    for (const auto& category : kCategories) {
        auto histogram = std::unique_ptr<TH1D>(static_cast<TH1D*>(
            source.at(0)->Clone(Form("%s_category%d", prefix, category.key))));
        histogram->Reset();
        histogram->SetDirectory(nullptr);
        histogram->SetTitle(category.label);
        categories[category.key] = std::move(histogram);
    }

    for (int fit_class : kFitClasses) {
        const int category_key =
            fit_class == 0 || fit_class == 4 || fit_class == 5 ? fit_class : 99;
        categories[category_key]->Add(source.at(fit_class).get());
    }
    return categories;
}

void draw(const EffHists& data_hists, const EffHists& mc_hists,
          TDirectory* canvas_dir, const std::string& output_dir) {
    TCanvas canvas("efficiency_by_data_fitclass",
                   "Data and MC efficiency by Data FitClass", 1800, 1100);
    canvas.Divide(4, 2);
    std::vector<std::unique_ptr<TH1D>> drawn;
    std::vector<std::unique_ptr<TLegend>> legends;

    for (size_t index = 0; index < kFitClasses.size(); ++index) {
        const int fit_class = kFitClasses[index];
        canvas.cd(index + 1);
        gPad->SetGrid();
        gPad->SetLogy();
        auto data = normalizedClone(*data_hists.at(fit_class),
                                    Form("draw_data_eff_class%d", fit_class),
                                    kBlack, 1);
        auto mc = normalizedClone(*mc_hists.at(fit_class),
                                  Form("draw_mc_eff_class%d", fit_class),
                                  kRed + 1, 1);
        const double maximum = std::max(data->GetMaximum(), mc->GetMaximum());
        data->SetMaximum(maximum > 0 ? maximum * 5.0 : 1.0);
        data->SetMinimum(1e-5);
        data->SetTitle(
            Form("Data FitClass %d: %s;Efficiency;Normalized channels",
                 fit_class, classLabel(fit_class).c_str()));
        data->Draw("HIST");
        mc->Draw("HIST SAME");
        auto legend = std::make_unique<TLegend>(0.50, 0.69, 0.89, 0.89);
        legend->AddEntry(
            data.get(),
            Form("Data (N=%.0f)", data_hists.at(fit_class)->GetEntries()), "l");
        legend->AddEntry(
            mc.get(), Form("MC (N=%.0f)", mc_hists.at(fit_class)->GetEntries()),
            "l");
        legend->Draw();
        drawn.push_back(std::move(data));
        drawn.push_back(std::move(mc));
        legends.push_back(std::move(legend));
    }
    canvas.Update();
    canvas_dir->cd();
    canvas.Write();
    canvas.SaveAs(Form("%s/efficiency_by_data_fitclass.pdf", output_dir.c_str()));
    canvas.SaveAs(Form("%s/efficiency_by_data_fitclass.png", output_dir.c_str()));
}

void drawRatio(const EffHists& ratio_hists, TDirectory* canvas_dir,
               const std::string& output_dir) {
    TCanvas canvas("data_over_mc_by_data_fitclass",
                   "Data / MC efficiency by Data FitClass", 1800, 1100);
    canvas.Divide(4, 2);
    std::vector<std::unique_ptr<TH1D>> drawn;

    for (size_t index = 0; index < kFitClasses.size(); ++index) {
        const int fit_class = kFitClasses[index];
        canvas.cd(index + 1);
        gPad->SetGrid();
        gPad->SetLogy();
        auto ratio = normalizedClone(
            *ratio_hists.at(fit_class), Form("draw_ratio_class%d", fit_class),
            kBlue + 1, 1);
        ratio->SetMinimum(1e-5);
        ratio->SetMaximum(ratio->GetMaximum() > 0 ? ratio->GetMaximum() * 5.0
                                                  : 1.0);
        ratio->SetTitle(
            Form("Data FitClass %d: %s;Data efficiency / MC efficiency;"
                 "Normalized channels",
                 fit_class, classLabel(fit_class).c_str()));
        ratio->Draw("HIST");
        drawn.push_back(std::move(ratio));
    }
    canvas.Update();
    canvas_dir->cd();
    canvas.Write();
    canvas.SaveAs(
        Form("%s/data_over_mc_by_data_fitclass.pdf", output_dir.c_str()));
    canvas.SaveAs(
        Form("%s/data_over_mc_by_data_fitclass.png", output_dir.c_str()));
}

void drawCategoryOverlay(const EffHists& categories, const char* canvas_name,
                         const char* title, const char* x_title,
                         TDirectory* canvas_dir,
                         const std::string& output_dir) {
    TCanvas canvas(canvas_name, title, 1200, 850);
    gPad->SetGrid();
    gPad->SetLogy();

    std::vector<std::unique_ptr<TH1D>> drawn;
    TLegend legend(0.62, 0.65, 0.89, 0.89);
    double maximum = 0.0;
    for (const auto& category : kCategories) {
        auto histogram = normalizedClone(
            *categories.at(category.key),
            Form("draw_%s_category%d", canvas_name, category.key),
            category.color, 1);
        maximum = std::max(maximum, histogram->GetMaximum());
        drawn.push_back(std::move(histogram));
    }

    for (size_t index = 0; index < kCategories.size(); ++index) {
        auto& histogram = drawn[index];
        histogram->SetTitle(Form("%s;%s;Normalized channels", title, x_title));
        histogram->SetMinimum(1e-5);
        histogram->SetMaximum(maximum > 0 ? maximum * 5.0 : 1.0);
        histogram->Draw(index == 0 ? "HIST" : "HIST SAME");
        legend.AddEntry(
            histogram.get(),
            Form("%s (N=%.0f)", kCategories[index].label,
                 categories.at(kCategories[index].key)->GetEntries()),
            "l");
    }
    legend.Draw();
    canvas.Update();
    canvas_dir->cd();
    canvas.Write();
    canvas.SaveAs(Form("%s/%s.pdf", output_dir.c_str(), canvas_name));
    canvas.SaveAs(Form("%s/%s.png", output_dir.c_str(), canvas_name));
}

}  // namespace fitclass_efficiency

void compare_fitclass(
    const char* data_fitted_file =
        "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/"
        "out_calibration/mip_calib/21987-22133/mip_fitted_noplot_7.root",
    const char* data_efficiency_file =
        "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/"
        "out_calibration/mip_calib/21987-22133/mip_neighborcheck_nofit.root",
    const char* mc_nofit_file =
        "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/"
        "out/SimMuon/mip_neighborcheck_nofit.root",
    const char* output_root = "fitclass_comparison.root",
    const char* output_dir = "fitclass_plots") {
    using namespace fitclass_efficiency;
    gStyle->SetOptStat(0);
    gSystem->mkdir(output_dir, true);

    std::map<int, DataChannel> data_channels;
    std::map<int, Counts> data_counts;
    std::map<int, Counts> mc_counts;
    if (!loadDataFitClasses(data_fitted_file, data_channels) ||
        !loadCounts(data_efficiency_file, "Data", data_counts) ||
        !loadCounts(mc_nofit_file, "MC", mc_counts)) {
        return;
    }
    auto data_hists = makeHists("data");
    auto mc_hists = makeHists("mc");
    auto ratio_hists = makeRatioHists();

    std::unique_ptr<TFile> output(TFile::Open(output_root, "RECREATE"));
    TDirectory* data_dir = output->mkdir("data");
    TDirectory* mc_dir = output->mkdir("mc");
    TDirectory* comparison_dir = output->mkdir("comparison");
    TDirectory* canvas_dir = output->mkdir("canvas");

    int cellid = -1;
    int data_fit_class = -999;
    Long64_t data_entries = 0, data_ntrack = -1;
    Long64_t mc_entries = 0, mc_ntrack = -1;
    double data_efficiency = -1.0, mc_efficiency = -1.0, data_over_mc = -1.0;
    TTree tree("fitclass_efficiency_comparison",
               "Efficiency comparison grouped by Data FitClass");
    tree.SetDirectory(comparison_dir);
    tree.Branch("cellid", &cellid);
    tree.Branch("data_fit_class", &data_fit_class);
    tree.Branch("data_entries", &data_entries);
    tree.Branch("data_ntrack", &data_ntrack);
    tree.Branch("mc_entries", &mc_entries);
    tree.Branch("mc_ntrack", &mc_ntrack);
    tree.Branch("data_efficiency", &data_efficiency);
    tree.Branch("mc_efficiency", &mc_efficiency);
    tree.Branch("data_over_mc", &data_over_mc);

    int missing_data = 0;
    int missing_mc = 0;
    for (const auto& item : data_channels) {
        cellid = item.first;
        data_fit_class = item.second.fit_class;
        if (kFitClassNames.find(data_fit_class) == kFitClassNames.end()) continue;
        const auto data_it = data_counts.find(cellid);
        if (data_it == data_counts.end()) {
            data_entries = 0;
            data_ntrack = -1;
            data_efficiency = -1.0;
            ++missing_data;
        } else {
            data_entries = data_it->second.entries;
            data_ntrack = data_it->second.ntrack;
            data_efficiency = getEfficiency(data_it->second);
            if (data_efficiency >= 0) {
                data_hists[data_fit_class]->Fill(data_efficiency);
            }
        }

        const auto mc_it = mc_counts.find(cellid);
        if (mc_it == mc_counts.end()) {
            mc_entries = 0;
            mc_ntrack = -1;
            mc_efficiency = -1.0;
            ++missing_mc;
        } else {
            mc_entries = mc_it->second.entries;
            mc_ntrack = mc_it->second.ntrack;
            mc_efficiency = getEfficiency(mc_it->second);
            if (mc_efficiency >= 0) mc_hists[data_fit_class]->Fill(mc_efficiency);
        }
        data_over_mc = data_efficiency > 0 && mc_efficiency > 0
                           ? data_efficiency / mc_efficiency
                           : -1.0;
        if (data_over_mc > 0) ratio_hists[data_fit_class]->Fill(data_over_mc);
        tree.Fill();
    }

    data_dir->cd();
    for (int fit_class : kFitClasses) data_hists[fit_class]->Write();
    mc_dir->cd();
    for (int fit_class : kFitClasses) mc_hists[fit_class]->Write();
    comparison_dir->cd();
    for (int fit_class : kFitClasses) ratio_hists[fit_class]->Write();
    tree.Write();
    auto data_categories = makeCategoryHists(data_hists, "data_efficiency");
    auto mc_categories = makeCategoryHists(mc_hists, "mc_efficiency");
    auto ratio_categories = makeCategoryHists(ratio_hists, "data_over_mc");
    data_dir->cd();
    for (const auto& category : kCategories) {
        data_categories[category.key]->Write();
    }
    mc_dir->cd();
    for (const auto& category : kCategories) {
        mc_categories[category.key]->Write();
    }
    comparison_dir->cd();
    for (const auto& category : kCategories) {
        ratio_categories[category.key]->Write();
    }
    draw(data_hists, mc_hists, canvas_dir, output_dir);
    drawRatio(ratio_hists, canvas_dir, output_dir);
    drawCategoryOverlay(data_categories, "data_efficiency_categories",
                        "Data efficiency by FitClass category", "Efficiency",
                        canvas_dir, output_dir);
    drawCategoryOverlay(mc_categories, "mc_efficiency_categories",
                        "MC efficiency by Data FitClass category", "Efficiency",
                        canvas_dir, output_dir);
    drawCategoryOverlay(ratio_categories, "data_over_mc_categories",
                        "Data / MC by Data FitClass category",
                        "Data efficiency / MC efficiency", canvas_dir,
                        output_dir);
    output->Close();

    std::cout << "Data channels=" << data_channels.size()
              << ", Data efficiency channels=" << data_counts.size()
              << ", MC channels=" << mc_counts.size()
              << ", missing Data efficiency cell IDs=" << missing_data
              << ", missing MC cell IDs=" << missing_mc << std::endl;
    for (int fit_class : kFitClasses) {
        std::cout << "  Data FitClass " << fit_class << ": Data="
                  << data_hists[fit_class]->GetEntries()
                  << ", MC=" << mc_hists[fit_class]->GetEntries() << std::endl;
    }
}
