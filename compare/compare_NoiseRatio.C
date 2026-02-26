#include <vector>
#include <string>
#include <algorithm>

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"
#include "TError.h"
#include "TString.h"

namespace {
struct Item {
	int run;      // set to -1 if you use explicit filePath
	std::string label;
	std::string filePath; // optional (overrides filePattern)
};

std::string makeFilePath(const Item& item, const char* filePattern)
{
	if (!item.filePath.empty()) return item.filePath;
	if (item.run < 0) return std::string(filePattern);
	return std::string(Form(filePattern, item.run));
}

std::string makeHistPath(const Item& item)
{
	// Directory: Layer12/Chip_0/HG/
	// Histogram: hittag_0_ADC_hg_12_0_0
	return std::string(Form(
		"NoiseMap/hNoiseRatio_Layer_Normalized"));
}
} // namespace

void compare_NoiseRatio()
{
	// --- User settings ---
	// If you want to compare across runs in this folder structure, keep this pattern.
	// Examples:
	//   "run%d/histogram.root"
	//   "run%d/hittag0_checker_histograms.root"
	// If you want a single file in the current directory, set e.g. "histogram.root" and set run=-1 in items.
	const char* filePattern = "../out/run%d/noise_hits.root";

	// (run, layer, chip, channel) を何個でも追加してOK
	// - run: filePatternに %d がある場合に使う（例: run21774/...）
	// - filePath: 個別にファイルを指定したい場合はここに入れる（filePatternより優先）
	const std::vector<Item> items = {
		// Example:
		// {21774, "50 Hz", ""},
        // {21773, "100 Hz", ""},
        // {21739, "300 Hz", ""},
        // {21729, "500 Hz", ""},
		{21936, "100 Hz (new)", ""},
		{21937, "300 Hz (new)", ""},
		{21938, "50 Hz (new)", ""}
	};

	// legend位置
	TLegend* leg = new TLegend(0.58, 0.20, 0.88, 0.38);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);

	// 色（足りなければ回す）
	const std::vector<int> colors = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 1, kBlack};

	// optional
	const bool logy = true;
	const char* outPng = "compare_NoiseRatio.png";
	const bool normalizeToUnitArea = false; // set to true if you want to compare shapes regardless of total counts
	const int referenceItemIndex = 0; // index in items used as denominator for ratio panel
	const double ratioMin = 0.5;
	const double ratioMax = 1.5;
	// --- Basic checks ---
	if (items.empty()) {
		::Error("compare_NoiseRatio", "items is empty. Please add at least one Item.");
		return;
	}

	gStyle->SetOptStat(0);

	TCanvas* c = new TCanvas("c_compare_NoiseRatio", "compare: hNoiseRatio_Layer_Normalized", 1000, 900);
	TPad* padTop = new TPad("padTop", "", 0.0, 0.30, 1.0, 1.0);
	TPad* padBottom = new TPad("padBottom", "", 0.0, 0.0, 1.0, 0.30);
	padTop->SetBottomMargin(0.02);
	padBottom->SetTopMargin(0.03);
	padBottom->SetBottomMargin(0.30);
	padTop->Draw();
	padBottom->Draw();

	bool firstDraw = true;
	double globalMax = 0.0;

	// First pass: load and detach all hists
	std::vector<TH1*> hists;
	hists.reserve(items.size());
	std::vector<std::string> loadedLabels;
	loadedLabels.reserve(items.size());
	std::vector<int> loadedItemIndices;
	loadedItemIndices.reserve(items.size());

	for (size_t i = 0; i < items.size(); ++i) {
		const Item& item = items[i];
		const std::string fpath = makeFilePath(item, filePattern);
		const std::string hpath = makeHistPath(item);

		TFile* f = TFile::Open(fpath.c_str(), "READ");
		if (!f || f->IsZombie()) {
			::Warning("compare_NoiseRatio", "Cannot open: %s (skip)", fpath.c_str());
			if (f) f->Close();
			continue;
		}

		TH1* h = dynamic_cast<TH1*>(f->Get(hpath.c_str()));
		if (!h) {
			::Warning("compare_NoiseRatio", "Histogram not found: %s in %s (skip)", hpath.c_str(), fpath.c_str());
			f->Close();
			continue;
		}
        h->GetXaxis()->SetRangeUser(300, 600); // zoom in
		h->SetDirectory(nullptr);
		f->Close();

		if (normalizeToUnitArea) {
			const double integral = h->Integral();
			if (integral > 0) h->Scale(1.0 / integral);
		}

		const int col = colors.empty() ? int(kBlack) : colors[i % colors.size()];
		h->SetLineColor(col);
		h->SetMarkerColor(col);
		h->SetLineWidth(2);

		globalMax = std::max(globalMax, h->GetMaximum());
		hists.push_back(h);
		loadedItemIndices.push_back(static_cast<int>(i));

		if (!item.label.empty()) {
			loadedLabels.push_back(item.label);
		} else if (item.run >= 0) {
			loadedLabels.push_back(Form("run%d", item.run));
		}
	}

	if (hists.empty()) {
		::Error("compare_NoiseRatio", "No histograms loaded. Check filePattern/items and histogram paths.");
		return;
	}

	int refLoadedIndex = -1;
	for (size_t i = 0; i < loadedItemIndices.size(); ++i) {
		if (loadedItemIndices[i] == referenceItemIndex) {
			refLoadedIndex = static_cast<int>(i);
			break;
		}
	}
	if (refLoadedIndex < 0) {
		::Error("compare_NoiseRatio", "Reference item index %d not loaded. Check items and files.", referenceItemIndex);
		return;
	}

	TH1* hRef = hists[refLoadedIndex];

	padTop->cd();
	if (logy) padTop->SetLogy();

	// Second pass: draw
	for (size_t i = 0; i < hists.size(); ++i) {
		TH1* h = hists[i];
		h->SetTitle(";Layer;Noise Ratio");
        if (logy && !normalizeToUnitArea) {
            // If log scale, set minimum to 0.1 or 1% of max to avoid issues with zero bins
            // h->SetMinimum(0.1);.
            // h->GetYaxis()->SetRangeUser(0.5, globalMax * 1.2);
        } else if (logy && normalizeToUnitArea) {
            // h->GetYaxis()->SetRangeUser(0.00005, globalMax * 1.2);
        } else {
            // h->GetYaxis()->SetRangeUser(0, globalMax * 1.2);
        }

		if (firstDraw) {
            h->SetTitle("hNoiseRatio_Layer_Normalized comparison");
			h->Draw("e");
			firstDraw = false;
		} else {
			h->Draw("e same");
		}
	}

	// legend entries (only the successfully loaded ones, in the same order)
	for (size_t i = 0; i < hists.size(); ++i) {
		const char* label = (i < loadedLabels.size()) ? loadedLabels[i].c_str() : hists[i]->GetName();
		leg->AddEntry(hists[i], label, "l");
	}

	leg->Draw();

	padBottom->cd();
	bool firstRatio = true;
	for (size_t i = 0; i < hists.size(); ++i) {
		TH1* h = hists[i];
		TH1* hRatio = dynamic_cast<TH1*>(h->Clone(Form("%s_ratio", h->GetName())));
		hRatio->SetDirectory(nullptr);
		hRatio->Divide(hRef);
		hRatio->SetTitle(";Layer;Ratio to ref");
		hRatio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
		hRatio->GetYaxis()->SetNdivisions(505);
		hRatio->GetYaxis()->SetTitleSize(0.10);
		hRatio->GetYaxis()->SetTitleOffset(0.45);
		hRatio->GetYaxis()->SetLabelSize(0.08);
		hRatio->GetXaxis()->SetTitleSize(0.10);
		hRatio->GetXaxis()->SetLabelSize(0.08);
		hRatio->GetXaxis()->SetTitleOffset(1.0);

		if (firstRatio) {
			hRatio->Draw("e");
			firstRatio = false;
		} else {
			hRatio->Draw("e same");
		}
	}

	c->Modified();
	c->Update();
	c->SaveAs(outPng);
}
