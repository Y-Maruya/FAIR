#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TPad.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>

#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {
constexpr int kNLayer = 40;
constexpr int kNChip = 9;
constexpr int kNChannel = 36;
constexpr int kNLayerChip = kNLayer * kNChip;
constexpr int kNBIN_XY = 18;
const double kXYMin = -40.3 * 18.0 / 2.0;
const double kXYMax = +40.3 * 18.0 / 2.0;

constexpr std::array<double, kNChannel> kPosX = {
    100.2411, 100.2411, 100.2411,  59.94146,  59.94146,  59.94146,
     19.64182,  19.64182,  19.64182,  19.64182,  59.94146, 100.2411,
    100.2411,  59.94146,  19.64182, 100.2411,  59.94146,  19.64182,
    -20.65782, -60.95746, -101.2571, -20.65782, -60.95746, -101.2571,
    -101.2571, -60.95746, -20.65782, -20.65782, -20.65782, -20.65782,
    -60.95746, -60.95746, -60.95746, -101.2571, -101.2571, -101.2571};

constexpr std::array<double, kNChannel> kPosY = {
    141.04874, 181.34838, 221.64802, 141.04874, 181.34838, 221.64802,
    141.04874, 181.34838, 221.64802, 261.94766, 261.94766, 261.94766,
    302.2473,  302.2473,  302.2473,  342.54694, 342.54694, 342.54694,
    342.54694, 342.54694, 342.54694, 302.2473,  302.2473,  302.2473,
    261.94766, 261.94766, 261.94766, 221.64802, 181.34838, 141.04874,
    221.64802, 181.34838, 141.04874, 221.64802, 181.34838, 141.04874};

inline double pos_x(int channel, int chip) {
  constexpr double chip_dis_y = 241.8;
  int chip_mod = chip % 3;
  int ch = channel;
  if (chip_mod != 0) {
    if (ch == 2) ch = 0;
    else if (ch == 0) ch = 2;
    if (ch == 33) ch = 35;
    else if (ch == 35) ch = 33;
  }
  return -(kPosY[ch] - chip_mod * chip_dis_y);
}

inline double pos_y(int channel, int chip) {
  constexpr double hbu_x = 239.3;
  int hbu_id = chip / 3;
  return -(-kPosX[channel] + (hbu_id - 1) * hbu_x);
}

void draw_label(const char* text, double x = 0.12, double y = 0.92, double size = 0.07) {
  TLatex l;
  l.SetNDC(true);
  l.SetTextSize(size);
  l.DrawLatex(x, y, text);
}

void draw_layerchip_boundaries(int bins_per_layer, double y_min, double y_max) {
  for (int layer = 1; layer < kNLayer; ++layer) {
    const double x = layer * bins_per_layer - 0.5;
    TLine* line = new TLine(x, y_min, x, y_max);
    line->SetLineColor(kRed + 1);
    line->SetLineStyle(2);
    line->Draw("SAME");
  }

  TLatex latex;
  latex.SetTextFont(42);
  latex.SetTextSize(0.025);
  latex.SetTextColor(kBlack);
  for (int layer = 0; layer < kNLayer; ++layer) {
    const double x = layer * bins_per_layer + bins_per_layer * 0.25;
    latex.DrawLatex(x, 1+ 0.02 * (y_max - y_min), Form("L%d", layer));
  }
}

void draw_chip_boundaries(double y_min, double y_max) {
  for (int chip = 1; chip < kNChip; ++chip) {
    const double x = chip - 0.5;
    TLine* line = new TLine(x, y_min, x, y_max);
    line->SetLineColor(kRed + 1);
    line->SetLineStyle(2);
    line->Draw("SAME");
  }

  TLatex latex;
  latex.SetTextFont(42);
  latex.SetTextSize(0.03);
  latex.SetTextColor(kBlack);
  for (int chip = 0; chip < kNChip; ++chip) {
    latex.DrawLatex(chip - 0.1, y_max + 0.05 * (y_max - y_min), Form("C%d", chip));
  }
}
}

void draw_mip_eff_th2d(
    const char* input_file = "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out_calibration/mip_calib/21987-22163/mip_neighborcheck_efficiency.root",
    const char* output_dir = "./mip_eff_th2d_plots") {
  gStyle->SetOptStat(0);
  gSystem->mkdir(output_dir, true);

  std::unique_ptr<TFile> fin(TFile::Open(input_file, "READ"));
  if (!fin || fin->IsZombie()) {
    std::cerr << "[ERROR] Cannot open input file: " << input_file << std::endl;
    return;
  }

  TTree* t = dynamic_cast<TTree*>(fin->Get("mip_efficiency"));
  if (!t) {
    std::cerr << "[ERROR] TTree 'mip_efficiency' not found in: " << input_file << std::endl;
    return;
  }

  Int_t layer = -1;
  Int_t chip = -1;
  Int_t channel = -1;
  Long64_t entries = 0;
  Int_t ntrack = 0;
  Double_t efficiency = 0.0;

  t->SetBranchAddress("layer", &layer);
  t->SetBranchAddress("chip", &chip);
  t->SetBranchAddress("channel", &channel);
  t->SetBranchAddress("entries", &entries);
  t->SetBranchAddress("ntrack_pass_through_channel", &ntrack);
  t->SetBranchAddress("efficiency", &efficiency);

  std::vector<std::unique_ptr<TH2D>> h_layer;
  h_layer.reserve(kNLayer);
  for (int l = 0; l < kNLayer; ++l) {
    auto h = std::make_unique<TH2D>(
        Form("h2_eff_layer%02d", l),
        Form("Efficiency map L%02d;x [mm];y [mm]", l),
        kNBIN_XY, kXYMin, kXYMax,
        kNBIN_XY, kXYMin, kXYMax);
    h->GetZaxis()->SetTitle("Efficiency");
    h->GetZaxis()->SetRangeUser(0.0, 1.0);
    h_layer.push_back(std::move(h));
  }

  std::vector<std::unique_ptr<TH2D>> h_chip_xy_avg;
  h_chip_xy_avg.reserve(kNChip);
  std::vector<std::unique_ptr<TH2D>> h_chip_xy_sum;
  h_chip_xy_sum.reserve(kNChip);
  std::vector<std::unique_ptr<TH2D>> h_chip_xy_cnt;
  h_chip_xy_cnt.reserve(kNChip);
  for (int c = 0; c < kNChip; ++c) {
    auto h_avg = std::make_unique<TH2D>(
        Form("h2_eff_chip%02d_xyavg", c),
        Form("Efficiency XY map (layer-avg) Chip %d;x [mm];y [mm]", c),
        kNBIN_XY, kXYMin, kXYMax,
        kNBIN_XY, kXYMin, kXYMax);
    h_avg->GetZaxis()->SetTitle("Efficiency");
    h_avg->GetZaxis()->SetRangeUser(0.0, 1.0);

    auto h_sum = std::make_unique<TH2D>(
        Form("h2_eff_chip%02d_xysum", c),
        "",
        kNBIN_XY, kXYMin, kXYMax,
        kNBIN_XY, kXYMin, kXYMax);
    auto h_cnt = std::make_unique<TH2D>(
        Form("h2_eff_chip%02d_xycnt", c),
        "",
        kNBIN_XY, kXYMin, kXYMax,
        kNBIN_XY, kXYMin, kXYMax);
    h_chip_xy_avg.push_back(std::move(h_avg));
    h_chip_xy_sum.push_back(std::move(h_sum));
    h_chip_xy_cnt.push_back(std::move(h_cnt));
  }

  auto h_eff_layerchip_channel = std::make_unique<TH2D>(
      "h2_eff_layerchip_channel",
      "Efficiency;chip + layer*9;channel",
      kNLayerChip, -0.5, kNLayerChip - 0.5,
      kNChannel, -0.5, kNChannel - 0.5);
  h_eff_layerchip_channel->GetZaxis()->SetTitle("Efficiency");
  h_eff_layerchip_channel->GetZaxis()->SetRangeUser(0.0, 1.0);

  auto h_one_minus_eff_layerchip_channel = std::make_unique<TH2D>(
      "h2_one_minus_eff_layerchip_channel",
      "1 - Efficiency;chip + layer*9;channel",
      kNLayerChip, -0.5, kNLayerChip - 0.5,
      kNChannel, -0.5, kNChannel - 0.5);
  h_one_minus_eff_layerchip_channel->GetZaxis()->SetTitle("1 - Efficiency");
  h_one_minus_eff_layerchip_channel->GetZaxis()->SetRangeUser(1e-4, 1.0);

  auto h_eff_layerchip = std::make_unique<TH1D>(
      "h1_eff_layerchip",
      "Mean efficiency per (layer,chip);chip + layer*9;Efficiency",
      kNLayerChip, -0.5, kNLayerChip - 0.5);
  h_eff_layerchip->GetYaxis()->SetRangeUser(0.0, 1.0);

  std::vector<double> entries_sum_layerchip(kNLayerChip, 0.0);
  std::vector<double> ntrack_sum_layerchip(kNLayerChip, 0.0);
  std::vector<double> entries_sum_chip(kNChip, 0.0);
  std::vector<double> ntrack_sum_chip(kNChip, 0.0);

  const Long64_t n = t->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    if (layer < 0 || layer >= kNLayer) continue;
    if (chip < 0 || chip >= kNChip) continue;
    if (channel < 0 || channel >= kNChannel) continue;

    const double x = pos_x(channel, chip);
    const double y = pos_y(channel, chip);
    const int layerchip = layer * kNChip + chip;
    const double one_minus_eff = std::max(1.0 - efficiency, 1e-6);

    h_layer[layer]->Fill(x, y, efficiency);
    h_chip_xy_sum[chip]->Fill(x, y, efficiency);
    h_chip_xy_cnt[chip]->Fill(x, y, 1.0);
    h_eff_layerchip_channel->SetBinContent(layerchip + 1, channel + 1,
                                           efficiency);
    h_one_minus_eff_layerchip_channel->SetBinContent(layerchip + 1,
                                                     channel + 1,
                                                     one_minus_eff);
    entries_sum_layerchip[layerchip] += static_cast<double>(entries);
    ntrack_sum_layerchip[layerchip] += static_cast<double>(ntrack);
    entries_sum_chip[chip] += static_cast<double>(entries);
    ntrack_sum_chip[chip] += static_cast<double>(ntrack);
  }

  for (int c = 0; c < kNChip; ++c) {
    for (int bx = 1; bx <= h_chip_xy_avg[c]->GetNbinsX(); ++bx) {
      for (int by = 1; by <= h_chip_xy_avg[c]->GetNbinsY(); ++by) {
        const double cnt = h_chip_xy_cnt[c]->GetBinContent(bx, by);
        const double sum = h_chip_xy_sum[c]->GetBinContent(bx, by);
        if (cnt > 0.0) {
          h_chip_xy_avg[c]->SetBinContent(bx, by, sum / cnt);
        }
      }
    }
  }

  for (int layerchip = 0; layerchip < kNLayerChip; ++layerchip) {
    if (ntrack_sum_layerchip[layerchip] > 0.0) {
      h_eff_layerchip->SetBinContent(layerchip + 1,
                                     entries_sum_layerchip[layerchip] /
                                         ntrack_sum_layerchip[layerchip]);
    }
  }

  std::unique_ptr<TFile> fout(TFile::Open(Form("%s/mip_eff_th2d.root", output_dir), "RECREATE"));
  if (!fout || fout->IsZombie()) {
    std::cerr << "[ERROR] Cannot create output ROOT file in: " << output_dir << std::endl;
    return;
  }

  for (int l = 0; l < kNLayer; ++l) {
    TCanvas c(Form("c_layer%02d", l), Form("Layer %02d", l), 900, 700);
    h_layer[l]->Draw("COLZ");
    draw_label(Form("Layer %d", l));
    c.SaveAs(Form("%s/eff_layer%02d_xy.png", output_dir, l));
    h_layer[l]->Write();
    c.Write();
  }

  TCanvas c_all("c_all_layers", "All layers", 5600, 4200);
  c_all.Divide(7, 6, 0.001, 0.001);
  for (int l = 0; l < kNLayer; ++l) {
    c_all.cd(l + 1);
    gPad->SetMargin(0.08, 0.14, 0.08, 0.10);
    h_layer[l]->Draw("COLZ");
    draw_label(Form("L%d", l), 0.10, 0.90, 0.08);
  }
  c_all.SaveAs(Form("%s/eff_all_layers_xy.png", output_dir));
  c_all.Write();

  TCanvas c_chip_row("c_chip_row", "Chips side-by-side (XY, layer-avg)", 5400, 750);
  c_chip_row.Divide(kNChip, 1, 0.001, 0.001);
  for (int c = 0; c < kNChip; ++c) {
    c_chip_row.cd(c + 1);
    gPad->SetMargin(0.08, 0.16, 0.10, 0.08);
    h_chip_xy_avg[c]->Draw("COLZ");
    draw_label(Form("Chip %d", c), 0.12, 0.90, 0.09);
    h_chip_xy_avg[c]->Write();
  }
  c_chip_row.SaveAs(Form("%s/eff_chip_side_by_side_xy_layeravg.png", output_dir));
  c_chip_row.Write();

    TCanvas c_eff_layerchip("c_eff_layerchip",
                "Efficiency vs chip+layer*9", 3800, 700);
    h_eff_layerchip->SetLineColor(kBlue + 1);
    h_eff_layerchip->SetMarkerColor(kBlue + 1);
    h_eff_layerchip->SetMarkerStyle(20);
    h_eff_layerchip->SetMarkerSize(0.7);
    h_eff_layerchip->Draw("EP");
    c_eff_layerchip.SaveAs(Form("%s/eff_chip_plus_layer9.png", output_dir));
    h_eff_layerchip->Write();
    c_eff_layerchip.Write();

    auto g_layerchip_eff = std::make_unique<TGraphErrors>();
    g_layerchip_eff->SetName("g_eff_layerchip_entries_over_ntrack");
    g_layerchip_eff->SetTitle(
        "Efficiency from entries/ntrack;chip + layer*9;Efficiency");
    g_layerchip_eff->SetLineColor(kMagenta + 2);
    g_layerchip_eff->SetMarkerColor(kMagenta + 2);
    g_layerchip_eff->SetMarkerStyle(20);
    g_layerchip_eff->SetMarkerSize(0.6);

    for (int layerchip = 0; layerchip < kNLayerChip; ++layerchip) {
      const double n = ntrack_sum_layerchip[layerchip];
      if (n <= 0.0) continue;
      const double p = entries_sum_layerchip[layerchip] / n;
      const double p_clamped = std::max(0.0, std::min(1.0, p));
      const double err = std::sqrt(p_clamped * (1.0 - p_clamped) / n);
      const int ip = g_layerchip_eff->GetN();
      g_layerchip_eff->SetPoint(ip, layerchip, p_clamped);
      g_layerchip_eff->SetPointError(ip, 0.0, err);
    }

    TCanvas c_layerchip_graph("c_layerchip_graph",
                              "LayerChip TGraph from entries/ntrack", 3800,
                              700);
    c_layerchip_graph.SetGrid();
    g_layerchip_eff->Draw("AP");
    g_layerchip_eff->GetYaxis()->SetRangeUser(0.0, 1.05);
    g_layerchip_eff->GetXaxis()->SetRangeUser(0.0-0.5, kNLayerChip+0.5);
    draw_layerchip_boundaries(kNChip, 0.0, 1.05);
    gPad->Update();
    c_layerchip_graph.SetGrid(0);
    c_layerchip_graph.SaveAs(
        Form("%s/eff_layerchip_tgraph_entries_over_ntrack.png", output_dir));
    g_layerchip_eff->Write();
    c_layerchip_graph.Write();

    auto g_chip_eff = std::make_unique<TGraphErrors>();
    g_chip_eff->SetName("g_chip_eff_entries_over_ntrack");
    g_chip_eff->SetTitle("Chip-by-chip efficiency from entries/ntrack;chip;Efficiency");
    g_chip_eff->SetLineColor(kBlue + 2);
    g_chip_eff->SetMarkerColor(kBlue + 2);
    g_chip_eff->SetMarkerStyle(20);
    g_chip_eff->SetMarkerSize(1.0);
    for (int c = 0; c < kNChip; ++c) {
      const double n = ntrack_sum_chip[c];
      if (n <= 0.0) continue;
      const double p = entries_sum_chip[c] / n;
      const double p_clamped = std::max(0.0, std::min(1.0, p));
      const double err = std::sqrt(p_clamped * (1.0 - p_clamped) / n);
      const int ip = g_chip_eff->GetN();
      g_chip_eff->SetPoint(ip, c, p_clamped);
      g_chip_eff->SetPointError(ip, 0.0, err);
    }

    TCanvas c_chip_graph("c_chip_graph", "Chip-by-chip TGraph", 1200, 700);
    c_chip_graph.SetGrid();
    g_chip_eff->Draw("AP");
    g_chip_eff->GetYaxis()->SetRangeUser(0.0, 1.05);
    g_chip_eff->GetXaxis()->SetRangeUser(-0.5, kNChip - 0.5);
    draw_chip_boundaries(0.0, 1.05);
    c_chip_graph.SaveAs(Form("%s/eff_chip_tgraph_entries_over_ntrack.png", output_dir));
    g_chip_eff->Write();
    c_chip_graph.Write();

    TCanvas c_eff_lc_ch("c_eff_layerchip_channel",
              "Efficiency map (layerchip-channel)", 3800, 700);
    h_eff_layerchip_channel->Draw("COLZ");
    c_eff_lc_ch.SaveAs(
      Form("%s/eff_chip_plus_layer9_channel_th2d.png", output_dir));
    h_eff_layerchip_channel->Write();
    c_eff_lc_ch.Write();

    TCanvas c_one_minus_eff_log("c_one_minus_eff_log",
                  "1-efficiency map logZ", 3800, 700);
    c_one_minus_eff_log.SetLogz();
    h_one_minus_eff_layerchip_channel->Draw("COLZ");
    c_one_minus_eff_log.SaveAs(
      Form("%s/one_minus_eff_chip_plus_layer9_channel_logz.png",
         output_dir));
    h_one_minus_eff_layerchip_channel->Write();
    c_one_minus_eff_log.Write();

  fout->Close();

  std::cout << "[INFO] Wrote outputs to: " << output_dir << std::endl;
  std::cout << "[INFO] ROOT file: " << output_dir << "/mip_eff_th2d.root" << std::endl;
}
