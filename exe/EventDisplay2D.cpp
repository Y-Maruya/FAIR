#include <TCanvas.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TString.h>
#include <TPad.h>
#include <TLine.h>

#include <iostream>
#include <vector>
#include <string>

#include <yaml-cpp/yaml.h>
#include "common/RunContext.hpp"
#include "common/config/ParseRunConfig.hpp"
#include "IO/reader/RootInput.hpp"
#include "IO/reader/ReaderRegistry.hpp"
#include "common/edm/EDM.hpp"
#include "common/AHCALGeometry.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config.yaml> [ievt]\n";
        return 1;
    }

    YAML::Node config = YAML::LoadFile(argv[1]);
    RunContext ctx;
    ctx.config = parse_run_config(config);

    const std::string fname = ctx.config.output;
    const long long ievt = (argc >= 3) ? std::stoll(argv[2]) : 0;

    FAIR::init_logger("AHCALApp", "eventdisplay/2D.log", spdlog::level::debug);
    RootInput in(fname, "events");

    ReaderRegistry rr;
    rr.register_vector_struct<AHCALRecoHit>("vector<AHCALRecoHit>");
    rr.register_struct<Track>("Track");

    if (!in.read_entry(ievt)) {
        std::cerr << "Event out of range: " << ievt << " (entries=" << in.entries() << ")\n";
        return 1;
    }

    // prefix は Writer の key（例: RecoHits）
    auto reco  = rr.read<std::vector<AHCALRecoHit>>("vector<AHCALRecoHit>", "RecoHits", in);
    auto track = rr.read<Track>("Track", "MuonKFTrack", in);

    std::cout << "Event " << ievt << ": " << reco.size()
              << " RecoHits, Track valid=" << track.valid << "\n";
    std::cout << "  Track: chi2=" << track.chi2 << ", ndof=" << track.ndof
              << ", nInTrackHits=" << track.nInTrackHits
              << ", nOutTrackHits=" << track.nOutTrackHits << "\n";

    if (reco.empty()) {
        std::cerr << "No RecoHits in event " << ievt << "\n";
        return 1;
    }

    // --- style ---
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(100);
    gStyle->SetPalette(kViridis);

    // --- ranges (from geometry) ---
    const double xMin = -AHCALGeometry::x_max ;
    const double xMax = +AHCALGeometry::x_max ;
    const double yMin = -AHCALGeometry::y_max ;
    const double yMax = +AHCALGeometry::y_max ;

    const double zMin = 0.0;
    const double zMax = AHCALGeometry::Pos_Z(AHCALGeometry::Layer_No - 1) + 30.0;

    // --- histograms: weight = Edep(MeV) ---
    TH2D hxy("hxy", Form("RecoHits XY (evt %lld);X [mm];Y [mm]", ievt),
             18, xMin, xMax, 18, yMin, yMax);
    TH2D hxz("hxz", Form("RecoHits XZ (evt %lld);Z [mm];X [mm]", ievt),
             120, zMin, zMax, 18, xMin, xMax);
    TH2D hyz("hyz", Form("RecoHits YZ (evt %lld);Z [mm];Y [mm]", ievt),
             120, zMin, zMax, 18, yMin, yMax);

    for (const auto& hit : reco) {
        const double x = hit.Xpos();
        const double y = hit.Ypos();
        const double z = hit.Zpos();
        if (hit.Nmip <0.5) continue; // skip low-energy hits
        const double w = hit.Edep; // color = deposited energy [MeV]
        // Nmip を使いたいなら: const double w = hit.Nmip;

        hxy.Fill(x, y, w);
        hxz.Fill(z, x, w);
        hyz.Fill(z, y, w);
    }

    hxy.SetMinimum(0.0);
    hxz.SetMinimum(0.0);
    hyz.SetMinimum(0.0);

    auto drawLabel = [&](TPad* pad, const char* title) {
        pad->cd();
        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextSize(0.04);
        lat.DrawLatex(0.12, 0.92, title);
        lat.SetTextSize(0.032);
        lat.DrawLatex(0.12, 0.87, Form("evt=%lld, nhit=%zu", ievt, reco.size()));
        lat.DrawLatex(0.12, 0.83, Form("track.valid=%d, chi2/ndof=%.2f/%d",
                                        (int)track.valid, track.chi2, track.ndof));
        if (track.valid) {
            lat.DrawLatex(0.12, 0.79, Form("state: x=%.1f y=%.1f z=%.1f  tx=%.4f ty=%.4f",
                                           track.x, track.y, track.z, track.tx, track.ty));
        } else {
            lat.DrawLatex(0.12, 0.79, "state: (invalid)");
        }
        lat.DrawLatex(0.12, 0.75, "color = Edep [MeV]");
    };

    // Track extrapolation: x(z)=x0+tx*(z-z0), y(z)=y0+ty*(z-z0)
    auto xAt = [&](double z) { return track.x + track.tx * (z - track.z); };
    auto yAt = [&](double z) { return track.y + track.ty * (z - track.z); };

    TCanvas c("c", "AHCAL EventDisplay2D", 1800, 900);
    c.Divide(2, 2);

    c.cd(1);
    gPad->SetRightMargin(0.14);
    hxy.Draw("COLZ");
    drawLabel((TPad*)gPad, "XY projection");

    c.cd(2);
    gPad->SetRightMargin(0.14);
    hxz.Draw("COLZ");
    drawLabel((TPad*)gPad, "XZ projection");

    c.cd(3);
    gPad->SetRightMargin(0.14);
    hyz.Draw("COLZ");
    drawLabel((TPad*)gPad, "YZ projection");

    // ---- overlay MuonKFTrack (red line) ----
    if (track.valid) {
        const double x1 = xAt(zMin);
        const double y1 = yAt(zMin);
        const double x2 = xAt(zMax);
        const double y2 = yAt(zMax);

        // XY: (x(z),y(z))
        c.cd(1);
        {
            auto* lxy = new TLine(x1, y1, x2, y2);
            lxy->SetLineColor(kRed + 1);
            lxy->SetLineWidth(3);
            lxy->Draw("SAME");
        }

        // XZ: histogram axes are (X, Z)
        c.cd(2);
        {
            auto* lxz = new TLine(x1, zMin, x2, zMax);
            lxz->SetLineColor(kRed + 1);
            lxz->SetLineWidth(3);
            lxz->Draw("SAME");
        }

        // YZ: histogram axes are (Y, Z)
        c.cd(3);
        {
            auto* lyz = new TLine(y1, zMin, y2, zMax);
            lyz->SetLineColor(kRed + 1);
            lyz->SetLineWidth(3);
            lyz->Draw("SAME");
        }
    }

    c.cd(4);
    TLatex info;
    info.SetNDC(true);
    info.SetTextSize(0.035);
    info.DrawLatex(0.10, 0.90, "RecoHit -> (X,Y,Z):");
    info.SetTextSize(0.03);
    info.DrawLatex(0.10, 0.84, "Xpos/Ypos/Zpos use AHCALGeometry");
    info.DrawLatex(0.10, 0.79, "cellID = layer*100000 + asic*10000 + channel");
    info.DrawLatex(0.10, 0.74, "MuonKFTrack: red line = (x,y) + (tx,ty)*(z-track.z)");

    c.Update();

    // const std::string outPdf = Form("eventdisplay/evt%lld.pdf", ievt);
    const std::string outPng = Form("eventdisplay/evt%lld.png", ievt);
    // c.SaveAs(outPdf.c_str());
    c.SaveAs(outPng.c_str());

    // std::cout << "Saved: " << outPdf << " and " << outPng << "\n";
    return 0;
}