#include <TCanvas.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TString.h>
#include <TPad.h>
#include <TLine.h>
#include <TBox.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <yaml-cpp/yaml.h>
#include "common/RunContext.hpp"
#include "common/config/ParseRunConfig.hpp"
#include "IO/reader/RootInput.hpp"
#include "IO/reader/ReaderRegistry.hpp"
#include "common/edm/EDM.hpp"
#include "common/AHCALGeometry.hpp"
    const double xMin = -AHCALGeometry::x_max ;
    const double xMax = +AHCALGeometry::x_max ;
    const double yMin = -AHCALGeometry::y_max ;
    const double yMax = +AHCALGeometry::y_max ;

    const double zMin = -9.0;
    const double zMax = 1191 + 15.0;
    // --- histograms: weight = Edep(MeV) ---
    TH2D hxy_tmp("hxy", "RecoHits XY;X [mm];Y [mm]", 18, xMin, xMax, 18, yMin, yMax);
    TH2D hxz_tmp("hxz", "RecoHits XZ;Z [mm];X [mm]", 135, zMin, zMax, 18, xMin, xMax);
    TH2D hyz_tmp("hyz", "RecoHits YZ;Z [mm];Y [mm]", 135, zMin, zMax, 18, yMin, yMax);

void draw_box(std::vector<AHCALRecoHit>& RmIsolatedHits,std::string option = "xy") {
    double xysize = AHCALGeometry::x_max*2/18; // cell size approx.
    double zsize = (1191+24)/135; // cell size approx.
    for (const auto& hit : RmIsolatedHits) {
        const double x = hxy_tmp.GetXaxis()->GetBinCenter(hxy_tmp.GetXaxis()->FindBin(hit.Xpos()));
        const double y = hxy_tmp.GetYaxis()->GetBinCenter(hxy_tmp.GetYaxis()->FindBin(hit.Ypos()));
        const double z = hxz_tmp.GetXaxis()->GetBinCenter(hxz_tmp.GetXaxis()->FindBin(hit.Zpos()));
        if (option == "xy") {
            TBox *box_xy = new TBox(x - xysize/2, y - xysize/2, x + xysize/2, y + xysize/2);
            box_xy->SetLineColor(kBlack);
            box_xy->SetLineWidth(1);
            box_xy->SetFillStyle(0); // transparent
            box_xy->Draw("SAME");
        }else if (option == "xz") {
            TBox *box_xz = new TBox(z - zsize/2, x - xysize/2, z + zsize/2, x + xysize/2);
            box_xz->SetLineColor(kBlack);
            box_xz->SetLineWidth(1);
            box_xz->SetFillStyle(0); // transparent
            box_xz->Draw("SAME");
        }else if (option == "yz") {
            TBox *box_yz = new TBox(z - zsize/2, y - xysize/2, z + zsize/2, y + xysize/2);
            box_yz->SetLineColor(kBlack);
            box_yz->SetLineWidth(1);
            box_yz->SetFillStyle(0); // transparent
            box_yz->Draw("SAME");
        }else{
            std::cerr << "Invalid option for draw_box: " << option << "\n";
        }
    }
}
int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <config.yaml> <text file with event list> -s [SamplingRate]\n";
        return 1;
    }
    std::ifstream event_list_file(argv[2]);
    if (!event_list_file.is_open()) {
        std::cerr << "Failed to open event list file: " << argv[2] << "\n";
        return 1;
    }
    int sampling_rate = 1; // default: process all events
    if (argc >= 5 && std::string(argv[3]) == "-s") {
        sampling_rate = std::stoi(argv[4]);
        if (sampling_rate <= 0) {
            std::cerr << "Invalid sampling rate: " << sampling_rate << "\n";
            return 1;
        }
    }
    YAML::Node config = YAML::LoadFile(argv[1]);
    RunContext ctx;
    ctx.config = parse_run_config(config);

    const std::string fname = ctx.config.output;

    FAIR::init_logger("AHCALApp", "eventdisplay/2D.log", spdlog::level::debug);
    RootInput in(fname, "events");
    ReaderRegistry rr;
    rr.register_vector_struct<AHCALRecoHit>("vector<AHCALRecoHit>");
    rr.register_struct<Track>("Track");
    rr.register_struct<SimpleFittedTrack>("SimpleFittedTrack");
    rr.register_struct<AHCALTLURawData>("AHCALTLURawData");
    std::string line;
    int event_count = 0;
    while (std::getline(event_list_file, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        int runNumber;
        int poolIndex;
        long long event_number;
        if (!(iss >> runNumber >> poolIndex >> event_number)) {
            std::cerr << "Invalid line in event list file: " << line << "\n";
            continue;
        }
        if (event_count % sampling_rate != 0) {
            event_count++;
            continue; // skip this event
        }
        event_count++;
        if (!in.read_entry(event_number)) {
            std::cerr << "Event out of range: " << event_number << " (entries=" << in.entries() << ")\n";
            return 1;
        }
        auto reco  = rr.read<std::vector<AHCALRecoHit>>("vector<AHCALRecoHit>", "RecoHits", in);
        auto isolated = rr.read<std::vector<AHCALRecoHit>>("vector<AHCALRecoHit>", "RmIsolatedHits", in);
        auto track = rr.read<Track>("Track", "MuonKFTrack", in);
        auto fitted_track = rr.read<SimpleFittedTrack>("SimpleFittedTrack", "FittedTrack", in);
        auto rawdata = rr.read<AHCALTLURawData>("AHCALTLURawData", "TLURawData", in);
        std::cout << "Event " << event_number << ": " << reco.size()
                << " RecoHits, " << isolated.size() << " IsolatedHits, Track valid=" << track.valid << "\n";
        std::cout << "  Track: chi2=" << track.chi2 << ", ndof=" << track.ndof
                << ", nInTrackHits=" << track.nInTrackHits
                << ", nOutTrackHits=" << track.nOutTrackHits << "\n";
        std::cout << "  FittedTrack: valid=" << fitted_track.valid
                << ", chi2_x/chi2_y=" << fitted_track.chi2_x << "/" << fitted_track.chi2_y << "\n";
        std::cout << "    nInTrackHits=" << fitted_track.inTrackHitsIndices.size()
                << ", nOutTrackHits=" << fitted_track.outTrackHitsIndices.size() << "\n";

        if (reco.empty()) {
            std::cerr << "No RecoHits in event " << event_number << "\n";
            return 1;
        }

        // --- style ---
        gStyle->SetOptStat(0);
        gStyle->SetNumberContours(100);
        gStyle->SetPalette(kViridis);

        TH2D hxy("hxy", Form("RecoHits XY (evt %lld, triggerID=%d);X [mm];Y [mm]", event_number, rawdata.TriggerID),
                18, xMin, xMax, 18, yMin, yMax);
        TH2D hxz("hxz", Form("RecoHits XZ (evt %lld, triggerID=%d);Z [mm];X [mm]", event_number, rawdata.TriggerID),
                135, zMin, zMax, 18, xMin, xMax);
        TH2D hyz("hyz", Form("RecoHits YZ (evt %lld, triggerID=%d);Z [mm];Y [mm]", event_number, rawdata.TriggerID),
                135, zMin, zMax, 18, yMin, yMax);
        for (const auto& hit : reco) {
            const double x = hit.Xpos();
            const double y = hit.Ypos();
            const double z = hit.Zpos();
            // if (hit.Nmip <0.5) continue; // skip low-energy hits
            const double w = hit.Edep; // color = deposited energy [MeV]

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
            lat.DrawLatex(0.12, 0.87, Form("evt=%lld, nhit=%zu", event_number, reco.size()));
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
        draw_box(isolated, "xy");
        drawLabel((TPad*)gPad, "XY projection");

        c.cd(2);
        gPad->SetRightMargin(0.14);
        hxz.Draw("COLZ");
        draw_box(isolated, "xz");
        drawLabel((TPad*)gPad, "XZ projection");

        c.cd(3);
        gPad->SetRightMargin(0.14);
        hyz.Draw("COLZ");
        draw_box(isolated, "yz");
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
                lxy->SetLineStyle(2); // dashed red for MuonKFTrack
                lxy->Draw("SAME");
            }

            // XZ: histogram axes are (Z, x)
            c.cd(2);
            {
                auto* lxz = new TLine(zMin, x1, zMax, x2);
                lxz->SetLineColor(kRed + 1);
                lxz->SetLineWidth(3);
                lxz->SetLineStyle(2); // dashed red for MuonKFTrack
                lxz->Draw("SAME");
            }

            // YZ: histogram axes are (Z, y)
            c.cd(3);
            {
                auto* lyz = new TLine(zMin, y1, zMax, y2);
                lyz->SetLineColor(kRed + 1);
                lyz->SetLineWidth(3);
                lyz->SetLineStyle(2); // dashed red for MuonKFTrack
                lyz->Draw("SAME");
            }
        }
        if (fitted_track.valid) {
            // ---- overlay SimpleFittedTrack (blue line) ----
            const double x1 = fitted_track.init_pos_x + fitted_track.direction_x * (zMin);
            const double y1 = fitted_track.init_pos_y + fitted_track.direction_y * (zMin);
            const double x2 = fitted_track.init_pos_x + fitted_track.direction_x * (zMax);
            const double y2 = fitted_track.init_pos_y + fitted_track.direction_y * (zMax);

            // XY: (x(z),y(z))
            c.cd(1);
            {
                auto* lxy = new TLine(x1, y1, x2, y2);
                lxy->SetLineColor(kBlue + 1);
                lxy->SetLineWidth(3);
                lxy->SetLineStyle(2); // dashed blue for SimpleFittedTrack
                lxy->Draw("SAME");
            }

            // XZ: histogram axes are (Z, x)
            c.cd(2);
            {
                auto* lxz = new TLine(zMin, x1, zMax, x2);
                lxz->SetLineColor(kBlue + 1);
                lxz->SetLineWidth(3);
                lxz->SetLineStyle(2); // dashed blue for SimpleFittedTrack
                lxz->Draw("SAME");
            }

            // YZ: histogram axes are (Z, y)
            c.cd(3);
            {
                auto* lyz = new TLine(zMin, y1, zMax, y2);
                lyz->SetLineColor(kBlue + 1);
                lyz->SetLineWidth(3);
                lyz->SetLineStyle(2); // dashed blue for SimpleFittedTrack
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
        info.DrawLatex(0.10, 0.74, "MuonKFTrack: red line");
        info.DrawLatex(0.10, 0.69, "SimpleFittedTrack: blue line");
        info.DrawLatex(0.10, 0.64, "Black boxes: RmIsolatedHits");

        c.Update();
        TLatex RunNum;
        RunNum.SetNDC();
        RunNum.SetTextSize(0.04);
        RunNum.DrawLatex(0.12, 0.96, Form("Run %d", runNumber));
        // const std::string outPdf = Form("eventdisplay/evt%lld.pdf", ievt);
        const std::string outPng = Form("eventdisplay_eff/Run%d_evt%lld.png", runNumber, event_number);
        // c.SaveAs(outPdf.c_str());
        c.SaveAs(outPng.c_str());
    }
    

    // std::cout << "Saved: " << outPdf << " and " << outPng << "\n";
    return 0;
}