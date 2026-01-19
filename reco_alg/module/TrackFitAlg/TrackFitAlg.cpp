#include "TrackFitAlg.hpp"
#include "TGraphErrors.h"
#include "TF1.h"
namespace AHCALRecoAlg {
    void TrackFitAlg::execute(EventStore& evt) {
        auto recohits = evt.get<std::vector<AHCALRecoAndRawHit>>(m_in_recohit_key);
        if (recohits.empty()) {
            throw std::runtime_error("TrackFitAlg: No input reco hits found");
        }
        std::unique_ptr<TGraphErrors> graph_xz = std::make_unique<TGraphErrors>();
        std::unique_ptr<TGraphErrors> graph_yz = std::make_unique<TGraphErrors>();
        Track track;
        for (const auto& hit : recohits) {
            if (hit.Nmip <0.5) continue; // MIP cut
            track.nTotalHits++;
            double x = hit.Xpos();
            double y = hit.Ypos();
            double z = hit.Zpos();
            int index = graph_xz->GetN();
            graph_xz->SetPoint(index, z, x);
            graph_yz->SetPoint(index, z, y);
            graph_xz->SetPointError(index, AHCALGeometry::z_size/2, AHCALGeometry::xy_size/2);
            graph_yz->SetPointError(index, AHCALGeometry::z_size/2, AHCALGeometry::xy_size/2);
        }
        // Fit lines
        graph_xz->Fit("pol1", "Q"); // "Q" for quiet mode
        graph_yz->Fit("pol1", "Q");
        // Extract fit parameters
        double p0_x = graph_xz->GetFunction("pol1")->GetParameter(0);
        double p1_x = graph_xz->GetFunction("pol1")->GetParameter(1);
        double p0_y = graph_yz->GetFunction("pol1")->GetParameter(0);
        double p1_y = graph_yz->GetFunction("pol1")->GetParameter(1);
        // Create Track object
        track.init_pos = {p0_x, p0_y};
        track.direction = {p1_x, p1_y};
        track.chi2 = {graph_xz->GetFunction("pol1")->GetChisquare(), graph_yz->GetFunction("pol1")->GetChisquare()};
        track.ndf = graph_xz->GetFunction("pol1")->GetNDF();
        // Classify hits
        for (const auto& hit : recohits) {
            double x_pred = p0_x + p1_x * hit.Zpos();
            double y_pred = p0_y + p1_y * hit.Zpos();
            double dx = hit.Xpos() - x_pred;
            double dy = hit.Ypos() - y_pred;
            if (abs(dx) < threshold_xy && abs(dy) < threshold_xy) {
                track.inTrackHits.push_back(hit);
            } else {
                track.outTrackHits.push_back(hit);
            }
        }
        // Store track
        evt.put(m_out_track_key, std::move(track));
    }
}

        
    