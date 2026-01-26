#include "TrackFitAlg.hpp"
#include <TGraphErrors.h>
#include <TF1.h>
#include "common/config/YAMLUtil.hpp"
#include "common/edm/EDM.hpp"
namespace AHCALRecoAlg {
    void TrackFitAlg::execute(EventStore& evt) {
        std::string m_in_recohit_key = m_cfg.in_recohit_key;
        std::string m_out_track_key = m_cfg.out_track_key;
        double threshold_xy = m_cfg.threshold_xy;
        auto recohits = evt.get<std::vector<AHCALRecoHit>>(m_in_recohit_key);
        if (recohits.empty()) {
            // throw std::runtime_error("TrackFitAlg: No input reco hits found");
            // LOG_WARN("TrackFitAlg: No input reco hits found.");
            LOG_DEBUG("TrackFitAlg: No input reco hits found.");
            SimpleFittedTrack empty_track;
            empty_track.valid = false;
            evt.put(m_out_track_key, std::move(empty_track));
            return;
        }
        std::unique_ptr<TGraphErrors> graph_xz = std::make_unique<TGraphErrors>();
        std::unique_ptr<TGraphErrors> graph_yz = std::make_unique<TGraphErrors>();
        SimpleFittedTrack track;
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
        if (track.nTotalHits < 3){
            LOG_DEBUG("TrackFitAlg: Not enough hits to fit a track. Hits found: {}", track.nTotalHits);
            evt.put(m_out_track_key, std::move(track));
            return; // Not enough hits to fit
        }
        track.valid = true;
        // Fit lines
        auto pol1_x = std::make_unique<TF1>("pol1_x", "[0] + [1]*x");
        auto pol1_y = std::make_unique<TF1>("pol1_y", "[0] + [1]*x");
        pol1_x->SetParLimits(1, -20, 20); // Limit slope to reasonable values
        pol1_y->SetParLimits(1, -20, 20);
        graph_xz->Fit("pol1_x", "Q"); // "Q" for quiet mode
        graph_yz->Fit("pol1_y", "Q");
        if (graph_xz->GetFunction("pol1_x") == nullptr || graph_yz->GetFunction("pol1_y") == nullptr) {
            LOG_WARN("TrackFitAlg: Fit failed.");
            track.valid = false;
            evt.put(m_out_track_key, std::move(track));
            return;
        }
            
        // Extract fit parameters
        double p0_x = graph_xz->GetFunction("pol1_x")->GetParameter(0);
        double p1_x = graph_xz->GetFunction("pol1_x")->GetParameter(1);
        double p0_y = graph_yz->GetFunction("pol1_y")->GetParameter(0);
        double p1_y = graph_yz->GetFunction("pol1_y")->GetParameter(1);
        // Create Track object
        track.init_pos_x = p0_x;
        track.init_pos_y = p0_y;
        track.direction_x = p1_x;
        track.direction_y = p1_y;
        track.chi2_x = graph_xz->GetFunction("pol1_x")->GetChisquare();
        track.chi2_y = graph_yz->GetFunction("pol1_y")->GetChisquare();
        track.ndf = graph_xz->GetFunction("pol1_x")->GetNDF();
        // Classify hits
        int index = 0;
        for (const auto& hit : recohits) {
            double x_pred = p0_x + p1_x * hit.Zpos();
            double y_pred = p0_y + p1_y * hit.Zpos();
            double dx = hit.Xpos() - x_pred;
            double dy = hit.Ypos() - y_pred;
            if (abs(dx) < threshold_xy && abs(dy) < threshold_xy) {
                track.inTrackHits.push_back(hit);
                track.inTrackHitsIndices.push_back(index);
            } else {
                track.outTrackHits.push_back(hit);
                track.outTrackHitsIndices.push_back(index);
            }
            index++;
        }
        // Store track
        evt.put(m_out_track_key, std::move(track));
    }
    void TrackFitAlg::parse_cfg(const YAML::Node& n) {
        m_cfg.in_recohit_key = get_or<std::string>(n, "in_recohit_key", m_cfg.in_recohit_key);
        m_cfg.out_track_key = get_or<std::string>(n, "out_track_key", m_cfg.out_track_key);
        m_cfg.threshold_xy = get_or<double>(n, "threshold_xy", m_cfg.threshold_xy);
    }
}

        
    