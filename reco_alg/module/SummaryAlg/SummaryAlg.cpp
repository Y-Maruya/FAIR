#include "SummaryAlg.hpp"
#include <TGraphErrors.h>
#include <TF1.h>
#include "common/config/YAMLUtil.hpp"
#include "common/edm/EDM.hpp"
#include "common/AlgRegistry.hpp"
#include <set>
AHCAL_REGISTER_ALG(AHCALRecoAlg::SummaryAlg, "SummaryAlg")
namespace AHCALRecoAlg {
    void SummaryAlg::execute(EventStore& evt) {
        std::string m_in_recohit_key = m_cfg.in_recohit_key;
        std::string m_out_summary_key = m_cfg.out_summary_key;
        auto recohits = evt.get<std::vector<AHCALRecoHit>>(m_in_recohit_key);
        EventSummary summary;
        summary.event_counter = evt.event_counter();
        if (recohits.empty()) {
            // throw std::runtime_error("SummaryAlg: No input reco hits found");
            // LOG_WARN("SummaryAlg: No input reco hits found.");
            evt.put(m_out_summary_key, std::move(summary));
            return;
        }
        summary.nHits = recohits.size();
        std::set<int> hitLayers;
        for (const auto& hit : recohits) {
            if (hit.Nmip > 0.5) {
                summary.nHitAbove0p5MIP++;
            }
            hitLayers.insert(hit.layer());
        }
        summary.nHitLayers = hitLayers.size();
        std::set<int> hitLayersAbove0p5MIP;
        for (const auto& hit : recohits) {
            if (hit.Nmip > 0.5) {
                hitLayersAbove0p5MIP.insert(hit.layer());
            }
        }
        summary.nHitLayersAbove0p5MIP = hitLayersAbove0p5MIP.size();
        for (const auto& hit : recohits) {
            summary.TotalEnergy += hit.Edep;
            summary.TotalnMIP += hit.Nmip;
            if (hit.Edep > summary.MaxEnergy) {
                summary.MaxEnergy = hit.Edep;
            }
            if (hit.Nmip > summary.MaxnMIP) {
                summary.MaxnMIP = hit.Nmip;
            }
            summary.CenterOfGravityX += hit.Xpos() * hit.Edep;
            summary.CenterOfGravityY += hit.Ypos() * hit.Edep;
            summary.CenterOfGravityZ += hit.Zpos() * hit.Edep;
        }
        if (summary.TotalEnergy > 0.0) {
            summary.CenterOfGravityX /= summary.TotalEnergy;
            summary.CenterOfGravityY /= summary.TotalEnergy;
            summary.CenterOfGravityZ /= summary.TotalEnergy;
            for (const auto& hit : recohits) {
                summary.RMSX += hit.Edep * std::pow(hit.Xpos() - summary.CenterOfGravityX, 2);
                summary.RMSY += hit.Edep * std::pow(hit.Ypos() - summary.CenterOfGravityY, 2);
                summary.RMSZ += hit.Edep * std::pow(hit.Zpos() - summary.CenterOfGravityZ, 2);
            }
            summary.RMSX = std::sqrt(summary.RMSX / summary.TotalEnergy);
            summary.RMSY = std::sqrt(summary.RMSY / summary.TotalEnergy);
            summary.RMSZ = std::sqrt(summary.RMSZ / summary.TotalEnergy);
        }
        if (summary.nHits > m_max_nhit_event.second) {
            m_max_nhit_event = {evt.event_counter(), summary.nHits};
        }
        if (summary.TotalEnergy > m_max_energy_event.second) {
            m_max_energy_event = {evt.event_counter(), summary.TotalEnergy};
        }
        evt.put(m_out_summary_key, std::move(summary));
    }
    void SummaryAlg::parse_cfg(const YAML::Node& n) {
        m_cfg.in_recohit_key = get_or<std::string>(n, "in_recohit_key", m_cfg.in_recohit_key);
        m_cfg.out_summary_key = get_or<std::string>(n, "out_summary_key", m_cfg.out_summary_key);
        m_cfg.output_summary_txt = get_or<bool>(n, "output_summary_txt", m_cfg.output_summary_txt);
        m_cfg.summary_txt_path = get_or<std::string>(n, "summary_txt_path", m_cfg.summary_txt_path);
    }
    void SummaryAlg::initialize() {
        if (m_cfg.output_summary_txt) {
            m_summary_txt_file.open(m_cfg.summary_txt_path);
            if (!m_summary_txt_file.is_open()) {
                throw std::runtime_error("SummaryAlg: Failed to open summary text file: " + m_cfg.summary_txt_path);
            }
            // Write header line
            // m_summary_txt_file << "nHits\tnHitAbove0p5MIP\tnHitLayers\tnHitLayersAbove0p5MIP\tTotalEnergy\tTotalnMIP\tMaxEnergy\tMaxnMIP\tCenterOfGravityX\tCenterOfGravityY\tCenterOfGravityZ\tRMSX\tRMSY\tRMSZ\n";
        }
    }
    void SummaryAlg::finalize() {
        if (m_cfg.output_summary_txt && m_summary_txt_file.is_open()) {
            m_summary_txt_file << "Event with max nHits: " << m_max_nhit_event.first << " (nHits = " << m_max_nhit_event.second << ")\n";
            m_summary_txt_file << "Event with max TotalEnergy: " << m_max_energy_event.first << " (TotalEnergy = " << m_max_energy_event.second << ")\n";
            m_summary_txt_file.close();
        }
        LOG_INFO("SummaryAlg: Event with max nHits: {}, nHits = {}", m_max_nhit_event.first, m_max_nhit_event.second);
        LOG_INFO("SummaryAlg: Event with max TotalEnergy: {}, TotalEnergy = {}", m_max_energy_event.first, m_max_energy_event.second);
    }
}

        
    