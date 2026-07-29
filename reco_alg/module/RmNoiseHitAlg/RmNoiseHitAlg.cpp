#include "RmNoiseHitAlg.hpp"

#include "common/edm/EDM.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"
#include <muParser.h>
AHCAL_REGISTER_ALG(AHCALRecoAlg::RmNoiseHitAlg, "RmNoiseHitAlg")
namespace AHCALRecoAlg {

void RmNoiseHitAlg::computeNoiseFlags(const std::vector<AHCALRawHit>& hits, std::vector<int>& isoFlag, bool hittag1_only) {
    for (size_t i = 0; i < hits.size(); ++i) {
        if (hittag1_only && hits[i].hittag != 1) {
            isoFlag[i] = 2; // Not Noise, but will be filtered out later
            continue;
        }
        const auto& h = hits[i];
        double hg_ped = 0.0;
        double lg_ped = 0.0;
        auto itp = ped_map_->find(h.cellID);
        if (itp != ped_map_->end()) {
            hg_ped = itp->second.HighGainPeak;
            lg_ped = itp->second.LowGainPeak;
        } else {
            isoFlag[i] = 2; // Not Noise, but will be filtered out later due to missing pedestal
            continue;
        }
        if (itp->second.HighGainStatus != 0 || itp->second.LowGainStatus != 0) {
            isoFlag[i] = 2; // Not Noise, but will be filtered out later due to bad pedestal
            continue;
        }
        const double hg = h.hg_adc - hg_ped;
        const double lg = h.lg_adc - lg_ped;
        if (cut && cut->eval(hg, lg)) {
            isoFlag[i] = 1; // Noise
        } else {
            isoFlag[i] = 0; // Not Noise
        }
    }
}

void RmNoiseHitAlg::execute(EventStore& evt) {
    auto& rawhits = evt.get<std::vector<AHCALRawHit>>(m_cfg.in_rawhit_key);
    RmNoiseSummary summary;
    std::vector<AHCALRawHit> out_rawhits;
    std::vector<AHCALRawHit> out_removed_hits; // for output if m_cfg.output_Noise_hits is true
    std::vector<int> isoFlag(rawhits.size(), 0);
    computeNoiseFlags(rawhits, isoFlag, m_cfg.hittag1_only);
    int total_hits = 0;
    int noise_hits = 0;
    for (size_t i = 0; i < rawhits.size(); ++i) {
        if (isoFlag[i] == 0) { // Not Noise
            out_rawhits.push_back(rawhits[i]);
            total_hits++;
        } else if (isoFlag[i] == 1) { // Noise
            noise_hits++;
            total_hits++;
            if (m_cfg.output_removed_hits) {
                out_removed_hits.push_back(rawhits[i]);
            }
        } else if (isoFlag[i] == 2) { // Hittag != 1, will be filtered out but not counted as noise
            // Do nothing, will be filtered out by hittag later
        }
    }
    summary.total_hits = total_hits;
    summary.noise_hits = noise_hits;
    summary.noise_fraction = (summary.total_hits > 0) ? static_cast<double>(noise_hits) / summary.total_hits : 0.0;
    LOG_DEBUG("RmNoiseHitAlg: total_hits={}, noise_hits={}, noise_fraction={}",
         summary.total_hits,
         summary.noise_hits,
         summary.noise_fraction);
    evt.put(m_cfg.out_rawhit_key, std::move(out_rawhits));
    evt.put(m_cfg.out_summary_key, std::move(summary));
    if (m_cfg.output_removed_hits) {
        evt.put(m_cfg.out_removedrawhit_key, std::move(out_removed_hits));
    }
}

void RmNoiseHitAlg::parse_cfg(const YAML::Node& n) {
    m_cfg.in_rawhit_key = get_or<std::string>(n, "in_rawhit_key", m_cfg.in_rawhit_key);
    m_cfg.out_rawhit_key = get_or<std::string>(n, "out_rawhit_key", m_cfg.out_rawhit_key);
    m_cfg.criteria = get_or<std::string>(n, "criteria", m_cfg.criteria);
    m_cfg.output_removed_hits = get_or<bool>(n, "output_removed_hits", m_cfg.output_removed_hits);
    m_cfg.out_removedrawhit_key = get_or<std::string>(n, "out_removedrawhit_key", m_cfg.out_removedrawhit_key);
    m_cfg.out_summary_key = get_or<std::string>(n, "out_summary_key", m_cfg.out_summary_key);
    m_cfg.hittag1_only = get_or<bool>(n, "hittag1_only", m_cfg.hittag1_only);
}
}
