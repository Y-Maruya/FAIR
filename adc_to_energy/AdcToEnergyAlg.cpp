#include "AdcToEnergyAlg.hpp"

#include "calibration/RefValues.hpp"
#include "common/AlgRegistry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/edm/EDM.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

AHCAL_REGISTER_ALG(AHCALRecoAlg::AdcToEnergyAlg, "AdcToEnergyAlg")

namespace AHCALRecoAlg {
namespace {

bool isFinitePositive(double value) {
  return std::isfinite(value) && value > 0.0;
}

} // namespace

AdcToEnergyAlg::~AdcToEnergyAlg() = default;

void AdcToEnergyAlg::initialize() {
  loadCalibrations();
}

void AdcToEnergyAlg::init_by_run() {
  loadCalibrations();
}

void AdcToEnergyAlg::loadCalibrations() {
  const int runNumber = ctx().config.runNumber;
  if (runNumber == loaded_run_) {
    return;
  }
  if (runNumber <= 0) {
    throw std::runtime_error("AdcToEnergyAlg requires a positive runNumber to read CalibDBIO");
  }

  initialize_pedestal(runNumber);
  initialize_mip(runNumber);
  initialize_intercalib(runNumber);
  loaded_run_ = runNumber;
  LOG_INFO("AdcToEnergyAlg: loaded CalibDBIO constants for run {}", runNumber);
}

void AdcToEnergyAlg::initialize_pedestal(int runNumber) {
  CalibDBIO::PedestalReader reader(runNumber);
  ped_map_ = reader.getPedestalMapPtr();
  LOG_INFO("AdcToEnergyAlg: loaded {} pedestal entries from CalibDBIO", ped_map_->size());
}

void AdcToEnergyAlg::initialize_mip(int runNumber) {
  CalibDBIO::MIPReader reader(runNumber);
  mip_map_ = reader.getMIPMapPtr();
  LOG_INFO("AdcToEnergyAlg: loaded {} MIP entries from CalibDBIO", mip_map_->size());
}

void AdcToEnergyAlg::initialize_intercalib(int runNumber) {
  CalibDBIO::InterCalibReader reader(runNumber);
  hglg_map_ = reader.getHGLGRatioMapPtr();
  hg_saturation_map_ = reader.getHGSaturationMapPtr();
  LOG_INFO("AdcToEnergyAlg: loaded {} HGLG entries and {} HG saturation entries from CalibDBIO",
           hglg_map_->size(), hg_saturation_map_->size());
}

void AdcToEnergyAlg::execute(EventStore &evt) { 
  const auto& raw_hits = evt.get<std::vector<AHCALRawHit>>(m_in_rawhit_key);

  std::vector<AHCALRecoHit> reco_hits;
  reco_hits.reserve(raw_hits.size());

  const double SwitchPoint = AHCALRefValues::SwitchPoint;
  std::size_t n_missing_ped = 0;
  std::size_t n_missing_mip = 0;
  std::size_t n_missing_hglg = 0;
  std::size_t n_missing_saturation = 0;
  std::size_t n_invalid_calib = 0;

  for (const auto &raw_hit : raw_hits) {
    if (raw_hit.hittag == 0) continue; // skip hittag==0
    double hg_ped = 0.0;
    double lg_ped = 0.0;
    double mpv = 0.0;
    double gain_ratio = 0.0;
    double gain_plat = 0.0;
    double hg_ped_error = 0.0;
    double lg_ped_error = 0.0;
    double mpv_error = 0.0;
    double gain_ratio_error = 0.0;
    const auto ped_it = ped_map_->find(raw_hit.cellID);
    if (ped_it == ped_map_->end()) {
      ++n_missing_ped;
      if (m_cfg.mask_level > 0){
        continue;
      }else{
        hg_ped = AHCALRefValues::ref_ped_highgain;
        lg_ped = AHCALRefValues::ref_ped_lowgain;
        hg_ped_error = AHCALRefValues::ref_ped_highgain_error; // default error value for HG pedestal
        lg_ped_error = AHCALRefValues::ref_ped_lowgain_error; // default error value for LG pedestal
      }
    }else{
      hg_ped = ped_it->second.HighGainPeak;
      lg_ped = ped_it->second.LowGainPeak;
      hg_ped_error = ped_it->second.HighGainPeakError;
      lg_ped_error = ped_it->second.LowGainPeakError;
      if (ped_it->second.HighGainUsable == 0 || (ped_it->second.HighGainStatus != 0 && ped_it->second.HighGainStatus != 999) || ped_it->second.HighGainPeak <= 0.0) {
        ++n_missing_ped;
        if (m_cfg.mask_level > 0){
          continue;
        }else{
          hg_ped = AHCALRefValues::ref_ped_highgain;
          hg_ped_error = AHCALRefValues::ref_ped_highgain_error;
        }
      } 
      if (ped_it->second.LowGainUsable == 0 || (ped_it->second.LowGainStatus != 0 && ped_it->second.LowGainStatus != 999) || ped_it->second.LowGainPeak <= 0.0) {
        ++n_missing_ped;
        if (m_cfg.mask_level > 0){
          continue;
        }else{
          lg_ped = AHCALRefValues::ref_ped_lowgain;
          lg_ped_error = AHCALRefValues::ref_ped_lowgain_error;
        }
      }
    }
    const auto mip_it = mip_map_->find(raw_hit.cellID);
    if (mip_it == mip_map_->end()) {
      ++n_missing_mip;
      continue;
    }
    if (mip_it->second.mpv <= 0.0) {
      ++n_missing_mip;
      continue;
    }
    mpv = mip_it->second.mpv;
    mpv_error = mip_it->second.mpverror;

    const auto hglg_it = hglg_map_->find(raw_hit.cellID);
    if (hglg_it == hglg_map_->end()) {
      ++n_missing_hglg;
      if (m_cfg.mask_level > 1){
        continue;
      }else{
        gain_ratio = AHCALRefValues::ref_gain_ratio;
        gain_ratio_error = AHCALRefValues::ref_gain_ratio_error; // default error value for HG/LG ratio
      }
    }else{
      gain_ratio = hglg_it->second.slope;
      gain_ratio_error = hglg_it->second.slopeerror;
      if (gain_ratio <= 0.0 || hglg_it->second.qualityflag == 1) {
        ++n_missing_hglg;
        if (m_cfg.mask_level > 1){
          continue;
        }else{
          gain_ratio = AHCALRefValues::ref_gain_ratio;
          gain_ratio_error = AHCALRefValues::ref_gain_ratio_error; // default error value for HG/LG ratio
        }
      }
    }

    const auto saturation_it = hg_saturation_map_->find(raw_hit.cellID);
    if (saturation_it == hg_saturation_map_->end()) {
      ++n_missing_saturation;
      continue;
    }
    if (saturation_it->second <= 0.0) {
      ++n_missing_saturation;
      continue;
    }
    gain_plat = saturation_it->second;

    AHCALRecoHit reco_hit;
    reco_hit.cellID = raw_hit.cellID;
    reco_hit.index = raw_hit.index;
    if (!isFinitePositive(mpv) ||
        !std::isfinite(hg_ped) ||
        !std::isfinite(lg_ped) ||
        !isFinitePositive(gain_ratio) ||
        !isFinitePositive(gain_plat) || 
        !isFinitePositive(hg_ped_error) ||
        !isFinitePositive(lg_ped_error) ||
        !isFinitePositive(mpv_error) ||
        !isFinitePositive(gain_ratio_error)) {
      ++n_invalid_calib;
      continue;
    }

    const double hg = static_cast<double>(raw_hit.hg_adc);
    const double lg = static_cast<double>(raw_hit.lg_adc);

    if ((hg - hg_ped) < (gain_plat - SwitchPoint)) {
      reco_hit.Nmip = (hg - hg_ped) / mpv;
      reco_hit.NmipError = std::sqrt(std::pow(hg_ped_error / mpv, 2) + std::pow((hg - hg_ped) * mpv_error / (mpv * mpv), 2));
      reco_hit.Edep = (hg - hg_ped) * AHCALRefValues::MIP_E / mpv;
      reco_hit.EdepError = std::sqrt(std::pow(hg_ped_error * AHCALRefValues::MIP_E / mpv, 2) + std::pow((hg - hg_ped) * AHCALRefValues::MIP_E * mpv_error / (mpv * mpv), 2));
    } else {
      reco_hit.Nmip = (lg - lg_ped) * gain_ratio / mpv;
      reco_hit.NmipError = std::sqrt(std::pow(lg_ped_error * gain_ratio / mpv, 2) + std::pow((lg - lg_ped) * gain_ratio_error / mpv, 2) + std::pow((lg - lg_ped) * gain_ratio * mpv_error / (mpv * mpv), 2));
      reco_hit.Edep = (lg - lg_ped) * gain_ratio * AHCALRefValues::MIP_E / mpv;
      reco_hit.EdepError = std::sqrt(std::pow(lg_ped_error * AHCALRefValues::MIP_E / mpv, 2) + std::pow((lg - lg_ped) * AHCALRefValues::MIP_E * mpv_error / (mpv * mpv), 2));
    }

    if (reco_hit.Edep < 0) {
      reco_hit.Edep = 0;
      // reco_hit.EdepError = 0;
      reco_hit.Nmip = 0;
      // reco_hit.NmipError = 0;
    }
    if (reco_hit.Nmip > 1e6) {
      LOG_DEBUG("Large Nmip detected: cellID={} Nmip={}", reco_hit.cellID, reco_hit.Nmip);
      LOG_DEBUG("  raw HG={} LG={} HG_ped={} LG_ped={} gain_ratio={} gain_plat={} mpv={}",
                hg, lg, hg_ped, lg_ped, gain_ratio, gain_plat, mpv);
    }
    reco_hits.push_back(reco_hit);
  }

  LOG_DEBUG("Converted {} raw hits to reco hits.", reco_hits.size());
  if (n_missing_ped || n_missing_mip || n_missing_hglg || n_missing_saturation || n_invalid_calib) {
    LOG_DEBUG("AdcToEnergyAlg: skipped hits due to missing/invalid calibration: pedestal={} mip={} hglg={} hg_saturation={} invalid={}",
              n_missing_ped, n_missing_mip, n_missing_hglg, n_missing_saturation, n_invalid_calib);
  }
  evt.put(m_out_recohit_key, std::move(reco_hits));
}

void AdcToEnergyAlg::parse_cfg(const YAML::Node& n) {
  m_cfg.in_rawhit_key = get_or<std::string>(n, "in_rawhit_key", m_cfg.in_rawhit_key);
  m_cfg.out_recohit_key = get_or<std::string>(n, "out_recohit_key", m_cfg.out_recohit_key);
  m_cfg.mask_level = get_or<int>(n, "mask_level", m_cfg.mask_level);
  if (m_cfg.mask_level < 0 || m_cfg.mask_level > 2) {
    throw std::runtime_error("Invalid mask_level in AdcToEnergyAlg configuration: "+ std::to_string(m_cfg.mask_level) + ". Valid values are 0, 1, or 2.");
  }
  m_in_rawhit_key = m_cfg.in_rawhit_key;
  m_out_recohit_key = m_cfg.out_recohit_key;
}

} // namespace AHCALRecoAlg
