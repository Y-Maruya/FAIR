#include "AdcToEnergyReadTTreeAlg.hpp"

#include "calibration/RefValues.hpp"
#include "common/edm/EDM.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>

#include <utility>
#include <vector>
AHCAL_REGISTER_ALG(AHCALRecoAlg::AdcToEnergyReadTTreeAlg, "AdcToEnergyReadTTreeAlg")
namespace AHCALRecoAlg {

bool AdcToEnergyReadTTreeAlg::initialize_mip(){
  std::string mip_file_name = m_cfg.mip_file;
  std::string cut_string = m_cfg.mip_cut_string;
  file_cellid_version = m_cfg.mip_cellid_version;

  m_in_file = std::make_unique<TFile>(mip_file_name.c_str(), "READ");
  if (!m_in_file || m_in_file->IsZombie()) {
    LOG_ERROR("Error: cannot open file: {}", mip_file_name);
    return false;
  }

  m_in_tree = dynamic_cast<TTree *>(m_in_file->Get("mip"));
  if (!m_in_tree) {
    LOG_ERROR("Error: cannot find TTree 'mip' in file: {}", mip_file_name);
    return false;
  }

  ROOT::RDataFrame df(*m_in_tree);

  ROOT::RDF::RNode df_cut = df;
  if (!cut_string.empty()) df_cut = df.Filter(cut_string);

  if (!cut_string.empty()) {
    if (df_cut.Count().GetValue() == 0) {
      LOG_ERROR("No entries left after applying cut: {}", cut_string);
      return false;
    }
    LOG_INFO("Applied cut string: {}", cut_string);
  }

  const auto cellids = df_cut.Take<int>("cellid").GetValue();
  const auto mpvs    = df_cut.Take<double>("MPV").GetValue();
  if (cellids.size() != mpvs.size()) {
    LOG_ERROR("MIP tree size mismatch: cellid={} MPV={}", cellids.size(),
              mpvs.size());
    return false;
  }

  for (size_t i = 0; i < cellids.size(); ++i) {
    const int conv_cellid = cellid_conversion(cellids[i]);
    const std::size_t idx = channel_index_from_cellid(conv_cellid);
    if (idx >= total_channel_count()) continue;
    mip_values_[idx] = mpvs[i];
  }

  int loaded = 0;
  LOG_INFO("Loaded {} MIP entries from {}", cellids.size(), mip_file_name);

  for (double &mip : mip_values_) {
    if (mip <= 100.0) {
      mip = AHCALRefValues::ref_MIP;
    } else {
      ++loaded;
    }
  }

  LOG_INFO("Reference MIP values assigned for cut channels.");
  LOG_INFO("Total channels with reference MIP: {}",
           static_cast<int>(mip_values_.size()) - loaded);
  
  // Clean up ROOT objects to avoid memory leaks
  if (m_in_file) {
    m_in_file->Close();
    m_in_file.reset();
  }
  m_in_tree = nullptr;
  
  return true;
}

bool AdcToEnergyReadTTreeAlg::initialize_ped( ){
  std::string ped_file_name = m_cfg.ped_file;
  std::string cut_string = m_cfg.ped_cut_string;
  file_cellid_version = m_cfg.ped_cellid_version;

  m_in_file = std::make_unique<TFile>(ped_file_name.c_str(), "READ");
  if (!m_in_file || m_in_file->IsZombie()) {
    LOG_ERROR("Error: cannot open file: {}", ped_file_name);
    return false;
  }

  m_in_tree = dynamic_cast<TTree *>(m_in_file->Get("pedestal"));
  if (!m_in_tree) {
    LOG_ERROR("Error: cannot find TTree 'pedestal' in file: {}", ped_file_name);
    return false;
  }

  ROOT::RDataFrame df(*m_in_tree);

  ROOT::RDF::RNode df_cut = df;
  if (!cut_string.empty()) df_cut = df.Filter(cut_string);

  if (!cut_string.empty()) {
    if (df_cut.Count().GetValue() == 0) {
      LOG_ERROR("No entries left after applying cut: {}", cut_string);
      return false;
    }
    LOG_INFO("Applied cut string: {}", cut_string);
  }

  const auto cellids = df_cut.Take<int>("cellid").GetValue();
  const auto hg      = df_cut.Take<double>("highgain_peak").GetValue();
  const auto lg      = df_cut.Take<double>("lowgain_peak").GetValue();

  if (cellids.size() != hg.size() || cellids.size() != lg.size()) {
    LOG_ERROR("Pedestal tree size mismatch: cellid={} HG={} LG={}",
              cellids.size(), hg.size(), lg.size());
    return false;
  }

  for (size_t i = 0; i < cellids.size(); ++i) {
    const int conv_cellid = cellid_conversion(cellids[i]);
    const std::size_t idx = channel_index_from_cellid(conv_cellid);
    if (idx >= total_channel_count()) continue;
    hg_ped_values_[idx] = hg[i];
    lg_ped_values_[idx] = lg[i];
  }

  LOG_INFO("Loaded {} pedestal entries from {}", cellids.size(), ped_file_name);
  LOG_INFO("Reference pedestal values assigned for cut channels.");

  // Clean up ROOT objects to avoid memory leaks
  if (m_in_file) {
    m_in_file->Close();
    m_in_file.reset();
  }
  m_in_tree = nullptr;

  return true;
}

bool AdcToEnergyReadTTreeAlg::initialize_dac( ){
  std::string dac_file_name = m_cfg.dac_file;
  std::string cut_string = m_cfg.dac_cut_string;
  file_cellid_version = m_cfg.dac_cellid_version;

  m_in_file = std::make_unique<TFile>(dac_file_name.c_str(), "READ");
  if (!m_in_file || m_in_file->IsZombie()) {
    LOG_ERROR("Error: cannot open file: {}", dac_file_name);
    return false;
  }

  m_in_tree = dynamic_cast<TTree *>(m_in_file->Get("dac"));
  if (!m_in_tree) {
    LOG_ERROR("Error: cannot find TTree 'dac' in file: {}", dac_file_name);
    return false;
  }

  ROOT::RDataFrame df(*m_in_tree);

  ROOT::RDF::RNode df_cut = df;
  if (!cut_string.empty()) df_cut = df.Filter(cut_string);

  if (!cut_string.empty()) {
    if (df_cut.Count().GetValue() == 0) {
      LOG_ERROR("No entries left after applying cut: {}", cut_string);
      return false;
    }
    LOG_INFO("Applied cut string: {}", cut_string);
  }

  const auto cellids = df_cut.Take<int>("cellid").GetValue();
  const auto slopes  = df_cut.Take<float>("slope").GetValue();
  const auto plats   = df_cut.Take<float>("plat").GetValue();

  if (cellids.size() != slopes.size() || cellids.size() != plats.size()) {
    LOG_ERROR("DAC tree size mismatch: cellid={} slope={} plat={}", cellids.size(),
              slopes.size(), plats.size());
    return false;
  }

  for (size_t i = 0; i < cellids.size(); ++i) {
    const int conv_cellid = cellid_conversion(cellids[i]);
    const std::size_t idx = channel_index_from_cellid(conv_cellid);
    if (idx >= total_channel_count()) continue;
    gainratio_values_[idx] = slopes[i];
    gainplat_values_[idx] = static_cast<int>(plats[i]);
  }

  LOG_INFO("Loaded {} DAC entries from {}", cellids.size(), dac_file_name);
  LOG_INFO("Reference DAC values assigned for cut channels.");

  // Clean up ROOT objects to avoid memory leaks
  if (m_in_file) {
    m_in_file->Close();
    m_in_file.reset();
  }
  m_in_tree = nullptr;

  return true;
}

AdcToEnergyReadTTreeAlg::~AdcToEnergyReadTTreeAlg() {
  if (m_in_file) {
    m_in_file->Close();
  }
  m_in_tree = nullptr;
}


std::size_t AdcToEnergyReadTTreeAlg::total_channel_count() {
  return static_cast<std::size_t>(AHCALGeometry::Layer_No) *
         static_cast<std::size_t>(AHCALGeometry::chip_No) *
         static_cast<std::size_t>(AHCALGeometry::channel_No);
}

std::size_t AdcToEnergyReadTTreeAlg::channel_index_from_cellid(int cellid) {
  const int layer = cellid / 100000;
  const int chip = (cellid / 10000) % 10;
  const int channel = cellid % 10000;
  if (layer < 0 || layer >= AHCALGeometry::Layer_No ||
      chip < 0 || chip >= AHCALGeometry::chip_No ||
      channel < 0 || channel >= AHCALGeometry::channel_No) {
    return total_channel_count();
  }
  return (static_cast<std::size_t>(layer) * AHCALGeometry::chip_No + static_cast<std::size_t>(chip)) *
         AHCALGeometry::channel_No + static_cast<std::size_t>(channel);
}

int AdcToEnergyReadTTreeAlg::cellid_conversion(int input_cellid) {
  return AHCALGeometry::cellid_conversion(input_cellid, file_cellid_version);
}

void AdcToEnergyReadTTreeAlg::execute(EventStore &evt) { 
  auto raw_hits = evt.get<std::vector<AHCALRawHit>>(m_in_rawhit_key);

  std::vector<AHCALRecoHit> reco_hits;
  reco_hits.reserve(raw_hits.size());

  const int SwitchPoint = AHCALRefValues::SwitchPoint;

  for (const auto &raw_hit : raw_hits) {
    if (raw_hit.hittag == 0) continue; // skip hittag==0
    AHCALRecoHit reco_hit;
    reco_hit.cellID = raw_hit.cellID;
    reco_hit.index = raw_hit.index;
    const std::size_t idx = channel_index_from_cellid(raw_hit.cellID);
    if (idx >= total_channel_count()) {
      LOG_DEBUG("Invalid cellID={}, skip hit", raw_hit.cellID);
      continue;
    }
    const double mpv        = mip_values_[idx];
    const double hg_ped     = hg_ped_values_[idx];
    const double lg_ped     = lg_ped_values_[idx];
    const double gain_ratio = gainratio_values_[idx];
    const int gain_plat     = gainplat_values_[idx];

    const double hg = static_cast<double>(raw_hit.hg_adc);
    const double lg = static_cast<double>(raw_hit.lg_adc);

    if ((hg - hg_ped) < (static_cast<double>(gain_plat) - SwitchPoint)) {
      reco_hit.Nmip = (hg - hg_ped) / mpv;
      reco_hit.Edep = (hg - hg_ped) * AHCALRefValues::MIP_E / mpv;
    } else {
      reco_hit.Nmip = (lg - lg_ped) * gain_ratio / mpv;
      reco_hit.Edep = (lg - lg_ped) * gain_ratio * AHCALRefValues::MIP_E / mpv;
    }
    if (reco_hit.Edep < 0) {
      reco_hit.Edep = 0;
      reco_hit.Nmip = 0;
    }
    if (reco_hit.Nmip >1e6){
      LOG_DEBUG("Large Nmip detected: cellID={} Nmip={}", reco_hit.cellID, reco_hit.Nmip);
      LOG_DEBUG("  raw HG={} LG={} HG_ped={} LG_ped={} gain_ratio={} gain_plat={} mpv={}",
                hg, lg, hg_ped, lg_ped, gain_ratio, gain_plat, mpv);
    }
    reco_hits.push_back(reco_hit);
  }
  LOG_DEBUG("Converted {} raw hits to reco hits.", reco_hits.size());
  evt.put(m_out_recohit_key, std::move(reco_hits));
}
void AdcToEnergyReadTTreeAlg::parse_cfg(const YAML::Node& n) {
    m_cfg.in_rawhit_key = get_or<std::string>(n, "in_rawhit_key", m_cfg.in_rawhit_key);
    m_cfg.out_recohit_key = get_or<std::string>(n, "out_recohit_key", m_cfg.out_recohit_key);
    const YAML::Node mip_node = n["mip"];
    if (mip_node) {
        m_cfg.mip_file = get_or<std::string>(mip_node, "file", m_cfg.mip_file);
        m_cfg.mip_cut_string = get_or<std::string>(mip_node, "cut", m_cfg.mip_cut_string);
        m_cfg.mip_cellid_version = get_or<int>(mip_node, "cellid_version", m_cfg.mip_cellid_version);
    }
    const YAML::Node ped_node = n["pedestal"];
    if (ped_node) {
        m_cfg.ped_file = get_or<std::string>(ped_node, "file", m_cfg.ped_file);
        m_cfg.ped_cut_string = get_or<std::string>(ped_node, "cut", m_cfg.ped_cut_string);
        m_cfg.ped_cellid_version = get_or<int>(ped_node, "cellid_version", m_cfg.ped_cellid_version);
    }
    const YAML::Node dac_node = n["dac"];
    if (dac_node) {
        m_cfg.dac_file = get_or<std::string>(dac_node, "file", m_cfg.dac_file);
        m_cfg.dac_cut_string = get_or<std::string>(dac_node, "cut", m_cfg.dac_cut_string);
        m_cfg.dac_cellid_version = get_or<int>(dac_node, "cellid_version", m_cfg.dac_cellid_version);
    }
    m_in_rawhit_key = m_cfg.in_rawhit_key;
    m_out_recohit_key = m_cfg.out_recohit_key;
  }
} // namespace AHCALRecoAlg
