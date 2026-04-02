#include "Geant4SimReader.hpp"

#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "simulation/geant4/PrimaryGeneratorAction.hpp"

#include <cmath>
#include <stdexcept>

Geant4SimReader::Geant4SimReader(const YAML::Node& cfg) {
  m_cfg.out_simhits_key = get_or<std::string>(cfg, "out_simhits_key", m_cfg.out_simhits_key);
  m_cfg.out_simdata_key = get_or<std::string>(cfg, "out_simdata_key", m_cfg.out_simdata_key);
  m_cfg.max_events = get_or<int>(cfg, "max_events", m_cfg.max_events);
  m_cfg.seed = get_or<int>(cfg, "seed", m_cfg.seed);

  if (cfg["detector"]) {
    m_cfg.g4.detector_model = get_or<std::string>(cfg["detector"], "model", m_cfg.g4.detector_model);
    m_cfg.g4.use_ahcal_exact_geometry = get_or<bool>(cfg["detector"], "use_ahcal_exact_geometry", m_cfg.g4.use_ahcal_exact_geometry);
    m_cfg.g4.ahcal_geometry_gdml_file = get_or<std::string>(cfg["detector"], "ahcal_geometry_gdml_file", m_cfg.g4.ahcal_geometry_gdml_file);
  }
  if (cfg["decayer"]) {
    m_cfg.g4.enable_pythia8_decayer = get_or<bool>(cfg["decayer"], "enable_pythia8", m_cfg.g4.enable_pythia8_decayer);
  }
  if (cfg["generator"]) {
    const auto& g = cfg["generator"];
    m_cfg.g4.generator_mode = get_or<std::string>(g, "mode", m_cfg.g4.generator_mode);
    const std::string p = get_or<std::string>(g, "particle", "mu-");
    m_cfg.g4.primary_pdg = get_or<int>(g, "pdg", pdg_from_name_(p));
    m_cfg.g4.primary_energy_GeV = get_or<double>(g, "energy_GeV", m_cfg.g4.primary_energy_GeV);

    if (g["position_mm"] && g["position_mm"].IsSequence() && g["position_mm"].size() == 3) {
      m_cfg.g4.primary_x_mm = g["position_mm"][0].as<double>();
      m_cfg.g4.primary_y_mm = g["position_mm"][1].as<double>();
      m_cfg.g4.primary_z_mm = g["position_mm"][2].as<double>();
    }
    if (g["direction"] && g["direction"].IsSequence() && g["direction"].size() == 3) {
      m_cfg.g4.dir_x = g["direction"][0].as<double>();
      m_cfg.g4.dir_y = g["direction"][1].as<double>();
      m_cfg.g4.dir_z = g["direction"][2].as<double>();
    }
    m_cfg.g4.input_root_file = get_or<std::string>(g, "input_root_file", m_cfg.g4.input_root_file);
    m_cfg.g4.input_tree_name = get_or<std::string>(g, "input_tree_name", m_cfg.g4.input_tree_name);
    m_cfg.g4.fixed_vertex = get_or<bool>(g, "fixed_vertex", m_cfg.g4.fixed_vertex);
  }

  if (cfg["physics"]) {
    m_cfg.g4.physics_list = get_or<std::string>(cfg["physics"], "list", m_cfg.g4.physics_list);
  }
  if (cfg["detector"]) {
    m_cfg.g4.birks_constant_mm_per_MeV = get_or<double>(cfg["detector"], "birks_constant_mm_per_MeV", m_cfg.g4.birks_constant_mm_per_MeV);
    m_cfg.g4.enable_trigger_component = get_or<bool>(cfg["detector"], "enable_trigger_component", m_cfg.g4.enable_trigger_component);
    m_cfg.g4.trigger_nplanes = get_or<int>(cfg["detector"], "trigger_nplanes", m_cfg.g4.trigger_nplanes);
    m_cfg.g4.hcal_step_time_limit_ns = get_or<double>(cfg["detector"], "hcal_step_time_limit_ns", m_cfg.g4.hcal_step_time_limit_ns);
    m_cfg.g4.passive_side_thickness_mm = get_or<double>(cfg["detector"], "passive_side_thickness_mm", m_cfg.g4.passive_side_thickness_mm);
    m_cfg.g4.passive_cover_thickness_mm = get_or<double>(cfg["detector"], "passive_cover_thickness_mm", m_cfg.g4.passive_cover_thickness_mm);
    m_cfg.g4.attach_thickness_mm = get_or<double>(cfg["detector"], "attach_thickness_mm", m_cfg.g4.attach_thickness_mm);
    m_cfg.g4.double_sided_readout = get_or<bool>(cfg["detector"], "double_sided_readout", m_cfg.g4.double_sided_readout);
    m_cfg.g4.house_pitch_mm = get_or<double>(cfg["detector"], "house_pitch_mm", m_cfg.g4.house_pitch_mm);
    m_cfg.g4.sensitive_dig_out_x_mm = get_or<double>(cfg["detector"], "sensitive_dig_out_x_mm", m_cfg.g4.sensitive_dig_out_x_mm);
    m_cfg.g4.sensitive_dig_out_y_mm = get_or<double>(cfg["detector"], "sensitive_dig_out_y_mm", m_cfg.g4.sensitive_dig_out_y_mm);
    m_cfg.g4.sensitive_dig_out_z_mm = get_or<double>(cfg["detector"], "sensitive_dig_out_z_mm", m_cfg.g4.sensitive_dig_out_z_mm);
  }

  if (std::abs(m_cfg.g4.dir_z) < 1e-9 && std::abs(m_cfg.g4.dir_x) < 1e-9 && std::abs(m_cfg.g4.dir_y) < 1e-9) {
    throw std::runtime_error("Geant4SimReader: invalid zero direction vector");
  }

  m_engine = std::make_unique<AHCALRecoAlg::Sim::SimulationEngine>(m_cfg.g4);
  m_engine->initialize();

  LOG_INFO("Geant4SimReader initialized (detector='{}', generator='{}', physics='{}', root='{}')",
           m_cfg.g4.detector_model, m_cfg.g4.generator_mode, m_cfg.g4.physics_list, m_cfg.g4.input_root_file);
}

bool Geant4SimReader::next(std::vector<AHCALSimHit>& out_hits, SimData& out_simdata) {
  if (m_cfg.max_events >= 0 && m_event >= m_cfg.max_events) {
    return false;
  }

  fill_simdata_(out_simdata);
  m_engine->run_one_event(out_hits);

  ++m_event;
  return true;
}

void Geant4SimReader::fill_simdata_(SimData& out) const {
  const auto& info = AHCALRecoAlg::Sim::SimReaderPrimaryGenerator::last_event_info();

  out.injected_pdgId = info.injected_pdg != 0 ? info.injected_pdg : m_cfg.g4.primary_pdg;
  out.injected_energy = info.injected_energy_GeV > 0 ? info.injected_energy_GeV : m_cfg.g4.primary_energy_GeV;
  out.injected_time = 0.0;
  out.injected_z = info.injected_z_mm;
  out.injected_x = info.injected_x_mm;
  out.injected_y = info.injected_y_mm;

  out.injected_px = info.injected_px_GeV;
  out.injected_py = info.injected_py_GeV;
  out.injected_pz = info.injected_pz_GeV;

  out.ftagNulabel = info.ftagNulabel;
  out.isNeutrinoEvent = (info.ftagNulabel >= 0 && info.ftagNulabel <= 3);
  out.isCC = (info.ftagNulabel == 0 || info.ftagNulabel == 1 || info.ftagNulabel == 2);
  out.Interaction_x = info.interaction_x_mm;
  out.Interaction_y = info.interaction_y_mm;
  out.Interaction_z = info.interaction_z_mm;
  out.Secondarypdgid = info.secondary_pdg;
  out.SecondaryEnergy = info.secondary_energy_GeV;
  out.Secondary_px = info.secondary_px_GeV;
  out.Secondary_py = info.secondary_py_GeV;
  out.Secondary_pz = info.secondary_pz_GeV;
}

int Geant4SimReader::pdg_from_name_(const std::string& name) {
  if (name == "mu-") return 13;
  if (name == "mu+") return -13;
  if (name == "e-") return 11;
  if (name == "e+") return -11;
  if (name == "pi+") return 211;
  if (name == "pi-") return -211;
  if (name == "proton") return 2212;
  if (name == "neutron") return 2112;
  return 0;
}
