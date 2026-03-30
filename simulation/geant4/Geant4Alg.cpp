#include "simulation/geant4/Geant4Alg.hpp"

#include "common/AlgRegistry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/edm/InjectedParticle.hpp"
#include "common/edm/SimData.hpp"
#include "common/edm/SimHit.hpp"
#include "simulation/geant4/PrimaryGeneratorAction.hpp"
#include "simulation/geant4/SimulationEngine.hpp"

AHCAL_REGISTER_ALG(AHCALRecoAlg::Sim::Geant4Alg, "Geant4Alg")

namespace AHCALRecoAlg::Sim {

Geant4Alg::Geant4Alg(RunContext& ctx, std::string name)
    : IAlg(ctx, std::move(name)) {}

Geant4Alg::~Geant4Alg() = default;

void Geant4Alg::parse_cfg(const YAML::Node& n) {
  m_cfg.in_injected_particles_key = get_or<std::string>(n, "in_injected_particles_key", m_cfg.in_injected_particles_key);
  m_cfg.out_simhits_key = get_or<std::string>(n, "out_simhits_key", m_cfg.out_simhits_key);
  m_cfg.out_simdata_key = get_or<std::string>(n, "out_simdata_key", m_cfg.out_simdata_key);
  m_cfg.require_input_particle = get_or<bool>(n, "require_input_particle", m_cfg.require_input_particle);

  m_cfg.g4.detector_model = get_or<std::string>(n, "detector_model", m_cfg.g4.detector_model);
  m_cfg.g4.generator_mode = get_or<std::string>(n, "generator_mode", "external_inject");
  m_cfg.g4.physics_list = get_or<std::string>(n, "physics_list", m_cfg.g4.physics_list);
  m_cfg.g4.birks_constant_mm_per_MeV = get_or<double>(n, "birks_constant_mm_per_MeV", m_cfg.g4.birks_constant_mm_per_MeV);
  m_cfg.g4.enable_trigger_component = get_or<bool>(n, "enable_trigger_component", m_cfg.g4.enable_trigger_component);
  m_cfg.g4.trigger_nplanes = get_or<int>(n, "trigger_nplanes", m_cfg.g4.trigger_nplanes);
  m_cfg.g4.hcal_step_time_limit_ns = get_or<double>(n, "hcal_step_time_limit_ns", m_cfg.g4.hcal_step_time_limit_ns);
  m_cfg.g4.enable_pythia8_decayer = get_or<bool>(n, "enable_pythia8_decayer", m_cfg.g4.enable_pythia8_decayer);
  m_cfg.g4.passive_side_thickness_mm = get_or<double>(n, "passive_side_thickness_mm", m_cfg.g4.passive_side_thickness_mm);
  m_cfg.g4.passive_cover_thickness_mm = get_or<double>(n, "passive_cover_thickness_mm", m_cfg.g4.passive_cover_thickness_mm);
  m_cfg.g4.attach_thickness_mm = get_or<double>(n, "attach_thickness_mm", m_cfg.g4.attach_thickness_mm);
  m_cfg.g4.double_sided_readout = get_or<bool>(n, "double_sided_readout", m_cfg.g4.double_sided_readout);
  m_cfg.g4.house_pitch_mm = get_or<double>(n, "house_pitch_mm", m_cfg.g4.house_pitch_mm);
  m_cfg.g4.sensitive_dig_out_x_mm = get_or<double>(n, "sensitive_dig_out_x_mm", m_cfg.g4.sensitive_dig_out_x_mm);
  m_cfg.g4.sensitive_dig_out_y_mm = get_or<double>(n, "sensitive_dig_out_y_mm", m_cfg.g4.sensitive_dig_out_y_mm);
  m_cfg.g4.sensitive_dig_out_z_mm = get_or<double>(n, "sensitive_dig_out_z_mm", m_cfg.g4.sensitive_dig_out_z_mm);
}

void Geant4Alg::initialize() {
  m_engine = std::make_unique<SimulationEngine>(m_cfg.g4);
  m_engine->initialize();
}

bool Geant4Alg::apply_external_primary_(const std::vector<InjectedParticle>& injected) {
  if (injected.empty()) {
    return false;
  }

  const auto& p = injected.front();
  ExternalPrimaryParticle external;
  external.pdg = p.pdgId;
  external.energy_GeV = p.energy_GeV;
  external.time_ns = p.time_ns;
  external.x_mm = p.x_mm;
  external.y_mm = p.y_mm;
  external.z_mm = p.z_mm;
  external.px_GeV = p.px_GeV;
  external.py_GeV = p.py_GeV;
  external.pz_GeV = p.pz_GeV;
  SimReaderPrimaryGenerator::set_external_primary(external);
  return true;
}

void Geant4Alg::fill_simdata_(SimData& out) const {
  const auto& info = SimReaderPrimaryGenerator::last_event_info();
  out.injected_pdgId = info.injected_pdg;
  out.injected_energy = info.injected_energy_GeV;
  out.injected_time = 0.0;
  out.injected_x = info.injected_x_mm;
  out.injected_y = info.injected_y_mm;
  out.injected_z = info.injected_z_mm;
  out.injected_px = info.injected_px_GeV;
  out.injected_py = info.injected_py_GeV;
  out.injected_pz = info.injected_pz_GeV;

  out.ftagNulabel = info.ftagNulabel;
  out.isCC = (out.ftagNulabel >= 0 && out.ftagNulabel <= 2);
  out.isNeutrinoEvent = (out.ftagNulabel >= 0 && out.ftagNulabel <= 3);
  out.Interaction_x = info.interaction_x_mm;
  out.Interaction_y = info.interaction_y_mm;
  out.Interaction_z = info.interaction_z_mm;

  out.Secondarypdgid = info.secondary_pdg;
  out.SecondaryEnergy = info.secondary_energy_GeV;
  out.Secondary_px = info.secondary_px_GeV;
  out.Secondary_py = info.secondary_py_GeV;
  out.Secondary_pz = info.secondary_pz_GeV;
}

void Geant4Alg::execute(EventStore& evt) {
  auto* injected = evt.try_get<std::vector<InjectedParticle>>(m_cfg.in_injected_particles_key);
  if ((!injected || injected->empty()) && m_cfg.require_input_particle) {
    throw std::runtime_error("Geant4Alg: required input key is missing or empty: " + m_cfg.in_injected_particles_key);
  }

  if (injected && !injected->empty()) {
    apply_external_primary_(*injected);
  }

  std::vector<AHCALSimHit> out_hits;
  SimData out_simdata;
  m_engine->run_one_event(out_hits);
  fill_simdata_(out_simdata);

  evt.set(m_cfg.out_simhits_key, std::move(out_hits));
  evt.set(m_cfg.out_simdata_key, std::move(out_simdata));
}

} // namespace AHCALRecoAlg::Sim
