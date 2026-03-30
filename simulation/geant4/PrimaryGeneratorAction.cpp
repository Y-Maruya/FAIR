#include "simulation/geant4/PrimaryGeneratorAction.hpp"

#include "common/Logger.hpp"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include "TFile.h"
#include "TTree.h"

#include <stdexcept>

namespace AHCALRecoAlg::Sim {
thread_local PrimaryEventInfo SimReaderPrimaryGenerator::s_lastInfo{};
thread_local bool SimReaderPrimaryGenerator::s_hasExternalPrimary = false;
thread_local ExternalPrimaryParticle SimReaderPrimaryGenerator::s_externalPrimary{};

SimReaderPrimaryGenerator::SimReaderPrimaryGenerator(const Geant4Config& cfg) : m_cfg(cfg) {
  m_gun = new G4ParticleGun(1);

  auto* p = G4ParticleTable::GetParticleTable()->FindParticle(m_cfg.primary_pdg);
  if (!p) p = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
  m_gun->SetParticleDefinition(p);
  m_gun->SetParticleEnergy(m_cfg.primary_energy_GeV * GeV);
  m_gun->SetParticleMomentumDirection(G4ThreeVector(m_cfg.dir_x, m_cfg.dir_y, m_cfg.dir_z).unit());
  m_gun->SetParticlePosition(G4ThreeVector(m_cfg.primary_x_mm * mm, m_cfg.primary_y_mm * mm, m_cfg.primary_z_mm * mm));

  if (m_cfg.generator_mode == "sim_reader") {
    setup_root_reader_();
  }
}

SimReaderPrimaryGenerator::~SimReaderPrimaryGenerator() { delete m_gun; }

void SimReaderPrimaryGenerator::setup_root_reader_() {
  if (m_cfg.input_root_file.empty()) {
    throw std::runtime_error("SimReaderPrimaryGenerator: generator_mode=sim_reader requires generator.input_root_file");
  }

  m_file.reset(TFile::Open(m_cfg.input_root_file.c_str(), "READ"));
  if (!m_file || m_file->IsZombie()) {
    throw std::runtime_error("SimReaderPrimaryGenerator: failed to open input ROOT file: " + m_cfg.input_root_file);
  }

  m_tree = dynamic_cast<TTree*>(m_file->Get(m_cfg.input_tree_name.c_str()));
  if (!m_tree) {
    throw std::runtime_error("SimReaderPrimaryGenerator: cannot find tree: " + m_cfg.input_tree_name);
  }

  m_tree->SetBranchStatus("*", 0);
  m_tree->SetBranchStatus("pdgc", 1);
  m_tree->SetBranchStatus("px", 1);
  m_tree->SetBranchStatus("py", 1);
  m_tree->SetBranchStatus("pz", 1);
  m_tree->SetBranchStatus("E", 1);
  m_tree->SetBranchStatus("status", 1);
  m_tree->SetBranchStatus("firstMother", 1);

  m_tree->SetBranchAddress("pdgc", &m_pdgc);
  m_tree->SetBranchAddress("px", &m_px);
  m_tree->SetBranchAddress("py", &m_py);
  m_tree->SetBranchAddress("pz", &m_pz);
  m_tree->SetBranchAddress("E", &m_E);
  m_tree->SetBranchAddress("status", &m_status);
  m_tree->SetBranchAddress("firstMother", &m_firstMother);

  if (!m_cfg.fixed_vertex) {
    m_tree->SetBranchStatus("vx", 1);
    m_tree->SetBranchStatus("vy", 1);
    m_tree->SetBranchStatus("vz", 1);
    m_tree->SetBranchAddress("vx", &m_vx);
    m_tree->SetBranchAddress("vy", &m_vy);
    m_tree->SetBranchAddress("vz", &m_vz);
  }
}

const PrimaryEventInfo& SimReaderPrimaryGenerator::last_event_info() { return s_lastInfo; }

void SimReaderPrimaryGenerator::set_external_primary(const ExternalPrimaryParticle& primary) {
  s_externalPrimary = primary;
  s_hasExternalPrimary = true;
}

void SimReaderPrimaryGenerator::clear_external_primary() {
  s_hasExternalPrimary = false;
}

void SimReaderPrimaryGenerator::update_interaction_info_() {
  s_lastInfo.ftagNulabel = 4;
  s_lastInfo.secondary_pdg = 0;
  s_lastInfo.secondary_energy_GeV = -999999;
  s_lastInfo.secondary_px_GeV = -999999;
  s_lastInfo.secondary_py_GeV = -999999;
  s_lastInfo.secondary_pz_GeV = -999999;
  if (!m_pdgc || !m_status || !m_firstMother || m_pdgc->empty()) return;

  const int nu = m_pdgc->at(0);
  for (size_t j = 2; j < m_pdgc->size(); ++j) {
    if (m_status->at(j) != 1 || m_firstMother->at(j) != 0) continue;
    if (std::abs(nu) == 12 && m_pdgc->at(j) == (nu > 0 ? 11 : -11)) {
      s_lastInfo.ftagNulabel = 0;
    } else if (std::abs(nu) == 14 && m_pdgc->at(j) == (nu > 0 ? 13 : -13)) {
      s_lastInfo.ftagNulabel = 1;
    } else if (std::abs(nu) == 16 && m_pdgc->at(j) == (nu > 0 ? 15 : -15)) {
      s_lastInfo.ftagNulabel = 2;
    } else if (m_pdgc->at(j) == nu) {
      s_lastInfo.ftagNulabel = 3;
    } else {
      continue;
    }
    s_lastInfo.secondary_pdg = m_pdgc->at(j);
    if (m_E && j < m_E->size()) s_lastInfo.secondary_energy_GeV = m_E->at(j);
    if (m_px && j < m_px->size()) s_lastInfo.secondary_px_GeV = m_px->at(j);
    if (m_py && j < m_py->size()) s_lastInfo.secondary_py_GeV = m_py->at(j);
    if (m_pz && j < m_pz->size()) s_lastInfo.secondary_pz_GeV = m_pz->at(j);
    return;
  }
}

bool SimReaderPrimaryGenerator::generate_from_root_(G4Event* event) {
  if (!m_tree) return false;
  if (m_evt >= m_tree->GetEntries()) return false;

  m_tree->GetEntry(m_evt);
  update_interaction_info_();

  G4ThreeVector vtx(m_cfg.primary_x_mm * mm, m_cfg.primary_y_mm * mm, m_cfg.primary_z_mm * mm);
  if (!m_cfg.fixed_vertex) {
    // AHCAL-simulation follows meters in gFaser for v*.
    vtx = G4ThreeVector(m_vx * m, m_vy * m, m_vz * m);
  }
  s_lastInfo.injected_x_mm = vtx.x() / mm;
  s_lastInfo.injected_y_mm = vtx.y() / mm;
  s_lastInfo.injected_z_mm = vtx.z() / mm;
  s_lastInfo.interaction_x_mm = s_lastInfo.injected_x_mm;
  s_lastInfo.interaction_y_mm = s_lastInfo.injected_y_mm;
  s_lastInfo.interaction_z_mm = s_lastInfo.injected_z_mm;
  if (m_pdgc && !m_pdgc->empty()) s_lastInfo.injected_pdg = m_pdgc->at(0);
  if (m_E && !m_E->empty()) s_lastInfo.injected_energy_GeV = m_E->at(0);

  auto* vertex = new G4PrimaryVertex(vtx, 0.0 * ns);

  if (m_pdgc && m_status && m_px && m_py && m_pz) {
    for (size_t i = 0; i < m_pdgc->size(); ++i) {
      if (m_status->at(i) != 1) continue;
      auto* p = new G4PrimaryParticle(m_pdgc->at(i), m_px->at(i) * GeV, m_py->at(i) * GeV, m_pz->at(i) * GeV);
      vertex->SetPrimary(p);
      if (i == 0) {
        s_lastInfo.injected_px_GeV = m_px->at(i);
        s_lastInfo.injected_py_GeV = m_py->at(i);
        s_lastInfo.injected_pz_GeV = m_pz->at(i);
      }
    }
  }

  event->AddPrimaryVertex(vertex);
  ++m_evt;
  return true;
}

void SimReaderPrimaryGenerator::generate_particle_gun_(G4Event* event) {
  s_lastInfo = {};
  s_lastInfo.injected_pdg = m_cfg.primary_pdg;
  s_lastInfo.injected_energy_GeV = m_cfg.primary_energy_GeV;
  s_lastInfo.injected_x_mm = m_cfg.primary_x_mm;
  s_lastInfo.injected_y_mm = m_cfg.primary_y_mm;
  s_lastInfo.injected_z_mm = m_cfg.primary_z_mm;
  s_lastInfo.interaction_x_mm = m_cfg.primary_x_mm;
  s_lastInfo.interaction_y_mm = m_cfg.primary_y_mm;
  s_lastInfo.interaction_z_mm = m_cfg.primary_z_mm;
  s_lastInfo.injected_px_GeV = m_cfg.primary_energy_GeV * m_cfg.dir_x;
  s_lastInfo.injected_py_GeV = m_cfg.primary_energy_GeV * m_cfg.dir_y;
  s_lastInfo.injected_pz_GeV = m_cfg.primary_energy_GeV * m_cfg.dir_z;
  s_lastInfo.ftagNulabel = -1;
  m_gun->GeneratePrimaryVertex(event);
  ++m_evt;
}

void SimReaderPrimaryGenerator::generate_external_primary_(G4Event* event) {
  if (!s_hasExternalPrimary) {
    throw std::runtime_error("SimReaderPrimaryGenerator: external_inject mode requires an injected primary");
  }

  const auto p = s_externalPrimary;
  s_lastInfo = {};
  s_lastInfo.injected_pdg = p.pdg;
  s_lastInfo.injected_energy_GeV = p.energy_GeV;
  s_lastInfo.injected_x_mm = p.x_mm;
  s_lastInfo.injected_y_mm = p.y_mm;
  s_lastInfo.injected_z_mm = p.z_mm;
  s_lastInfo.interaction_x_mm = p.x_mm;
  s_lastInfo.interaction_y_mm = p.y_mm;
  s_lastInfo.interaction_z_mm = p.z_mm;
  s_lastInfo.injected_px_GeV = p.px_GeV;
  s_lastInfo.injected_py_GeV = p.py_GeV;
  s_lastInfo.injected_pz_GeV = p.pz_GeV;
  s_lastInfo.ftagNulabel = -1;

  auto* vertex = new G4PrimaryVertex(G4ThreeVector(p.x_mm * mm, p.y_mm * mm, p.z_mm * mm), p.time_ns * ns);
  auto* primary = new G4PrimaryParticle(p.pdg, p.px_GeV * GeV, p.py_GeV * GeV, p.pz_GeV * GeV);
  vertex->SetPrimary(primary);
  event->AddPrimaryVertex(vertex);

  s_hasExternalPrimary = false;
  ++m_evt;
}

void SimReaderPrimaryGenerator::GeneratePrimaries(G4Event* event) {
  if (m_cfg.generator_mode == "external_inject") {
    generate_external_primary_(event);
    return;
  }

  if (m_cfg.generator_mode == "sim_reader") {
    if (generate_from_root_(event)) return;
    throw std::runtime_error("SimReaderPrimaryGenerator: ROOT input exhausted");
  }
  generate_particle_gun_(event);
}

} // namespace AHCALRecoAlg::Sim
