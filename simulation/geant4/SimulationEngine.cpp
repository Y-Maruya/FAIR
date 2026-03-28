#include "simulation/geant4/SimulationEngine.hpp"

#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "simulation/geant4/ActionInitialization.hpp"
#include "simulation/geant4/DetectorConstruction.hpp"
#include "simulation/geant4/Pythia8DecayerPhysics.hpp"

#include "G4PhysListFactory.hh"
#include "G4RunManagerFactory.hh"
#include <algorithm>

namespace AHCALRecoAlg::Sim {

SimulationEngine::SimulationEngine(const Geant4Config& cfg) : m_cfg(cfg) {}

SimulationEngine::~SimulationEngine() = default;

void SimulationEngine::initialize() {
  if (m_initialized) return;

  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
  m_runManager.reset(rm);
  m_hitMap.reserve(AHCALGeometry::Layer_No * AHCALGeometry::chip_No * AHCALGeometry::channel_No);

  m_runManager->SetUserInitialization(new DetectorConstruction(m_cfg, &m_hitMap));

  G4PhysListFactory factory;
  auto* phys = factory.GetReferencePhysList(m_cfg.physics_list);
  if (!phys) {
    LOG_WARN("Unknown physics list '{}', fallback to FTFP_BERT", m_cfg.physics_list);
    phys = factory.GetReferencePhysList("FTFP_BERT");
  }
  m_runManager->SetUserInitialization(phys);
  if (m_cfg.enable_pythia8_decayer) {
#if FAIR_HAS_PYTHIA8
    phys->RegisterPhysics(new Pythia8DecayerPhysics());
#else
    LOG_WARN("enable_pythia8_decayer=true was requested, but FAIR_HAS_PYTHIA8 is not enabled.");
#endif
  }
  m_runManager->SetUserInitialization(new ActionInitialization(m_cfg));

  m_runManager->Initialize();
  m_initialized = true;
}

void SimulationEngine::run_one_event(std::vector<AHCALSimHit>& out_hits) {
  if (!m_initialized) {
    initialize();
  }

  m_hitMap.clear();
  m_runManager->BeamOn(1);

  out_hits.clear();
  out_hits.reserve(m_hitMap.size());
  for (const auto& kv : m_hitMap) {
    AHCALSimHit hit;
    hit.cellID = kv.first;
    hit.Edep = kv.second.edepMeV;
    hit.Nmip = hit.Edep / AHCALGeometry::MIPEnergy;
    hit.HitTime = kv.second.hitTimeNs;
    hit.TimeOfArrival = kv.second.toaNs < 1e17 ? kv.second.toaNs : kv.second.hitTimeNs;
    out_hits.push_back(hit);
  }
  std::sort(out_hits.begin(), out_hits.end(),
            [](const AHCALSimHit& a, const AHCALSimHit& b) { return a.cellID < b.cellID; });
}

} // namespace AHCALRecoAlg::Sim
