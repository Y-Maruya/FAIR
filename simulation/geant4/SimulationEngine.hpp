#pragma once

#include "common/edm/SimHit.hpp"
#include "simulation/geant4/Geant4Config.hpp"
#include "simulation/geant4/HitRecord.hpp"

#include <memory>
#include <unordered_map>
#include <vector>

class G4RunManager;

namespace AHCALRecoAlg::Sim {

class SimulationEngine {
public:
  explicit SimulationEngine(const Geant4Config& cfg);
  ~SimulationEngine();

  void initialize();
  void run_one_event(std::vector<AHCALSimHit>& out_hits);

private:
  Geant4Config m_cfg;
  std::unique_ptr<G4RunManager> m_runManager;
  std::unordered_map<int, HitRecord> m_hitMap;
  bool m_initialized = false;
};

} // namespace AHCALRecoAlg::Sim
