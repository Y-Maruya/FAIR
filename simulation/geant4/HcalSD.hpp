#pragma once

#include "simulation/geant4/HitRecord.hpp"
#include <unordered_map>

#include "G4VSensitiveDetector.hh"

namespace AHCALRecoAlg::Sim {

class HcalSD final : public G4VSensitiveDetector {
public:
  explicit HcalSD(std::unordered_map<int, HitRecord>* hitMap);

  G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

private:
  std::unordered_map<int, HitRecord>* m_hitMap;
};

} // namespace AHCALRecoAlg::Sim
