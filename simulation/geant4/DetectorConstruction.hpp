#pragma once

#include "simulation/geant4/Geant4Config.hpp"
#include "simulation/geant4/HitRecord.hpp"

#include "G4VUserDetectorConstruction.hh"

#include <unordered_map>

class G4VPhysicalVolume;

namespace AHCALRecoAlg::Sim {

class DetectorConstruction final : public G4VUserDetectorConstruction {
public:
  DetectorConstruction(const Geant4Config& cfg, std::unordered_map<int, HitRecord>* hitMap);

  G4VPhysicalVolume* Construct() override;

private:
  Geant4Config m_cfg;
  std::unordered_map<int, HitRecord>* m_hitMap;
};

} // namespace AHCALRecoAlg::Sim
