#pragma once

#include "simulation/geant4/Geant4Config.hpp"

#include "G4VUserActionInitialization.hh"

namespace AHCALRecoAlg::Sim {

class ActionInitialization final : public G4VUserActionInitialization {
public:
  explicit ActionInitialization(const Geant4Config& cfg);

  void BuildForMaster() const override;
  void Build() const override;

private:
  Geant4Config m_cfg;
};

} // namespace AHCALRecoAlg::Sim
