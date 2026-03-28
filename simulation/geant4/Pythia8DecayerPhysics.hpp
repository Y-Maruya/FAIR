#pragma once

#include "G4VPhysicsConstructor.hh"

namespace AHCALRecoAlg::Sim {

class Pythia8DecayerPhysics final : public G4VPhysicsConstructor {
public:
  explicit Pythia8DecayerPhysics(G4int ver = 1);
  ~Pythia8DecayerPhysics() override;

  void ConstructParticle() override;
  void ConstructProcess() override;
};

} // namespace AHCALRecoAlg::Sim
