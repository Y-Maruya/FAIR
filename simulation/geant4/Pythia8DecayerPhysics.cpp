#include "simulation/geant4/Pythia8DecayerPhysics.hpp"

#include "simulation/geant4/Pythia8Decayer.hpp"

#include "G4Decay.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include <algorithm>
#include <vector>

namespace {
bool is_target_particle(int pdgCode) {
  static const std::vector<int> targetPDGCodes = {
      15, 521, 511, 541, 531, 411, 421, 431, 443, 553,
      5122, 4122, 5332, 4322, 4324, 5232, 5132, 4232, 4132};
  return std::find(targetPDGCodes.begin(), targetPDGCodes.end(), std::abs(pdgCode)) != targetPDGCodes.end();
}
} // namespace

namespace AHCALRecoAlg::Sim {

Pythia8DecayerPhysics::Pythia8DecayerPhysics(G4int)
    : G4VPhysicsConstructor("Pythia8DecayerPhysics") {}

Pythia8DecayerPhysics::~Pythia8DecayerPhysics() = default;

void Pythia8DecayerPhysics::ConstructParticle() {}

void Pythia8DecayerPhysics::ConstructProcess() {
#if FAIR_HAS_PYTHIA8
  auto* extDecayer = new Pythia8Decayer("Pythia8DecayerPhysics");

  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)()) {
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (!pmanager) continue;

    G4ProcessVector* processVector = pmanager->GetProcessList();
    for (size_t i = 0; i < processVector->length(); ++i) {
      auto* decay = dynamic_cast<G4Decay*>((*processVector)[i]);
      if (!decay) continue;

      if (is_target_particle(std::abs(particle->GetPDGEncoding()))) {
        if (particle->GetDecayTable()) {
          delete particle->GetDecayTable();
          particle->SetDecayTable(nullptr);
        }
      }
      if (!particle->GetDecayTable()) {
        decay->SetExtDecayer(extDecayer);
      }
    }
  }
#endif
}

} // namespace AHCALRecoAlg::Sim
