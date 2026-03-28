#include "simulation/geant4/ActionInitialization.hpp"

#include "simulation/geant4/PrimaryGeneratorAction.hpp"

namespace AHCALRecoAlg::Sim {

ActionInitialization::ActionInitialization(const Geant4Config& cfg) : m_cfg(cfg) {}

void ActionInitialization::BuildForMaster() const {}

void ActionInitialization::Build() const {
  SetUserAction(new SimReaderPrimaryGenerator(m_cfg));
}

} // namespace AHCALRecoAlg::Sim
