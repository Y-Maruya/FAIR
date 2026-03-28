#pragma once

#include "G4VExtDecayer.hh"

#if FAIR_HAS_PYTHIA8
#include "Pythia8/Pythia.h"
#include <memory>
#endif

namespace AHCALRecoAlg::Sim {

class Pythia8Decayer final : public G4VExtDecayer {
public:
  explicit Pythia8Decayer(const std::string& name);
  ~Pythia8Decayer() override = default;

  G4DecayProducts* ImportDecayProducts(const G4Track& aTrack) override;

private:
#if FAIR_HAS_PYTHIA8
  std::unique_ptr<Pythia8::Pythia> m_decayer;
#endif
};

} // namespace AHCALRecoAlg::Sim
