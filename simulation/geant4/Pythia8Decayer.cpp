#include "simulation/geant4/Pythia8Decayer.hpp"

#include "common/Logger.hpp"

#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

#include <cmath>
#include <stdexcept>

namespace AHCALRecoAlg::Sim {

Pythia8Decayer::Pythia8Decayer(const std::string& name)
    : G4VExtDecayer(name) {
#if FAIR_HAS_PYTHIA8
  std::string docstring = "../share/Pythia8/xmldoc";
  m_decayer = std::make_unique<Pythia8::Pythia>(docstring);

  m_decayer->readString("ProcessLevel:all = off");
  m_decayer->readString("ProcessLevel:resonanceDecays=on");
  m_decayer->readString("Init:showAllSettings=false");
  m_decayer->readString("Init:showChangedSettings=false");
  m_decayer->readString("Init:showAllParticleData=false");
  m_decayer->readString("Init:showChangedParticleData=false");
  m_decayer->readString("Next:numberShowProcess = 0");
  m_decayer->readString("Next:numberShowEvent = 0");
  m_decayer->init();
  m_decayer->readString("111:onMode = off");
#else
  (void)name;
  throw std::runtime_error("Pythia8Decayer was constructed but FAIR_HAS_PYTHIA8 is not enabled.");
#endif
}

G4DecayProducts* Pythia8Decayer::ImportDecayProducts(const G4Track& aTrack) {
#if FAIR_HAS_PYTHIA8
  m_decayer->event.reset();

  G4ParticleDefinition* pd = aTrack.GetDefinition();
  const int pdgid = pd->GetPDGEncoding();

  if (!m_decayer->particleData.findParticle(pdgid)) {
    return nullptr;
  }
  if (!m_decayer->particleData.canDecay(pdgid)) {
    return nullptr;
  }

  m_decayer->event.append(
      pdgid, 1, 0, 0,
      aTrack.GetMomentum().x() / CLHEP::GeV,
      aTrack.GetMomentum().y() / CLHEP::GeV,
      aTrack.GetMomentum().z() / CLHEP::GeV,
      aTrack.GetDynamicParticle()->GetTotalEnergy() / CLHEP::GeV,
      pd->GetPDGMass() / CLHEP::GeV);

  const double spinup = std::round(std::cos(aTrack.GetPolarization().angle(aTrack.GetMomentumDirection())));
  m_decayer->event.back().pol(spinup);

  const int before = m_decayer->event.size();
  m_decayer->next();
  const int after = m_decayer->event.size();

  auto* dproducts = new G4DecayProducts(*(aTrack.GetDynamicParticle()));
  for (int ip = before; ip < after; ++ip) {
    if (m_decayer->event[ip].status() < 0) continue;

    G4ParticleDefinition* pddec = G4ParticleTable::GetParticleTable()->FindParticle(m_decayer->event[ip].id());
    if (!pddec) continue;

    G4ThreeVector mom(
        m_decayer->event[ip].px() * CLHEP::GeV,
        m_decayer->event[ip].py() * CLHEP::GeV,
        m_decayer->event[ip].pz() * CLHEP::GeV);
    dproducts->PushProducts(new G4DynamicParticle(pddec, mom));
  }
  return dproducts;
#else
  (void)aTrack;
  return nullptr;
#endif
}

} // namespace AHCALRecoAlg::Sim
