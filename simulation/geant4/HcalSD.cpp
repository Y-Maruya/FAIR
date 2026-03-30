#include "simulation/geant4/HcalSD.hpp"

#include "common/AHCALGeometry.hpp"

#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"

namespace AHCALRecoAlg::Sim {

HcalSD::HcalSD(std::unordered_map<int, HitRecord>* hitMap, double timeLimitNs)
    : G4VSensitiveDetector("FairHcalSD"), m_hitMap(hitMap), m_timeLimitNs(timeLimitNs) {
}

G4bool HcalSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
  if (!m_hitMap) return false;

  const double rawEdepMeV = step->GetTotalEnergyDeposit() / MeV;
  if (rawEdepMeV <= 0.0) return false;

  const auto* pre = step->GetPreStepPoint();
  const auto* post = step->GetPostStepPoint();
  if (!pre || !post) return false;
  if (step->GetTrack()->GetGlobalTime() / ns > m_timeLimitNs) return false;

  const auto* pv = pre->GetPhysicalVolume();
  if (!pv) return false;

  int layer = 0;
  const auto name = pv->GetName();
  if (name == "HCALtriggerPhysicalFront" || name == "HCALDownstreamPhysical") {
    layer = pv->GetCopyNo();
  } else {
    const auto* touch = pre->GetTouchable();
    if (!touch) return false;
    layer = touch->GetCopyNumber(2);
  }
  if (layer < 0 || layer >= AHCALGeometry::Layer_No) return false;

  const double x = pre->GetPosition().x() / mm;
  const double y = pre->GetPosition().y() / mm;
  int chip = 0;
  int channel = 0;
  AHCALGeometry::inverse(x, y, chip, channel);
  if (chip < 0 || chip >= AHCALGeometry::chip_No || channel < 0 || channel >= AHCALGeometry::channel_No) {
    return false;
  }

  const int cellID = AHCALGeometry::CellID(layer, chip, channel);
  auto& rec = (*m_hitMap)[cellID];

  const double visibleEdepMeV = birks_attenuation_(step);
  if (visibleEdepMeV <= 0.0) return false;

  rec.edepMeV += visibleEdepMeV;
  rec.nSteps += 1;

  const double midTimeNs = 0.5 * (pre->GetGlobalTime() + post->GetGlobalTime()) / ns;
  if (rawEdepMeV > rec.maxStepEdepMeV) {
    rec.maxStepEdepMeV = rawEdepMeV;
    rec.hitTimeNs = midTimeNs;
  }

  const double toaNs = step->GetTrack()->GetGlobalTime() / ns;
  if (toaNs < rec.toaNs) rec.toaNs = toaNs;

  return true;
}

double HcalSD::birks_attenuation_(const G4Step* step) const {
  auto* material = step->GetTrack()->GetMaterial();
  if (!material) return step->GetTotalEnergyDeposit() / MeV;

  const double birk1 = material->GetIonisation()->GetBirksConstant();
  const double destep = step->GetTotalEnergyDeposit();
  const double stepl = step->GetStepLength();
  const double charge = step->GetTrack()->GetDefinition()->GetPDGCharge();

  if (birk1 == 0.0 || destep <= 0.0 || stepl <= 0.0 || charge == 0.0) {
    return destep / MeV;
  }
  return (destep / (1.0 + birk1 * destep / stepl)) / MeV;
}

} // namespace AHCALRecoAlg::Sim
