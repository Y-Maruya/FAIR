#include "simulation/geant4/HcalSD.hpp"

#include "common/AHCALGeometry.hpp"

#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"

namespace AHCALRecoAlg::Sim {

HcalSD::HcalSD(std::unordered_map<int, HitRecord>* hitMap)
    : G4VSensitiveDetector("FairHcalSD"), m_hitMap(hitMap) {}

G4bool HcalSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
  if (!m_hitMap) return false;

  const auto* pre = step->GetPreStepPoint();
  const auto* touch = pre->GetTouchable();
  if (!touch) return false;

  const int layer = touch->GetCopyNumber();
  if (layer < 0 || layer >= AHCALGeometry::Layer_No) return false;

  const double x = pre->GetPosition().x() / mm;
  const double y = pre->GetPosition().y() / mm;
  int chip = 0;
  int channel = 0;
  AHCALGeometry::inverse(x, y, chip, channel);

  if (chip < 0 || chip >= AHCALGeometry::chip_No || channel < 0 || channel >= AHCALGeometry::channel_No) {
    return false;
  }

  const double edep = step->GetTotalEnergyDeposit() / MeV;
  if (edep <= 0.0) return false;

  const int cellID = AHCALGeometry::CellID(layer, chip, channel);
  const double t = pre->GetGlobalTime() / ns;

  auto& rec = (*m_hitMap)[cellID];
  rec.edepMeV += edep;
  if (t > rec.hitTimeNs) rec.hitTimeNs = t;
  if (t < rec.toaNs) rec.toaNs = t;

  return true;
}

} // namespace AHCALRecoAlg::Sim
