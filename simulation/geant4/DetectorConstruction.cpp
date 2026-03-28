#include "simulation/geant4/DetectorConstruction.hpp"

#include "common/AHCALGeometry.hpp"
#include "simulation/geant4/HcalSD.hpp"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

namespace AHCALRecoAlg::Sim {

DetectorConstruction::DetectorConstruction(const Geant4Config& cfg, std::unordered_map<int, HitRecord>* hitMap)
    : m_cfg(cfg), m_hitMap(hitMap) {}

G4VPhysicalVolume* DetectorConstruction::Construct() {
  auto* nist = G4NistManager::Instance();
  auto* worldMat = nist->FindOrBuildMaterial("G4_AIR");
  auto* sciMat = nist->BuildMaterialWithNewDensity("FAIR_POLY", "G4_POLYSTYRENE", 1.032 * g / cm3);
  sciMat->GetIonisation()->SetBirksConstant(m_cfg.birks_constant_mm_per_MeV * mm / MeV);
  auto* steelMat = nist->FindOrBuildMaterial("G4_Fe");

  const double worldXY = 2600.0 * mm;
  const double worldZ = 6000.0 * mm;
  auto* worldSolid = new G4Box("World", worldXY / 2, worldXY / 2, worldZ / 2);
  auto* worldLV = new G4LogicalVolume(worldSolid, worldMat, "WorldLV");
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "WorldPV", nullptr, false, 0, false);

  const double slabXY = AHCALGeometry::xy_size * 18.0 * mm;
  const double activeThick = 3.0 * mm;
  const double absorberThick = 20.0 * mm;

  auto* activeSolid = new G4Box("HcalActive", slabXY / 2, slabXY / 2, activeThick / 2);
  auto* absSolid = new G4Box("HcalAbs", slabXY / 2, slabXY / 2, absorberThick / 2);
  auto* activeLV = new G4LogicalVolume(activeSolid, sciMat, "HcalActiveLV");
  auto* absLV = new G4LogicalVolume(absSolid, steelMat, "HcalAbsLV");

  for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer) {
    const double z_act = AHCALGeometry::Pos_Z(layer) * mm;
    const double z_abs = z_act + (activeThick + absorberThick) / 2.0;
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z_act), activeLV, "HcalActivePV", worldLV, false, layer, false);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z_abs), absLV, "HcalAbsPV", worldLV, false, layer, false);
  }

  if (m_cfg.enable_trigger_component && m_cfg.trigger_nplanes > 0) {
    auto* trgSolid = new G4Box("HcalTrigger", slabXY / 2, slabXY / 2, activeThick / 2);
    auto* trgLV = new G4LogicalVolume(trgSolid, sciMat, "HcalTriggerLV");
    for (int i = 0; i < m_cfg.trigger_nplanes; ++i) {
      const double z = AHCALGeometry::Pos_Z(AHCALGeometry::Layer_No - 1) * mm + (i + 1) * 40.0 * mm;
      new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z), trgLV, "HcalTriggerPV", worldLV, false, i, false);
    }
  }

  auto* sdManager = G4SDManager::GetSDMpointer();
  auto* sd = new HcalSD(m_hitMap);
  sdManager->AddNewDetector(sd);
  activeLV->SetSensitiveDetector(sd);

  return worldPV;
}

} // namespace AHCALRecoAlg::Sim
