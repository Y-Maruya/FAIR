#include "simulation/geant4/DetectorConstruction.hpp"

#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "simulation/geant4/HcalSD.hpp"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <algorithm>

namespace AHCALRecoAlg::Sim {

DetectorConstruction::DetectorConstruction(const Geant4Config& cfg, std::unordered_map<int, HitRecord>* hitMap)
    : m_cfg(cfg), m_hitMap(hitMap) {}

G4VPhysicalVolume* DetectorConstruction::Construct() {
  auto* nist = G4NistManager::Instance();
  auto* worldMat = nist->FindOrBuildMaterial("G4_AIR");
  auto* sciMat = nist->BuildMaterialWithNewDensity("FAIR_POLY", "G4_POLYSTYRENE", 1.032 * g / cm3);
  sciMat->GetIonisation()->SetBirksConstant(m_cfg.birks_constant_mm_per_MeV * mm / MeV);
  auto* steelMat = nist->FindOrBuildMaterial("G4_Fe");
  auto* passiveMat = nist->FindOrBuildMaterial("G4_MYLAR");
  auto* attachMat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  auto* galactic = nist->FindOrBuildMaterial("G4_Galactic");

  const double worldXY = 2600.0 * mm;
  const double worldZ = 6000.0 * mm;
  auto* worldSolid = new G4Box("World", worldXY / 2, worldXY / 2, worldZ / 2);
  auto* worldLV = new G4LogicalVolume(worldSolid, worldMat, "WorldLV");
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "WorldPV", nullptr, false, 0, false);

  // AHCAL-simulation CaloUnitVolume-inspired parameters.
  const double sensitiveXY = AHCALGeometry::xy_size * 18.0 * mm;
  const double houseXY = std::max(0.0, m_cfg.house_pitch_mm) * 18.0 * mm;
  const double sensitiveZ = 3.0 * mm;
  const double absorberThick = 20.0 * mm;

  const double passiveSide = std::max(0.0, m_cfg.passive_side_thickness_mm) * mm;
  const double passiveCover = std::max(0.0, m_cfg.passive_cover_thickness_mm) * mm;
  const double attachThick = std::max(0.0, m_cfg.attach_thickness_mm) * mm;
  const bool doubleSided = m_cfg.double_sided_readout;

  const double passiveX = sensitiveXY / 2.0 + passiveSide;
  const double passiveY = sensitiveXY / 2.0 + passiveSide;
  const double passiveZ = sensitiveZ / 2.0 + passiveCover;
  const double housingX = houseXY / 2.0;
  const double housingY = houseXY / 2.0;
  const double housingZ = passiveZ + (doubleSided ? 1.0 : 0.5) * attachThick;

  const double sensitiveDigOutX = std::max(0.0, m_cfg.sensitive_dig_out_x_mm) * mm;
  const double sensitiveDigOutY = std::max(0.0, m_cfg.sensitive_dig_out_y_mm) * mm;
  const double sensitiveDigOutZ = std::max(0.0, m_cfg.sensitive_dig_out_z_mm) * mm;

  auto* unitHousingSolid = new G4Box("HcalUnitHousing", housingX, housingY, housingZ);
  auto* unitHousingLV = new G4LogicalVolume(unitHousingSolid, galactic, "HcalUnitHousingLV");

  auto* passiveSolid = new G4Box("HcalUnitPassive", passiveX, passiveY, passiveZ);
  auto* passiveLV = new G4LogicalVolume(passiveSolid, passiveMat, "HcalUnitPassiveLV");
  const double passiveZOffset = (!doubleSided) ? attachThick / 2.0 : 0.0;
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, passiveZOffset), passiveLV, "HcalUnitPassivePV", unitHousingLV, false, 0, false);

  auto* sensitiveSolid = new G4Box("HcalUnitSensitive", sensitiveXY / 2.0, sensitiveXY / 2.0, sensitiveZ / 2.0);
  auto* sensitiveLV = new G4LogicalVolume(sensitiveSolid, sciMat, "HcalUnitSensitiveLV");
  new G4PVPlacement(nullptr, G4ThreeVector(), sensitiveLV, "HcalUnitSensitivePV", passiveLV, false, 0, false);

  if (sensitiveDigOutX > 0.0 && sensitiveDigOutY > 0.0 && sensitiveDigOutZ > 0.0) {
    auto* digOutSolid = new G4Box("HcalUnitSensitiveDigOut", sensitiveDigOutX / 2.0, sensitiveDigOutY / 2.0, sensitiveDigOutZ / 2.0);
    auto* digOutLV = new G4LogicalVolume(digOutSolid, worldMat, "HcalUnitSensitiveDigOutLV");
    new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, -sensitiveZ / 2.0 + sensitiveDigOutZ / 2.0), digOutLV,
                      "HcalUnitSensitiveDigOutPV", sensitiveLV, false, 0, false);
  }

  if (sensitiveDigOutX > 0.0 && sensitiveDigOutY > 0.0 && passiveCover > 0.0) {
    auto* digOutEsrSolid = new G4Box("HcalUnitSensitiveDigOutESR", sensitiveDigOutX / 2.0, sensitiveDigOutY / 2.0, passiveCover / 2.0);
    auto* digOutEsrLV = new G4LogicalVolume(digOutEsrSolid, worldMat, "HcalUnitSensitiveDigOutESRLV");
    new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, -sensitiveZ / 2.0 - passiveCover / 2.0), digOutEsrLV,
                      "HcalUnitSensitiveDigOutESRPV", passiveLV, false, 0, false);
  }

  if (attachThick > 0.0) {
    auto* attachSolid = new G4Box("HcalUnitAttach", passiveX, passiveY, attachThick / 2.0);
    auto* attachLV = new G4LogicalVolume(attachSolid, attachMat, "HcalUnitAttachLV");
    if (doubleSided) {
      new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -(passiveZ + attachThick / 2.0)), attachLV,
                        "HcalUnitAttachPV", unitHousingLV, false, 0, false);
      new G4PVPlacement(nullptr, G4ThreeVector(0, 0, +(passiveZ + attachThick / 2.0)), attachLV,
                        "HcalUnitAttachPV", unitHousingLV, false, 1, false);
    } else {
      new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -(passiveZ)), attachLV,
                        "HcalUnitAttachPV", unitHousingLV, false, 0, false);
    }
  }

  auto* absSolid = new G4Box("HcalAbs", houseXY / 2.0, houseXY / 2.0, absorberThick / 2);
  auto* absLV = new G4LogicalVolume(absSolid, steelMat, "HcalAbsLV");

  for (int layer = 0; layer < AHCALGeometry::Layer_No; ++layer) {
    const double z_act = AHCALGeometry::Pos_Z(layer) * mm;
    const double z_abs = z_act + (sensitiveZ + absorberThick) / 2.0;
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z_act), unitHousingLV, "HcalUnitHousingPV", worldLV, false, layer, false);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z_abs), absLV, "HcalAbsPV", worldLV, false, layer, false);
  }

  if (m_cfg.enable_trigger_component && m_cfg.trigger_nplanes > 0) {
    LOG_WARN("Trigger planes are currently passive geometry in Geant4 output (no AHCAL channel mapping).");
    auto* trgSolid = new G4Box("HcalTrigger", sensitiveXY / 2, sensitiveXY / 2, sensitiveZ / 2);
    auto* trgLV = new G4LogicalVolume(trgSolid, sciMat, "HcalTriggerLV");
    for (int i = 0; i < m_cfg.trigger_nplanes; ++i) {
      const double z = AHCALGeometry::Pos_Z(AHCALGeometry::Layer_No - 1) * mm + (i + 1) * 40.0 * mm;
      new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z), trgLV, "HcalTriggerPV", worldLV, false, i, false);
    }
  }

  auto* sdManager = G4SDManager::GetSDMpointer();
  auto* sd = new HcalSD(m_hitMap, m_cfg.hcal_step_time_limit_ns);
  sdManager->AddNewDetector(sd);
  sensitiveLV->SetSensitiveDetector(sd);

  return worldPV;
}

} // namespace AHCALRecoAlg::Sim
