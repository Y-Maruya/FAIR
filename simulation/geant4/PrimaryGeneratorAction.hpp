#pragma once

#include "simulation/geant4/Geant4Config.hpp"

#include "G4VUserPrimaryGeneratorAction.hh"

#include <memory>
#include <string>
#include <vector>

class G4ParticleGun;
class G4Event;
class TFile;
class TTree;

namespace AHCALRecoAlg::Sim {

struct PrimaryEventInfo {
  int injected_pdg = 0;
  double injected_energy_GeV = 0.0;
  double injected_px_GeV = 0.0;
  double injected_py_GeV = 0.0;
  double injected_pz_GeV = 0.0;
  double injected_x_mm = -999999;
  double injected_y_mm = -999999;
  double injected_z_mm = -999999;
  double interaction_x_mm = -999999;
  double interaction_y_mm = -999999;
  double interaction_z_mm = -999999;

  int ftagNulabel = -1; // 0:nueCC 1:numuCC 2:nutauCC 3:NC 4:other
  int secondary_pdg = 0;
  double secondary_energy_GeV = -999999;
  double secondary_px_GeV = -999999;
  double secondary_py_GeV = -999999;
  double secondary_pz_GeV = -999999;
};

class SimReaderPrimaryGenerator final : public G4VUserPrimaryGeneratorAction {
public:
  explicit SimReaderPrimaryGenerator(const Geant4Config& cfg);
  ~SimReaderPrimaryGenerator() override;

  void GeneratePrimaries(G4Event* event) override;
  static const PrimaryEventInfo& last_event_info();

private:
  void setup_root_reader_();
  bool generate_from_root_(G4Event* event);
  void generate_particle_gun_(G4Event* event);
  void update_interaction_info_();

  Geant4Config m_cfg;
  G4ParticleGun* m_gun = nullptr;

  std::unique_ptr<TFile> m_file;
  TTree* m_tree = nullptr;
  long long m_evt = 0;

  std::vector<int>* m_pdgc = nullptr;
  std::vector<double>* m_px = nullptr;
  std::vector<double>* m_py = nullptr;
  std::vector<double>* m_pz = nullptr;
  std::vector<double>* m_E = nullptr;
  std::vector<int>* m_status = nullptr;
  std::vector<int>* m_firstMother = nullptr;
  double m_vx = 0.0;
  double m_vy = 0.0;
  double m_vz = 0.0;

  static thread_local PrimaryEventInfo s_lastInfo;
};

} // namespace AHCALRecoAlg::Sim
