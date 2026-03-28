#pragma once

#include <string>

namespace AHCALRecoAlg::Sim {
struct Geant4Config {
  std::string detector_model = "hcal_option3";
  std::string generator_mode = "sim_reader"; // sim_reader | particle_gun
  std::string physics_list = "FTFP_BERT";

  // SimReader mode
  std::string input_root_file;
  std::string input_tree_name = "gFaser";
  bool fixed_vertex = false;

  int primary_pdg = 13;
  double primary_energy_GeV = 100.0;
  double primary_x_mm = 0.0;
  double primary_y_mm = 0.0;
  double primary_z_mm = -2000.0;
  double dir_x = 0.0;
  double dir_y = 0.0;
  double dir_z = 1.0;

  // detector details
  double birks_constant_mm_per_MeV = 0.126;
  bool enable_trigger_component = true;
  int trigger_nplanes = 1;
  bool enable_pythia8_decayer = false;
};
} // namespace AHCALRecoAlg::Sim
