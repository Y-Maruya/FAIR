#pragma once

#include <string>

namespace AHCALRecoAlg::Sim {
struct Geant4Config {
  std::string detector_model = "hcal_option3";
  std::string generator_mode = "sim_reader"; // sim_reader | particle_gun | external_inject
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
  double hcal_step_time_limit_ns = 150.0;
  bool enable_pythia8_decayer = false;

  // CaloUnit-like mechanical structure (AHCAL-simulation-inspired)
  double passive_side_thickness_mm = 0.065;
  double passive_cover_thickness_mm = 0.065;
  double attach_thickness_mm = 0.0;
  bool double_sided_readout = false;
  double house_pitch_mm = 40.3;
  double sensitive_dig_out_x_mm = 5.5;
  double sensitive_dig_out_y_mm = 5.5;
  double sensitive_dig_out_z_mm = 1.1;
};
} // namespace AHCALRecoAlg::Sim
