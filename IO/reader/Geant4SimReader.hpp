#pragma once

#include "common/edm/SimData.hpp"
#include "common/edm/SimHit.hpp"
#include "simulation/geant4/Geant4Config.hpp"
#include "simulation/geant4/SimulationEngine.hpp"

#include <memory>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

class Geant4SimReader {
public:
  struct Config {
    std::string out_simhits_key = "SimHits";
    std::string out_simdata_key = "SimData";
    int max_events = -1;
    int seed = 12345;
    AHCALRecoAlg::Sim::Geant4Config g4;
  };

  explicit Geant4SimReader(const YAML::Node& cfg);

  bool next(std::vector<AHCALSimHit>& out_hits, SimData& out_simdata);

  long long generated_events() const { return m_event; }
  const Config& config() const { return m_cfg; }

private:
  Config m_cfg;
  long long m_event = 0;

  std::unique_ptr<AHCALRecoAlg::Sim::SimulationEngine> m_engine;

  void fill_simdata_(SimData& out) const;
  static int pdg_from_name_(const std::string& name);
};
