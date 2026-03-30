#pragma once

#include "common/IAlg.hpp"
#include "simulation/geant4/Geant4Config.hpp"

#include <memory>
#include <string>
#include <vector>

class AHCALSimHit;
class SimData;
class InjectedParticle;

namespace AHCALRecoAlg::Sim {
class SimulationEngine;

struct Geant4AlgCfg {
  std::string in_injected_particles_key = "InjectedParticles";
  std::string out_simhits_key = "SimHits";
  std::string out_simdata_key = "SimData";
  bool require_input_particle = true;
  Geant4Config g4;
};

class Geant4Alg final : public IAlg {
public:
  Geant4Alg(RunContext& ctx, std::string name);
  ~Geant4Alg() override;

  void parse_cfg(const YAML::Node& n) override;
  void initialize() override;
  void execute(EventStore& evt) override;

private:
  Geant4AlgCfg m_cfg;
  std::unique_ptr<SimulationEngine> m_engine;

  void fill_simdata_(SimData& out) const;
  bool apply_external_primary_(const std::vector<InjectedParticle>& injected);
};

} // namespace AHCALRecoAlg::Sim
