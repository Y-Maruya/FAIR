#pragma once
#include <yaml-cpp/yaml.h>
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>

#include "common/config/YAMLUtil.hpp"
#include "common/RunContext.hpp"
#include "common/IAlg.hpp"
#include "common/config/ParseRunConfig.hpp"

// ---- your alg headers ----
#include "adc_to_energy/AdcToEnergyReadTTreeAlg.hpp"
#include "reco_alg/module/TrackFitAlg/TrackFitAlg.hpp"
#include "reco_alg/module/MuonKFAlg/MuonKFAlg.hpp"
// -------------------------------------
// writer algs
#include "IO/writer/RootWriterAlg.hpp"
// ---------------------------------------------


inline WriterRegistry parse_writer_registry(const YAML::Node& n) {
  WriterRegistry reg;
  const auto out = require_node(n, "outputlist");
  if (!out.IsSequence()) throw std::runtime_error("RootWriterAlg.cfg.outputlist must be a sequence");

  for (const auto& x : out) {
    const auto type = x.as<std::string>();
    if (type == "vector<AHCALRawHit>") {
      reg.register_vector_struct<AHCALRawHit>();
    } else if (type == "vector<AHCALRecoHit>") {
      reg.register_vector_struct<AHCALRecoHit>();
    } else if (type == "SimpleFittedTrack") {
      reg.register_struct<SimpleFittedTrack>();
    } else if (type == "AHCALRawHit") {
      reg.register_struct<AHCALRawHit>();
    } else if (type == "AHCALRecoHit") {
      reg.register_struct<AHCALRecoHit>();
    } else if (type == "vector<Track>") {
      reg.register_vector_struct<Track>();
    } else if (type == "Track") {
      reg.register_struct<Track>();
    } else {
        LOG_ERROR("AlgFactory::parse_writer_registry: Unknown writer type '{}'", type);
        LOG_ERROR("IF YOU HAVE ADDED A NEW STRUCT, MAKE SURE TO REGISTER IT HERE.");
        throw std::runtime_error("Unknown writer type: " + type);
    }
  }
  return reg;
}
inline std::unique_ptr<IAlg> make_alg(RunContext& ctx, const YAML::Node& alg_node) {
  const std::string type = require_string(alg_node, "type");
  const YAML::Node cfg = alg_node["cfg"] ? alg_node["cfg"] : YAML::Node(YAML::NodeType::Map);

  if (type == "AdcToEnergyReadTTreeAlg") {   
    auto alg = std::make_unique<AHCALRecoAlg::AdcToEnergyReadTTreeAlg>(ctx, "AdcToEnergyReadTTreeAlg");
    alg->parse_cfg(cfg);
    alg->initialize_mip();
    alg->initialize_ped();
    alg->initialize_dac();
    return alg;
  }

  if (type == "TrackFitAlg") {
    auto alg = std::make_unique<AHCALRecoAlg::TrackFitAlg>(ctx, "TrackFitAlg");
    alg->parse_cfg(cfg);
    return alg;
  }

  // Add other alg types here
  if (type == "MuonKFAlg") {
    auto alg = std::make_unique<AHCALRecoAlg::MuonKFAlg>(ctx, "MuonKFAlg");
    alg->parse_cfg(cfg);
    return alg;
  }


  // Writer algs
  if (type == "RootWriterAlg") {
    return std::make_unique<RootWriterAlg>(
        ctx,
        "RootWriterAlg",
        ctx.config.output,
        parse_writer_registry(cfg));
  }

  throw std::runtime_error("Unknown alg type: " + type);
}

inline std::vector<std::unique_ptr<IAlg>> build_pipeline(RunContext& ctx, const YAML::Node& root) {
  const auto& algs = require_node(root, "algs");
  if (!algs.IsSequence()) throw std::runtime_error("algs must be a YAML sequence");

  std::vector<std::unique_ptr<IAlg>> v;
  v.reserve(algs.size());
  for (const auto& a : algs) v.emplace_back(make_alg(ctx, a));
  return v;
}


