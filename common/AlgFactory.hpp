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
#include "common/AlgRegistry.hpp"
// ---- calibration structs ----
#include "calibration/module/pedestal/PedestalAlg.hpp"
// ---- your alg headers ----
#include "adc_to_energy/AdcToEnergyReadTTreeAlg.hpp"
#include "reco_alg/module/TrackFitAlg/TrackFitAlg.hpp"
#include "reco_alg/module/MuonKFAlg/MuonKFAlg.hpp"
// -------------------------------------
// reader algs
#include "IO/reader/ReaderRegistry.hpp"
// writer algs
#include "IO/writer/RootWriterAlg.hpp"
// ---------------------------------------------


inline WriterRegistry parse_writer_registry(const YAML::Node& n) {
  WriterRegistry reg;
  const auto out = require_node(n, "outputlist");
  if (!out.IsSequence()) throw std::runtime_error("RootWriterAlg.cfg.outputlist must be a sequence");

  for (const auto& x : out) {
    const auto type = x.as<std::string>();

    if (auto fn = IOTypeRegistry::instance().get_writer(type)) {
      fn(reg);
      continue;
    }

    LOG_ERROR("AlgFactory::parse_writer_registry: Unknown writer type '{}'", type);
    LOG_ERROR("IF YOU HAVE ADDED A NEW STRUCT, MAKE SURE TO REGISTER IT (e.g. in its EDM header).");
    throw std::runtime_error("Unknown writer type: " + type);
  }
  return reg;
}

inline ReaderRegistry parse_reader_registry(const YAML::Node& n) {
  ReaderRegistry reg;
  const auto out = require_node(n, "inputlist");
  if (!out.IsSequence()) throw std::runtime_error("RootInput.cfg.inputlist must be a sequence");

  for (const auto& x : out) {
    const auto type_array = x.as<std::vector<std::string>>();
    if (type_array.size() != 2) {
      LOG_ERROR("AlgFactory::parse_reader_registry: Each inputlist entry must be a pair of [type, key]");
      throw std::runtime_error("Invalid inputlist entry size");
    }
    const auto& type = type_array[0];
    const auto& key  = type_array[1];

    if (auto fn = IOTypeRegistry::instance().get_reader(type)) {
      fn(reg, key);
      continue;
    }

    LOG_ERROR("AlgFactory::parse_reader_registry: Unknown reader type '{}'", type);
    LOG_ERROR("IF YOU HAVE ADDED A NEW STRUCT, MAKE SURE TO REGISTER IT (e.g. in its EDM header).");
    throw std::runtime_error("Unknown reader type: " + type);
  }
  return reg;
}

inline void readandput(const YAML::Node& n, EventStore& store, ReaderRegistry& rr, RootInput& in) {
  const auto out = require_node(n, "inputlist");
  if (!out.IsSequence()) throw std::runtime_error("RootInput.cfg.inputlist must be a sequence");

  for (const auto& x : out) {
    const auto type_array = x.as<std::vector<std::string>>();
    if (type_array.size() != 2) {
      LOG_ERROR("AlgFactory::readandput: Each inputlist entry must be a pair of [type, key]");
      throw std::runtime_error("Invalid inputlist entry size");
    }

    const auto& type = type_array[0];
    const auto& key  = type_array[1];

    if (auto fn = IOTypeRegistry::instance().get_readput(type)) {
      fn(store, rr, in, type, key);
      continue;
    }

    LOG_ERROR("AlgFactory::readandput: Unknown reader type '{}'", type);
    LOG_ERROR("IF YOU HAVE ADDED A NEW STRUCT, MAKE SURE TO REGISTER IT (e.g. in its EDM header).");
    throw std::runtime_error("Unknown reader type: " + type);
  }
}

inline std::unique_ptr<IAlg> make_alg(RunContext& ctx, const YAML::Node& alg_node) {
  const std::string type = require_string(alg_node, "type");
  const YAML::Node cfg = alg_node["cfg"] ? alg_node["cfg"] : YAML::Node(YAML::NodeType::Map);

    if (type == "RootWriterAlg") {
      LOG_INFO("Creating RootWriterAlg");
      return std::make_unique<RootWriterAlg>(
          ctx,
          "RootWriterAlg",
          ctx.config.output,
          parse_writer_registry(cfg));
    }

    try {
        auto alg = AlgRegistry::instance().create(type, ctx, cfg);
        LOG_INFO("{}: initialized", type);
        return alg;
    } catch (const std::exception& e) {
        LOG_ERROR("AlgFactory::make_alg: {}", e.what());
        throw;
    }
}

inline std::vector<std::unique_ptr<IAlg>> build_pipeline(RunContext& ctx, const YAML::Node& root) {
  const auto& algs = require_node(root, "algs");
  if (!algs.IsSequence()) throw std::runtime_error("algs must be a YAML sequence");

  std::vector<std::unique_ptr<IAlg>> v;
  v.reserve(algs.size());
  for (const auto& a : algs) v.emplace_back(make_alg(ctx, a));
  return v;
}


