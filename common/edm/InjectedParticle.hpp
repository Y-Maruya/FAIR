#pragma once

#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"

#include <vector>

class InjectedParticle {
public:
  int pdgId = 13;
  double energy_GeV = 0.0;
  double time_ns = 0.0;
  double x_mm = 0.0;
  double y_mm = 0.0;
  double z_mm = -2000.0;
  double px_GeV = 0.0;
  double py_GeV = 0.0;
  double pz_GeV = 1.0;
};

inline std::vector<FieldDesc> describe(const InjectedParticle*) {
  return {
      field("pdgId", &InjectedParticle::pdgId),
      field("energy_GeV", &InjectedParticle::energy_GeV),
      field("time_ns", &InjectedParticle::time_ns),
      field("x_mm", &InjectedParticle::x_mm),
      field("y_mm", &InjectedParticle::y_mm),
      field("z_mm", &InjectedParticle::z_mm),
      field("px_GeV", &InjectedParticle::px_GeV),
      field("py_GeV", &InjectedParticle::py_GeV),
      field("pz_GeV", &InjectedParticle::pz_GeV),
  };
}

inline std::vector<FieldDescVector> describe_vector(const InjectedParticle*) {
  return {
      field_vector("v.pdgId", &InjectedParticle::pdgId),
      field_vector("v.energy_GeV", &InjectedParticle::energy_GeV),
      field_vector("v.time_ns", &InjectedParticle::time_ns),
      field_vector("v.x_mm", &InjectedParticle::x_mm),
      field_vector("v.y_mm", &InjectedParticle::y_mm),
      field_vector("v.z_mm", &InjectedParticle::z_mm),
      field_vector("v.px_GeV", &InjectedParticle::px_GeV),
      field_vector("v.py_GeV", &InjectedParticle::py_GeV),
      field_vector("v.pz_GeV", &InjectedParticle::pz_GeV),
  };
}

AHCAL_REGISTER_IO_STRUCT(InjectedParticle, "InjectedParticle");
AHCAL_REGISTER_IO_STRUCT_VECTOR(InjectedParticle, "vector<InjectedParticle>");
