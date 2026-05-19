#pragma once
#include "common/AHCALGeometry.hpp"
#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"
#include <vector> 

struct LHCData {
    uint64_t since = 0;
    uint64_t until = 0;
    int fillNumber = 0;
    bool stableBeams = false;
    std::string beamMode = "";
    std::string machineMode = "";
    int beamType1 = 0;
    int beamType2 = 0;
    int beamEnergyGeV = 0;
    double betaStar = 0.0;
    double crossingAngle = 0.0;
    std::string injectionScheme = "";
};

inline std::vector<FieldDesc> describe(const LHCData*) {
    return {
        field("since", &LHCData::since),
        field("until", &LHCData::until),
        field("fillNumber", &LHCData::fillNumber),
        field("stableBeams", &LHCData::stableBeams),
        field("beamMode", &LHCData::beamMode),
        field("machineMode", &LHCData::machineMode),
        field("beamType1", &LHCData::beamType1),
        field("beamType2", &LHCData::beamType2),
        field("beamEnergyGeV", &LHCData::beamEnergyGeV),
        field("betaStar", &LHCData::betaStar),
        field("crossingAngle", &LHCData::crossingAngle),
        field("injectionScheme", &LHCData::injectionScheme)
    };
}
AHCAL_REGISTER_IO_STRUCT(LHCData, "LHCData");