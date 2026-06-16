#pragma once
#include "common/AHCALGeometry.hpp"
#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"
#include <vector> 

struct LumiRecord {
    uint64_t since = 0;
    uint64_t until = 0;
    double instLumi = 0.0;
    double AvrageEventsPerBunchCrossing = 0.0;
    Long64_t runNumber = 0;
    Long64_t LumiBlockNumber = 0;
    Long64_t runlb = 0;
    int status = 0;
};

inline std::vector<FieldDesc> describe(const LumiRecord*) {
    return {
        field("since", &LumiRecord::since),
        field("until", &LumiRecord::until),
        field("instLumi", &LumiRecord::instLumi),
        field("AvrageEventsPerBunchCrossing", &LumiRecord::AvrageEventsPerBunchCrossing),
        field("runNumber", &LumiRecord::runNumber),
        field("LumiBlockNumber", &LumiRecord::LumiBlockNumber),
        field("runlb", &LumiRecord::runlb),
        field("status", &LumiRecord::status)
    };
}
AHCAL_REGISTER_IO_STRUCT(LumiRecord, "LumiRecord");