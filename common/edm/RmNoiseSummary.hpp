#pragma once
#include "common/AHCALGeometry.hpp"
#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"
#include <vector> 
struct RmNoiseSummary {
    int total_hits = 0;
    int noise_hits = 0;
    double noise_fraction = 0.0;
};

inline std::vector<FieldDesc> describe(const RmNoiseSummary*) {
    return {
        field("total_hits", &RmNoiseSummary::total_hits),
        field("noise_hits", &RmNoiseSummary::noise_hits),
        field("noise_fraction", &RmNoiseSummary::noise_fraction)
    };
}
AHCAL_REGISTER_IO_STRUCT(RmNoiseSummary, "RmNoiseSummary");