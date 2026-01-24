#pragma once
#include "common/edm/RecoHit.hpp"
#include <vector>
#include "IO/Descriptor.hpp"


struct SimpleFittedTrack {
    double init_pos_x; // (x0, y0) at z=0
    double init_pos_y;
    double direction_x; // (dx/dz, dy/dz)
    double direction_y;
    double chi2_x; // (chi2_x, chi2_y)
    double chi2_y;
    int ndf = 0;
    std::vector<int> inTrackHitsIndices;
    std::vector<int> outTrackHitsIndices;
    int nTotalHits = 0;
    bool valid = false;
    std::vector<AHCALRecoHit> inTrackHits; // for FAIR internal use
    std::vector<AHCALRecoHit> outTrackHits; // for FAIR internal use
};

inline std::vector<FieldDesc> describe(const SimpleFittedTrack*) {
    return {
        field("init_pos_x", &SimpleFittedTrack::init_pos_x),
        field("init_pos_y", &SimpleFittedTrack::init_pos_y),
        field("direction_x", &SimpleFittedTrack::direction_x),
        field("direction_y", &SimpleFittedTrack::direction_y),
        field("chi2_x", &SimpleFittedTrack::chi2_x),
        field("chi2_y", &SimpleFittedTrack::chi2_y),
        field("ndf", &SimpleFittedTrack::ndf),
        field("inTrackHitsIndices", &SimpleFittedTrack::inTrackHitsIndices),
        field("outTrackHitsIndices", &SimpleFittedTrack::outTrackHitsIndices),
        field("nTotalHits", &SimpleFittedTrack::nTotalHits),
        field("valid", &SimpleFittedTrack::valid),
    };
}