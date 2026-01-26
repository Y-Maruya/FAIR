#pragma once
#include <vector>
#include "common/edm/RecoHit.hpp"
struct Track {
    // state at the last updated z in the fit loop: (x,y,tx,ty)
    double x = 0.0;
    double y = 0.0;
    double tx = 0.0;
    double ty = 0.0;
    double z = 0.0;
    double chi2 = 0.0;
    int ndof = 0;
    int consecutive_skips = 0;
    bool valid = false;
    int nInTrackHits = 0;
    int nOutTrackHits = 0;
    std::vector<int> inTrackHitsIndices;
    std::vector<int> outTrackHitsIndices;
    std::vector<AHCALRecoHit> inTrackHits; // for FAIR internal use
    std::vector<AHCALRecoHit> outTrackHits; // for FAIR internal use
};

inline std::vector<FieldDesc> describe(const Track*) {
    return {
        field("x", &Track::x),
        field("y", &Track::y),
        field("tx", &Track::tx),
        field("ty", &Track::ty),
        field("z", &Track::z),
        field("chi2", &Track::chi2),
        field("ndof", &Track::ndof),
        field("consecutive_skips", &Track::consecutive_skips),
        field("nInTrackHits", &Track::nInTrackHits),
        field("nOutTrackHits", &Track::nOutTrackHits),
        field("inTrackHitsIndices", &Track::inTrackHitsIndices),
        field("outTrackHitsIndices", &Track::outTrackHitsIndices),
        field("valid", &Track::valid)
    };
}
inline std::vector<FieldDescVector> describe_vector(const Track*) {
    return {
        field_vector("v.x", &Track::x),
        field_vector("v.y", &Track::y),
        field_vector("v.tx", &Track::tx),
        field_vector("v.ty", &Track::ty),
        field_vector("v.z", &Track::z),
        field_vector("v.chi2", &Track::chi2),
        field_vector("v.ndof", &Track::ndof),
        field_vector("v.consecutive_skips", &Track::consecutive_skips),
        field_vector("v.nInTrackHits", &Track::nInTrackHits),
        field_vector("v.nOutTrackHits", &Track::nOutTrackHits),
        field_vector("v.inTrackHitsIndices", &Track::inTrackHitsIndices),
        field_vector("v.outTrackHitsIndices", &Track::outTrackHitsIndices),
        field_vector("v.valid", &Track::valid)
    };
}