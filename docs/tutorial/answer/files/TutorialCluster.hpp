#pragma once

#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"
#include "common/edm/RecoHit.hpp"
#include <vector>

namespace AHCALRecoAlg {

struct TutorialCluster {
    std::vector<AHCALRecoHit> hits;
    std::vector<int> index;
    int nHitInCluster = 0;
    double totalEdep = 0.0;
    double avgEdep = 0.0;
    void updateSummary();
};

inline std::vector<FieldDesc> describe(const TutorialCluster*) {
    return {
        field("nHitInCluster", &TutorialCluster::nHitInCluster),
        field("totalEdep", &TutorialCluster::totalEdep),
        field("avgEdep", &TutorialCluster::avgEdep),
    };
}

inline std::vector<FieldDescVector> describe_vector(const TutorialCluster*) {
    return {
        field_vector("v.nHitInCluster", &TutorialCluster::nHitInCluster),
        field_vector("v.totalEdep", &TutorialCluster::totalEdep),
        field_vector("v.avgEdep", &TutorialCluster::avgEdep),
    };
}

AHCAL_REGISTER_IO_STRUCT(TutorialCluster, "TutorialCluster");
AHCAL_REGISTER_IO_STRUCT_VECTOR(TutorialCluster, "vector<TutorialCluster>");

} // namespace AHCALRecoAlg
