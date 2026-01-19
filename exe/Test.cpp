#include "TrackFitAlg.hpp"
#include <iostream>
using namespace AHCALRecoAlg;
int main() {
    std::vector<std::unique_ptr<IAlg>> algs;
    TrackFitAlg trackFitAlg("RecoAndRawHits", "FittedTrack");
    trackFitAlg.setThresholdXY(20.0);
    algs.emplace_back(std::make_unique<TrackFitAlg>(trackFitAlg));
    std::cout << "TrackFitAlg created successfully." << std::endl;
    return 0;
}