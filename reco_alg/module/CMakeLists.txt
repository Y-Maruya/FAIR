# reco_alg/module/CMakeLists.txt
add_subdirectory(TrackFitAlg)
add_subdirectory(MuonKFAlg)
add_subdirectory(RmIsolatedHitAlg)
add_subdirectory(SummaryAlg)
add_subdirectory(RmNoiseHitAlg)

# Collect all reco module algorithms into a single interface library
add_library(reco_module_algs INTERFACE)
target_link_libraries(reco_module_algs INTERFACE
  TrackFitAlg MuonKFAlg RmIsolatedHitAlg SummaryAlg RmNoiseHitAlg
)