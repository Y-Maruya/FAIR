add_subdirectory(MuonEffAlg)
add_subdirectory(VetoEffAlg)
add_subdirectory(InputEffAlg)
add_subdirectory(NoiseHitAlg)
add_subdirectory(RateAnaAlg)
add_subdirectory(VetoAnaAlg)

# Collect all analysis algorithms into a single interface library
add_library(analysis_algs INTERFACE)
target_link_libraries(analysis_algs INTERFACE
  MuonEffAlg VetoEffAlg InputEffAlg NoiseHitAlg RateAnaAlg VetoAnaAlg
)