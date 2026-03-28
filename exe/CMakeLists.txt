add_executable(trackfit_test
  Test.cpp
)
add_executable(calib_tester CalibITester.cpp)
add_executable(EventDisplay2D EventDisplay2D.cpp)
add_executable(EventDisplay2D_iso EventDisplay2D_iso.cpp)
add_executable(EventDisplay2D_multi EventDisplay2D_multi.cpp)
add_executable(fair_multi MultiInOut.cpp)
add_executable(fair_single MultiInOneOut.cpp)
if(TARGET fair_options)
  target_link_libraries(trackfit_test 
    PRIVATE 
      fair_options
      -Wl,--whole-archive
      AdcToEnergyReadTTreeAlg
      TrackFitAlg
      MuonKFAlg
      RootRawHitReader
      BinaryRawHitReader
      PedestalAlg
      RmIsolatedHitAlg
      -Wl,--no-whole-archive
  )
  target_link_libraries(calib_tester
    PRIVATE
      fair_options
      PedestalReader
  )
  target_link_libraries(fair_multi 
    PRIVATE
      fair_options
      -Wl,--whole-archive
      AdcToEnergyReadTTreeAlg
      TrackFitAlg
      MuonKFAlg
      RootRawHitReader
      BinaryRawHitReader
      PedestalAlg
      RmIsolatedHitAlg
      MIPAlg
      MuonEffAlg
      VetoEffAlg
      InputEffAlg
      NoiseHitAlg
      SimHitReader
      RateAnaAlg
      -Wl,--no-whole-archive
  )
  target_link_libraries(fair_single 
    PRIVATE
      fair_options
      -Wl,--whole-archive
      AdcToEnergyReadTTreeAlg
      TrackFitAlg
      MuonKFAlg
      RootRawHitReader
      BinaryRawHitReader
      PedestalAlg
      RmIsolatedHitAlg
      MIPAlg
      MuonEffAlg
      VetoEffAlg
      InputEffAlg
      NoiseHitAlg
      SimHitReader
      RateAnaAlg
      -Wl,--no-whole-archive
  )
  if(FAIR_WITH_GEANT4 AND Geant4_FOUND)
    target_link_libraries(fair_multi PRIVATE Geant4SimReader)
    target_link_libraries(fair_single PRIVATE Geant4SimReader)
    target_compile_definitions(fair_multi PRIVATE FAIR_HAS_GEANT4=1)
    target_compile_definitions(fair_single PRIVATE FAIR_HAS_GEANT4=1)
  endif()
  target_link_libraries(EventDisplay2D
    PRIVATE
      fair_options
      AdcToEnergyReadTTreeAlg
      MuonKFAlg
      RmIsolatedHitAlg
  )
  target_link_libraries(EventDisplay2D_iso 
    PRIVATE
      fair_options
      AdcToEnergyReadTTreeAlg
      MuonKFAlg
      RmIsolatedHitAlg
  )
  target_link_libraries(EventDisplay2D_multi
    PRIVATE
      fair_options
      AdcToEnergyReadTTreeAlg
      MuonKFAlg
      RmIsolatedHitAlg
  )
endif()

target_include_directories(trackfit_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

target_include_directories(calib_tester
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

target_include_directories(EventDisplay2D
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)
target_include_directories(EventDisplay2D_iso
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)
target_include_directories(EventDisplay2D_multi
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

target_include_directories(fair_multi
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)
target_include_directories(fair_single
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)
