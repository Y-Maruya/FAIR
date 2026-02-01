add_executable(trackfit_test
  Test.cpp
)
add_executable(eventdisplay2D EventDisplay2D.cpp)
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
      -Wl,--no-whole-archive
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
      -Wl,--no-whole-archive
  )
  target_link_libraries(eventdisplay2D 
    PRIVATE
      fair_options
      AdcToEnergyReadTTreeAlg
      MuonKFAlg
  )
endif()

target_include_directories(trackfit_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

target_include_directories(eventdisplay2D
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
