add_executable(trackfit_test
  Test.cpp
)
#add_executable(muonkf_test
#  TestMuonKF.cpp
#)
if(TARGET fair_options)
  target_link_libraries(trackfit_test 
    PRIVATE 
      fair_options
      AdcToEnergyReadTTreeAlg
      TrackFitAlg
      MuonKFAlg
      RootRawHitReader
  )
#  target_link_libraries(muonkf_test PRIVATE fair_options)
endif()

target_include_directories(trackfit_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

#target_include_directories(muonkf_test
#  PRIVATE
#    ${CMAKE_SOURCE_DIR}
#)
