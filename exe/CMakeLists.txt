add_executable(trackfit_test
  Test.cpp
)

target_link_libraries(trackfit_test
  PRIVATE
    TrackFitAlg
)

target_include_directories(trackfit_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)
