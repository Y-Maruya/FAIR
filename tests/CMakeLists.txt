add_executable(test_expand_config
  test_expand_config.cpp
)

target_link_libraries(test_expand_config
  PRIVATE
    fair_options
)

target_include_directories(test_expand_config
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

add_test(NAME test_expand_config COMMAND test_expand_config)

add_executable(test_resolved_config
  test_resolved_config.cpp
)

target_link_libraries(test_resolved_config
  PRIVATE
    fair_options
)

target_include_directories(test_resolved_config
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

add_test(NAME test_resolved_config COMMAND test_resolved_config)

add_executable(test_runinfo_snapshot
  test_runinfo_snapshot.cpp
)

target_link_libraries(test_runinfo_snapshot
  PRIVATE
    fair_options
)

target_include_directories(test_runinfo_snapshot
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

add_test(NAME test_runinfo_snapshot COMMAND test_runinfo_snapshot)

find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_test(
  NAME test_build_runinfo_snapshot
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_build_runinfo_snapshot.py
)

add_executable(test_mip_crossing_selection
  test_mip_crossing_selection.cpp
)

target_link_libraries(test_mip_crossing_selection
  PRIVATE
    MIPAlg
    PedestalReader
    fair_options
)

target_include_directories(test_mip_crossing_selection
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)

add_test(NAME test_mip_crossing_selection COMMAND test_mip_crossing_selection)
