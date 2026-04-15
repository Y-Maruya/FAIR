add_executable(trackfit_test
  Test.cpp
)
add_executable(calib_tester CalibITester.cpp)
add_executable(EventDisplay2D EventDisplay2D.cpp)
add_executable(EventDisplay2D_iso EventDisplay2D_iso.cpp)
add_executable(EventDisplay2D_multi EventDisplay2D_multi.cpp)
add_executable(fair_multi MultiInOut.cpp)
add_executable(fair_single MultiInOneOut.cpp)
add_executable(fair_run ExecRunbyRun.cpp)
add_executable(fair_calib ExecCalib.cpp)

if(TARGET fair_options)
  # Define library groups for reuse.
  add_library(core_reco_libs INTERFACE)
  target_link_libraries(core_reco_libs INTERFACE
    fair_options
    AdcToEnergyReadTTreeAlg
    reco_module_algs
    io_readers
    calib_algs
  )

  add_library(full_analysis_libs INTERFACE)
  target_link_libraries(full_analysis_libs INTERFACE
    core_reco_libs
    analysis_algs
    calib_io
  )

  function(collect_static_libs_from_target root_target out_var)
    set(_queue "${root_target}")
    set(_seen "")
    set(_static_libs "")

    while(_queue)
      list(POP_FRONT _queue _cur)

      list(FIND _seen "${_cur}" _seen_idx)
      if(NOT _seen_idx EQUAL -1)
        continue()
      endif()
      list(APPEND _seen "${_cur}")

      if(NOT TARGET ${_cur})
        continue()
      endif()

      get_target_property(_type ${_cur} TYPE)
      if(_type STREQUAL "STATIC_LIBRARY")
        list(APPEND _static_libs "${_cur}")
        continue()
      endif()

      if(_type STREQUAL "INTERFACE_LIBRARY")
        get_target_property(_deps ${_cur} INTERFACE_LINK_LIBRARIES)
      else()
        get_target_property(_deps ${_cur} LINK_LIBRARIES)
      endif()

      if(NOT _deps)
        continue()
      endif()

      foreach(_dep IN LISTS _deps)
        if(_dep MATCHES "^\\$<LINK_ONLY:([^>]+)>$")
          set(_dep "${CMAKE_MATCH_1}")
        elseif(_dep MATCHES "^\\$<")
          continue()
        endif()
        if(_dep MATCHES "^-")
          continue()
        endif()
        if(TARGET ${_dep})
          list(APPEND _queue "${_dep}")
        endif()
      endforeach()
    endwhile()

    list(REMOVE_DUPLICATES _static_libs)
    set(${out_var} "${_static_libs}" PARENT_SCOPE)
  endfunction()

  function(link_factory_executable target_name root_target)
    collect_static_libs_from_target(${root_target} _factory_static_libs)
    target_link_libraries(${target_name}
      PRIVATE
        fair_options
        -Wl,--whole-archive
        ${_factory_static_libs}
        -Wl,--no-whole-archive
        ${root_target}
    )
  endfunction()

  # Link trackfit_test with core reco libraries.
  link_factory_executable(trackfit_test core_reco_libs)
  target_link_libraries(trackfit_test
    PRIVATE
      calib_io
  )

  # Link calib_tester with fair_options and calib I/O
  target_link_libraries(calib_tester
    PRIVATE
      fair_options
      calib_io
  )

  # Link factory-based executables.
  link_factory_executable(fair_multi full_analysis_libs)
  link_factory_executable(fair_single full_analysis_libs)
  link_factory_executable(fair_run full_analysis_libs)
  link_factory_executable(fair_calib full_analysis_libs)
  # Link EventDisplay variants with core reco libraries (no whole-archive)
  target_link_libraries(EventDisplay2D
    PRIVATE
      fair_options
      AdcToEnergyReadTTreeAlg
      reco_module_algs
  )
  target_link_libraries(EventDisplay2D_iso 
    PRIVATE
      fair_options
      AdcToEnergyReadTTreeAlg
      reco_module_algs
  )
  target_link_libraries(EventDisplay2D_multi
    PRIVATE
      fair_options
      AdcToEnergyReadTTreeAlg
      reco_module_algs
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
target_include_directories(fair_run
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)
target_include_directories(fair_calib
  PRIVATE
    ${CMAKE_SOURCE_DIR}
)