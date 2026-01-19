# CMakeLists.txt (top)
cmake_minimum_required(VERSION 3.16)

project(FAIR
  VERSION 0.1.0
  DESCRIPTION "FAIR offline framework (AHCAL etc.)"
  LANGUAGES CXX
)

# ---------------------------
# Global build settings
# ---------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default build type (if single-config generator)
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

# Output dirs (nice on lxplus)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

option(FAIR_BUILD_TESTS "Build tests" OFF)
option(FAIR_ENABLE_WARNINGS "Enable extra compiler warnings" ON)

if(FAIR_ENABLE_WARNINGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
  endif()
endif()

# Helpful for shared libs / plugins
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ---------------------------
# Dependencies
# ---------------------------
# ROOT is optional at top-level (enable when needed in subdirs)
option(FAIR_WITH_ROOT "Enable ROOT-dependent components" ON)
if(FAIR_WITH_ROOT)
  find_package(ROOT QUIET)
  if(ROOT_FOUND)
    message(STATUS "ROOT found: ${ROOT_VERSION}")
    # ROOT::Core etc. will be available to subdirectories
  else()
    message(STATUS "ROOT not found (FAIR_WITH_ROOT=ON). ROOT-dependent targets may be skipped.")
  endif()
endif()

# ---------------------------
# Common include / interface target
# ---------------------------
add_library(fair_options INTERFACE)
target_compile_features(fair_options INTERFACE cxx_std_17)

# Common include dirs (your headers)
target_include_directories(fair_options INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/common
  ${CMAKE_CURRENT_SOURCE_DIR}/IO
  ${CMAKE_CURRENT_SOURCE_DIR}/calibration
  ${CMAKE_CURRENT_SOURCE_DIR}/adc_to_energy
  ${CMAKE_CURRENT_SOURCE_DIR}/reco_alg
)

# If ROOT exists, expose include dirs via interface as well (optional)
if(FAIR_WITH_ROOT AND ROOT_FOUND)
  target_link_libraries(fair_options INTERFACE ROOT::Core)
endif()

# ---------------------------
# Subprojects
# ---------------------------

# faser-common (has its own CMakeLists.txt)
# This defines e.g. Logging, Exceptions, EventFormats targets (depending on its CMake)
add_subdirectory(faser-common)

# Your modules (add CMakeLists.txt in each as you implement)
# You can comment out directories that don't have CMakeLists yet.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/common/CMakeLists.txt")
  add_subdirectory(common)
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/reco_alg/CMakeLists.txt")
  add_subdirectory(reco_alg)
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/IO/CMakeLists.txt")
  add_subdirectory(IO)
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/calibration/CMakeLists.txt")
  add_subdirectory(calibration)
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/adc_to_energy/CMakeLists.txt")
  add_subdirectory(adc_to_energy)
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/simulation/CMakeLists.txt")
  add_subdirectory(simulation)
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/analysis/CMakeLists.txt")
  add_subdirectory(analysis)
endif()

# Executables (driver programs)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/exe/CMakeLists.txt")
  add_subdirectory(exe)
endif()

# QA/QC (optional)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/qaqc/CMakeLists.txt")
  add_subdirectory(qaqc)
endif()

# ---------------------------
# Tests (optional)
# ---------------------------
if(FAIR_BUILD_TESTS)
  enable_testing()
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
    add_subdirectory(tests)
  endif()
endif()

## install rules (optional)
install(TARGETS fair_options
  EXPORT fairOptionsTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

message(STATUS "FAIR configured. Build dir: ${CMAKE_BINARY_DIR}")
