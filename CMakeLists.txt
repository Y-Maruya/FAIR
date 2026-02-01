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
option(FAIR_WITH_SPDLOG "Use spdlog for logging" ON)
option(FAIR_WITH_YAML_CPP "Use yaml-cpp for configuration files" ON)

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
  # Prefer ROOT from the current environment (e.g. LCG view) over /usr
  if(DEFINED ENV{ROOTSYS} AND EXISTS "$ENV{ROOTSYS}")
    list(PREPEND CMAKE_PREFIX_PATH "$ENV{ROOTSYS}")
  endif()

  # Use CONFIG mode so we get imported targets like ROOT::Core, ROOT::ROOTDataFrame, etc.
  find_package(ROOT CONFIG QUIET COMPONENTS
    Core RIO Tree Hist Graf Gpad Graf3d MathCore
    ROOTDataFrame TreePlayer
  )

  if(ROOT_FOUND)
    message(STATUS "ROOT found: ${ROOT_VERSION} (ROOTSYS=$ENV{ROOTSYS})")
  else()
    message(FATAL_ERROR "ROOT not found or wrong ROOT selected. Please 'source' the desired ROOT (e.g. LCG view) and reconfigure.")
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
  ${CMAKE_CURRENT_SOURCE_DIR}/faser-common/EventFormats
  ${CMAKE_CURRENT_SOURCE_DIR}/faser-common/Exceptions
  ${CMAKE_CURRENT_SOURCE_DIR}/faser-common/Logging/include
)

# If ROOT exists, expose include dirs via interface as well (optional)
if(FAIR_WITH_ROOT AND ROOT_FOUND)
  target_link_libraries(fair_options INTERFACE 
  ROOT::Core
  ROOT::Graf
  ROOT::Graf3d
  ROOT::Hist
  ROOT::MathCore
  ROOT::RIO
  ROOT::Tree 
  ROOT::ROOTDataFrame
  )
endif()

# Use spdlog for logging (optional)
if(FAIR_WITH_SPDLOG)
  find_package(spdlog REQUIRED)

  target_link_libraries(fair_options
    INTERFACE
      spdlog::spdlog
  )

  target_compile_definitions(fair_options
    INTERFACE
      FAIR_USE_SPDLOG
  )
endif()

# Use YAML-CPP for configuration files (optional)
if(FAIR_WITH_YAML_CPP)
  find_package(yaml-cpp REQUIRED)

  target_link_libraries(fair_options
    INTERFACE
      yaml-cpp
  )

  target_compile_definitions(fair_options
    INTERFACE
      FAIR_USE_YAML_CPP
  )
endif()

# ---------------------------
# Subprojects
# ---------------------------

# faser-common (has its own CMakeLists.txt)
# This defines e.g. Logging, Exceptions, EventFormats targets (depending on its CMake)
# add_subdirectory(faser-common)
# Your modules (add CMakeLists.txt in each as you implement)
# You can comment out directories that don't have CMakeLists yet.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/common/CMakeLists.txt")
  add_subdirectory(common)
endif()
# --- add: zlib (needed by BinaryRawHitReader) ---
find_package(ZLIB REQUIRED)
target_link_libraries(fair_options INTERFACE ZLIB::ZLIB)
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
