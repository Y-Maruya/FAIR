# reco_alg/module/LHCDataReaderAlg/CMakeLists.txt
cmake_minimum_required(VERSION 3.16)


set(TDAQ_EXTERNAL
    "/cvmfs/atlas.cern.ch/repo/sw/tdaq/offline/tdaq-common/tdaq-common-11-02-01/installed/external/x86_64-el9-gcc13-opt"
    CACHE PATH "TDAQ external installation containing COOL/CORAL")
set(BOOST_ROOT
    "/cvmfs/atlas.cern.ch/repo/sw/software/24.0/sw/lcg/releases/LCG_104d_ATLAS_7/Boost/1.82.0/x86_64-el9-gcc13-opt"
    CACHE PATH "LCG Boost installation used by COOL/CORAL")
set(GCC_ROOT
    "/cvmfs/sft.cern.ch/lcg/releases/gcc/13.1.0-b3d18/x86_64-el9"
    CACHE PATH "GCC runtime installation")


# Build as a library (STATIC is simplest; change to SHARED if you want plugins)
add_library(LHCDataReaderAlg STATIC
  LHCDataReaderAlg.cpp
)
add_library(LumiRecordReaderAlg STATIC
  LumiRecordReaderAlg.cpp
)

# Public include paths for headers used by this target
target_include_directories(LHCDataReaderAlg
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FAIR_ROOT_DIR}
    "${TDAQ_EXTERNAL}/include"
    "${BOOST_ROOT}/include"
)

target_include_directories(LumiRecordReaderAlg
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FAIR_ROOT_DIR}
    "${TDAQ_EXTERNAL}/include"
    "${BOOST_ROOT}/include"
)

# Link directories for COOL/CORAL libraries
target_link_directories(LHCDataReaderAlg
  PUBLIC
    "${TDAQ_EXTERNAL}/lib"
    "${BOOST_ROOT}/lib"
    "${GCC_ROOT}/lib64"
)

target_link_directories(LumiRecordReaderAlg
  PUBLIC
    "${TDAQ_EXTERNAL}/lib"
    "${BOOST_ROOT}/lib"
    "${GCC_ROOT}/lib64"
)

# Link COOL/CORAL and Boost libraries
target_link_libraries(LHCDataReaderAlg
  PUBLIC
    lcg_CoolApplication
    lcg_CoolKernel
    lcg_RelationalCool
    lcg_CoralBase
    lcg_CoralKernel
    lcg_RelationalAccess
    boost_thread
    boost_system
    boost_chrono
    pthread
)

target_link_libraries(LumiRecordReaderAlg
  PUBLIC
    lcg_CoolApplication
    lcg_CoolKernel
    lcg_RelationalCool
    lcg_CoralBase
    lcg_CoralKernel
    lcg_RelationalAccess
    boost_thread
    boost_system
    boost_chrono
    pthread
)

# Use common compile options / include dirs from top-level interface lib (if present)
if(TARGET fair_options)
  target_link_libraries(LHCDataReaderAlg PUBLIC fair_options)
  target_link_libraries(LumiRecordReaderAlg PUBLIC fair_options)
endif()

# Optional: nice warnings locally
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(LHCDataReaderAlg PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(LumiRecordReaderAlg PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Set runtime library paths
set_target_properties(LHCDataReaderAlg
  PROPERTIES
    BUILD_RPATH "${TDAQ_EXTERNAL}/lib;${BOOST_ROOT}/lib;${GCC_ROOT}/lib64"
    INSTALL_RPATH "${TDAQ_EXTERNAL}/lib;${BOOST_ROOT}/lib;${GCC_ROOT}/lib64"
)

set_target_properties(LumiRecordReaderAlg
  PROPERTIES
    BUILD_RPATH "${TDAQ_EXTERNAL}/lib;${BOOST_ROOT}/lib;${GCC_ROOT}/lib64"
    INSTALL_RPATH "${TDAQ_EXTERNAL}/lib;${BOOST_ROOT}/lib;${GCC_ROOT}/lib64"
)

# Export an alias target name (clean usage)
add_library(FAIR::LHCDataReaderAlg ALIAS LHCDataReaderAlg)
add_library(FAIR::LumiRecordReaderAlg ALIAS LumiRecordReaderAlg)

target_link_libraries(calib_io INTERFACE 
    LHCDataReaderAlg
    LumiRecordReaderAlg
)
