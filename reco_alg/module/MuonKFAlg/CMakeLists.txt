# reco_alg/module/MuonKFAlg/CMakeLists.txt
cmake_minimum_required(VERSION 3.16)

# Build as a library (STATIC is simplest; change to SHARED if you want plugins)
add_library(MuonKFAlg STATIC
  MuonKFAlg.cpp
)

# Public include paths for headers used by this target
# - Project root/common is already exposed by fair_options in top-level
# - But we add it here as well to be safe if someone builds this target alone.
target_include_directories(MuonKFAlg
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}         # so "common/..." resolves
)


# Use common compile options / include dirs from top-level interface lib (if present)
if(TARGET fair_options)
  target_link_libraries(MuonKFAlg PUBLIC fair_options)
endif()

# Optional: nice warnings locally
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(MuonKFAlg PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Export an alias target name (clean usage)
add_library(FAIR::MuonKFAlg ALIAS MuonKFAlg)