# reco_alg/module/TrackFitAlg/CMakeLists.txt
cmake_minimum_required(VERSION 3.16)

# Build as a library (STATIC is simplest; change to SHARED if you want plugins)
add_library(TrackFitAlg STATIC
  TrackFitAlg.cpp
)

# Public include paths for headers used by this target
# - Project root/common is already exposed by fair_options in top-level
# - But we add it here as well to be safe if someone builds this target alone.
target_include_directories(TrackFitAlg
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}         # so "common/..." resolves
)

find_package(ROOT REQUIRED COMPONENTS Core Graf Graf3d Hist MathCore RIO Tree)
target_link_libraries(TrackFitAlg PUBLIC
  ROOT::Core
  ROOT::Graf
  ROOT::Graf3d
  ROOT::Hist
  ROOT::MathCore
  ROOT::RIO
  ROOT::Tree
)

# Use common compile options / include dirs from top-level interface lib (if present)
if(TARGET fair_options)
  target_link_libraries(TrackFitAlg PUBLIC fair_options)
endif()

# If you use any faser-common targets inside TrackFitAlg.cpp, link them here.
# (Uncomment once you know which are needed.)
# Example possibilities:
# if(TARGET Logging)
#   target_link_libraries(TrackFitAlg PUBLIC Logging)
# endif()
# if(TARGET Exceptions)
#   target_link_libraries(TrackFitAlg PUBLIC Exceptions)
# endif()
# if(TARGET EventFormats)
#   target_link_libraries(TrackFitAlg PUBLIC EventFormats)
# endif()

# Optional: nice warnings locally
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(TrackFitAlg PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Export an alias target name (clean usage)
add_library(FAIR::TrackFitAlg ALIAS TrackFitAlg)