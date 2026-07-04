# IO/
cmake_minimum_required(VERSION 3.16)

# Build as a library
add_library(AdcToEnergyReadTTreeAlg STATIC
  AdcToEnergyReadTTreeAlg.cpp
)

add_library(AdcToEnergyAlg STATIC
  AdcToEnergyAlg.cpp
)

target_include_directories(AdcToEnergyReadTTreeAlg
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FAIR_ROOT_DIR}         # so "common/..." resolves
)

target_include_directories(AdcToEnergyAlg
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FAIR_ROOT_DIR}
)

# Use common compile options / include dirs from top-level interface lib (if present)
if(TARGET fair_options)
  target_link_libraries(AdcToEnergyReadTTreeAlg PUBLIC fair_options)
  target_link_libraries(AdcToEnergyAlg PUBLIC fair_options)
endif()

# Optional: nice warnings locally
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(AdcToEnergyReadTTreeAlg PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(AdcToEnergyAlg PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Export an alias target name (clean usage)
add_library(FAIR::AdcToEnergyReadTTreeAlg ALIAS AdcToEnergyReadTTreeAlg)
add_library(FAIR::AdcToEnergyAlg ALIAS AdcToEnergyAlg)
