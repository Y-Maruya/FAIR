# IO/
cmake_minimum_required(VERSION 3.16)

# Build as a library
add_library(PedestalAlg STATIC
  PedestalAlg.cpp
)

target_include_directories(PedestalAlg
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}         # so "common/..." resolves
)

# Use common compile options / include dirs from top-level interface lib (if present)
if(TARGET fair_options)
  target_link_libraries(PedestalAlg PUBLIC fair_options)
endif()

# Optional: nice warnings locally
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(PedestalAlg PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Export an alias target name (clean usage)
add_library(FAIR::PedestalAlg ALIAS PedestalAlg)