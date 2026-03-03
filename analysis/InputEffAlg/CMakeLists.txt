# IO/
cmake_minimum_required(VERSION 3.16)

# Build as a library
add_library(InputEffAlg STATIC
  InputEffAlg.cpp
)

target_include_directories(InputEffAlg
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FAIR_ROOT_DIR}         # so "common/..." resolves
)

# Use common compile options / include dirs from top-level interface lib (if present)
if(TARGET fair_options)
  target_link_libraries(InputEffAlg PUBLIC fair_options)
endif()

# Optional: nice warnings locally
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(InputEffAlg PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Export an alias target name (clean usage)
add_library(FAIR::InputEffAlg ALIAS InputEffAlg)