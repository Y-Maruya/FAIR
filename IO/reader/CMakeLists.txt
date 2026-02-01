# IO/
cmake_minimum_required(VERSION 3.16)

# Build as a library
add_library(RootRawHitReader STATIC
  RootRawHitReader.cpp
)
add_library(BinaryRawHitReader STATIC
  BinaryRawHitReader.cpp
)

target_include_directories(RootRawHitReader
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}         # so "common/..." resolves
)

target_include_directories(BinaryRawHitReader
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}         # so "common/..." resolves
)

# Use common compile options / include dirs from top-level interface lib (if present)
if(TARGET fair_options)
  target_link_libraries(RootRawHitReader PUBLIC fair_options)
  target_link_libraries(BinaryRawHitReader PUBLIC fair_options)
endif()

# Optional: nice warnings locally
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(RootRawHitReader PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(BinaryRawHitReader PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Export an alias target name (clean usage)
add_library(FAIR::RootRawHitReader ALIAS RootRawHitReader)
add_library(FAIR::BinaryRawHitReader ALIAS BinaryRawHitReader)