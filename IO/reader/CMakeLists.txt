# IO/
cmake_minimum_required(VERSION 3.16)

# Build as a library
add_library(RootRawHitReader STATIC
  RootRawHitReader.cpp
)
add_library(BinaryRawHitReader STATIC
  BinaryRawHitReader.cpp
)
add_library(SimHitReader STATIC
  SimHitReader.cpp
)

if(FAIR_WITH_GEANT4 AND Geant4_FOUND)
  add_library(Geant4SimReader STATIC
    Geant4SimReader.cpp
    ${FAIR_ROOT_DIR}/simulation/geant4/ActionInitialization.cpp
    ${FAIR_ROOT_DIR}/simulation/geant4/DetectorConstruction.cpp
    ${FAIR_ROOT_DIR}/simulation/geant4/HcalSD.cpp
    ${FAIR_ROOT_DIR}/simulation/geant4/Pythia8Decayer.cpp
    ${FAIR_ROOT_DIR}/simulation/geant4/Pythia8DecayerPhysics.cpp
    ${FAIR_ROOT_DIR}/simulation/geant4/PrimaryGeneratorAction.cpp
    ${FAIR_ROOT_DIR}/simulation/geant4/SimulationEngine.cpp
  )
endif()

target_include_directories(RootRawHitReader
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FAIR_ROOT_DIR}         # so "common/..." resolves
)

target_include_directories(BinaryRawHitReader
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FAIR_ROOT_DIR}         # so "common/..." resolves
)

target_include_directories(SimHitReader
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FAIR_ROOT_DIR}         # so "common/..." resolves
)

if(FAIR_WITH_GEANT4 AND Geant4_FOUND)
  target_include_directories(Geant4SimReader
    PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${FAIR_ROOT_DIR}
      ${Geant4_INCLUDE_DIRS}
  )
endif()

# Use common compile options / include dirs from top-level interface lib (if present)
if(TARGET fair_options)
  target_link_libraries(RootRawHitReader PUBLIC fair_options)
  target_link_libraries(BinaryRawHitReader PUBLIC fair_options)
  target_link_libraries(SimHitReader PUBLIC fair_options)
  if(FAIR_WITH_GEANT4 AND Geant4_FOUND)
    target_link_libraries(Geant4SimReader PUBLIC fair_options ${Geant4_LIBRARIES})
  endif()
endif()

# Optional: nice warnings locally
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(RootRawHitReader PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(BinaryRawHitReader PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(SimHitReader PRIVATE -Wall -Wextra -Wpedantic)
  if(FAIR_WITH_GEANT4 AND Geant4_FOUND)
    target_compile_options(Geant4SimReader PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endif()

# Export an alias target name (clean usage)
add_library(FAIR::RootRawHitReader ALIAS RootRawHitReader)
add_library(FAIR::BinaryRawHitReader ALIAS BinaryRawHitReader)
add_library(FAIR::SimHitReader ALIAS SimHitReader)
if(FAIR_WITH_GEANT4 AND Geant4_FOUND)
  add_library(FAIR::Geant4SimReader ALIAS Geant4SimReader)
endif()
