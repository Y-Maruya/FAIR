add_library(FairGeant4Core STATIC
  ActionInitialization.cpp
  Geant4Alg.cpp
  DetectorConstruction.cpp
  HcalSD.cpp
  PrimaryGeneratorAction.cpp
  Pythia8Decayer.cpp
  Pythia8DecayerPhysics.cpp
  SimulationEngine.cpp
)

target_include_directories(FairGeant4Core
  PUBLIC
    ${FAIR_ROOT_DIR}
)

if(TARGET fair_options)
  target_link_libraries(FairGeant4Core PUBLIC fair_options)
endif()

target_link_libraries(FairGeant4Core PUBLIC ${Geant4_LIBRARIES})
target_include_directories(FairGeant4Core PUBLIC ${Geant4_INCLUDE_DIRS})

add_library(FAIR::Geant4Core ALIAS FairGeant4Core)
