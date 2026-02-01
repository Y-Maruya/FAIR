#pragma once
#include "common/AHCALGeometry.hpp"
#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"
#include <vector>

class MCInfo {
    public:
    //injected particle information
    int injected_pdgId; // PDG ID of the injected particle
    double injected_energy; // Energy of the injected particle in GeV
    double injected_time; // Time of the injected particle in ns
    double injected_z; // z position of the injection point in mm
    double injected_x; // x position of the injection point in mm
    double injected_y; // y position of the injection point in mm
    double injected_px; // x component of momentum of the injected particle in GeV/c
    double injected_py; // y component of momentum of the injected particle in GeV/c
    double injected_pz; // z component of momentum of the injected particle in GeV/c
    //neutrino interaction information
    bool isNeutrinoEvent; // true if the event contains a neutrino interaction
    int ftagNulabel; // neutrino flavor label (0: nueCC, 1: numuCC, 2: nutauCC, 3: NC, 4: other)
    bool isCC; // true if charged current interaction, false if neutral current interaction
    double Interaction_x; // x position of the neutrino interaction vertex in mm
    double Interaction_y; // y position of the neutrino interaction vertex in mm
    double Interaction_z; // z position of the neutrino interaction vertex in mm
    double SecondaryEnergy; // the particle energy from neutrino interaction in GeV
};

inline std::vector<FieldDesc> describe(const MCInfo&) {
    return {
        field("injected_pdgId", &MCInfo::injected_pdgId),
        field("injected_energy", &MCInfo::injected_energy),
        field("injected_time", &MCInfo::injected_time),
        field("injected_z", &MCInfo::injected_z),
        field("injected_x", &MCInfo::injected_x),
        field("injected_y", &MCInfo::injected_y),
        field("injected_px", &MCInfo::injected_px),
        field("injected_py", &MCInfo::injected_py),
        field("injected_pz", &MCInfo::injected_pz),
        field("isNeutrinoEvent", &MCInfo::isNeutrinoEvent),
        field("ftagNulabel", &MCInfo::ftagNulabel),
        field("isCC", &MCInfo::isCC),
        field("Interaction_x", &MCInfo::Interaction_x),
        field("Interaction_y", &MCInfo::Interaction_y),
        field("Interaction_z", &MCInfo::Interaction_z),
        field("SecondaryEnergy", &MCInfo::SecondaryEnergy)
    };
}
AHCAL_REGISTER_IO_STRUCT(MCInfo, "MCInfo");


