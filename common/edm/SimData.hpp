#pragma once
#include "common/AHCALGeometry.hpp"
#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"
#include <vector>

class SimData {
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
    int Secondarypdgid; // PDG ID of the secondary particle from neutrino interaction
    double SecondaryEnergy; // the particle energy from neutrino interaction in GeV
    double Secondary_px; // x component of momentum of the particle from neutrino interaction in GeV/c
    double Secondary_py; // y component of momentum of the particle from neutrino interaction in GeV/c
    double Secondary_pz; // z component of momentum of the particle from neutrino
};

inline std::vector<FieldDesc> describe(const SimData*) {
    return {
        field("injected_pdgId", &SimData::injected_pdgId),
        field("injected_energy", &SimData::injected_energy),
        field("injected_time", &SimData::injected_time),
        field("injected_z", &SimData::injected_z),
        field("injected_x", &SimData::injected_x),
        field("injected_y", &SimData::injected_y),
        field("injected_px", &SimData::injected_px),
        field("injected_py", &SimData::injected_py),
        field("injected_pz", &SimData::injected_pz),
        field("isNeutrinoEvent", &SimData::isNeutrinoEvent),
        field("ftagNulabel", &SimData::ftagNulabel),
        field("isCC", &SimData::isCC),
        field("Interaction_x", &SimData::Interaction_x),
        field("Interaction_y", &SimData::Interaction_y),
        field("Interaction_z", &SimData::Interaction_z),
        field("SecondaryEnergy", &SimData::SecondaryEnergy),
        field("Secondarypdgid", &SimData::Secondarypdgid),
        field("Secondary_px", &SimData::Secondary_px),
        field("Secondary_py", &SimData::Secondary_py),
        field("Secondary_pz", &SimData::Secondary_pz)
    };
}
AHCAL_REGISTER_IO_STRUCT(SimData, "SimData");


