#pragma once
#include "common/AHCALGeometry.hpp"
#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"
#include <vector> 
class AHCALSimHit {
public:
    // int layer;      // 0..39
    // int asic;       // 0..8
    // int channel;    // 0..35
    int cellID;    // layer*100000 + asic*10000 + channel
    int layer() const {
        return cellID / 100000;
    }
    int asic() const {
        return (cellID / 10000) % 10;
    }
    int channel() const {
        return cellID % 10000;
    }
    int chip() const {
        return asic();
    }
    double Edep = 0.0;      // in MeV
    double Nmip = 0.0;      // in MIP
    double HitTime = 0.0;      // in ns
    double TimeOfArrival = 0.0;   // in ns
    double Xpos() const {
        return AHCALGeometry::Pos_X(this->channel(),this->asic());
    }
    double Ypos() const {
        return AHCALGeometry::Pos_Y(this->channel(),this->asic());
    }
    double Zpos() const {
        return AHCALGeometry::Pos_Z(this->layer());
    }
    int Xindex() const {
        return static_cast<int>(Xpos() / 40.3 + 9.0);
    }
    int Yindex() const {
        return static_cast<int>(Ypos() / 40.3 + 9.0);
    }
    int index = 0; // for internal use
};
inline std::vector<FieldDesc> describe(const AHCALSimHit*) {
    return {
        field("cellID", &AHCALSimHit::cellID),
        field("Edep", &AHCALSimHit::Edep),
        field("Nmip", &AHCALSimHit::Nmip),
        field("HitTime", &AHCALSimHit::HitTime),
        field("TimeOfArrival", &AHCALSimHit::TimeOfArrival)
    };
}
inline std::vector<FieldDescVector> describe_vector(const AHCALSimHit*) {
    return {
        field_vector("v.cellID", &AHCALSimHit::cellID),
        field_vector("v.Edep", &AHCALSimHit::Edep),
        field_vector("v.Nmip", &AHCALSimHit::Nmip),
        field_vector("v.HitTime", &AHCALSimHit::HitTime),
        field_vector("v.TimeOfArrival", &AHCALSimHit::TimeOfArrival)
    };
}
AHCAL_REGISTER_IO_STRUCT(AHCALSimHit, "AHCALSimHit");
AHCAL_REGISTER_IO_STRUCT_VECTOR(AHCALSimHit, "vector<AHCALSimHit>");