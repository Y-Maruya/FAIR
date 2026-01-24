#pragma once
#include "common/AHCALGeometry.hpp"
#include "IO/Descriptor.hpp"
class AHCALRecoHit {
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
};
inline std::vector<FieldDesc> describe(const AHCALRecoHit*) {
    return {
        field("cellID", &AHCALRecoHit::cellID),
        field("Edep", &AHCALRecoHit::Edep),
        field("Nmip", &AHCALRecoHit::Nmip),
    };
}
inline std::vector<FieldDescVector> describe_vector(const AHCALRecoHit*) {
    return {
        field_vector("v.cellID", &AHCALRecoHit::cellID),
        field_vector("v.Edep", &AHCALRecoHit::Edep),
        field_vector("v.Nmip", &AHCALRecoHit::Nmip),
    };
}