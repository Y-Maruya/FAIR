#pragma once
#include "common/AHCALGeometry.hpp"
#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"
#include <vector> 

class AHCALTrueHit {
public:
    int cellID;
    int MCcellID;
    float energy;
    float time;
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

inline std::vector<FieldDesc> describe(const AHCALTrueHit*) {
    return {
        field("cellID", &AHCALTrueHit::cellID),
        field("MCcellID", &AHCALTrueHit::MCcellID),
        field("energy", &AHCALTrueHit::energy),
        field("time", &AHCALTrueHit::time),
    };
}

inline std::vector<FieldDescVector> describe_vector(const AHCALTrueHit*) {
    return {
        field_vector("v.cellID", &AHCALTrueHit::cellID),
        field_vector("v.MCcellID", &AHCALTrueHit::MCcellID),
        field_vector("v.energy", &AHCALTrueHit::energy),
        field_vector("v.time", &AHCALTrueHit::time),
    };
}
AHCAL_REGISTER_IO_STRUCT(AHCALTrueHit, "AHCALTrueHit");
AHCAL_REGISTER_IO_STRUCT_VECTOR(AHCALTrueHit, "vector<AHCALTrueHit>");