//
#pragma once
#include "IO/Descriptor.hpp"
class AHCALRawHit {
public:
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

    int hg_adc;        // 12-bit
    int lg_adc;        // 12-bit
    int hittag = 0;
    int bcid = 0;
    int index = 0; // for internal use
};

inline std::vector<FieldDesc> describe(const AHCALRawHit*) {
    return {
        field("cellID", &AHCALRawHit::cellID),
        field("hg_adc", &AHCALRawHit::hg_adc),
        field("lg_adc", &AHCALRawHit::lg_adc),
        field("hittag", &AHCALRawHit::hittag),
        field("bcid", &AHCALRawHit::bcid),
    };
}
inline std::vector<FieldDescVector> describe_vector(const AHCALRawHit*) {
    return {
        field_vector("v.cellID", &AHCALRawHit::cellID),
        field_vector("v.hg_adc", &AHCALRawHit::hg_adc),
        field_vector("v.lg_adc", &AHCALRawHit::lg_adc),
        field_vector("v.hittag", &AHCALRawHit::hittag),
        field_vector("v.bcid", &AHCALRawHit::bcid),
    };
}