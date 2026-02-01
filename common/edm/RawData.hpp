#pragma once
#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"


struct AHCALTLURawData{
    int Timestamp;
    int BCID_TLU;
    std::vector<int> Inputs = std::vector<int>(6,0);
    std::vector<int> FineTimestamps = std::vector<int>(6,0);
    int RunNo;
    int CycleID;
    int TriggerID;
    int Event_Time;
};

inline std::vector<FieldDesc> describe(const AHCALTLURawData*){
    return {
        field("Timestamp", &AHCALTLURawData::Timestamp),
        field("BCID_TLU", &AHCALTLURawData::BCID_TLU),
        field("Inputs", &AHCALTLURawData::Inputs),
        field("FineTimestamps", &AHCALTLURawData::FineTimestamps),
        field("RunNo", &AHCALTLURawData::RunNo),
        field("CycleID", &AHCALTLURawData::CycleID),
        field("TriggerID", &AHCALTLURawData::TriggerID),
        field("Event_Time", &AHCALTLURawData::Event_Time),
    };
}

AHCAL_REGISTER_IO_STRUCT(AHCALTLURawData, "AHCALTLURawData");