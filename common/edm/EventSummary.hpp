#pragma once
#include "common/AHCALGeometry.hpp"
#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"
#include <vector> 
struct EventSummary{
    int event_counter = 0;
    int nHits = 0;
    int nHitAbove0p5MIP = 0;
    int nHitLayers = 0;
    int nHitLayersAbove0p5MIP = 0;
    double TotalEnergy = 0.0;
    double TotalnMIP = 0.0;
    double MaxEnergy = 0.0;
    double MaxnMIP = 0.0;
    double CenterOfGravityX = 0.0;
    double CenterOfGravityY = 0.0;
    double CenterOfGravityZ = 0.0;
    double RMSX = 0.0;
    double RMSY = 0.0;
    double RMSZ = 0.0;
};

inline std::vector<FieldDesc> describe(const EventSummary*) {
    return {
        field("event_counter", &EventSummary::event_counter),
        field("nHits", &EventSummary::nHits),
        field("nHitAbove0p5MIP", &EventSummary::nHitAbove0p5MIP),
        field("nHitLayers", &EventSummary::nHitLayers),
        field("nHitLayersAbove0p5MIP", &EventSummary::nHitLayersAbove0p5MIP),
        field("TotalEnergy", &EventSummary::TotalEnergy),
        field("TotalnMIP", &EventSummary::TotalnMIP),
        field("MaxEnergy", &EventSummary::MaxEnergy),
        field("MaxnMIP", &EventSummary::MaxnMIP),
        field("CenterOfGravityX", &EventSummary::CenterOfGravityX),
        field("CenterOfGravityY", &EventSummary::CenterOfGravityY),
        field("CenterOfGravityZ", &EventSummary::CenterOfGravityZ),
        field("RMSX", &EventSummary::RMSX),
        field("RMSY", &EventSummary::RMSY),
        field("RMSZ", &EventSummary::RMSZ)
    };
}
AHCAL_REGISTER_IO_STRUCT(EventSummary, "EventSummary");