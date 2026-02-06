#pragma once
#include "common/edm/RecoHit.hpp"
#include <vector>
#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"
#include <muParser.h> // for track selection string parsing
struct SimpleFittedTrack {
    double init_pos_x; // (x0, y0) at z=0
    double init_pos_y;
    double direction_x; // (dx/dz, dy/dz)
    double direction_y;
    double chi2_x; // (chi2_x, chi2_y)
    double chi2_y;
    int ndf = 0;
    std::vector<int> inTrackHitsIndices;
    std::vector<int> outTrackHitsIndices;
    int nTotalHits = 0;
    bool valid = false;
    std::vector<AHCALRecoHit> inTrackHits; // for FAIR internal use
    std::vector<AHCALRecoHit> outTrackHits; // for FAIR internal use
};

inline std::vector<FieldDesc> describe(const SimpleFittedTrack*) {
    return {
        field("init_pos_x", &SimpleFittedTrack::init_pos_x),
        field("init_pos_y", &SimpleFittedTrack::init_pos_y),
        field("direction_x", &SimpleFittedTrack::direction_x),
        field("direction_y", &SimpleFittedTrack::direction_y),
        field("chi2_x", &SimpleFittedTrack::chi2_x),
        field("chi2_y", &SimpleFittedTrack::chi2_y),
        field("ndf", &SimpleFittedTrack::ndf),
        field("inTrackHitsIndices", &SimpleFittedTrack::inTrackHitsIndices),
        field("outTrackHitsIndices", &SimpleFittedTrack::outTrackHitsIndices),
        field("nTotalHits", &SimpleFittedTrack::nTotalHits),
        field("valid", &SimpleFittedTrack::valid)
    };
}
AHCAL_REGISTER_IO_STRUCT(SimpleFittedTrack, "SimpleFittedTrack");


class SimpleFittedTrackCut {
public:
    explicit SimpleFittedTrackCut(const std::string& expr) {
        std::string e = expr;
        replace_all(e, "&&", " and ");
        replace_all(e, "||", " or ");
        replace_all(e, "!=",  " not "); 

        parser_.SetExpr(e);
        parser_.DefineVar("init_pos_x", &init_pos_x_);
        parser_.DefineVar("init_pos_y", &init_pos_y_);
        parser_.DefineVar("direction_x", &direction_x_);
        parser_.DefineVar("direction_y", &direction_y_);
        parser_.DefineVar("chi2_x", &chi2_x_);
        parser_.DefineVar("chi2_y", &chi2_y_);
        parser_.DefineVar("ndf", &ndf_);
        parser_.DefineVar("nTotalHits", &nTotalHits_);
        parser_.DefineVar("valid", &valid_);
    }

    bool eval(const SimpleFittedTrack& t) {
        init_pos_x_  = t.init_pos_x;
        init_pos_y_  = t.init_pos_y;
        direction_x_ = t.direction_x;
        direction_y_ = t.direction_y;
        chi2_x_      = t.chi2_x;
        chi2_y_      = t.chi2_y;
        ndf_         = static_cast<double>(t.ndf);
        nTotalHits_  = static_cast<double>(t.nTotalHits);
        valid_       = t.valid ? 1.0 : 0.0;
        return parser_.Eval() != 0.0;
    }

private:
    mu::Parser parser_;
    double init_pos_x_{0}, init_pos_y_{0}, direction_x_{0}, direction_y_{0};
    double chi2_x_{0}, chi2_y_{0}, ndf_{0}, nTotalHits_{0}, valid_{0};

    static void replace_all(std::string& s, const std::string& from, const std::string& to) {
        if (from.empty()) return;
        size_t pos = 0;
        while ((pos = s.find(from, pos)) != std::string::npos) {
            s.replace(pos, from.size(), to);
            pos += to.size();
        }
    }
};