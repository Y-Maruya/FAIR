#pragma once
#include <vector>
#include "common/edm/RecoHit.hpp"
#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"
#include <muParser.h> // for track selection string parsing
struct Track {
    // state at the last updated z in the fit loop: (x,y,tx,ty)
    double x = 0.0;
    double y = 0.0;
    double tx = 0.0;
    double ty = 0.0;
    double z = 0.0;
    double chi2 = 0.0;
    int ndof = 0;
    int consecutive_skips = 0;
    bool valid = false;
    int nInTrackHits = 0;
    int nOutTrackHits = 0;
    std::vector<int> inTrackHitsIndices;
    std::vector<int> outTrackHitsIndices;
    std::vector<AHCALRecoHit> inTrackHits; // for FAIR internal use
    std::vector<AHCALRecoHit> outTrackHits; // for FAIR internal use
};

inline std::vector<FieldDesc> describe(const Track*) {
    return {
        field("x", &Track::x),
        field("y", &Track::y),
        field("tx", &Track::tx),
        field("ty", &Track::ty),
        field("z", &Track::z),
        field("chi2", &Track::chi2),
        field("ndof", &Track::ndof),
        field("consecutive_skips", &Track::consecutive_skips),
        field("nInTrackHits", &Track::nInTrackHits),
        field("nOutTrackHits", &Track::nOutTrackHits),
        field("inTrackHitsIndices", &Track::inTrackHitsIndices),
        field("outTrackHitsIndices", &Track::outTrackHitsIndices),
        field("valid", &Track::valid)
    };
}
inline std::vector<FieldDescVector> describe_vector(const Track*) {
    return {
        field_vector("v.x", &Track::x),
        field_vector("v.y", &Track::y),
        field_vector("v.tx", &Track::tx),
        field_vector("v.ty", &Track::ty),
        field_vector("v.z", &Track::z),
        field_vector("v.chi2", &Track::chi2),
        field_vector("v.ndof", &Track::ndof),
        field_vector("v.consecutive_skips", &Track::consecutive_skips),
        field_vector("v.nInTrackHits", &Track::nInTrackHits),
        field_vector("v.nOutTrackHits", &Track::nOutTrackHits),
        field_vector("v.inTrackHitsIndices", &Track::inTrackHitsIndices),
        field_vector("v.outTrackHitsIndices", &Track::outTrackHitsIndices),
        field_vector("v.valid", &Track::valid)
    };
}
AHCAL_REGISTER_IO_STRUCT(Track, "Track");
AHCAL_REGISTER_IO_STRUCT_VECTOR(Track, "vector<Track>");

class TrackCut {
public:
    explicit TrackCut(const std::string& expr) {
        std::string e = expr;
        replace_all(e, "&&", " and ");
        replace_all(e, "||", " or ");
        replace_all(e, "!=",  " not "); 

        parser_.SetExpr(e);
        parser_.DefineVar("x", &x_);
        parser_.DefineVar("y", &y_);
        parser_.DefineVar("tx", &tx_);
        parser_.DefineVar("ty", &ty_);
        parser_.DefineVar("z", &z_);
        parser_.DefineVar("chi2", &chi2_);
        parser_.DefineVar("ndof", &ndof_);
        parser_.DefineVar("consecutive_skips", &consecutive_skips_);
        parser_.DefineVar("nInTrackHits", &nInTrackHits_);
        parser_.DefineVar("nOutTrackHits", &nOutTrackHits_);
        parser_.DefineVar("valid", &valid_);
    }

    bool eval(const Track& t) {
        x_ = t.x;
        y_ = t.y;
        tx_ = t.tx;
        ty_ = t.ty;
        z_ = t.z;
        chi2_ = t.chi2;
        ndof_ = static_cast<double>(t.ndof);
        consecutive_skips_ = static_cast<double>(t.consecutive_skips);
        nInTrackHits_ = static_cast<double>(t.nInTrackHits);
        nOutTrackHits_ = static_cast<double>(t.nOutTrackHits);
        valid_ = t.valid ? 1.0 : 0.0;
        return parser_.Eval() != 0.0;
    }

private:
    mu::Parser parser_;
    double x_{0}, y_{0}, tx_{0}, ty_{0};
    double z_{0}, chi2_{0}, ndof_{0}, consecutive_skips_{0}, nInTrackHits_{0}, nOutTrackHits_{0}, valid_{0};

    static void replace_all(std::string& s, const std::string& from, const std::string& to) {
        if (from.empty()) return;
        size_t pos = 0;
        while ((pos = s.find(from, pos)) != std::string::npos) {
            s.replace(pos, from.size(), to);
            pos += to.size();
        }
    }
};