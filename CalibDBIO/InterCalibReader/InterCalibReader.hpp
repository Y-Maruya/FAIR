#pragma once
#include <memory>
#include <string>
#include <unordered_map>

namespace CalibDBIO {

struct HGLGRatio{ // Check the qualityflag != 1 before using the HGLG values
    double slope;
    double slopeerror;
    double intercept;
    double intercepterror;
    int qualityflag;
};

using HGLGRatioMap = std::unordered_map<int, HGLGRatio>;
using HGSaturationMap = std::unordered_map<int, double>;

class InterCalibReader {
public:
    InterCalibReader(int runNumber);
    ~InterCalibReader();

    void readHGLGRatios(int runNumber);
    HGLGRatio getHGLGRatio(int cellID) const;
    double getHGSaturationPoint(int cellID) const;
    const HGLGRatioMap& getHGLGRatioMap() const { return *hglgMap_; }
    const HGSaturationMap& getHGSaturationMap() const { return *hgSaturationMap_; } // Check the saturation point > 0 before using the saturation point values (not filled, so not needed to check the qualityflag)
    std::shared_ptr<const HGLGRatioMap> getHGLGRatioMapPtr() const { return hglgMap_; }
    std::shared_ptr<const HGSaturationMap> getHGSaturationMapPtr() const { return hgSaturationMap_; }

private:
    std::string runinfo_url = "https://faser-runinfo.app.cern.ch/cgibin/getRunInfo.py?runno=";
    std::shared_ptr<const HGSaturationMap> hgSaturationMap_; // cellID to HG saturation point
    std::shared_ptr<const HGLGRatioMap> hglgMap_; // cellID to HGLGRatio
};

} // namespace CalibDBIO
