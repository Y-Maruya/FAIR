#pragma once
#include <memory>
#include <string>
#include <unordered_map>

namespace CalibDBIO {
struct Pedestal{
    double HighGainPeak;
    double HighGainPeakError;
    double HighGainSigma;
    double HighGainSigmaError;
    double LowGainPeak;
    double LowGainPeakError;
    double LowGainSigma;
    double LowGainSigmaError;
    int HighGainStatus;
    int LowGainStatus;
    int HighGainUsable;
    int LowGainUsable;
};

using PedestalMap = std::unordered_map<int, Pedestal>;

class PedestalReader {
public:
    PedestalReader(int runNumber);
    ~PedestalReader();

    void readPedestals(int runNumber);
    Pedestal getPedestal(int cellID) const;
    const PedestalMap& getPedestalMap() const { return *pedestalMap_; }
    std::shared_ptr<const PedestalMap> getPedestalMapPtr() const { return pedestalMap_; }

private:
    std::string runinfo_url = "https://faser-runinfo.app.cern.ch/cgibin/getRunInfo.py?runno=";
    std::shared_ptr<const PedestalMap> pedestalMap_; // cellID to Pedestal
    std::string pedestal_runs_file_ = "pedestal_runs.txt";
    int selectNearestPedestalRun(int runNumber);
};

} // namespace CalibDBIO
