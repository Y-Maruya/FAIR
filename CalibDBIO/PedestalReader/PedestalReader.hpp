#pragma once
#include <string>
#include <unordered_map>

namespace CalibDBIO {
struct Pedestal{
    double HighGainPeak;
    double HighGainSigma;
    double LowGainPeak;
    double LowGainSigma;
    int HighGainStatus;
    int LowGainStatus;
};
class PedestalReader {
public:
    PedestalReader(int runNumber);
    ~PedestalReader();

    void readPedestals(int runNumber);
    Pedestal getPedestal(int cellID);
    const std::unordered_map<int, Pedestal>& getPedestalMap() const { return pedestalMap_; }

private:
    std::string runinfo_url = "https://faser-runinfo.app.cern.ch/cgibin/getRunInfo.py?runno=";
    std::unordered_map<int, Pedestal> pedestalMap_; // cellID to Pedestal
    std::string pedestal_runs_file_ = "pedestal_runs.txt";
    int selectNearestPedestalRun(int runNumber);
};

} // namespace CalibDBIO