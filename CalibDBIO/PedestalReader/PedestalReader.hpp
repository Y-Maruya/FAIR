#pragma once
#include <string>
#include <unordered_map>

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

private:
    std::string runinfo_url = "https://faser-runinfo.app.cern.ch/cgibin/getRunInfo.py?runno=";
    std::unordered_map<int, Pedestal> pedestalMap_; // cellID to Pedestal
    std::string pedestal_runs_file_ = "pedestal_runs.txt";
    int selectNearestPedestalRun(int runNumber);
};