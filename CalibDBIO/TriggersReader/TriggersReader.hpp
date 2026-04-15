#pragma once
#include <string>
#include <unordered_map>

namespace CalibDBIO {
struct Triggers{
    int AnyVetoed = -1;
    int NonVetoed = -1;
    int NonPhysical = -1;
    int Physical = -1;
    int PreVeto = -1;
    int PostVeto = -1;
    int Recorded = -1;
};
class TriggersReader {
public:
    TriggersReader(int runNumber);
    ~TriggersReader();

    void readTriggers(int runNumber);
    Triggers getTriggers();

private:
    int m_runNumber;
    std::string runinfo_url = "https://faser-runinfo.app.cern.ch/cgibin/getRunInfo.py?runno=";
    Triggers triggers_;
};

} // namespace CalibDBIO