#pragma once
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

namespace CalibDBIO {
enum DecisionClass {
  kNoData = -1,
  kInsufficientStatistics = 0,
  kFirstFitSuccess = 1,
  kThresholdFitSuccess = 2,
  kGoodFitWithoutThreshold = 3,
  kSecondandFirstFitFailed = 4,
  kSecondFitFailedFirstFitSuccess = 5,
  kThresholdFitFailedButNeeded = 6,
  kNDecisionClasses = 7
};

enum StateContainer {
  kInsufficientStatisticsState = 1<<0,
  kFirstFitFailed_BigChi2Ndf = 1<<1,
  kFirstFitFailed_ParameterAtLimit = 1<<2,
  kFirstFitFailed_FitStatusNotZero = 1<<3,
  kFirstFitFailed_Others = 1<<4,
  kFirstFitSuccess_LowGausSigma = 1<<5,
  kFirstFitSuccess_KnownEfficiencyDegraded = 1<<6,
  kFirstFitSuccessState = 1<<7,
  kSecondFitFailed_RefitStatusNotZero = 1<<8,
  kSecondFitFailed_BigChi2Ndf = 1<<9,
  kSecondFitFailed_ParameterAtLimit = 1<<10,
  kSecondFitFailed_Others = 1<<11,
  kSecondFitSuccess_LowMPV = 1<<12,
  kSecondFitSuccess_LowGausSigma = 1<<13,
  kSecondFitSuccess_LowWidth = 1<<14,
  kSecondFitSuccess_HighMPVError = 1<<15,
  kSecondFitSuccessState = 1<<16,
  kThresholdFitFailed_BigChi2Ndf = 1<<29,
  kThresholdFitFailed_BigWidth = 1<<17,
  kThresholdFitFailed_BigWidthRatio = 1<<28,
  kThresholdFitFailed_LowLeftMean = 1<<18,
  kThresholdFitFailed_LowRightMean = 1<<19,
  kThresholdFitFailed_ParameterAtLimit = 1<<20,
  kThresholdFitFailed_BigUncertainty = 1<<30,
  kThresholdFitFailed_Others = 1<<21,
  kThresholdFitFailed_and_GoodEffMPVStability = 1<<22,
  kThresholdFitFailed_and_BadEffMPVStability = 1<<23,
  kThresholdFitSuccessState = 1<<24,
  kEfficiencyRatioOver1p0 = 1<<25,
  kEfficiencyRatioOver0p98 = 1<<26,
  kEfficiencyRatioUnder0p98 = 1<<27,
  kNStateContainers = 31
};


struct MIP{ // Check the decision != kNoData before using the MIP values
    double mpv; // most probable value
    double width; // width of the Landau distribution
    double gaussigma; // sigma of the Gaussian convolution
    double mpverror; // error on the most probable value
    double widtherror; // error on the width
    double gaussigmaerror; // error on the Gaussian sigma
    int decision; // decision on the MIP quality 
    int state; // state of the MIP fitting
    int fullfit; // flag indicating if the MIP value was obtained from a full fit
    int imputed; // flag indicating if the MIP value was averaged of other runs
    int ref; // layer avg is used as reference for the MIP value
};
struct Threshold{ // Check the threshold > 0 before using the threshold values
    double threshold; // threshold value
    double thresholderror; // error on the threshold value
    double thresholdwidth; // width of the threshold distribution
    double thresholdwidtherror; // error on the threshold width
    int decision; // decision on the threshold quality
    int state; // state of the threshold fitting
    int fullfit; // flag indicating if the threshold value was obtained from a full fit
    int imputed; // flag indicating if the threshold value was averaged of other runs
    int ref; // layer avg is used as reference for the threshold value
};

using MIPMap = std::unordered_map<int, MIP>;
using ThresholdMap = std::unordered_map<int, Threshold>;
using MIPMPVMap = std::unordered_map<int, std::pair<double, double>>;

class MIPReader {
public:
    MIPReader(int runNumber);
    ~MIPReader();

    void readMIPs(int runNumber);
    MIP getMIP(int cellID) const;
    Threshold getThreshold(int cellID) const;
    const MIPMap& getMIPMap() const { return *mipMap_; }
    const ThresholdMap& getThresholdMap() const { return *thresholdMap_; }
    const MIPMPVMap& getMIPMPVMap() const { return *mipMPVMap_; }
    std::shared_ptr<const MIPMap> getMIPMapPtr() const { return mipMap_; }
    std::shared_ptr<const ThresholdMap> getThresholdMapPtr() const { return thresholdMap_; }
    std::shared_ptr<const MIPMPVMap> getMIPMPVMapPtr() const { return mipMPVMap_; }

private:
    std::string runinfo_url = "https://faser-runinfo.app.cern.ch/cgibin/getRunInfo.py?runno=";
    std::shared_ptr<const MIPMap> mipMap_; // cellID to MIP
    std::shared_ptr<const ThresholdMap> thresholdMap_; // cellID to Threshold
    std::shared_ptr<const MIPMPVMap> mipMPVMap_; // cellID to MIP MPV
};

} // namespace CalibDBIO
