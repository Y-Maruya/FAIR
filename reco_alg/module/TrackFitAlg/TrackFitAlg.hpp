#include "common/EventStore.hpp"
#include "common/Hit.hpp"
#include <string>
#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace AHCALRecoAlg {
    struct Track {
        std::pair <double, double> init_pos; // (x0, y0) at z=0
        std::pair <double, double> direction; // (dx/dz, dy/dz)
        std::pair <double, double> chi2; // (chi2_x, chi2_y)
        int ndf = 0;
        std::vector<AHCALRecoAndRawHit> inTrackHits;
        std::vector<AHCALRecoAndRawHit> outTrackHits;
        int nTotalHits = 0;
    };

    class TrackFitAlg final : public IAlg { // final to prevent inheritance
    public:
        TrackFitAlg(std::string in_recohit_key, std::string out_track_key)
            : m_in_recohit_key(std::move(in_recohit_key)),
            m_out_track_key(std::move(out_track_key)) {}

        void execute(EventStore& evt) override;
        void setThresholdXY(double threshold) {
            threshold_xy = threshold;
        }
    private:
        std::string m_in_recohit_key;
        std::string m_out_track_key;
        double threshold_xy = 40./2;
    };
} // namespace AHCALRecoAlg