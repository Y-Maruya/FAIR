#ifndef AHCAL_RECO_ALG_TRACK_FIT_ALG_HPP
#define AHCAL_RECO_ALG_TRACK_FIT_ALG_HPP
#include "common/EventStore.hpp"
#include "common/edm/EDM.hpp"
#include <string>
#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace AHCALRecoAlg {
    class TrackFitAlg final : public IAlg { // final to prevent inheritance
    public:
        TrackFitAlg(std::string in_recohit_key,
            std::string out_track_key)
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

#endif // AHCAL_RECO_ALG_TRACK_FIT_ALG_HPP