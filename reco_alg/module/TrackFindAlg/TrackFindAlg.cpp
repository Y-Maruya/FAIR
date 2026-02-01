#include "TrackFindAlg.hpp"

#include <algorithm>
#include <set>

namespace AHCALRecoAlg {

namespace {
inline bool isSkippedLayer(int layer, const std::vector<int>& skip) {
    return std::find(skip.begin(), skip.end(), layer) != skip.end();
}
} // namespace

void TrackFindAlg::execute(EventStore& evt) {
    // NOTE: EventStore API はプロジェクト側の実装に合わせて調整してください。
    // 想定: evt.get<T>(key) -> T& , evt.put<T>(key, value)

    auto& hits = evt.get<std::vector<AHCALRecoHit>>(m_in_recohit_key);

    std::vector<int> isoFlag;
    computeIsolatedFlags(hits, isoFlag);

    std::vector<int> trackLabel, trackLength;
    int numTrack = 0;
    buildTracks(hits, isoFlag, trackLabel, trackLength, numTrack);

    std::vector<TrackSummary> sum;
    summarizeTracks(hits, trackLabel, numTrack, sum);

    // pick main muon track (macro 同等: Loose の中で nLayers 最大、同点なら chi2 最小)
    int mainId = 0;
    int bestLayers = -1;
    double bestChi = 1e30;
    bool mainIsTight = false;

    for (const auto& ts : sum) {
        if (!ts.isLoose) continue;
        const double chi = ts.chi2x + ts.chi2y;
        if (ts.nLayers > bestLayers || (ts.nLayers == bestLayers && chi < bestChi)) {
            bestLayers = ts.nLayers;
            bestChi = chi;
            mainId = ts.label;
            mainIsTight = ts.isTight;
        }
    }

    find = (mainId > 0);

    std::vector<AHCALRecoHit> muonHits;
    if (mainId > 0) {
        muonHits.reserve(hits.size());
        for (size_t i = 0; i < hits.size(); ++i) {
            const bool inMain = (trackLabel[i] == mainId);
            const bool okForMuon =
                m_use_tight_only ? (mainIsTight && inMain) : inMain; // macro の useTightOnly 相当（最小版）

            if (okForMuon) muonHits.push_back(hits[i]);
        }
    }

    evt.put<std::vector<AHCALRecoHit>>(m_out_track_key, std::move(muonHits));
}

void TrackFindAlg::computeIsolatedFlags(const std::vector<AHCALRecoHit>& hits,
                                       std::vector<int>& isoFlag) const
{
    // macro: cone=42mm, 同レイヤ内に近傍が無ければ isolated=1
    const double cone = 42.0;

    const int n = static_cast<int>(hits.size());
    isoFlag.assign(n, 0);

    for (int i = 0; i < n; ++i) {
        if (hits[i].Nmip <= m_thr_mip_frac) continue;

        const int li = hits[i].layer();
        if (isSkippedLayer(li, m_skiplayers)) continue;

        bool hasNeighbor = false;
        const double xi = hits[i].Xpos();
        const double yi = hits[i].Ypos();

        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            if (hits[j].Nmip <= m_thr_mip_frac) continue;
            if (hits[j].layer() != li) continue;

            const double dx = hits[j].Xpos() - xi;
            const double dy = hits[j].Ypos() - yi;
            if (std::fabs(dx) < cone && std::fabs(dy) < cone) { hasNeighbor = true; break; }
        }

        isoFlag[i] = hasNeighbor ? 0 : 1;
    }
}

void TrackFindAlg::buildTracks(const std::vector<AHCALRecoHit>& hits,
                               const std::vector<int>& isoFlag,
                               std::vector<int>& trackLabel,
                               std::vector<int>& trackLength,
                               int& numTrack) const
{
    // macro FollowTrackFromSeed をそのまま C++ 化（前進: layer+1..+3）
    const double c1 = 42.0, c2 = 84.0, c3 = 126.0;

    const int n = static_cast<int>(hits.size());
    trackLabel.assign(n, 0);
    trackLength.assign(n, 0);

    std::vector<int> scanned(n, 0);
    std::vector<int> order(n);
    for (int i = 0; i < n; ++i) order[i] = i;

    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return hits[a].layer() > hits[b].layer(); });

    auto followFromSeed = [&](int seedIdx, int label) {
        int cur = seedIdx;

        while (true) {
            int bestK = -1;
            int bestDL = 9999;
            double bestR2 = 1e30;

            const int ls = hits[cur].layer();
            const double xs = hits[cur].Xpos();
            const double ys = hits[cur].Ypos();

            for (int k = 0; k < n; ++k) {
                if (k == cur) continue;
                if (hits[k].Nmip <= m_thr_mip_frac) continue;
                if (isoFlag[k] == 0) continue;
                if (scanned[k] != 0) continue;

                const int lk = hits[k].layer();
                if (isSkippedLayer(lk, m_skiplayers)) continue;

                const int dL = lk - ls;
                if (dL <= 0 || dL > 3) continue;

                const double dx = hits[k].Xpos() - xs;
                const double dy = hits[k].Ypos() - ys;

                bool in = false;
                if (dL == 1) in = (std::fabs(dx) < c1 && std::fabs(dy) < c1);
                if (dL == 2) in = (std::fabs(dx) < c2 && std::fabs(dy) < c2);
                if (dL == 3) in = (std::fabs(dx) < c3 && std::fabs(dy) < c3);
                if (!in) continue;

                const double r2 = dx * dx + dy * dy;
                if (dL < bestDL || (dL == bestDL && r2 < bestR2)) {
                    bestDL = dL;
                    bestR2 = r2;
                    bestK = k;
                }
            }

            if (bestK < 0) break;

            scanned[bestK] = 1;
            trackLabel[bestK] = label;
            trackLength[bestK] = trackLength[cur] + bestDL;
            cur = bestK;
        }
    };

    numTrack = 0;
    for (int idx : order) {
        if (hits[idx].Nmip <= m_thr_mip_frac) continue;
        if (isoFlag[idx] == 0) continue;
        if (isSkippedLayer(hits[idx].layer(), m_skiplayers)) continue;

        if (scanned[idx] == 0 && trackLabel[idx] == 0) {
            ++numTrack;
            const int label = numTrack;

            scanned[idx] = 1;
            trackLabel[idx] = label;
            trackLength[idx] = 0;

            followFromSeed(idx, label);
        }
    }
}

void TrackFindAlg::summarizeTracks(const std::vector<AHCALRecoHit>& hits,
                                  const std::vector<int>& trackLabel,
                                  int numTrack,
                                  std::vector<TrackSummary>& out) const
{
    out.clear();
    if (numTrack <= 0) return;

    std::vector<std::vector<int>> idxs(numTrack + 1);
    idxs.shrink_to_fit(); // no-op-ish; keeps behavior explicit

    for (int i = 0; i < static_cast<int>(hits.size()); ++i) {
        const int lab = trackLabel[i];
        if (lab > 0 && lab <= numTrack) idxs[lab].push_back(i);
    }

    for (int lab = 1; lab <= numTrack; ++lab) {
        if (idxs[lab].empty()) continue;

        TrackSummary ts;
        ts.label = lab;
        ts.nHits = static_cast<int>(idxs[lab].size());

        std::set<int> layers;
        double sumZ = 0, sumZZ = 0, sumX = 0, sumZX = 0, sumY = 0, sumZY = 0;
        double sumQ = 0, sumQ2 = 0;

        for (int idx : idxs[lab]) {
            const auto& h = hits[idx];
            layers.insert(h.layer());

            const double z = h.Zpos();
            const double x = h.Xpos();
            const double y = h.Ypos();

            sumZ += z; sumZZ += z * z;
            sumX += x; sumZX += z * x;
            sumY += y; sumZY += z * y;

            const double q = h.Nmip; // macro の (HG-Ped)/MIP に相当する量として Nmip を使用
            sumQ += q; sumQ2 += q * q;
        }

        ts.nLayers = static_cast<int>(layers.size());
        const int N = ts.nHits;

        if (N >= 2) {
            const double denom = N * sumZZ - sumZ * sumZ;
            double ax = 0, bx = 0, ay = 0, by = 0;

            if (std::fabs(denom) > 0) {
                ax = (N * sumZX - sumZ * sumX) / denom;
                bx = (sumX - ax * sumZ) / N;
                ay = (N * sumZY - sumZ * sumY) / denom;
                by = (sumY - ay * sumZ) / N;
            }

            ts.tanTheta = std::sqrt(ax * ax + ay * ay);

            double s2x = 0, s2y = 0;
            for (int idx : idxs[lab]) {
                const auto& h = hits[idx];
                const double z = h.Zpos();
                const double xfit = ax * z + bx;
                const double yfit = ay * z + by;
                const double dx = h.Xpos() - xfit;
                const double dy = h.Ypos() - yfit;
                s2x += dx * dx;
                s2y += dy * dy;
            }

            if (N > 2) {
                const double tile = 30.0;
                ts.chi2x = (s2x / (tile * tile)) / (N - 2);
                ts.chi2y = (s2y / (tile * tile)) / (N - 2);
            } else {
                ts.chi2x = 0;
                ts.chi2y = 0;
            }
        }

        ts.meanQ = (N > 0 ? sumQ / N : 0.0);
        double var = (N > 0 ? (sumQ2 / N - ts.meanQ * ts.meanQ) : 0.0);
        if (var < 0) var = 0;
        ts.rmsQ = std::sqrt(var);

        // macro と同じカット
        ts.isLoose =
            (ts.nLayers >= 8) && (ts.nHits >= 8) &&
            (ts.tanTheta < 0.9) &&
            (ts.meanQ > 0.5 && ts.meanQ < 1.5) &&
            (ts.rmsQ < 0.8);

        ts.isTight =
            ts.isLoose &&
            (ts.nLayers >= 15) &&
            (ts.tanTheta < 0.5) &&
            (ts.chi2x < 5.0) && (ts.chi2y < 5.0) &&
            (ts.meanQ > 0.8 && ts.meanQ < 1.2) &&
            (ts.rmsQ < 0.5);

        out.push_back(ts);
    }
}

} // namespace AHCALRecoAlg

