#include "MuonKFAlg.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"

#include <array>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>

AHCAL_REGISTER_ALG(AHCALRecoAlg::MuonKFAlg, "MuonKFAlg")
// using namespace AHCALRecoAlg;
namespace AHCALRecoAlg{

// --------------------------------------------
// parse_skip_layers
// --------------------------------------------
std::bitset<40> parse_skip_layers(std::vector<int> skipLayers) {
  std::bitset<40> mask;
  for (int L : skipLayers) {
    if (L >= 0 && L < 40) {
      mask.set(L);
    }
  }
  return mask;
}


// --------------------------------------------
// Minimal linear algebra helpers (4D state, 2D meas)
// state: (x, y, tx, ty)
// meas : (x, y)
// --------------------------------------------
struct TrackInternal {
  double xv[4] = {0,0,0,0};      // x,y,tx,ty
  double C[4][4] = {{0}};        // covariance
  double z = 0.0;
  double chi2 = 0.0;
  int ndof = 0;
  int consecutive_skips = 0;
  std::vector<const AHCALRecoHit*> used;
};

static inline double default_sigma_mm() {
  // pitch ~ 40.3 mm from your Xindex/Yindex conversion; uniform -> sigma = pitch/sqrt(12)
  return AHCALGeometry::xy_size / std::sqrt(12.0);
}

static inline void mat4_eye(double M[4][4]) {
  for(int i=0;i<4;i++) for(int j=0;j<4;j++) M[i][j] = (i==j)?1.0:0.0;
}

static inline void mat4_mul(const double A[4][4], const double B[4][4], double R[4][4]) {
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      double s=0;
      for(int k=0;k<4;k++) s += A[i][k]*B[k][j];
      R[i][j]=s;
    }
  }
}

static inline void mat4_transpose(const double A[4][4], double R[4][4]) {
  for(int i=0;i<4;i++) for(int j=0;j<4;j++) R[i][j]=A[j][i];
}

static inline void mat4_add(const double A[4][4], const double B[4][4], double R[4][4]) {
  for(int i=0;i<4;i++) for(int j=0;j<4;j++) R[i][j]=A[i][j]+B[i][j];
}

static inline void mat4_sub(const double A[4][4], const double B[4][4], double R[4][4]) {
  for(int i=0;i<4;i++) for(int j=0;j<4;j++) R[i][j]=A[i][j]-B[i][j];
}

static inline void mat4_vec4_mul(const double A[4][4], const double x[4], double r[4]) {
  for(int i=0;i<4;i++){
    double s=0;
    for(int j=0;j<4;j++) s += A[i][j]*x[j];
    r[i]=s;
  }
}

static inline void propagate(TrackInternal& trk, double z_to, double sigmaTheta) {
  const double dz = z_to - trk.z;

  // F
  double F[4][4]; mat4_eye(F);
  F[0][2]=dz;
  F[1][3]=dz;

  // x <- F x
  double xnew[4];
  mat4_vec4_mul(F, trk.xv, xnew);
  for(int i=0;i<4;i++) trk.xv[i]=xnew[i];

  // C <- F C F^T + Q
  double Ft[4][4]; mat4_transpose(F, Ft);
  double FC[4][4]; mat4_mul(F, trk.C, FC);
  double FCFt[4][4]; mat4_mul(FC, Ft, FCFt);

  double Q[4][4] = {{0}};
  Q[2][2] = sigmaTheta*sigmaTheta;
  Q[3][3] = sigmaTheta*sigmaTheta;

  mat4_add(FCFt, Q, trk.C);

  trk.z = z_to;
}

// compute d2 = r^T S^-1 r for 2x2 S
static inline double d2_2d(double rx, double ry,
                           double S00, double S01, double S10, double S11) {
  const double det = S00*S11 - S01*S10;
  const double inv00 =  S11 / det;
  const double inv01 = -S01 / det;
  const double inv10 = -S10 / det;
  const double inv11 =  S00 / det;
  return rx*(inv00*rx + inv01*ry) + ry*(inv10*rx + inv11*ry);
}

static inline bool update_with_hit(TrackInternal& trk, const AHCALRecoHit& h,
                                   double sigma_xy_mm, double gateD2) {
  const double mx = h.Xpos();
  const double my = h.Ypos();

  const double hx = trk.xv[0];
  const double hy = trk.xv[1];

  const double rx = mx - hx;
  const double ry = my - hy;

  // S = top-left 2x2 of C + R
  double S00 = trk.C[0][0] + sigma_xy_mm*sigma_xy_mm;
  double S01 = trk.C[0][1];
  double S10 = trk.C[1][0];
  double S11 = trk.C[1][1] + sigma_xy_mm*sigma_xy_mm;

  const double det = S00*S11 - S01*S10;
  if (std::abs(det) < 1e-24) return false;

  const double inv00 =  S11 / det;
  const double inv01 = -S01 / det;
  const double inv10 = -S10 / det;
  const double inv11 =  S00 / det;

  const double d2 = rx*(inv00*rx + inv01*ry) + ry*(inv10*rx + inv11*ry);
  if (d2 > gateD2) return false;

  // K = C H^T S^-1  (4x2) where H picks (x,y)
  // First take C[:,0..1] then multiply by S^-1
  double K[4][2];
  for(int i=0;i<4;i++){
    const double c0 = trk.C[i][0];
    const double c1 = trk.C[i][1];
    K[i][0] = c0*inv00 + c1*inv10;
    K[i][1] = c0*inv01 + c1*inv11;
  }

  // x <- x + K r
  trk.xv[0] += K[0][0]*rx + K[0][1]*ry;
  trk.xv[1] += K[1][0]*rx + K[1][1]*ry;
  trk.xv[2] += K[2][0]*rx + K[2][1]*ry;
  trk.xv[3] += K[3][0]*rx + K[3][1]*ry;

  // C <- (I - K H) C
  // (I - K H) only affects columns 0,1:
  // M = I; M(i,0)-=K(i,0); M(i,1)-=K(i,1)
  double M[4][4]; mat4_eye(M);
  for(int i=0;i<4;i++){
    M[i][0] -= K[i][0];
    M[i][1] -= K[i][1];
  }
  double MC[4][4];
  mat4_mul(M, trk.C, MC);
  for(int i=0;i<4;i++) for(int j=0;j<4;j++) trk.C[i][j]=MC[i][j];

  trk.chi2 += d2;
  trk.ndof += 2;
  trk.consecutive_skips = 0;
  trk.used.push_back(&h);
  return true;
}

static inline const AHCALRecoHit* pick_nearest_hit(
    const std::vector<const AHCALRecoHit*>& hits,
    double xpred, double ypred,
    const MuonKFAlgCfg& cfg) {
  const AHCALRecoHit* best = nullptr;
  double best_d2 = std::numeric_limits<double>::infinity();
  for (const auto* h : hits) {
    if (cfg.useNmipWindow) {
      if (h->Nmip < cfg.nmipMin || h->Nmip > cfg.nmipMax) continue;
    }
    const double dx = h->Xpos() - xpred;
    const double dy = h->Ypos() - ypred;
    const double d2 = dx*dx + dy*dy;
    if (d2 < best_d2) {
      best_d2 = d2;
      best = h;
    }
  }
  return best;
}

static inline std::vector<const AHCALRecoHit*> topK_for_seed(
    std::vector<const AHCALRecoHit*> hits,
    int K,
    const MuonKFAlgCfg& cfg) {
  if (cfg.useNmipWindow) {
    hits.erase(std::remove_if(hits.begin(), hits.end(), [&](auto* h){
      return (h->Nmip < cfg.nmipMin || h->Nmip > cfg.nmipMax);
    }), hits.end());
  }
  std::sort(hits.begin(), hits.end(), [](auto* a, auto* b){
    return std::abs(a->Nmip - 1.0) < std::abs(b->Nmip - 1.0);
  });
  if ((int)hits.size() > K) hits.resize(K);
  return hits;
}

// --------------------------------------------
// Public function: find_muon_track_kf
// --------------------------------------------
bool find_muon_track_kf(const std::vector<AHCALRecoHit>& recoHits,
                        Track& bestOut,
                        const MuonKFAlgCfg& cfg) {
  if (recoHits.empty()) return false;

  // group by layer, respecting skipLayer
  std::array<std::vector<const AHCALRecoHit*>, 40> byLayer;
  for (auto& v : byLayer) v.clear();

  int maxLayer = -1;
  for (const auto& h : recoHits) {
    const int L = h.layer();
    if (L < 0 || L >= 40) continue;
    if (cfg.skipLayer.test(L)) continue;
    byLayer[L].push_back(&h);
    maxLayer = std::max(maxLayer, L);
  }
  if (maxLayer < 0) return false;

  const int Lend   = maxLayer;
  const int Lstart = std::max(0, Lend - cfg.lastNLayers + 1);

  // build active layer list (hit exists + not skipped)
  std::vector<int> layers;
  layers.reserve(cfg.lastNLayers);
  for (int L = Lstart; L <= Lend; ++L) {
    if (cfg.skipLayer.test(L)) continue;
    if (!byLayer[L].empty()) layers.push_back(L);
  }
  if (layers.size() < 3) return false;

  const double sigma_xy = (cfg.measSigmaXY_mm > 0.0) ? cfg.measSigmaXY_mm : default_sigma_mm();

  // choose seed back layer (last layer with hits)
  const int L2 = layers.back();
  const auto hitsL2 = topK_for_seed(byLayer[L2], cfg.maxSeedHitsPerLayer, cfg);
  if (hitsL2.empty()) return false;

  // choose candidate earlier layers for seeds with sufficient gap
  std::vector<int> seedL1s;
  for (int i = (int)layers.size() - 2; i >= 0; --i) {
    if (L2 - layers[i] >= cfg.seedLayerGap) seedL1s.push_back(layers[i]);
    if ((int)seedL1s.size() >= 4) break;
  }
  if (seedL1s.empty()) return false;

  bool found = false;
  double bestScore = std::numeric_limits<double>::infinity();
  TrackInternal bestInt;

  const double z2 = AHCALGeometry::Pos_Z(L2);

  for (int L1 : seedL1s) {
    const auto hitsL1 = topK_for_seed(byLayer[L1], cfg.maxSeedHitsPerLayer, cfg);
    if (hitsL1.empty()) continue;

    const double z1 = AHCALGeometry::Pos_Z(L1);
    const double dz = z2 - z1;
    if (std::abs(dz) < 1e-6) continue;

    for (const auto* h1 : hitsL1) {
      for (const auto* h2 : hitsL2) {

        TrackInternal trk;
        trk.z = z2;
        trk.xv[0] = h2->Xpos();
        trk.xv[1] = h2->Ypos();
        trk.xv[2] = (h2->Xpos() - h1->Xpos()) / dz;
        trk.xv[3] = (h2->Ypos() - h1->Ypos()) / dz;

        // init covariance
        for(int i=0;i<4;i++) for(int j=0;j<4;j++) trk.C[i][j]=0.0;
        trk.C[0][0] = sigma_xy*sigma_xy;
        trk.C[1][1] = sigma_xy*sigma_xy;
        const double slope0 = 0.05; // 50 mrad prior (loose)
        trk.C[2][2] = slope0*slope0;
        trk.C[3][3] = slope0*slope0;

        trk.used.clear();
        trk.used.push_back(h2);
        trk.chi2 = 0.0;
        trk.ndof = 0;
        trk.consecutive_skips = 0;

        // iterate backward over active layers (excluding L2 itself at the end)
        for (int idx = (int)layers.size() - 2; idx >= 0; --idx) {
          const int L = layers[idx];
          const double z = AHCALGeometry::Pos_Z(L);

          propagate(trk, z, cfg.sigmaTheta);

          const double xpred = trk.xv[0];
          const double ypred = trk.xv[1];

          const AHCALRecoHit* hbest = pick_nearest_hit(byLayer[L], xpred, ypred, cfg);
          if (!hbest) {
            trk.consecutive_skips++;
            if (trk.consecutive_skips > cfg.maxConsecutiveSkips) break;
            continue;
          }

          const bool ok = update_with_hit(trk, *hbest, sigma_xy, cfg.gateD2);
          if (!ok) {
            trk.consecutive_skips++;
            if (trk.consecutive_skips > cfg.maxConsecutiveSkips) break;
            continue;
          }
        }

        const int nUsed = (int)trk.used.size();
        if (nUsed < cfg.minUsedLayers) continue;

        const double chi2ndof = (trk.ndof > 0) ? (trk.chi2 / trk.ndof) : 1e9;
        const double score = chi2ndof + 2.0 / (double)nUsed;

        if (score < bestScore) {
          bestScore = score;
          bestInt = trk;
          found = true;
        }
      }
    }
  }

  if (!found) return false;

  // export to public Track
  bestOut = Track{};
  bestOut.x  = bestInt.xv[0];
  bestOut.y  = bestInt.xv[1];
  bestOut.tx = bestInt.xv[2];
  bestOut.ty = bestInt.xv[3];
  bestOut.z  = bestInt.z;
  bestOut.chi2 = bestInt.chi2;
  bestOut.ndof = bestInt.ndof;
  bestOut.consecutive_skips = bestInt.consecutive_skips;
  bestOut.inTrackHitsIndices.clear();
  bestOut.outTrackHitsIndices.clear();
  bestOut.valid = true;
  for (const auto* h : bestInt.used) {
    bestOut.inTrackHits.push_back(*h);
    bestOut.inTrackHitsIndices.push_back(h->index);
    bestOut.nInTrackHits++;
  }
  for (const auto& h : recoHits) {
    bool used = false;
    for (const auto* uh : bestInt.used) {
      if (h.index == uh->index) {
        used = true;
        break;
      }
    }
    if (!used) {
      bestOut.outTrackHits.push_back(h);
      bestOut.outTrackHitsIndices.push_back(h.index);
      bestOut.nOutTrackHits++;
    }
  }

  return true;
}

// --------------------------------------------
// MuonKFAlg (RecoAlg-style wrapper)
// --------------------------------------------
void MuonKFAlg::execute(EventStore & evt) {
  auto recohits = evt.get<std::vector<AHCALRecoHit>>(m_cfg.in_recohit_key);
  m_last_found = false;
  m_last = Track{};

  if (recohits.empty()) {
    evt.put<Track>(m_cfg.out_track_key, Track{});
    LOG_DEBUG("MuonKFAlg: No input reco hits found.");
    return;
  }

  Track trk;
  if (find_muon_track_kf(recohits, trk, m_cfg)) {
    m_last_found = true;
    m_last = std::move(trk);
    evt.put<Track>(m_cfg.out_track_key, m_last);
  }else{
    trk = Track{};
    m_last_found = false;
    evt.put<Track>(m_cfg.out_track_key, Track{});
  }
}

void MuonKFAlg::parse_cfg(const YAML::Node& n) {
  m_cfg.in_recohit_key = get_or<std::string>(n, "in_recohit_key", m_cfg.in_recohit_key);
  m_cfg.out_track_key = get_or<std::string>(n, "out_track_key", m_cfg.out_track_key);
  m_cfg.lastNLayers = get_or<int>(n, "lastNLayers", m_cfg.lastNLayers);
  m_cfg.minUsedLayers = get_or<int>(n, "minUsedLayers", m_cfg.minUsedLayers);
  m_cfg.maxConsecutiveSkips = get_or<int>(n, "maxConsecutiveSkips", m_cfg.maxConsecutiveSkips);
  m_cfg.useNmipWindow = get_or<bool>(n, "useNmipWindow", m_cfg.useNmipWindow);
  m_cfg.nmipMin = get_or<double>(n, "nmipMin", m_cfg.nmipMin);
  m_cfg.nmipMax = get_or<double>(n, "nmipMax", m_cfg.nmipMax);
  m_cfg.measSigmaXY_mm = get_or<double>(n, "measSigmaXY_mm", m_cfg.measSigmaXY_mm);
  m_cfg.sigmaTheta = get_or<double>(n, "sigmaTheta", m_cfg.sigmaTheta);
  m_cfg.gateD2 = get_or<double>(n, "gateD2", m_cfg.gateD2);
  m_cfg.seedLayerGap = get_or<int>(n, "seedLayerGap", m_cfg.seedLayerGap);
  m_cfg.maxSeedHitsPerLayer = get_or<int>(n, "maxSeedHitsPerLayer", m_cfg.maxSeedHitsPerLayer);
  std::vector<int> skipLayers = get_or<std::vector<int>>(n, "skipLayers", {0,2,14});
  m_cfg.skipLayer = parse_skip_layers(skipLayers);
}
} // namespace AHCALRecoAlg
