#include "PCAClusterAlg.hpp"
#include "common/AlgRegistry.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

AHCAL_REGISTER_ALG(AHCALRecoAlg::PCAClusterAlg, "PCAClusterAlg")

namespace AHCALRecoAlg {

namespace {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

static Vec3 matVec(
    double a00, double a01, double a02,
    double a11, double a12,
    double a22,
    const Vec3& v
) {
    return {
        a00 * v.x + a01 * v.y + a02 * v.z,
        a01 * v.x + a11 * v.y + a12 * v.z,
        a02 * v.x + a12 * v.y + a22 * v.z
    };
}

static double norm(const Vec3& v) {
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

static void normalize(Vec3& v) {
    const double n = norm(v);
    if (n <= 0.0) return;
    v.x /= n;
    v.y /= n;
    v.z /= n;
}

static double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static Vec3 powerIteration(
    double a00, double a01, double a02,
    double a11, double a12,
    double a22
) {
    Vec3 v{1.0, 1.0, 1.0};
    normalize(v);

    for (int iter = 0; iter < 50; ++iter) {
        Vec3 av = matVec(a00, a01, a02, a11, a12, a22, v);
        normalize(av);
        v = av;
    }

    return v;
}

static double rayleighQuotient(
    double a00, double a01, double a02,
    double a11, double a12,
    double a22,
    const Vec3& v
) {
    const Vec3 av = matVec(a00, a01, a02, a11, a12, a22, v);
    return dot(v, av);
}

} // namespace

void PCAClusterAlg::initialize() {}

void PCAClusterAlg::finalize() {}

void PCAClusterAlg::execute(EventStore& evt) {
    auto clusters = evt.get<std::vector<Cluster>>(cfg_.in_cluster_key);

    for (auto& c : clusters) {
        computePCA(c);
    }

    evt.put(cfg_.out_cluster_key, std::move(clusters));
}

void PCAClusterAlg::computePCA(Cluster& c) const {
    const int nHits = static_cast<int>(c.hits.size());

    if (nHits < cfg_.min_hits) {
        c.pca_lambda1 = 0.0;
        c.pca_lambda2 = 0.0;
        c.pca_lambda3 = 0.0;
        c.pca_linearity = 0.0;
        c.pca_width = 0.0;
        c.is_track_like = false;
        return;
    }

    double sum_w = 0.0;
    double mean_x = 0.0;
    double mean_y = 0.0;
    double mean_z = 0.0;

    for (const auto& h : c.hits) {
        const double w = cfg_.energy_weighted ? std::max(h.Edep, 0.0) : 1.0;

        sum_w += w;
        mean_x += w * h.Xpos();
        mean_y += w * h.Ypos();
        mean_z += w * h.Zpos();
    }

    if (sum_w <= 0.0) {
        c.is_track_like = false;
        return;
    }

    mean_x /= sum_w;
    mean_y /= sum_w;
    mean_z /= sum_w;

    c.pca_mean_x = mean_x;
    c.pca_mean_y = mean_y;
    c.pca_mean_z = mean_z;

    double sxx = 0.0;
    double syy = 0.0;
    double szz = 0.0;
    double sxy = 0.0;
    double sxz = 0.0;
    double syz = 0.0;

    for (const auto& h : c.hits) {
        const double w = cfg_.energy_weighted ? std::max(h.Edep, 0.0) : 1.0;

        const double dx = h.Xpos() - mean_x;
        const double dy = h.Ypos() - mean_y;
        const double dz = h.Zpos() - mean_z;

        sxx += w * dx * dx;
        syy += w * dy * dy;
        szz += w * dz * dz;
        sxy += w * dx * dy;
        sxz += w * dx * dz;
        syz += w * dy * dz;
    }

    sxx /= sum_w;
    syy /= sum_w;
    szz /= sum_w;
    sxy /= sum_w;
    sxz /= sum_w;
    syz /= sum_w;

    Vec3 v1 = powerIteration(sxx, sxy, sxz, syy, syz, szz);
    double lambda1 = rayleighQuotient(sxx, sxy, sxz, syy, syz, szz, v1);

    // Deflation: A2 = A - lambda1 * v1 v1^T
    const double bxx = sxx - lambda1 * v1.x * v1.x;
    const double byy = syy - lambda1 * v1.y * v1.y;
    const double bzz = szz - lambda1 * v1.z * v1.z;
    const double bxy = sxy - lambda1 * v1.x * v1.y;
    const double bxz = sxz - lambda1 * v1.x * v1.z;
    const double byz = syz - lambda1 * v1.y * v1.z;

    Vec3 v2 = powerIteration(bxx, bxy, bxz, byy, byz, bzz);
    double lambda2 = rayleighQuotient(sxx, sxy, sxz, syy, syz, szz, v2);

    const double trace = sxx + syy + szz;
    double lambda3 = trace - lambda1 - lambda2;

    if (lambda1 < lambda2) std::swap(lambda1, lambda2);
    if (lambda2 < lambda3) std::swap(lambda2, lambda3);
    if (lambda1 < lambda2) std::swap(lambda1, lambda2);

    lambda1 = std::max(lambda1, 0.0);
    lambda2 = std::max(lambda2, 0.0);
    lambda3 = std::max(lambda3, 0.0);

    c.pca_axis_x = v1.x;
    c.pca_axis_y = v1.y;
    c.pca_axis_z = v1.z;

    c.pca_lambda1 = lambda1;
    c.pca_lambda2 = lambda2;
    c.pca_lambda3 = lambda3;

    if (lambda1 > 0.0) {
        c.pca_linearity = 1.0 - lambda2 / lambda1;
    } else {
        c.pca_linearity = 0.0;
    }

    c.pca_width = std::sqrt(lambda2 + lambda3);

    c.is_track_like =
        c.pca_linearity > cfg_.track_linearity_threshold &&
        c.pca_width < cfg_.track_width_threshold_mm;
}

void PCAClusterAlg::parse_cfg(const YAML::Node& n) {
    cfg_.in_cluster_key =
        n["in_cluster_key"].as<std::string>(cfg_.in_cluster_key);

    cfg_.out_cluster_key =
        n["out_cluster_key"].as<std::string>(cfg_.out_cluster_key);

    cfg_.min_hits =
        n["min_hits"].as<int>(cfg_.min_hits);

    cfg_.track_linearity_threshold =
        n["track_linearity_threshold"].as<double>(cfg_.track_linearity_threshold);

    cfg_.track_width_threshold_mm =
        n["track_width_threshold_mm"].as<double>(cfg_.track_width_threshold_mm);

    cfg_.energy_weighted =
        n["energy_weighted"].as<bool>(cfg_.energy_weighted);
}

} // namespace AHCALRecoAlg