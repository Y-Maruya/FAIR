#pragma once

#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"
#include "common/edm/RecoHit.hpp"
#include <vector>

namespace AHCALRecoAlg {

struct Cluster {
    int cluster_id = -1;
    std::vector<AHCALRecoHit> hits;
    std::vector<int> index;
    int nHitsInCluster = 0;
    int nHitsUnder0p5MIPInCluster = 0;
    int nHitsUnder1p5MIPInCluster = 0;
    int nHitsOver2p5MIPInCluster = 0;
    double totalEdep = 0.0;
    double totalEdepError = 0.0;
    double avgEdep = 0.0;
    double avgEdepError = 0.0;
    double pca_mean_x = 0.0;
    double pca_mean_y = 0.0;
    double pca_mean_z = 0.0;

    double pca_axis_x = 0.0;
    double pca_axis_y = 0.0;
    double pca_axis_z = 1.0;

    double pca_lambda1 = 0.0;
    double pca_lambda2 = 0.0;
    double pca_lambda3 = 0.0;

    double pca_linearity = 0.0; // 1 - lambda2/lambda1
    double pca_width = 0.0;     // sqrt(lambda2 + lambda3)
    bool is_track_like = false;
    void updateSummary();
};

inline void Cluster::updateSummary() {
    nHitsInCluster = static_cast<int>(hits.size());

    nHitsUnder0p5MIPInCluster = 0;
    nHitsUnder1p5MIPInCluster = 0;
    nHitsOver2p5MIPInCluster = 0;

    totalEdep = 0.0;
    totalEdepError = 0.0;

    for (const auto& h : hits) {
        totalEdep += h.Edep;
        totalEdepError += h.EdepError * h.EdepError;

        if (h.Nmip < 0.5) nHitsUnder0p5MIPInCluster++;
        if (h.Nmip < 1.5) nHitsUnder1p5MIPInCluster++;
        if (h.Nmip > 2.5) nHitsOver2p5MIPInCluster++;
    }

    totalEdepError = std::sqrt(totalEdepError);
    avgEdep = nHitsInCluster > 0 ? totalEdep / nHitsInCluster : 0.0;
    avgEdepError = nHitsInCluster > 0 ? totalEdepError / nHitsInCluster : 0.0;
}
inline std::vector<FieldDesc> describe(const Cluster*) {
    return {
        optional_field("index", &Cluster::index),
        field("nHitsInCluster", &Cluster::nHitsInCluster),
        field("totalEdep", &Cluster::totalEdep),
        field("avgEdep", &Cluster::avgEdep),
        field("totalEdepError", &Cluster::totalEdepError),
        field("avgEdepError", &Cluster::avgEdepError),
        field("nHitsUnder0p5MIPInCluster", &Cluster::nHitsUnder0p5MIPInCluster),
        field("nHitsUnder1p5MIPInCluster", &Cluster::nHitsUnder1p5MIPInCluster),
        field("nHitsOver2p5MIPInCluster", &Cluster::nHitsOver2p5MIPInCluster),
        field("cluster_id", &Cluster::cluster_id),
        field("pca_mean_x", &Cluster::pca_mean_x),
        field("pca_mean_y", &Cluster::pca_mean_y),
        field("pca_mean_z", &Cluster::pca_mean_z),
        field("pca_axis_x", &Cluster::pca_axis_x),
        field("pca_axis_y", &Cluster::pca_axis_y),
        field("pca_axis_z", &Cluster::pca_axis_z),
        field("pca_lambda1", &Cluster::pca_lambda1),
        field("pca_lambda2", &Cluster::pca_lambda2),
        field("pca_lambda3", &Cluster::pca_lambda3),
        field("pca_linearity", &Cluster::pca_linearity),
        field("pca_width", &Cluster::pca_width),
        field("is_track_like", &Cluster::is_track_like),
    };
}

inline std::vector<FieldDescVector> describe_vector(const Cluster*) {
    return {
        field_vector("v.cluster_id", &Cluster::cluster_id),
        optional_field_vector("v.index", &Cluster::index),
        field_vector("v.nHitsInCluster", &Cluster::nHitsInCluster),
        field_vector("v.totalEdep", &Cluster::totalEdep),
        field_vector("v.avgEdep", &Cluster::avgEdep),
        field_vector("v.totalEdepError", &Cluster::totalEdepError),
        field_vector("v.avgEdepError", &Cluster::avgEdepError),
        field_vector("v.nHitsUnder0p5MIPInCluster", &Cluster::nHitsUnder0p5MIPInCluster),
        field_vector("v.nHitsUnder1p5MIPInCluster", &Cluster::nHitsUnder1p5MIPInCluster),
        field_vector("v.nHitsOver2p5MIPInCluster", &Cluster::nHitsOver2p5MIPInCluster),
        field_vector("v.pca_mean_x", &Cluster::pca_mean_x),
        field_vector("v.pca_mean_y", &Cluster::pca_mean_y),
        field_vector("v.pca_mean_z", &Cluster::pca_mean_z),
        field_vector("v.pca_axis_x", &Cluster::pca_axis_x),
        field_vector("v.pca_axis_y", &Cluster::pca_axis_y),
        field_vector("v.pca_axis_z", &Cluster::pca_axis_z),
        field_vector("v.pca_lambda1", &Cluster::pca_lambda1),
        field_vector("v.pca_lambda2", &Cluster::pca_lambda2),
        field_vector("v.pca_lambda3", &Cluster::pca_lambda3),
        field_vector("v.pca_linearity", &Cluster::pca_linearity),
        field_vector("v.pca_width", &Cluster::pca_width),
        field_vector("v.is_track_like", &Cluster::is_track_like),
    };
}

AHCAL_REGISTER_IO_STRUCT(Cluster, "Cluster");
AHCAL_REGISTER_IO_STRUCT_VECTOR(Cluster, "vector<Cluster>");

} // namespace AHCALRecoAlg
