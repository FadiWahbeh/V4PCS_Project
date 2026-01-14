//
// Created by Fabien WAHBEH on 22/12/2025.
// Mise à jour : Implémentation conforme à V4PCS (Base 4 points + Invariants)
//

#ifndef V4PCS_ALGO_H
#define V4PCS_ALGO_H

#include <vector>
#include <unordered_map>
#include <Eigen/Core>
#include "patch_builder.h"

struct Base4PCS {
    const patch_planaire* p1;
    const patch_planaire* p2;
    const patch_planaire* p3;
    const patch_planaire* p4;

    // ratios d'intersection
    // r1 = ||p1 - e|| / ||p1 - p2||
    // r2 = ||p3 - e|| / ||p3 - p4||
    double r1;
    double r2;

    // Distances pour filtrage rapide
    // d1 = ||p1 - p2||, d2 = ||p3 - p4||
    double dist1;
    double dist2;

    Eigen::Vector3d intersection_point;
};

// Clé de hachage pour l'indexation des paires cibles par distance
struct DistanceKey {
    int dist_bin;

    bool operator==(const DistanceKey& other) const {
        return dist_bin == other.dist_bin;
    }

    struct Hash {
        size_t operator()(const DistanceKey& k) const {
            return std::hash<int>()(k.dist_bin);
        }
    };
};

struct IndexedPair {
    const patch_planaire* p1;
    const patch_planaire* p2;
    double distance;
};

class v4pcs_algo {
public:
    // delta_dist : Tolérance de distance
    // delta_angle_deg : Tolérance d'angle normale
    v4pcs_algo(double delta_dist, double delta_angle_deg);

    Eigen::Matrix4d aligner(
        const std::vector<patch_planaire>& source_patches,
        const std::vector<patch_planaire>& target_patches,
        int max_bases = 1000,
        int max_seconds = 30
    );

private:
    double m_delta_dist;
    double m_delta_angle_cos;
    std::unordered_map<DistanceKey, std::vector<IndexedPair>, DistanceKey::Hash> m_target_index;

    // Etape 1 : Indexation des paires de la cible
    void construire_index_cible(const std::vector<patch_planaire>& patches);

    // Etape 2 : Sélectionner une base "large" et valide dans la source
    bool selectionner_base_aleatoire(
        const std::vector<patch_planaire>& patches,
        Base4PCS& out_base) const;

    // Etape 3 : Trouver les ensembles congruents dans la cible correspondant à la base
    std::vector<Eigen::Matrix4d> trouver_transformations_candidates(
        const Base4PCS& base) const;

    Eigen::Matrix4d calculer_transform_4points(
        const Base4PCS& base,
        const IndexedPair& pair1,
        const IndexedPair& pair2) const;

    double verifier_LCP(
        const Eigen::Matrix4d& T,
        const std::vector<patch_planaire>& source,
        const std::vector<patch_planaire>& target) const;
};

#endif