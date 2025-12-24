//
// Created by Fabien WAHBEH on 22/12/2025..
//

#ifndef V4PCS_ALGO_H
#define V4PCS_ALGO_H

#include <vector>
#include <unordered_map>
#include <Eigen/Core>
#include "patch_pair_generator.h"
#include "patch_builder.h"

struct PairKey {
    int dist_bin;
    int angle_bin;

    bool operator==(const PairKey& other) const {
        return dist_bin == other.dist_bin && angle_bin == other.angle_bin;
    }

    struct Hash {
        size_t operator()(const PairKey& k) const {
            return std::hash<int>()(k.dist_bin) ^ (std::hash<int>()(k.angle_bin) << 1);
        }
    };
};

class v4pcs_algo {
public:
    v4pcs_algo(double delta_dist, double delta_angle_deg);

    Eigen::Matrix4d aligner(
        const std::vector<patch_pair>& source_pairs,
        const std::vector<patch_pair>& target_pairs,
        const std::vector<patch_planaire>& source_patches,
        const std::vector<patch_planaire>& target_patches,
        int max_iterations = 1000,
        double overlap_dist = 0.20
    );

private:
    double m_delta_dist;
    double m_delta_angle_cos;
    double m_delta_angle_deg;

    std::unordered_map<PairKey, std::vector<const patch_pair*>, PairKey::Hash> m_target_index;

    void construire_index(const std::vector<patch_pair>& target_pairs);
    std::vector<const patch_pair*> trouver_candidats(const patch_pair& src_pair) const;

    Eigen::Matrix4d estimer_transform(const patch_pair& src, const patch_pair& tgt) const;

    double calculer_score_LCP(
        const Eigen::Matrix4d& transform,
        const std::vector<patch_planaire>& source,
        const std::vector<patch_planaire>& target,
        double dist_thresh) const;
};

#endif