//
// Created by Fabien WAHBEH on 22/12/2025.
//

#include "v4pcs_algo.h"
#include <iostream>
#include <random>
#include <cmath>
#include <Eigen/SVD>
#include <Eigen/Dense>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

v4pcs_algo::v4pcs_algo(double delta_dist, double delta_angle_deg)
    : m_delta_dist(delta_dist), m_delta_angle_deg(delta_angle_deg)
{
    m_delta_angle_cos = std::cos(delta_angle_deg * M_PI / 180.0);
}

void v4pcs_algo::construire_index(const std::vector<patch_pair>& target_pairs) {
    m_target_index.clear();
    for(const auto& p : target_pairs) {
        int d_bin = static_cast<int>(p.distance / m_delta_dist);

        double angle = std::acos(std::min(1.0, std::abs(p.angle_dot))) * 180.0 / M_PI;
        int a_bin = static_cast<int>(angle / m_delta_angle_deg);

        m_target_index[{d_bin, a_bin}].push_back(&p);
    }
}

std::vector<const patch_pair*> v4pcs_algo::trouver_candidats(const patch_pair& src) const {
    int d_bin = static_cast<int>(src.distance / m_delta_dist);
    double angle = std::acos(std::min(1.0, std::abs(src.angle_dot))) * 180.0 / M_PI;
    int a_bin = static_cast<int>(angle / m_delta_angle_deg);

    std::vector<const patch_pair*> candidats;

    for(int dd = -1; dd <= 1; ++dd) {
        for(int da = -1; da <= 1; ++da) {
            auto it = m_target_index.find({d_bin + dd, a_bin + da});
            if(it != m_target_index.end()) {
                for(const auto* p_tgt : it->second) {
                    if(std::abs(p_tgt->distance - src.distance) < m_delta_dist &&
                       std::abs(p_tgt->angle_dot - src.angle_dot) < 0.1) { // 0.1 approx tolerance cos
                        candidats.push_back(p_tgt);
                    }
                }
            }
        }
    }
    return candidats;
}

Eigen::Matrix4d v4pcs_algo::estimer_transform(const patch_pair& src, const patch_pair& tgt) const {
    Eigen::Matrix<double, 3, 4> P, Q;

    P.col(0) = src.p1->center;
    P.col(1) = src.p2->center;
    P.col(2) = src.p1->center + src.p1->normal;
    P.col(3) = src.p2->center + src.p2->normal;

    Q.col(0) = tgt.p1->center;
    Q.col(1) = tgt.p2->center;
    Q.col(2) = tgt.p1->center + tgt.p1->normal;
    Q.col(3) = tgt.p2->center + tgt.p2->normal;

    return Eigen::umeyama(P, Q, false);
}

double v4pcs_algo::calculer_score_LCP(
    const Eigen::Matrix4d& T,
    const std::vector<patch_planaire>& source,
    const std::vector<patch_planaire>& target,
    double dist_thresh) const
{
    double score = 0;
    int step = 10;
    double sq_dist = dist_thresh * dist_thresh;

    for(size_t i = 0; i < source.size(); i += step) {
        Eigen::Vector3d p_trans = (T * source[i].center.homogeneous()).head<3>();
        Eigen::Vector3d n_trans = T.block<3,3>(0,0) * source[i].normal;

        for(const auto& t : target) {
            if((t.center - p_trans).squaredNorm() < sq_dist) {
                if(t.normal.dot(n_trans) > 0.8) {
                    score += 1.0;
                    break;
                }
            }
        }
    }
    return score;
}

Eigen::Matrix4d v4pcs_algo::aligner(
    const std::vector<patch_pair>& source_pairs,
    const std::vector<patch_pair>& target_pairs,
    const std::vector<patch_planaire>& source_patches,
    const std::vector<patch_planaire>& target_patches,
    int max_iterations,
    double overlap_dist)
{
    std::cout << "[V4PCS] Construction de l'index geometrique..." << std::endl;
    construire_index(target_pairs);
    std::cout << "[V4PCS] Index pret. Demarrage RANSAC." << std::endl;

    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist_src(0, source_pairs.size() - 1);

    Eigen::Matrix4d best_T = Eigen::Matrix4d::Identity();
    double best_score = -1.0;

    for(int i=0; i<max_iterations; ++i) {
        const auto& src_pair = source_pairs[dist_src(rng)];
        auto candidats = trouver_candidats(src_pair);

        for(const auto* tgt_pair : candidats) {
            Eigen::Matrix4d T = estimer_transform(src_pair, *tgt_pair);

            double score = calculer_score_LCP(T, source_patches, target_patches, overlap_dist);

            if(score > best_score) {
                best_score = score;
                best_T = T;
                std::cout << " Iter " << i << " | Score: " << score << " | Candidats: " << candidats.size() << std::endl;
            }
        }

        if(i % 1000 == 0) std::cout << "." << std::flush;
    }
    std::cout << "\n[V4PCS] Termine." << std::endl;
    return best_T;
}