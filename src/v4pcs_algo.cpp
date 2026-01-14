//
// Created by Fabien WAHBEH on 22/12/2025.
// Mise à jour : VERSION PERMUTATIONS (Teste tous les sens d'alignement)
//

#include "v4pcs_algo.h"
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <Eigen/Geometry>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

v4pcs_algo::v4pcs_algo(double delta_dist, double delta_angle_deg)
    : m_delta_dist(delta_dist)
{
    double rad = delta_angle_deg * M_PI / 180.0;
    m_delta_angle_cos = std::cos(rad);
}

void v4pcs_algo::construire_index_cible(const std::vector<patch_planaire>& patches) {
    m_target_index.clear();
    const size_t n = patches.size();

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double d = (patches[i].center - patches[j].center).norm();
            if (d < 0.5) continue;

            IndexedPair pair;
            pair.p1 = &patches[i];
            pair.p2 = &patches[j];
            pair.distance = d;

            DistanceKey key;
            key.dist_bin = static_cast<int>(d / m_delta_dist);
            m_target_index[key].push_back(pair);
        }
    }
}

static bool intersection_segments(
    const Eigen::Vector3d& p1, const Eigen::Vector3d& p2,
    const Eigen::Vector3d& q1, const Eigen::Vector3d& q2,
    Eigen::Vector3d& out_intersection,
    double& out_r1, double& out_r2)
{
    Eigen::Vector3d u = p2 - p1;
    Eigen::Vector3d v = q2 - q1;
    Eigen::Vector3d w = p1 - q1;

    double a = u.dot(u);
    double b = u.dot(v);
    double c = v.dot(v);
    double d = u.dot(w);
    double e = v.dot(w);
    double D = a * c - b * b;

    if (D < 1e-8) return false;

    double sc = (b * e - c * d) / D;
    double tc = (a * e - b * d) / D;

    Eigen::Vector3d intersect_p = p1 + sc * u;
    Eigen::Vector3d intersect_q = q1 + tc * v;

    if ((intersect_p - intersect_q).squaredNorm() > 0.5 * 0.5) return false;
    if (sc < 0.1 || sc > 0.9 || tc < 0.1 || tc > 0.9) return false;

    out_intersection = (intersect_p + intersect_q) * 0.5;
    out_r1 = sc;
    out_r2 = tc;
    return true;
}

bool v4pcs_algo::selectionner_base_aleatoire(
    const std::vector<patch_planaire>& patches,
    Base4PCS& out_base) const
{
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist(0, patches.size() - 1);

    for (int k = 0; k < 500; ++k) {
        const auto& p1 = patches[dist(rng)];
        const auto& p2 = patches[dist(rng)];
        const auto& p3 = patches[dist(rng)];
        const auto& p4 = patches[dist(rng)];

        if (&p1 == &p2 || &p3 == &p4) continue;

        double d1 = (p1.center - p2.center).norm();
        double d2 = (p3.center - p4.center).norm();

        // On garde une base assez large pour la stabilité
        if (d1 < 1.0 || d2 < 1.0) continue;

        Eigen::Vector3d e;
        double r1, r2;
        if (intersection_segments(p1.center, p2.center, p3.center, p4.center, e, r1, r2)) {
            out_base.p1 = &p1; out_base.p2 = &p2;
            out_base.p3 = &p3; out_base.p4 = &p4;
            out_base.r1 = r1;
            out_base.r2 = r2;
            out_base.dist1 = d1;
            out_base.dist2 = d2;
            out_base.intersection_point = e;
            return true;
        }
    }
    return false;
}

static Eigen::Matrix4d solve_umeyama(
    const std::vector<Eigen::Vector3d>& src_pts,
    const std::vector<Eigen::Vector3d>& tgt_pts)
{
    Eigen::Matrix<double, 3, 4> src, tgt;
    for(int i=0; i<4; ++i) {
        src.col(i) = src_pts[i];
        tgt.col(i) = tgt_pts[i];
    }
    return Eigen::umeyama(src, tgt, false);
}

std::vector<Eigen::Matrix4d> v4pcs_algo::trouver_transformations_candidates(
    const Base4PCS& base) const
{
    std::vector<Eigen::Matrix4d> transforms;
    transforms.reserve(100);

    DistanceKey key1{static_cast<int>(base.dist1 / m_delta_dist)};
    DistanceKey key2{static_cast<int>(base.dist2 / m_delta_dist)};

    std::vector<const IndexedPair*> cands1, cands2;

    // Récupération large
    for(int i=-1; i<=1; ++i) {
        DistanceKey k = key1; k.dist_bin += i;
        auto it = m_target_index.find(k);
        if(it != m_target_index.end()) {
            for(const auto& p : it->second) {
                if(std::abs(p.distance - base.dist1) < m_delta_dist) cands1.push_back(&p);
            }
        }
    }
    for(int i=-1; i<=1; ++i) {
        DistanceKey k = key2; k.dist_bin += i;
        auto it = m_target_index.find(k);
        if(it != m_target_index.end()) {
            for(const auto& p : it->second) {
                if(std::abs(p.distance - base.dist2) < m_delta_dist) cands2.push_back(&p);
            }
        }
    }

    // Sampling pour éviter l'explosion
    static std::mt19937 rng(std::random_device{}());
    if (cands1.size() > 50) { // Réduit pour tester les permutations
        std::shuffle(cands1.begin(), cands1.end(), rng);
        cands1.resize(50);
    }
    if (cands2.size() > 50) {
        std::shuffle(cands2.begin(), cands2.end(), rng);
        cands2.resize(50);
    }

    double sq_delta = m_delta_dist * m_delta_dist;

    std::vector<Eigen::Vector3d> src_pts = {
        base.p1->center, base.p2->center, base.p3->center, base.p4->center
    };

    for(const auto* pair1 : cands1) {
        for(const auto* pair2 : cands2) {

            // On construit les 2 variantes pour la première paire
            // Variante A: (p1->p1, p2->p2)
            // Variante B: (p1->p2, p2->p1)
            struct Config { Eigen::Vector3d a; Eigen::Vector3d b; };
            Config cfg1_A = { pair1->p1->center, pair1->p2->center };
            Config cfg1_B = { pair1->p2->center, pair1->p1->center };

            Config cfg2_A = { pair2->p1->center, pair2->p2->center };
            Config cfg2_B = { pair2->p2->center, pair2->p1->center };

            // On teste les 4 combinaisons de permutations (AA, AB, BA, BB)
            Config configs1[2] = {cfg1_A, cfg1_B};
            Config configs2[2] = {cfg2_A, cfg2_B};

            for(int i=0; i<2; ++i) {
                for(int j=0; j<2; ++j) {
                    const auto& c1 = configs1[i];
                    const auto& c2 = configs2[j];

                    // Vérification de l'intersection pour cette permutation spécifique
                    // e1 = c1.a + r1 * (c1.b - c1.a)
                    Eigen::Vector3d vec1 = c1.b - c1.a;
                    Eigen::Vector3d e1 = c1.a + base.r1 * vec1;

                    Eigen::Vector3d vec2 = c2.b - c2.a;
                    Eigen::Vector3d e2 = c2.a + base.r2 * vec2;

                    // Si les points d'intersection ne matchent pas, cette permutation est fausse
                    if ((e1 - e2).squaredNorm() > sq_delta) continue;

                    std::vector<Eigen::Vector3d> tgt_pts = { c1.a, c1.b, c2.a, c2.b };
                    transforms.push_back(solve_umeyama(src_pts, tgt_pts));
                }
            }

            if (transforms.size() >= 50) return transforms;
        }
    }
    return transforms;
}

Eigen::Matrix4d v4pcs_algo::calculer_transform_4points(
    const Base4PCS&, const IndexedPair&, const IndexedPair&) const
{
    return Eigen::Matrix4d::Identity();
}

double v4pcs_algo::verifier_LCP(
    const Eigen::Matrix4d& T,
    const std::vector<patch_planaire>& source,
    const std::vector<patch_planaire>& target) const
{
    double score = 0;
    // Vérification rapide (1 point sur 20)
    int step = std::max(1, (int)source.size() / 20);

    double sq_dist = m_delta_dist * m_delta_dist;

    for(size_t i = 0; i < source.size(); i += step) {
        Eigen::Vector3d p_world = (T * source[i].center.homogeneous()).head<3>();

        // On vérifie juste la distance pour commencer (plus robuste)
        for(const auto& t : target) {
            if((t.center - p_world).squaredNorm() < sq_dist) {
                score += 1.0;
                break;
            }
        }
    }
    return score;
}

Eigen::Matrix4d v4pcs_algo::aligner(
    const std::vector<patch_planaire>& source_patches,
    const std::vector<patch_planaire>& target_patches,
    int max_bases,
    int max_seconds)
{
    std::cout << "[V4PCS] Indexation cible..." << std::endl;
    construire_index_cible(target_patches);

    std::cout << "[V4PCS] Recherche RANSAC avec Permutations..." << std::endl;

    Eigen::Matrix4d best_T = Eigen::Matrix4d::Identity();
    double best_score = -1.0;

    auto start_time = std::chrono::high_resolution_clock::now();
    int valid_bases = 0;

    for(int i = 0; i < max_bases; ++i) {
        auto now = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(now - start_time).count();
        if (elapsed > max_seconds) {
            std::cout << "\n[TIMEOUT] Arret." << std::endl;
            break;
        }

        Base4PCS base;
        if(!selectionner_base_aleatoire(source_patches, base)) continue;

        auto candidates = trouver_transformations_candidates(base);
        if(!candidates.empty()) valid_bases++;

        if (i % 200 == 0) {
            std::cout << " Iter " << i << " | Valid: " << valid_bases
                      << " | Best: " << best_score << " | T=" << (int)elapsed << "s \r" << std::flush;
        }

        for(const auto& T : candidates) {
            double score = verifier_LCP(T, source_patches, target_patches);
            if(score > best_score) {
                best_score = score;
                best_T = T;
                std::cout << " * New Best Score: " << score << " (Iter " << i << ")" << std::endl;
            }
        }
    }
    std::cout << "\n";
    return best_T;
}