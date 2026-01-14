//
// Created by Enzo Gallet on 09/12/2025.
//

#include "patch_builder.h"

#include <cmath>
#include <fstream>
#include <queue>
#include <random>
#include <stdexcept>
#include <unordered_set>

#include <Eigen/Eigenvalues>

static void generer_couleur_aleatoire(std::uint8_t& r, std::uint8_t& g, std::uint8_t& b) {
    static std::mt19937 rng(std::random_device{}());
    static std::uniform_int_distribution<int> dist(50, 255);
    r = static_cast<std::uint8_t>(dist(rng));
    g = static_cast<std::uint8_t>(dist(rng));
    b = static_cast<std::uint8_t>(dist(rng));
}

patch_builder::patch_builder(double angle_thresh_deg, double distance_thresh)
{
    const double rad = angle_thresh_deg * (3.14159265358979323846 / 180.0);
    m_angle_thresh_cos = std::cos(rad);
    m_distance_thresh  = distance_thresh;
}

bool patch_builder::est_connexion_lisse_patch(const Eigen::Vector3d& patch_center,
                                              const Eigen::Vector3d& patch_normal,
                                              const voxel_planaire& candidate) const
{
    const double dot = std::abs(patch_normal.dot(candidate.normale));
    if (dot < m_angle_thresh_cos) return false;

    const double dist = std::abs(patch_normal.dot(candidate.barycentre - patch_center));
    if (dist > m_distance_thresh) return false;

    return true;
}

void patch_builder::refit_patch_plane_and_prune(patch_planaire& patch) const
{
    const std::size_t n = patch.voxel_centers.size();
    if (n < 3) return;

    Eigen::Vector3d mean = Eigen::Vector3d::Zero();
    for (const auto& c : patch.voxel_centers) mean += c;
    mean /= static_cast<double>(n);

    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    for (const auto& c : patch.voxel_centers) {
        const Eigen::Vector3d d = c - mean;
        cov += d * d.transpose();
    }
    cov /= static_cast<double>(n);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
    if (es.info() != Eigen::Success) return;

    Eigen::Vector3d normal = es.eigenvectors().col(0);
    if (normal.squaredNorm() < 1e-12) return;
    normal.normalize();

    std::vector<voxel_clé> new_keys;
    std::vector<Eigen::Vector3d> new_centers;
    new_keys.reserve(n);
    new_centers.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
        const double dist = std::abs(normal.dot(patch.voxel_centers[i] - mean));
        if (dist <= m_distance_thresh) {
            new_keys.push_back(patch.voxel_keys[i]);
            new_centers.push_back(patch.voxel_centers[i]);
        }
    }

    if (new_centers.size() >= 3 && new_centers.size() < n) {
        patch.voxel_keys.swap(new_keys);
        patch.voxel_centers.swap(new_centers);

        const std::size_t n2 = patch.voxel_centers.size();
        Eigen::Vector3d mean2 = Eigen::Vector3d::Zero();
        for (const auto& c : patch.voxel_centers) mean2 += c;
        mean2 /= static_cast<double>(n2);

        Eigen::Matrix3d cov2 = Eigen::Matrix3d::Zero();
        for (const auto& c : patch.voxel_centers) {
            const Eigen::Vector3d d = c - mean2;
            cov2 += d * d.transpose();
        }
        cov2 /= static_cast<double>(n2);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es2(cov2);
        if (es2.info() == Eigen::Success) {
            Eigen::Vector3d n2v = es2.eigenvectors().col(0);
            if (n2v.squaredNorm() > 1e-12) {
                n2v.normalize();
                patch.center = mean2;
                patch.normal = n2v;
                return;
            }
        }

        patch.center = mean2;
        patch.normal = normal;
        return;
    }

    patch.center = mean;
    patch.normal = normal;
}

void patch_builder::construire_patches(
    const std::unordered_map<voxel_clé, voxel_planaire, voxel_clé::Hash>& grille_planaire,
    const std::unordered_map<voxel_clé, voxel_data, voxel_clé::Hash>& /*grille*/)
{
    m_patches.clear();
    if (grille_planaire.empty()) return;

    // 6-voisinage
    const std::vector<Eigen::Vector3i> neighbors_offsets = {
        { 1, 0, 0}, {-1, 0, 0},
        { 0, 1, 0}, { 0,-1, 0},
        { 0, 0, 1}, { 0, 0,-1}
    };

    std::unordered_set<voxel_clé, voxel_clé::Hash> visited;
    visited.reserve(grille_planaire.size());

    for (const auto& [start_key, start_voxel] : grille_planaire)
    {
        if (visited.contains(start_key)) continue;

        patch_planaire patch;
        generer_couleur_aleatoire(patch.r, patch.g, patch.b);

        std::queue<voxel_clé> q;
        q.push(start_key);
        visited.insert(start_key);

        Eigen::Vector3d sum_centers = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_normals = Eigen::Vector3d::Zero();
        std::size_t count = 0;

        Eigen::Vector3d ref_normal = start_voxel.normale;
        if (ref_normal.squaredNorm() < 1e-12) ref_normal = Eigen::Vector3d(0,0,1);
        ref_normal.normalize();

        while (!q.empty())
        {
            const voxel_clé current_key = q.front();
            q.pop();

            const auto it = grille_planaire.find(current_key);
            if (it == grille_planaire.end()) continue;

            const voxel_planaire& v = it->second;

            patch.voxel_keys.push_back(current_key);
            patch.voxel_centers.push_back(v.barycentre);

            Eigen::Vector3d cur_n = v.normale;
            if (cur_n.squaredNorm() < 1e-12) cur_n = ref_normal;
            else cur_n.normalize();
            if (cur_n.dot(ref_normal) < 0.0) cur_n = -cur_n;

            sum_centers += v.barycentre;
            sum_normals += cur_n;
            ++count;

            Eigen::Vector3d patch_center = sum_centers / static_cast<double>(count);
            Eigen::Vector3d patch_normal = sum_normals;
            if (patch_normal.squaredNorm() < 1e-12) patch_normal = ref_normal;
            else patch_normal.normalize();

            for (const auto& off : neighbors_offsets)
            {
                voxel_clé nk;
                nk.index = current_key.index + off;

                if (visited.contains(nk)) continue;

                const auto itn = grille_planaire.find(nk);
                if (itn == grille_planaire.end()) continue;

                if (est_connexion_lisse_patch(patch_center, patch_normal, itn->second)) {
                    visited.insert(nk);
                    q.push(nk);
                }
            }
        }

        if (!patch.voxel_centers.empty())
        {
            patch.center = sum_centers / static_cast<double>(patch.voxel_centers.size());
            patch.normal = sum_normals;
            if (patch.normal.squaredNorm() < 1e-12) patch.normal = ref_normal;
            else patch.normal.normalize();
            if (patch.normal.dot(ref_normal) < 0.0) patch.normal = -patch.normal;

            refit_patch_plane_and_prune(patch);
            if (patch.normal.dot(ref_normal) < 0.0) patch.normal = -patch.normal;

            m_patches.push_back(std::move(patch));
        }
    }
}

void patch_builder::exporter_patches_ply(const std::filesystem::path& nom_fichier) const
{
    // Patch centers
    {
        std::filesystem::path patches_path = nom_fichier.string() + ".patch_centers.ply";
        std::ofstream ofs(patches_path);
        if (!ofs) throw std::runtime_error("Erreur ouverture : " + patches_path.string());

        ofs << "ply\nformat ascii 1.0\n"
            << "element vertex " << m_patches.size() << "\n"
            << "property float x\nproperty float y\nproperty float z\n"
            << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
            << "end_header\n";

        for (const auto& p : m_patches) {
            ofs << (float)p.center.x() << " " << (float)p.center.y() << " " << (float)p.center.z() << " "
                << (int)p.r << " " << (int)p.g << " " << (int)p.b << "\n";
        }
    }

    // Patch voxels (centers)
    {
        std::filesystem::path voxels_path = nom_fichier.string() + ".patch_voxels.ply";
        std::ofstream ofs(voxels_path);
        if (!ofs) throw std::runtime_error("Erreur ouverture : " + voxels_path.string());

        std::size_t total = 0;
        for (const auto& p : m_patches) total += p.voxel_centers.size();

        ofs << "ply\nformat ascii 1.0\n"
            << "element vertex " << total << "\n"
            << "property float x\nproperty float y\nproperty float z\n"
            << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
            << "end_header\n";

        for (const auto& p : m_patches) {
            for (const auto& c : p.voxel_centers) {
                ofs << (float)c.x() << " " << (float)c.y() << " " << (float)c.z() << " "
                    << (int)p.r << " " << (int)p.g << " " << (int)p.b << "\n";
            }
        }
    }
}