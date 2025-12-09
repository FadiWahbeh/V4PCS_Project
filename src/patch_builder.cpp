//
// Created by Enzo Gallet on 09/12/2025.
//

#include "patch_builder.h"

#include <vector>

#include <fstream>
#include <queue>
#include <random>
#include <set>
#include <unordered_set>

void patch_builder::construire_patches(
    const std::unordered_map<voxel_clé, voxel_planaire, voxel_clé::Hash>& grille_planaire,
    const std::unordered_map<voxel_clé, voxel_data, voxel_clé::Hash>& grille)
{
    m_patches.clear();
    if (grille_planaire.empty())
        return;

    // 26-voisinage
    std::vector<Eigen::Vector3i> neighbors_offsets;
    neighbors_offsets.reserve(26);
    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz)
                if (dx != 0 || dy != 0 || dz != 0)
                    neighbors_offsets.emplace_back(dx, dy, dz);

    // Ensemble des clés déjà visitées
    std::unordered_set<voxel_clé, voxel_clé::Hash> visited;
    visited.reserve(grille_planaire.size());

    // Générateur de couleurs
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> color_dist(50, 255);

    // Parcours de tous les voxels planaires
    for (const auto& kv : grille_planaire)
    {
        const voxel_clé& start_key = kv.first;

        if (visited.find(start_key) != visited.end())
            continue;

        patch_planaire patch;
        patch.center.setZero();
        patch.normal.setZero();
        patch.voxels.clear();

        patch.r = static_cast<std::uint8_t>(color_dist(rng));
        patch.g = static_cast<std::uint8_t>(color_dist(rng));
        patch.b = static_cast<std::uint8_t>(color_dist(rng));

        Eigen::Vector3d sum_centers = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_normals = Eigen::Vector3d::Zero();

        std::queue<voxel_clé> q;
        q.push(start_key);
        visited.insert(start_key);

        while (!q.empty())
        {
            voxel_clé current_key = q.front();
            q.pop();

            auto it_current = grille_planaire.find(current_key);
            if (it_current == grille_planaire.end())
                continue; // ne devrait pas arriver, par sécurité

            const voxel_planaire& current_voxel = it_current->second;

            // Ajouter ce voxel au patch
            patch.voxels.push_back(&current_voxel);

            // Accumulation du centre
            sum_centers += current_voxel.barycentre;

            // Accumulation de la normale en gérant l'orientation
            if (sum_normals.isZero(1e-8))
            {
                sum_normals += current_voxel.normale;
            }
            else
            {
                if (current_voxel.normale.dot(sum_normals) < 0.0)
                    sum_normals -= current_voxel.normale;
                else
                    sum_normals += current_voxel.normale;
            }

            // Exploration des voisins dans la grille
            const Eigen::Vector3i& idx = current_key.index;
            for (const auto& off : neighbors_offsets)
            {
                auto res = idx + off;
                voxel_clé neighbor_key;
                neighbor_key.index = idx + off;


                if (visited.find(neighbor_key) != visited.end())
                    continue;

                auto it_neighbor = grille_planaire.find(neighbor_key);
                if (it_neighbor == grille_planaire.end())
                    continue; // ce voxel n'existe pas dans la grille planaire

                const voxel_planaire& neighbor_voxel = it_neighbor->second;

                if (est_connexion_lisse(current_voxel, neighbor_voxel))
                {
                    visited.insert(neighbor_key);
                    q.push(neighbor_key);
                }
            }
        }

        if (!patch.voxels.empty())
        {
            const double n = static_cast<double>(patch.voxels.size());
            patch.center = sum_centers / n;

            if (!sum_normals.isZero(1e-8))
                patch.normal = sum_normals.normalized();
            else
                patch.normal = patch.voxels.front()->normale;

            m_patches.push_back(std::move(patch));
        }
    }


}


void patch_builder::exporter_patches_ply(const std::filesystem::path& nom_fichier) const
{
    // Construction des chemins de sortie
    const std::filesystem::path centers_path =
        std::filesystem::path(DATA_DIR_OUTPUT) /
        (nom_fichier.filename().string() + ".patch_centers.ply");

    const std::filesystem::path voxels_path =
        std::filesystem::path(DATA_DIR_OUTPUT) /
        (nom_fichier.filename().string() + ".patch_voxels.ply");

    // ============================
    // Export des centres
    // ============================

    {
        std::ofstream ofs(centers_path, std::ios::binary);
        if (!ofs)
            throw std::runtime_error("Impossible d'ouvrir le fichier : " + centers_path.string());

        const std::size_t n_points = m_patches.size();

        ofs << "ply\n"
               "format ascii 1.0\n"
               "element vertex " << n_points << "\n"
               "property float x\n"
               "property float y\n"
               "property float z\n"
               "property uchar r\n"
               "property uchar g\n"
               "property uchar b\n"
               "end_header\n";

        for (const patch_planaire& p : m_patches)
        {
            ofs << static_cast<float>(p.center.x()) << " "
                << static_cast<float>(p.center.y()) << " "
                << static_cast<float>(p.center.z()) << " "
                << static_cast<int>(p.r) << " "
                << static_cast<int>(p.g) << " "
                << static_cast<int>(p.b) << "\n";
        }
    }

    // ============================
    // Export des voxels
    // ============================

    {
        // Comptage du nombre total de voxels
        std::size_t n_voxels = 0;
        for (const patch_planaire& p : m_patches)
            n_voxels += p.voxels.size();

        std::ofstream ofs(voxels_path, std::ios::binary);
        if (!ofs)
            throw std::runtime_error("Impossible d'ouvrir le fichier : " + voxels_path.string());

        ofs << "ply\n"
               "format ascii 1.0\n"
               "element vertex " << n_voxels << "\n"
               "property float x\n"
               "property float y\n"
               "property float z\n"
               "property uchar r\n"
               "property uchar g\n"
               "property uchar b\n"
               "end_header\n";

        for (const patch_planaire& p : m_patches)
        {
            for (const voxel_planaire* v : p.voxels)
            {
                ofs << static_cast<float>(v->barycentre.x()) << " "
                    << static_cast<float>(v->barycentre.y()) << " "
                    << static_cast<float>(v->barycentre.z()) << " "
                    << static_cast<int>(p.r) << " "
                    << static_cast<int>(p.g) << " "
                    << static_cast<int>(p.b) << "\n";
            }
        }
    }
}


bool patch_builder::est_connexion_lisse(const voxel_planaire &a, const voxel_planaire &b) const {
    const double dot = std::abs(a.normale.dot(b.normale));
    if (dot < m_angle_thresh_cos) return false;

    const double dist1 = std::abs(a.normale.dot(a.barycentre - b.barycentre));
    const double dist2 = std::abs(b.normale.dot(b.barycentre - a.barycentre));

    if (dist1 > m_distance_thresh || dist2 > m_distance_thresh) return false;
    return true;
}

