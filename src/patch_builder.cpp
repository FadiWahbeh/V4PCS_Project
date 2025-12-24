//
// Created by Enzo Gallet on 09/12/2025.
// Correction appliquée : Nommage des propriétés PLY (r->red) pour compatibilité CloudCompare.
//

#include "patch_builder.h"

#include <vector>
#include <fstream>
#include <queue>
#include <random>
#include <unordered_set>

static void generer_couleur_aleatoire(std::uint8_t &r, std::uint8_t &g, std::uint8_t &b) {
    static std::mt19937 rng(std::random_device{}());
    static std::uniform_int_distribution<int> dist(50, 255); // Eviter les couleurs trop sombres
    r = static_cast<std::uint8_t>(dist(rng));
    g = static_cast<std::uint8_t>(dist(rng));
    b = static_cast<std::uint8_t>(dist(rng));
}

void patch_builder::construire_patches(
    const std::unordered_map<voxel_clé, voxel_planaire, voxel_clé::Hash>& grille_planaire,
    const std::unordered_map<voxel_clé, voxel_data, voxel_clé::Hash>& grille)
{
    m_patches.clear();
    if (grille_planaire.empty())
        return;

    std::vector<Eigen::Vector3i> neighbors_offsets;
    neighbors_offsets.reserve(26);
    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz)
                if (dx != 0 || dy != 0 || dz != 0)
                    neighbors_offsets.emplace_back(dx, dy, dz);

    std::unordered_set<voxel_clé, voxel_clé::Hash> visited;
    visited.reserve(grille_planaire.size());

    for (const auto& [start_key, start_voxel] : grille_planaire)
    {
        if (visited.contains(start_key))
            continue;

        patch_planaire patch;
        generer_couleur_aleatoire(patch.r, patch.g, patch.b);
        std::queue<voxel_clé> file;
        file.push(start_key);
        visited.insert(start_key);

        Eigen::Vector3d sum_normals = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_centers = Eigen::Vector3d::Zero();
        Eigen::Vector3d ref_normal = start_voxel.normale;

        while (!file.empty())
        {
            voxel_clé current_key = file.front();
            file.pop();

            const auto& current_voxel = grille_planaire.at(current_key);
            patch.voxels.push_back(&current_voxel);

            Eigen::Vector3d cur_n = current_voxel.normale;
            if (cur_n.dot(ref_normal) < 0) cur_n = -cur_n;

            sum_centers += current_voxel.barycentre;
            sum_normals += cur_n;

            for (const auto& offset : neighbors_offsets)
            {
                voxel_clé neighbor_key;
                neighbor_key.index = current_key.index + offset;

                if (visited.contains(neighbor_key)) continue;
                if (!grille_planaire.contains(neighbor_key)) continue;

                const auto& neighbor_voxel = grille_planaire.at(neighbor_key);

                if (est_connexion_lisse(current_voxel, neighbor_voxel))
                {
                    visited.insert(neighbor_key);
                    file.push(neighbor_key);
                }
            }
        }

        if (!patch.voxels.empty())
        {
            double n = static_cast<double>(patch.voxels.size());
            patch.center = sum_centers / n;
            patch.normal = sum_normals.normalized();
            m_patches.push_back(std::move(patch));
        }
    }
}

bool patch_builder::est_connexion_lisse(const voxel_planaire &a, const voxel_planaire &b) const {
    double dot = std::abs(a.normale.dot(b.normale));
    if (dot < m_angle_thresh_cos)
        return false;

    Eigen::Vector3d vec = b.barycentre - a.barycentre;
    double dist = std::abs(a.normale.dot(vec));

    if (dist > m_distance_thresh)
        return false;

    return true;
}

void patch_builder::exporter_patches_ply(const std::filesystem::path& nom_fichier) const
{
    {
        std::filesystem::path patches_path = nom_fichier.string() + ".patch_centers.ply";
        std::ofstream ofs(patches_path);
        if (!ofs) throw std::runtime_error("Erreur ouverture : " + patches_path.string());

        ofs << "ply\n"
               "format ascii 1.0\n"
               "element vertex " << m_patches.size() << "\n"
               "property float x\n"
               "property float y\n"
               "property float z\n"
               "property uchar red\n"
               "property uchar green\n"
               "property uchar blue\n"
               "end_header\n";

        for (const auto& p : m_patches) {
            ofs << p.center.x() << " " << p.center.y() << " " << p.center.z() << " "
                << (int)p.r << " " << (int)p.g << " " << (int)p.b << "\n";
        }
    }

    {
        std::filesystem::path voxels_path = nom_fichier.string() + ".patch_voxels.ply";
        std::ofstream ofs(voxels_path);
        if (!ofs) throw std::runtime_error("Erreur ouverture : " + voxels_path.string());

        std::size_t n_voxels = 0;
        for (const auto& p : m_patches) n_voxels += p.voxels.size();

        ofs << "ply\n"
               "format ascii 1.0\n"
               "element vertex " << n_voxels << "\n"
               "property float x\n"
               "property float y\n"
               "property float z\n"
               "property uchar red\n"    // <--- CORRECTION: red au lieu de r
               "property uchar green\n"  // <--- CORRECTION: green au lieu de g
               "property uchar blue\n"   // <--- CORRECTION: blue au lieu de b
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