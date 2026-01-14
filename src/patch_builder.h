//
// Created by Enzo Gallet on 09/12/2025.
//

#ifndef PATCH_BUILDER_H
#define PATCH_BUILDER_H

#include <cstdint>
#include <filesystem>
#include <vector>
#include <unordered_map>
#include <Eigen/Core>

#include "voxel_planaire_extraction.h"

// il ne faut pas stocker des pointeurs vers des éléments d'une unordered_map.
//  il faut plutot stocker les clés plus une copie des barycentres. 'plus rapide'
struct patch_planaire
{
    Eigen::Vector3d center;
    Eigen::Vector3d normal;

    std::vector<voxel_clé> voxel_keys;
    std::vector<Eigen::Vector3d> voxel_centers;

    std::uint8_t r = 200, g = 200, b = 200;
};

class patch_builder {
public:
    patch_builder(double angle_thresh_deg, double distance_thresh);

    void construire_patches(
        const std::unordered_map<voxel_clé, voxel_planaire, voxel_clé::Hash>& grille_planaire,
        const std::unordered_map<voxel_clé, voxel_data, voxel_clé::Hash>& grille);

    void exporter_patches_ply(const std::filesystem::path& nom_fichier) const;

    int get_patches_taille() const { return static_cast<int>(m_patches.size()); }

    const std::vector<patch_planaire>& get_patches() const { return m_patches; }

private:
    double m_angle_thresh_cos;
    double m_distance_thresh;

    std::vector<patch_planaire> m_patches;

    bool est_connexion_lisse_patch(const Eigen::Vector3d& patch_center,
                                  const Eigen::Vector3d& patch_normal,
                                  const voxel_planaire& candidate) const;

    void refit_patch_plane_and_prune(patch_planaire& patch) const;
};

#endif