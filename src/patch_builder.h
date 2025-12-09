//
// Created by Enzo Gallet on 09/12/2025.
//

#ifndef PATCH_BUILDER_H
#define PATCH_BUILDER_H
#include <vector>

#include <Eigen/Core>

#include "voxel_planaire_extraction.h"

struct patch_planaire
{
    Eigen::Vector3d center;                    // centre du patch
    Eigen::Vector3d normal;                    // normale moyenne
    std::vector<const voxel_planaire*> voxels; // pointeurs vers les voxels planaires membres

    std::uint8_t r, g, b; // couleur pour visualisation
};

class patch_builder {
public:
    patch_builder(const double angle_thresh_deg, const double distance_thresh)
        : m_angle_thresh_cos(std::cos(angle_thresh_deg * M_PI / 180.0)),
          m_distance_thresh(distance_thresh)
    {}

    void construire_patches(const std::unordered_map<voxel_clé, voxel_planaire, voxel_clé::Hash>& grille_planaire,
                            const std::unordered_map<voxel_clé, voxel_data, voxel_clé::Hash>& grille);

    void exporter_patches_ply(const std::filesystem::path& nom_fichier) const;

    int get_patches_taille() const {return m_patches.size();}


private:
    double m_angle_thresh_cos;
    double m_distance_thresh;

    std::vector<patch_planaire> m_patches;

    bool sont_voisins(const voxel_planaire& a, const voxel_planaire& b) const;

    bool est_connexion_lisse(const voxel_planaire& a, const voxel_planaire& b) const;

};



#endif //PATCH_BUILDER_H
