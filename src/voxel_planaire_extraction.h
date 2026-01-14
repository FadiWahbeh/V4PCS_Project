//
// Created by Enzo Gallet on 09/12/2025.
//

#ifndef VOXEL_PLANAIRE_EXTRACTION_H
#define VOXEL_PLANAIRE_EXTRACTION_H

#include <cstddef>
#include <unordered_map>
#include <Eigen/Core>

#include "grille_voxel.h"

// PCA eigenvalues λ0 ≤ λ1 ≤ λ2
// courbure  = λ0 / (λ0+λ1+λ2)  -> petit pour un plan
// linéaire  = (λ2-λ1) / λ2     -> grand pour une ligne
// planaire  = (λ1-λ0) / λ2     -> grand pour un plan
struct voxel_planaire
{
    Eigen::Vector3d barycentre;
    Eigen::Vector3d normale;
    double          linearity;
    double          planarity;
    double          curvature;
    std::size_t     count;
};

class voxel_planaire_extraction {
public:
    // linearity_max : max autorisé rejette les formes "ligne"
    // curvature_max : max autorisé rejette voxels trop bruités / non-planaires
    // planarity_min : min exigé garde seulement les vrais plans
    // min_points    : nb minimal de points pour PCA stable
    voxel_planaire_extraction(double linearity_max,
                              double curvature_max,
                              double planarity_min = 0.40,
                              std::size_t min_points = 10)
        : m_linearity_max(linearity_max),
          m_curvature_max(curvature_max),
          m_planarity_min(planarity_min),
          m_min_points(min_points)
    {}

    void extraire(const grille_voxel& gv);

    const std::unordered_map<voxel_clé, voxel_planaire, voxel_clé::Hash>& get_grille_planaire() const noexcept
    {
        return m_grille_planaire;
    }

private:
    double m_linearity_max;
    double m_curvature_max;
    double m_planarity_min;
    std::size_t m_min_points;

    std::unordered_map<voxel_clé, voxel_planaire, voxel_clé::Hash> m_grille_planaire;
};

#endif