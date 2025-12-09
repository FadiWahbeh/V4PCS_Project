//
// Created by Enzo Gallet on 09/12/2025.
//

#ifndef VOXEL_PLANAIRE_EXTRACTION_H
#define VOXEL_PLANAIRE_EXTRACTION_H
#include "grille_voxel.h"
#include <Eigen/Core>

struct voxel_planaire
{
    //voxel_clé       key;        // indice (ix, iy, iz)
    Eigen::Vector3d barycentre; // barycentre des points du voxel
    Eigen::Vector3d normale;    // normale estimée
    double          linearity;  // mesure de "planarité"
    double          curvature;  // idem
    std::size_t     count;      // nb de points
};


class voxel_planaire_extraction {
public:
    voxel_planaire_extraction(const double linearity_thresh,
                              const double curvature_thresh)
        : m_linearity_thresh(linearity_thresh),
          m_curvature_thresh(curvature_thresh)
    {};

    void extraire(const grille_voxel &gv);

    const std::unordered_map<voxel_clé, voxel_planaire, voxel_clé::Hash>& get_grille_planaire() const noexcept
    {
        return m_grille_planaire;
    }


private:
    double m_linearity_thresh;
    double m_curvature_thresh;

    std::unordered_map<voxel_clé, voxel_planaire, voxel_clé::Hash> m_grille_planaire;

};



#endif //VOXEL_PLANAIRE_EXTRACTION_H
