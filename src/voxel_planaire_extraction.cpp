//
// Created by Enzo Gallet on 09/12/2025.
// Correction appliquée : Seuil de points réduit pour éviter de vider le nuage.
//

#include "voxel_planaire_extraction.h"
#include <Eigen/Eigenvalues>

void voxel_planaire_extraction::extraire(const grille_voxel& gv)
{
    const auto& grille = gv.get_grille();

    m_grille_planaire.clear();
    m_grille_planaire.reserve(grille.size());

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;

    for (const auto& [cle, data] : grille)
    {
        if (data.compteur < static_cast<std::uint64_t>(m_min_points))
            continue;

        es.compute(data.covariance, Eigen::ComputeEigenvectors);
        if (es.info() != Eigen::Success)
            continue;

        const auto& evals = es.eigenvalues();
        const auto& evecs = es.eigenvectors();

        const double l0  = evals[0];
        const double l1  = evals[1];
        const double l2  = evals[2];
        const double sum = l0 + l1 + l2;

        if (sum <= 1e-12 || l2 <= 1e-12)
            continue;

        const double curvature = l0 / sum;
        const double linearity = (l2 - l1) / l2;
        const double planarity = (l1 - l0) / l2;

        if (curvature > m_curvature_max)  continue;
        if (linearity > m_linearity_max)  continue;
        if (planarity < m_planarity_min)  continue;

        Eigen::Vector3d normal = evecs.col(0);
        if (normal.squaredNorm() < 1e-12)
            continue;
        normal.normalize();

        voxel_planaire pv;
        pv.barycentre = data.barycentre;
        pv.normale    = normal;
        pv.linearity  = linearity;
        pv.planarity  = planarity;
        pv.curvature  = curvature;
        pv.count      = static_cast<std::size_t>(data.compteur);

        m_grille_planaire.emplace(cle, pv);
    }
}