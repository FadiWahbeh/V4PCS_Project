//
// Created by Enzo Gallet on 09/12/2025.
//

#include "voxel_planaire_extraction.h"
#include <Eigen/Eigenvalues>


void voxel_planaire_extraction::extraire(const grille_voxel &gv)
{
    const auto& grille = gv.get_grille();

    m_grille_planaire.clear();
    m_grille_planaire.reserve(grille.size());

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;

    for (const auto& [clé, data] : grille)
    {
        if (data.compteur < 3)
            continue;

        const Eigen::Vector3d& barycentre = data.barycentre;
        const Eigen::Matrix3d& covariance = data.covariance;

        es.compute(covariance, Eigen::ComputeEigenvectors);
        if (es.info() != Eigen::Success)
            continue;

        const auto& evals = es.eigenvalues();
        const auto& evecs = es.eigenvectors();

        const double l0  = evals[0]; // plus petite
        const double l1  = evals[1];
        const double l2  = evals[2]; // plus grande
        const double sum = l0 + l1 + l2;

        if (sum <= 0.0)
            continue;

        const double curvature = l0 / sum;
        const double linearity = (l2 - l1) / l2;

        // Critères (adaptables selon besoin)
        if (curvature > m_curvature_thresh)
            continue;
        if (linearity < m_linearity_thresh)
            continue;

        Eigen::Vector3d normal = evecs.col(0); // déjà normalisé

        voxel_planaire pv;
        pv.barycentre = barycentre;
        pv.normale    = normal;
        pv.linearity  = linearity;
        pv.curvature  = curvature;
        pv.count      = static_cast<std::size_t>(data.compteur);

        m_grille_planaire.emplace(clé, pv);
    }
}


