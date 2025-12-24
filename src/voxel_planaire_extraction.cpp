//
// Created by Enzo Gallet on 09/12/2025.
// Correction appliquée : Seuil de points réduit pour éviter de vider le nuage.
//

#include "voxel_planaire_extraction.h"
#include <Eigen/Eigenvalues>
#include <iostream> // Ajout pour debug si besoin

void voxel_planaire_extraction::extraire(const grille_voxel &gv)
{
    const auto& grille = gv.get_grille();

    m_grille_planaire.clear();
    m_grille_planaire.reserve(grille.size());

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;

    for (const auto& [clé, data] : grille)
    {
        // Correction :
        // avant y avais "const std::size_t minPts = 20;"
        // y avais un problème avec des voxels de 20cm, beaucoup de murs valides ont seulement 5-10 points.
        // 20 points supprime 80% des données. On passe à 4 c'est mieux mais à tester avec d'autre valeurs.
        const std::size_t minPts = 4;

        if (data.compteur < minPts)
            continue;

        const Eigen::Vector3d& barycentre = data.barycentre;
        const Eigen::Matrix3d& covariance = data.covariance;

        es.compute(covariance, Eigen::ComputeEigenvectors);
        if (es.info() != Eigen::Success)
            continue;

        const auto& evals = es.eigenvalues();
        const auto& evecs = es.eigenvectors();
        const double l0  = evals[0];
        const double l1  = evals[1];
        const double l2  = evals[2];
        const double sum = l0 + l1 + l2;

        if (sum <= 1e-9)
            continue;

        const double curvature = l0 / sum;
        const double linearity = (l2 - l1) / l2;

        if (curvature > m_curvature_thresh)
            continue;

        Eigen::Vector3d normal = evecs.col(0);

        voxel_planaire pv;
        pv.barycentre = barycentre;
        pv.normale    = normal;
        pv.linearity  = linearity;
        pv.curvature  = curvature;
        pv.count      = static_cast<std::size_t>(data.compteur);

        m_grille_planaire.emplace(clé, pv);
    }
}