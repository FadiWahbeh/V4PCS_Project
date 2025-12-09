//
// Created by Enzo Gallet on 09/12/2025.
//

#include "grille_voxel.h"

#include <fstream>
#include <iostream>

void grille_voxel::charger_points() {
    std::ifstream in(m_fichier);
    if (!in)
        throw std::runtime_error("Impossible d'ouvrir le fichier : " + m_fichier.string());

    m_points.clear();


    double x, y, z;
    while (in >> x >> y >> z)
    {
        m_points.emplace_back(x, y, z);
    }

    if (m_points.empty())
        throw std::runtime_error("Fichier vide : " + m_fichier.string());
}

void grille_voxel::estimer_voxel_size(const double sommets_par_voxel,
                                      const double min_size,
                                      const double max_size)
{
    double min_x =  std::numeric_limits<double>::infinity();
    double min_y =  std::numeric_limits<double>::infinity();
    double min_z =  std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();

    const std::size_t count = m_points.size();

    for (const auto& p : m_points)
    {
        const double x = p.x();
        const double y = p.y();
        const double z = p.z();

        if (x < min_x) min_x = x;
        if (y < min_y) min_y = y;
        if (z < min_z) min_z = z;
        if (x > max_x) max_x = x;
        if (y > max_y) max_y = y;
        if (z > max_z) max_z = z;
    }

    const double dx = max_x - min_x;
    const double dy = max_y - min_y;
    const double dz = max_z - min_z;

    double volume = dx * dy * dz;
    if (volume <= 0.0)
        volume = 1e-9;

    const double densité = static_cast<double>(count) / volume;

    double voxel_volume = sommets_par_voxel / densité;
    double voxel_taille = std::cbrt(voxel_volume);

    if (voxel_taille < min_size) voxel_taille = min_size;
    if (voxel_taille > max_size) voxel_taille = max_size;

    m_voxel_taille = voxel_taille;
}

void grille_voxel::construire_grille()
{
    m_grille.clear();

    for (const auto& sommet : m_points)
    {
        voxel_clé clé = calcule_clé(sommet);

        auto [it, inserted] = m_grille.try_emplace(clé);
        it->second.add(sommet);
    }

    for (auto& [clé, data] : m_grille)
    {
        data.calculer_barycentre();
        data.calculer_covariance();
        data.calculer_couleur(clé);
    }
}

grille_voxel::grille_voxel(const std::filesystem::path &fichier,
                           const double voxel_taille,
                           const double sommets_par_voxel,
                           const double min_size,
                           const double max_size)
    : m_fichier(fichier)
{
    if (!exists(m_fichier))
        throw std::runtime_error("Fichier introuvable : " + m_fichier.string());

    charger_points();

    if (voxel_taille <= 0.0)
        estimer_voxel_size(sommets_par_voxel, min_size, max_size);
    else
        m_voxel_taille = voxel_taille;

    construire_grille();
}


voxel_clé grille_voxel::calcule_clé(const Eigen::Vector3d &sommet) const noexcept
{
    const double inv = 1.0 / m_voxel_taille;
    const int ix = static_cast<int>(std::floor(sommet.x() * inv));
    const int iy = static_cast<int>(std::floor(sommet.y() * inv));
    const int iz = static_cast<int>(std::floor(sommet.z() * inv));

    voxel_clé clé;
    clé.index = Eigen::Vector3i(ix, iy, iz);
    return clé;
}

