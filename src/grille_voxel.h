//
// Created by Enzo Gallet on 09/12/2025.
//

#ifndef GRILLE_VOXEL_H
#define GRILLE_VOXEL_H
#include <cstdint>
#include <functional>
#include <filesystem>
#include <iostream>
#include <Eigen/Core>


struct voxel_clé
{
    Eigen::Vector3i index;

    bool operator==(const voxel_clé& other) const noexcept
    {
        return index == other.index;
    }

    struct Hash
    {
        std::size_t operator()(const voxel_clé& k) const noexcept
        {
            const int ix = k.index.x();
            const int iy = k.index.y();
            const int iz = k.index.z();

            const std::size_t h1 = std::hash<int>{}(ix);
            const std::size_t h2 = std::hash<int>{}(iy);
            const std::size_t h3 = std::hash<int>{}(iz);

            // Combinaison des trois hachages pour former la clé finale
            return h1 ^ (h2 * 1315423911u) ^ (h3 * 2654435761u);
        }
    };
};

struct voxel_data
{
    Eigen::Vector3d somme_position = Eigen::Vector3d::Zero();
    Eigen::Matrix3d somme_produit_externe = Eigen::Matrix3d::Zero();
    std::uint64_t compteur = 0;
    Eigen::Vector3d barycentre = Eigen::Vector3d::Zero();
    Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();

    std::uint8_t r = 0;
    std::uint8_t g = 0;
    std::uint8_t b = 0;

    void add(const Eigen::Vector3d& s) noexcept
    {
        somme_position += s;

        // s*s^T explicite pour limiter les temporaires
        const double x = s.x();
        const double y = s.y();
        const double z = s.z();

        somme_produit_externe(0,0) += x*x;
        somme_produit_externe(0,1) += x*y;
        somme_produit_externe(0,2) += x*z;

        somme_produit_externe(1,0) += y*x;
        somme_produit_externe(1,1) += y*y;
        somme_produit_externe(1,2) += y*z;

        somme_produit_externe(2,0) += z*x;
        somme_produit_externe(2,1) += z*y;
        somme_produit_externe(2,2) += z*z;

        ++compteur;
    }


    void calculer_barycentre() noexcept
    {
        if (compteur != 0) {
            barycentre = somme_position / static_cast<double>(compteur);
        }
    }

    void calculer_couleur(const voxel_clé &idx) {
        constexpr voxel_clé::Hash hasher;
        const std::size_t h = hasher(idx);

        r = static_cast<std::uint8_t>(h & 0xFF);
        g = static_cast<std::uint8_t>((h >> 8) & 0xFF);
        b = static_cast<std::uint8_t>((h >> 16) & 0xFF);

        // Ajustement pour éviter les couleurs trop sombres
        if (r < 30 && g < 30 && b < 30) {
            r = static_cast<std::uint8_t>(r + 50);
            g = static_cast<std::uint8_t>(g + 50);
            b = static_cast<std::uint8_t>(b + 50);
        }
    }

    void calculer_covariance() noexcept
    {
        // 1. Calculer le barycentre (déjà fait)

        // 2. Calcul de la matrice de covariance C
        // C = (1/N) * [ (somme(P * P^T)) - N * (barycentre * barycentre^T) ]
        const double N = static_cast<double>(compteur);

        const auto E_XXt = somme_produit_externe / N;
        const auto mu_muT = barycentre * barycentre.transpose();


        // Covariance = E[X X^T] - E[X] E[X]^T
        covariance = E_XXt - mu_muT;

    }
};

class grille_voxel {
public:
    explicit grille_voxel(const std::filesystem::path& fichier,
                          double voxel_taille = -1.0,
                          double sommets_par_voxel = 150.0,
                          double min_size = 0.03,
                          double max_size = 0.30);



    //void exporter_en_ply(const std::filesystem::path& nom_fichier) const;

    double get_voxel_taille() const {
        return  m_voxel_taille;
    }

    const std::unordered_map<voxel_clé, voxel_data, voxel_clé::Hash>& get_grille() const noexcept
    {
        return m_grille;
    }

private:
    std::filesystem::path m_fichier;
    double m_voxel_taille = -1.0;

    std::vector<Eigen::Vector3d> m_points; // nouveau
    std::unordered_map<voxel_clé, voxel_data, voxel_clé::Hash> m_grille;

    void charger_points();
    void estimer_voxel_size(double sommets_par_voxel, double min_size, double max_size);
    void construire_grille();
    voxel_clé calcule_clé(const Eigen::Vector3d &sommet) const noexcept;
};





#endif //GRILLE_VOXEL_H
