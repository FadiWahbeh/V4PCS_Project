//
// Created by Enzo Gallet on 09/12/2025.
//

#ifndef GRILLE_VOXEL_H
#define GRILLE_VOXEL_H

#include <cstdint>
#include <functional>
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <cctype>
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

            return h1 ^ (h2 * 1315423911u) ^ (h3 * 2654435761u);
        }
    };
};

struct voxel_data
{
    Eigen::Vector3d somme_position = Eigen::Vector3d::Zero();
    Eigen::Matrix3d somme_produit_externe = Eigen::Matrix3d::Zero();
    Eigen::Vector3d somme_couleur = Eigen::Vector3d::Zero();

    std::uint64_t compteur = 0;

    Eigen::Vector3d barycentre = Eigen::Vector3d::Zero();
    Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();

    std::uint8_t r = 200;
    std::uint8_t g = 200;
    std::uint8_t b = 200;

    void add(const Eigen::Vector3d& s, const Eigen::Vector3d& c) noexcept
    {
        somme_position += s;
        somme_couleur += c;

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
            double div = static_cast<double>(compteur);
            barycentre = somme_position / div;

            Eigen::Vector3d moy_col = somme_couleur / div;

            auto clamp255 = [](double v) -> std::uint8_t {
                if (v < 0.0) v = 0.0;
                if (v > 255.0) v = 255.0;
                return static_cast<std::uint8_t>(v);
            };

            r = clamp255(moy_col.x());
            g = clamp255(moy_col.y());
            b = clamp255(moy_col.z());
        }
    }

    void calculer_covariance() noexcept
    {
        const double N = static_cast<double>(compteur);
        if (N <= 1e-12) { covariance.setZero(); return; }
        const auto E_XXt = somme_produit_externe / N;
        const auto mu_muT = barycentre * barycentre.transpose();
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

    double get_voxel_taille() const {
        return m_voxel_taille;
    }

    const std::unordered_map<voxel_clé, voxel_data, voxel_clé::Hash>& get_grille() const noexcept
    {
        return m_grille;
    }

    static inline bool is_comment_or_empty(const std::string& s) {
        for (char c : s) {
            if (!std::isspace(static_cast<unsigned char>(c))) {
                return (c == '#' || c == '/' || c == ';');
            }
        }
        return true;
    }

private:
    std::filesystem::path m_fichier;
    double m_voxel_taille = -1.0;

    std::unordered_map<voxel_clé, voxel_data, voxel_clé::Hash> m_grille;

    // Optimisations
    // 1- passe si voxel_size connu
    void build_grid_streaming();
    void scan_bbox_and_count(Eigen::Vector3d& mn, Eigen::Vector3d& mx, std::uint64_t& count);

    void estimer_voxel_size(double sommets_par_voxel, double min_size, double max_size);

    voxel_clé calcule_clé(const Eigen::Vector3d &sommet) const noexcept;

    // 2- Parser rapide : lit xyz + garde last3 comme couleur si >=6 colonnes
    static bool parse_xyz_and_color(const std::string& line, Eigen::Vector3d& pos, Eigen::Vector3d& col);
};

#endif