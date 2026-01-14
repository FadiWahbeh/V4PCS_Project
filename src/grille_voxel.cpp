//
// Created by Enzo Gallet on 09/12/2025.
//

#include "grille_voxel.h"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdlib>

static inline void skip_spaces(const char*& p) {
    while (*p && std::isspace(static_cast<unsigned char>(*p))) ++p;
}

static bool read_double_token(const char*& p, double& out) {
    skip_spaces(p);
    if (!*p) return false;

    char buf[64];
    int i = 0;

    // copie token
    while (*p && !std::isspace(static_cast<unsigned char>(*p))) {
        char c = *p++;
        if (c == ',') c = '.';
        if (i < 63) buf[i++] = c;
    }
    buf[i] = '\0';

    char* end = nullptr;
    out = std::strtod(buf, &end);
    return end && end != buf;
}

bool grille_voxel::parse_xyz_and_color(const std::string& line, Eigen::Vector3d& pos, Eigen::Vector3d& col) {
    if (line.empty() || is_comment_or_empty(line)) return false;

    const char* p = line.c_str();

    double x, y, z;
    if (!read_double_token(p, x)) return false;
    if (!read_double_token(p, y)) return false;
    if (!read_double_token(p, z)) return false;

    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) return false;

    pos = Eigen::Vector3d(x, y, z);
    col = Eigen::Vector3d(200.0, 200.0, 200.0);

    // On lit le reste des colonnes et on garde les 3 dernières
    // total_cols = 3 + rest_cols
    double last3[3] = {200.0, 200.0, 200.0};
    int total_cols = 3;
    double v;

    while (total_cols < 32 && read_double_token(p, v)) {
        // shift last3
        last3[0] = last3[1];
        last3[1] = last3[2];
        last3[2] = v;
        ++total_cols;
    }

    if (total_cols >= 6) {
        col = Eigen::Vector3d(last3[0], last3[1], last3[2]);
    }

    return true;
}

void grille_voxel::scan_bbox_and_count(Eigen::Vector3d& mn, Eigen::Vector3d& mx, std::uint64_t& count)
{
    std::ifstream in(m_fichier, std::ios::in | std::ios::binary);
    if (!in) throw std::runtime_error("Impossible d'ouvrir le fichier : " + m_fichier.string());

    // buffer plus gros
    std::vector<char> buffer(4 * 1024 * 1024);
    in.rdbuf()->pubsetbuf(buffer.data(), static_cast<std::streamsize>(buffer.size()));

    mn = Eigen::Vector3d( std::numeric_limits<double>::infinity(),
                          std::numeric_limits<double>::infinity(),
                          std::numeric_limits<double>::infinity() );
    mx = Eigen::Vector3d(-std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity() );

    count = 0;

    std::string line;
    line.reserve(128);

    Eigen::Vector3d p, c;
    std::size_t skipped = 0;

    while (std::getline(in, line)) {
        if (!parse_xyz_and_color(line, p, c)) { ++skipped; continue; }

        mn = mn.cwiseMin(p);
        mx = mx.cwiseMax(p);
        ++count;
    }

    if (count == 0) {
        throw std::runtime_error("Fichier vide ou illisible : " + m_fichier.string());
    }

    std::cout << "[scan_bbox] " << m_fichier.filename().string()
              << " -> points=" << count
              << ", lignes ignorees=" << skipped << "\n";
}

void grille_voxel::estimer_voxel_size(double sommets_par_voxel, double min_size, double max_size)
{
    Eigen::Vector3d mn, mx;
    std::uint64_t count = 0;
    scan_bbox_and_count(mn, mx, count);

    const double dx = mx.x() - mn.x();
    const double dy = mx.y() - mn.y();
    const double dz = mx.z() - mn.z();

    double volume = dx * dy * dz;
    if (volume <= 0.0) volume = 1e-9;

    const double densite = static_cast<double>(count) / volume;
    double voxel_volume = sommets_par_voxel / densite;
    double voxel_taille = std::cbrt(voxel_volume);

    if (voxel_taille < min_size) voxel_taille = min_size;
    if (voxel_taille > max_size) voxel_taille = max_size;

    m_voxel_taille = voxel_taille;
}

void grille_voxel::build_grid_streaming()
{
    std::ifstream in(m_fichier, std::ios::in | std::ios::binary);
    if (!in) throw std::runtime_error("Impossible d'ouvrir le fichier : " + m_fichier.string());

    std::vector<char> buffer(4 * 1024 * 1024);
    in.rdbuf()->pubsetbuf(buffer.data(), static_cast<std::streamsize>(buffer.size()));

    m_grille.clear();
    // éviter des rehashs multiples
    m_grille.reserve(1'800'000);

    std::string line;
    line.reserve(128);

    std::size_t ok = 0;
    std::size_t skipped = 0;

    Eigen::Vector3d pos, col;

    while (std::getline(in, line)) {
        if (!parse_xyz_and_color(line, pos, col)) { ++skipped; continue; }

        voxel_clé clé = calcule_clé(pos);

        auto [it, inserted] = m_grille.try_emplace(clé);
        it->second.add(pos, col);

        ++ok;
    }

    std::cout << "[charger_points+grille] " << m_fichier.filename().string()
              << " -> points OK=" << ok
              << ", lignes ignorees=" << skipped
              << ", voxels=" << m_grille.size() << "\n";

    if (ok == 0) {
        throw std::runtime_error("Fichier vide ou illisible : " + m_fichier.string());
    }

    for (auto& [clé, data] : m_grille) {
        data.calculer_barycentre();
        data.calculer_covariance();
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

    if (voxel_taille > 0.0) {
        // 1- passe
        m_voxel_taille = voxel_taille;
        build_grid_streaming();
    } else {
        // 2- passes mais sans stocker 40M points
        estimer_voxel_size(sommets_par_voxel, min_size, max_size);
        build_grid_streaming();
    }
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