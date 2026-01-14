//
// Created by Fabien WAHBEH on 22/12/2025.
//

#include "patch_pair_generator.h"
#include <cmath>
#include <iostream>
#include <fstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

patch_pair_generator::patch_pair_generator(double min_distance, double min_angle_deg)
    : m_min_distance(min_distance)
{
    m_max_dot_product = std::cos(min_angle_deg * M_PI / 180.0);
}

void patch_pair_generator::generer_paires(const std::vector<patch_planaire>& patches) {
    m_paires.clear();
    const size_t n = patches.size();

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const auto& p1 = patches[i];
            const auto& p2 = patches[j];

            double dist = (p1.center - p2.center).norm();
            if (dist < m_min_distance) {
                continue;
            }

            double dot = std::abs(p1.normal.dot(p2.normal));

            if (dot > m_max_dot_product) {
                continue;
            }

            patch_pair paire;
            paire.p1 = &p1;
            paire.p2 = &p2;
            paire.distance = dist;
            paire.angle_dot = dot;

            m_paires.push_back(paire);
        }
    }
}

void patch_pair_generator::exporter_paires_ply(const std::filesystem::path& chemin_base) const {
    std::string path_str = chemin_base.string() + ".pairs.ply";
    std::ofstream ofs(path_str);

    if (!ofs) {
        std::cerr << "Erreur export paires : " << path_str << std::endl;
        return;
    }
    ofs << "ply\n"
           "format ascii 1.0\n"
           "element vertex " << m_paires.size() * 2 << "\n"
           "property float x\n"
           "property float y\n"
           "property float z\n"
           "property uchar red\n"
           "property uchar green\n"
           "property uchar blue\n"
           "element edge " << m_paires.size() << "\n"
           "property int vertex1\n"
           "property int vertex2\n"
           "end_header\n";

    for (const auto& paire : m_paires) {
        ofs << paire.p1->center.x() << " " << paire.p1->center.y() << " " << paire.p1->center.z()
            << " 255 255 0\n";

        ofs << paire.p2->center.x() << " " << paire.p2->center.y() << " " << paire.p2->center.z()
            << " 0 255 255\n";
    }

    for (size_t i = 0; i < m_paires.size(); ++i) {
        ofs << (i * 2) << " " << (i * 2 + 1) << "\n";
    }

    std::cout << " -> Export Paires VISUEL OK : " << path_str << std::endl;
}