//
//
//

#ifndef PATCH_PAIR_GENERATOR_H
#define PATCH_PAIR_GENERATOR_H

#include <vector>
#include <filesystem>
#include "patch_builder.h" 

struct patch_pair {
    const patch_planaire* p1;
    const patch_planaire* p2;
    
    double distance;    // d in the paper
    double angle_dot;   // normal dot product
};

class patch_pair_generator {
public:
    patch_pair_generator(double min_distance, double min_angle_deg);

    void generer_paires(const std::vector<patch_planaire>& patches);

    const std::vector<patch_pair>& get_paires() const { return m_paires; }

    void exporter_paires_ply(const std::filesystem::path& nom_fichier) const;

private:
    double m_min_distance;
    double m_max_dot_product;

    std::vector<patch_pair> m_paires;
};

#endif