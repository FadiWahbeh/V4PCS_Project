//
// Created by Enzo Gallet on 09/12/2025.
//

#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <chrono>

#include "fichier_manager.h"
#include "grille_voxel.h"
#include "patch_builder.h"
#include "patch_pair_generator.h"
#include "voxel_planaire_extraction.h"
#include "v4pcs_algo.h"

namespace fs = std::filesystem;

struct ResultatNuage {
    std::string nom;
    std::vector<patch_planaire> patchs;
    std::vector<patch_pair> paires;
};

ResultatNuage traiter_nuage(const fs::path& chemin_fichier, const fs::path& dossier_sortie, const std::string& suffixe_debug) {
    ResultatNuage res;
    res.nom = chemin_fichier.stem().string();
    fs::path chemin_base = dossier_sortie / (chemin_fichier.filename().string() + suffixe_debug);

    std::cout << "\n------------------------------------------------" << std::endl;
    std::cout << "Traitement : " << res.nom << " (" << suffixe_debug << ")" << std::endl;

    auto t_start = std::chrono::high_resolution_clock::now();
    double forced_vs = 0.05;
    auto grille = grille_voxel(chemin_fichier, forced_vs);
    auto t_load = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d_load = t_load - t_start;
    std::cout << " [Time] Voxelisation (Lecture fichier) : " << d_load.count() << " sec" << std::endl;

    {
        std::cout << " -> Export des Voxels Bruts (Couleurs)..." << std::endl;
        std::ofstream ofs(chemin_base.string() + ".raw_voxels.ply");
        auto raw_grid = grille.get_grille();

        ofs << "ply\nformat ascii 1.0\n"
            << "element vertex " << raw_grid.size() << "\n"
            << "property float x\nproperty float y\nproperty float z\n"
            << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
            << "end_header\n";

        for(const auto& [k, v] : raw_grid) {
            ofs << v.barycentre.x() << " " << v.barycentre.y() << " " << v.barycentre.z() << " "
                << (int)v.r << " " << (int)v.g << " " << (int)v.b << "\n";
        }
    }

    double lin_thresh = 0.4;
    double curv_thresh = 0.35;
    voxel_planaire_extraction extraction(lin_thresh, curv_thresh);
    extraction.extraire(grille);

    auto g_simple = grille.get_grille();
    auto g_plan = extraction.get_grille_planaire();
    std::cout << "Voxels plans extraits : " << g_plan.size() << std::endl;

    double angle_deg = 15.0;
    double dist_thr = forced_vs * 3.0;
    patch_builder builder(angle_deg, dist_thr);
    builder.construire_patches(g_plan, g_simple);

    builder.exporter_patches_ply(chemin_base);
    std::cout << "Patchs construits : " << builder.get_patches_taille() << std::endl;

    std::vector<patch_planaire> tous = builder.get_patches();
    std::sort(tous.begin(), tous.end(), [](const auto& a, const auto& b){
        return a.voxels.size() > b.voxels.size();
    });

    int max_p = 200;
    int min_v = 50;
    for(const auto& p : tous) {
        if(p.voxels.size() >= min_v) res.patchs.push_back(p);
        if(res.patchs.size() >= max_p) break;
    }
    std::cout << "Patchs filtres conserves : " << res.patchs.size() << " (Top " << max_p << ")" << std::endl;

    double pair_dist = 3.0;
    double pair_ang = 15.0;
    patch_pair_generator gen(pair_dist, pair_ang);
    gen.generer_paires(res.patchs);
    res.paires = gen.get_paires();

    gen.exporter_paires_ply(chemin_base);
    std::cout << "Paires valides generees : " << res.paires.size() << std::endl;

    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d_total = t_end - t_start;
    std::cout << "Total preparation : " << d_total.count() << " sec" << std::endl;

    return res;
}

int main() {
    auto start_global = std::chrono::high_resolution_clock::now();

    std::cout << "V4PCS PROJECT MAIN" << std::endl;
    std::cout << "Input  : " << DATA_DIR_INPUT << std::endl;
    std::cout << "Output : " << DATA_DIR_OUTPUT << std::endl;

    auto manager = fichier_manager(DATA_DIR_INPUT);
    auto fichiers = manager.lister_fichiers(".txt");

    fs::path path_target;
    fs::path path_source;

    if (fichiers.empty()) {
        std::cerr << "Aucun fichier trouvé " << std::endl;
        return -1;
    }
    else if (fichiers.size() == 1) {
        std::cout << "\nMode SELF-REGISTRATION" << std::endl;
        path_target = fichiers[0];
        path_source = fichiers[0];
    }
    else {
        std::cout << "\nMode RECALAGE NORMAL" << std::endl;
        path_target = fichiers[0];
        path_source = fichiers[1];
    }

    ResultatNuage target_data = traiter_nuage(path_target, DATA_DIR_OUTPUT, "_TARGET");
    ResultatNuage source_data = traiter_nuage(path_source, DATA_DIR_OUTPUT, "_SOURCE");

    std::cout << "\n----------------------------------" << std::endl;
    std::cout << "Lancement de l'algorithme V4PCS..." << std::endl;

    auto start_algo = std::chrono::high_resolution_clock::now();

    if (source_data.paires.empty() || target_data.paires.empty()) {
        std::cerr << "Pas assez de paires." << std::endl;
        return -1;
    }

    v4pcs_algo algo(0.20, 5.0);

    Eigen::Matrix4d transformation = algo.aligner(
        source_data.paires,
        target_data.paires,
        source_data.patchs,
        target_data.patchs,
        5000,
        0.20
    );

    auto end_algo = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_algo = end_algo - start_algo;

    std::cout << "\n Temps algorithme RANSAC : " << diff_algo.count() << " sec" << std::endl;

    std::cout << "\n Transformation finale" << std::endl;
    std::cout << transformation << std::endl;

    auto end_global = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_total = end_global - start_global;
    std::cout << "\n[Time] TEMPS TOTAL : " << diff_total.count() << " sec" << std::endl;

    return 0;
}