#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <random>
#include <Eigen/Core>

#include "fichier_manager.h"
#include "grille_voxel.h"
#include "patch_builder.h"
#include "voxel_planaire_extraction.h"
#include "v4pcs_algo.h"

namespace fs = std::filesystem;

struct ColoredPoint {
    Eigen::Vector3d p;
    std::uint8_t r, g, b;
};

// Fonction pour écrire un PLY simple
static void write_ply_colored(const fs::path& path, const std::vector<ColoredPoint>& pts) {
    std::ofstream ofs(path);
    if (!ofs) {
        std::cerr << "Erreur ecriture : " << path.string() << "\n";
        return;
    }
    ofs << "ply\nformat ascii 1.0\n"
        << "element vertex " << pts.size() << "\n"
        << "property float x\nproperty float y\nproperty float z\n"
        << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
        << "end_header\n";
    for (const auto& cp : pts) {
        ofs << (float)cp.p.x() << " " << (float)cp.p.y() << " " << (float)cp.p.z() << " "
            << (int)cp.r << " " << (int)cp.g << " " << (int)cp.b << "\n";
    }
    std::cout << " -> Fichier genere : " << path.filename() << "\n";
}

static inline Eigen::Vector3d apply_T(const Eigen::Matrix4d& T, const Eigen::Vector3d& p) {
    Eigen::Vector4d hp(p.x(), p.y(), p.z(), 1.0);
    return (T * hp).head<3>();
}

// Données intermédiaires
struct ResultatNuage {
    std::string nom;
    std::vector<patch_planaire> patchs;
    std::vector<Eigen::Vector3d> planar_voxels;
};

static ResultatNuage traiter_nuage(const fs::path& chemin_fichier,
                                  const fs::path& dossier_sortie,
                                  const std::string& suffixe_debug)
{
    ResultatNuage res;
    res.nom = chemin_fichier.stem().string();
    fs::path chemin_base = dossier_sortie / (res.nom + suffixe_debug);

    std::cout << "\n--- Traitement : " << res.nom << " ---\n";

    // 1. Voxelisation
    double voxel_size = 0.20;
    std::cout << " -> Voxelisation (grid=" << voxel_size << ")..." << std::flush;
    grille_voxel grille(chemin_fichier, voxel_size);
    std::cout << " OK (" << grille.get_grille().size() << " voxels)\n";

    // 2. Extraction Planarité
    std::cout << " -> Extraction Planarite..." << std::flush;
    voxel_planaire_extraction extraction(0.4, 0.35);
    extraction.extraire(grille);
    const auto& g_plan = extraction.get_grille_planaire();
    std::cout << " OK (" << g_plan.size() << " voxels plans)\n";

    std::vector<ColoredPoint> debug_planar;
    debug_planar.reserve(g_plan.size());
    for (const auto& [k, v] : g_plan) {
        res.planar_voxels.push_back(v.barycentre);
        debug_planar.push_back({v.barycentre, 0, 255, 0});
    }
    write_ply_colored(chemin_base.string() + ".planar_voxels.ply", debug_planar);

    // 3. Construction Patchs
    std::cout << " -> Construction Patchs..." << std::flush;
    patch_builder builder(15.0, voxel_size * 3.0);
    builder.construire_patches(g_plan, grille.get_grille());
    std::cout << " OK (" << builder.get_patches_taille() << " patchs)\n";

    builder.exporter_patches_ply(chemin_base);

    // 4. Filtrage et Sélection
    std::vector<patch_planaire> tous = builder.get_patches();
    std::vector<patch_planaire> valides;
    valides.reserve(tous.size());

    // suppression du bruit
    for(const auto& p : tous) {
        if(p.voxel_keys.size() >= 5) {
            valides.push_back(p);
        }
    }

    std::cout << " -> Melange aleatoire de " << valides.size() << " patchs..." << std::endl;
    std::shuffle(valides.begin(), valides.end(), std::mt19937(std::random_device{}()));

    // On garde 1200 patchs pour donner un maximum de chances à l'algo
    int max_p = 1200;
    for (const auto& p : valides) {
        res.patchs.push_back(p);
        if (res.patchs.size() >= static_cast<size_t>(max_p)) break;
    }
    std::cout << " -> Patchs Conserves pour l'algo : " << res.patchs.size() << "\n";

    return res;
}

int main()
{
    std::cout << "=== V4PCS PROJECT : RECALAGE REEL (FULL DEBUG) ===\n";

    const std::string INPUT_DIR = "../data/input";
    const std::string OUTPUT_DIR = "../data/output";

    if (!fs::exists(INPUT_DIR)) {
        fs::create_directories(INPUT_DIR);
        std::cerr << "Dossier input cree. Ajoutez vos fichiers.\n";
        return -1;
    }
    if (!fs::exists(OUTPUT_DIR)) fs::create_directories(OUTPUT_DIR);

    fichier_manager manager(INPUT_DIR);
    auto fichiers = manager.lister_fichiers("");

    if (fichiers.size() < 2) {
        std::cerr << "Erreur: Il faut 2 fichiers dans input.\n";
        return -1;
    }

    std::sort(fichiers.begin(), fichiers.end());

    fs::path path_target = fichiers[0]; // Station 1 (Fixe)
    fs::path path_source = fichiers[1]; // Station 3 (Mobile)

    std::cout << " TARGET (Fixe)   : " << path_target.filename() << "\n";
    std::cout << " SOURCE (Mobile) : " << path_source.filename() << "\n";

    ResultatNuage data_tgt = traiter_nuage(path_target, OUTPUT_DIR, "_TARGET");
    ResultatNuage data_src = traiter_nuage(path_source, OUTPUT_DIR, "_SOURCE");

    if (data_src.patchs.empty() || data_tgt.patchs.empty()) return -1;

    std::cout << "\n--- Lancement V4PCS ---\n";

    // Distance 1.0m : Large tolérance car peu de recouvrement
    // Angle 45.0 deg : Très permissif pour les normales
    v4pcs_algo algo(1.0, 45.0);

    auto t1 = std::chrono::high_resolution_clock::now();

    // 30000 essais ou 300 secondes (5 minutes) pour trouver la bonne combinaison
    Eigen::Matrix4d T = algo.aligner(data_src.patchs, data_tgt.patchs, 30000, 300);

    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "\nTemps V4PCS : " << std::chrono::duration<double>(t2-t1).count() << "s\n";
    std::cout << "Transformation Finale :\n" << T << "\n";

    std::cout << "Generation du fichier ALIGNEMENT...\n";
    std::vector<ColoredPoint> visual;

    // Cible en ROUGE
    for(const auto& p : data_tgt.planar_voxels) visual.push_back({p, 255, 0, 0});

    // Source recalée en CYAN (Bleu-Vert)
    for(const auto& p : data_src.planar_voxels) visual.push_back({apply_T(T, p), 0, 255, 255});

    fs::path out_path = fs::path(OUTPUT_DIR) / ("ALIGNEMENT_" + data_src.nom + "_SUR_" + data_tgt.nom + ".ply");
    write_ply_colored(out_path, visual);

    std::cout << "Termine. Verifiez le dossier output.\n";
    return 0;
}