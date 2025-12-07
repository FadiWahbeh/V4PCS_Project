#include "IO/FileManager.h"
#include "Geometry/VoxelGrid.h"
#include "Geometry/PatchMerger.h"
#include "Geometry/PlaneExtractor.h"

#include <iostream>
#include <filesystem>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <limits>
#include <algorithm>
#include <cmath> 

namespace fs = std::filesystem;

bool isSupportedExtension(const fs::path& p) {
    std::string ext = p.extension().string();
    for (auto& c : ext) c = static_cast<char>(std::tolower(c));
    return (ext == ".ply" || ext == ".pcd" || ext == ".txt" || ext == ".xyz");
}

// Essaie de trouver un dossier data/source autour de l'exécutable
bool findDataDirs(fs::path& sourceDir, fs::path& outputDir) {
    fs::path exeDir = fs::current_path();

    std::vector<fs::path> candidateRoots = {
        exeDir,
        exeDir / "..",
        exeDir / ".." / "..",
        exeDir / ".." / ".." / ".."
    };

    for (const auto& root : candidateRoots) {
        fs::path s = root / "data" / "source";
        fs::path o = root / "data" / "output";
        if (fs::exists(s) && fs::is_directory(s)) {
            if (!fs::exists(o)) {
                fs::create_directories(o);
            }
            sourceDir = fs::canonical(s);
            outputDir = fs::canonical(o);
            return true;
        }
    }
    return false;
}

// Estimation automatique d'un voxel_size "raisonnable"

double estimateOptimalVoxelSize(const V4PCS::PointCloud& cloud,
    double targetPointsPerVoxel,
    double minSize,
    double maxSize)
{
    if (cloud.empty()) {
        return 0.1;
    }

    Eigen::Vector3d minPt = cloud[0].position;
    Eigen::Vector3d maxPt = cloud[0].position;

    for (const auto& p : cloud) {
        minPt = minPt.cwiseMin(p.position);
        maxPt = maxPt.cwiseMax(p.position);
    }

    Eigen::Vector3d diag = maxPt - minPt;
    double volume = std::max(1e-9,
        diag.x() * diag.y() * diag.z());

    double density = static_cast<double>(cloud.size()) / volume;
    double voxelVolume = targetPointsPerVoxel / density;
    double voxelSize = std::cbrt(voxelVolume);
    voxelSize = std::clamp(voxelSize, minSize, maxSize);

    return voxelSize;
}

void processOneFile(const fs::path& inputPath,
    const fs::path& outputDir,
    double voxel_size) {

    if (!fs::exists(inputPath)) {
        std::cerr << " Fichier inexistant : " << inputPath << std::endl;
        return;
    }

    if (!isSupportedExtension(inputPath)) {
        std::cerr << "Extension non supportee : " << inputPath << std::endl;
        return;
    }

    std::cout << "\n=============================\n";
    std::cout << "Traitement de : " << inputPath << "\n";

    // Étape 1 : Lecture

    std::cout << "Etape 1 : Lecture (" << inputPath.extension().string() << ")\n";
    V4PCS::PointCloud cloud = V4PCS::IO::FileManager::loadPointCloudAuto(inputPath.string());

    if (cloud.empty()) {
        std::cerr << "[ERREUR] Nuage vide ou illisible : " << inputPath << std::endl;
        return;
    }

    std::cout << "Nombre de points lus : " << cloud.size() << "\n";

    // Choix du voxel_size 

    double voxel_size_local = voxel_size;

    if (voxel_size_local <= 0.0) {
        double targetPointsPerVoxel = 150.0;
        double minSize = 0.03;
        double maxSize = 0.30;

        voxel_size_local = estimateOptimalVoxelSize(cloud,
            targetPointsPerVoxel,
            minSize,
            maxSize);

        std::cout << "[AUTO] voxel_size estime = " << voxel_size_local << " m\n";
    }
    else {
        std::cout << "[MANUEL] voxel_size = " << voxel_size_local << " m\n";
    }

    std::string stem = inputPath.stem().string();
    std::string vs_str = std::to_string(voxel_size_local);

    fs::path outputVoxelPath = outputDir / ("voxelized_vs" + vs_str + "m_" + stem + ".ply");
    fs::path outputPatchCentersPath = outputDir / ("patches_vs" + vs_str + "m_" + stem + ".ply");
    fs::path outputPatchVoxelsPath = outputDir / ("patch_voxels_vs" + vs_str + "m_" + stem + ".ply");

    // Étape 2 : Voxelisation

    std::cout << "Etape 2 : Voxelisation (taille: " << voxel_size_local << " m)\n";
    V4PCS::Geometry::VoxelGrid grid(voxel_size_local);
    grid.build(cloud);

    V4PCS::PointCloud voxel_centroids = grid.getVoxelCentroids();
    std::cout << "voxel_centroids.size() = " << voxel_centroids.size() << "\n";
    std::cout << "Nombre de voxels crees : " << voxel_centroids.size() << "\n";
    std::cout << "Ratio de reduction : "
        << (1.0 - (double)voxel_centroids.size() / cloud.size()) * 100.0
        << " %\n";

    std::cout << "--- Etape 3 : Sauvegarde des voxels ---\n";
    V4PCS::IO::FileManager::savePLY(outputVoxelPath.string(), voxel_centroids);
    std::cout << "[OK] Voxelisation : " << outputVoxelPath << "\n";

    // Étape 4 : Extraction de plans par voxel

    std::cout << "Etape 4 : Extraction des plans par voxel\n";

    std::vector<V4PCS::Voxel> voxels = grid.getVoxelsAsVector();

    double linearity_thresh = 0.8;
    double curvature_thresh = 0.1;

    V4PCS::PlaneExtractor extractor(voxel_size_local, linearity_thresh, curvature_thresh);
    std::vector<V4PCS::Plane> planes = extractor.extractPlanes(voxels);

    std::cout << "Nombre de voxels planaires (plans extraits) : "
        << planes.size() << "\n";

    if (planes.empty()) {
        std::cout << "[INFO] Aucun plan extrait, on s'arrête ici pour ce fichier.\n";
        return;
    }

    // Étape 5 : Construction du grid_map

    std::cout << "Etape 5 : Construction du grid_map pour les plans\n";

    std::unordered_map<Eigen::Vector3i, int, V4PCS::Vector3iHash> grid_map;
    grid_map.reserve(planes.size());

    for (int i = 0; i < static_cast<int>(planes.size()); ++i) {
        grid_map[planes[i].grid_index] = i;
    }

    // Étape 6 : Fusion des plans en patches

    std::cout << "Etape 6 : Fusion des plans en patches\n";

    double angle_thresh_deg = 10.0;
    double distance_thresh = voxel_size_local * 1.5;

    V4PCS::Geometry::PatchMerger merger(angle_thresh_deg, distance_thresh);
    std::vector<V4PCS::PlanarPatch> patches = merger.merge(planes, grid_map);

    std::cout << "Nombre de patches obtenus : " << patches.size() << "\n";

    if (patches.empty()) {
        std::cout << "[INFO] Aucun patch fusionné, pas de PLY patches.\n";
        return;
    }

    // Debug : stats sur la taille des patches

    {
        size_t min_sz = std::numeric_limits<size_t>::max();
        size_t max_sz = 0;
        size_t sum_sz = 0;

        for (const auto& patch : patches) {
            size_t s = patch.component_planes.size();
            min_sz = std::min(min_sz, s);
            max_sz = std::max(max_sz, s);
            sum_sz += s;
        }

        double mean_sz = patches.empty()
            ? 0.0
            : static_cast<double>(sum_sz) / patches.size();

        std::cout << "Patch size (nb de voxels-plan par patch) : "
            << "min=" << min_sz
            << ", max=" << max_sz
            << ", mean=" << mean_sz << "\n";
    }

    // Étape 7 : Nuage de centres de patches

    std::cout << "Etape 7 : Construction du nuage de centres de patches\n";

    V4PCS::PointCloud patchCentersCloud;
    patchCentersCloud.reserve(patches.size());

    for (const auto& patch : patches) {
        V4PCS::Point p;
        p.position = patch.center;
        p.r = patch.r;
        p.g = patch.g;
        p.b = patch.b;
        patchCentersCloud.push_back(p);
    }

    // Étape 7 bis : Nuage de voxels planaires colorés par patch

    std::cout << "Etape 7 bis : Construction du nuage de voxels colorés par patch\n";

    V4PCS::PointCloud patchVoxelsCloud;
    patchVoxelsCloud.reserve(planes.size());

    for (const auto& patch : patches) {
        for (const auto* planePtr : patch.component_planes) {
            if (!planePtr) continue;

            V4PCS::Point p;
            p.position = planePtr->centroid;
            p.r = patch.r;
            p.g = patch.g;
            p.b = patch.b;
            patchVoxelsCloud.push_back(p);
        }
    }

    // Étape 8 : Sauvegardes

    std::cout << "--- Etape 8 : Sauvegarde des patches (centres) ---\n";
    V4PCS::IO::FileManager::savePLY(outputPatchCentersPath.string(), patchCentersCloud);
    std::cout << "[OK] Patches (centres) sauvegardes : " << outputPatchCentersPath << "\n";

    std::cout << "--- Etape 9 : Sauvegarde des voxels colorés par patch ---\n";
    V4PCS::IO::FileManager::savePLY(outputPatchVoxelsPath.string(), patchVoxelsCloud);
    std::cout << "[OK] Voxels par patch sauvegardes : " << outputPatchVoxelsPath << "\n";
}

int main(int argc, char** argv) {

    fs::path sourceDir, outputDir;
    if (!findDataDirs(sourceDir, outputDir)) {
        std::cerr << "[ERREUR] Impossible de trouver un dossier data/source.\n";
        std::cerr << "Je cherche autour de : " << fs::current_path() << "\n";
        std::cerr << "Assure-toi d'avoir une arborescence du type :\n"
            << "  <racine_projet>/data/source\n"
            << "  <racine_projet>/data/output\n";
        return -1;
    }

    std::cout << "Répertoire de travail  : " << fs::current_path() << "\n";
    std::cout << "Dossier source choisi  : " << sourceDir << "\n";
    std::cout << "Dossier output choisi  : " << outputDir << "\n";

    // 1. Lecture du voxel_size

    double voxel_size = -1.0;

    int firstFileArg = 1;
    if (argc >= 2) {

        char* endPtr = nullptr;
        double vs_candidate = std::strtod(argv[1], &endPtr);

        if (endPtr != argv[1] && *endPtr == '\0') {
            voxel_size = vs_candidate;
            firstFileArg = 2;
        }
        else {
            firstFileArg = 1;
        }
    }

    if (voxel_size <= 0.0) {
        std::cout << "voxel_size = AUTO (sera estime par fichier)\n";
    }
    else {
        std::cout << "voxel_size fixe (manuel) = " << voxel_size << " m\n";
    }

    // 2) Traitement des fichiers

    // CAS A : des fichiers sont passés en argument
    if (argc > firstFileArg) {
        for (int i = firstFileArg; i < argc; ++i) {
            fs::path p(argv[i]);

            if (p.is_relative()) {
                p = sourceDir / p;
            }

            processOneFile(p, outputDir, voxel_size);
        }
    }
    // CAS B : aucun fichier, on traite tout ce qu'il y a dans sourceDir
    else {
        bool anyFile = false;
        for (const auto& entry : fs::directory_iterator(sourceDir)) {
            if (!entry.is_regular_file()) continue;
            const fs::path& p = entry.path();
            if (!isSupportedExtension(p)) continue;

            anyFile = true;
            processOneFile(p, outputDir, voxel_size);
        }

        if (!anyFile) {
            std::cerr << "[INFO] Aucun fichier .ply/.pcd/.txt/.xyz trouve dans "
                << sourceDir << "\n";
        }
    }

    std::cout << "\n=== Termine ===\n";
    return 0;
}
