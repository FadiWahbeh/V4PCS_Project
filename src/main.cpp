#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <unordered_set>
#include <Eigen/Core>
#include "fichier_manager.h"
#include "grille_voxel.h"
#include "patch_builder.h"
#include "patch_pair_generator.h"
#include "voxel_planaire_extraction.h"
#include "v4pcs_algo.h"

namespace fs = std::filesystem;

// Helpers
template <class Patch>
static std::size_t patch_voxel_count(const Patch& p) {
    if constexpr (requires { p.voxel_keys.size(); }) {
        return p.voxel_keys.size();
    } else {
        return p.voxels.size();
    }
}

struct ColoredPoint {
    Eigen::Vector3d p;
    std::uint8_t r, g, b;
};

static void write_ply_colored(const fs::path& path, const std::vector<ColoredPoint>& pts) {
    std::ofstream ofs(path);
    if (!ofs) {
        std::cerr << "Impossible d'ouvrir : " << path.string() << "\n";
        return;
    }
    ofs << "ply\nformat ascii 1.0\n";
    ofs << "element vertex " << pts.size() << "\n";
    ofs << "property float x\nproperty float y\nproperty float z\n";
    ofs << "property uchar red\nproperty uchar green\nproperty uchar blue\n";
    ofs << "end_header\n";
    for (const auto& cp : pts) {
        ofs << (float)cp.p.x() << " " << (float)cp.p.y() << " " << (float)cp.p.z() << " "
            << (int)cp.r << " " << (int)cp.g << " " << (int)cp.b << "\n";
    }
}

static inline Eigen::Vector3d apply_T(const Eigen::Matrix4d& T, const Eigen::Vector3d& p) {
    Eigen::Vector4d hp(p.x(), p.y(), p.z(), 1.0);
    Eigen::Vector4d out = T * hp;
    return out.head<3>();
}

// Overlap
struct QKey { int x,y,z; };
static inline bool operator==(const QKey& a, const QKey& b){ return a.x==b.x && a.y==b.y && a.z==b.z; }
struct QKeyHash {
    std::size_t operator()(const QKey& k) const {
        std::size_t h1 = std::hash<int>{}(k.x);
        std::size_t h2 = std::hash<int>{}(k.y);
        std::size_t h3 = std::hash<int>{}(k.z);
        return h1 ^ (h2 * 1315423911u) ^ (h3 * 2654435761u);
    }
};
static inline QKey quantize(const Eigen::Vector3d& p, double eps) {
    return QKey{
        (int)llround(p.x() / eps),
        (int)llround(p.y() / eps),
        (int)llround(p.z() / eps)
    };
}

// Data
struct ResultatNuage {
    std::string nom;

    std::vector<patch_planaire> patchs;
    std::vector<patch_pair> paires;

    double voxel_size = 0.05;
    std::vector<Eigen::Vector3d> planar_voxels;
    std::vector<Eigen::Vector3d> patch_centers_all;
};

static ResultatNuage traiter_nuage(const fs::path& chemin_fichier,
                                  const fs::path& dossier_sortie,
                                  const std::string& suffixe_debug)
{
    ResultatNuage res;
    res.nom = chemin_fichier.stem().string();
    fs::path chemin_base = dossier_sortie / (chemin_fichier.filename().string() + suffixe_debug);

    std::cout << "\n------------------------------------------------\n";
    std::cout << "Traitement : " << res.nom << " (" << suffixe_debug << ")\n";

    auto t_start = std::chrono::high_resolution_clock::now();

    // Voxel size
    double voxel_size = 0.05;
    grille_voxel grille(chemin_fichier, voxel_size);
    res.voxel_size = grille.get_voxel_taille();

    auto t_load = std::chrono::high_resolution_clock::now();
    std::cout << "Voxelisation : "
              << std::chrono::duration<double>(t_load - t_start).count() << " sec\n";

    {
        std::cout << " -> Export des Voxels Bruts (Couleurs)...\n";
        std::ofstream ofs(chemin_base.string() + ".raw_voxels.ply");

        const auto& raw_grid = grille.get_grille(); // PAS de copie
        ofs << "ply\nformat ascii 1.0\n"
            << "element vertex " << raw_grid.size() << "\n"
            << "property float x\nproperty float y\nproperty float z\n"
            << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
            << "end_header\n";

        for (const auto& [k, v] : raw_grid) {
            ofs << (float)v.barycentre.x() << " " << (float)v.barycentre.y() << " " << (float)v.barycentre.z() << " "
                << (int)v.r << " " << (int)v.g << " " << (int)v.b << "\n";
        }
    }

    // Extraction voxels planaires
    double lin_thresh  = 0.4;
    double curv_thresh = 0.35;
    voxel_planaire_extraction extraction(lin_thresh, curv_thresh);
    extraction.extraire(grille);

    const auto& g_simple = grille.get_grille();
    const auto& g_plan   = extraction.get_grille_planaire();

    std::cout << "Voxels plans extraits : " << g_plan.size() << "\n";

    // Stock debug planar voxels
    res.planar_voxels.reserve(g_plan.size());
    for (const auto& [k, v] : g_plan) res.planar_voxels.push_back(v.barycentre);

    // Patches
    double angle_deg = 15.0;
    double dist_thr  = voxel_size * 3.0;

    patch_builder builder(angle_deg, dist_thr);
    builder.construire_patches(g_plan, g_simple);
    builder.exporter_patches_ply(chemin_base);

    std::cout << "Patchs construits : " << builder.get_patches_taille() << "\n";

    // Stock debug patch centers
    res.patch_centers_all.reserve(builder.get_patches().size());
    for (const auto& p : builder.get_patches()) res.patch_centers_all.push_back(p.center);

    // Tri + filtre Top patches
    std::vector<patch_planaire> tous = builder.get_patches(); // copie volontaire (tri)
    std::sort(tous.begin(), tous.end(), [](const auto& a, const auto& b){
        return patch_voxel_count(a) > patch_voxel_count(b);
    });

    int max_p = 200;
    int min_v = 50;
    for (const auto& p : tous) {
        if ((int)patch_voxel_count(p) >= min_v) res.patchs.push_back(p);
        if ((int)res.patchs.size() >= max_p) break;
    }
    std::cout << "Patchs filtres conserves : " << res.patchs.size() << " (Top " << max_p << ")\n";

    double pair_dist = 3.0;
    double pair_ang  = 15.0;

    patch_pair_generator gen(pair_dist, pair_ang);
    gen.generer_paires(res.patchs);
    res.paires = gen.get_paires();
    gen.exporter_paires_ply(chemin_base);

    std::cout << "Paires valides generees : " << res.paires.size() << "\n";

    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Total preparation : "
              << std::chrono::duration<double>(t_end - t_start).count() << " sec\n";

    return res;
}

static void export_and_overlap_debug(const ResultatNuage& target,
                                     const ResultatNuage& source,
                                     const Eigen::Matrix4d& T,
                                     const fs::path& out_dir)
{
    // 1- Compare patch centers : TARGET en rouge et SOURCE aligné en bleu
    std::vector<ColoredPoint> compare_centers;
    compare_centers.reserve(target.patch_centers_all.size() + source.patch_centers_all.size());

    for (const auto& p : target.patch_centers_all) compare_centers.push_back({p, 255, 0, 0});
    for (const auto& p : source.patch_centers_all) compare_centers.push_back({apply_T(T, p), 0, 0, 255});

    write_ply_colored(out_dir / "COMPARE_patch_centers_TARGET_RED__SOURCE_BLUE.ply", compare_centers);

    // 2- Export voxels planaires SOURCE alignés
    std::vector<ColoredPoint> src_planar_aligned;
    src_planar_aligned.reserve(source.planar_voxels.size());
    for (const auto& p : source.planar_voxels) src_planar_aligned.push_back({apply_T(T, p), 0, 0, 255});
    write_ply_colored(out_dir / "SOURCE_planar_voxels_ALIGNED_BLUE.ply", src_planar_aligned);

    // 3- Overlap : pourcentatges des voxels planaires SOURCE alignés proches des voxels planaires TARGET
    const double eps = target.voxel_size;
    const int radius = 1;

    std::unordered_set<QKey, QKeyHash> tgt_set;
    tgt_set.reserve(target.planar_voxels.size() * 2);
    for (const auto& p : target.planar_voxels) tgt_set.insert(quantize(p, eps));

    std::size_t ok = 0;
    for (const auto& p : source.planar_voxels) {
        QKey q = quantize(apply_T(T, p), eps);
        bool found = false;
        for (int dx=-radius; dx<=radius && !found; ++dx)
        for (int dy=-radius; dy<=radius && !found; ++dy)
        for (int dz=-radius; dz<=radius && !found; ++dz) {
            if (tgt_set.find(QKey{q.x+dx, q.y+dy, q.z+dz}) != tgt_set.end()) found = true;
        }
        if (found) ++ok;
    }

    const double ratio = source.planar_voxels.empty() ? 0.0 : (double)ok / (double)source.planar_voxels.size();
    std::cout << " planar voxels overlap (eps=" << eps << ", radius=" << radius << ") = "
              << ratio * 100.0 << "%\n";
}

static void generer_rot180_decale(const fs::path& in_path,
                                 const fs::path& out_path,
                                 char axis,
                                 const Eigen::Vector3d& t)
{
    std::ifstream in(in_path, std::ios::in | std::ios::binary);
    if (!in) throw std::runtime_error("Impossible d'ouvrir: " + in_path.string());

    std::ofstream out(out_path, std::ios::out | std::ios::binary);
    if (!out) throw std::runtime_error("Impossible d'ecrire: " + out_path.string());

    std::vector<char> buffer(4 * 1024 * 1024);
    in.rdbuf()->pubsetbuf(buffer.data(), (std::streamsize)buffer.size());

    std::string line;
    line.reserve(128);

    auto skip_spaces_local = [](const char*& p){
        while (*p && std::isspace((unsigned char)*p)) ++p;
    };

    auto read_tok = [&](const char*& p, double& v)->bool{
        skip_spaces_local(p);
        if (!*p) return false;
        char buf[64]; int i=0;
        while (*p && !std::isspace((unsigned char)*p)) {
            char c=*p++;
            if (c==',') c='.';
            if (i<63) buf[i++]=c;
        }
        buf[i]='\0';
        char* end=nullptr;
        v = std::strtod(buf, &end);
        return end && end!=buf;
    };

    std::uint64_t ok=0, skipped=0;

    while (std::getline(in, line)) {
        if (line.empty() || grille_voxel::is_comment_or_empty(line)) {
            out << line << "\n";
            ++skipped;
            continue;
        }

        const char* p = line.c_str();
        double x,y,z;
        if (!read_tok(p,x) || !read_tok(p,y) || !read_tok(p,z) ||
            !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
            out << line << "\n";
            ++skipped;
            continue;
        }

        Eigen::Vector3d P(x,y,z), Q;

        // rotation 180°
        switch (axis) {
            case 'X': case 'x': Q = Eigen::Vector3d( P.x(), -P.y(), -P.z() ); break;
            case 'Y': case 'y': Q = Eigen::Vector3d(-P.x(),  P.y(), -P.z() ); break;
            default:            Q = Eigen::Vector3d(-P.x(), -P.y(),  P.z() ); break;
        }
        Q += t;

        skip_spaces_local(p);
        out.setf(std::ios::fixed);
        out.precision(9);
        out << Q.x() << " " << Q.y() << " " << Q.z();
        if (*p) out << " " << p;
        out << "\n";

        ++ok;
    }

    std::cout << " OK=" << ok << " skipped=" << skipped
              << " -> " << out_path.filename().string() << "\n";
}

int main()
{

    // Attention : mettre à false après 1 run (si n'est pas déjà généré)
    const bool generate_test_file = false;

    if (generate_test_file) {
        fs::path in  = fs::path(DATA_DIR_INPUT) / "birdfountain_station1_xyz_intensity_rgb.txt";
        fs::path out = fs::path(DATA_DIR_INPUT) / "birdfountain_Rot180Decale.txt";

        if (!fs::exists(out)) {
            generer_rot180_decale(in, out, 'Z', Eigen::Vector3d(0.8, 0.2, 0.0));
        }
    }

    auto start_global = std::chrono::high_resolution_clock::now();

    std::cout << "V4PCS PROJECT MAIN\n";
    std::cout << "Input  : " << DATA_DIR_INPUT << "\n";
    std::cout << "Output : " << DATA_DIR_OUTPUT << "\n";

    fichier_manager manager(DATA_DIR_INPUT);
    auto fichiers = manager.lister_fichiers(".txt");

    if (fichiers.empty()) {
        std::cerr << "Aucun fichier trouvé\n";
        return -1;
    }

    auto contains = [](const fs::path& p, const std::string& s) {
        return p.filename().string().find(s) != std::string::npos;
    };

    fs::path path_target, path_source;

    if (fichiers.size() == 1) {
        std::cout << "\nMode SELF-REGISTRATION\n";
        path_target = fichiers[0];
        path_source = fichiers[0];
    } else {
        std::cout << "\nMode RECALAGE NORMAL\n";

        for (const auto& f : fichiers) {
            if (contains(f, "station1_xyz_intensity_rgb")) path_target = f;
            if (contains(f, "Decale"))                    path_source = f;
        }

        // fallback si pas trouvé
        if (path_target.empty() || path_source.empty()) {
            std::sort(fichiers.begin(), fichiers.end());
            path_target = fichiers[0];
            path_source = fichiers[1];
        }
    }

    std::cout << "\nTARGET = " << path_target.filename().string() << "\n";
    std::cout << "SOURCE = " << path_source.filename().string() << "\n";

    const fs::path out_dir = fs::path(DATA_DIR_OUTPUT);

    ResultatNuage target_data = traiter_nuage(path_target, out_dir, "_TARGET");
    ResultatNuage source_data = traiter_nuage(path_source, out_dir, "_SOURCE");

    std::cout << "\n----------------------------------\n";
    std::cout << "Lancement de l'algorithme V4PCS...\n";

    auto start_algo = std::chrono::high_resolution_clock::now();

    if (source_data.paires.empty() || target_data.paires.empty()) {
        std::cerr << "Pas assez de paires.\n";
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
    std::cout << "\n Temps algorithme RANSAC : "
              << std::chrono::duration<double>(end_algo - start_algo).count() << " sec\n";

    std::cout << "\n Transformation finale\n" << transformation << "\n";

    // Debug superposition
    const bool enable_align_debug = true;
    if (enable_align_debug) {
        export_and_overlap_debug(target_data, source_data, transformation, out_dir);
    }

    auto end_global = std::chrono::high_resolution_clock::now();
    std::cout << "\n TEMPS TOTAL : "
              << std::chrono::duration<double>(end_global - start_global).count() << " sec\n";

    return 0;
}