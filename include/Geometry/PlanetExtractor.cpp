#include "Geometry/PlaneExtractor.h"
#include "../libs/eigen/Eigen/Eigenvalues"

namespace V4PCS {

    PlaneExtractor::PlaneExtractor(double voxel_size,
        double linearity_thresh,
        double curvature_thresh)
        : voxel_size_(voxel_size),
        linearity_thresh_(linearity_thresh),
        curvature_thresh_(curvature_thresh) {
    }

    // Analyse un voxel : PCA + test de planéité
    bool PlaneExtractor::analyzeVoxel(const PointCloud& voxel_points, Plane& out_plane) {

        // Il faut au moins 3 points pour définir un plan
        if (voxel_points.size() < 3) return false;

        // 1) Barycentre
        Vector3 mean = Vector3::Zero();
        for (const auto& p : voxel_points) {
            mean += p.position;
        }
        mean /= static_cast<double>(voxel_points.size());

        // 2) Matrice de covariance
        Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
        for (const auto& p : voxel_points) {
            Eigen::Vector3d d = p.position - mean;
            cov += d * d.transpose();
        }
        cov /= static_cast<double>(voxel_points.size());

        // 3) Décomposition en valeurs propres
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
        if (solver.info() != Eigen::Success) return false;

        // SelfAdjointEigenSolver renvoie les valeurs propres triées par ordre croissant
        Eigen::Vector3d evals = solver.eigenvalues();
        double l0 = evals(0);   // plus petite
        double l1 = evals(1);
        double l2 = evals(2);   // plus grande

        double sum = l0 + l1 + l2;
        if (sum <= 0.0) return false;

        // mesure de courbure (classique en PCA)
        double curvature = l0 / sum;

        // (optionnel) mesure de "linéarité / planéité"
        // on pourrait ajouter un critère sur (l1 - l0) / l2 ou autre.
        // pour l'instant on ne se sert que de curvature_thresh_
        if (curvature > curvature_thresh_) {
            return false;   // pas assez planaire
        }

        // La normale du plan : vecteur propre associé à la plus petite valeur propre
        Eigen::Vector3d normal = solver.eigenvectors().col(0);
        normal.normalize();

        out_plane.centroid = mean;
        out_plane.normal = normal;
        // grid_index et id_voxel seront remplis dans extractPlanes()
        out_plane.grid_index = Eigen::Vector3i::Zero();
        out_plane.id_voxel = -1;

        return true;
    }

    // Itère sur les voxels et garde ceux qui sont suffisamment planaires
    std::vector<Plane> PlaneExtractor::extractPlanes(const std::vector<Voxel>& voxels) {

        std::vector<Plane> planes;
        planes.reserve(voxels.size());

        for (const auto& v : voxels) {
            Plane pl;
            if (analyzeVoxel(v.points, pl)) {
                // On recopie les infos de voxel
                pl.grid_index = v.grid_index;
                pl.id_voxel = v.id;
                planes.push_back(pl);
            }
        }

        return planes;
    }

} // namespace V4PCS
