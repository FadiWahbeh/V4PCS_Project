#include "Geometry/VoxelGrid.h"
#include <cmath>

namespace V4PCS {
    namespace Geometry {

        VoxelGrid::VoxelGrid(double voxel_size)
            : voxel_size_(voxel_size) {
        }

        void VoxelGrid::build(const PointCloud& cloud) {
            voxels_.clear();
            if (cloud.empty()) return;

            // Calculer les index de grille pour chaque point
            for (const auto& pt : cloud) {
                Eigen::Vector3i idx = getGridIndex(pt.position);

                auto it = voxels_.find(idx);
                if (it == voxels_.end()) {
                    Voxel v;
                    v.grid_index = idx;
                    it = voxels_.emplace(idx, std::move(v)).first;
                }
                it->second.points.push_back(pt);
            }

            // Calculer le barycentre pour chaque voxel
            for (auto& kv : voxels_) {
                computeCentroid(kv.second);
            }
        }

        PointCloud VoxelGrid::getVoxelCentroids() const {
            PointCloud centroids;
            centroids.reserve(voxels_.size());

            for (const auto& kv : voxels_) {
                const Eigen::Vector3i& idx = kv.first;
                const Voxel& voxel = kv.second;

                Point p;
                p.position = voxel.centroid;

                std::uint8_t r, g, b;
                computeColorFromIndex(idx, r, g, b);
                p.r = r;
                p.g = g;
                p.b = b;

                centroids.push_back(p);
            }
            return centroids;
        }

        const std::unordered_map<Eigen::Vector3i, Voxel, Vector3iHash>&
            VoxelGrid::getVoxels() const {
            return voxels_;
        }

        Eigen::Vector3i VoxelGrid::getGridIndex(const Vector3& pt) const {
            return Eigen::Vector3i(
                static_cast<int>(std::floor(pt.x() / voxel_size_)),
                static_cast<int>(std::floor(pt.y() / voxel_size_)),
                static_cast<int>(std::floor(pt.z() / voxel_size_))
            );
        }

        void VoxelGrid::computeCentroid(Voxel& v) const {
            Vector3 sum = Vector3::Zero();
            for (const auto& p : v.points) {
                sum += p.position;
            }
            // Calcul du barycentre
            v.centroid = sum / static_cast<double>(v.points.size());
        }

        void VoxelGrid::computeColorFromIndex(const Eigen::Vector3i& idx,
            std::uint8_t& r,
            std::uint8_t& g,
            std::uint8_t& b) {
            Vector3iHash hasher;
            std::size_t h = hasher(idx);

            r = static_cast<std::uint8_t>(h & 0xFF);
            g = static_cast<std::uint8_t>((h >> 8) & 0xFF);
            b = static_cast<std::uint8_t>((h >> 16) & 0xFF);

            // pour éviter des couleurs trop sombres
            if (r < 30 && g < 30 && b < 30) {
                r = static_cast<std::uint8_t>(r + 50);
                g = static_cast<std::uint8_t>(g + 50);
                b = static_cast<std::uint8_t>(b + 50);
            }
        }

        std::vector<Voxel> VoxelGrid::getVoxelsAsVector() const {
            std::vector<Voxel> list;
            list.reserve(voxels_.size());
            for (const auto& kv : voxels_) {
                list.push_back(kv.second);
            }
            return list;
        }

    }
}
