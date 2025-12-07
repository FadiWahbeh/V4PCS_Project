#pragma once
#include "Core/Types.h"
#include <unordered_map>
#include <cstdint>

namespace V4PCS {
    namespace Geometry {

        class VoxelGrid {
        public:
            // Construction de la grille
            explicit VoxelGrid(double voxel_size);

            // Exécution de la voxelisation
            void build(const PointCloud& cloud);

            // Récupérer les centroïdes pour la visualisation
            PointCloud getVoxelCentroids() const;

            // Récupérer la map complète
            const std::unordered_map<Eigen::Vector3i, Voxel, Vector3iHash>& getVoxels() const;

            std::vector<Voxel> getVoxelsAsVector() const;

        private:
            double voxel_size_;
            std::unordered_map<Eigen::Vector3i, Voxel, Vector3iHash> voxels_;

            Eigen::Vector3i getGridIndex(const Vector3& pt) const;
            void computeCentroid(Voxel& v) const;

            static void computeColorFromIndex(const Eigen::Vector3i& idx,
                std::uint8_t& r,
                std::uint8_t& g,
                std::uint8_t& b);
        };
    } 
} 
