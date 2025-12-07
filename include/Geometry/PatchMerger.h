#pragma once

#include "Core/Types.h"
#include "Geometry/VoxelGrid.h"   // pour Vector3iHash
#include "../libs/eigen/Eigen/Dense"
#include <vector>
#include <unordered_map>

namespace V4PCS {
    namespace Geometry {

        class PatchMerger {
        public:
            PatchMerger(double angle_thresh_deg, double distance_thresh);

            // Fonction principale
            std::vector<PlanarPatch> merge(
                const std::vector<Plane>& voxel_planes,
                const std::unordered_map<Eigen::Vector3i, int, Vector3iHash>& grid_map
            );

        private:
            double angle_thresh_cos_;
            double distance_thresh_;

            bool isSmoothConnection(const Plane& p1, const Plane& p2) const;
        };

    } // namespace Geometry
} // namespace V4PCS
